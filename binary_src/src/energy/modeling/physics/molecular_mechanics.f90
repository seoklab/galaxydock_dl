!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/physics/molecular_mechanics.f90
!
! Description: Molecular mechanics energy terms
!
! TODO add comment
! Available topology and parameter types : Amber ()
!                                          Charmm ()
!
!-------------------------------------------------------------------------------
MODULE MOLECULAR_MECHANICS
!-------------------------------------------------------------------------------
use globals
use logger, only: log_p
use in_out_utils, only: find_atom_idx
use mathfunctions, only: bound_ang, cosine_potential
use geometry, only: calc_ang_and_grad, calc_tor_and_grad
use energy_vars, only: R, i_R, i_L, ii_R, ii_L, i_P, atm_parm, &
    rdie, die_const_pro, die_const_slv, vdw_r_scale, rattle_type, &
    use_bond, use_elec, use_LJ, use_softsphere, soften_short, &
    vdw_soft_k, LJ_att_w, LJ_rep_w, Roff_sqr, truncate_by_shift, &
    use_vdw_scwrl4
use cutoff_methods

implicit none
save
private

! Variables for bonded energy
integer :: n_ang_no, n_bnd_no     ! no. of bonds & angles

type bond_energy_para
   logical :: disabled
   integer :: ia(2)
   real(dp) :: b0, kb
end type bond_energy_para

type angle_energy_para
   logical :: disabled
   integer :: ang_type
   integer :: ia(4)
   integer :: n_ang_fold
   real(dp) :: w0(6), kw(6), n(6)
end type angle_energy_para

type(bond_energy_para),  allocatable :: bnd_para(:), bnd_para_bck(:)
type(angle_energy_para), allocatable :: ang_para(:), ang_para_bck(:)

! Variables for nonbonded energy
real(dp) :: perm_vac_inv  ! Unit conversion factor for Coulomb interaction
!TODO: merge to some value and use that? (when calculating GB)
real(dp) :: Eqq_fac       ! Scaling factor of Coulomb
real(dp) :: Egb_fac       ! Scaling factor of GB
integer, allocatable :: LJ_type(:), QQ_type(:) ! Atomic type of LJ/Coulomb
logical, allocatable :: backbone_hbpair(:,:)
real(dp), allocatable :: qiqj(:,:) ! charge multiplication of i-th&j-th coulombic types
real(dp) :: LJ_soften_a, LJ_soften_b

public :: perm_vac_inv
public :: Egb_fac
public :: LJ_type
public :: QQ_type
public :: qiqj
public :: initialize_MM
public :: finalize_MM
public :: calc_MM_energy
public :: E_vdw_single
public :: update_bonded_para
public :: restore_bonded_para

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_MM(molecule)
!-------------------------------------------------------------------------------
! Initialize molecular mechanics potential evaluation
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule
integer :: i_atm, atmno, resno, ref_res_no
integer :: n_res, n_het, i_res, i_het, n_link, n_lig, i_lig
type(residue_type), allocatable :: residues(:)
type(ligand_residue_type), allocatable :: ligands(:)
type(link_type), allocatable :: links(:)

n_res = molecule%n_res
n_het = molecule%n_het
n_lig = molecule%n_lig
n_link = molecule%n_link
  
allocate(residues(n_res+n_het))
allocate(ligands(max(1,n_lig)))

! Copy into 'residues' to be recognized by following subroutines
do i_res = 1, molecule%n_res
    residues(i_res) = molecule%residue(i_res)
end do
do i_het = 1, molecule%n_het
    residues(n_res+i_het) = molecule%hetmol(i_het)
end do
do i_lig = 1, molecule%n_lig
    ligands(i_lig) = molecule%ligand(i_lig)
end do

! Set nonbonded parameters
call setup_nonbonded_energy(residues(:), ligands(:), n_res, n_het, n_lig)

! Set bond energy indices
allocate(links(0:n_link))
links(1:n_link) = molecule%link(1:n_link)
call setup_bonded_energy(residues(:), ligands(:), n_res, n_het, n_lig, &
    links(0:n_link), n_link)

! Set atomic mass
do i_atm = 1, tn%nonligatm
    resno = i_R(1,i_atm)
    atmno = i_R(2,i_atm)
    ref_res_no = residues(resno)%res_type
    atm_parm(i_atm)%mass = eng_para%mass(ref_res_eng(ref_res_no)%atm_cls(atmno))     
end do
do i_atm = 1, tn%ligatm
    resno = i_L(1, i_atm)
    atmno = i_L(2, i_atm)
    ref_res_no = ligands(resno)%lig_type
    atm_parm(i_atm+tn%nonligatm)%mass = eng_para%mass(ref_lig(ref_res_no)%atm_cls(atmno))
end do

deallocate(residues)
deallocate(ligands)
deallocate(links)

end subroutine initialize_MM
!-------------------------------------------------------------------------------
subroutine finalize_MM()
!------------------------------------------------------------------------------- 
deallocate(bnd_para)
deallocate(ang_para)
deallocate(LJ_type, QQ_type, qiqj)
deallocate(backbone_hbpair)
  
end subroutine finalize_MM
!-------------------------------------------------------------------------------
subroutine setup_bonded_energy(residues, ligands, n_res, n_het, n_lig, links, n_link)
!-------------------------------------------------------------------------------
! Setup for bonded energy arrays
!-------------------------------------------------------------------------------
integer, intent(in) :: n_res, n_het, n_lig
type(residue_type), intent(in) :: residues(n_res+n_het)
type(ligand_residue_type), intent(in) :: ligands(n_lig)
integer, intent(in) :: n_link
type(link_type), intent(in) :: links(0:n_link)
type(ref_res_eng_type) :: eref
integer :: i, k, kk, atm_no, prm_no
integer :: i_res, i_het, i_lig, ang_no, bnd_no, ang_type
integer :: res_no, ref_res_no, i_link
integer :: i_ang_no, i_bnd_no, imax
character(len=4) :: atom_name
character(len=6) :: atom_error_mode = ''
character(len=4) :: next_resname, prev_resname

! 1. Allocate indexing arrays for bond/angle 
n_ang_no = 0
n_bnd_no = 0
  
! Count no. of bond and angles
do i_res = 1, n_res + n_het
    ref_res_no = residues(i_res)%res_type
    n_bnd_no = n_bnd_no + ref_res_eng(ref_res_no)%n_bnd_E
    n_ang_no = n_ang_no + ref_res_eng(ref_res_no)%n_ang_E
end do
do i_lig = 1, n_lig
    ref_res_no = ligands(i_lig)%lig_type
    n_bnd_no = n_bnd_no + ref_lig(ref_res_no)%n_bnd
    n_ang_no = n_ang_no + ref_lig(ref_res_no)%n_ang
    n_ang_no = n_ang_no + ref_lig(ref_res_no)%n_dih
end do

if (n_link > 0) then
    do i_link = 1, n_link
        ref_res_no = links(i_link)%link_type
        n_bnd_no = n_bnd_no + ref_res_eng(ref_res_no)%n_bnd_E
        n_ang_no = n_ang_no + ref_res_eng(ref_res_no)%n_ang_E
    end do
end if

allocate(bnd_para(n_bnd_no))
allocate(ang_para(n_ang_no))

do i = 1, n_bnd_no
    bnd_para(i)%disabled = .false.
end do
do i = 1, n_ang_no
    ang_para(i)%disabled = .false.
end do

i_bnd_no = 0
i_ang_no = 0

! Standard protein residues
do i_res = 1, n_res
    ref_res_no = residues(i_res)%res_type
    eref = ref_res_eng(ref_res_no)
    !to check prev/next residue type (whether it's PRO or GLY)
    if (i_res == n_res) then
        prev_resname = ref_res(residues(n_res-1)%res_type)%res_name
        next_resname = 'XXXX'  !to exclude
    else if (i_res == 1) then
        prev_resname = 'XXXX'  !to exclude
        next_resname = ref_res(residues(2)%res_type)%res_name
    else
        prev_resname = ref_res(residues(i_res-1)%res_type)%res_name
        next_resname = ref_res(residues(i_res+1)%res_type)%res_name
    end if
    
    do bnd_no = 1, eref%n_bnd_E
        i_bnd_no = i_bnd_no + 1
        prm_no = ref_res_eng(ref_res_no)%bnd_prm(bnd_no)
        bnd_para(i_bnd_no)%ia(1:2) = res_index(i_res)%loc_no(eref%atm_in_bnd_E(1:2,bnd_no))
        bnd_para(i_bnd_no)%kb = eng_para%bnd_para(1, prm_no)
        bnd_para(i_bnd_no)%b0 = eng_para%bnd_para(2, prm_no)
        if (eref%atom_in_bnd_E(2,bnd_no) == '+N  ' .and. &
             (next_resname == 'PRO' .or. next_resname == 'CPRO')) then
           bnd_para(i_bnd_no)%disabled = .true.
       else if (eref%atom_in_bnd_E(2,bnd_no) == '+NP ' .and. &
           next_resname /= 'PRO' .and. next_resname /= 'CPRO') then
           bnd_para(i_bnd_no)%disabled = .true.
       end if
   end do
   
   do ang_no = 1, eref%n_ang_E
       ang_type = eref%ang_E_type(ang_no)
       imax = eref%n_atm_in_ang(ang_no)
       i_ang_no = i_ang_no + 1
       prm_no = ref_res_eng(ref_res_no)%ang_prm(ang_no)
       ang_para(i_ang_no)%ang_type = ref_res_eng(ref_res_no)%ang_E_type(ang_no)
       ang_para(i_ang_no)%ia(1:imax) = res_index(i_res)%loc_no(eref%atm_in_ang_E(1:imax,ang_no))
       
       ang_para(i_ang_no)%n_ang_fold = eng_para%n_ang_fold(prm_no)
       do kk = 1, eng_para%n_ang_fold(prm_no)
           ang_para(i_ang_no)%kw(kk) = eng_para%ang_para(1, kk, prm_no)
           ang_para(i_ang_no)%w0(kk) = eng_para%ang_para(2, kk, prm_no)
           ang_para(i_ang_no)%n(kk)  = eng_para%ang_para(3, kk, prm_no)
       end do
       if (ang_type == 2 .or. ang_type == 3) then
           if (eref%atom_in_ang_E(3,ang_no) == '+N  ' .and. &
               (next_resname == 'PRO' .or. next_resname == 'CPRO')) then
               ang_para(i_ang_no)%disabled = .true.
           else if (eref%atom_in_ang_E(3,ang_no) == '+NP ' .and. &
               next_resname /='PRO' .and. next_resname /= 'CPRO') then
               ang_para(i_ang_no)%disabled = .true.
           end if
       else if (ang_type == 1) then
           if (eref%atom_in_ang_E(4,ang_no) == '+N  ' .and. &
               (next_resname == 'PRO' .or. next_resname == 'CPRO')) then
               ang_para(i_ang_no)%disabled = .true.
           else if (eref%atom_in_ang_E(4,ang_no) == '+NP ' .and. &
               (next_resname /= 'PRO' .and. next_resname /= 'CPRO')) then
               ang_para(i_ang_no)%disabled = .true.
           else if (eref%atom_in_ang_E(1,ang_no) == '-CA ' .and. &
               (prev_resname == 'PRO' .or. prev_resname == 'NPRO' .or. &
               prev_resname == 'GLY' .or. prev_resname == 'NGLY')) then
               ang_para(i_ang_no)%disabled = .true.
           else if (eref%atom_in_ang_E(1,ang_no) == '-CAP' .and. &
               prev_resname /= 'PRO' .and. prev_resname /= 'NPRO') then
               ang_para(i_ang_no)%disabled = .true.
           else if (eref%atom_in_ang_E(1,ang_no) == '-CAG' .and. &
               prev_resname /= 'GLY' .and. prev_resname /= 'NGLY') then
               ang_para(i_ang_no)%disabled = .true.
           end if
       end if
   end do
end do

! Hetero molecules
do i_het = 1, n_het
    ref_res_no = residues(n_res+i_het)%res_type
    eref = ref_res_eng(ref_res_no)
    
    do bnd_no = 1, eref%n_bnd_E
        i_bnd_no = i_bnd_no + 1
        prm_no = ref_res_eng(ref_res_no)%bnd_prm(bnd_no)
        bnd_para(i_bnd_no)%ia(1:2) = res_index(n_res+i_het)%loc_no(eref%atm_in_bnd_E(1:2,bnd_no))
        bnd_para(i_bnd_no)%kb = eng_para%bnd_para(1, prm_no)
        bnd_para(i_bnd_no)%b0 = eng_para%bnd_para(2, prm_no)
    end do

    do ang_no = 1, eref%n_ang_E
        imax = eref%n_atm_in_ang(ang_no)
        i_ang_no = i_ang_no + 1
        
        prm_no = ref_res_eng(ref_res_no)%ang_prm(ang_no)
        ang_para(i_ang_no)%ang_type = ref_res_eng(ref_res_no)%ang_E_type(ang_no)
        ang_para(i_ang_no)%ia(1:imax) = res_index(n_res+i_het)%loc_no(eref%atm_in_ang_E(1:imax,ang_no))
        ang_para(i_ang_no)%n_ang_fold = eng_para%n_ang_fold(prm_no)
        do kk = 1, eng_para%n_ang_fold(prm_no)
            ang_para(i_ang_no)%kw(kk) = eng_para%ang_para(1, kk, prm_no)
            ang_para(i_ang_no)%w0(kk) = eng_para%ang_para(2, kk, prm_no)
            ang_para(i_ang_no)%n(kk)  = eng_para%ang_para(3, kk, prm_no)
        end do
    end do
end do

do i_lig = 1, n_lig
    ref_res_no = ligands(i_lig)%lig_type

    do bnd_no = 1, ref_lig(ref_res_no)%n_bnd
        i_bnd_no = i_bnd_no + 1
        prm_no = ref_lig(ref_res_no)%eng_para(bnd_no, 1)
        bnd_para(i_bnd_no)%ia(1:2) = ii_L(ref_lig(i_lig)%bnd(2:3,bnd_no), i_lig)
        bnd_para(i_bnd_no)%b0 = ref_lig(ref_res_no)%b_len0(bnd_no)
        bnd_para(i_bnd_no)%kb = eng_para%bnd_para(1, prm_no)
    end do

    do ang_no = 1, ref_lig(ref_res_no)%n_ang
        i_ang_no = i_ang_no + 1
        prm_no = ref_lig(ref_res_no)%eng_para(ang_no, 2)
        ang_para(i_ang_no)%ang_type = eng_para%ang_type(prm_no)
        ang_para(i_ang_no)%n_ang_fold = eng_para%n_ang_fold(prm_no)
        ang_para(i_ang_no)%ia(1:3) = ii_L(ref_lig(i_lig)%ang(2:4,ang_no), i_lig)
        ang_para(i_ang_no)%w0(1) = ref_lig(ref_res_no)%b_ang0(ang_no)
        ang_para(i_ang_no)%kw(1) = eng_para%ang_para(1, 1, prm_no)
        ang_para(i_ang_no)%n(1)  = eng_para%ang_para(3, 1, prm_no)
    end do

    do ang_no = 1, ref_lig(ref_res_no)%n_dih
        i_ang_no = i_ang_no + 1
        prm_no = ref_lig(ref_res_no)%eng_para(ang_no, 3)
        ang_para(i_ang_no)%ang_type = eng_para%ang_type(prm_no)
        ang_para(i_ang_no)%n_ang_fold = eng_para%n_ang_fold(prm_no)
        ang_para(i_ang_no)%ia(1:4) = ii_L(ref_lig(i_lig)%dih(2:5,ang_no), i_lig)
        ang_para(i_ang_no)%w0(1) = pi-ref_lig(ref_res_no)%d_ang0(ang_no)
        ang_para(i_ang_no)%kw(1) = eng_para%ang_para(1, 1, prm_no)
        ang_para(i_ang_no)%n(1)  = eng_para%ang_para(3, 1, prm_no)
    end do
end do

if (n_link > 0) then
    do i_link = 1, n_link
        ref_res_no = links(i_link)%link_type
        eref = ref_res_eng(ref_res_no)
        
        do bnd_no = 1, eref%n_bnd_E
            i_bnd_no = i_bnd_no + 1
            prm_no = ref_res_eng(ref_res_no)%bnd_prm(bnd_no)
            bnd_para(i_bnd_no)%kb = eng_para%bnd_para(1, prm_no)
            bnd_para(i_bnd_no)%b0 = eng_para%bnd_para(2, prm_no)
            do i = 1, 2
                read(eref%atom_in_bnd_E(i,bnd_no)(1:1),"(I1)") k
                read(eref%atom_in_bnd_E(i,bnd_no)(2:4),"(A4)") atom_name
                res_no = links(i_link)%link_res_no(k)
                call find_atom_idx(res_no, residues(res_no)%res_type, atom_name, &
                    atm_no, atom_error_mode)
                bnd_para(i_bnd_no)%ia(i) = ii_R(atm_no,res_no)
            end do
        end do
        
        do ang_no = 1, eref%n_ang_E
            imax = eref%n_atm_in_ang(ang_no)
            i_ang_no = i_ang_no + 1
            prm_no = ref_res_eng(ref_res_no)%ang_prm(ang_no)
            ang_para(i_ang_no)%ang_type = ref_res_eng(ref_res_no)%ang_E_type(ang_no)
            ang_para(i_ang_no)%n_ang_fold = eng_para%n_ang_fold(prm_no)
            do kk = 1, eng_para%n_ang_fold(prm_no)
                ang_para(i_ang_no)%kw(kk) = eng_para%ang_para(1, kk, prm_no)
                ang_para(i_ang_no)%w0(kk) = eng_para%ang_para(2, kk, prm_no)
                ang_para(i_ang_no)%n(kk)  = eng_para%ang_para(3, kk, prm_no)
            end do
            do i = 1, imax
                read(eref%atom_in_ang_E(i,ang_no)(1:1),"(I1)") k
                read(eref%atom_in_ang_E(i,ang_no)(2:4),"(A4)") atom_name
                res_no = links(i_link)%link_res_no(k)
                call find_atom_idx(res_no, residues(res_no)%res_type, atom_name, &
                    atm_no, atom_error_mode)
                ang_para(i_ang_no)%ia(i) = ii_R(atm_no, res_no)
            end do
        end do
    end do
end if

end subroutine setup_bonded_energy
!-------------------------------------------------------------------------------
subroutine setup_nonbonded_energy(residues, ligands, n_res, n_het, n_lig)
!-------------------------------------------------------------------------------
! Setup for nonbonded energy 
!-------------------------------------------------------------------------------
integer, intent(in) :: n_res, n_het, n_lig
type(residue_type), intent(in) :: residues(n_res+n_het)
type(ligand_residue_type), intent(in) :: ligands(n_lig)
type(ref_res_eng_type) :: eref
integer :: i, j, i_q, j_q
integer :: res_no, atm_no, ref_res_no
real(dp) :: q_product, s2, qi
! Parameters for LJ softening
real(dp), parameter :: LJ_soften_a1 = 12.0d0*(1.0d0/0.6d0)**13 ! Coeff. for repulsive term 
real(dp), parameter :: LJ_soften_a2 = -12.0d0*(1.0d0/0.6d0)**7 ! Coeff. for attraction term
real(dp), parameter :: LJ_soften_b1 = (1.0d0/0.6d0)**12 ! Coeff. for repulsive term 
real(dp), parameter :: LJ_soften_b2 = -2.0d0/(0.6d0)**6 ! Coeff. for attraction term

! Set conversion factors
if (force_field_type == 'CHARMM') then
    perm_vac_inv = 332.0716d0
else if (force_field_type == 'AMBER') then
    perm_vac_inv = 1.0d0
end if

if (rdie) then
    Eqq_fac = perm_vac_inv
else
    Eqq_fac = perm_vac_inv/die_const_pro
    Egb_fac = perm_vac_inv*(1.0d0/die_const_pro - 1.0d0/die_const_slv)
end if

allocate(LJ_type(tn%atom), QQ_type(tn%atom))
allocate(backbone_hbpair(tn%atom,2))
allocate(qiqj(eng_para%n_q_typ, eng_para%n_q_typ))

! Set auxiliary parameters for fast evaluation
s2 = vdw_r_scale*vdw_r_scale
do i = 1, eng_para%n_atom_cls
    do j = 1, eng_para%n_atom_cls
        eng_para%aux_para(2,j,i) = eng_para%aux_para(2,j,i)*s2
        eng_para%aux_para(4,j,i) = eng_para%aux_para(4,j,i)*s2
    end do
end do

LJ_soften_a = LJ_soften_a1*LJ_rep_w + LJ_soften_a2*LJ_att_w
LJ_soften_b = LJ_soften_a*0.6d0 + LJ_soften_b1*LJ_rep_w + LJ_soften_b2*LJ_att_w

! Set non_zero charges and qi*qj
eng_para%non0q(0) = .false.
do i_q = 1, eng_para%n_q_typ
    if (abs(eng_para%charge(i_q)) < small_real) then
        eng_para%non0q(i_q) = .false.
    else
        eng_para%non0q(i_q) = .true.
    end if
    qi = eng_para%charge(i_q)
    do j_q = i_q, eng_para%n_q_typ
        q_product = qi*eng_para%charge(j_q)
        qiqj(j_q,i_q) = q_product
        qiqj(i_q,j_q) = q_product
    end do
end do

! Fill atomic type indices
backbone_hbpair(:,:) = .false.
do i = 1, tn%nonligatm
    res_no = i_R(1,i)
    atm_no = i_R(2,i)
    ref_res_no = residues(res_no)%res_type
    eref = ref_res_eng(ref_res_no)
    LJ_type(i) = eref%atm_cls(atm_no)
    QQ_type(i) = eref%qq_idx(atm_no)
    atm_parm(i)%vdwradii = eng_para%LJ_para(2,LJ_type(i))
    if (trim(ref_res(ref_res_no)%atom_name(atm_no)) == 'H') then
        backbone_hbpair(i,1) = .true.
    else if (trim(ref_res(ref_res_no)%atom_name(atm_no)) == 'O') then
        backbone_hbpair(i,2) = .true.
    end if
end do
do i = 1, tn%ligatm
    res_no = i_L(1,i)
    atm_no = i_L(2,i)
    ref_res_no = ligands(res_no)%lig_type
    LJ_type(tn%nonligatm+i) = ref_lig(ref_res_no)%atm_cls(atm_no)
    QQ_type(tn%nonligatm+i) = ref_lig(ref_res_no)%q_typ_no(atm_no)
    atm_parm(tn%nonligatm+i)%vdwradii = eng_para%LJ_para(2,LJ_type(atm_no))
end do

end subroutine setup_nonbonded_energy
!-------------------------------------------------------------------------------
subroutine update_bonded_para(protein)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer :: i_bnd_no, i_ang_no, idx1, idx2, i_res, j_res, i_atm, j_atm
integer :: n_indx(3), ang_type, na, ia, ja
real(dp) :: dr(3), dist

if (.not. allocated(bnd_para_bck)) allocate(bnd_para_bck(n_bnd_no))
bnd_para_bck(:) = bnd_para(:)

do i_bnd_no = 1, n_bnd_no
    if (bnd_para(i_bnd_no)%disabled) cycle
    idx1 = bnd_para(i_bnd_no)%ia(1)
    idx2 = bnd_para(i_bnd_no)%ia(2)
    !
    i_res = i_R(1, idx1)
    i_atm = i_R(2, idx1)
    j_res = i_R(1, idx2)
    j_atm = i_R(2, idx2)
    !
    if ((i_res > tn%stdres) .or. (j_res > tn%stdres)) exit
    if (i_res == j_res) cycle
    if (protein%residue(i_res)%is_broken .and. protein%residue(j_res)%is_broken) then
        dr(:) = protein%residue(i_res)%R(:,i_atm) - protein%residue(j_res)%R(:,j_atm)
        dist = sqrt(dot_product(dr,dr))
        bnd_para(i_bnd_no)%b0 = dist
        bnd_para(i_bnd_no)%kb = bnd_para(i_bnd_no)%kb * 0.01d0
    end if
end do

if (.not. allocated(ang_para_bck)) allocate(ang_para_bck(n_ang_no))
ang_para_bck(:) = ang_para(:)
n_indx(:) = (/ 4, 3, 4 /)

do i_ang_no = 1, n_ang_no
    if (ang_para(i_ang_no)%disabled) cycle
    !
    ang_type = ang_para(i_ang_no)%ang_type
    na = n_indx(ang_type)
    !
    i_res = i_R(1, ang_para(i_ang_no)%ia(1))
    if (i_res > tn%stdres) exit
    !
    do ia = 1, na-1
        if (ang_para(i_ang_no)%disabled) cycle
        do ja = ia+1, na
            idx1 = ang_para(i_ang_no)%ia(ia)
            idx2 = ang_para(i_ang_no)%ia(ja)
            !
            i_res = i_R(1, idx1)
            i_atm = i_R(2, idx1)
            j_res = i_R(1, idx2)
            j_atm = i_R(2, idx2)
            !
            if (i_res == j_res) cycle
            if (protein%residue(i_res)%is_broken .and. &
                protein%residue(j_res)%is_broken) then
                ang_para(i_ang_no)%disabled = .true.
                exit
            end if
        end do
    end do
end do

end subroutine update_bonded_para
!-------------------------------------------------------------------------------
subroutine restore_bonded_para()
!-------------------------------------------------------------------------------

bnd_para(:) = bnd_para_bck(:)
ang_para(:) = ang_para_bck(:)

deallocate(bnd_para_bck)
deallocate(ang_para_bck)

end subroutine restore_bonded_para
!-------------------------------------------------------------------------------
subroutine calc_MM_energy(f, g, appl_pair, skip_fixed, calc_g)
!-------------------------------------------------------------------------------
! Main subroutine for MM energy calculation
! Evaluate all the bond/angles and non-bonded for appl_pairs
! skip_fixed: Whether to skip fixed atomic pairs
!             This is useful if atoms are fixed but their interaction still should
!             be evaluated (e.g. rigid body minimization in PPDOCK)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f(1:6), g(3,tn%atom,6)
logical, intent(in) :: calc_g, appl_pair(tn%residue,tn%residue), skip_fixed
real(dp) :: f_tmp(3), g_tmp(3,tn%atom,3)

f(:) = 0.0d0
g(:,:,:) = 0.0d0

! Bonded terms
if (use_bond) then
    !when all bonds are constrained do not calculate
    if (rattle_type /= 2) then
        call bond_energy(f(1), g(:,:,1), calc_g)
    end if
    call angle_energy(f_tmp(1:3), g_tmp(:,:,1:3), calc_g)
    !Swap location
    f(2) = f_tmp(2) ! Bond ang
    f(3) = f_tmp(1) ! Proper torsion
    f(4) = f_tmp(3) ! !Improper torsion
    g(:,:,2) = g_tmp(:,:,2)
    g(:,:,3) = g_tmp(:,:,1)
    g(:,:,4) = g_tmp(:,:,3)
end if
  
! Non-bonded terms (only for residue pairs in appl_pair)
if (use_softsphere .or. use_LJ .or. use_vdw_scwrl4 .or. use_elec) then
    call non_bonded_energy(f(5), f(6), g(:,:,5:6), appl_pair(:,:), &
        skip_fixed, calc_g)
end if

end subroutine calc_MM_energy
!-------------------------------------------------------------------------------
subroutine bond_energy(Ebnd, g, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: Ebnd, g(3,tn%atom)
logical, intent(in) :: calc_g
integer  :: idx1, idx2, bnd_no
real(dp) :: b(3), dEbnd(3)
real(dp) :: b_len, db, b0, kb, kb2, C1

g(:,:) = 0.0d0
Ebnd = 0.0d0
do bnd_no = 1, n_bnd_no
    idx1 = bnd_para(bnd_no)%ia(1)
    idx2 = bnd_para(bnd_no)%ia(2)

    if (bnd_para(bnd_no)%disabled) cycle
    if (atm_parm(idx1)%fixed .and. atm_parm(idx2)%fixed) cycle
    if (rattle_type == 1) then
        if (atm_parm(idx1)%is_H .or. atm_parm(idx2)%is_H) cycle
    end if
     
    ! parameters
    b0 = bnd_para(bnd_no)%b0
    kb = bnd_para(bnd_no)%kb

    b(:) = R(:,idx2) - R(:,idx1)
    b_len = sqrt(dot_product(b(:),b(:)))
    db = b_len - b0
    Ebnd = Ebnd + kb*db*db

    if (calc_g) then
        kb2 = 2.0d0*kb
        C1 = kb2*db/b_len
        dEbnd(:) = C1*b(:)
        g(:,idx1) = g(:,idx1) - dEbnd(:)
        g(:,idx2) = g(:,idx2) + dEbnd(:)
    end if
end do

end subroutine bond_energy
!-------------------------------------------------------------------------------
subroutine angle_energy(Eang, g, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: Eang(3), g(3,tn%atom,3)
logical, intent(in) :: calc_g
logical :: ang_fix
integer :: ang_type, na, ang_no
integer :: i, k, idx, kk
integer :: n_indx(3)
integer :: n
real(dp) :: angle, fint, df_int(3,4), dadr(3,4), ra(3,4)
real(dp) :: dw, kw, tmp, w0, df0, fv(2)

Eang = 0.0d0
g(:,:,:) = 0.0d0
n_indx(:) = (/ 4, 3, 4 /)

do ang_no = 1, n_ang_no
    if (ang_para(ang_no)%disabled) cycle
    
    ! Get index, type
    ang_type = ang_para(ang_no)%ang_type
    na = n_indx(ang_type)
     
    ang_fix = .true.
    do i = 1, na
        ra(:,i) = R(:,ang_para(ang_no)%ia(i))
        if (.not. atm_parm(ang_para(ang_no)%ia(i))%fixed) ang_fix = .false.
    end do
    if (ang_fix) cycle

    ! Calc angle
    if (ang_type == 1 .or. ang_type == 3) then ! torsion angle
        call calc_tor_and_grad(ra(:,1:4), angle, dadr(:,:), calc_g)
    else if (ang_type == 2) then ! bond angle 
        call calc_ang_and_grad(ra(:,1:3), angle, dadr(:,1:3), calc_g)
    end if

    ! Parameter index
    ! Iter over n_ang_fold (usually 1, but sometimes > 1)
    do kk = 1, ang_para(ang_no)%n_ang_fold
        w0 = ang_para(ang_no)%w0(kk)
        kw = ang_para(ang_no)%kw(kk)
        if (force_field_type == 'CHARMM') then
            if (ang_type == 1) then ! proper torsion
                n = ang_para(ang_no)%n(kk)
                fv(:) = cosine_potential(angle, w0-pi, kw, n, calc_g)
                fint = fv(1)
                df0 = fv(2)
            else  ! bnd ang & improper torsion
                if (ang_type == 3) w0 = -w0 ! Convention is opposite
                dw = bound_ang(angle - w0)
                tmp = kw*dw
                fint = tmp*dw
                df0 = 2.0d0*tmp
            end if

        else if (force_field_type == 'AMBER') then
            if (ang_type == 1.or. ang_type == 3) then ! proper & improper torsion
                n = ang_para(ang_no)%n(kk)
                fv(:) = cosine_potential(angle, w0-pi, kw, n, calc_g)
                fint = fv(1)
                df0 = fv(2)
            else  ! bnd ang
                dw = bound_ang(angle - w0)
                tmp = kw*dw
                fint = tmp*dw
                df0 = 2.0d0*tmp
            end if
        end if
        
        ! Derivatives
        if (calc_g) then
            df_int(:,1:na) = df0*dadr(:,1:na)
            do k = 1, na
                idx = ang_para(ang_no)%ia(k)
                g(:,idx,ang_type) = g(:,idx,ang_type) + df_int(:,k)
            end do
        end if

        Eang(ang_type) = Eang(ang_type) + fint
    end do
end do

end subroutine angle_energy
!-------------------------------------------------------------------------------
subroutine non_bonded_energy(Evdw, Eqq, g, appl_pair, skip_fixed, calc_g)
!-------------------------------------------------------------------------------
! This calculates non-bonded energy (VdW & electrostatics) only for
! residue pairs having true values of appl_pair(:,:).
!-------------------------------------------------------------------------------
real(dp), intent(out) :: Evdw, Eqq, g(3,tn%atom,2)
logical, intent(in) :: calc_g, appl_pair(tn%residue,tn%residue), skip_fixed
logical :: lj_hb, i_fixed, j_fixed
logical :: vdw_on, qq_on, non0qi
integer :: i, j, np
integer :: LJi, LJj, QQi, QQj
real(dp) :: Ri(3), Rij(3)
real(dp) :: fvdw, fqq
real(dp) :: r_sqr, r_abs
real(dp) :: df0v, df0q
real(dp) :: ftr, gtr
real(dp) :: q2

Evdw = 0.0d0
Eqq = 0.0d0
g(:,:,:) = 0.0d0
lj_hb = .false.

do i = 1, tn%atom
    i_fixed = atm_parm(i)%fixed
    Ri = R(:,i)
    LJi = LJ_type(i)
    QQi = QQ_type(i)
    non0qi = eng_para%non0q(QQi)

    !TODO
    !if (LJi == 0) cycle to skip calculation when type is not assigned?

    ! 1-4 Interactions
    do np = i_P(i)%pair_end_index(2)+1, i_P(i)%pair_end_index(3)
        j = i_P(i)%i_pair(np)
        vdw_on = (use_LJ .or. use_softsphere .or. use_vdw_scwrl4) .and. appl_pair(i_R(1,i),i_R(1,j))
        qq_on  = use_elec .and. appl_pair(i_R(1,i),i_R(1,j))
        
        j_fixed = atm_parm(j)%fixed
        if (i_fixed .and. j_fixed .and. skip_fixed) cycle
        
        r_abs = i_P(i)%d(np)
        r_sqr = r_abs*r_abs
        Rij = R(:,j) - Ri(:)
        
        if (vdw_on) then
            LJj = LJ_type(j)
            if (use_softsphere) then
                call softsphere(vdw_r_scale*eng_para%vdwsum(LJi,LJj), &
                    r_abs, r_sqr, .true., fvdw, df0v, calc_g)
            else if (use_LJ) then
                lj_hb = ((backbone_hbpair(i,1) .and. backbone_hbpair(j,2)) .or. &
                         (backbone_hbpair(j,1) .and. backbone_hbpair(i,2)))
                call LennardJones(r_abs, r_sqr, LJi, LJj, .true., lj_hb, &
                    fvdw, df0v, calc_g)
            else if (use_vdw_scwrl4) then
                lj_hb = ((backbone_hbpair(i,1) .and. backbone_hbpair(j,2)) .or. &
                         (backbone_hbpair(j,1) .and. backbone_hbpair(i,2)))
                call vdw_scwrl4(r_abs, r_sqr, LJi, LJj, .true., lj_hb, &
                    fvdw, df0v, calc_g)
            end if

            call nonbonded_switch_function(r_sqr, ftr, gtr, calc_g)
            Evdw = Evdw + fvdw*ftr
            
            if (calc_g) then
                if ((.not. i_fixed .or. (.not. skip_fixed))) then
                    g(:,i,1) = g(:,i,1) - fvdw*gtr*Rij(:) - df0v*ftr*Rij(:)/r_sqr
                end if
                if ((.not. j_fixed .or. (.not. skip_fixed))) then
                    g(:,j,1) = g(:,j,1) + fvdw*gtr*Rij(:) + df0v*ftr*Rij(:)/r_sqr
                end if
            end if
        end if

        if (qq_on) then
            QQj = QQ_type(j)
            if (non0qi .and. eng_para%non0q(QQj)) then ! non zero charge
                q2 = qiqj(QQi,QQj)
                call Coulomb(r_abs, r_sqr, q2, .true., fqq, df0q, calc_g)
            else
                cycle
            end if

            if (truncate_by_shift) then
                call nonbonded_shift_function(r_sqr, Roff_sqr, ftr, gtr, calc_g)
            else if (.not. vdw_on) then
                call nonbonded_switch_function(r_sqr, ftr, gtr, calc_g)
            end if

            Eqq = Eqq + fqq*ftr
            
            if (calc_g) then
                if ((.not. i_fixed .or. (.not. skip_fixed))) then
                    g(:,i,2) = g(:,i,2) - fqq*gtr*Rij(:) - df0q*ftr*Rij(:)/r_sqr
                end if
                if ((.not. j_fixed .or. (.not. skip_fixed))) then
                    g(:,j,2) = g(:,j,2) + fqq*gtr*Rij(:) + df0q*ftr*Rij(:)/r_sqr
                end if
            end if
        end if
    end do
    
    ! Non-bonded interactions
    do np = i_P(i)%pair_end_index(3)+1, i_P(i)%n_pair
        j = i_P(i)%i_pair(np)
        vdw_on = (use_LJ .or. use_softsphere .or. use_vdw_scwrl4) .and. appl_pair(i_R(1,i),i_R(1,j))
        qq_on  = use_elec .and. appl_pair(i_R(1,i),i_R(1,j))
        
        j_fixed = atm_parm(j)%fixed
        if (i_fixed .and. j_fixed .and. skip_fixed) cycle
        
        r_abs = i_P(i)%d(np)
        r_sqr = r_abs*r_abs
        Rij = R(:,j) - Ri(:)

        if (vdw_on) then
            LJj = LJ_type(j)
            if (use_softsphere) then
                call softsphere(vdw_r_scale*eng_para%vdwsum(LJi,LJj), &
                    r_abs, r_sqr, .false., fvdw, df0v, calc_g)
            else if (use_LJ) then
                lj_hb = ((backbone_hbpair(i,1) .and. backbone_hbpair(j,2)) .or. &
                         (backbone_hbpair(j,1) .and. backbone_hbpair(i,2)))
                call LennardJones(r_abs, r_sqr, LJi, LJj, .false., lj_hb, &
                    fvdw, df0v, calc_g)
            else if (use_vdw_scwrl4) then
                lj_hb = ((backbone_hbpair(i,1) .and. backbone_hbpair(j,2)) .or. &
                         (backbone_hbpair(j,1) .and. backbone_hbpair(i,2)))
                call vdw_scwrl4(r_abs, r_sqr, LJi, LJj, .false., lj_hb, &
                    fvdw, df0v, calc_g)
            end if

            call nonbonded_switch_function(r_sqr, ftr, gtr, calc_g)
            Evdw = Evdw + fvdw*ftr
            
            if (calc_g) then
                if ((.not. i_fixed .or. (.not. skip_fixed))) then
                    g(:,i,1) = g(:,i,1) - fvdw*gtr*Rij(:) - df0v*ftr*Rij(:)/r_sqr
                end if
                if ((.not. j_fixed .or. (.not. skip_fixed))) then
                    g(:,j,1) = g(:,j,1) + fvdw*gtr*Rij(:) + df0v*ftr*Rij(:)/r_sqr
                end if
            end if
        end if

        if (qq_on) then
            QQj = QQ_type(j)
            if (non0qi .and. eng_para%non0q(QQj)) then ! non zero charge
                q2 = qiqj(QQi,QQj)
                call Coulomb(r_abs, r_sqr, q2, .false., fqq, df0q, calc_g)
            else
                cycle
            end if
            
            if (truncate_by_shift) then
                call nonbonded_shift_function(r_sqr, Roff_sqr, ftr, gtr, calc_g)
            else if (.not. vdw_on) then
                call nonbonded_switch_function(r_sqr, ftr, gtr, calc_g)
            end if
            
            Eqq = Eqq + fqq*ftr

            if (calc_g) then
                if ((.not. i_fixed .or. (.not. skip_fixed))) then
                    g(:,i,2) = g(:,i,2) - fqq*gtr*Rij(:) - df0q*ftr*Rij(:)/r_sqr
               end if
               if ((.not. j_fixed .or. (.not. skip_fixed))) then
                   g(:,j,2) = g(:,j,2) + fqq*gtr*Rij(:) + df0q*ftr*Rij(:)/r_sqr
               end if
           end if
       end if
   end do
end do

end subroutine non_bonded_energy
!-------------------------------------------------------------------------------
subroutine softsphere(r_vdwsum, r_abs, r_sqr, e14, fvdw, df0, calc_g)
!-------------------------------------------------------------------------------
! Calculate softsphere clash score for vdW interaction between i-j pair
!-------------------------------------------------------------------------------
real(dp), intent(in) :: r_abs, r_sqr, r_vdwsum
logical, intent(in) :: calc_g, e14
real(dp), intent(out) :: fvdw, df0
real(dp) :: r_vdwsumsqr, r_diff
  
fvdw = 0.0d0
df0 = 0.0d0

r_vdwsumsqr = r_vdwsum**2
if (r_sqr < r_vdwsumsqr .and. .not. e14) then
    r_diff = r_vdwsum - r_abs
    
    fvdw = vdw_soft_k*r_diff*r_diff ! force constant = 200.0 (Default)
     
    if (calc_g) then
        df0 = -2.0d0*vdw_soft_k*r_diff*r_abs
    end if
end if

end subroutine softsphere
!-------------------------------------------------------------------------------
subroutine vdw_scwrl4(r_abs, r_sqr, LJi, LJj, e14, lj_hb, fvdw, df0, calc_g)
!-------------------------------------------------------------------------------
! Based on the vdw potential used by SCWRL4 but modified in a simple form
! to calculate gradients.
! reference: G. G. Krivov et al., Improved prediction of protein side-chain
! conformations with SCWRL4, Proteins, 2009.
! The original parameters were trained by using Emin values from charmm19.
! So I encourage you to use with charmm19 ff.
!-------------------------------------------------------------------------------
real(dp), intent(in) :: r_abs, r_sqr
integer, intent(in) :: LJi, LJj
logical, intent(in) :: calc_g, e14, lj_hb
real(dp), intent(out) :: fvdw, df0
integer :: eps_idx, sig_idx
real(dp) :: sig2, eps_2
real(dp) :: d_sig_1, d_sig_2
real(dp) :: x, x_2, arg1, arg2

!Range parameters for (d/sig)^2
d_sig_1 = (0.8254d0)**2.0d0
d_sig_2 = 100.0d0/81.0d0

fvdw = 0.0d0
df0 = 0.0d0

! Set index considering if i-j is 1-4 pair
if (e14) then ! if 1-4
    eps_idx = 3
    sig_idx = 4 
else
    eps_idx = 1
    sig_idx = 2
end if

sig2 = eng_para%aux_para(sig_idx,LJj,LJi)
x_2 = r_sqr/sig2 
x = x_2**0.5d0

if (x_2 < d_sig_1) then
    fvdw = 10.0d0
else if (x_2 < d_sig_2) then
    !using cosine function upto (10.0/9.0)
    arg1 = pi/(10.0d0/9.0d0-0.8254d0)
    arg2 = cos(arg1*(x-0.8254d0))
    eps_2 = eng_para%aux_para(eps_idx,LJj,LJi)*0.5d0
    fvdw = (5.0d0+eps_2)*arg2 + 5.0d0 - eps_2
    if (calc_g) then
        df0 = (-5.0d0-eps_2)*sqrt(1.0d0 - arg2**2.0d0)*pi/ &
            (10.0d0/9.0d0-0.8254d0)
    end if
else
    call LennardJones(r_abs, r_sqr, LJi, LJj, e14, lj_hb, fvdw, df0, calc_g)
end if

end subroutine vdw_scwrl4
!-------------------------------------------------------------------------------
subroutine LennardJones(r_abs, r_sqr, LJi, LJj, e14, lj_hb, fvdw, df0, calc_g)
!-------------------------------------------------------------------------------
! Calculate Lennard-Jones potential for vdW interaction between i-j pair
!-------------------------------------------------------------------------------
real(dp), intent(in) :: r_sqr, r_abs
integer, intent(in) :: LJi, LJj
logical, intent(in) :: calc_g, e14, lj_hb
real(dp), intent(out) :: fvdw, df0
integer :: eps_idx, sig_idx
real(dp) :: sig2, epsilon_r6, eps, tmp, r2, r6

fvdw = 0.0d0
df0 = 0.0d0

! Set index considering if i-j is 1-4 pair
if (e14) then ! if 1-4
    eps_idx = 3
    sig_idx = 4 
else
    eps_idx = 1
    sig_idx = 2
end if

sig2 = eng_para%aux_para(sig_idx,LJj,LJi)
  
! Modification by softening for singular point at Origin
if (r_sqr/sig2 < 0.36d0 .and. soften_short) then
    eps = eng_para%aux_para(eps_idx,LJj,LJi)
    
    tmp = eps*LJ_soften_a/sqrt(sig2)
    fvdw = tmp*r_abs - LJ_soften_b*eps
    
    if (calc_g) then
        df0 = tmp*r_abs
    end if
! Regular cases
else
    r2 = sig2/r_sqr 
    r6 = r2*r2*r2
    epsilon_r6 = r6*eng_para%aux_para(eps_idx,LJj,LJi)
    if (lj_hb) then
        fvdw = epsilon_r6*(LJ_att_w*2.0d0 - LJ_rep_w*r2*r2)
        if (calc_g) then
            df0 = -10.0d0*(fvdw - 0.8d0*LJ_att_w*epsilon_r6)
        end if
    else
        fvdw = epsilon_r6*(LJ_att_w*2.0d0 - LJ_rep_w*r6) ! epsilone is negative in general
        if (calc_g) then
            df0 = -12.0d0*(fvdw - LJ_att_w*epsilon_r6)
        end if
    end if
end if
  
end subroutine LennardJones
!-------------------------------------------------------------------------------
subroutine Coulomb(r_abs, r_sqr, q2, e14, fqq, df0, calc_g)
!-------------------------------------------------------------------------------
! Calculate Coulomb potential for electrostatic interaction between i-j pair
!-------------------------------------------------------------------------------
real(dp), intent(in) :: r_abs, r_sqr, q2
real(dp), intent(out) :: fqq, df0
logical, intent(in) :: e14, calc_g
real(dp) :: qq2
real(dp), parameter :: a = 0.85d0

qq2 = q2*Eqq_fac
if (e14) qq2 = qq2*eng_para%E14fac ! Take care if 1-4

! Modification for singular point at Origin
if (r_sqr <= 1.0d0 .and. soften_short) then
    fqq = qq2*(-2.0d0*r_abs + 3.0d0)
    if (calc_g) df0 = - 2.0d0*qq2*r_abs
    
else if (rdie) then ! distance dependent dielectric
    fqq = qq2 / r_sqr
    if (calc_g) df0 = -2.0d0*fqq
    
else
    fqq = qq2 / r_abs
    if (calc_g) df0 = - fqq
     
end if
  
end subroutine Coulomb
!-------------------------------------------------------------------------------
subroutine E_vdw_single(i_atm, j_atm, scale, fvdw, soften, pair_type)
!-------------------------------------------------------------------------------
! TODO
! This subroutine is quite redundant to 'LennardJones'.Can be modified later.
! Lennard-Jones energy for single atomic pair
! Used for checking clash
!-------------------------------------------------------------------------------
integer, intent(in) :: i_atm, j_atm, pair_type
real(dp), intent(in) :: scale
real(dp), intent(out) :: fvdw
logical, intent(in) :: soften
integer :: sig_idx, eps_idx
integer :: LJi, LJj
real(dp) :: r_sqr, r_abs, Rij(3)
real(dp) :: sig2, eps, arg
real(dp) :: r2, r6, epsilon_r6

Rij(:) = R(:,j_atm) - R(:,i_atm)
r_sqr = dot_product(Rij(:), Rij(:))
r_abs = sqrt(r_sqr)
           
LJi = LJ_type(i_atm)
LJj = LJ_type(j_atm)

! Set index considering if i-j is 1-4 pair
if (pair_type == 3) then ! if 1-4
    eps_idx = 3
    sig_idx = 4 
else
    eps_idx = 1
    sig_idx = 2
end if

sig2 = scale*scale*eng_para%aux_para(sig_idx,LJj,LJi)

! Modification for singular point at Origin
if (r_sqr/sig2 < 0.36d0 .and. soften) then
    eps = eng_para%aux_para(eps_idx,LJj,LJi)
    arg = eps*LJ_soften_a/sqrt(sig2)
    fvdw =  arg*r_abs - LJ_soften_b*eps
else
    r2 = sig2/r_sqr 
    r6 = r2*r2*r2
    epsilon_r6 = r6 * eng_para%aux_para(eps_idx,LJj,LJi)
    fvdw = epsilon_r6 * (2.0d0 - r6) ! epsilone is negative
end if

end subroutine E_vdw_single
!-------------------------------------------------------------------------------
END MODULE MOLECULAR_MECHANICS
!-------------------------------------------------------------------------------
