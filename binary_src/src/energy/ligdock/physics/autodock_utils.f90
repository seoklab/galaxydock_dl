!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ligdock/physics/autodock_utils.f90
!
! Description:
!   This module contains subroutines can be shared only between autodock energy
!   related modules.
!-------------------------------------------------------------------------------
MODULE AUTODOCK_UTILS
!-------------------------------------------------------------------------------

use globals
use logger
use mathfunctions, only: cross, v_norm
!
use energy_vars
use autodock_vars

implicit none
save
private

public :: initialize_atdk_para
public :: set_nb_matrix
public :: set_energy_table
public :: get_dist_idx
public :: get_charges
public :: HD_hbond_parameter
public :: NA_hbond_parameter
public :: OA_SA_hbond_parameter
public :: calc_ramp
public :: racc_rdon_H
public :: racc_rdon_N
public :: racc_rdon_O
public :: get_lig_atom_type

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_atdk_para(para)
!-------------------------------------------------------------------------------
type(atdk_eng_para_type), intent(inout) :: para
integer :: i

para%n_atm_types(:) = 0
para%lig_atom_types(:) = 0
para%type_idx(:) = 0
para%type_name(:) = ''
para%is_Hdon(:) = .false.
para%is_Hacc_N(:) = .false.
para%is_Hacc_O(:) = .false.
para%solpar(:) = 0.0d0
para%vol(:) = 0.0d0

do i = 1, max_atdk_atm_type
    para%ljrij(:,i) = 0.0d0
    para%ljeij(:,i) = 0.0d0
    para%hbond_read(:,i) = .false.
    para%hbrij(:,i) = 0.0d0
    para%hbeij(:,i) = 0.0d0
end do

para%solcn(:) = 0.0d0

end subroutine initialize_atdk_para
!-------------------------------------------------------------------------------
subroutine set_nb_matrix(ligand)
!-------------------------------------------------------------------------------
! nb_matrix(i_atm, j_atm) = 0 ; this pair is excluded in internal E calculation
!                           1 ; this pair is included in internal E calculation
!                           4 ; interaction E between this pair will be rescaled.
!                               (1-4 pair)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
integer :: atm_no, i_atm, j_atm, partner, ii_atm, jj_atm
integer :: i_bnd, i_br, j_br

allocate(nb_matrix(ligand%n_atm, ligand%n_atm))

nb_matrix(:,:) = 1

! 1. exclude 1-2, 1-3 pair, and mark 1-4 pair
do i_atm = 1, ligand%n_atm
    nb_matrix(i_atm, i_atm) = 0
    atm_no = ii_L(i_atm, ligand%lig_no)

    ! 1-2 pair
    do i_bnd = 1, i_P(atm_no)%n_bnd(1)
        j_atm = i_P(atm_no)%i_bnd(i_bnd,1) - tn%nonligatm
        partner = i_L(2, j_atm)
        nb_matrix(i_atm, partner) = 0
    end do

    ! 1-3 pair
    do i_bnd = 1, i_P(atm_no)%n_bnd(2)
        j_atm = i_P(atm_no)%i_bnd(i_bnd,2) - tn%nonligatm
        partner = i_L(2, j_atm)
        nb_matrix(i_atm, partner) = 0
    end do

    ! 1-4 pair
    do i_bnd = 1, i_P(atm_no)%n_bnd(3)
        j_atm = i_P(atm_no)%i_bnd(i_bnd,3) - tn%nonligatm
        partner = i_L(2, j_atm)
        nb_matrix(i_atm, partner) = 4
    end do
end do

! 2. exclude nonbonded interaction within same piece
do i_atm = 1, ligand%n_atm - 1
    do j_atm = i_atm + 1, ligand%n_atm
        if (ligand%piece(i_atm) == ligand%piece(j_atm)) then
            nb_matrix(i_atm, j_atm) = 0
            nb_matrix(j_atm, i_atm) = 0
        end if
    end do
end do

! 3. exclude nonbonded interaction between bridging atom and bridged branch
do i_br = 2, ligand%n_br
    i_atm = ligand%bridge(1,i_br)
    j_atm = ligand%bridge(2,i_br)

    ! exclude bridge_atom - rigid piece - bridge_atom
    do j_br = 2, ligand%n_br
        ii_atm = ligand%bridge(1, j_br)
        jj_atm = ligand%bridge(2, j_br)

        ! j_atm - (i_atm =piece= ii_atm) - jj_atm
        if (ligand%piece(i_atm) == ligand%piece(ii_atm)) then
            nb_matrix(j_atm, jj_atm) = 0
            nb_matrix(jj_atm, j_atm) = 0

        ! j_atm - (i_atm =piece= jj_atm) - ii_atm
        else if (ligand%piece(i_atm) == ligand%piece(jj_atm)) then
            nb_matrix(j_atm, ii_atm) = 0
            nb_matrix(ii_atm, j_atm) = 0

        ! i_atm - (j_atm =piece= ii_atm) - jj_atm
        else if (ligand%piece(j_atm) == ligand%piece(ii_atm)) then
            nb_matrix(i_atm, jj_atm) = 0
            nb_matrix(jj_atm, i_atm) = 0

        ! i_atm - (j_atm =piece= jj_atm) - ii_atm
        else if (ligand%piece(j_atm) == ligand%piece(jj_atm)) then
            nb_matrix(i_atm, ii_atm) = 0
            nb_matrix(ii_atm, i_atm) = 0
        end if
    end do
    
    do atm_no = 1, ligand%n_atm
        ! exclude i_atm - atoms in same piece to j_atm
        if (ligand%piece(atm_no) == ligand%piece(j_atm)) then
            nb_matrix(i_atm, atm_no) = 0
            nb_matrix(atm_no, i_atm) = 0
        ! exclude j_atm - atoms in same piece to i_atm
        else if (ligand%piece(atm_no) == ligand%piece(i_atm)) then
            nb_matrix(j_atm, atm_no) = 0
            nb_matrix(atm_no, j_atm) = 0
        end if
    end do
end do

end subroutine set_nb_matrix
!-------------------------------------------------------------------------------
subroutine set_energy_table(n_atm_type, n_lig)
!-------------------------------------------------------------------------------
! Generate energy table to speed up grid construction
!-------------------------------------------------------------------------------
integer, intent(in) :: n_atm_type, n_lig

allocate(T_vdw(max_dist_short, n_lig, n_atm_type))
allocate(T_solv(max_dist_short))
allocate(T_diel(max_dist_long))
allocate(int_E_table(3, max_dist_long, n_lig, n_lig))

call set_vdw_table(T_vdw, n_atm_type, n_lig)
call set_solv_table(T_solv)
call set_diel_table(T_diel)
call set_int_E_table(int_E_table) 

end subroutine set_energy_table
!-------------------------------------------------------------------------------
subroutine set_vdw_table(T_vdw, n_atm_type, n_lig)
!-------------------------------------------------------------------------------
! Generate van der waals term table
! index of table : T_vdw(i_pro,i_lig,distance)
! smoothing is done
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: T_vdw(:,:,:)
integer, intent(in) :: n_atm_type, n_lig
real(dp) :: Rij, Eij, cA, cB, r, rA, rB, tmp, soften_w
integer :: expA, expB, dist, i_lig, i_pro, lig_type
real(dp), parameter :: sm_width = smooth/2.0

soften_w = 1.0
if (soften_dock_E) soften_w = 0.0001

do i_pro = 1, n_atm_type
    do i_lig = 1, n_lig
        T_vdw(1,i_lig,i_pro) = max_E * soften_w
        lig_type = atdk_para%lig_atom_types(i_lig)
        if (atdk_para%hbond_read(lig_type,i_pro)) then ! H-bond, 12-10 potential
            Rij = atdk_para%hbrij(lig_type,i_pro)
            Eij = atdk_para%hbeij(lig_type,i_pro) * atdk_hbond_w
            expA = 12
            expB = 10
        else  ! vdW, 12-6 potential
            Rij = atdk_para%ljrij(lig_type,i_pro)
            Eij = atdk_para%ljeij(lig_type,i_pro) * atdk_vdw_w
            expA = 12
            expB = 6
        end if
        
        tmp = Eij/dble(expA - expB)
        cA = tmp * (Rij**expA) * expB
        cB = tmp * (Rij**expB) * expA
        do dist = 2, max_dist_short
            r = dble(dist-1)/divisor

!            ! vdw & h_bond
!            if (r < Rij - sm_width) then
!                r = r + sm_width
!            else if (r > Rij + sm_width) then
!                r = r - sm_width
!            else 
!                r = Rij
!            end if

            rB = r**expB
            rA = r**expA
            T_vdw(dist,i_lig,i_pro) = min(max_E*soften_w, cA/rA - cB/rB)
        end do
    end do
end do

call smooth_vdw(T_vdw, n_atm_type, n_lig, max_dist_short)

end subroutine set_vdw_table
!-------------------------------------------------------------------------------
subroutine smooth_vdw(T_vdw, n_atm_type, n_lig, n_dist)
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: T_vdw(:,:,:)
integer, intent(in) :: n_atm_type, n_lig, n_dist
real(dp) :: sm_T_vdw(n_dist,n_lig,n_atm_type), sm_E, soften_w
integer :: i_pro, i_lig, dist, sm_dist
integer :: start, finish, sm_range

sm_range = int(smooth*divisor/2.0)
soften_w = 1.0
if (soften_dock_E) soften_w = 0.0001

do i_pro = 1, n_atm_type
    do i_lig = 1, n_lig
        do dist = 1, n_dist
            sm_E = max_E*soften_w
            start = max(1,dist-sm_range)
            finish = min(n_dist,dist+sm_range-1)
            do sm_dist = start,finish
                sm_E = min(sm_E,T_vdw(sm_dist,i_lig,i_pro))
            end do
            sm_T_vdw(dist,i_lig,i_pro) = sm_E
        end do
    end do
end do

T_vdw = sm_T_vdw
!-------------------------------------------------------------------------------
end subroutine smooth_vdw
!-------------------------------------------------------------------------------
subroutine set_solv_table(T_solv)
!-------------------------------------------------------------------------------
! Generate solvation term table
! dG = -exp(-(r/sig)^2 /2) * S * V
! r is distance btw protein atom and ligand atom(without H-H solv term)
! sig = 3.6 angstrom
! S is ligand atom solvation parameter from autodock parameter
! V is protein atom solvation volume from pdbqs file
! In here, calculate exponential term and returns in table
!-------------------------------------------------------------------------------
real(dp), intent(out) :: T_solv(:)
real(dp) :: const, rij
integer :: i
const = -1.0d0/(2.0d0*(sol_sig**2))

do i = 1, max_dist_short
    rij = dble(i-1)/divisor
    T_solv(i) = exp((rij ** 2) * const) * atdk_solv_w
end do

end subroutine set_solv_table
!-------------------------------------------------------------------------------
subroutine set_diel_table(T_diel)
!-------------------------------------------------------------------------------
! Generate dielectric constant table
! e(r) = (A + B/(1.0 + rk*exp(lamB*r)))
! and it stores in kcal unit.
! Therefore dG = e'(r)*q(lig)*q(prot)/r
! where e'(r) = 332.0/e(r) * atdk_elec_w
!-------------------------------------------------------------------------------
real(dp), intent(out) :: T_diel(:)
real(dp) :: lam, epsij0, A, B, rk, lamB, e, d
integer :: i

lam = 0.003627d0
epsij0 = 78.4d0
A = -8.5525d0
B = epsij0 - A
rk = 7.7839d0
lamB = -lam * B

T_diel(1) = 1.0

!   using distance dependent diel const
do i = 2,max_dist_long
    d = dble(i-1)/divisor
    e = A + B/(1.0d0 + rk*exp(lamB*d))
    e = 332.0d0/e
    e = e * atdk_elec_w
    T_diel(i) = e
end do

end subroutine set_diel_table
!-------------------------------------------------------------------------------
subroutine set_int_E_table(int_E_table)
!-------------------------------------------------------------------------------
! Read internal energy parameters and construct internal energy tables
! Global variables :
!   to store parameters : intern_para
!   Internal energy tables : int_E_table
!-------------------------------------------------------------------------------
real(dp), intent(out) :: int_E_table(:,:,:,:)
real(dp) :: Rij, Eij, cA, cB, rA, rB, tmp
integer :: type1, type2, i1, i2, dist, exp1, exp2
real(dp) :: r_dist, r_sqr, r, sm_width
real(dp) :: lam, epsij0, A, B, rk, lamB, e

! calc diel_fn
lam = 0.003627d0
epsij0 = 78.4d0
A = -8.5525d0
B = epsij0 - A
rk = 7.7839d0
lamB = -lam * B

sm_width = smooth/2.0

do type1 = 1, atdk_para%n_atm_types(2)
    i1 = atdk_para%lig_atom_types(type1)
    do type2 = 1, atdk_para%n_atm_types(2)
        i2 = atdk_para%lig_atom_types(type2)
        if (atdk_para%hbond_read(i2,i1)) then ! H-bond
            Rij = atdk_para%hbrij(i2,i1)
            Eij = atdk_para%hbeij(i2,i1) * atdk_int_h_w
            exp1 = 12
            exp2 = 10
        else
            Rij = atdk_para%ljrij(i2,i1)
            Eij = atdk_para%ljeij(i2,i1) * atdk_int_v_w
            exp1 = 12
            exp2 = 6 
        end if
        tmp = Eij/(exp1-exp2)
        cA = tmp*(Rij**exp1)*exp2
        cB = tmp*(Rij**exp2)*exp1
        ! For NON-SQRT VERSION
        do dist = 1,max_dist_long
            r_dist = dble(dist)
            r_sqr = r_dist/32.0d0
            r = sqrt(r_sqr)
            
            ! diel_fn
            e = A + B / (1.0d0 + rk * exp(lamB * r))
            e = 332.d0 / e
            int_E_table(1,dist,type2,type1) = (e / r) * atdk_int_e_w
         
            ! sol_fn
            int_E_table(2,dist,type2,type1) = exp(sol_const * r_sqr) * atdk_int_s_w
        
!            ! vdw & h_bond
!            if (r < Rij - sm_width) then
!                r = r + sm_width
!            else if (r > Rij + sm_width) then
!                r = r - sm_width
!            else 
!                r = Rij
!            end if
!   
            rB = r**exp2
            rA = r**exp1
            int_E_table(3,dist,type2,type1) = min(max_E,cA/rA-cB/rB)
        end do
    end do
end do

end subroutine set_int_E_table
!-------------------------------------------------------------------------------
subroutine  get_dist_idx(dist, idx_long, idx_short)
!-------------------------------------------------------------------------------
! find index for T_vdw, T_diel, T_solv
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: dist
integer, intent(out) :: idx_long, idx_short
integer :: idx

idx = int(dist*divisor) + 1
if (idx > max_dist_long) then
    idx_long = max_dist_long
    idx_short = max_dist_short
else if (idx > max_dist_short) then
    idx_long = idx
    idx_short = max_dist_short
else
    idx_long = idx
    idx_short = idx
end if

if (dist < 0.5) dist = 0.50
!-------------------------------------------------------------------------------
end subroutine get_dist_idx
!-------------------------------------------------------------------------------
subroutine get_charges(protein, q_s)
!-------------------------------------------------------------------------------
! get charge information of atoms
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
real(dp), intent(out) :: q_s(:)
integer :: i_atm, res_no, atm_no, ref_res_no
integer :: qq_type

do i_atm = 1, tn%nonligatm
    res_no = i_R(1, i_atm)
    atm_no = i_R(2, i_atm)
    if (res_no > protein%n_res) then
        ref_res_no = protein%hetmol(res_no-protein%n_res)%res_type
    else
        ref_res_no = protein%residue(res_no)%res_type
    end if
    qq_type = ref_res_eng(ref_res_no)%qq_idx(atm_no)
    if (force_field_type(1:5) == 'AMBER') then 
        q_s(i_atm) = eng_para%charge(qq_type) / electron2kcal
    else
        q_s(i_atm) = eng_para%charge(qq_type)
    end if
end do

do i_atm = 1, tn%ligatm
    res_no = i_L(1, i_atm)
    atm_no = i_L(2, i_atm)
    ref_res_no = protein%ligand(res_no)%lig_type
    q_s(i_atm + tn%nonligatm) = ref_lig(ref_res_no)%charge(atm_no)
end do

end subroutine get_charges
!-------------------------------------------------------------------------------
subroutine HD_hbond_parameter(hbond, atm_idx)
!-------------------------------------------------------------------------------
! Set basic parameters(rexp, rvec) of hydrogen atom needed to calculate Hbond term.
!         
!   H        
!    \  rvec(:,1) = r_H(:) - r_stem(:).    If stem == 'N', rexp = 2. 
!     stem                                 Otherwise, rexp = 4.
!
!-------------------------------------------------------------------------------
type(atdk_hbond_type), intent(inout) :: hbond
integer, intent(in) :: atm_idx
integer :: stem_atm, type_idx
character(len=2) :: atm_type
real(dp) :: r_H(3), r_stem(3), dr(3)

! find connected atom
stem_atm = i_P(atm_idx)%i_bnd(1,1)
type_idx = atdk_para_idx(stem_atm)

atm_type = atdk_para%type_name(type_idx)(1:1)

if (atm_type == 'N') then
    hbond%rexp = 2
else
    hbond%rexp = 4
end if

r_H(:) = R(:,atm_idx)
r_stem(:) = R(:,stem_atm)
dr(:) = r_H(:) - r_stem(:)
call v_norm(dr)
hbond%rvec(:,1) = dr(:)

end subroutine HD_hbond_parameter
!-------------------------------------------------------------------------------
subroutine NA_hbond_parameter(hbond, atm_idx)
!-------------------------------------------------------------------------------
! Set basic parameters(rvec) of nitrogen atom needed to calculate Hbond term.
!-------------------------------------------------------------------------------
type(atdk_hbond_type), intent(inout) :: hbond
integer, intent(in) :: atm_idx

if (i_P(atm_idx)%n_bnd(1) == 1) then
    call set_azide_hbond(hbond, atm_idx)
else if (i_P(atm_idx)%n_bnd(1) == 2) then
    call set_nitrogen_two_bonds(hbond, atm_idx)
else if (i_P(atm_idx)%n_bnd(1) == 3) then
    call set_nitrogen_three_bonds(hbond, atm_idx)
end if
!-------------------------------------------------------------------------------
end subroutine NA_hbond_parameter
!-------------------------------------------------------------------------------
subroutine set_azide_hbond(hbond, atm_idx)
!-------------------------------------------------------------------------------
! set hbond parameter for azide: :N = stem - X
!                                  rvec(:,1) = r_N(:) - r_stem(:)
!-------------------------------------------------------------------------------
type(atdk_hbond_type), intent(inout) :: hbond
integer, intent(in) :: atm_idx
integer :: stem_idx
real(dp) :: r_N(3), r_stem(3), dr(3)

r_N(1:3) = R(1:3, atm_idx)

stem_idx = i_P(atm_idx)%i_bnd(1,1)
r_stem(1:3) = R(1:3, stem_idx)

dr(1:3) = r_N(1:3) - r_stem(1:3)
call v_norm(dr)

hbond%rvec(1:3,1) = dr(1:3)

end subroutine set_azide_hbond
!-------------------------------------------------------------------------------
subroutine set_nitrogen_two_bonds(hbond, atm_idx)
!-------------------------------------------------------------------------------
! set hbond parameter for nitrogen with two bonds: X1-N=X2
!    
!                 stem1
!       rvec(:,1) /
!         <--- : N          rvec(:,1) = r_N - r_mid(middle point of stem1,stem2)
!                 \           
!                 stem2
!      
!-------------------------------------------------------------------------------
type(atdk_hbond_type), intent(inout) :: hbond
integer, intent(in) :: atm_idx
integer :: stem_1, stem_2
real(dp), dimension(3) :: r_N, r_stem_1, r_stem_2, r_mid, dr

r_N(1:3) = R(1:3, atm_idx)

stem_1 = i_P(atm_idx)%i_bnd(1,1)
stem_2 = i_P(atm_idx)%i_bnd(2,1)

r_stem_1(1:3) = R(1:3, stem_1)
r_stem_2(1:3) = R(1:3, stem_2)

r_mid(1:3) = (r_stem_1(1:3) + r_stem_2(1:3)) / 2.0d0
dr(1:3) = r_N(1:3) - r_mid(1:3)
call v_norm(dr)

hbond%rvec(1:3,1) = dr(1:3) 

end subroutine set_nitrogen_two_bonds
!-------------------------------------------------------------------------------
subroutine set_nitrogen_three_bonds(hbond, atm_idx)
!-------------------------------------------------------------------------------
! set nitrogen hbond parameter with three bonds: 
!
!                   stem1
!       rvec(:,1)  /  
!          <--- : N - stem2   rvec(:,1) = r_N - r_mid 
!                  \                               
!                   stem3
!
!-------------------------------------------------------------------------------
type(atdk_hbond_type), intent(inout) :: hbond
integer, intent(in) :: atm_idx
integer :: stem_1, stem_2, stem_3
real(dp), dimension(3) :: r_N, r_stem_1, r_stem_2, r_stem_3, r_mid, dr

r_N(1:3) = R(1:3, atm_idx)

stem_1 = i_P(atm_idx)%i_bnd(1,1)
stem_2 = i_P(atm_idx)%i_bnd(2,1)
stem_3 = i_P(atm_idx)%i_bnd(3,1)

r_stem_1(1:3) = R(1:3, stem_1)
r_stem_2(1:3) = R(1:3, stem_2)
r_stem_3(1:3) = R(1:3, stem_3)

r_mid(1:3) = (r_stem_1(1:3) + r_stem_2(1:3) + r_stem_3(1:3)) / 3.0d0

dr(1:3) = r_N(1:3) - r_mid(1:3)
call v_norm(dr)

hbond%rvec(1:3,1) = dr(1:3)

end subroutine set_nitrogen_three_bonds
!-------------------------------------------------------------------------------
subroutine OA_SA_hbond_parameter(hbond, atm_idx)
!-------------------------------------------------------------------------------
! set parameters for directional hydrogen bond acceptor
!-------------------------------------------------------------------------------
type(atdk_hbond_type), intent(inout) :: hbond
integer, intent(in) :: atm_idx

if (i_P(atm_idx)%n_bnd(1) == 1) then ! Carbonyl group
    call set_carbonyl_hbond(hbond, atm_idx)
else if (i_P(atm_idx)%n_bnd(1) == 2) then ! Hydroxyl or Ether
    call set_hydroxyl_ether_hbond(hbond, atm_idx)
end if
!-------------------------------------------------------------------------------
end subroutine OA_SA_hbond_parameter
!-------------------------------------------------------------------------------
subroutine set_carbonyl_hbond(hbond, atm_idx)
!-------------------------------------------------------------------------------
! set hbond paramters for S/O = C - X
!       
!         <--- rvec(:,1)
!       S/O = stem          @ rvec(:,2) = vector perpendicular to this face
!              \   dr(:)           i.e.   rvec(:,1) X dr(:)
!               sec
!
!-------------------------------------------------------------------------------
type(atdk_hbond_type), intent(inout) :: hbond
integer, intent(in) :: atm_idx
integer :: stem_idx, sec_atm_idx
real(dp), dimension(3) :: dr, r_O, r_stem, r_sec, rvec2, rvec1

! calculate r vector 1 (O=C vector)
r_O(1:3) = R(1:3, atm_idx)

stem_idx = i_P(atm_idx)%i_bnd(1,1)
r_stem(1:3) = R(1:3, stem_idx)

rvec1(1:3) = r_O(1:3) - r_stem(1:3)
call v_norm(rvec1)

! calculate r vector 2
sec_atm_idx = i_P(atm_idx)%i_bnd(1,2)  ! 1-3 pair
r_sec(1:3) = R(1:3, sec_atm_idx)
dr(1:3) = r_sec(1:3) - r_stem(1:3)
call v_norm(dr)

call cross(rvec1(1:3), dr(1:3), rvec2(1:3))
call v_norm(rvec2)

hbond%rvec(1:3,1) = rvec1(1:3)
hbond%rvec(1:3,2) = rvec2(1:3)

end subroutine set_carbonyl_hbond
!-------------------------------------------------------------------------------
subroutine set_hydroxyl_ether_hbond(hbond, atm_idx)
!-------------------------------------------------------------------------------
! set hbond parameters for hydroxyl and ether group: X1-O-X2
!-------------------------------------------------------------------------------
type(atdk_hbond_type), intent(inout) :: hbond
integer, intent(in) :: atm_idx
integer :: stem_1, stem_2
real(dp), dimension(3) :: dr, r_O, r_stem_1, r_stem_2, rvec1, rvec2
real(dp) :: rdot
real(dp) :: dr1(3), dr2(3), dist1, dist2

r_O(1:3) = R(1:3, atm_idx)

! Set stem1, stem2
dr1(1:3) = r_O(1:3) - R(1:3, i_P(atm_idx)%i_bnd(1,1))
dr2(1:3) = r_O(1:3) - R(1:3, i_P(atm_idx)%i_bnd(2,1))
dist1 = dot_product(dr1,dr1)
dist2 = dot_product(dr2,dr2)
if (dist1 < dist2) then
    stem_1 = i_P(atm_idx)%i_bnd(1,1)
    stem_2 = i_P(atm_idx)%i_bnd(2,1)
else
    stem_2 = i_P(atm_idx)%i_bnd(1,1)
    stem_1 = i_P(atm_idx)%i_bnd(2,1)
end if
r_stem_1(1:3) = R(1:3, stem_1)
r_stem_2(1:3) = R(1:3, stem_2)

rvec2(1:3) = r_stem_1(1:3) - r_stem_2(1:3)
call v_norm(rvec2)

dr(1:3) = r_O(1:3) - r_stem_2(1:3)
rdot = dot_product(dr(1:3), rvec2(1:3))
rvec1(1:3) = r_O(1:3) - ((rdot * rvec2(1:3)) + r_stem_2(1:3))
call v_norm(rvec1)

hbond%rvec(1:3,1) = rvec1(1:3)
hbond%rvec(1:3,2) = rvec2(1:3)

end subroutine set_hydroxyl_ether_hbond
!-------------------------------------------------------------------------------
subroutine calc_ramp(atm_idx, close_atm, hbond_curr, hbond_close, Hramp)
!-------------------------------------------------------------------------------
! TODO: comment
!-------------------------------------------------------------------------------
integer, intent(in) :: atm_idx, close_atm
type(atdk_hbond_type), intent(in) :: hbond_curr, hbond_close
real(dp), intent(out) :: Hramp
real(dp) :: cos_theta, theta

if (atm_idx == close_atm) then
    Hramp = 1.0d0
else
    cos_theta = dot_product(hbond_curr%rvec(1:3, 1), hbond_close%rvec(1:3, 1))
    cos_theta = min(cos_theta, 1.0d0)
    cos_theta = max(cos_theta, -1.0d0)
    theta = acos(cos_theta)
    Hramp = 0.5d0 - 0.5d0 * cos(theta*120.0d0/90.0d0)
end if

end subroutine calc_ramp
!-------------------------------------------------------------------------------
subroutine racc_rdon_H(hbond_para, d)
!-------------------------------------------------------------------------------
! Calculate hbond factor to consider Hbond directionality of hydrogen atom.
!
!        d(:)
!      X --- H
!             \  hbond_para%rvec(:) = r_H(:) - r_stem(:)
!            stem 
!
! If angle theta between d and rvec is smaller than 90 degree, racc = 0.0
! Otherwise, racc = [-cos(theta)]**2 (when stem is Nitrogen) 
!         or racc = [-cos(theta)]**4 (when stem is other hetero atom)
! If theta -> 180, racc -> 1 (more favorable Hbond direction).
!-------------------------------------------------------------------------------
type(atdk_hbond_type), intent(inout) :: hbond_para
real(dp), intent(in) :: d(:)
real(dp) :: cos_theta

hbond_para%racc = 1.0d0
hbond_para%rdon = 1.0d0

cos_theta = -dot_product(d,hbond_para%rvec(:,1)) ! It's actually -cosine(theta)
if (cos_theta < 0.0d0) then ! theta < 90
    hbond_para%racc = 0.0d0
else
    hbond_para%racc = cos_theta**(hbond_para%rexp) ! theta -> 180, racc -> 1
endif

end subroutine racc_rdon_H
!-------------------------------------------------------------------------------
subroutine racc_rdon_O(hbond_para, d)
!-------------------------------------------------------------------------------
! Calculate hbond factor to consider Hbond directionality of oxygen atom.
! It's very hard to describe this term graphically as a comment...
!-------------------------------------------------------------------------------
type(atdk_hbond_type), intent(inout) :: hbond_para
real(dp), intent(in) :: d(:)
real(dp) :: cos_theta, vecX(3), t0, t1, dp_d(3)

hbond_para%racc = 1.0d0
hbond_para%rdon = 0.0d0

cos_theta = -dot_product(d,hbond_para%rvec(:,1))
t0 = dot_product(d,hbond_para%rvec(:,2))
if (t0 > 1.0) then
    t0 = 1.0
else if (t0 < -1.0) then
    t0 = -1.0
end if

! Get deviation from ideal angle
t0 = pi_half - acos(t0) ! Ideal value of t0 = pi_half

if (cos_theta > 0.0d0) then ! 90 < theta < 180
    dp_d = dble(d)
    call cross(dp_d(:),hbond_para%rvec(:,2),vecX(:)) ! vecX is vector in lone-pair plane
    call v_norm(vecX(:))
    t1 = dot_product(vecX(:),hbond_para%rvec(:,1)) 
    if (t1 > 1.0d0) then
        t1 = 1.0d0
    else if (t1 < -1.0d0) then
        t1 = -1.0d0
    end if
    t1 = -pi_half + acos(t1)
    if (t1 < 0.0d0) then
        t1 = -t1  ! range of t1 is set to 0 ~ 90 deg
    end if 
    ! rdon = 1 when t0 = 0, t1 = 45 deg
    hbond_para%rdon = (0.9d0 + 0.1d0*sin(2.0d0*t1))*cos(t0)
else if (cos_theta >= -0.34202d0) then  ! 0.34202 = cos(100 deg)
    hbond_para%rdon = 562.25d0 * (0.116978d0 - cos_theta**2)**3 * cos(t0)
end if
!-------------------------------------------------------------------------------
end subroutine racc_rdon_O
!-------------------------------------------------------------------------------
subroutine racc_rdon_N(hbond_para, dr)
!-------------------------------------------------------------------------------
! Calculate hbond factor to consider Hbond directionality of nitrogen atom.
! It's very hard to describe this term graphically as a comment...
!-------------------------------------------------------------------------------
type(atdk_hbond_type), intent(inout) :: hbond_para
real(dp), intent(in) :: dr(:)
real(dp) :: cos_theta

hbond_para%racc = 1.0d0
hbond_para%rdon = 1.0d0

cos_theta = -dot_product(dr(1:3), hbond_para%rvec(1:3, 1))

if(cos_theta <= 0) then
    hbond_para%rdon = 0.0d0
else
    hbond_para%rdon = cos_theta * cos_theta
end if
!-------------------------------------------------------------------------------
end subroutine racc_rdon_N
!-------------------------------------------------------------------------------
subroutine get_lig_atom_type(type_idx, lig_type)
!-------------------------------------------------------------------------------
integer, intent(in) :: type_idx
integer, intent(out) :: lig_type
integer :: i_lig

do i_lig = 1, atdk_para%n_atm_types(2)
    if (atdk_para%lig_atom_types(i_lig) == type_idx) then
        lig_type = i_lig
        exit
    end if
end do

end subroutine get_lig_atom_type
!-------------------------------------------------------------------------------
END MODULE AUTODOCK_UTILS
!-------------------------------------------------------------------------------
