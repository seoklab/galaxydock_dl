!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/knowledge/hbond_energy.f90
!
! Description: This module is for evaluating Rosetta H-bond energy.
!
!  1. Consideration of Solvation state on hydrogen bonding is held by 
!    setting [use_SAHbond] as
!    .false. (default) : no consideration
!    .true. : consider, using parameters defined in Rosetta
! 
!  2. Assignment for Acceptor & Donor and hybridization state can be slightly 
!    different from Rosetta.
!
! Reference : Kortemme, T. et al. An orientation-dependent Hydrogen bonding 
!             potential improves prediction of specificity and structure 
!             for protein and protein-protein complexes, 
!             J. Mol. Biol. 2003, 326, 1239-1259.
!
!-------------------------------------------------------------------------------
MODULE HBOND_ENERGY
!-------------------------------------------------------------------------------

use globals
use logger
!
use mathfunctions, only: fade_function, poly3, poly5, poly7, poly8, &
                         v_norm, sigmoidal1, sigmoidal2, cosine_well
use convert_res_name, only: convert_to_stdres
use geometry,         only: calc_bnd_angle, calc_tor_angle
!
use energy_vars,      only: R, i_R, ii_R, use_coHbond, use_SAHbond,i_P
use in_out_utils,     only: find_atom_idx
implicit none
save
private

real(dp), parameter :: hb_weight(1:16) = (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0, 0.49d0, &
       0.49d0, 0.49d0, 0.49d0, 0.49d0, 0.49d0, 0.49d0,&
       0.49d0, 0.49d0, 0.49d0, 0.49d0, 0.49d0 /) ! Weight more for backbone-backbone Hbond

! rosetta weight parameter.
!  real(dp), parameter :: hb_weight(1:16) = (/ 1.17d0, 1.17d0, 1.17d0, 1.17d0, 1.170d0, &
!       1.17d0, 1.10d0, 1.17d0, 1.10d0, 1.17d0, 1.10d0,&
!       1.17d0, 1.10d0, 1.10d0, 1.10d0, 1.10d0 /) !
  
! Parameters for fade function (fade away for less important geometry)
real(dp), parameter :: min_r = 1.4d0
real(dp), parameter :: max_r = 3.0d0
real(dp), parameter :: min_r2 = min_r*min_r
real(dp), parameter :: max_r2 = max_r*max_r
real(dp), parameter :: r_interp_edge = 2.1d0
real(dp), parameter :: angle_interp_edge = 0.05d0
real(dp), parameter :: interp_range = 0.2d0
real(dp), parameter :: min_xD = 0.0d0
real(dp), parameter :: min_xH = 0.0d0
real(dp), parameter :: max_xD = 1.0d0
real(dp), parameter :: max_xH = 1.0d0
real(dp), parameter :: switch_dis = 2.1d0
real(dp), parameter :: interp_min = switch_dis - interp_range
real(dp), parameter :: interp_max = switch_dis + interp_range
real(dp), parameter :: MAX_HB_ENERGY = 0.0
!Coefficients for H-bond calculation
real(dp), parameter :: c8_AHdisBBhelix(1:8) = (/ 12.93768086d0, -221.0155722d0, &
       1604.391304d0, -6409.335773d0, 15200.86425d0, -21375.00216d0, 16475.98811d0,-5361.556440d0 /)
real(dp), parameter :: c8_AHdisBBOther(1:8) = (/ 13.58980244d0, -224.0452428d0, &
       1568.933094d0, -6044.257847d0, 13820.14980d0, -18730.96076d0, 13912.92238d0,-4361.995425d0 /)
real(dp), parameter :: c5_AHdisSP2(1:5)     = (/ 10.98727738d0, -100.2401419d0, &
       340.9733405d0, -511.6111233d0, 285.0061262d0 /)
real(dp), parameter :: c5_AHdisSP3(1:5)     = (/ 7.011735538d0, -68.99968829d0, &
       251.8209310d0, -403.3593133d0, 238.7378958d0 /)
real(dp), parameter :: c8_xDBBHelix(1:8)    = (/ 223.5268153d0, -757.7254095d0, &
       1019.593508d0, -689.2232431d0, 240.1436064d0, -37.84119583d0, 0.858689040d0, 0.278181985d0 /)
real(dp), parameter :: c8_xDBBOther(1:8)    = (/ 111.9877946d0, -380.3066184d0, &
       514.7650204d0, -352.4092342d0, 124.6219703d0, -19.94401946d0, 0.149314979d0, 0.635771774d0 /)
real(dp), parameter :: c3_xDSP2short(1:3)   = (/-0.562582503d0, -0.746682668d0, 0.809265171d0 /)
real(dp), parameter :: c3_xDSP2long(1:3)    = (/ 0.094962885d0, -0.254313172d0, 0.0d0 /)
real(dp), parameter :: c3_xDSP3short(1:3)   = (/-0.100140144d0, -1.139139041d0, 0.739279186d0 /)
real(dp), parameter :: c3_xDSP3long(1:3)    = (/ 0.089380221d0, -0.207503776d0, 0.0d0 /)
real(dp), parameter :: c8_xHBBHelix(1:8)    = (/ 54.80664331d0, -196.8196655d0, &
       295.9418886d0, -232.1056020d0, 96.99124565d0, -20.60918361d0, 1.573169816d0, 0.000745458d0 /)
real(dp), parameter :: c8_xHBBOther(1:8)    = (/ 43.94483847d0, -144.3836033d0, &
       193.5865176d0, -132.4469355d0, 47.28137288d0, -8.945888012d0,-0.227035135d0, 0.791902995d0 /)
real(dp), parameter :: c3_xHSP2short(1:3)   = (/ 1.720984644d0, -1.855254573d0, 0.0d0 /)
real(dp), parameter :: c3_xHSP2long(1:3)    = (/ 0.439598249d0, -0.444673076d0, 0.0d0 /)
real(dp), parameter :: c3_xHSP3(1:3)        = (/ 1.761487842d0, -1.876959406d0, 0.0d0 /)
real(dp), parameter :: c7_xHRing(1:7)       = (/ 37.74431600d0, -117.7316740d0, &
       143.0759275d0, -86.22588350d0, 26.74481750d0, -4.469970500d0, 0.365845500d0 /)
!Parameters for co-operativity
real(dp), parameter :: dipole_max0 = -0.5d0
real(dp), parameter :: dipole_fmax = -1.0d0
real(dp), parameter :: coop_factor_antiparallel = 0.5d0
real(dp), parameter :: coop_factor_parallel = 0.2d0

real(dp), parameter :: buriedness_exposed      = 0.16d0/0.49d0
real(dp), parameter :: buriedness_intermediate = 0.44d0/0.49d0
real(dp), parameter :: buriedness_buried       = 0.94d0/0.49d0

! Variables
integer :: max_hbpair ! Max available Hbonds
integer :: nH, nA     ! # of Hydrogens & Acceptors
integer :: max_hbADH  ! Maximum number of hydrogen bonds
type hbond_acceptor_type
    integer :: hybrid_state !1: sp2, 2:sp3, 3: ring(EX: ND,NE in HIS)
    integer :: n_B_atm
    integer :: A_atm
    integer :: B_atm(2) !1: B_atm in rosetta hbond, 2: B2_atm in rosetta hbond
    integer :: res
end type hbond_acceptor_type

type hbond_donor_type
    integer :: H_atm
    integer :: D_atm
    integer :: res
end type hbond_donor_type

type hbond_pair_type
    integer :: hbond_type
    real(dp) :: AHdist
    real(dp) :: xD
    real(dp) :: xH
    type(hbond_acceptor_type) :: acc
    type(hbond_donor_type) :: don

end type hbond_pair_type

integer, allocatable :: acc_atom_index(:)   ! tn%atom   -> 0 for disabled
integer, allocatable :: don_atom_index(:)   ! tn%atom   -> 0 for disabled
logical, allocatable :: is_BB(:)   
logical, allocatable :: is_HOH(:)   
type(hbond_acceptor_type), allocatable :: hbond_acc(:)
type(hbond_donor_type), allocatable :: hbond_don(:)
!type(hbond_pair_type) :: hbond_pair(max_hbpair)
type(hbond_pair_type), allocatable :: hbond_pair(:)
public :: initialize_hbond
public :: finalize_hbond
!
public :: calc_hbond_energy
!
public :: calc_hbond_energy_single
!
public :: max_hbpair
public :: get_hbond_pair

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_hbond(protein)
!-------------------------------------------------------------------------------
! Initial setup subroutine of Hbond energy for modeling
! Hbonds will be enumerated only for residue pairs defined at appl_pair
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein

max_hbADH = 3*tn%residue
max_hbpair = 30*tn%residue

allocate(is_BB(tn%atom),is_HOH(tn%atom),acc_atom_index(tn%atom), don_atom_index(tn%atom))
allocate(hbond_acc(max_hbADH), hbond_don(max_hbADH))
allocate(hbond_pair(max_hbpair))

call setup_hbond_misc(protein)

end subroutine initialize_hbond
!-------------------------------------------------------------------------------
subroutine finalize_hbond()
!-------------------------------------------------------------------------------
deallocate(is_BB,is_HOH,acc_atom_index, don_atom_index)
deallocate(hbond_acc, hbond_don)
deallocate(hbond_pair)

end subroutine finalize_hbond
!-------------------------------------------------------------------------------
subroutine setup_hbond_misc(protein)
!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer :: i_atm, resno, atmno, ref_res_no
character(len=4) :: atmname, resname

is_BB(:) = .false.
is_HOH(:) = .false.
acc_atom_index(:) = 0
don_atom_index(:) = 0
nH = 0
nA = 0
  
do i_atm = 1, tn%stdatm
    resno = i_R(1,i_atm)
    atmno = i_R(2,i_atm)
    ref_res_no = protein%residue(resno)%res_type

    ! Correct residue name
    atmname = trim(ref_res(ref_res_no)%atom_name(atmno))
    resname = ref_res(ref_res_no)%res_name
    if ((.not. resname == 'PTR') .and. &
        (.not. resname == 'SEP') .and. &
        (.not. resname == 'TPO')) &
        call convert_to_stdres(resname)

    ! Set arrays about Acceptor / Donor candidates
    call setup_hbond_index(resname, atmname, i_atm, resno, atmno, ref_res_no)
end do
!for water
do i_atm = tn%stdatm + 1 , tn%stdatm + tn%hetatm!for water atm
   resno = i_R(1,i_atm)
   atmno = i_R(2,i_atm)
   ref_res_no = protein%hetmol(resno - protein%n_res)%res_type

   ! Correct residue name
   atmname = trim(ref_res(ref_res_no)%atom_name(atmno))
   resname = protein%hetmol(resno - protein%n_res)%res_name

   ! Set arrays about Acceptor / Donor candidates
   if (resname == 'HOH') then 
      call setup_hbond_index(resname, atmname, i_atm, resno, atmno, ref_res_no)
   endif
enddo

end subroutine setup_hbond_misc
!-------------------------------------------------------------------------------
subroutine setup_hbond_index(restype, atmtype, i_atm, resno, atmno, ref_res_no)
!-------------------------------------------------------------------------------
! Assign 'all available' candidates for acceptor/donor and their types
!-------------------------------------------------------------------------------
character(len=4), intent(in) :: restype, atmtype
integer, intent(in) :: ref_res_no, i_atm, resno, atmno
integer :: ac_in_res, ac2_in_res
! Setup for Hydrogen & Acceptor Indices
! 1. Donor & Hydrogen, consider only polar-hydrogen
if (trim(atmtype) == 'H' .or. &
       trim(atmtype) == 'H1' .or. trim(atmtype) == 'H2' .or. trim(atmtype) == 'H3' .or. &
       (restype(1:3) == 'LYS' .and. atmtype(1:2) == 'HZ') .or. &
       (restype(1:3) == 'ASN' .and. atmtype(1:3) == 'HD2') .or. &
       (restype(1:3) == 'GLN' .and. atmtype(1:3) == 'HE2') .or. &
       (restype(1:3) == 'ARG' .and. atmtype(1:2) == 'HH') .or. &
       (restype(1:3) == 'ARG' .and. atmtype(1:2) == 'HE') .or. &
       (restype(1:3) == 'SER' .and. atmtype(1:2) == 'HG') .or. &
       (restype(1:3) == 'THR' .and. atmtype(1:2) == 'HG') .or. &
       (restype(1:3) == 'TYR' .and. atmtype(1:2) == 'HH') .or. &
       (restype(1:3) == 'TRP' .and. atmtype(1:3) == 'HE1')) then
    nH = nH + 1
    don_atom_index(i_atm) = nH
    hbond_don(nH)%H_atm = i_atm
    hbond_don(nH)%D_atm = ii_R(ref_res(ref_res_no)%atm_in_bnd(1,atmno), resno)
    hbond_don(nH)%res = resno
    if (trim(restype(1:3)) == 'HOH') is_HOH(i_atm) = .true.

else if (trim(atmtype(1:3)) == 'N') then
    is_BB(i_atm) = .true.

! 2. Acceptor & Acceptor-connected-atom, SP2 Hybridization state
! Backbone Nitrogen excluded
else if (trim(atmtype(1:3)) == 'OXT' &
        .or. ((restype(1:3) /= 'HOH') .and. trim(atmtype(1:3)) == 'O') &
           .or. (restype(1:3) == 'ASP' .and. atmtype(1:2) == 'OD') &
           .or. (restype(1:3) == 'GLU' .and. atmtype(1:2) == 'OE') &
           .or. (restype(1:3) == 'ASN' .and. atmtype(1:2) == 'OD') & 
           .or. (restype(1:3) == 'GLN' .and. atmtype(1:2) == 'OE') & 
           .or. (restype(1:3) == 'PTR' .and. atmtype(1:3) == 'O1P') &
           .or. (restype(1:3) == 'PTR' .and. atmtype(1:3) == 'O2P') &
           .or. (restype(1:3) == 'PTR' .and. atmtype(1:3) == 'O3P') &
           .or. (restype(1:3) == 'SEP' .and. atmtype(1:3) == 'O1P') &
           .or. (restype(1:3) == 'SEP' .and. atmtype(1:3) == 'O2P') &
           .or. (restype(1:3) == 'SEP' .and. atmtype(1:3) == 'O3P') &
           .or. (restype(1:3) == 'TPO' .and. atmtype(1:3) == 'O1P') &
           .or. (restype(1:3) == 'TPO' .and. atmtype(1:3) == 'O2P') &
           .or. (restype(1:3) == 'TPO' .and. atmtype(1:3) == 'O3P')) then
    if (trim(atmtype(1:3)) == 'O' .or. trim(atmtype(1:3)) == 'OXT') then
        is_BB(i_atm) = .true.
    end if

    nA = nA + 1
    acc_atom_index(i_atm) = nA
    hbond_acc(nA)%hybrid_state = 1
    hbond_acc(nA)%n_B_atm = 1
    hbond_acc(nA)%A_atm = i_atm
    hbond_acc(nA)%B_atm(1) = ii_R(ref_res(ref_res_no)%atm_in_bnd(1,atmno), resno)
    hbond_acc(nA)%res   = resno

! 3. Acceptor & Acceptor-connected-atom, SP3 Hybridization state
else if ((restype(1:3) == 'SER' .and. atmtype(1:2) == 'OG') .or. &
         (restype(1:3) == 'SEP' .and. atmtype(1:2) == 'OG') .or. &
         (restype(1:3) == 'THR' .and. atmtype(1:2) == 'OG') .or. &
         (restype(1:3) == 'TPO' .and. atmtype(1:2) == 'OG') .or. &
         (restype(1:3) == 'TYR' .and. atmtype(1:2) == 'OH') .or. &
         (restype(1:3) == 'PTR' .and. atmtype(1:2) == 'OH') .or. &
         (restype(1:3) == 'HOH' .and. trim(atmtype(1:3)) == 'O')) then
    nA = nA + 1
    acc_atom_index(i_atm) = nA
    hbond_acc(nA)%hybrid_state = 2
    hbond_acc(nA)%n_B_atm = 2
    hbond_acc(nA)%A_atm = i_atm
    hbond_acc(nA)%res   = resno
    if (restype(1:3) == 'HOH') then
        call find_atom_idx(resno, ref_res_no,'H1  ',ac_in_res,'NORMAL')
        call find_atom_idx(resno, ref_res_no,'H2  ',ac2_in_res,'NORMAL')
        hbond_acc(nA)%B_atm(1) =  ii_R(ac_in_res, resno) !H1
        hbond_acc(nA)%B_atm(2) =  ii_R(ac2_in_res, resno) !H2
        is_HOH(i_atm) = .true.
    else if (restype(1:3) == 'SER') then 
        call find_atom_idx(resno, ref_res_no,'HG  ',ac2_in_res,'NORMAL')
        hbond_acc(nA)%B_atm(1) =  ii_R(ref_res(ref_res_no)%atm_in_bnd(1,atmno), resno)
        hbond_acc(nA)%B_atm(2) =  ii_R(ac2_in_res, resno) !HG
    else if (restype(1:3) == 'THR') then 
        call find_atom_idx(resno, ref_res_no,'HG1 ',ac2_in_res,'NORMAL')
        hbond_acc(nA)%B_atm(1) =  ii_R(ref_res(ref_res_no)%atm_in_bnd(1,atmno), resno)
        hbond_acc(nA)%B_atm(2) =  ii_R(ac2_in_res, resno) !HG1
    else if (restype(1:3) == 'TYR') then 
        call find_atom_idx(resno, ref_res_no,'HH  ',ac2_in_res,'NORMAL')
        hbond_acc(nA)%B_atm(1) =  ii_R(ref_res(ref_res_no)%atm_in_bnd(1,atmno), resno)
        hbond_acc(nA)%B_atm(2) =  ii_R(ac2_in_res, resno) !HH
    else !TPO,PTR,SEP
        hbond_acc(nA)%n_B_atm = 1
        hbond_acc(nA)%B_atm(1) =  ii_R(ref_res(ref_res_no)%atm_in_bnd(1,atmno), resno)
    endif
     
! 4. Acceptor & Acceptor-connected-atom, Deallocated-aromatic state
else if ((restype(1:3) == 'TRP' .and. atmtype(1:2) == 'NE') .or. &
         (restype(1:3) == 'HIS' .and. atmtype(1:2) == 'NE') .or. &
         (restype(1:3) == 'HIS' .and. atmtype(1:2) == 'ND')) then
    nA = nA + 1
    acc_atom_index(i_atm) = nA
    hbond_acc(nA)%hybrid_state = 3
    hbond_acc(nA)%n_B_atm = 2
    hbond_acc(nA)%A_atm = i_atm
    hbond_acc(nA)%res   = resno
    hbond_acc(nA)%B_atm(1) =  ii_R(ref_res(ref_res_no)%atm_in_bnd(1,atmno), resno)
    if (restype(1:3) == 'HIS' .and. atmtype(1:2) == 'ND') then
        call find_atom_idx(resno, ref_res_no,'CE1 ',ac2_in_res,'NORMAL')
        hbond_acc(nA)%B_atm(2) =  ii_R(ac2_in_res, resno) !CE1
    else if (restype(1:3) == 'HIS' .and. atmtype(1:2) == 'NE') then
        call find_atom_idx(resno, ref_res_no,'CD2 ',ac2_in_res,'NORMAL')
        hbond_acc(nA)%B_atm(2) =  ii_R(ac2_in_res, resno) !CD2
    else ! TRP - NE
        call find_atom_idx(resno, ref_res_no,'CE2 ',ac2_in_res,'NORMAL')
        hbond_acc(nA)%B_atm(2) =  ii_R(ac2_in_res, resno) !CE2
    end if
end if

end subroutine setup_hbond_index
!-------------------------------------------------------------------------------
subroutine calc_hbond_energy(f, g, appl_pair, calc_g)
!-------------------------------------------------------------------------------
! Total Hbond energy calculation of corresponding conformation
! This has same structure with ROSETTA++ function 'get_hbE'
!
! Weights and parameters are different for each H-bond type
!
! H-bond Types are: (Acceptor(hybridization) - Donor)
! 1 : BB(sp2)-BB for helix-helix
! 2 : BB(sp2)-BB for the others         
! 3 : BB(sp2)-BB for beta sheet (may be?)
! 4 : BB(sp2)-BB for turn
! 5 : BB(sp2)-SC general
! 6 : SC(sp2)-BB
! 7 : SC(sp2)-SC
! 8 : SC(sp3)-BB
! 9 : SC(sp3)-SC
!10 : SC(ring)-BB
!11 : SC(ring)-SC
!12 : BB(sp2)-H2O
!13 : SC(sp2)-H2O
!14 : SC(sp3)-H2O
!15 : SC(ring)-H2O
!16 : H2O - BB/SC/H2O general
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f, g(3,tn%atom)
logical, intent(in) :: appl_pair(tn%residue,tn%residue), calc_g
integer :: idA, idD, idH, idAC,idAC2
integer :: i_pair, npair
real(dp) :: f_hb, dfdr, dfdxD, dfdxH
real(dp) :: R_HD(3), R_HA(3), R_AH(3), R_AAC(3), norm_HD, norm_HA, norm_AAC
real(dp) :: drdR(3), dxDdR(3,3), dxHdR(3,3), arg
real(dp) :: buriedness
! Variables for co-operativity
real(dp) :: coop_factor(max_hbpair)
real(dp) :: dffdr(3,2,2,max_hbpair), dffdr1_s(3), dffdr2_s(3)
integer :: n_copair(max_hbpair), copair_id(2,max_hbpair), i_co, co_A, co_H
logical :: exist_pair(tn%residue,tn%residue)

f = 0.0d0
g(:,:) = 0.0d0
! First get hbond pairlist from current conformation
! Cooperativity is also considered if turned on
call get_hbond_pair(appl_pair, npair, n_copair, copair_id, & 
                    coop_factor, exist_pair, dffdr, calc_g)

! Iter over pairs obtained above
do i_pair = 1, npair
    ! Atomic indices for i_pair
    idD   =hbond_pair(i_pair)%don%D_atm 
    idH   =hbond_pair(i_pair)%don%H_atm 
    idA   =hbond_pair(i_pair)%acc%A_atm 
    idAC  =hbond_pair(i_pair)%acc%B_atm(1)
    if (hbond_pair(i_pair)%acc%n_B_atm == 2) then
        idAC2 =hbond_pair(i_pair)%acc%B_atm(2) 
    else
        idAC2 =0
    endif
    ! Get Hbond energy and derivatives
    call hbond_compute_energy(hbond_pair(i_pair)%hbond_type, &
                              hbond_pair(i_pair)%AHdist, &
                              hbond_pair(i_pair)%xD, &
                              hbond_pair(i_pair)%xH,f_hb,dfdr,dfdxD,dfdxH)
    
    if (f_hb>= MAX_HB_ENERGY) cycle
    
    f_hb = f_hb*hb_weight(hbond_pair(i_pair)%hbond_type)
    if (calc_g) then
        dfdr = hb_weight(hbond_pair(i_pair)%hbond_type)*dfdr
        dfdxD = hb_weight(hbond_pair(i_pair)%hbond_type)*dfdxD
        dfdxH = hb_weight(hbond_pair(i_pair)%hbond_type)*dfdxH
    end if
    ! Consider solvation state of a hydrogen bond
    !disabled for water. water doesn't have Ca
    if (use_SAHbond) then        
        if (.not.(is_HOH(idH) .or. is_HOH(idA))) then
            call get_solvation_state(i_R(1,idH), i_R(1,idA), buriedness)

            f_hb = buriedness*f_hb
            if (calc_g) then
                dfdr = buriedness*dfdr
                dfdxD = buriedness*dfdxD
                dfdxH = buriedness*dfdxH
            end if
        end if
    end if

    ! Consider cooperativity
    if (use_coHbond .and. n_copair(i_pair) /= 0) then
        f = f + coop_factor(i_pair)*f_hb
        ! Add gradients
        if (calc_g) then
            do i_co = 1, n_copair(i_pair)
                dffdr1_s(:) = dffdr(:,1,i_co,i_pair)*f_hb
                g(:,idH) = g(:,idH) - dffdr1_s(:)
                g(:,idA) = g(:,idA) + dffdr1_s(:)
                co_H =hbond_pair(copair_id(i_co,i_pair))%don%H_atm 
                co_A =hbond_pair(copair_id(i_co,i_pair))%acc%A_atm 
                dffdr2_s(:) = dffdr(:,2,i_co,i_pair)*f_hb
                g(:,co_H) = g(:,co_H) - dffdr2_s(:)
                g(:,co_A) = g(:,co_A) + dffdr2_s(:)
            end do
        end if
    else
        f = f + f_hb
    end if
    ! Below is only for gradient calculation
    if (.not. calc_g) cycle
    ! Calculate gradient
    R_HD(:) = R(:,idD) - R(:,idH)
    R_HA(:) = R(:,idA) - R(:,idH)
    R_AH(:) = -R_HA(:)
    if (hbond_pair(i_pair)%acc%hybrid_state == 3) then !ring type
        R_AAC(:) = 0.5*(R(:,idAC)+R(:,idAC2)) - R(:,idA)
    else if (hbond_pair(i_pair)%acc%hybrid_state == 2 .and. &
             hbond_pair(i_pair)%acc%n_B_atm == 2) then !sp3,hydrogen exists
        R_AAC(:) = R(:,idAC2) - R(:,idA)
    else    
        R_AAC(:) = R(:,idAC) - R(:,idA)
    endif
    norm_HD = sqrt(dot_product(R_HD(:),R_HD(:)))
    norm_HA = hbond_pair(i_pair)%AHdist
    norm_AAC = sqrt(dot_product(R_AAC(:),R_AAC(:)))

    ! 1. Gradient from distance term
    drdR(:) = -R_HA(:)/hbond_pair(i_pair)%AHdist

    ! 2. Gradient from xD angle
    arg = -1.0d0*hbond_pair(i_pair)%xD
    dxDdR(:,1) = (R_HD(:)*arg/norm_HD - R_HA(:)/norm_HA) / norm_HD
    dxDdR(:,3) = (R_HA(:)*arg/norm_HA - R_HD(:)/norm_HD) / norm_HA
    dxDdR(:,2) = -dxDdR(:,1) - dxDdR(:,3)
     
    ! 3. Gradient from xH angle
    arg = -1.0d0*hbond_pair(i_pair)%xH
    dxHdR(:,1) = (R_AH(:)*arg/norm_HA - R_AAC(:)/norm_AAC) / norm_HA
    dxHdR(:,3) = (R_AAC(:)*arg/norm_AAC - R_AH(:)/norm_HA) / norm_AAC
    dxHdR(:,2) = -dxHdR(:,1) - dxHdR(:,3)

    ! Sum
    g(:,idD) = g(:,idD)                  + coop_factor(i_pair)*dfdxD*dxDdR(:,1)
    g(:,idH) = g(:,idH) + coop_factor(i_pair)*( dfdr*drdR(:) + dfdxD*dxDdR(:,2) + dfdxH*dxHdR(:,1))
    g(:,idA) = g(:,idA) + coop_factor(i_pair)*(-dfdr*drdR(:) + dfdxD*dxDdR(:,3) + dfdxH*dxHdR(:,2))
    
    if (hbond_pair(i_pair)%acc%hybrid_state == 3) then !ring type
        g(:,idAC)  = g(:,idAC)  + 0.5*coop_factor(i_pair)*dfdxH*dxHdR(:,3)
        g(:,idAC2) = g(:,idAC2) + 0.5*coop_factor(i_pair)*dfdxH*dxHdR(:,3)
    
    else if (hbond_pair(i_pair)%acc%hybrid_state == 2 .and. &
             hbond_pair(i_pair)%acc%n_B_atm == 2) then !sp3,hydrogen exists
        g(:,idAC2) = g(:,idAC2) + coop_factor(i_pair)*dfdxH*dxHdR(:,3)
    
    else    
        g(:,idAC)  = g(:,idAC)  + coop_factor(i_pair)*dfdxH*dxHdR(:,3)
    endif
end do

end subroutine calc_hbond_energy
!-------------------------------------------------------------------------------
subroutine calc_hbond_energy_single(f, idH, hb_type) 
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f
integer, intent(in) :: idH, hb_type
integer :: idA, idD, idAC,idAC2
integer :: i_A, i_H
real(dp) :: d, xD, xH, dsqr
real(dp) :: f_hb, dfdr, dfdxD, dfdxH
real(dp) :: Rbnd(3,3)
real(dp) :: buriedness

f = 0.0d0
idD = 0

! Find Hindex
do i_H = 1, nH
    if (hbond_don(i_H)%H_atm == idH) then
        idD = hbond_don(i_H)%D_atm
        exit
    end if
end do

! Bypass if Hindex is errorneous
if (idD == 0) return

! Calculate

do i_A = 1, nA
    idA = hbond_acc(i_A)%A_atm
    idAC = hbond_acc(i_A)%B_atm(1)
    if (hbond_acc(i_A)%n_B_atm == 2) then
        idAC2 =hbond_acc(i_A)%B_atm(2) 
    else
        idAC2 =0
    endif
    
    if (i_R(1,idH) == i_R(1,idA)) cycle
    Rbnd(:,1) = R(:,idD) - R(:,idH)
    Rbnd(:,2) = R(:,idH) - R(:,idA)
    if (hbond_acc(i_A)%hybrid_state == 3) then !ring type
        Rbnd(:,3)= R(:,idA) - 0.5*(R(:,idAC)+R(:,idAC2))
    else if (hbond_acc(i_A)%hybrid_state == 2 .and. &
             hbond_acc(i_A)%n_B_atm == 2) then !sp3,hydrogen exists
        Rbnd(:,3) = R(:,idA) - R(:,idAC2)
    else    
        Rbnd(:,3) = R(:,idA) - R(:,idAC)
    endif

    dsqr = dot_product(Rbnd(:,2),Rbnd(:,2))
    if (dsqr > 9.0d0) cycle

    d = sqrt(dsqr)
    xD = dot_product(Rbnd(:,3),Rbnd(:,2))/sqrt(dot_product(Rbnd(:,3),Rbnd(:,3))*dsqr)
    xH = dot_product(Rbnd(:,2),Rbnd(:,1))/sqrt(dot_product(Rbnd(:,1),Rbnd(:,1))*dsqr)

    ! Get Hbond energy and derivatives
    call hbond_compute_energy(1, d, xD, xH, f_hb, dfdr, dfdxD, dfdxH)
    f_hb = f_hb*hb_weight(hb_type)

    ! Consider solvation state of a hydrogen bond
    if (use_SAHbond) then
        !pair(1) = idH
        !pair(2) = idA
        !call get_solvation_state(pair, i_R(1,idH), i_R(1,idA), buriedness)
        call get_solvation_state(i_R(1,idH), i_R(1,idA), buriedness)
        f_hb = buriedness*f_hb
    end if
    f = f + f_hb
end do

end subroutine calc_hbond_energy_single
!-------------------------------------------------------------------------------
subroutine calc_perres_hbE(f, g) 
!-------------------------------------------------------------------------------
!  Report residue-wise docomposed H-bond energy
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f, g(3,tn%atom)
integer :: hb_type(max_hbpair), idA, idD, idH, idAC, i_pair, npair
real(dp) :: f_hb, dfdr, dfdxD, dfdxH
real(dp) :: perres_energy(2,tn%residue)
integer :: i_res, res1, res2
! Variables for co-operativity

f = 0.0d0
g(:,:) = 0.0d0
perres_energy(:,:) = 0.0d0
! First get hbond pairlist from current conformation
! Cooperativity is also considered if turned on
                    
do i_pair = 1, npair
    idD   =hbond_pair(i_pair)%don%D_atm 
    idH   =hbond_pair(i_pair)%don%H_atm 
    idA   =hbond_pair(i_pair)%acc%A_atm 
    idAC  =hbond_pair(i_pair)%acc%B_atm(1)
    res1 = i_R(1,idA)
    res2 = i_R(1,idD)
    call hbond_compute_energy(hbond_pair(i_pair)%hbond_type, &
                              hbond_pair(i_pair)%AHdist, &
                              hbond_pair(i_pair)%xD, &
                              hbond_pair(i_pair)%xH,f_hb,dfdr,dfdxD,dfdxH)

    f_hb = f_hb*hb_weight(hb_type(i_pair))
    !For pairwise decomposition
    if (hb_type(i_pair) < 5) then
        perres_energy(1,res1) = perres_energy(1,res1) + 0.5d0*f_hb
        perres_energy(1,res2) = perres_energy(1,res2) + 0.5d0*f_hb
    else
        perres_energy(2,res1) = perres_energy(2,res1) + 0.5d0*f_hb
        perres_energy(2,res2) = perres_energy(2,res2) + 0.5d0*f_hb
    end if
    f = f + f_hb
end do

do i_res = 1, tn%residue
    write(log_msg,"(I5,2F10.3)") i_res, perres_energy(:,i_res)
    call log_p(log_msg, me=me, level=30)
end do

end subroutine calc_perres_hbE
!-------------------------------------------------------------------------------
subroutine hbond_compute_energy(hb_type, AHdist, xD, xH, f, dfdr, dfdxD, dfdxH)
!-------------------------------------------------------------------------------
! Single H-bond energy calculation (most fundamental level)
! This has been written following the structure implemented in ROSETTA++ 
! function name 'hbond_compute_energy'
!-------------------------------------------------------------------------------
integer, intent(in) :: hb_type
real(dp), intent(in) :: AHdist, xD, xH
real(dp), intent(out) :: f, dfdr, dfdxD, dfdxH
real(dp) :: f1, f2, f3
real(dp) :: FLr, FSr, FxD, FxH, Pr, PSxD, PSxH, PLxD, PLxH
real(dp) :: dFLr, dFSr, dFxD, dFxH, dPr, dPSxD, dPSxH, dPLxD, dPLxH
PLxD = 0.0d0
PLxH = 0.0d0
dPLxD = 0.0d0
dPLxH = 0.0d0

! Call fade_function: (x, min0, fmin, fmax, max0)
call fade_function(xD, min_xD, angle_interp_edge, 1.0d0, 1.0d0, FxD, dFxD)
call fade_function(xH, min_xH, angle_interp_edge, 1.0d0, 1.0d0, FxH, dFxH)

if (hb_type <= 4) then !BB-BB interaction
    FLr = 0.0d0
    dFLr = 0.0d0
    call fade_function(AHdist, min_r, min_r, r_interp_edge, max_r, FSr, dFSr)
else !Else (Including BB-SC, SC-SC)
    call fade_function(AHdist, interp_min, interp_max, interp_max, max_r, FLr, dFLr)
    call fade_function(AHdist, min_r, min_r, interp_min, interp_max, FSr, dFSr)
end if

select case(hb_type)
case(1) !HelixBB-HelixBB
    call poly8(AHdist, c8_AHdisBBhelix, min_r, 2.8d0, Pr, dPr)
    call poly8(xD, c8_xDBBHelix, min_xD, max_xD, PSxD, dPSxD)
    call poly8(xH, c8_xHBBHelix, min_xH, max_xH, PSxH, dPSxH)
case(2) !BB-BB
    call poly8(AHdist, c8_AHdisBBOther, min_r, 2.745d0, Pr, dPr)
    call poly8(xD, c8_xDBBOther, min_xD, max_xD, PSxD, dPSxD)
    call poly8(xH, c8_xHBBOther, min_xH, max_xH, PSxH, dPSxH)
case(3:7) !Other sp2-H
    call poly5(AHdist, c5_AHdisSP2, min_r, 2.5d0, Pr, dPr)
    call poly3(xD, c3_xDSP2short, min_xD, max_xD, PSxD, dPSxD)
    call poly3(xH, c3_xHSP2short, min_xH, max_xH, PSxH, dPSxH)
    call poly3(xD, c3_xDSP2long, min_xD, max_xD, PLxD, dPLxD)
    call poly3(xH, c3_xHSP2long, min_xH, max_xH, PLxH, dPLxH)
case(8:9) !Other sp3-H
    call poly5(AHdist, c5_AHdisSP3, min_r, 2.5d0, Pr, dPr)
    call poly3(xD, c3_xDSP3short, min_xD, max_xD, PSxD, dPSxD)
    call poly3(xH, c3_xHSP3, min_xH, max_xH, PSxH, dPSxH)
    call poly3(xD, c3_xDSP3long, min_xD, max_xD, PLxD, dPLxD)
    PLxH = PSxH
    dPLxH = dPSxH
     
case(10:11) !Other ring-H
    call poly5(AHdist, c5_AHdisSP2, min_r, 2.5d0, Pr, dPr)
    call poly3(xD, c3_xDSP2short, min_xD, max_xD, PSxD, dPSxD)
    call poly7(xH, c7_xHRing, min_xH, max_xH, PSxH, dPSxH)
    call poly3(xD, c3_xDSP2long, min_xD, max_xD, PLxD, dPLxD)
    PLxH = PSxH
    dPLxH = dPSxH
!added!   
case(12:13) !H2O - BB / SC(sp2)
   call poly5(AHdist, c5_AHdisSP2, min_r, 2.5d0, Pr, dPr)
   call poly3(xD, c3_xDSP2short, min_xD, max_xD, PSxD, dPSxD)
   call poly3(xH, c3_xHSP2short, min_xH, max_xH, PSxH, dPSxH)
   call poly3(xD, c3_xDSP2long, min_xD, max_xD, PLxD, dPLxD)
   call poly3(xH, c3_xHSP2long, min_xH, max_xH, PLxH, dPLxH)
case(14) !H2O - SC(sp3)
   call poly5(AHdist, c5_AHdisSP3, min_r, 2.5d0, Pr, dPr)
   call poly3(xD, c3_xDSP3short, min_xD, max_xD, PSxD, dPSxD)
   call poly3(xH, c3_xHSP3, min_xH, max_xH, PSxH, dPSxH)
   call poly3(xD, c3_xDSP3long, min_xD, max_xD, PLxD, dPLxD)
   PLxH = PSxH
   dPLxH = dPSxH
case(15) !H2O - SC(ring) 
   call poly5(AHdist, c5_AHdisSP2, min_r, 2.5d0, Pr, dPr)
   call poly3(xD, c3_xDSP2short, min_xD, max_xD, PSxD, dPSxD)
   call poly7(xH, c7_xHRing, min_xH, max_xH, PSxH, dPSxH)
   call poly3(xD, c3_xDSP2long, min_xD, max_xD, PLxD, dPLxD)
   PLxH = PSxH
   dPLxH = dPSxH
case(16) !Other - H2O    
   call poly5(AHdist, c5_AHdisSP3, min_r, 2.5d0, Pr, dPr)
   call poly3(xD, c3_xDSP3short, min_xD, max_xD, PSxD, dPSxD)
   call poly3(xH, c3_xHSP3, min_xH, max_xH, PSxH, dPSxH)
   call poly3(xD, c3_xDSP3long, min_xD, max_xD, PLxD, dPLxD)
   PLxH = PSxH
   dPLxH = dPSxH
end select

f1 = Pr*FxD*FxH
f2 = FSr*(PSxD*FxH + FxD*PSxH)
f3 = FLr*(PLxD*FxH + FxD*PLxH)
f = f1 + f2 + f3

dfdr = dPr*FxD*FxH + dFSr*(PSxD*FxH + FxD*PSxH) + dFLr*(PLxD*FxH + FxD*PLxH)
dfdxD = dFxD*(Pr*FxH + FLr*PLxH + FSr*PSxH) + FxH*(FSr*dPSxD + FLr*dPLxD)
dfdxH = dFxH*(Pr*FxD + FSr*PSxD + FLr*PLxD) + FxD*(FSr*dPSxH + FLr*dPLxH)

end subroutine hbond_compute_energy
!-------------------------------------------------------------------------------
subroutine get_hbond_pair(appl_pair,npair, n_copair,copair_id, &
                          coop_factor, exist_pair, dffdr, calc_g)
!-------------------------------------------------------------------------------
! Return pairlist information (atmno/distance/angle/type) 
! from current conformation R
! Cooperativity information is also extracted here if turned on
!-------------------------------------------------------------------------------
logical, intent(in) :: appl_pair(tn%residue,tn%residue), calc_g
integer, intent(out) :: npair !id_hb, hb_type -> hbond_pair
type(hbond_acceptor_type) :: acc
type(hbond_donor_type) :: don
integer :: i_atm,ia,j_atm
real(dp) :: Ri(3), Rac(3), Rij(3), r_sqr, arg
real(dp) :: p(3), q(3), qq
integer, intent(out) :: copair_id(2,max_hbpair), n_copair(max_hbpair)
real(dp), intent(out) :: coop_factor(max_hbpair), dffdr(3,2,2,max_hbpair)
integer :: pair_id(tn%residue,tn%residue), i_pairid(2,max_hbpair)
! Added variables for co-operativity
real(dp) :: d_v(3,max_hbpair), dipole_sinval, d_v1(3), d_v2(3)
real(dp) :: ff, dcf, r1, r2, factor
integer :: i_p1, i_p2, pairing_type
integer :: i_res,j_res,k_res,l_res
logical,intent(out) :: exist_pair(tn%residue,tn%residue)
npair = 0
pair_id(:,:) = 0
i_pairid(:,:) = 0
coop_factor(:) = 1.0d0
dffdr(:,:,:,:) = 0.0d0
n_copair(:) = 0
copair_id(:,:) = 0
d_v(:,:) = 0.0d0
exist_pair(:,:) = .false.
do i_atm = 1, tn%atom
    if (don_atom_index(i_atm) == 0 .and. acc_atom_index(i_atm) == 0) cycle
    
    do ia =1, i_P(i_atm)%n_pair
        j_atm = i_P(i_atm)%i_pair(ia)
        if (don_atom_index(i_atm) >0 .and. acc_atom_index(j_atm)>0 ) then
            don = hbond_don(don_atom_index(i_atm))
            acc = hbond_acc(acc_atom_index(j_atm))
        else if (acc_atom_index(i_atm) >0 .and. don_atom_index(j_atm)>0 ) then
            don = hbond_don(don_atom_index(j_atm))
            acc = hbond_acc(acc_atom_index(i_atm))
        else
            cycle
        endif

        ! if same residue or the pair is not for calculation, cycle
        if (.not. appl_pair(don%res ,acc%res) .or. don%res == acc%res) cycle
        
        !if HA dist >max_r(default: 3.0A), cycle
        Ri = R(:,acc%A_atm) !vector of A_atm
        Rij(:) = R(:,don%H_atm) - Ri(:)
        r_sqr = dot_product(Rij(:), Rij(:))
        if (r_sqr > max_r2 .or. r_sqr < min_r2) cycle
        
        npair = npair + 1
        if (npair > max_hbpair) then
            npair = npair - 1
            call log_p('Warning: Number of Hbond pair exceeds {max_hbpair}.',me=me,level=30)
            call log_p('         Further enumeration will be ignored.',me=me,level=30)
            exit
        end if
    
        hbond_pair(npair)%acc = acc
        hbond_pair(npair)%don = don
        ! consideration for ring type, sp3 type(in this case,B atom changes to H atom)
        if (acc%hybrid_state == 3) then !ring type
            Rac = 0.5d0*(R(:,acc%B_atm(1)) + R(:,acc%B_atm(2)))
        else if (acc%hybrid_state == 2 .and. acc%n_B_atm == 2) then !sp3,hydrogen exists
            Rac = R(:,acc%B_atm(2))
        else    
            Rac = R(:,acc%B_atm(1))
        endif
        !calcuate AHdist, cos(BAH) (=xH) ,cos(AHD) (=xD)
        hbond_pair(npair)%AHdist = sqrt(r_sqr)
        p(:) = Ri(:) - Rac(:)
        q(:) = Rij(:)
        qq = dot_product(q,q)
        arg = dot_product(p,q)/sqrt(dot_product(p,p)*qq)
        hbond_pair(npair)%xH = sign(min(abs(arg),1.0d0),arg) ! make sure abs(arg)<=1

        p(:) = q(:)
        q(:) = R(:,don%D_atm) - R(:,don%H_atm)
        arg = dot_product(p,q)/sqrt(qq*dot_product(q,q))
        hbond_pair(npair)%xD = sign(min(abs(arg),1.0d0),arg) ! make sure abs(arg)<=1
        
        d_v(:,npair) = Rij(:)
        
        ! Assign H-bond type
        ! Backbone-Backbone 
        if (is_BB(acc%A_atm) .and. is_BB(don%D_atm)) then
            ! Save bb-interaction info (for cooperativity consideration at next step)
            exist_pair(acc%res,don%res) = .true.
            pair_id(acc%res,don%res) = npair
            i_pairid(1:2,npair) = (/ acc%res, don%res /)
            
            if( (acc%res - don%res)==4 .or. (acc%res - don%res)==-4) then
                hbond_pair(npair)%hbond_type = 1 !Helix-helix
            else if( (acc%res - don%res)<=3 .and. (acc%res - don%res)>=-3) then
                hbond_pair(npair)%hbond_type = 3 !Strand-strand
            else
                hbond_pair(npair)%hbond_type = 2 !Others
            endif
        ! Backbone-Sidechain
        else if (is_BB(acc%A_atm) .and. .not.is_BB(don%D_atm) .and. &
        (.not. is_HOH(don%D_atm))) then
            hbond_pair(npair)%hbond_type = 5 
        
        ! Sidechain-Backbone
        else if ((.not. is_BB(acc%A_atm) .and. is_BB(don%D_atm)) .and. &
        (.not. is_HOH(acc%A_atm))) then
            if (acc%hybrid_state == 1) then
                hbond_pair(npair)%hbond_type = 6 
            else if (acc%hybrid_state == 2) then
                hbond_pair(npair)%hbond_type = 8 
            else if (acc%hybrid_state == 3) then
                hbond_pair(npair)%hbond_type = 10 
            endif
           
        ! Sidechain-Sidechain
        else if ((.not. is_BB(acc%A_atm) .and. .not. is_BB(don%D_atm)) .and. &
        (.not. is_HOH(acc%A_atm) .and. .not. is_HOH(don%D_atm))) then
            if (acc%hybrid_state == 1) then
                hbond_pair(npair)%hbond_type = 7 
            else if (acc%hybrid_state == 2) then
                hbond_pair(npair)%hbond_type = 9 
            else if (acc%hybrid_state == 3) then
                hbond_pair(npair)%hbond_type = 11 
            endif

        ! H2O - ...
        else if (is_HOH(don%D_atm)) then
            if (is_BB(acc%A_atm)) then
                hbond_pair(npair)%hbond_type = 12 
            else if (is_HOH(acc%A_atm)) then
                hbond_pair(npair)%hbond_type = 16 
            else if (acc%hybrid_state == 1) then
                hbond_pair(npair)%hbond_type = 13 
            else if (acc%hybrid_state == 2) then
                hbond_pair(npair)%hbond_type = 14 
            else if (acc%hybrid_state == 3) then
                hbond_pair(npair)%hbond_type = 15 
            endif
        ! ... - H2O
        else if (is_HOH(acc%A_atm)) then
            hbond_pair(npair)%hbond_type = 16 
        endif

    enddo
    if (npair >= max_hbpair) exit
enddo
! Below are optional procedure for cooperativity consideration
if (.not. use_coHbond) return 
! Co-operativity factor (for interaction of backbone-backbone hbonds)
do i_p1 = 1, npair ! Iter over i-th hbond
    i_res = i_pairid(1,i_p1) !acceptor resno of i-th hbond (bb type only)
    j_res = i_pairid(2,i_p1) !donor resno of i-th hbond (bb type only)
    if (i_res == 0 .or. j_res == 0) cycle ! Consider only if bb hbond exist

    do i_p2 = 1, npair ! Iter over j-th hbond
        k_res = i_pairid(1,i_p2)
        l_res = i_pairid(2,i_p2)
        if (k_res == 0 .or. l_res == 0) cycle ! Consider only if bb hbond exist
        pairing_type = 0

        ! Anti-parallel: Single count only when i < j
        if (i_res == l_res .and. j_res == k_res) then
            factor = coop_factor_antiparallel !0.5
            pairing_type = 1
        ! Parallel: Count both available co-pairs
        else if ((i_res == l_res .and. j_res == k_res+2) .or. &
                 (j_res == k_res .and. i_res == l_res-2)) then
            factor = coop_factor_parallel !0.2
            pairing_type = 2
        end if
        ! Otherwise pass
        if (pairing_type == 0) cycle
        
        d_v1(:) = d_v(:,i_p1) !AHvec of ith hb
        d_v2(:) = d_v(:,i_p2)
        r1 = sqrt(dot_product(d_v1(:), d_v1(:)))
        r2 = sqrt(dot_product(d_v2(:), d_v2(:)))
        dipole_sinval = dot_product(d_v1(:),d_v2(:))/(r1*r2)

        ! Consider only if angle between dipole is proper
        if (dipole_sinval < dipole_max0) then
            call fade_function(dipole_sinval, -2.0d0, -2.0d0, dipole_fmax, dipole_max0, ff, dcf)
            n_copair(i_p1) = n_copair(i_p1) + 1
            copair_id(n_copair(i_p1),i_p1) = i_p2

            ff = ff*factor
            coop_factor(i_p1) = coop_factor(i_p1) + ff
            if (calc_g) then
                dffdr(:,1,n_copair(i_p1),i_p1) = factor*dcf*(dipole_sinval*d_v1(:)/(r1*r1) - d_v2(:)/(r1*r2))
                dffdr(:,2,n_copair(i_p1),i_p1) = factor*dcf*(dipole_sinval*d_v2(:)/(r2*r2) - d_v1(:)/(r1*r2))
            end if
        end if
        if (n_copair(i_p1) == 2) exit ! To prevent weird case
    end do
end do

end subroutine get_hbond_pair
!-------------------------------------------------------------------------------
subroutine get_solvation_state(res1, res2, buriedness) 
!-------------------------------------------------------------------------------
! Calculate solvation state of H-bonding pair residues
!-------------------------------------------------------------------------------
!integer, intent(in) :: idH,idA, res1, res2
integer, intent(in) :: res1, res2
real(dp), intent(out) :: buriedness
integer :: i_res
integer :: n_contact1, n_contact2, state
real(dp) :: r1(3), r2(3), r1_sqr, r2_sqr
real(dp), parameter :: r_contact_sqr = 64.0d0

n_contact1 = 0
n_contact2 = 0

! Count number of C-beta contacts

do i_res = 1, tn%stdres
    if (i_res == res1 .or. i_res == res2) cycle

    r1(:) = R(:,res_index(res1)%Cb_id(1)) - R(:,res_index(i_res)%Cb_id(1))
    r2(:) = R(:,res_index(res2)%Cb_id(1)) - R(:,res_index(i_res)%Cb_id(1))
    r1_sqr = dot_product(r1(:),r1(:))
    r2_sqr = dot_product(r2(:),r2(:))

    if (r1_sqr < r_contact_sqr)  n_contact1 = n_contact1 + 1
    if (r2_sqr < r_contact_sqr)  n_contact2 = n_contact2 + 1
end do

! Determine solvation state of Hbond by combining states from both partner
state = contactno_to_state(n_contact1) + contactno_to_state(n_contact2)
if (state < 2) then !Exposed
    buriedness = buriedness_exposed
else if (state == 2) then !Intermediate
    buriedness = buriedness_intermediate
else !Buried
    buriedness = buriedness_buried
end if

end subroutine get_solvation_state
!-------------------------------------------------------------------------------
function contactno_to_state(n_contact) 
!-------------------------------------------------------------------------------
! Simple solvation state used in Rosetta
!-------------------------------------------------------------------------------
integer, intent(in) :: n_contact
integer :: contactno_to_state
  
if (n_contact < 9) then ! Solvated
    contactno_to_state = 0
else if (n_contact >= 9 .and. n_contact < 15) then ! Medium buried
    contactno_to_state = 1
else ! Fully buried
    contactno_to_state = 2
end if

end function contactno_to_state
!-------------------------------------------------------------------------------
END MODULE HBOND_ENERGY
!-------------------------------------------------------------------------------
