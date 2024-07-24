!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/physics/surface_area.f90
!
! Description: Nonpolar solvation energy calculation by using surface area
!              Contains Haber and Hasel type SA calculations
! Note: SA calculation using FACTS is implemented in FACTS for programming and
!       using efficiency
!
! TODO add comment
! Reference
!
!-------------------------------------------------------------------------------
MODULE SURFACE_AREA

use globals
use logger, only: terminate_with_error
use in_out_vars, only: infile_parameter
use energy_vars, only: R, i_R, i_P, appl_res, SA_type, ASP_type, r_probe, max_neigh, &
    SOLV_HASEL, SOLV_HABER, max_energy
use molecular_mechanics, only: LJ_type

implicit none
save
private

integer :: n_grp
!
integer, allocatable :: neighno_SA(:,:), n_neigh_SA(:) ! Neighbor atom no. index (w/ excluded)
real(dp), allocatable :: rabs_neigh_SA(:,:), i_rabs_neigh_SA(:,:)
!
integer, allocatable :: bnd_type(:,:) ! Bond connection type , used for Hasel method
real(dp), allocatable :: area_ref(:)  ! Reference area to be subtracted (e.g. area at monomer for PPDOCK)
real(dp) :: d_probe          ! radius, diameter of probe
real(dp), allocatable :: p_type(:), SA_vdwr(:) ! connection type & atom-type vdW radius
real(dp), allocatable :: ASP_parm(:)           ! Atomic solvation parameter
real(dp) :: p_con(3)                  ! weighting factors for bond connection type
real(dp), allocatable :: SA_max(:)
real(dp), allocatable :: SA_g_ref(:,:) ! Reference gradient to be subtracted

public :: initialize_SA
public :: finalize_SA
public :: calc_Hasel_SA
public :: refresh_SA_ref
public :: SA_g_ref

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_SA(molecule)
!-------------------------------------------------------------------------------
! TODO comment
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule

allocate(n_neigh_SA(tn%stdatm), neighno_SA(max_neigh,tn%stdatm))
allocate(rabs_neigh_SA(max_neigh,tn%stdatm), i_rabs_neigh_SA(max_neigh,tn%stdatm))
allocate(bnd_type(max_neigh,tn%stdatm))
allocate(area_ref(tn%stdatm))
allocate(SA_g_ref(3,tn%stdatm))

d_probe = 2.0d0*r_probe
area_ref(:) = 0.0d0
n_neigh_SA(:) = 0
neighno_SA(:,:) = 0
bnd_type(:,:) = 0
rabs_neigh_SA(:,:) = 0.0d0
i_rabs_neigh_SA(:,:) = 0.0d0
SA_g_ref(:,:) = 0.0d0

call setup_SA_parm()

end subroutine initialize_SA
!-------------------------------------------------------------------------------
subroutine finalize_SA
!-------------------------------------------------------------------------------
! TODO comment
!-------------------------------------------------------------------------------

deallocate(n_neigh_SA, neighno_SA)
deallocate(rabs_neigh_SA, i_rabs_neigh_SA)
deallocate(bnd_type)
deallocate(area_ref)
!TODO added by curie
deallocate(p_type)
deallocate(SA_vdwr)
deallocate(ASP_parm)
deallocate(SA_max)
deallocate(SA_g_ref)

end subroutine finalize_SA
!-------------------------------------------------------------------------------
subroutine setup_SA_parm()
!-------------------------------------------------------------------------------
! Assign SA parameters for 'Hasel' type SA calculations
! Parameters for FACTS type is not assigned here.
!-------------------------------------------------------------------------------
  integer :: i_grp

if (infile_parameter == trim(data_dir)//'parameters_charmm22.in') then
    n_grp = eng_para%n_atom_cls
    allocate(ASP_parm(n_grp))
else
    n_grp = eng_para%n_atom_cls
    allocate(p_type(n_grp))
    allocate(SA_vdwr(n_grp))
    allocate(ASP_parm(n_grp))
    allocate(SA_max(n_grp))
end if

! If using CHARMM22 polarh parameters
if (infile_parameter == trim(data_dir)//'parameters_charmm22ph.in') then
    ! SA type
    if (SA_type == SOLV_HASEL) then ! Hasel
! Order Atom types defined in CHARMM22ph
!             H         HC          C         CA        CT1
!           CT2        CT3       CPH1       CPH2        CPT
!            CY        CP1        CP2        CP3         CC 
!             N        NR1        NR2        NR3        NH1
!           NH2        NH3        NC2         NY         NP
!             O         OC        OH1          S         OT(water)
!            CP       O2NB        ON3        ON4          P(->S)

        SA_vdwr(:) = &
       (/    1.10d0,    1.10d0,    1.72d0,    1.72d0,    1.80d0, &
             1.90d0,    2.00d0,    1.80d0,    1.80d0,    1.72d0, &
             1.72d0,    1.80d0,    1.90d0,    1.90d0,    1.72d0, &
             1.55d0,    1.55d0,    1.55d0,    1.55d0,    1.55d0, &
             1.60d0,    1.60d0,    1.60d0,    1.55d0,    0.00d0, &
             1.50d0,    1.70d0,    1.52d0,    1.80d0,    0.00d0, &
             1.80d0,    1.52d0,    1.52d0,    1.52d0,    1.80d0  /)

        p_type(:) = &
        (/   1.128d0,   1.128d0,   1.554d0,   1.554d0,   1.276d0, &
             1.045d0,   0.880d0,   1.073d0,   1.073d0,   1.554d0, &
             1.554d0,   1.276d0,   1.276d0,   1.045d0,   1.554d0, &
             1.028d0,   1.028d0,   1.028d0,   1.028d0,   1.028d0, &
             1.215d0,   1.215d0,   1.215d0,   1.028d0,   0.000d0, &
             0.926d0,   0.922d0,   1.080d0,   1.121d0,   0.000d0, &
             1.073d0,   1.080d0,   1.080d0,   1.080d0,   1.121d0  /)

        p_con(:) = (/ 0.8875d0, 0.3516d0, 0.3516d0 /)

    else if (SA_type == SOLV_HABER) then ! Haber
        SA_vdwr(:) = &
        (/ 0.10000d0,0.491489d0, 1.36946d0, 1.57867d0, 1.89590d0, &
            2.28648d0, 2.65143d0, 2.07567d0, 2.07567d0, 1.36946d0, &
            1.36946d0, 1.89590d0, 2.28648d0, 2.28648d0, 1.36946d0, &
          0.323677d0, 1.39086d0, 1.39086d0, 1.39086d0, 0.10000d0, &
            1.46960d0, 1.54421d0, 1.55222d0, 0.10000d0, 0.00000d0, &
            1.68185d0, 1.65213d0, 1.54896d0, 2.18892d0, 0.00000d0, &
            2.07567d0, 1.54896d0, 1.54896d0, 1.54896d0, 2.18892d0  /)
        
        p_type(:) = &
         (/0.988957d0, 1.77757d0, 1.30420d0, 1.27625d0, 1.31643d0, &
            1.18278d0, 1.08302d0, 1.13437d0, 1.13437d0, 1.30420d0, &
            1.30420d0, 1.31643d0, 1.18278d0, 1.18278d0, 1.30420d0, &
            1.13022d0, 1.51070d0, 1.51070d0, 1.51070d0, 1.59251d0, &
            1.26162d0, 1.18710d0, 1.05387d0, 1.59251d0, 0.00000d0, &
            1.03612d0, 1.13999d0, 1.09665d0, 1.68013d0, 0.00000d0, &
            1.13437d0, 1.09665d0, 1.09665d0, 1.09665d0, 1.68013d0  /)

        p_con(:) = (/ 0.3430d0, 0.5075d0, 0.5075d0 /)
    end if

    do i_grp = 1, n_grp
        SA_max(i_grp) = 4.0d0*pi*(SA_vdwr(i_grp)+r_probe)**2
    end do

    ! ASP type
    if (ASP_type == 0) then
        ASP_parm(:) = 1.0d0

    else if (ASP_type == 1) then
        ASP_parm(:) = &
   (/ 0.00000d0, 0.00000d0, 0.01200d0, 0.01200d0, 0.01200d0, &
      0.01200d0, 0.01200d0, 0.01200d0, 0.01200d0, 0.01200d0, &
      0.01200d0, 0.01200d0, 0.01200d0, 0.01200d0, 0.01200d0, &
    -0.06000d0,-0.06000d0,-0.06000d0,-0.06000d0,-0.06000d0, &
    -0.06000d0,-0.06000d0,-0.06000d0,-0.06000d0,-0.06000d0, &
    -0.06000d0,-0.06000d0,-0.06000d0, 0.01200d0, 0.00000d0, &
      0.01200d0,-0.06000d0,-0.06000d0,-0.06000d0, 0.01200d0  /)

    else if (ASP_type == 2) then
        ASP_parm(:) = &
   (/     0.0d0,     0.0d0,   0.025d0,   0.025d0,   0.025d0, &
        0.025d0,   0.025d0,   0.025d0,   0.025d0,   0.025d0, &
        0.025d0,   0.025d0,   0.025d0,   0.025d0,   0.025d0, &
        0.012d0,   0.012d0,   0.012d0,   0.012d0,   0.012d0, &
        0.012d0,   0.012d0,   0.028d0,   0.012d0,   0.012d0, &
        0.012d0,   0.014d0,   0.012d0,   0.018d0,     0.0d0, &
        0.025d0,   0.012d0,   0.014d0,   0.012d0,   0.018d0  /)

    else if (ASP_type == 3) then
        ASP_parm(:) = &
   (/     0.0d0,     0.0d0,   0.018d0,   0.018d0,   0.018d0, &
        0.018d0,   0.018d0,   0.018d0,   0.018d0,   0.018d0, &
        0.018d0,   0.018d0,   0.018d0,   0.018d0,   0.018d0, &
       -0.009d0,  -0.038d0,  -0.038d0,  -0.038d0,  -0.009d0, &
       -0.009d0,  -0.038d0,  -0.038d0,  -0.009d0,  -0.009d0, &
       -0.009d0,  -0.037d0,  -0.009d0,  -0.005d0,     0.0d0, &
        0.018d0,  -0.009d0,  -0.037d0,  -0.009d0,  -0.005d0  /)

    else if (ASP_type == 4) then
        ASP_parm(:) = 0.025d0
        ASP_parm(4) = 0.04d0
    end if

! If using CHARMM19 parameters
else if (infile_parameter == trim(data_dir)//'parameters_charmm19.in') then
! Order of Atom types defined in CHARMM19
!             H         HC         HA         HT          C
!          CH1E       CH2E       CH3E       CR1E         CT
!            CM         CR          N        NC2        NH1
!           NH2        NH3         NP         NR          O
!           OH1         OM         OS         OC         OT
!           OH2         LP         FE          S       SH1E
    if (SA_type == SOLV_HASEL) then ! Hasel
        SA_vdwr(:) = &
    (/   1.128d0,   1.128d0,   0.000d0,   0.000d0,   1.554d0, &
          1.276d0,   1.045d0,   0.880d0,   1.073d0,   0.000d0, &
          0.000d0,   1.554d0,   1.028d0,   1.215d0,   1.028d0, &
          1.215d0,   1.215d0,   0.000d0,   1.028d0,   0.926d0, &
          1.080d0,   0.000d0,   0.000d0,   0.922d0,   0.000d0, &
          0.000d0,   0.000d0,   0.000d0,   1.121d0,   1.121d0 /)

        p_type(:) = &
    (/   1.128d0,   1.128d0,   0.000d0,   0.000d0,   1.554d0, &
          1.276d0,   1.045d0,   0.880d0,   1.073d0,   0.000d0, &
          0.000d0,   1.554d0,   1.028d0,   1.215d0,   1.028d0, &
          1.215d0,   1.215d0,   0.000d0,   1.028d0,   0.926d0, &
          1.080d0,   0.000d0,   0.000d0,   0.922d0,   0.000d0, &
          0.000d0,   0.000d0,   0.000d0,   1.121d0,   1.121d0 /)

        p_con(:) = (/ 0.8875d0, 0.3516d0, 0.3516d0 /)

    else if (SA_type == SOLV_HABER) then ! Haber
        SA_vdwr(:) = &
    (/ 0.10000d0,0.491498d0, 0.00000d0, 0.00000d0, 1.36946d0, &
        1.89590d0, 2.28648d0, 2.65143d0, 2.07567d0, 0.00000d0, &
        0.00000d0, 1.57867d0,0.323677d0, 1.55222d0,0.100001d0, &
        1.46960d0, 1.55421d0, 0.00000d0, 1.39086d0, 1.68185d0, &
        1.54896d0, 0.00000d0, 0.00000d0, 1.65213d0, 0.00000d0, &
        0.00000d0, 0.00000d0, 0.00000d0, 2.18892d0, 1.87829d0 /)

        p_type(:) = &
    (/0.988957d0, 1.77757d0, 0.00000d0, 0.00000d0, 1.30420d0, &
        1.31643d0, 1.18278d0, 1.08302d0, 1.13437d0, 0.00000d0, &
        0.00000d0, 1.27625d0, 1.13022d0, 1.05387d0, 1.59251d0, &
        1.26162d0, 1.18710d0, 0.00000d0, 1.51070d0, 1.03612d0, &
        1.09665d0, 0.00000d0, 0.00000d0, 1.13999d0, 0.00000d0, &
        0.00000d0, 0.00000d0, 0.00000d0, 1.68013d0,0.907302d0 /)

        p_con(:) = (/ 0.3430d0, 0.5075d0, 0.5075d0 /)
    end if

    do i_grp = 1, n_grp
        SA_max(i_grp) = 4.0d0*pi*(SA_vdwr(i_grp)+r_probe)**2
    end do

    ! ASP type
    if (ASP_type == 0) then
        ASP_parm(:) = 1.0d0

    else if (ASP_type == 1) then
        ASP_parm(:) = &
   (/ 0.00000d0,0.000000d0, 0.00000d0, 0.00000d0, 0.01200d0, &
      0.01200d0, 0.01200d0, 0.01200d0, 0.01200d0, 0.00000d0, &
      0.00000d0, 0.01200d0,-0.06000d0,-0.06000d0,-0.06000d0, &
    -0.06000d0,-0.06000d0, 0.00000d0,-0.06000d0,-0.06000d0, &
    -0.06000d0, 0.00000d0, 0.00000d0,-0.06000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.01200d0, 0.01200d0 /)
        
    else if (ASP_type == 2) then
        ASP_parm(:) = &
   (/     0.0d0,     0.0d0,     0.0d0,     0.0d0,   0.025d0, &
        0.025d0,   0.025d0,   0.025d0,   0.025d0,     0.0d0, &
          0.0d0,   0.025d0,   0.012d0,   0.028d0,   0.012d0, &
        0.012d0,   0.012d0,     0.0d0,     0.0d0,   0.012d0, &
        0.012d0,     0.0d0,     0.0d0,   0.014d0,     0.0d0, &
          0.0d0,     0.0d0,     0.0d0,   0.018d0,   0.018d0 /)

    else if (ASP_type == 3) then
        ASP_parm(:) = &
   (/     0.0d0,     0.0d0,     0.0d0,     0.0d0,   0.018d0, &
        0.018d0,   0.018d0,   0.018d0,   0.018d0,   0.018d0, &
        0.018d0,   0.018d0,  -0.009d0,  -0.038d0,  -0.009d0, &
       -0.009d0,  -0.038d0,  -0.009d0,  -0.038d0,  -0.009d0, &
       -0.009d0,  -0.009d0,  -0.009d0,  -0.037d0,  -0.037d0, &
       -0.009d0,     0.0d0,     0.0d0,  -0.005d0,  -0.005d0  /)
    end if

else if (infile_parameter == trim(data_dir)//'parameters_charmm22.in') then
    !use with ASP_type 0 only
    if (ASP_type == 0) then
        ASP_parm(:) = 1.0d0
    end if

else if (infile_parameter == trim(data_dir)//'parameters_martini.in') then
! Order of Atom types defined in Martini
!             N          C      1  P5      2  P4         P3
!            P2      3  P1        Nda         Nd      4  Na
!         5  N0      6 C5M      6 C5C         C4      7  C3
!         8  C2      9  C1        Qda     10  Qd     11 QaD
!        11 QaE         Q0        SP4     12 SP1        SNd
!           SNa        SC5     13 SC4        SC3        SC2

    SA_vdwr(:) = &
    (/   1.500d0,   1.500d0,   2.600d0,   2.600d0,   0.000d0, & ! N/C: For C-alpha correction
          0.000d0,   2.600d0,   0.000d0,   0.000d0,   2.600d0, &
          2.600d0,   2.600d0,   2.600d0,   0.000d0,   2.600d0, &
          2.600d0,   2.600d0,   0.000d0,   2.600d0,   2.600d0, &
          2.600d0,   0.000d0,   0.000d0,   2.400d0,   0.000d0, &
          0.000d0,   0.000d0,   2.400d0,   0.000d0,   0.000d0 /)

    p_type(:) = &
    ! Modification by PHB: 2012/4/4
    ! Naive parameters + Manual correction for aromatic & glutamate
    (/   0.00000d0, 0.00000d0, 0.90000d0, 0.90000d0, 0.00000d0, &
          0.00000d0, 0.90000d0, 0.00000d0, 0.00000d0, 0.90000d0, &
          0.90000d0, 0.90000d0, 0.90000d0, 0.00000d0, 0.90000d0, &
          0.90000d0, 0.90000d0, 0.00000d0, 0.90000d0, 0.90000d0, &
          1.20000d0, 0.00000d0, 0.00000d0, 0.75000d0, 0.00000d0, &
          0.00000d0, 0.00000d0, 0.65000d0, 0.00000d0, 0.00000d0 /)

    p_con(:) = (/ 1.11602d0, 0.00000d0, 0.75840d0 /)

    ! ASP type
    if (ASP_type == 0) then
        ASP_parm(:) = 1.0d0
    else if (ASP_type == 1) then
        ! Derived by parameter fitting to decoy, but seems to be poor...
        ASP_parm(:) = &
    (/ 0.00000d0, 0.00000d0, 0.50735d0, 0.45520d0, 0.00000d0, &
        0.00000d0, 0.35393d0, 0.00000d0, 0.00000d0, 0.33457d0, &
        0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.30858d0, &
        0.38240d0, 0.05966d0, 0.00000d0, 0.03861d0, 0.12917d0, &
        0.12917d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
        0.00000d0, 0.00000d0, 0.18708d0, 0.00000d0, 0.00000d0 /)

        SA_max(:) = &
    (/   0.00000d0,   0.00000d0, 204.01848d0, 262.51540d0,   0.00000d0, &
          0.00000d0, 229.41247d0,   0.00000d0,   0.00000d0, 195.70571d0, & 
        227.81550d0, 229.35880d0, 229.35880d0,   0.00000d0, 228.73088d0, &
        233.25004d0, 230.93507d0,   0.00000d0, 235.44605d0, 263.78993d0, &
        263.78993d0,   0.00000d0,   0.00000d0, 192.06289d0,   0.00000d0, &
          0.00000d0,   0.00000d0, 202.53189d0,   0.00000d0,   0.00000d0 /)

    else if (ASP_type == 3) then
        ! Derived by combining atomic-ASP: 2012/4/4
        ASP_parm(:) = &
    (/ 0.00000d0, 0.00000d0, 0.05000d0, 0.08000d0, 0.00000d0, &
        0.00000d0, 0.05000d0, 0.00000d0, 0.00000d0, 0.05000d0, &
        0.10000d0, 0.12500d0, 0.10000d0, 0.00000d0, 0.10000d0, &
        0.20000d0, 0.25000d0, 0.00000d0,-0.20000d0,-0.10000d0, &
       -0.10000d0, 0.00000d0, 0.00000d0,-0.15000d0, 0.00000d0, &
        0.00000d0, 0.00000d0, 0.20000d0, 0.00000d0, 0.00000d0 /)
 
        SA_max(:) = &
    (/   0.00000d0,   0.00000d0, 200.00000d0, 250.00000d0,   0.00000d0, &
          0.00000d0, 150.00000d0,   0.00000d0,   0.00000d0, 200.00000d0, & 
        200.00000d0, 200.00000d0, 150.00000d0,   0.00000d0, 150.00000d0, &
        150.00000d0, 200.00000d0,   0.00000d0, 150.00000d0, 200.00000d0, &
        250.00000d0,   0.00000d0,   0.00000d0, 125.00000d0,   0.00000d0, &
          0.00000d0,   0.00000d0, 125.00000d0,   0.00000d0,   0.00000d0 /)
    end if

else
    call terminate_with_error('Error: No SA parameters exists for '//trim(infile_parameter))
end if

end subroutine setup_SA_parm
!-------------------------------------------------------------------------------
subroutine calc_Hasel_SA(f, g, area, appl_res, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f, g(3,tn%stdatm),area(tn%stdatm)
logical, intent(in) :: calc_g, appl_res(tn%residue)
real(dp) :: darea(3,0:max_neigh)
integer :: i_atm, j_atm, i_type, i_j
logical :: status

f = 0.0d0
g(:,:) = 0.0d0
call set_neighbor_SA(status)
if (.not. status) then
    f = max_energy
    return
end if

do i_atm = 1, tn%stdatm
    i_type = LJ_type(i_atm)
     
    call calc_Hasel_SA_atm(i_atm, area(i_atm), darea(:,0:max_neigh), appl_res, calc_g)
    darea(:,0:) = area(i_atm)*ASP_parm(i_type)*darea(:,0:)

    f = f + (area(i_atm)-area_ref(i_atm))*ASP_parm(i_type)

    if (calc_g) then
        g(:,i_atm) = g(:,i_atm) + darea(:,0)
        do i_j = 1, n_neigh_SA(i_atm)
            j_atm = neighno_SA(i_j,i_atm)
            g(:,j_atm) = g(:,j_atm) + darea(:,i_j)
        end do
    end if
end do

end subroutine calc_Hasel_SA
!-------------------------------------------------------------------------------
subroutine set_neighbor_SA(status)
!-------------------------------------------------------------------------------
! neighbors required for hasel-type SA calculation
!-------------------------------------------------------------------------------
logical, intent(out) :: status
integer :: i_atm, j_atm, i_j
real(dp) :: d
  
status = .true.
do i_atm = 1, tn%stdatm
    n_neigh_SA(i_atm) = 0
    do i_j = 1, i_P(i_atm)%n_pair
        j_atm = i_P(i_atm)%i_pair(i_j)
        d = i_P(i_atm)%d(i_j)

        n_neigh_SA(i_atm) = n_neigh_SA(i_atm) + 1
        ! Debug mode: errorneous structure
        if (n_neigh_SA(i_atm) > max_neigh) then
            status = .false.
            return
        end if
        neighno_SA(n_neigh_SA(i_atm),i_atm) = j_atm
        rabs_neigh_SA(n_neigh_SA(i_atm),i_atm) = d
        i_rabs_neigh_SA(n_neigh_SA(i_atm),i_atm) = 1.0d0/d
    end do
    bnd_type(1:i_P(i_atm)%pair_end_index(1), i_atm) = 1
    bnd_type(i_P(i_atm)%pair_end_index(1)+1:i_P(i_atm)%pair_end_index(2), i_atm) = 2
    bnd_type(i_P(i_atm)%pair_end_index(2)+1:n_neigh_SA(i_atm), i_atm) = 3
end do
  
do i_atm = 1, tn%stdatm
    do i_j = 1, i_P(i_atm)%n_pair
        j_atm = i_P(i_atm)%i_pair(i_j)
        d = i_P(i_atm)%d(i_j)

        n_neigh_SA(j_atm) = n_neigh_SA(j_atm) + 1
        if (n_neigh_SA(j_atm) > max_neigh) then
            status = .false.
            return
        end if
        neighno_SA(n_neigh_SA(j_atm),j_atm) = i_atm
        rabs_neigh_SA(n_neigh_SA(j_atm),j_atm) = d
        i_rabs_neigh_SA(n_neigh_SA(j_atm),j_atm) = 1.0d0/d
        if (i_j <= i_P(i_atm)%pair_end_index(1)) then
            bnd_type(n_neigh_SA(j_atm), j_atm) = 1
        else if (i_j <= i_P(i_atm)%pair_end_index(2)) then
            bnd_type(n_neigh_SA(j_atm), j_atm) = 2
        else
            bnd_type(n_neigh_SA(j_atm), j_atm) = 3
        end if
    end do
end do
  
end subroutine set_neighbor_SA
!-------------------------------------------------------------------------------
subroutine calc_Hasel_SA_atm(i_atm, area, darea, appl_res, calc_g)
!-------------------------------------------------------------------------------
! TODO comment
!-------------------------------------------------------------------------------
integer, intent(in) :: i_atm
real(dp), intent(out) :: area, darea(3,0:max_neigh)
logical, intent(in) :: calc_g, appl_res(tn%residue)
integer :: j_atm, i_j, i_type
real(dp) :: probability
real(dp) :: p_i, cb, s, b_ij, dbdd
real(dp) :: r_ip, rij_d, r_ijp
real(dp) :: R_ij(3), r_i, r_j, dij, i_dij
real(dp) :: dbdr(3)

area = 0.0d0
darea(:,0:n_neigh_SA(i_atm)) = 0.0d0
if (.not. appl_res(i_R(1,i_atm))) return

i_type = LJ_type(i_atm)
r_i = SA_vdwr(i_type)
r_ip = r_i + r_probe
s = SA_max(i_type) !4.0d0*pi*r_ip*r_ip
p_i = p_type(i_type)

if (s < small_real) return
  
probability = 1.0d0

do i_j = 1, n_neigh_SA(i_atm)
    j_atm = neighno_SA(i_j,i_atm)
    r_j = SA_vdwr(LJ_type(j_atm))
    if (.not. appl_res(i_R(1,j_atm))) cycle

    cb = p_i*p_con(bnd_type(i_j,i_atm))/s
    R_ij(:) = R(:,j_atm) - R(:,i_atm)
    dij = rabs_neigh_SA(i_j,i_atm)
    i_dij = i_rabs_neigh_SA(i_j,i_atm)
     
    r_ijp = r_i + r_j + d_probe

    if (r_ijp < dij) cycle
     
    rij_d = (r_j - r_i)*i_dij
    b_ij = pi*r_ip*(r_ijp - dij)*(1.0d0 + rij_d)
    if (calc_g) then
        dbdd = -pi*r_ip*(1.0d0 + r_ijp*rij_d*i_dij)
        dbdr(:) = dbdd/(1.0d0/cb - b_ij) * R_ij(:)*i_dij
    end if

    probability = min(1.0d0,probability*(1.0d0 - cb*b_ij))

    if (calc_g) then
        darea(:,0) = darea(:,0) + dbdr(:)
        darea(:,i_j) = -dbdr(:)
    end if
end do

area = s*max(0.0d0,probability)

end subroutine calc_Hasel_SA_atm
!-------------------------------------------------------------------------------
subroutine refresh_SA_ref(calc_g)
!-------------------------------------------------------------------------------
  logical, intent(in) :: calc_g
  logical :: is_interres(tn%stdres)
  real(dp) :: area(tn%stdatm), f_tmp, g_tmp(3,tn%stdatm)

  SA_g_ref(:,:) = 0.0d0
  ! Get initial SA of unbound subunits
  is_interres(1:tn%recres) = .true. !Setup only receptor protein
  is_interres(tn%recres+1:tn%stdres) = .false. 
  call calc_Hasel_SA(f_tmp, g_tmp, area, is_interres(:), calc_g)
  area_ref(1:tn%recatm) = area(1:tn%recatm) 
  if (calc_g) SA_g_ref(:,1:tn%recatm) = g_tmp(:,1:tn%recatm)

  is_interres(1:tn%recres) = .false.
  is_interres(tn%recres+1:tn%stdres) = .true. !Setup only ligand protein
  call calc_Hasel_SA(f_tmp, g_tmp, area, is_interres(:), calc_g)
  area_ref(tn%recatm+1:tn%stdatm) = area(tn%recatm+1:tn%stdatm)
  if (calc_g) SA_g_ref(:,tn%recatm+1:tn%stdatm) = g_tmp(:,tn%recatm+1:tn%stdatm)

end subroutine refresh_SA_ref
!-------------------------------------------------------------------------------
END MODULE SURFACE_AREA
