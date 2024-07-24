!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ligdock/physics/autodock_vars.f90
!
! Description:
!   variables/derived types needed to autodock energy calculation
!-------------------------------------------------------------------------------
MODULE AUTODOCK_VARS
!-------------------------------------------------------------------------------

use globals

implicit none
save
public

!===============================================================================
! PARAMETERS
!===============================================================================
integer, parameter :: max_atdk_atm_type = 35  ! maximum # of autodock atomtypes
! related to energy table (T_vdw, T_diel, T_solv, int_E_table
integer, parameter :: max_dist_long = 5000    ! for coulomb, internal energy
integer, parameter :: max_dist_short = 1000   ! for vdw/hbond/desolv
real(dp), parameter :: divisor = 100.0d0      ! for distance binning.
real(dp), parameter :: max_E = 100000.0d0     ! upper cutoff of energy
real(dp), parameter :: smooth = 0.500d0       ! range of smoothing
! parameter used in desolvation
real(dp), parameter :: sol_sig = 3.6d0        
real(dp), parameter :: sol_const = -1.0d0 / (2.0d0 * (sol_sig**2))
real(dp), parameter :: solpar_q = 0.01097d0
! energy scale for 1-4 interaction pair
real(dp), parameter :: scale_1_4 = 0.50d0
! Cutoff of energy calculation(vdw/hbond/desolvation)
real(dp), parameter :: NBcutoff = 8.0d0
real(dp), parameter :: NBcutoff2 = NBcutoff**2

!===============================================================================
! DERIVED TYPES
!===============================================================================
type ad4_param_type
!-------------------------------------------------------------------------------
! Collection of autodock4 parameters
!-------------------------------------------------------------------------------
character(len=2), dimension(max_atdk_atm_type) :: atm_type
real(dp), dimension(max_atdk_atm_type) :: Rii, epsii ! Non H-bond
real(dp), dimension(max_atdk_atm_type) :: Rij_hb, epsij_hb ! H-bond
integer, dimension(max_atdk_atm_type) :: hbond, bond_index
real(dp), dimension(max_atdk_atm_type) :: vol, solpar ! solvation
!-------------------------------------------------------------------------------
end type ad4_param_type
!-------------------------------------------------------------------------------
type atdk_eng_para_type 
!-------------------------------------------------------------------------------
! Collection of autodock4 parameters which are used in current job.
!-------------------------------------------------------------------------------
integer :: n_atm_types(2) ! Number of atom types 
                          ! 1 -> total ; 2 -> lig
integer :: lig_atom_types(max_atdk_atm_type) ! list of ligand_atm_type idx
integer, dimension(max_atdk_atm_type) :: type_idx          ! Atom type index
character(len=2), dimension(max_atdk_atm_type) :: type_name ! Atom type name
logical, dimension(max_atdk_atm_type) :: is_Hdon, is_Hacc_N, is_Hacc_O
! parameters related to autodock energy.
! below arrays have dimension of (ligand_type, prot_type)
real(dp), dimension(max_atdk_atm_type,max_atdk_atm_type) :: ljrij, ljeij ! L-J
logical, dimension(max_atdk_atm_type,max_atdk_atm_type):: hbond_read ! use hbond parameter or not
real(dp), dimension(max_atdk_atm_type,max_atdk_atm_type):: hbrij, hbeij  ! Hbond 
real(dp), dimension(max_atdk_atm_type) :: solpar, vol ! solvation
                                                        ! solpar = sol_S,
                                                        ! vol = sol_V
real(dp), dimension(max_atdk_atm_type) :: solcn ! for autodock3 only
!-------------------------------------------------------------------------------
end type atdk_eng_para_type
!-------------------------------------------------------------------------------
type atdk_hbond_type    
!-------------------------------------------------------------------------------
! Parameters related to calculate autodock H-bond term
!-------------------------------------------------------------------------------
integer :: rexp              ! power of cosine term
real(dp), dimension(3,2) :: rvec   ! bond direction : 1 is for rvec1, 2 is for rvec2
real(dp) :: racc, rdon             !  constant for hydrogen bond 
!-------------------------------------------------------------------------------
end type atdk_hbond_type
!-------------------------------------------------------------------------------
type(atdk_eng_para_type) :: atdk_para
type(atdk_hbond_type), allocatable :: atdk_hbond(:)

! Energy table
real(dp), allocatable, dimension(:,:,:) :: T_vdw  
real(dp), allocatable, dimension(:) :: T_diel, T_solv
real(dp), allocatable, dimension(:,:,:,:) :: int_E_table

! To store atom types
integer, allocatable :: atdk_para_idx(:)

! To flag target ligand
logical, allocatable :: is_ligand(:)
integer :: target_lig_no        ! ligand%lig_no + tn%nonlig

! for atomic charges
real(dp), allocatable :: q_s(:)

! nonbonded_matrix to calculate internal energy
integer, allocatable :: nb_matrix(:,:)
real(dp), allocatable :: qiqj(:), sum_SiVj(:)

!-------------------------------------------------------------------------------
END MODULE AUTODOCK_VARS
