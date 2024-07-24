!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/energy_vars.f90
!
! Description:
!  This file defines parameters and derived types related to energy calculation.
!-------------------------------------------------------------------------------
MODULE ENERGY_VARS
!-------------------------------------------------------------------------------
use globals

implicit none
save
public

! Coordinate indices
real(dp), dimension(:,:), allocatable :: R ! Cartesian coordinate of the system
integer, dimension(:,:), allocatable :: i_R, ii_R ! index & inverse index of R
integer, dimension(:,:), allocatable :: i_L, ii_L ! index & inverse index of R for Ligand

! Distance cutoff related to switch/shift function and pairlist update
real(dp) :: Ron                            ! Distance cut for switch/shift function
real(dp) :: Roff                           ! Distance cut for pairlist
real(dp) :: LRoff                          ! Distance cut for pairlist
real(dp) :: Ron_sqr, Roff_sqr              ! Dummy variables for fast calculation
real(dp) :: LRoff_sqr                      ! Dummy variables for fast calculation
real(dp) :: R3on_off, Ron_off_3            ! Dummy variables for fast calculation

! Index of atoms in order of degree of freedom
integer, allocatable :: dof_id(:)          

! Energy evaluation residues/pairs
logical, allocatable :: appl_res(:), appl_respair(:,:), appl_symm_rsr(:)

!-------------------------------------------------------------------------------
type atm_parm_type
!-------------------------------------------------------------------------------
! TODO: Comments
!-------------------------------------------------------------------------------
logical :: fixed              ! Whether fixed by constraints
logical :: excluded           ! Whether excluded by cut-off
real(dp) :: mass              ! Atomic mass
real(dp) :: vdwradii          ! van der Waals radii
logical :: is_H               ! Whether is hydrogen

end type atm_parm_type
!-------------------------------------------------------------------------------
type(atm_parm_type), allocatable :: atm_parm(:)
!
real(dp), parameter :: max_energy = 9.99999d6
integer, parameter :: max_bonded_pair = 50 ! max no of bonded atoms; 1-2,3,4 pairs
!integer, parameter :: max_neigh = 1000     ! max no of neighboring atoms at 1st shell 
integer, parameter :: max_neigh = 2000     ! max no of neighboring atoms at 1st shell 
integer, parameter :: max_neigh2 = max_neigh*4 ! max no of neighboring atoms at 2nd shell
!-------------------------------------------------------------------------------
type pair_index_type
!-------------------------------------------------------------------------------
! TODO: Comments
!-------------------------------------------------------------------------------
! Bonded pairs
integer :: n_bnd(3)                 ! No. bonded pairs:
                                    ! 1-> 1-2; 2-> 1-3; 3-> 1-4
integer :: i_bnd(max_bonded_pair,3) ! Index bonded pairs
                                    ! 1-> 1-2; 2-> 1-3; 3-> 1-4
integer :: pair_end_index(3)        ! Indices for each interaction types: 
                                    ! 1-> 1-2; 2-> 1-3; 3-> 1-4
!
! Non-bonded pairs
integer :: n_pair, n_Lpair          ! Number of neighbor pairs within Roff/LRoff
integer :: i_pair(max_neigh)        ! Index of neighboring atom with Roff
integer :: i_Lpair(max_neigh2)      ! Index of neighboring atom with LRoff
!
real(dp) :: d(max_neigh)   ! Pre-calculated distances with Roff

end type pair_index_type
!-------------------------------------------------------------------------------
type(pair_index_type), allocatable :: i_P(:)

!-------------------------------------------------------------------------------
integer, parameter :: n_E_component_modeling = 23
integer, parameter :: n_E_component_ppdock   = 4
integer, parameter :: n_E_component_ligdock  = 15

character(len=len_fname), parameter :: energy_header_modeling = ''
character(len=len_fname), parameter :: energy_header_ppdock   = ''
character(len=len_fname), parameter :: energy_header_ligdock  = ''
!-------------------------------------------------------------------------------
type energy_type
!-------------------------------------------------------------------------------
! TODO: Comment
!-------------------------------------------------------------------------------
real(dp) :: total
real(dp) :: modeling(0:n_E_component_modeling)
real(dp) :: ppdock  (0:n_E_component_ppdock)
real(dp) :: ligdock (0:n_E_component_ligdock)
!-------------------------------------------------------------------------------
end type energy_type
!-------------------------------------------------------------------------------

integer :: energy_print_level
logical :: report_pairwise
logical :: energy_on_full ! report full-atom energy
logical :: reinitialize

! Flags for energy calculation
logical :: use_modeling_E, use_ppdock_E, use_ligdock_E

! Energy scale
type(energy_type) :: Escale

! Molecular Mechanics 
logical :: use_bond, use_LJ, use_softsphere, use_elec
logical :: use_vdw_scwrl4 ! (for scwrl4 vdw)
integer :: rattle_type ! 0: not applied 1: with hydrogen 2: all
!
real(dp) :: bond_w, vdw_w, elec_w
real(dp) :: LJ_att_w, LJ_rep_w ! weight for attractive or repulsive VDW
!
logical :: rdie ! Whether use distance dependent dielectric
real(dp) :: die_const_pro, die_const_slv ! Dielectric constant
!
logical :: soften_short
character(len=len_fname) :: vdw_type
real(dp) :: vdw_soft_k ! Soft-sphere force constant
real(dp) :: vdw_r_scale, trclash_scale, scclash_scale ! radius_scale
logical :: skip_fixed
integer :: pair_type_cut
!
logical :: truncate_by_shift ! Use shift_function

! Solvation
logical :: use_solv, use_SA, use_factsmem, use_SApp
!
real(dp) :: solv_w, SA_w
!
real(dp) :: r_probe ! Diameter of probe for solvation
integer :: solv_type, SA_type, ASP_type
integer, parameter :: SOLV_FACTS = 1
integer, parameter :: SOLV_EEF1 = 2
integer, parameter :: SOLV_HABER = 4
integer, parameter :: SOLV_HASEL = 5
logical :: update_reff

! Parameters for dDFIRE/DFIRE
logical :: use_ddfire
!
real(dp) :: dfire_w
real(dp) :: ddfire_add_scale
!
character(len=len_fname) :: dfire_score_file, dfire_atype_file
character(len=len_fname) :: kgb_water_file
logical :: use_water_kgb, use_water_kgb_only

! Parameters for KGB
logical :: use_kgb
real(dp) :: kgb_w, wkgb_w

! GOAP
logical :: use_goap
real(dp) :: goap_w
logical :: use_goap_discrete

! Restraints
logical :: use_rsr
logical :: use_meld
integer, parameter :: max_n_features = 11
integer, parameter :: max_n_grp = 5
real(dp) :: meld_active(max_n_grp, max_n_features)
!
real(dp) :: rsr_w
!
real(dp) :: update_rsr_dcut !Dcut for rsr update
real(dp) :: update_rsr_sig, update_lig_rsr_sig ! sigma for rsr update E
real(dp) :: segrsr_scale 
logical :: use_update_rsr, update_NO_rsr, update_lig_rsr
character(len=len_fname) :: rsr_file, rsr_pdb, infile_update_rsr_pdblist

! Distogram
logical :: use_distogram
real(dp) :: distogram_w
character(len=len_fname) :: distogram_file, distogram_pdb

! Torsion
logical :: use_bb_torsion
real(dp) :: bb_torsion_w
character(len=len_fname) :: bb_torsion_file

! Machine Learning score
logical :: use_ml
real(dp) :: ml_w

! Hydrogen bond
logical :: use_hbond
!
real(dp) :: hbond_w
!
logical :: use_coHbond, use_SAHbond

! Ramachandran score
logical :: use_rama_score
!
real(dp) :: rama_w

! Rotamer score
logical :: use_rotamer_score
!
real(dp) ::  rot_w 

! Conservation Score
logical :: use_conserve
!
real(dp) :: conserve_w
!
logical :: blastfile_defined 
character(len=len_fname) :: blastfile

! InterEV score
logical :: use_interev
real(dp) :: interev_w

! symmetric restraint energy
logical :: use_symm_rsr
real(dp) :: symm_rsr_w

! Autodock score
logical :: use_atdk3, use_atdk4
!
! TODO: In atdk_*_w, FE_*_coeff should be included.
real(dp) :: atdk_w
real(dp) :: atdk_vdw_w, atdk_hbond_w, atdk_elec_w, atdk_solv_w
real(dp) :: atdk_int_v_w, atdk_int_h_w, atdk_int_e_w, atdk_int_s_w, atdk_int_w
!
character(len=len_fname) :: infile_atdk_prm, infile_atdk3_sol_prm

! DrugScore
logical :: use_drugscore
!
real(dp) :: drugscore_w
!
character(len=len_fname) :: infile_drugscore_prm

! X-score
logical :: use_Xscore
real(dp) :: X_HM_w, conf_w

! PLP_tor score
logical :: use_PLP
real(dp) :: PLP_w

! dSAS
logical :: use_dSAS
real(dp) :: dSAS_w

! ROTA score
logical :: use_ROTA
real(dp) :: ROTA_w, prot_w

! For ligand docking only. docking_grid related variables
!-------------------------------------------------------------------------------
type docking_grid_param_type
!-------------------------------------------------------------------------------
integer :: n_elem(3)
real(dp) :: grid_cntr(3), grid_width
!-------------------------------------------------------------------------------
end type docking_grid_param_type
!-------------------------------------------------------------------------------
type docking_grid_type
!-------------------------------------------------------------------------------
real(dp), allocatable :: atdk_grid(:,:,:)
real(dp), allocatable :: drugscore_grid(:,:,:)
real(dp), allocatable :: HM_grid(:,:)
real(dp), allocatable :: rsr_grid(:,:,:)
!-------------------------------------------------------------------------------
end type docking_grid_type
!-------------------------------------------------------------------------------
type docking_FFT_grid_type
!-------------------------------------------------------------------------------
integer :: n_atdk, n_drugscore, n_rsr
Complex, allocatable :: atdk_grid(:,:)
Complex, allocatable :: drugscore_grid(:,:)
Complex, allocatable :: HM_grid(:)
Complex, allocatable :: rsr_grid(:,:)
!-------------------------------------------------------------------------------
end type docking_FFT_grid_type
!-------------------------------------------------------------------------------

type(docking_grid_param_type) :: dock_grid_info
type(docking_grid_type) :: dock_grid
type(docking_FFT_grid_type) :: r_grid, l_grid
logical :: use_input_cntr
logical :: soften_dock_E

!-------------------------------------------------------------------------------
END MODULE ENERGY_VARS
!-------------------------------------------------------------------------------
