!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/sampling/optimizers/ligdock_csa/ligdock_csa_vars.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE LIGDOCK_CSA_VARS
!-------------------------------------------------------------------------------

use globals
use energy_vars

implicit none
save
public

integer, parameter :: max_gene_dim = 6 + (max_lig_br-1)
integer, parameter :: max_usc_for_dock = 20

!-------------------------------------------------------------------------------
type flex_sc_type
!-------------------------------------------------------------------------------
integer :: n_atm
real(dp), allocatable :: R(:,:)
!-------------------------------------------------------------------------------
end type flex_sc_type
!-------------------------------------------------------------------------------
type ligdock_bank_type
!-------------------------------------------------------------------------------
! ligand coordinates
!real(dp) :: R(3,max_lig_atom)
real(dp), allocatable :: R(:,:)
! flexible sc
integer, allocatable :: rot_idx(:)
type(flex_sc_type), allocatable :: flex_sc(:)
! gene information
real(dp), dimension(max_gene_dim) :: gene
integer :: gene_dim
! energy
real(dp) :: E(0:n_E_component_ligdock)
real(dp) :: E_ML
integer :: nft
! flags for seed selection
logical :: used
!-------------------------------------------------------------------------------
end type ligdock_bank_type
!-------------------------------------------------------------------------------

! variables for csa cycles
integer :: n_csa_cycle      ! No. of independent CSA runs
integer :: n_csa_iter       ! No. of iteration for one csa runs
                            ! n_csa_iter = int(max_bank/n_bank_add)
integer :: n_seed_cycle     ! No. of seed cycle for one csa iteration
integer :: max_opt_cycle    ! Maximum number of optimization to achieve
                            ! n_seed_cycle

! variables related to bank
integer :: max_bank         ! Maximum bank size
integer :: n_bank_add       ! # of conformation to add to bank in each CSA
                            ! iteration
integer :: max_trial        ! maximum number of trials to generate first bank
real(dp) :: e0max, e1max    ! energy cutoff used in first bank generation

! variables related to seed selection
integer :: n_seed             ! # of seeds to generate new conformations
integer :: cut_pos_seed_size  ! cutoff to terminate SEED_CYCLE
integer :: irr                ! TODO: comment
real(dp) :: rep_conf_seed_cut ! TODO: comment

! variables related to generate new conformation
integer :: n_opr_1, n_opr_2, n_opr_3, n_opr_4, n_opr_5, n_opr_6, n_opr_7,n_opr_8
integer :: min_tor_select
real(dp) :: max_tor_select_ratio
real(dp) :: mutation_cut

! variables related to distance cutoff in CSA
real(dp) :: D_ave, D_cut, D_min
real(dp) :: factor_init_D_cut, factor_min_D_cut, xctdif
integer :: n_opt_to_D_min

! variables related to write bank
character(len=len_fname) :: ligdock_prefix
logical :: print_bank_evol, print_curr_bank, print_bank

! For CSA log (to estimate amount of calculation)
integer :: n_gen_str
integer :: total_nft

integer :: n_unused           ! # of conformations in bank which are not used as
                              ! seed.
integer, allocatable :: idx_seed(:)

integer :: n_new_conf

! variables related to initial docking
integer :: lig_conf                      ! lig_conf = 10 -> random 
                                         !          = 20 -> BetaDock (Voro)
integer :: n_gen_conf, n_voro_per_conf   ! max_trial = n_gen_conf * n_voro_per_conf
character(len=len_fname) :: voro_prefix

! variable related to USC information
integer :: usc_list(max_usc_for_dock)

! variable rerlated to ML_rerank
logical :: ML_rerank
!-------------------------------------------------------------------------------
END MODULE LIGDOCK_CSA_VARS
!-------------------------------------------------------------------------------
