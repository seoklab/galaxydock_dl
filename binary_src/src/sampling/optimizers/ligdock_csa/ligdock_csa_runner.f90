!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/sampling/optimizers/ligdock_csa/ligdock_csa_runner.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE LIGDOCK_CSA_RUNNER
!-------------------------------------------------------------------------------

use globals
use string
!
use logger
use geometry, only: cartesian2internal
use clustering, only: NMRclust
!
use energy_utils, only: protein_to_R
!
use ligand_operator, only: construct_ligand
!
use simplex, only: initialize_simplex, finalize_simplex
!
use ligdock_csa_vars
use ligdock_csa_dist, only: initialize_distance, update_distance, setup_distance, &
                            finalize_distance
use ligdock_csa_log_out
use ligdock_csa_gen_bank, only: generate_first_bank_pool, read_first_bank, &
                                read_add_conf_to_bank
use ligdock_csa_update_bank
use ligdock_csa_new_conf, only: make_new_conf, optimize_new_conf_rigid, &
                                optimize_seed_rigid, optimize_new_conf_flex_sc,&
                                optimize_seed_flex_sc, optimize_bank_steepest
use ligdock_csa_seed, only: select_seeds, done_seed_cycle,&
                            bank_used_reset, bank_all_used

implicit none
save
public

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_ligdock_csa(gene_dim, csa_level)
!-------------------------------------------------------------------------------
integer, intent(in) :: gene_dim
integer, optional, intent(in) :: csa_level
integer :: i_res, i_usc

call log_p('Initializing LigDock CSA...')

if (present(csa_level)) then
    call set_default_ligdock_csa(csa_level)
end if
allocate(idx_seed(n_seed))
call initialize_simplex(gene_dim)

if (n_usc > 0) then
    i_usc = 0
    do i_res = 1, tn%stdres
        if (is_usc(i_res)) then
            i_usc = i_usc + 1
            usc_list(i_usc) = i_res
        end if
    end do
end if

call log_p('Done.')

end subroutine initialize_ligdock_csa
!-------------------------------------------------------------------------------
subroutine finalize_ligdock_csa()
!-------------------------------------------------------------------------------
deallocate(idx_seed)
call finalize_simplex()

end subroutine finalize_ligdock_csa
!-------------------------------------------------------------------------------
subroutine set_default_ligdock_csa(csa_level)
!-------------------------------------------------------------------------------
integer, intent(in) :: csa_level
integer, parameter :: EASY = 1
integer, parameter :: MEDIUM = 2
integer, parameter :: HARD = 3
integer, parameter :: APPS_RIGID = 4
integer, parameter :: APPS_FLEX = 5
integer, parameter :: APPS_BA = 6

n_csa_cycle = 1
n_seed_cycle = 2
max_opt_cycle = 10000

max_trial = 10000
e0max = 1000.0d0
e1max = 1000000.0d0

irr = 1
rep_conf_seed_cut = 0.001d0

min_tor_select = 1
max_tor_select_ratio = 0.4d0
mutation_cut = 0.10d0

factor_init_D_cut = 2.0d0
factor_min_D_cut = 5.0d0
n_opt_to_D_min = 5000

ligdock_prefix = 'GalaxyDock'
print_bank_evol = .false.
print_curr_bank = .false.

if (csa_level == EASY) then
    max_bank = 30
    n_bank_add = 30
    n_opr_1 = 2
    n_opr_2 = 3
    n_opr_3 = 2
    n_opr_4 = 3
    n_opr_5 = 0
    n_opr_6 = 0
    n_opr_7 = 0
    n_opr_8 = 0
    n_seed = 15
    cut_pos_seed_size = 2
    print_bank = .false.
else if (csa_level == MEDIUM) then
    max_bank = 50
    n_bank_add = 50
    n_opr_1 = 2
    n_opr_2 = 3
    n_opr_3 = 2
    n_opr_4 = 3
    n_opr_5 = 0
    n_opr_6 = 0
    n_opr_7 = 0
    n_opr_8 = 0
    n_seed = 25
    cut_pos_seed_size = 2
    print_bank = .false.
else if (csa_level == HARD) then
    max_bank = 100
    n_bank_add = 50
    n_opr_1 = 2
    n_opr_2 = 3
    n_opr_3 = 2
    n_opr_4 = 3
    n_opr_5 = 0
    n_opr_6 = 0
    n_opr_7 = 0
    n_opr_8 = 0
    n_seed = 25
    cut_pos_seed_size = 2
    print_bank = .false.
else if (csa_level == APPS_RIGID) then
    max_bank = 100
    n_bank_add = 100
    n_opr_1 = 2
    n_opr_2 = 3
    n_opr_3 = 2
    n_opr_4 = 3
    n_opr_5 = 0
    n_opr_6 = 0
    n_opr_7 = 0
    n_opr_8 = 0
    n_seed = 25
    cut_pos_seed_size = 2
    print_bank = .true.
else if (csa_level == APPS_FLEX) then
    max_bank = 30
    n_bank_add = 30
    n_opr_1 = 1
    n_opr_2 = 2
    n_opr_3 = 1
    n_opr_4 = 2
    n_opr_5 = 0
    n_opr_6 = 0
    n_opr_7 = 1
    n_opr_8 = 2
    n_seed = 15
    cut_pos_seed_size = 2
    print_bank = .true.
    rep_conf_seed_cut = 0.005d0
else if (csa_level == APPS_BA) then
    max_bank = 30
    n_bank_add = 30
    n_opr_1 = 0
    n_opr_2 = 0
    n_opr_3 = 2
    n_opr_4 = 3
    n_opr_5 = 0
    n_opr_6 = 0
    n_opr_7 = 0
    n_opr_8 = 0
    n_seed = 15
    cut_pos_seed_size = 2
    print_bank = .false.
else 
    write(log_msg,'(A,1X,I1)') 'ERROR: Given CSA level is not pre-defiend.',&
                                csa_level
    call terminate_with_error(log_msg)
end if

n_csa_iter = max_bank/n_bank_add
n_new_conf = n_seed * (n_opr_1 + n_opr_2 + n_opr_3 + n_opr_4 + n_opr_5 + n_opr_6)

n_gen_str = 0
total_nft = 0

lig_conf = 10

end subroutine set_default_ligdock_csa
!-------------------------------------------------------------------------------
subroutine run_ligdock_csa_cycle(protein, ligand, final_E, include_curr_conf, &
                                 i_csa, time_start)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(in) :: ligand
real(dp), intent(out) :: final_E
logical, intent(in) :: include_curr_conf
integer, intent(in) :: i_csa
real(dp), intent(in) :: time_start
!
integer :: i_csa_iter, i_seed_cycle, i_opt_cycle
!
integer :: n_bank           ! current bank size
type(ligdock_bank_type), allocatable :: bank(:)
type(ligdock_bank_type), allocatable :: bank_1st(:)
type(ligdock_bank_type), allocatable :: bank_new(:)
type(ligdock_bank_type), allocatable :: bank_1st_pool(:)
integer :: i_Emax_bank_u, i_Emin_bank_u
!
real(dp) :: Dij(max_bank, max_bank)
!
character(len=len_fname) :: bank_prefix
character(len=3) :: cycle_num
!
real(dp) :: E_bank(max_bank)
integer :: cl_id(max_bank,max_bank), cl_size(max_bank), n_cl, i_cl
integer :: i_bank, i_res, i_usc
type(ligdock_bank_type), allocatable :: rep_bank(:)

allocate(bank(max_bank))
allocate(bank_1st(max_bank))
allocate(bank_new(n_new_conf))
allocate(bank_1st_pool(max_bank))

!call rotate_ligand_random(ref_lig(ligand%lig_type),protein,ligand)
call initialize_distance(ligand)
call generate_first_bank_pool(protein, ligand, include_curr_conf, bank_1st_pool)
!sumin : call ML energy 
call calc_ML_energy_ligand(bank_1st_pool, max_bank, protein, ligand, .false.)
!end sumin
call read_first_bank(bank_1st_pool, bank_1st)

bank = bank_1st
n_bank = n_bank_add
call bank_used_reset(bank, n_bank)

call update_min_max_bank_unit(bank(1:n_bank), n_bank, i_Emax_bank_u, i_Emin_bank_u)
call update_distance(bank(1:n_bank), n_bank, ligand%n_atm, Dij)
call setup_distance(n_bank_add, n_bank, Dij)
call optimize_bank_steepest(protein, ligand, bank, n_bank)

call print_ligdock_csa_log_header()
call log_p(' - First bank info.', level=30)
call print_bank_info(bank_1st, n_bank)

if (print_bank) then
    write(bank_prefix,'(A,A)') trim(ligdock_prefix), "_ib"
    call write_bank(bank_prefix, bank_1st, n_bank, protein, ligand)
end if

CSA_ITER: do i_csa_iter = 1, n_csa_iter
    i_seed_cycle = 0
    SEED_CYCLE: do i_opt_cycle = 1, max_opt_cycle
        call print_seed_cycle_log(i_csa, i_csa_iter,i_opt_cycle)
        call print_ligdock_csa_log(bank, n_bank, i_csa, i_seed_cycle, &
                            i_Emin_bank_u, i_Emax_bank_u, time_start)
        !
        call select_seeds(bank, n_bank, Dij)
        call make_new_conf(protein, ligand, bank, bank_1st, bank_new,&
                           n_bank, n_bank_add, i_seed_cycle)
        if (n_usc > 0) then
            call optimize_new_conf_flex_sc(protein, ligand, bank_new, n_new_conf)
            call optimize_seed_flex_sc(protein, ligand, bank)
        else
            call optimize_new_conf_rigid(protein, ligand, bank_new, n_new_conf, use_grad=.true.)
        end if
        !
        call calc_ML_energy_ligand(bank_new, n_new_conf, protein, ligand, .false.)
        call update_bank(bank, bank_new, n_bank, n_new_conf, ligand%n_atm, &
                         i_Emax_bank_u, i_Emin_bank_u)
        call update_distance(bank(1:n_bank), n_bank, ligand%n_atm, Dij)
        call print_bank_info(bank, n_bank)

        D_cut = D_cut * xctdif
        if (D_cut < D_min) D_cut = D_min

        if (done_seed_cycle(bank, n_bank)) then
            i_seed_cycle = i_seed_cycle + 1
            call bank_used_reset(bank, n_bank)
        end if

        if (i_seed_cycle == n_seed_cycle) exit

        if (print_curr_bank) then
            bank_prefix = 'curr_bank'
            call write_bank(bank_prefix, bank, n_bank, protein, ligand)
        end if
        if (print_bank_evol) then
            call num2str_999(i_opt_cycle, cycle_num)
            bank_prefix = 'bank' // cycle_num
            call write_bank(bank_prefix, bank, n_bank, protein, ligand)
        end if
    end do SEED_CYCLE

    if (n_bank == max_bank) exit

    call bank_all_used(bank, n_bank)
    call read_add_conf_to_bank(bank_1st_pool, bank_1st, bank, n_bank, n_bank_add)
    call update_min_max_bank_unit(bank(1:n_bank), n_bank, i_Emax_bank_u, i_Emin_bank_u)
    call update_distance(bank(1:n_bank), n_bank, ligand%n_atm, Dij)
    call setup_distance(n_bank_add, n_bank, Dij)
end do CSA_ITER

if (ML_rerank) call calc_ML_energy_ligand(bank, n_bank, protein, ligand, ML_rerank)
call sort_csa_bank(bank, n_bank)
call construct_ligand(protein, ligand, bank(1)%gene)
if (n_usc > 0) then
    do i_usc = 1, n_usc
        i_res = usc_list(i_usc)
        protein%residue(i_res)%R = bank(1)%flex_sc(i_usc)%R
    end do
    call cartesian2internal(1, protein%n_res, protein%residue(1:protein%n_res))
    call protein_to_R(protein)
end if
final_E = bank(1)%E_ML

if (print_bank) then
    ! write final bank
    write(bank_prefix,'(A,A)') trim(ligdock_prefix), '_fb'
    call write_bank(bank_prefix, bank, n_bank, protein, ligand)
    ! do clustering
    call update_distance(bank(1:n_bank), n_bank, ligand%n_atm, Dij)
    do i_bank = 1, n_bank
        E_bank(i_bank) = bank(i_bank)%E_ML
    end do
    call NMRclust(Dij, n_bank, E_bank, cl_id, cl_size, n_cl)
    allocate(rep_bank(n_cl))
    do i_cl = 1, n_cl
        rep_bank(i_cl) = bank(cl_id(1,i_cl))
    end do
    write(bank_prefix,'(A,A)') trim(ligdock_prefix), '_cl'
    call write_bank(bank_prefix, rep_bank, n_cl, protein, ligand)
    call write_clust_info(bank_prefix, cl_id, cl_size, n_cl)
    deallocate(rep_bank)
    ! calc colony energy
    write(bank_prefix,'(A,A)') trim(ligdock_prefix), '_co'
    call write_colony_info(bank_prefix, bank, n_bank, ligand)
end if

deallocate(bank)
deallocate(bank_1st)
deallocate(bank_new)
deallocate(bank_1st_pool)
call finalize_distance()

end subroutine run_ligdock_csa_cycle
!-------------------------------------------------------------------------------
END MODULE LIGDOCK_CSA_RUNNER
!-------------------------------------------------------------------------------
