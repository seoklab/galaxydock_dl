!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/sampling/optimizers/ligdock_csa/ligdock_csa_gen_bank.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE LIGDOCK_CSA_GEN_BANK
!-------------------------------------------------------------------------------

use globals
use logger
use ran
use string, only: num2str
use mathfunctions, only: v_norm, bound_ang
use allocate_molecule, only: deallocate_molecule_type
use geometry
use rotamer, only: identify_rotamer_state_single, perturb_rotamer_state, &
                   get_rotamer_index, place_rotamer
use rmsd, only: prep_RMSD_PM_calculation, calc_RMSD_PM
!
use energy_vars
use energy_utils,   only: update_R_for_ligand, update_R_for_sc, protein_to_R
use ligdock_energy, only: check_clash_within_ligand, ligdock_energy_using_grid,&
                          set_atdk_hbond_flex, update_prot_prot_intrxn
use xlogp_Score_m, only: get_xlogP_atom_types
!
use ligand_operator, only: construct_ligand, get_initial_ligand_gene
use simplex
use gradient_minimization, only: opt_ligand_steepest_descent_intcrd
!
use ligdock_csa_vars
use ligdock_csa_utils
use transrot_operator, only: translate_ligand, rotate_ligand
use ligdock_csa_fragFFT, only: hotspot_graph_type, rec_graph, generate_lig_graph, &
                               do_pharmdock, generate_fragment_binding_map
!use ligdock_csa_pharmFFT

implicit none
save
private


public :: generate_first_bank_pool
public :: read_first_bank
public :: read_add_conf_to_bank

CONTAINS
!-------------------------------------------------------------------------------
subroutine generate_first_bank_pool(protein, ligand, include_curr_conf, &
                                    bank_1st_pool)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(in) :: ligand
logical, intent(in) :: include_curr_conf
type(ligdock_bank_type), intent(out) :: bank_1st_pool(:)
integer :: i_iter, bank_i, bank_f
real(dp) :: time1

call log_p('Generating first bank...')
call my_timer(time1)
if (lig_conf == 30) then ! using FFT
    call generate_fragment_binding_map(protein, ligand)
end if
call finish_timer(time1)
do i_iter = 1, n_csa_iter
    bank_i = n_bank_add*(i_iter-1) + 1
    bank_f = n_bank_add*i_iter
    if (lig_conf == 10) then
        call log_p('First bank is generated randomly', level=20)
        if (n_usc > 0) then
            call generate_bank_flex_sc(protein, ligand, n_bank_add, &
                                       bank_1st_pool(bank_i:bank_f),&
                                       include_curr_conf)
        else
            call generate_bank_rigid(protein, ligand, n_bank_add, &
                                     bank_1st_pool(bank_i:bank_f), &
                                     include_curr_conf)
        end if
    else if (lig_conf == 20) then
        call log_p('First bank is generated using pre-dock', level=20)
        if (n_usc > 0) then
            call generate_bank_VD_flex_sc(protein, ligand, n_bank_add, &
                                          bank_1st_pool(bank_i:bank_f))
        else
            call generate_bank_VD_rigid(protein, ligand, n_bank_add, &
                                        bank_1st_pool(bank_i:bank_f))
        end if
    else if (lig_conf == 30) then ! using FFT
        call log_p('First bank is generated using FFT', level=20)
        call generate_bank_FFT(protein, ligand, n_bank_add, &
                               bank_1st_pool(bank_i:bank_f), include_curr_conf)
    end if
end do
call log_p('Done.')

end subroutine generate_first_bank_pool
!-------------------------------------------------------------------------------
subroutine read_first_bank(bank_1st_pool, bank_1st)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: bank_1st_pool(:)
type(ligdock_bank_type), intent(out) :: bank_1st(:)
integer :: i
   
do i = 1, n_bank_add
    bank_1st(i) = bank_1st_pool(i)
    n_gen_str = n_gen_str + 1
    total_nft = total_nft + bank_1st(i)%nft
end do

end subroutine read_first_bank
!-------------------------------------------------------------------------------
subroutine read_add_conf_to_bank(bank_1st_pool, bank_1st, bank, n_bank, &
                                 n_bank_add)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: bank_1st_pool(:)
type(ligdock_bank_type), intent(inout) :: bank_1st(:), bank(:)
integer, intent(inout) :: n_bank
integer, intent(in) :: n_bank_add
integer :: i

do i = n_bank + 1, n_bank + n_bank_add
    bank_1st(i) = bank_1st_pool(i)
    bank(i) = bank_1st(i)
    n_gen_str = n_gen_str + 1
end do

n_bank = n_bank + n_bank_add

end subroutine read_add_conf_to_bank
!-------------------------------------------------------------------------------
subroutine generate_bank_FFT(protein, ligand, n_bank, bank, include_curr_conf)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type), intent(out) :: bank(:)
integer, intent(in) :: n_bank
logical, intent(in) :: include_curr_conf
!
integer, parameter :: max_gene = 500
!
type(molecule_type) :: tmp_prot
type(ligdock_bank_type) :: tmp_bank, soln_pool(max_gene)
integer :: i_trial, gene_dim, i_gene, n_gene
!
integer :: n_succ, n_pool, n_soln, n_pool_e0, n_soln_e0, i_soln
type(ligdock_bank_type), allocatable :: bank_pool(:), conf_bank(:)
type(energy_type) :: ff
real(dp) :: g(3,ligand%n_atm), Ecut
integer :: status, ierr
!
integer :: i_frag
integer :: i, i_bank
!
type(hotspot_graph_type) :: lig_graph
real(dp), allocatable :: gene_s(:,:)
!
integer, parameter :: cand_per_trial = 50
integer :: tnf
integer :: max_candidate_pool
! DEBUG
integer :: i_rot, n_rot
real(dp) :: rmsd, gene(5+ligand%n_br), diff(ligand%n_br-1), prev_E, min_prev_E
real(dp) :: trans, dr(3), rot
integer :: n_type_d, n_atm_s_d(ligand%n_atm), atm_id_s_d(ligand%n_atm, ligand%n_atm)
character(len=20) :: atom_types_d(ligand%n_atm)

tmp_prot = protein
gene_dim = 6 + (ligand%n_br-1)
!
n_succ = 0
n_pool = 0
n_pool_e0 = 0
max_candidate_pool = min(n_bank*50, max_trial)
allocate(conf_bank(n_bank))
allocate(bank_pool(max_trial))
allocate(gene_s(gene_dim, max_gene))
!
call allocate_ligdock_bank_type(tmp_bank, protein, ligand)
if (include_curr_conf) then
    write(log_msg, '(A, I5)') 'Include input conformation in first bank.'
    call log_p(log_msg, level=30)
    call get_initial_ligand_gene(tmp_prot%ligand(ligand%lig_no), ligand, &
                                 tmp_bank%gene)
    tmp_bank%nft = 1
    call do_simplex(tmp_prot, ligand, tmp_bank%gene, tmp_bank%E, tmp_bank%nft)
    call ligand2bank(tmp_prot%ligand(ligand%lig_no), tmp_bank)

    n_succ = n_succ + 1
    bank(n_succ) = tmp_bank
end if

do i_trial = 1, max_trial
    tmp_prot = protein ! TODO: remove this
    call generate_conf_lib(tmp_prot, ligand, conf_bank, n_bank, gene_dim)
    
    do i_bank = 1, n_bank
        tmp_bank = conf_bank(i_bank)
        call construct_ligand(tmp_prot, ligand, tmp_bank%gene(1:gene_dim))
        call generate_lig_graph(tmp_prot%ligand(ligand%lig_no), lig_graph)
        do i_gene = 1, max_gene
            gene_s(1:gene_dim, i_gene) = tmp_bank%gene(1:gene_dim)
        end do
        call do_pharmdock(tmp_prot, ligand, rec_graph, lig_graph, &
                          tmp_bank%gene(1:3), gene_s, n_gene, e1max, max_gene)
        if (n_gene < 1) cycle
        !
        n_soln = 0
        do i_gene = 1, n_gene
            tmp_bank%gene(1:gene_dim) = gene_s(1:gene_dim, i_gene)
            call construct_ligand(tmp_prot, ligand, tmp_bank%gene(1:gene_dim))
            call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
            call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
            if (ff%ligdock(0) > e1max) cycle
            !
            call do_simplex(tmp_prot, ligand, tmp_bank%gene, &
                            tmp_bank%E, tmp_bank%nft) 
            call ligand2bank(tmp_prot%ligand(ligand%lig_no), tmp_bank)
            tmp_bank%E = ff%ligdock
            n_soln = n_soln + 1
            n_pool = n_pool + 1
            bank_pool(n_pool) = tmp_bank
            if (tmp_bank%E(0) < e0max) then
                n_pool_e0 = n_pool_e0 + 1
                !print*, i_gene, n_pool, n_pool_e0, tmp_bank%E(0)
            end if
            if (n_soln == cand_per_trial) exit
            if (n_pool == max_candidate_pool) exit
        end do
        if (n_pool == max_candidate_pool) exit
    end do
    ! 
    if (n_pool == max_candidate_pool) exit
end do
deallocate(gene_s)
deallocate(conf_bank)

call sort_csa_bank(bank_pool(1:n_pool), n_pool)

if (n_pool_e0 > n_bank) then
    n_pool_e0 = min(n_pool_e0, n_bank*20)
    write(log_msg, '(A,I5)') ' - Selecting initial bank from conformation pool', &
                             n_pool_e0
    call log_p(log_msg, level=30)
    call select_diverse_bank(ligand, bank_pool(1:n_pool_e0), &
                             bank(n_succ+1:n_bank), n_pool_e0, n_bank-n_succ)
else
    write(log_msg, '(A)') ' - Fill initial bank by energy'
    call log_p(log_msg, level=30)
    do i = 1, n_pool
        n_succ = n_succ + 1
        bank(n_succ) = bank_pool(i)
        if (n_succ == n_bank) exit
    end do
        
    if (n_succ /= n_bank) then
        write(log_msg, '(A, 2I3)') 'Cannot make first bank', n_succ, n_bank
        call terminate_with_error(log_msg)
    end if
end if

call protein_to_R(protein)
    
deallocate(bank_pool)

end subroutine generate_bank_FFT
!-------------------------------------------------------------------------------
!subroutine generate_bank_FFT(protein, ligand)
!!-------------------------------------------------------------------------------
!type(molecule_type), intent(in) :: protein
!type(ligand_type), intent(in) :: ligand
!!
!integer, parameter :: max_frag = 4
!!
!type(molecule_type) :: tmp_prot
!type(ligdock_bank_type) :: tmp_bank
!integer :: i_trial, gene_dim
!!
!integer :: fragment(max_lig_atom, max_frag), n_atm_frag(max_frag), n_frag
!integer :: i_frag, max_size
!integer :: i
!!
!type(docking_grid_param_type) :: FFT_grid_info
!type(pharm_rec_type) :: pharm_rec(max_frag)
!type(pharm_graph_type) :: rec_graph, lig_graph
!
!
!call ligand_fragmentation(ligand, fragment, n_atm_frag, n_frag, &
!                          max_frag, max_size)
!call initialize_FFT_dock(FFT_grid_info, max_size)
!do i_frag = 1, n_frag
!    call do_FFT_docking(FFT_grid_info, ligand, fragment(:,i_frag), &
!                        n_atm_frag(i_frag), max_size, pharm_rec(i_frag))
!    write(*,'(A,I1)') '(frag_idx) ', i_frag
!    do i = 1, n_atm_frag(i_frag)
!        write(*,*) '(frag_mem) ', fragment(i,i_frag)
!    end do
!    write(*,*) '(pharm_info) ', pharm_rec(i_frag)%n_pharm
!    do i = 1, pharm_rec(i_frag)%n_pharm
!        write(*,'(A,3F10.4)') '(pharm) ', pharm_rec(i_frag)%R(:,i)
!    end do
!end do
!call finalize_FFT_dock()
!!
!tmp_prot = protein
!gene_dim = 7 + (ligand%n_br-1)
!call generate_rec_graph(pharm_rec, rec_graph, n_frag)
!do i_trial = 1, max_trial
!    call generate_conf_wo_clash(tmp_prot, ligand, tmp_bank, gene_dim)
!    call generate_lig_graph(tmp_prot%ligand(ligand%lig_no), &
!                            fragment, n_atm_frag, n_frag, lig_graph)
!    call do_pharmdock(rec_graph, lig_graph) 
!end do
!stop
!
!end subroutine generate_bank_FFT
!-------------------------------------------------------------------------------
subroutine generate_bank_rigid(protein, ligand, n_bank, bank, &
                               include_curr_conf)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type), intent(out) :: bank(:)
integer, intent(in) :: n_bank
logical, intent(in) :: include_curr_conf
!
integer :: gene_dim
integer :: i_trial, i_fail
integer :: n_succ, n_fail
type(molecule_type) :: tmp_prot
type(ligdock_bank_type) :: tmp_bank
type(ligdock_bank_type), allocatable :: bank_pool(:)
type(energy_type) :: ff
logical :: status
real(dp) :: g(3,ligand%n_atm)

gene_dim = 6 + (ligand%n_br-1)

n_succ = 0
n_fail = 0

! backup initial conformation
tmp_prot = protein
allocate(bank_pool(max_trial))

if (include_curr_conf) then
    write(log_msg, '(A, I5)') 'Include input conformation in first bank.'
    call log_p(log_msg, level=30)
    call allocate_ligdock_bank_type(tmp_bank, protein, ligand)
    call get_initial_ligand_gene(tmp_prot%ligand(ligand%lig_no), ligand, &
                                 tmp_bank%gene)
    tmp_bank%nft = 1
    call do_simplex(tmp_prot, ligand, tmp_bank%gene, tmp_bank%E, tmp_bank%nft)
    call ligand2bank(tmp_prot%ligand(ligand%lig_no), tmp_bank)

    n_succ = n_succ + 1
    bank(n_succ) = tmp_bank
    write(log_msg, '(A, I3)') 'optimize first bank: ', n_succ
    call log_p(log_msg, level=40)
end if

do i_trial = 1, max_trial
    call generate_conf_wo_clash(tmp_prot, ligand, tmp_bank, gene_dim)
    call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
    tmp_bank%E(:) = ff%ligdock(:)
    !print*, tmp_bank%E(:)

    if (tmp_bank%E(0) < e1max) then
        tmp_bank%nft = 1
        write(log_msg, '(A, I5)') 'make bank pool       :   ', i_trial
        call log_p(log_msg, level=40)

        call do_simplex(tmp_prot, ligand, tmp_bank%gene, tmp_bank%E, tmp_bank%nft, pert_full=.true.)
        call opt_ligand_steepest_descent_intcrd(tmp_bank%E, status, tmp_bank%gene, &
                                                tmp_prot, ligand, tmp_bank%nft)
        call ligand2bank(tmp_prot%ligand(ligand%lig_no), tmp_bank)

        if (tmp_bank%E(0) < e0max) then
            n_succ = n_succ + 1
            bank(n_succ) = tmp_bank
            write(log_msg, '(A, I3)') 'optimize first bank: ', n_succ
            call log_p(log_msg, level=40)
        else
            n_fail = n_fail + 1
            bank_pool(n_fail) = tmp_bank
        end if
        if (n_succ == n_bank) exit
    end if
end do

if (n_succ < n_bank) then
    call log_p('Failed to generate first bank.', level=30)
    call log_p('Fill up the rest of first bank from the failed units.', level=30)
    
    call sort_bank_to_fill(bank_pool, n_fail)
    
    do i_fail = 1, n_fail
        n_succ = n_succ + 1
        bank(n_succ) = bank_pool(i_fail)
        write(log_msg, '(A, I3)') 'Fill out first bank: ', n_succ
        call log_p(log_msg, level=30)
        if (n_succ == n_bank) exit
    end do

    if (n_succ /= n_bank) then
        write(log_msg, '(A, 2I3)') 'ERROR: Cannot make first bank', &
                                    n_succ, n_bank
        call terminate_with_error(log_msg)
    end if
end if

call protein_to_R(protein)
deallocate(bank_pool)

end subroutine generate_bank_rigid
!-------------------------------------------------------------------------------
subroutine generate_bank_flex_sc(protein, ligand, n_bank, bank,&
                                 include_curr_conf)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type), intent(out) :: bank(:)
integer, intent(in) :: n_bank
logical, intent(in) :: include_curr_conf
!
integer :: gene_dim
integer :: i_trial, i_fail
integer :: i_res, i_usc
integer :: n_succ, n_fail
integer :: i_rot, i_rot_full
type(molecule_type) :: tmp_prot
type(ligdock_bank_type) :: tmp_bank
type(ligdock_bank_type), allocatable :: bank_pool(:)
type(energy_type) :: ff
real(dp) :: g(3, ligand%n_atm)

gene_dim = 6 + (ligand%n_br-1)

n_succ = 0
n_fail = 0

! backup initial conformation
tmp_prot = protein
allocate(bank_pool(max_trial))

if (include_curr_conf) then
    write(log_msg, '(A, I5)') 'Include input conformation in first bank.'
    call log_p(log_msg, level=30)
    call allocate_ligdock_bank_type(tmp_bank, protein, ligand)
    call get_initial_ligand_gene(tmp_prot%ligand(ligand%lig_no), ligand, &
                                 tmp_bank%gene)
    do i_usc = 1, n_usc
        i_res = usc_list(i_usc)
        call identify_rotamer_state_single(tmp_prot, i_res, i_rot, &
                                           i_rot_full)
        tmp_bank%rot_idx(i_usc) = i_rot
        tmp_bank%flex_sc(i_usc)%R = tmp_prot%residue(i_res)%R
        tmp_bank%flex_sc(i_usc)%n_atm = tmp_prot%residue(i_res)%n_atm
    end do
    !
    call protein_to_R(protein)
    call set_atdk_hbond_flex(usc_list, n_usc)
    !
    tmp_bank%nft = 1
    call do_simplex(tmp_prot, ligand, tmp_bank%gene, tmp_bank%E, tmp_bank%nft)
    call update_prot_prot_intrxn(tmp_bank%rot_idx(1:n_usc), tmp_bank%E)
    call ligand2bank(tmp_prot%ligand(ligand%lig_no), tmp_bank)

    n_succ = n_succ + 1
    bank(n_succ) = tmp_bank
    write(log_msg, '(A, I3)') 'optimize first bank: ', n_succ
    call log_p(log_msg, level=40)
end if

do i_trial = 1, max_trial
    call generate_conf_wo_clash(tmp_prot, ligand, tmp_bank, gene_dim)
    call generate_prot_conf(tmp_bank, tmp_prot, ligand)
    !
    call protein_to_R(tmp_prot)
    call set_atdk_hbond_flex(usc_list, n_usc)
    call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
    call update_prot_prot_intrxn(tmp_bank%rot_idx(1:n_usc), ff%ligdock)
    tmp_bank%E(:) = ff%ligdock(:)

    if (tmp_bank%E(0) < e1max) then
        tmp_bank%nft = 1
        write(log_msg, '(A, I5)') 'make bank pool       :   ', i_trial
        call log_p(log_msg, level=40)

        call do_simplex(tmp_prot, ligand, tmp_bank%gene, tmp_bank%E, tmp_bank%nft)
        call update_prot_prot_intrxn(tmp_bank%rot_idx(1:n_usc), tmp_bank%E)
        call ligand2bank(tmp_prot%ligand(ligand%lig_no), tmp_bank)

        if (tmp_bank%E(0) < e0max) then
            n_succ = n_succ + 1
            bank(n_succ) = tmp_bank
            write(log_msg, '(A, I3)') 'optimize first bank: ', n_succ
            call log_p(log_msg, level=40)
        else
            n_fail = n_fail + 1
            bank_pool(n_fail) = tmp_bank
        end if
        if (n_succ == n_bank) exit
    end if
end do

if (n_succ < n_bank) then
    call log_p('Failed to generate first bank.', level=30)
    call log_p('Fill up the rest of first bank from the failed units.', level=30)

    call sort_bank_to_fill(bank_pool, n_fail)
    
    do i_fail = 1, n_fail
        n_succ = n_succ + 1
        bank(n_succ) = bank_pool(i_fail)
        write(log_msg, '(A, I3)') 'Fill out first bank: ', n_succ
        call log_p(log_msg, level=30)
        if (n_succ == n_bank) exit
    end do

    if (n_succ /= n_bank) then
        write(log_msg, '(A, 2I3)') 'ERROR: Cannot make first bank', &
                                    n_succ, n_bank
        call terminate_with_error(log_msg)
    end if
end if

call protein_to_R(protein)
deallocate(bank_pool)

end subroutine generate_bank_flex_sc
!-------------------------------------------------------------------------------
subroutine generate_bank_VD_rigid(protein, ligand, n_bank, bank)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type), intent(out) :: bank(:)
integer, intent(in) :: n_bank
!
type(molecule_type) :: tmp_prot
type(ligdock_bank_type), allocatable :: bank_pool(:)
type(ligdock_bank_type), allocatable :: bank_pool_rand(:)
type(energy_type) :: ff
integer :: i, n_pool, n_pool_e0, i_succ, n_left
real(dp) :: g(3,ligand%n_atm)

! Read initial docking results
allocate(bank_pool(max_trial))
do i = 1, max_trial
    call allocate_ligdock_bank_type(bank_pool(i), protein, ligand)
end do
tmp_prot = protein
call read_VD_conformations(tmp_prot, ligand, bank_pool)

! Run energy minimization
n_pool = 0
n_pool_e0 = 0
tmp_prot = protein

do i = 1, max_trial
    call bank2ligand(bank_pool(i), tmp_prot%ligand(ligand%lig_no))
    call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
    call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)

    bank_pool(i)%E(0) = ff%ligdock(0)
    if (bank_pool(i)%E(0) < e1max) then
        n_pool = n_pool + 1
        bank_pool(i)%nft = 1
        write(log_msg, '(A, I5)') 'make bank pool       :   ', i
        call log_p(log_msg, level=40)

        call do_simplex(tmp_prot, ligand, bank_pool(i)%gene, &
                        bank_pool(i)%E, bank_pool(i)%nft)
        call ligand2bank(tmp_prot%ligand(ligand%lig_no), bank_pool(i))
        if (bank_pool(i)%E(0) < e0max) then
            n_pool_e0 = n_pool_e0 + 1
        end if
    end if
end do

call sort_bank_to_fill(bank_pool, max_trial)

if (n_pool_e0 > n_bank) then
    n_pool_e0 = min(n_pool_e0, n_bank*10)
    call select_diverse_bank(ligand, bank_pool(1:n_pool_e0), &
                             bank(1:n_bank), n_pool_e0, n_bank)
else
    i_succ = 0
    do i = 1, n_pool
        i_succ = i_succ + 1
        bank(i_succ) = bank_pool(i)
        if (i_succ == n_bank) exit
    end do
        
    if (i_succ /= n_bank) then
        write(log_msg, '(A, 2I3)') 'WARNING: Faild to fill first bank with VD results, ', i_succ, n_bank
        call log_p(log_msg, level=30)
        !
        write(log_msg, '(A, 2I3)') 'WARNING: Try to fill using random sampling.'
        call log_p(log_msg, level=30)
        !
        ! try to fill by random generation
        n_left = n_bank - i_succ
        allocate(bank_pool_rand(n_left))
        do i = 1, n_left
            call allocate_ligdock_bank_type(bank_pool_rand(i), protein, ligand)
        end do
        call generate_bank_rigid(protein, ligand, n_left, &
                                 bank_pool_rand(:), &
                                 .false.)
        bank(i_succ+1:n_bank) = bank_pool_rand(1:n_left)
        do i = 1, n_left
            call deallocate_ligdock_bank_type(bank_pool_rand(i))
        end do
        deallocate(bank_pool_rand)
    end if
end if

do i = 1, max_trial
    call deallocate_ligdock_bank_type(bank_pool(i))
end do
deallocate(bank_pool)

end subroutine generate_bank_VD_rigid
!-------------------------------------------------------------------------------
subroutine generate_bank_VD_flex_sc(protein, ligand, n_bank, bank)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type), intent(out) :: bank(:)
integer, intent(in) :: n_bank
!
type(molecule_type) :: tmp_prot
type(ligdock_bank_type), allocatable :: bank_pool(:)
type(energy_type) :: ff
integer :: i, n_pool, i_succ
real(dp) :: g(3, ligand%n_atm)

! Read initial docking results
allocate(bank_pool(max_trial))
do i = 1, max_trial
    call allocate_ligdock_bank_type(bank_pool(i), protein, ligand)
end do
tmp_prot = protein
call read_VD_conformations(tmp_prot, ligand, bank_pool)

! Run energy minimization
n_pool = 0
tmp_prot = protein
do i = 1, max_trial
    call bank2ligand(bank_pool(i), tmp_prot%ligand(ligand%lig_no))
    call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
    !
    call generate_prot_conf(bank_pool(i), tmp_prot, ligand)
    call protein_to_R(tmp_prot)
    call set_atdk_hbond_flex(usc_list, n_usc)
    call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
    call update_prot_prot_intrxn(bank_pool(i)%rot_idx(1:n_usc), ff%ligdock)
    bank_pool(i)%E = ff%ligdock
    !
    if (bank_pool(i)%E(0) < e1max) then
        n_pool = n_pool + 1
        bank_pool%nft = 1
        write(log_msg, '(A, I5)') 'make bank pool       :   ', i
        call log_p(log_msg, level=40)
        !
        call do_simplex(tmp_prot, ligand, bank_pool(i)%gene, &
                        bank_pool(i)%E, bank_pool(i)%nft)
        call update_prot_prot_intrxn(bank_pool(i)%rot_idx(1:n_usc), bank_pool(i)%E)
        call ligand2bank(tmp_prot%ligand(ligand%lig_no), bank_pool(i))
    end if
end do

call sort_bank_to_fill(bank_pool, max_trial)

i_succ = 0
do i = 1, n_pool
    i_succ = i_succ + 1
    bank(i_succ) = bank_pool(i)
    if (i_succ == n_bank) exit
end do
    
if (i_succ /= n_bank) then
    write(log_msg, '(A, 2I3)') 'Cannot make first bank', i_succ, n_bank
    call terminate_with_error(log_msg)
end if

call deallocate_molecule_type(tmp_prot)
do i = 1, max_trial
    call deallocate_ligdock_bank_type(bank_pool(i))
end do
deallocate(bank_pool)

end subroutine generate_bank_VD_flex_sc
!-------------------------------------------------------------------------------
subroutine read_VD_conformations(protein, ligand, bank_pool)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type), intent(inout) :: bank_pool(:)
integer :: i_bank, i_conf, i_voro
character(len=4) :: voro_num
character(len=len_fname) :: file_name

file_name = trim(voro_prefix) // '.pdbq'
call get_predock_conf(file_name, bank_pool, protein, ligand)

end subroutine read_VD_conformations
!-------------------------------------------------------------------------------
subroutine get_predock_conf(file_name, bank, protein, ligand)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: file_name
type(ligdock_bank_type), intent(inout) :: bank(:)
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(in) :: ligand
!
integer :: i_unit, ioerr
character(len=len_fname) :: line
!
logical :: trans_read, rot_read
real(dp) :: tor_ang(ligand%n_br-1)
integer :: i_atm, i_tor, cntr_atom, i_bank
real(dp) :: trans(3), euler_ang(3), angle, axis(3)

i_unit = 19
open(unit=i_unit, file=trim(file_name), status='old', action='read', &
     iostat=ioerr)
if (ioerr /= 0) then
    write(log_msg, '(A,A)') 'Cannot open pdbq file: ', trim(file_name)
    call terminate_with_error(log_msg)
end if

i_bank = 0

do
    read(i_unit, '(A120)', iostat=ioerr) line
    if (ioerr/=0) exit
    if (line(1:5) == "MODEL") then
        i_bank = i_bank + 1
        trans_read = .false.
        rot_read = .false.
        i_tor = 0
        bank(i_bank)%gene_dim = 7 + (ligand%n_br-1)
        bank(i_bank)%nft = 0
    else if (line(8:12) == 'Trans' .and. (.not. trans_read)) then
        read(line(14:),*) trans(1:3)
        trans_read = .true.
        !
        cntr_atom = ref_lig(ligand%lig_type)%cntr_atm
        bank(i_bank)%gene(1:3) = trans(1:3) - ref_lig(ligand%lig_type)%R(:,cntr_atom)
    else if (line(8:10) == 'Rot' .and. (.not. rot_read)) then
        read(line(14:),*) euler_ang(1:3)
        rot_read = .true.
        call euler_to_quaternion(euler_ang(1:3), angle, axis(1:3))
        call v_norm(axis(1:3))
        bank(i_bank)%gene(4) = angle
        bank(i_bank)%gene(5:7) = axis(1:3)
    else if (line(8:10) == 'Tor') then
        i_tor = i_tor + 1
        read(line(18:),*) tor_ang(i_tor)
        bank(i_bank)%gene(i_tor+7) = bound_ang(tor_ang(i_tor))
    !else if (line(1:4) == 'ATOM') then
    !    read(line(7:11), *) i_atm
    !    if (line(36:38) == '-na') then
    !        bank(i_bank)%R(:,i_atm) = 0.0d0
    !    else
    !        read(line(31:38), '(F8.3)') bank(i_bank)%R(1,i_atm)
    !        read(line(39:46), '(F8.3)') bank(i_bank)%R(2,i_atm)
    !        read(line(47:54), '(F8.3)') bank(i_bank)%R(3,i_atm)
    !    end if
    end if
end do
close(i_unit)
do i_bank = 1, n_gen_conf*n_voro_per_conf
    if (isnan(bank(i_bank)%gene(1))) then
        print*, 'nan occured, re-sample gene'
        call construct_gene_random(ref_lig(ligand%lig_type), bank(i_bank)%gene_dim, &
                                   bank(i_bank)%gene(1:bank(i_bank)%gene_dim))
    end if
    call construct_ligand(protein, ligand, bank(i_bank)%gene)
    call ligand2bank(protein%ligand(ligand%lig_no), bank(i_bank))
end do

end subroutine get_predock_conf
!-------------------------------------------------------------------------------
subroutine generate_conf_wo_clash(protein, ligand, conf_bank, gene_dim)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type), intent(out) :: conf_bank
integer, intent(in) :: gene_dim
!
integer :: clash
integer, parameter :: allow_bump = 1

call allocate_ligdock_bank_type(conf_bank, protein, ligand)

do
    call construct_gene_random(ref_lig(ligand%lig_type), gene_dim, &
                               conf_bank%gene(1:gene_dim))
    call construct_ligand(protein, ligand, conf_bank%gene(1:gene_dim))
    call update_R_for_ligand(protein%ligand(ligand%lig_no), ligand%lig_no)
    call check_clash_within_ligand(ligand, clash)
    if (clash <= allow_bump) exit
end do

end subroutine generate_conf_wo_clash
!-------------------------------------------------------------------------------
subroutine generate_conf_lib(protein, ligand, conf_bank, n_bank, gene_dim)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type), intent(out) :: conf_bank(:)
integer, intent(in) :: n_bank, gene_dim
!
type(ligdock_bank_type), allocatable :: conf_pool(:), failed_bank(:)
type(ligdock_bank_type) :: tmp_bank, min_bank
integer :: n_conf_pool, n_trial, i_trial, i_dim, i_iter, i_conf
integer :: n_succ, n_fail
real(dp) :: pert_ang
!
real(dp) :: atdk_w0
logical :: use_drugscore0, use_Xscore0, use_ROTA0, prot_w0
integer :: clash, status, ierr
integer, parameter :: allow_bump = 0
real(dp), parameter :: int_E_cut = 20.0d0
real(dp), parameter :: int_E_cut_2 = 50.0d0*int_E_cut
type(energy_type) :: ff
real(dp) :: g(3,ligand%n_atm)
!
type(molecule_type) :: tmp_prot
integer :: cntr_idx

!print*, 'Generate conformation library...'
atdk_w0 = atdk_w
use_drugscore0 = use_drugscore
use_Xscore0 = use_Xscore
use_ROTA0 = use_ROTA
prot_w0 = prot_w
!
atdk_w = 0.0d0
prot_w = 0.0d0
use_drugscore = .false.
use_Xscore = .false.
use_ROTA = .false.

n_conf_pool = n_bank*3
n_trial = n_bank*20
allocate(conf_pool(n_conf_pool))
allocate(failed_bank(n_trial))

n_succ = 0
n_fail = 0
call allocate_ligdock_bank_type(tmp_bank, protein, ligand)
call allocate_ligdock_bank_type(min_bank, protein, ligand)
tmp_prot = protein
do i_trial = 1, n_trial
    call sample_torsion_wo_clash(tmp_prot, ligand, tmp_bank, gene_dim)
    !
    !call get_initial_ligand_gene(tmp_prot%ligand(ligand%lig_no), ligand, &
    !                             tmp_bank%gene)
    cntr_idx = ref_lig(ligand%lig_type)%cntr_atm
    do i_dim = 1, 3
        tmp_bank%gene(i_dim) = dock_grid_info%grid_cntr(i_dim) - ref_lig(ligand%lig_type)%R(i_dim, cntr_idx)
    end do
    tmp_bank%gene(4:6) = (/0.0d0, 0.0d0, 0.0d0/)
    !do i_dim = 8, gene_dim
    !    pert_ang = (10.0d0 * random() - 5.0d0) * deg2rad
    !    tmp_bank%gene(i_dim) = tmp_bank%gene(i_dim) + pert_ang
    !    tmp_bank%gene(i_dim) = bound_ang(tmp_bank%gene(i_dim))
    !end do
    call construct_ligand(tmp_prot, ligand, tmp_bank%gene(1:gene_dim))
    call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
    !
    call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
    !print*, 'check: ', i_trial, ff%ligdock(0), int_E_cut_2
    if (ff%ligdock(0) < int_E_cut_2) then
        min_bank%gene(:) = tmp_bank%gene(:)
        min_bank%E(:) = ff%ligdock(:)
        min_bank%R(:, 1:ligand%n_atm) = &
            tmp_prot%ligand(ligand%lig_no)%R(:,1:ligand%n_atm)
        do i_dim = 7, gene_dim ! for torsion angles
            tmp_bank = min_bank
            do i_iter = 1, 20
                pert_ang = (10.0d0 * random() - 5.0d0) * deg2rad
                tmp_bank%gene(i_dim) = min_bank%gene(i_dim) + pert_ang
                tmp_bank%gene(i_dim) = bound_ang(tmp_bank%gene(i_dim))
                call construct_ligand(tmp_prot, ligand, tmp_bank%gene(1:gene_dim))
                call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
                call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
                if (ff%ligdock(0) < min_bank%E(0)) then
                    min_bank%E(:) = ff%ligdock(:)
                    min_bank%R(:, 1:ligand%n_atm) = &
                        protein%ligand(ligand%lig_no)%R(:,1:ligand%n_atm)
                    min_bank%gene(:) = tmp_bank%gene(:)
                end if
            end do
        end do
        if (min_bank%E(0) < int_E_cut) then
            n_succ = n_succ + 1
            conf_pool(n_succ) = min_bank
        else
            n_fail = n_fail + 1
            failed_bank(n_fail) = min_bank
        end if
        !
        if (n_succ == n_conf_pool) exit
    end if
end do

if (n_succ < n_bank) then
    call sort_csa_bank(failed_bank, n_fail)
    do i_conf = 1, n_fail
        n_succ = n_succ + 1
        conf_pool(n_succ) = failed_bank(i_conf)
        if (n_succ == n_bank) exit
    end do
    if (n_succ < n_bank) then
        write(log_msg, '(A, 2I4)') 'Cannot make enough ligand conformations', &
                                    n_succ, n_conf_pool
        call terminate_with_error(log_msg)
    end if
    return
end if

call select_diverse_bank(ligand, conf_pool(1:n_succ), &
                         conf_bank(1:n_bank), n_succ, n_bank)

deallocate(conf_pool, failed_bank)
!
atdk_w = atdk_w0
use_drugscore = use_drugscore0
use_Xscore = use_Xscore0
use_ROTA = use_ROTA0
prot_w = prot_w0

call deallocate_ligdock_bank_type(tmp_bank)
call deallocate_ligdock_bank_type(min_bank)

!print*, 'Done.'
end subroutine generate_conf_lib
!-------------------------------------------------------------------------------
subroutine construct_gene_random(ref_lig, gene_dim, gene)
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(in) :: ref_lig
integer, intent(in) :: gene_dim
real(dp), intent(out) :: gene(:)
!
integer :: i
real(dp) :: half_box(3), bound_1(3), bound_2(3)
real(dp) :: angle, axis(3)

! set boundary to locate ligand within the grid box
half_box(:) = dble((dock_grid_info%n_elem(:)-1)/2) * dock_grid_info%grid_width
if (use_rsr) then
    half_box(:) = half_box(:) * 0.10d0
end if

do i = 1, 3
    bound_1(i) = dock_grid_info%grid_cntr(i) + half_box(i)
    bound_2(i) = dock_grid_info%grid_cntr(i) - half_box(i)
end do

! generate random population
! 1. generate gene for translation
do i = 1,3
    gene(i) = (bound_1(i) - bound_2(i))*random() + bound_2(i)
    gene(i) = gene(i) - ref_lig%R(i,ref_lig%cntr_atm)
end do

! 2. generate gene for rotation
angle = two_pi * random()
call gen_random_axis(axis)
gene(4:6) = axis(1:3)*angle

! 3. generate gene for torsions
do i = 7, gene_dim
    gene(i) = two_pi*random() - pi
end do

!-------------------------------------------------------------------------------
end subroutine construct_gene_random
!-------------------------------------------------------------------------------
subroutine generate_prot_conf(bank, protein, ligand)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(inout) :: bank
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(in) :: ligand
!
integer :: i_trial, n_trial, i_res, i_usc, i_rot, i_rot_full
integer :: i_start, i_end, n_rot
integer :: rot_idx(n_usc)
type(molecule_type) :: tmp_prot
type(energy_type) :: ff
type(ligdock_bank_type) :: min_bank
integer :: rotamer_start_num(2,n_usc)
real(dp) :: g(3, ligand%n_atm)

do i_usc = 1, n_usc
    i_res = usc_list(i_usc)
    call get_rotamer_index(protein%residue(i_res:i_res+1), i_start, i_end)
    rotamer_start_num(1,i_usc) = i_start
    rotamer_start_num(2,i_usc) = i_end - i_start + 1
end do

n_trial = n_usc * 10

tmp_prot = protein
min_bank = bank
do i_trial = 1, n_trial
    do i_usc = 1, n_usc
        i_res = usc_list(i_usc)
        n_rot = rotamer_start_num(2,i_usc)
        i_rot = int(random()*n_rot)
        i_rot_full = rotamer_start_num(1,i_usc) + i_rot
        call place_rotamer(tmp_prot, i_res, i_rot_full)
        rot_idx(i_usc) = i_rot + 1
        call update_R_for_sc(tmp_prot%residue(i_res), i_res)
    end do
    call set_atdk_hbond_flex(usc_list, n_usc)
    call ligdock_energy_using_grid(tmp_prot, ligand, ff,g, .false.)
    call update_prot_prot_intrxn(rot_idx, ff%ligdock)
    if (i_trial == 1) then
        min_bank%E(:) = ff%ligdock(:)
        do i_usc = 1, n_usc
            i_res = usc_list(i_usc)
            min_bank%rot_idx(i_usc) = rot_idx(i_usc)
            min_bank%flex_sc(i_usc)%R = tmp_prot%residue(i_res)%R
            min_bank%flex_sc(i_usc)%n_atm = tmp_prot%residue(i_res)%n_atm
            protein%residue(i_res) = tmp_prot%residue(i_res)
        end do
        cycle
    end if
    if (ff%ligdock(0) < min_bank%E(0)) then
        min_bank%E(:) = ff%ligdock(:)
        do i_usc = 1, n_usc
            i_res = usc_list(i_usc)
            min_bank%rot_idx(i_usc) = rot_idx(i_usc)
            min_bank%flex_sc(i_usc)%R = tmp_prot%residue(i_res)%R
            min_bank%flex_sc(i_usc)%n_atm = tmp_prot%residue(i_res)%n_atm
            protein%residue(i_res) = tmp_prot%residue(i_res)
        end do
    end if
end do
    
bank = min_bank

end subroutine generate_prot_conf
!-------------------------------------------------------------------------------
subroutine sort_bank_to_fill(bank, n_bank)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(inout) :: bank(:)
integer, intent(in) :: n_bank
integer :: i_bank, j_bank
integer :: min_loc
real(dp) :: min_E
type(ligdock_bank_type) :: first_bank, tmp_bank

do i_bank = 1, n_bank
    first_bank = bank(i_bank)
    min_E = bank(i_bank)%E(0)
    min_loc = i_bank

    do j_bank = i_bank, n_bank
        if(bank(j_bank)%E(0) < min_E) then
            min_loc = j_bank
            min_E = bank(j_bank)%E(0)
        end if
    end do
    tmp_bank = bank(min_loc)
    bank(min_loc) = first_bank
    bank(i_bank) = tmp_bank
end do

end subroutine sort_bank_to_fill
!-------------------------------------------------------------------------------
subroutine select_diverse_bank(ligand, bank_pool, bank, n_pool, n_bank)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type), intent(in) :: bank_pool(:)
type(ligdock_bank_type), intent(out) :: bank(:)
integer, intent(in) :: n_pool, n_bank
!
integer :: i_conf, j_conf
!
character(len=20) :: atom_types(ligand%n_atm)
integer :: n_type, n_atm_s(ligand%n_atm), atm_id_s(ligand%n_atm, ligand%n_atm)
real(dp) :: R_mat(n_pool, n_pool), rmsd
!
real(dp) :: energy(n_pool), min_E
integer :: used_list(n_pool), n_used
integer :: min_loc(1), selected_conf

! Calc distance matrix
call get_xlogP_atom_types(ref_lig(ligand%lig_type), atom_types)
call prep_RMSD_PM_calculation(atom_types, ligand%n_atm, n_type, n_atm_s, atm_id_s)
do i_conf = 1, n_pool - 1
    do j_conf = i_conf + 1, n_pool
        call calc_RMSD_PM(bank_pool(i_conf)%R, bank_pool(j_conf)%R, &
                          n_type, n_atm_s, atm_id_s, rmsd)
        R_mat(i_conf, j_conf) = rmsd
        R_mat(j_conf, i_conf) = rmsd
    end do
end do

! Fill energy list & find minimum E conf
do i_conf = 1, n_pool
    energy(i_conf) = bank_pool(i_conf)%E(0)
end do

min_E = minval(energy(1:n_pool))
min_loc = minloc(energy(1:n_pool))

! Fill bank with lowest E conformation
bank(1) = bank_pool(min_loc(1))
used_list(1) = min_loc(1)
n_used = 1

! Find farthest conformation from selected conformations and select it.
do i_conf = 2, n_bank
    call find_farthest_member(R_mat, energy, used_list, n_used, n_pool, &
                              selected_conf)
    used_list(i_conf) = selected_conf
    bank(i_conf) = bank_pool(selected_conf)
    n_used = n_used + 1
end do

end subroutine select_diverse_bank
!-------------------------------------------------------------------------------
subroutine find_farthest_member(R_mat, energy, used_list, n_used, n_pool, &
                                selected_idx)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: R_mat(:,:), energy(:)
integer, intent(in) :: used_list(:)
integer, intent(in) :: n_used, n_pool
integer, intent(out) :: selected_idx
!
logical :: is_used(n_pool)
integer :: i_conf, i_used
real(dp) :: max_dist, max_conf_E, dist
integer :: max_loc

is_used(:) = .false.
do i_conf = 1, n_used
    is_used(used_list(i_conf)) = .true.
end do

max_dist = 0.0
max_loc = 0
max_conf_E = e0max

do i_conf = 1, n_pool
    if (is_used(i_conf)) cycle

    dist = 0.0
    do i_used = 1, n_used
        dist = dist + R_mat(used_list(i_used), i_conf)
    end do

    if (dist > max_dist) then
        max_dist = dist
        max_loc = i_conf
        max_conf_E = energy(i_conf)
    else if (dist == max_dist) then
        if (energy(i_conf) < max_conf_E) then
            max_dist = dist
            max_loc = i_conf
            max_conf_E = energy(i_conf)
        end if
    end if
end do

selected_idx = max_loc

end subroutine find_farthest_member
!-------------------------------------------------------------------------------
END MODULE LIGDOCK_CSA_GEN_BANK
!-------------------------------------------------------------------------------
