!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/sampling/optimizers/ligdock_csa/ligdock_csa_new_conf.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE LIGDOCK_CSA_NEW_CONF
!-------------------------------------------------------------------------------

use globals
use logger
use ran
use geometry, only: cartesian2internal, pert_rotation, expmap_quat, inv_expmap_quat
use mathfunctions, only: v_norm, bound_ang
use rmsd, only: calc_rmsd_simple, calc_rmsd_PM
use rotamer, only: get_rotamer_index, place_rotamer
!
use energy_utils, only: protein_to_R, update_R_for_sc, update_R_for_ligand
use ligdock_energy, only: ligdock_energy_using_grid, set_atdk_hbond_flex, &
                          update_prot_prot_intrxn
use gradient_minimization, only: opt_ligand_steepest_descent_intcrd
!
use simplex, only: do_simplex
!
use ligdock_csa_vars
use ligdock_csa_utils
use ligdock_csa_dist, only: atom_types, n_atm_s, n_type, atm_id_s

use ligdock_csa_fragFFT

implicit none
save
private

public :: make_new_conf
public :: optimize_new_conf_rigid
public :: optimize_seed_rigid
public :: optimize_new_conf_flex_sc
public :: optimize_seed_flex_sc
public :: optimize_bank_steepest

CONTAINS
!-------------------------------------------------------------------------------
subroutine make_new_conf(protein, ligand, bank, bank_1st, bank_new, &
                         n_bank, n_bank_add, i_seed_cycle)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type), intent(in) :: bank(:), bank_1st(:)
type(ligdock_bank_type), intent(out) :: bank_new(:)
integer, intent(in) :: n_bank, n_bank_add, i_seed_cycle
!
integer :: i_seed, i_opr, i_seed_selected, i_partner
integer :: i_new, n_pert

call log_divider()
call log_p("Make new conformation", level = 20)

310 format (6X,A,2X,I4,2X,I4)

i_new = 0
n_pert = 0
if (n_opr_1 > 0) then
    do i_seed = 1, n_seed
        do i_opr = 1, n_opr_1
            call select_seed_and_partner(i_seed, i_seed_selected, i_partner, &
                                         n_bank, n_bank_add, i_seed_cycle)
            write(log_msg, 310) "Operator 1 (Seed B / Partner R) =", &
                                i_seed_selected, i_partner
            call log_p(log_msg, level=40)
            i_new = i_new + 1
            !call operator_1(bank(i_seed_selected), bank_1st(i_partner), &
            !                bank_new(i_new), n_pert, ligand%n_atm)
            call operator_1_1(bank(i_seed_selected), bank_1st(i_partner), &
                            bank_new(i_new), protein, ligand)
        end do
    end do
end if

if (n_opr_2 > 0) then
    do i_seed = 1, n_seed
        do i_opr = 1, n_opr_2
            call select_seed_and_partner(i_seed, i_seed_selected, i_partner, &
                                         n_bank, n_bank_add, i_seed_cycle)
            write(log_msg, 310) "Operator 2 (Seed B / Partner R) =", &
                                i_seed_selected, i_partner
            call log_p(log_msg, level=40)
            i_new = i_new + 1
            !call operator_2(bank(i_seed_selected), bank(i_partner), &
            !                bank_new(i_new), n_pert, ligand%n_atm)
            call operator_2_1(bank(i_seed_selected), bank(i_partner), &
                            bank_new(i_new), protein, ligand)
        end do
    end do
end if

if (n_opr_3 > 0) then
    do i_seed = 1, n_seed
        do i_opr = 1, n_opr_3
            call select_seed_and_partner(i_seed, i_seed_selected, i_partner, &
                                         n_bank, n_bank_add, i_seed_cycle)
            write(log_msg, 310) "Operator 3 (Seed B / Partner R) =", &
                                i_seed_selected, i_partner
            call log_p(log_msg, level=40)
            i_new = i_new + 1
            call operator_3(bank(i_seed_selected), bank_1st(i_partner), &
                            bank_new(i_new), ligand%n_br-1)
            !call operator_3_1(bank(i_seed_selected), bank_1st(i_partner), &
            !                bank_new(i_new), protein, ligand, i_seed_cycle)
            !call operator_3_2(bank(i_seed_selected), bank_1st(i_partner), &
            !                bank_new(i_new), ligand%n_br-1, protein, ligand)
        end do
    end do
end if

if (n_opr_4 > 0) then
    do i_seed = 1, n_seed
        do i_opr = 1, n_opr_4
            call select_seed_and_partner(i_seed, i_seed_selected, i_partner, &
                                         n_bank, n_bank_add, i_seed_cycle)
            write(log_msg, 310) "Operator 4 (Seed B / Partner R) =", &
                                i_seed_selected, i_partner
            call log_p(log_msg, level=40)
            i_new = i_new + 1
            call operator_4(bank(i_seed_selected), bank(i_partner), &
                            bank_new(i_new), ligand%n_br-1)
            !call operator_4_1(bank(i_seed_selected), bank(i_partner), &
            !                bank_new(i_new), protein, ligand, i_seed_cycle)
            !call operator_4_2(bank(i_seed_selected), bank(i_partner), &
            !                bank_new(i_new), ligand%n_br-1, protein, ligand)
        end do
    end do
end if

if (n_opr_5 > 0) then
    do i_seed = 1, n_seed
        do i_opr = 1, n_opr_5
            call select_seed_and_partner(i_seed, i_seed_selected, i_partner, &
                                         n_bank, n_bank_add, i_seed_cycle)
            write(log_msg, 310) "Operator 5 (Seed B / Partner R) =", &
                                i_seed_selected, i_partner
            call log_p(log_msg, level=40)
            i_new = i_new + 1
            call operator_5(bank(i_seed_selected), bank_1st(i_partner), &
                            bank_new(i_new), protein, ligand, i_seed_cycle)
        end do
    end do
end if

if (n_opr_6 > 0) then
    do i_seed = 1, n_seed
        do i_opr = 1, n_opr_6
            call select_seed_and_partner(i_seed, i_seed_selected, i_partner, &
                                         n_bank, n_bank_add, i_seed_cycle)
            write(log_msg, 310) "Operator 6 (Seed B / Partner R) =", &
                                i_seed_selected, i_partner
            call log_p(log_msg, level=40)
            i_new = i_new + 1
            call operator_6(bank(i_seed_selected), bank(i_partner), &
                            bank_new(i_new), protein, ligand, i_seed_cycle)
        end do
    end do
end if

if (n_opr_7 > 0) then
    do i_seed = 1, n_seed
        do i_opr = 1, n_opr_7
            call select_seed_and_partner(i_seed, i_seed_selected, i_partner, &
                                         n_bank, n_bank_add, i_seed_cycle)
            write(log_msg, 310) "Operator 6 (Seed B / Partner R) =", &
                                i_seed_selected, i_partner
            call log_p(log_msg, level=40)
            i_new = i_new + 1
            call operator_7(bank(i_seed_selected), bank_1st(i_partner), &
                            bank_new(i_new), protein)
        end do
    end do
end if

if (n_opr_8 > 0) then
    do i_seed = 1, n_seed
        do i_opr = 1, n_opr_8
            call select_seed_and_partner(i_seed, i_seed_selected, i_partner, &
                                         n_bank, n_bank_add, i_seed_cycle)
            write(log_msg, 310) "Operator 7 (Seed B / Partner R) =", &
                                i_seed_selected, i_partner
            call log_p(log_msg, level=40)
            i_new = i_new + 1
            call operator_8(bank(i_seed_selected), bank(i_partner), &
                            bank_new(i_new), protein)
        end do
    end do
end if

if (n_opr_1 /= 0 .or. n_opr_2 /= 0) then
    write(log_msg, '(A, f8.3)')  'perturbation rate      :', &
          dble(n_pert)/(dble(n_seed*(n_opr_1+n_opr_2)))
    call log_p(log_msg, level=30)
end if

end subroutine make_new_conf
!-------------------------------------------------------------------------------
subroutine select_seed_and_partner(i_seed, i_seed_selected, i_partner, &
                                   n_bank, n_bank_add, i_seed_cycle)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_seed, n_bank, n_bank_add, i_seed_cycle
integer, intent(out) :: i_seed_selected, i_partner

! select seed
i_seed_selected = idx_seed(i_seed)
if (i_seed_cycle == 0 .and.  n_unused > n_bank_add-irr) then
    i_seed_selected = n_bank - n_bank_add + random()*n_bank_add + 1
    if (i_seed_selected > n_bank) i_seed_selected = n_bank
end if

! select partner        
do 
    i_partner = n_bank*random()+1
    if (i_partner > n_bank) i_partner = n_bank
    if (i_seed_selected /= i_partner) exit
end do

end subroutine select_seed_and_partner
!-------------------------------------------------------------------------------
subroutine operator_1(seed_u, partner_u, new_u, n_pert, n_lig_atm)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: seed_u, partner_u
type(ligdock_bank_type), intent(out) :: new_u
integer, intent(inout) :: n_pert
integer, intent(in) :: n_lig_atm
!
real(dp), parameter :: len_tr = 0.5d0
real(dp), parameter :: len_ang = 15.0d0 * deg2rad
real(dp) :: dist

new_u = seed_u

! initialize nft
new_u%nft = 0
! trans/rot cross over
new_u%gene(1:6) = partner_u%gene(1:6)

! If seed_u & partner_u aren't similar, perturb trans/rot
!call calc_rmsd_simple(seed_u%R(:,1:n_lig_atm), partner_u%R(:,1:n_lig_atm), &
!                      n_lig_atm, dist)
call calc_rmsd_PM(seed_u%R(:,1:n_lig_atm), partner_u%R(:,1:n_lig_atm), &
                  n_type, n_atm_s, atm_id_s, dist)
if (dist > mutation_cut * D_ave) then
    n_pert = n_pert + 1
    call perturb_trans_rot_gene(new_u%gene, len_tr, len_ang)
end if

end subroutine operator_1
!-------------------------------------------------------------------------------
subroutine operator_1_1(seed_u, partner_u, new_u, protein, ligand)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: seed_u, partner_u
type(ligdock_bank_type), intent(out) :: new_u
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
!
type(ligdock_bank_type) :: tmp_u
type(molecule_type) :: tmp_prot
type(energy_type) :: ff
real(dp) :: g(3,ligand%n_atm)
real(dp), parameter :: len_tr = 0.5d0
real(dp), parameter :: len_ang = 15.0d0 * deg2rad
real(dp) :: dist
integer :: n_lig_atm

new_u = seed_u
tmp_prot = protein

! initialize nft
new_u%nft = 0
! trans/rot cross over
new_u%gene(1:6) = partner_u%gene(1:6)
!
call construct_ligand(tmp_prot, ligand, new_u%gene)
call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
new_u%E(:) = ff%ligdock(:)

!! try T_crossover
tmp_u = seed_u
tmp_u%gene(1:3) = partner_u%gene(1:3)
call construct_ligand(tmp_prot, ligand, tmp_u%gene)
call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
tmp_u%E(:) = ff%ligdock(:)
if (tmp_u%E(0) < new_u%E(0)) then
    new_u = tmp_u
end if

! try R_crossover
tmp_u = seed_u
tmp_u%gene(4:6) = partner_u%gene(4:6)
call construct_ligand(tmp_prot, ligand, tmp_u%gene)
call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
tmp_u%E(:) = ff%ligdock(:)
if (tmp_u%E(0) < new_u%E(0)) then
    new_u = tmp_u
end if
!
! If seed_u & partner_u aren't similar, perturb trans/rot
n_lig_atm = ligand%n_atm
!call calc_rmsd_simple(seed_u%R(:,1:n_lig_atm), partner_u%R(:,1:n_lig_atm), &
!                      n_lig_atm, dist)
call calc_rmsd_PM(seed_u%R(:,1:n_lig_atm), partner_u%R(:,1:n_lig_atm), &
                  n_type, n_atm_s, atm_id_s, dist)
if (dist > mutation_cut * D_ave) then
    call perturb_trans_rot_gene(new_u%gene, len_tr, len_ang)
end if

end subroutine operator_1_1
!-------------------------------------------------------------------------------
subroutine operator_2(seed_u, partner_u, new_u, n_pert, n_lig_atm)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: seed_u, partner_u
type(ligdock_bank_type), intent(out) :: new_u
integer, intent(inout) :: n_pert
integer, intent(in) :: n_lig_atm
!
!real(dp), parameter :: len_tr = 1.0d0
!real(dp), parameter :: len_ang = 30.0d0 * deg2rad
real(dp), parameter :: len_tr = 0.5d0
real(dp), parameter :: len_ang = 15.0d0 * deg2rad
real(dp) :: dist

new_u = seed_u

! initialize nft
new_u%nft = 0
! trans/rot cross over
new_u%gene(1:6) = partner_u%gene(1:6)

! If seed_u & partner_u aren't similar, perturb trans/rot
!call calc_rmsd_simple(seed_u%R(:,1:n_lig_atm), partner_u%R(:,1:n_lig_atm), &
!                      n_lig_atm, dist)
call calc_rmsd_PM(seed_u%R(:,1:n_lig_atm), partner_u%R(:,1:n_lig_atm), &
                  n_type, n_atm_s, atm_id_s, dist)
if (dist > mutation_cut * D_ave) then
    n_pert = n_pert + 1
    call perturb_trans_rot_gene(new_u%gene, len_tr, len_ang)
end if

end subroutine operator_2
!-------------------------------------------------------------------------------
subroutine operator_2_1(seed_u, partner_u, new_u, protein, ligand)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: seed_u, partner_u
type(ligdock_bank_type), intent(out) :: new_u
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
!
type(ligdock_bank_type) :: tmp_u
type(molecule_type) :: tmp_prot
type(energy_type) :: ff
real(dp) :: g(3,ligand%n_atm)
real(dp), parameter :: len_tr = 0.5d0
real(dp), parameter :: len_ang = 15.0d0 * deg2rad
!real(dp), parameter :: len_tr = 1.0d0
!real(dp), parameter :: len_ang = 30.0d0 * deg2rad
!real(dp), parameter :: len_axis = 0.3d0
real(dp) :: dist
integer :: n_lig_atm

new_u = seed_u
tmp_prot = protein

! initialize nft
new_u%nft = 0
! trans/rot cross over
new_u%gene(1:6) = partner_u%gene(1:6)
!
call construct_ligand(tmp_prot, ligand, new_u%gene)
call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
new_u%E(:) = ff%ligdock(:)

!! try T_crossover
tmp_u = seed_u
tmp_u%gene(1:3) = partner_u%gene(1:3)
call construct_ligand(tmp_prot, ligand, tmp_u%gene)
call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
tmp_u%E(:) = ff%ligdock(:)
if (tmp_u%E(0) < new_u%E(0)) then
    new_u = tmp_u
end if
!
!! try R_crossover
tmp_u = seed_u
tmp_u%gene(4:6) = partner_u%gene(4:6)
call construct_ligand(tmp_prot, ligand, tmp_u%gene)
call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
tmp_u%E(:) = ff%ligdock(:)
if (tmp_u%E(0) < new_u%E(0)) then
    new_u = tmp_u
end if

! If seed_u & partner_u aren't similar, perturb trans/rot
n_lig_atm = ligand%n_atm
!call calc_rmsd_simple(seed_u%R(:,1:n_lig_atm), partner_u%R(:,1:n_lig_atm), &
!                      n_lig_atm, dist)
call calc_rmsd_PM(seed_u%R(:,1:n_lig_atm), partner_u%R(:,1:n_lig_atm), &
                  n_type, n_atm_s, atm_id_s, dist)
if (dist > mutation_cut * D_ave) then
    call perturb_trans_rot_gene(new_u%gene, len_tr, len_ang)
end if

end subroutine operator_2_1
!-------------------------------------------------------------------------------
subroutine operator_3(seed_u, partner_u, new_u, n_tor)
!-------------------------------------------------------------------------------
! Swapping selected torsion angles of ligand
! Partners are selected from 1st bank
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: seed_u, partner_u
type(ligdock_bank_type), intent(out) :: new_u
integer, intent(in) :: n_tor
!
integer :: n_swp_tor, i_tor, tor_idx
logical :: is_used(n_tor)
real(dp), parameter :: pert_tor = 20.0d0 * deg2rad

new_u = seed_u

! Initialize NFT
new_u%nft = 0

if (n_tor == 0) return

n_swp_tor = int(max_tor_select_ratio*n_tor)
if (n_swp_tor < min_tor_select) n_swp_tor = min_tor_select

is_used(:) = .false.
do i_tor = 1, n_swp_tor
    do
        tor_idx = int(n_tor*random()) + 1
        if (.not. is_used(tor_idx)) exit
    end do
    is_used(tor_idx) = .true.

    new_u%gene(tor_idx+6) = partner_u%gene(tor_idx+6) &
                            +(2.0d0*random() - 1.0d0)*pert_tor
    new_u%gene(tor_idx+6) = bound_ang(new_u%gene(tor_idx+6))
end do

end subroutine operator_3
!-------------------------------------------------------------------------------
subroutine operator_3_2(seed_u, partner_u, new_u, n_tor, protein, ligand)
!-------------------------------------------------------------------------------
! Swapping selected torsion angles of ligand
! Partners are selected from 1st bank
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: seed_u, partner_u
type(ligdock_bank_type), intent(out) :: new_u
integer, intent(in) :: n_tor
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type) :: tmp_u
type(molecule_type) :: tmp_prot
type(energy_type) :: ff
real(dp) :: g(3,ligand%n_atm)
!
integer :: n_swp_tor, i_tor, tor_idx
logical :: is_used(n_tor)
real(dp), parameter :: pert_tor = 20.0d0 * deg2rad

new_u = seed_u
tmp_prot = protein

! Initialize NFT
new_u%nft = 0

if (n_tor == 0) return

n_swp_tor = int(max_tor_select_ratio*n_tor)
if (n_swp_tor < min_tor_select) n_swp_tor = min_tor_select

is_used(:) = .false.
do i_tor = 1, n_swp_tor
    do
        tor_idx = int(n_tor*random()) + 1
        if (.not. is_used(tor_idx)) exit
    end do
    is_used(tor_idx) = .true.

    new_u%gene(tor_idx+6) = partner_u%gene(tor_idx+6) &
                            +(2.0d0*random() - 1.0d0)*pert_tor
    new_u%gene(tor_idx+6) = bound_ang(new_u%gene(tor_idx+6))
end do

call construct_ligand(tmp_prot, ligand, new_u%gene)
call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
new_u%E(:) = ff%ligdock(:)

!! try T_crossover
tmp_u = new_u
tmp_u%gene(1:3) = partner_u%gene(1:3)
call construct_ligand(tmp_prot, ligand, tmp_u%gene)
call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
tmp_u%E(:) = ff%ligdock(:)
if (tmp_u%E(0) < new_u%E(0)) then
    new_u = tmp_u
end if
!
!! try R_crossover
tmp_u = new_u
tmp_u%gene(4:6) = partner_u%gene(4:6)
call construct_ligand(tmp_prot, ligand, tmp_u%gene)
call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
tmp_u%E(:) = ff%ligdock(:)
if (tmp_u%E(0) < new_u%E(0)) then
    new_u = tmp_u
end if

end subroutine operator_3_2
!-------------------------------------------------------------------------------
subroutine operator_4(seed_u, partner_u, new_u, n_tor)
!-------------------------------------------------------------------------------
! Swapping selected torsion angles of ligand
! Partners are selected from current bank
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: seed_u, partner_u
type(ligdock_bank_type), intent(out) :: new_u
integer, intent(in) :: n_tor
!
integer :: n_swp_tor, i_tor, tor_idx
logical :: is_used(n_tor)
real(dp), parameter :: pert_tor = 20.0d0 * deg2rad

new_u = seed_u

! Initialize NFT
new_u%nft = 0

if (n_tor == 0) return

n_swp_tor = int(max_tor_select_ratio*n_tor)
if (n_swp_tor < min_tor_select) n_swp_tor = min_tor_select

is_used(:) = .false.
do i_tor = 1, n_swp_tor
    do
        tor_idx = int(n_tor*random()) + 1
        if (.not. is_used(tor_idx)) exit
    end do
    is_used(tor_idx) = .true.

    new_u%gene(tor_idx+6) = partner_u%gene(tor_idx+6) &
                            +(2.0d0*random() - 1.0d0)*pert_tor
    new_u%gene(tor_idx+6) = bound_ang(new_u%gene(tor_idx+6))
end do

end subroutine operator_4
!-------------------------------------------------------------------------------
subroutine operator_4_2(seed_u, partner_u, new_u, n_tor,protein, ligand)
!-------------------------------------------------------------------------------
! Swapping selected torsion angles of ligand
! Partners are selected from current bank
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: seed_u, partner_u
type(ligdock_bank_type), intent(out) :: new_u
integer, intent(in) :: n_tor
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type) :: tmp_u
type(molecule_type) :: tmp_prot
type(energy_type) :: ff
real(dp) :: g(3,ligand%n_atm)
!
integer :: n_swp_tor, i_tor, tor_idx
logical :: is_used(n_tor)
real(dp), parameter :: pert_tor = 20.0d0 * deg2rad

new_u = seed_u
tmp_prot = protein

! Initialize NFT
new_u%nft = 0

if (n_tor == 0) return

n_swp_tor = int(max_tor_select_ratio*n_tor)
if (n_swp_tor < min_tor_select) n_swp_tor = min_tor_select

is_used(:) = .false.
do i_tor = 1, n_swp_tor
    do
        tor_idx = int(n_tor*random()) + 1
        if (.not. is_used(tor_idx)) exit
    end do
    is_used(tor_idx) = .true.

    new_u%gene(tor_idx+6) = partner_u%gene(tor_idx+6) &
                            +(2.0d0*random() - 1.0d0)*pert_tor
    new_u%gene(tor_idx+6) = bound_ang(new_u%gene(tor_idx+6))
end do

call construct_ligand(tmp_prot, ligand, new_u%gene)
call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
new_u%E(:) = ff%ligdock(:)

!! try T_crossover
tmp_u = new_u
tmp_u%gene(1:3) = partner_u%gene(1:3)
call construct_ligand(tmp_prot, ligand, tmp_u%gene)
call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
tmp_u%E(:) = ff%ligdock(:)
if (tmp_u%E(0) < new_u%E(0)) then
    new_u = tmp_u
end if
!
!! try R_crossover
tmp_u = new_u
tmp_u%gene(4:6) = partner_u%gene(4:6)
call construct_ligand(tmp_prot, ligand, tmp_u%gene)
call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
tmp_u%E(:) = ff%ligdock(:)
if (tmp_u%E(0) < new_u%E(0)) then
    new_u = tmp_u
end if

end subroutine operator_4_2
!-------------------------------------------------------------------------------
subroutine operator_3_1(seed_u, partner_u, new_u, protein, ligand, i_cyc)
!-------------------------------------------------------------------------------
! Swapping selected torsion angles of ligand
! Partners are selected from current bank
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: seed_u, partner_u
type(ligdock_bank_type), intent(inout) :: new_u
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: i_cyc
!
type(molecule_type) :: tmp_prot
integer :: n_swp_tor, i_tor, tor_idx, selected, selected_tor, n_tor
integer :: n_tot, n_acc_atm(ligand%n_br)
logical :: is_used(ligand%n_br-1)
real(dp), parameter :: pert_tor = 15.0d0 * deg2rad
integer, parameter :: max_weight = 20
integer, parameter :: reduce_weight = 1
integer :: weight
integer :: clash, n_fail, n_add
integer :: allow_bump

if (ligand%n_br == 1) return

tmp_prot = protein

n_tor = ligand%n_br - 1
n_swp_tor = int(max_tor_select_ratio*n_tor)
if (n_swp_tor < min_tor_select) n_swp_tor = min_tor_select
is_used(:) = .false.

weight = max(max_weight - reduce_weight*(i_cyc-1),1)
n_fail = 0
n_add = 1
allow_bump = 0
do
    new_u = seed_u
    new_u%nft = 0
    do i_tor = 1, n_swp_tor
        n_tot = 0
        n_acc_atm(:) = 0
        do tor_idx = 1, ligand%n_br-1
            n_acc_atm(tor_idx) = n_tot
            if (is_used(tor_idx)) cycle
            n_tot = n_tot + min(ligand%n_atm_br(tor_idx+1), weight)
        end do
        n_acc_atm(ligand%n_br) = n_tot
        !
        selected = int(random()*n_tot) + 1
        do tor_idx = 1, ligand%n_br-1
            if (is_used(tor_idx)) cycle
            if (selected > n_acc_atm(tor_idx) .and. &
                selected <= n_acc_atm(tor_idx+1)) then
                selected_tor = tor_idx
                exit
            end if
        end do
        is_used(selected_tor) = .true.
        !
        new_u%gene(selected_tor+6) = partner_u%gene(selected_tor+6)
        if (random() > 0.5d0) then
            new_u%gene(selected_tor+6) = new_u%gene(selected_tor+6) &
                                 + (2.0d0*random() - 1.0d0)*pert_tor
        end if
    end do
    !
    call construct_ligand(tmp_prot, ligand, new_u%gene)
    call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
    call check_clash_within_ligand(ligand, clash)
    !
    if (clash <= allow_bump) exit
    n_fail = n_fail + 1
    if (n_fail > 10*n_add) then
        n_add = n_add + 1
        allow_bump = allow_bump + 1
    end if
end do

end subroutine operator_3_1
!-------------------------------------------------------------------------------
subroutine operator_4_1(seed_u, partner_u, new_u, protein, ligand, i_cyc)
!-------------------------------------------------------------------------------
! Swapping selected torsion angles of ligand
! Partners are selected from current bank
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: seed_u, partner_u
type(ligdock_bank_type), intent(inout) :: new_u
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: i_cyc
!
type(molecule_type) :: tmp_prot
integer :: n_swp_tor, i_tor, tor_idx, selected, selected_tor, n_tor
integer :: n_tot, n_acc_atm(ligand%n_br)
logical :: is_used(ligand%n_br-1)
real(dp), parameter :: pert_tor = 15.0d0 * deg2rad
integer, parameter :: max_weight = 20
integer, parameter :: reduce_weight = 1
integer :: weight
integer :: clash, n_fail, n_add
integer :: allow_bump

if (ligand%n_br == 1) return

tmp_prot = protein

n_tor = ligand%n_br - 1
n_swp_tor = int(max_tor_select_ratio*n_tor)
if (n_swp_tor < min_tor_select) n_swp_tor = min_tor_select
is_used(:) = .false.

weight = max(max_weight - reduce_weight*(i_cyc-1),1)
n_fail = 0
n_add = 1
allow_bump = 0
do
    new_u = seed_u
    new_u%nft = 0
    do i_tor = 1, n_swp_tor
        n_tot = 0
        n_acc_atm(:) = 0
        do tor_idx = 1, ligand%n_br-1
            n_acc_atm(tor_idx) = n_tot
            if (is_used(tor_idx)) cycle
            n_tot = n_tot + min(ligand%n_atm_br(tor_idx+1), weight)
        end do
        n_acc_atm(ligand%n_br) = n_tot
        !
        selected = int(random()*n_tot) + 1
        do tor_idx = 1, ligand%n_br-1
            if (is_used(tor_idx)) cycle
            if (selected > n_acc_atm(tor_idx) .and. &
                selected <= n_acc_atm(tor_idx+1)) then
                selected_tor = tor_idx
                exit
            end if
        end do
        is_used(selected_tor) = .true.
        !
        new_u%gene(selected_tor+6) = partner_u%gene(selected_tor+6)
        if (random() > 0.5d0) then
            new_u%gene(selected_tor+6) = new_u%gene(selected_tor+6) &
                                 + (2.0d0*random() - 1.0d0)*pert_tor
        end if
    end do
    !
    call construct_ligand(tmp_prot, ligand, new_u%gene)
    call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
    call check_clash_within_ligand(ligand, clash)
    !
    if (clash <= allow_bump) exit
    n_fail = n_fail + 1
    if (n_fail > 10*n_add) then
        n_add = n_add + 1
        allow_bump = allow_bump + 1
    end if
end do


end subroutine operator_4_1
!-------------------------------------------------------------------------------
subroutine operator_5(seed_u, partner_u, new_u, protein, ligand, i_cyc)
!-------------------------------------------------------------------------------
! crossover torsion & re-orient using fftmap
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: seed_u, partner_u
type(ligdock_bank_type), intent(inout) :: new_u
type(ligand_type), intent(in) :: ligand
type(molecule_type), intent(in) :: protein
integer, intent(in) :: i_cyc
!
integer, parameter :: max_gene = 50
type(molecule_type) :: tmp_prot
type(ligdock_bank_type) :: tmp_bank, soln_pool(max_gene)
type(hotspot_graph_type) :: lig_graph
integer :: n_swp_tor, i_tor, tor_idx, selected, selected_tor, n_tor
integer :: n_tot, n_acc_atm(ligand%n_br)
logical :: is_used(ligand%n_br-1)
real(dp), parameter :: pert_tor = 15.0d0 * deg2rad
integer, parameter :: max_weight = 20
integer, parameter :: reduce_weight = 1
integer :: weight
!
integer :: i_gene, n_gene, gene_dim, n_soln, i_soln
real(dp) :: gene_s(ligand%n_br+6, max_gene)
type(energy_type) :: ff
real(dp) :: g(3, ligand%n_atm)
integer :: status, ierr
!
real(dp), parameter :: len_tr = 0.5d0
real(dp), parameter :: len_ang = 15.0d0 * deg2rad
real(dp), parameter :: len_axis = 0.15d0
integer :: clash, n_fail, n_add
integer :: allow_bump

if (ligand%n_br == 1) return

tmp_prot = protein
n_tor = ligand%n_br - 1
!n_swp_tor = int(max_tor_select_ratio*n_tor)
!if (n_swp_tor < min_tor_select) n_swp_tor = min_tor_select
is_used(:) = .false.

weight = max(max_weight - reduce_weight*(i_cyc-1),1)
n_fail = 0
n_add = 1
allow_bump = 0
do
    new_u = seed_u
    new_u%nft = 0
    n_swp_tor = int(random()*n_tor) + 1
    do i_tor = 1, n_swp_tor
        n_tot = 0
        n_acc_atm(:) = 0
        do tor_idx = 1, ligand%n_br-1
            n_acc_atm(tor_idx) = n_tot
            if (is_used(tor_idx)) cycle
            n_tot = n_tot + min(ligand%n_atm_br(tor_idx+1), weight)
        end do
        n_acc_atm(ligand%n_br) = n_tot
        !
        selected = int(random()*n_tot) + 1
        do tor_idx = 1, ligand%n_br-1
            if (is_used(tor_idx)) cycle
            if (selected > n_acc_atm(tor_idx) .and. &
                selected <= n_acc_atm(tor_idx+1)) then
                selected_tor = tor_idx
                exit
            end if
        end do
        is_used(selected_tor) = .true.
        !
        new_u%gene(selected_tor+6) = partner_u%gene(selected_tor+6)
        if (random() > 0.5d0) then
            new_u%gene(selected_tor+6) = new_u%gene(selected_tor+6) &
                                 + (2.0d0*random() - 1.0d0)*pert_tor
        end if
    end do
    !
    call construct_ligand(tmp_prot, ligand, new_u%gene)
    call update_r_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
    call check_clash_within_ligand(ligand, clash)
    !
    if (clash <= allow_bump) exit
    n_fail = n_fail + 1
    if (n_fail > 10*n_add) then
        n_add = n_add + 1
        allow_bump = allow_bump + 1
    end if
end do

!new_u = seed_u
!! Initialize NFT
!new_u%nft = 0
!if (n_tor == 0) return
!n_swp_tor = int(max_tor_select_ratio*n_tor)
!if (n_swp_tor < min_tor_select) n_swp_tor = min_tor_select
!
!is_used(:) = .false.
!do i_tor = 1, n_swp_tor
!    do
!        tor_idx = int(n_tor*random()) + 1
!        if (.not. is_used(tor_idx)) exit
!    end do
!    is_used(tor_idx) = .true.
!    new_u%gene(tor_idx+6) = partner_u%gene(tor_idx+6) &
!                            +(2.0d0*random() - 1.0d0)*pert_tor
!    new_u%gene(tor_idx+6) = bound_ang(new_u%gene(tor_idx+6))
!end do

gene_dim = ligand%n_br + 5
tmp_bank = new_u
call construct_ligand(tmp_prot, ligand, tmp_bank%gene(1:gene_dim))
call generate_lig_graph(tmp_prot%ligand(ligand%lig_no), lig_graph)
do i_gene = 1, max_gene
    gene_s(1:gene_dim, i_gene) = tmp_bank%gene(1:gene_dim)
end do
call do_pharmdock(tmp_prot, ligand, rec_graph, lig_graph,&
                  tmp_bank%gene(1:3), gene_s, n_gene, e1max*5.0d0, max_gene)
do i_gene = 1, n_gene
    tmp_bank%gene(1:gene_dim) = gene_s(1:gene_dim,i_gene)
    soln_pool(i_gene) = tmp_bank
end do

! minimize trial bank made by torsion crossover
call do_simplex(tmp_prot, ligand, new_u%gene, &
                new_u%e, new_u%nft, max_iter_in=10)
!
do i_soln = 1, min(3,n_gene)
    call do_simplex(tmp_prot, ligand, soln_pool(i_soln)%gene, &
                    soln_pool(i_soln)%e, soln_pool(i_soln)%nft, max_iter_in=10)
    if (soln_pool(i_soln)%e(0) < new_u%e(0)) then
        new_u = soln_pool(i_soln)
    end if
end do
!
!if (random() > 0.5) then
!    call perturb_trans_rot_gene(new_u%gene, len_tr, len_ang, len_axis)
!end if

end subroutine operator_5
!-------------------------------------------------------------------------------
subroutine operator_6(seed_u, partner_u, new_u, protein, ligand, i_cyc)
!-------------------------------------------------------------------------------
! crossover torsion & re-orient using fftmap
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: seed_u, partner_u
type(ligdock_bank_type), intent(inout) :: new_u
type(ligand_type), intent(in) :: ligand
type(molecule_type), intent(in) :: protein
integer, intent(in) :: i_cyc
!
integer, parameter :: max_gene = 50
type(molecule_type) :: tmp_prot
type(ligdock_bank_type) :: tmp_bank, soln_pool(max_gene)
type(hotspot_graph_type) :: lig_graph
integer :: n_swp_tor, i_tor, tor_idx, selected, selected_tor, n_tor
integer :: n_tot, n_acc_atm(ligand%n_br)
logical :: is_used(ligand%n_br-1)
real(dp), parameter :: pert_tor = 15.0d0 * deg2rad
integer, parameter :: max_weight = 20
integer, parameter :: reduce_weight = 1
integer :: weight
!
integer :: i_gene, n_gene, gene_dim, n_soln, i_soln
real(dp) :: gene_s(ligand%n_br+6, max_gene)
type(energy_type) :: ff
real(dp) :: g(3, ligand%n_atm)
integer :: status, ierr
!
real(dp), parameter :: len_tr = 0.5d0
real(dp), parameter :: len_ang = 15.0d0 * deg2rad
real(dp), parameter :: len_axis = 0.15d0
integer :: clash, n_fail, n_add
integer :: allow_bump

if (ligand%n_br == 1) return
!
tmp_prot = protein
n_tor = ligand%n_br - 1
n_swp_tor = int(max_tor_select_ratio*n_tor)
if (n_swp_tor < min_tor_select) n_swp_tor = min_tor_select
is_used(:) = .false.

weight = max(max_weight - reduce_weight*(i_cyc-1),1)
n_fail = 0
n_add = 1
allow_bump = 0
do
    new_u = seed_u
    new_u%nft = 0
    n_swp_tor = int(random()*n_tor) + 1
    do i_tor = 1, n_swp_tor
        n_tot = 0
        n_acc_atm(:) = 0
        do tor_idx = 1, ligand%n_br-1
            n_acc_atm(tor_idx) = n_tot
            if (is_used(tor_idx)) cycle
            n_tot = n_tot + min(ligand%n_atm_br(tor_idx+1), weight)
        end do
        n_acc_atm(ligand%n_br) = n_tot
        !
        selected = int(random()*n_tot) + 1
        do tor_idx = 1, ligand%n_br-1
            if (is_used(tor_idx)) cycle
            if (selected > n_acc_atm(tor_idx) .and. &
                selected <= n_acc_atm(tor_idx+1)) then
                selected_tor = tor_idx
                exit
            end if
        end do
        is_used(selected_tor) = .true.
        !
        new_u%gene(selected_tor+6) = partner_u%gene(selected_tor+6)
        if (random() > 0.5d0) then
            new_u%gene(selected_tor+6) = new_u%gene(selected_tor+6) &
                                 + (2.0d0*random() - 1.0d0)*pert_tor
        end if
    end do
    !
    call construct_ligand(tmp_prot, ligand, new_u%gene)
    call update_r_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
    call check_clash_within_ligand(ligand, clash)
    !
    if (clash <= allow_bump) exit
    n_fail = n_fail + 1
    if (n_fail > 10*n_add) then
        n_add = n_add + 1
        allow_bump = allow_bump + 1
    end if
end do

!new_u = seed_u
!! Initialize NFT
!new_u%nft = 0
!if (n_tor == 0) return
!n_swp_tor = int(max_tor_select_ratio*n_tor)
!if (n_swp_tor < min_tor_select) n_swp_tor = min_tor_select
!
!is_used(:) = .false.
!do i_tor = 1, n_swp_tor
!    do
!        tor_idx = int(n_tor*random()) + 1
!        if (.not. is_used(tor_idx)) exit
!    end do
!    is_used(tor_idx) = .true.
!    new_u%gene(tor_idx+6) = partner_u%gene(tor_idx+6) &
!                            +(2.0d0*random() - 1.0d0)*pert_tor
!    new_u%gene(tor_idx+6) = bound_ang(new_u%gene(tor_idx+6))
!end do


gene_dim = ligand%n_br + 5
tmp_bank = new_u
call construct_ligand(tmp_prot, ligand, tmp_bank%gene(1:gene_dim))
call generate_lig_graph(tmp_prot%ligand(ligand%lig_no), lig_graph)
do i_gene = 1, max_gene
    gene_s(1:gene_dim, i_gene) = tmp_bank%gene(1:gene_dim)
end do
call do_pharmdock(tmp_prot, ligand, rec_graph, lig_graph,&
                  tmp_bank%gene(1:3), gene_s, n_gene, e1max*5.0d0, max_gene)
do i_gene = 1, n_gene
    tmp_bank%gene(1:gene_dim) = gene_s(1:gene_dim,i_gene)
    soln_pool(i_gene) = tmp_bank
end do

! minimize trial bank made by torsion crossover
call do_simplex(tmp_prot, ligand, new_u%gene, &
                new_u%e, new_u%nft, max_iter_in=10)

do i_soln = 1, min(3,n_gene)
    call do_simplex(tmp_prot, ligand, soln_pool(i_soln)%gene, &
                    soln_pool(i_soln)%e, soln_pool(i_soln)%nft, max_iter_in=10)
    if (soln_pool(i_soln)%e(0) < new_u%e(0)) then
        new_u = soln_pool(i_soln)
    end if
end do

end subroutine operator_6
!-------------------------------------------------------------------------------
subroutine operator_7(seed_u, partner_u, new_u, protein)
!-------------------------------------------------------------------------------
! Swapping chi angles of flexible sc.
! partner comes from first bank
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: seed_u, partner_u
type(ligdock_bank_type), intent(out) :: new_u
type(molecule_type), intent(in) :: protein
!
type(molecule_type) :: tmp_prot
integer :: i_usc, i_res, i_rot_start, i_rot_end, n_rot, rot_idx
real(dp) :: swap_check, mut_check

new_u = seed_u
new_u%nft = 0

do i_usc = 1, n_usc
    swap_check = random()
    if (swap_check > 0.5d0) then
        new_u%rot_idx(i_usc) = partner_u%rot_idx(i_usc)
        new_u%flex_sc(i_usc) = partner_u%flex_sc(i_usc)
    end if
end do

tmp_prot = protein
do i_usc = 1, n_usc
    i_res = usc_list(i_usc)
    mut_check = random()
    if (mut_check > 0.5d0) then
        call get_rotamer_index(tmp_prot%residue(i_res:i_res+1), i_rot_start, &
                               i_rot_end)
        n_rot = i_rot_end - i_rot_start + 1
        rot_idx = i_rot_start + int(random()*n_rot)
        call place_rotamer(tmp_prot, i_res, rot_idx)
        new_u%rot_idx(i_usc) = rot_idx - i_rot_start + 1
        new_u%flex_sc(i_usc)%R = tmp_prot%residue(i_res)%R
    end if
end do

end subroutine operator_7
!-------------------------------------------------------------------------------
subroutine operator_8(seed_u, partner_u, new_u, protein)
!-------------------------------------------------------------------------------
! Swapping chi angles of flexible sc.
! partner comes from current bank
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: seed_u, partner_u
type(ligdock_bank_type), intent(out) :: new_u
type(molecule_type), intent(in) :: protein
!
type(molecule_type) :: tmp_prot
integer :: i_usc, i_res, i_rot_start, i_rot_end, n_rot, rot_idx
real(dp) :: swap_check, mut_check

new_u = seed_u
new_u%nft = 0

do i_usc = 1, n_usc
    swap_check = random()
    if (swap_check > 0.5d0) then
        new_u%rot_idx(i_usc) = partner_u%rot_idx(i_usc)
        new_u%flex_sc(i_usc) = partner_u%flex_sc(i_usc)
    end if
end do

tmp_prot = protein
do i_usc = 1, n_usc
    i_res = usc_list(i_usc)
    mut_check = random()
    if (mut_check > 0.5d0) then
        call get_rotamer_index(tmp_prot%residue(i_res:i_res+1), i_rot_start, &
                               i_rot_end)
        n_rot = i_rot_end - i_rot_start + 1
        rot_idx = i_rot_start + int(random()*n_rot)
        call place_rotamer(tmp_prot, i_res, rot_idx)
        new_u%rot_idx(i_usc) = rot_idx - i_rot_start + 1
        new_u%flex_sc(i_usc)%R = tmp_prot%residue(i_res)%R
    end if
end do

end subroutine operator_8
!-------------------------------------------------------------------------------
subroutine perturb_trans_rot_gene(gene, len_trans, len_ang)
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: gene(:)
real(dp), intent(in) :: len_trans, len_ang
!real(dp) :: mut_check, q0(4), q1(4)
real(dp) :: mut_check, p(3), q(4)
integer :: i

do i = 1, 3
   mut_check = random()
   if (mut_check > 0.5d0) then
      gene(i) = gene(i) + (2.0d0*random() - 1.0d0)*len_trans
   end if
end do

p = gene(4:6)
call expmap_quat(p,q)
call pert_rotation(len_ang, q)
call inv_expmap_quat(q,p)
gene(4:6) = p

end subroutine perturb_trans_rot_gene
!-------------------------------------------------------------------------------
subroutine optimize_new_conf_rigid(protein, ligand, bank, n_bank, use_grad)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type), intent(inout) :: bank(:)
logical, optional, intent(in) :: use_grad
integer, intent(in) :: n_bank
type(molecule_type) :: tmp_prot
integer :: i_bank
logical :: status

tmp_prot = protein
do i_bank = 1, n_bank
    write(log_msg, '(A, I4)') 'Minimizing new conformation: ', i_bank
    call log_p(log_msg, level=40)
   
    if (present(use_grad) .and. use_grad) then
        call do_simplex(tmp_prot, ligand, bank(i_bank)%gene, bank(i_bank)%E,&
                        bank(i_bank)%nft, max_iter_in=10, pert_full=.true.)
        call opt_ligand_steepest_descent_intcrd(bank(i_bank)%E, status, bank(i_bank)%gene, &
                                                tmp_prot, ligand, bank(i_bank)%nft)
    else
        call do_simplex(tmp_prot, ligand, bank(i_bank)%gene, bank(i_bank)%E,&
                        bank(i_bank)%nft)
    end if
    call ligand2bank(tmp_prot%ligand(ligand%lig_no), bank(i_bank))
end do

n_gen_str = n_gen_str + n_bank
do i_bank = 1, n_bank
    total_nft = total_nft + bank(i_bank)%nft
end do
call update_R_for_ligand(protein%ligand(ligand%lig_no), ligand%lig_no)

end subroutine optimize_new_conf_rigid
!-------------------------------------------------------------------------------
subroutine optimize_seed_rigid(protein, ligand, bank)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type), intent(inout) :: bank(:)
type(molecule_type) :: tmp_prot
integer :: i_bank, i_seed

tmp_prot = protein
do i_seed = 1, n_seed
    i_bank = idx_seed(i_seed)
    write(log_msg, '(A, I4)') 'Re-minimizing seed: ', i_bank
    call log_p(log_msg, level=40)
    
    call do_simplex(tmp_prot, ligand, bank(i_bank)%gene, bank(i_bank)%E,&
                    bank(i_bank)%nft)
    call ligand2bank(tmp_prot%ligand(ligand%lig_no), bank(i_bank))
end do

call update_R_for_ligand(protein%ligand(ligand%lig_no), ligand%lig_no)

end subroutine optimize_seed_rigid
!-------------------------------------------------------------------------------
subroutine optimize_new_conf_flex_sc(protein, ligand, bank, n_bank)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type), intent(inout) :: bank(:)
integer, intent(in) :: n_bank
type(molecule_type) :: tmp_prot
integer :: i_bank, i_res, i_usc

tmp_prot = protein
do i_bank = 1, n_bank
    write(log_msg, '(A, I4)') 'Minimizing new conformation: ', i_bank
    call log_p(log_msg, level=40)
   
    do i_usc = 1, n_usc
        i_res = usc_list(i_usc)
        tmp_prot%residue(i_res)%R = bank(i_bank)%flex_sc(i_usc)%R
        call update_R_for_sc(tmp_prot%residue(i_res), i_res)
    end do
    call set_atdk_hbond_flex(usc_list, n_usc)
    call do_simplex(tmp_prot, ligand, bank(i_bank)%gene, bank(i_bank)%E,&
                    bank(i_bank)%nft)
    call update_prot_prot_intrxn(bank(i_bank)%rot_idx(1:n_usc), bank(i_bank)%E)
    call ligand2bank(tmp_prot%ligand(ligand%lig_no), bank(i_bank))
end do

n_gen_str = n_gen_str + n_bank
do i_bank = 1, n_bank
    total_nft = total_nft + bank(i_bank)%nft
end do
call protein_to_R(protein)

end subroutine optimize_new_conf_flex_sc
!-------------------------------------------------------------------------------
subroutine optimize_seed_flex_sc(protein, ligand, bank)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type), intent(inout) :: bank(:)
type(molecule_type) :: tmp_prot
integer :: i_bank, i_seed, i_res, i_usc

tmp_prot = protein
do i_seed = 1, n_seed
    i_bank = idx_seed(i_seed)
    write(log_msg, '(A, I4)') 'Re-minimizing seed: ', i_bank
    call log_p(log_msg, level=40)
        
    do i_usc = 1, n_usc
        i_res = usc_list(i_usc)
        tmp_prot%residue(i_res)%R = bank(i_bank)%flex_sc(i_usc)%R
        call update_R_for_sc(tmp_prot%residue(i_res), i_res)
    end do
    call set_atdk_hbond_flex(usc_list, n_usc)
    call do_simplex(tmp_prot, ligand, bank(i_bank)%gene, bank(i_bank)%E,&
                    bank(i_bank)%nft)
    call update_prot_prot_intrxn(bank(i_bank)%rot_idx(1:n_usc), bank(i_bank)%E)
    call ligand2bank(tmp_prot%ligand(ligand%lig_no), bank(i_bank))
end do

call protein_to_R(protein)

end subroutine optimize_seed_flex_sc
!-------------------------------------------------------------------------------
subroutine optimize_bank_steepest(protein, ligand, bank, n_bank)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: n_bank
type(ligdock_bank_type), intent(inout) :: bank(:)
type(molecule_type) :: tmp_prot
integer :: i_bank
logical :: status
type(energy_type) :: ff
real(dp) :: g(3,ligand%n_atm)

!$OMP PARALLEL DO default(shared) private(status, tmp_prot, ff, g) COPYIN(R) schedule(dynamic)
tmp_prot = protein
do i_bank = 1, n_bank
    write(log_msg, '(A, I4)') 'Re-minimizing bank: ', i_bank
    call log_p(log_msg, level=40)
    call opt_ligand_steepest_descent_intcrd(bank(i_bank)%E, status, bank(i_bank)%gene, &
                                            tmp_prot,ligand, bank(i_bank)%nft)
    call ligand2bank(tmp_prot%ligand(ligand%lig_no), bank(i_bank))
end do
!$OMP END PARALLEL DO

call update_R_for_ligand(protein%ligand(ligand%lig_no), ligand%lig_no)
do i_bank = 1, n_bank
    total_nft = total_nft + bank(i_bank)%nft
end do

end subroutine optimize_bank_steepest
!-------------------------------------------------------------------------------
END MODULE LIGDOCK_CSA_NEW_CONF
!-------------------------------------------------------------------------------
