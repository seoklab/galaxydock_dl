!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/sampling/optimizers/ligdock_csa/ligdock_csa_utils.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE LIGDOCK_CSA_UTILS
!-------------------------------------------------------------------------------

use globals
use ran
!
use energy_utils, only: update_R_for_ligand
use ligdock_energy, only: check_clash_within_ligand
use ligand_operator, only: construct_ligand
!
use ligdock_csa_vars

implicit none
save
public

CONTAINS
!-------------------------------------------------------------------------------
subroutine sort_csa_bank(bank, n_bank)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(inout) :: bank(:)
integer, intent(in) :: n_bank
integer :: i_bank, j_bank
integer :: min_loc
real :: min_E
real :: first_E, tmp_E
type(ligdock_bank_type) :: start_bank, tmp_bank

do i_bank = 1, n_bank
    start_bank = bank(i_bank)
    min_E = bank(i_bank)%E(0)
    min_loc = i_bank
    do j_bank = i_bank, n_bank
        if(bank(j_bank)%E(0) < min_E) then
            min_loc = j_bank
            min_E = bank(j_bank)%E(0)
        end if
    end do
    tmp_bank = bank(min_loc)
    bank(min_loc) = start_bank
    bank(i_bank) = tmp_bank
end do

end subroutine sort_csa_bank
!-------------------------------------------------------------------------------
subroutine bank2ligand(bank, ligand_residue)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: bank
type(ligand_residue_type), intent(inout) :: ligand_residue
integer :: i_atm

do i_atm = 1, ligand_residue%n_atm
    ligand_residue%R(:,i_atm) = bank%R(:,i_atm)
end do

end subroutine bank2ligand
!-------------------------------------------------------------------------------
subroutine ligand2bank(ligand_residue, bank)
!-------------------------------------------------------------------------------
type(ligand_residue_type), intent(in) :: ligand_residue
type(ligdock_bank_type), intent(inout) :: bank
integer :: i_atm

do i_atm = 1, ligand_residue%n_atm
    bank%R(:,i_atm) = ligand_residue%R(:,i_atm)
end do

end subroutine ligand2bank
!-------------------------------------------------------------------------------
subroutine allocate_ligdock_bank_type(bank, protein, ligand)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(inout) :: bank
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
integer :: i_res, i_usc, n_atm

if (allocated(bank%R)) return
!
allocate(bank%R(3, ligand%n_atm))

if (n_usc == 0) return
!
allocate(bank%rot_idx(n_usc))
allocate(bank%flex_sc(n_usc))
!
do i_usc = 1, n_usc
    i_res = usc_list(i_usc)
    n_atm = protein%residue(i_res)%n_atm
    !
    bank%flex_sc(i_usc)%n_atm = n_atm
    allocate(bank%flex_sc(i_usc)%R(3, 0:n_atm))
end do

end subroutine allocate_ligdock_bank_type
!-------------------------------------------------------------------------------
subroutine deallocate_ligdock_bank_type(bank)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(inout) :: bank

deallocate(bank%R)
!
if (allocated(bank%rot_idx)) then
    deallocate(bank%rot_idx)
    deallocate(bank%flex_sc)
end if

end subroutine deallocate_ligdock_bank_type
!-------------------------------------------------------------------------------
subroutine sample_torsion_wo_clash(protein, ligand, conf_bank, gene_dim)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(in) :: ligand
type(ligdock_bank_type), intent(inout) :: conf_bank
integer, intent(in) :: gene_dim
!
integer :: n_fail
integer :: clash, tmp_allow_bump
integer, parameter :: allow_bump = 0
!integer, parameter :: allow_bump = 1
!integer, parameter :: allow_bump = 1000

tmp_allow_bump = allow_bump
n_fail = 0
do
    call sample_torsion_random(ligand, gene_dim, &
                               conf_bank%gene(1:gene_dim))
    call construct_ligand(protein, ligand, conf_bank%gene(1:gene_dim))
    call update_R_for_ligand(protein%ligand(ligand%lig_no), ligand%lig_no)
    call check_clash_within_ligand(ligand, clash)
    if (clash <= tmp_allow_bump) exit
    if (n_fail > 1000) then
        n_fail = 0
        tmp_allow_bump = tmp_allow_bump + 1
        cycle
    endif
    n_fail = n_fail + 1
end do

end subroutine sample_torsion_wo_clash
!-------------------------------------------------------------------------------
subroutine sample_torsion_random(ligand, gene_dim, gene)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: gene_dim
real(dp), intent(out) :: gene(:)
!
integer :: i, cntr_idx

!cntr_idx = ligand%atm_in_br(1,1) 
cntr_idx = ref_lig(ligand%lig_type)%cntr_atm
do i = 1, 3
    gene(i) = dock_grid_info%grid_cntr(i) - ref_lig(ligand%lig_type)%R(i, cntr_idx)
end do
gene(4:6) = (/0.0d0, 0.0d0, 0.0d0/)

! 3. generate gene for torsions
do i = 7, gene_dim
    gene(i) = two_pi*random() - pi
end do

end subroutine sample_torsion_random
!-------------------------------------------------------------------------------
END MODULE LIGDOCK_CSA_UTILS
!-------------------------------------------------------------------------------
