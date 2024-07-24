!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/apps/ligand_docking/ligand_docking.f90
!
! Description: 
!
!-------------------------------------------------------------------------------
MODULE LIGAND_DOCKING
!-------------------------------------------------------------------------------
use globals
use logger, only: log_p, my_timer
use in_out_ligand, only: write_mol2
!
use ligdock_csa_vars
use ligdock_csa_runner, only: initialize_ligdock_csa, finalize_ligdock_csa, &
                              run_ligdock_csa_cycle

implicit none
save
private

public :: ligdock_runner

CONTAINS
!-------------------------------------------------------------------------------
subroutine ligdock_runner(protein, ligand, include_curr_conf)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(inout) :: ligand
logical, intent(in) :: include_curr_conf
!
integer :: gene_dim, i_csa
real(dp) :: csa_start, opt_E
!
character(len=len_fname) :: mol2_file

mol2_file = 'merged_ligand.mol2'
call write_mol2(mol2_file, ref_lig(ligand%lig_type)%R, ref_lig(ligand%lig_type))

write(log_msg, '(A,I5)') '  Running ligand docking...'
call log_p(log_msg, level=20)

call my_timer(csa_start)
gene_dim = 6 + (ligand%n_br - 1)
call initialize_ligdock_csa(gene_dim)

do i_csa = 1, n_csa_cycle
    call run_ligdock_csa_cycle(protein, ligand, opt_E, include_curr_conf, &
                               i_csa, csa_start)
end do

call finalize_ligdock_csa()

end subroutine ligdock_runner
!-------------------------------------------------------------------------------
END MODULE LIGAND_DOCKING
!-------------------------------------------------------------------------------
