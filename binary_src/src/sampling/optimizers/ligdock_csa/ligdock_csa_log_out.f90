!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/sampling/optimizers/ligdock_csa/ligdock_csa_log_out.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE LIGDOCK_CSA_LOG_OUT
!-------------------------------------------------------------------------------

use globals
!
use logger
use in_out_structure, only: open_write_pdb, close_write_pdb, write_pdb_model
use in_out_ligand, only: write_multiple_mol2
!
use sort, only: sort2
use rmsd, only: prep_RMSD_PM_calculation, calc_RMSD_PM
use rotamer, only: place_rotamer
!
use xlogp_Score_m, only: get_xlogP_atom_types
!
use ligdock_csa_vars

implicit none
save
public

CONTAINS
!===============================================================================
! Print CSA log
!===============================================================================
subroutine print_ligdock_csa_log_header()
!-------------------------------------------------------------------------------
write(log_msg,"(A,A)") "LCSA PLD: ",&
   "icsa icy  nstep   cutdif  imn imx    ebmin     ebmax        nft   ius  nbk  time"
call log_p(log_msg, level=20)

end subroutine print_ligdock_csa_log_header
!-------------------------------------------------------------------------------
subroutine print_bank_info(bank, n_bank)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: bank(:)
integer, intent(in) :: n_bank
integer :: i_bank

call log_p('    Bank No.:   Energy', level=30)
do i_bank = 1, n_bank
    write(log_msg, '(6X,A,1X,I3,A,F10.3)') 'Bank', i_bank, ': ', bank(i_bank)%E(0)
    call log_p(log_msg, level = 30)
end do

end subroutine print_bank_info
!-------------------------------------------------------------------------------
subroutine print_seed_cycle_log(i_csa, i_iter, i_seed)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_csa, i_iter, i_seed
call log_divider(level = 20)
write(log_msg,"(4X,A,I4,1X,1A,1X,I4,1X,1A,1X,I4)") &
               "CSA Run Number / Interation / SEED CYCLE = ", &
                i_csa, '/', i_iter, '/', i_seed
call log_p(log_msg, level=30)
write(log_msg,"(6X,A,F8.3,1X,1A,1X,F8.3,1X,1A,1X,F8.3)") &
              "D_ave / D_cut / D_min = ", D_ave,'/', D_cut,'/', D_min
call log_p(log_msg, level=20)

end subroutine print_seed_cycle_log
!-------------------------------------------------------------------------------
subroutine print_ligdock_csa_log(bank, n_bank, i_csa, i_seed_cycle, &
                                 i_Emin_bank_u, i_Emax_bank_u, time_start)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: bank(:)
integer, intent(in) :: n_bank
integer, intent(in) :: i_csa, i_seed_cycle
integer, intent(in) :: i_Emin_bank_u, i_Emax_bank_u
real(dp), intent(in) :: time_start
!
character(len=len_fname) :: time_message
integer :: n_unused_seeds
integer :: i
real(dp) :: time_cycle, time_diff

time_message = 'Used time for CSA cycle'

n_unused_seeds=0
do i = 1, n_bank
   if (.not. bank(i)%used) n_unused_seeds = n_unused_seeds + 1
end do

call my_timer(time_cycle)
call report_time_diff(time_message, time_start, time_cycle)

time_diff = time_cycle - time_start

330 format (A,1X,I3,1X,I3,1X,I6,1X,F8.3,1X,I3,1X,I3,1X,2F10.3,1X,I10,1X,I4,1X,I4,1X,F10.3,A)

write(log_msg,330) &
     "LCSA PLD:", i_csa, i_seed_cycle, n_gen_str, D_cut, i_Emin_bank_u, i_Emax_bank_u, &
     bank(i_Emin_bank_u)%E(0), bank(i_Emax_bank_u)%E(0), &
     total_nft, n_unused_seeds, n_bank, time_diff,' sec'

call log_p(log_msg, level=20)
  
end subroutine print_ligdock_csa_log
!===============================================================================
! write CSA outputs
!===============================================================================
subroutine write_bank(bank_prefix, bank, n_bank, protein0, ligand)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: bank_prefix
type(ligdock_bank_type), intent(in) :: bank(:)
integer, intent(in) :: n_bank
type(molecule_type), intent(in) :: protein0
type(ligand_type), intent(in) :: ligand
!
integer :: i_bank, i_res, i_atm, i_usc
real(dp) :: lig_crd_s(3,ligand%n_atm,n_bank)
type(molecule_type) :: protein
character(len=len_fname) :: out_file
integer, parameter :: pdb_unit=74

do i_bank = 1, n_bank
    do i_atm = 1, ligand%n_atm
        lig_crd_s(:,i_atm,i_bank) = bank(i_bank)%R(:,i_atm)
    end do
end do

out_file = trim(bank_prefix) // '.mol2'
call write_multiple_mol2(out_file, lig_crd_s, ref_lig(ligand%lig_type), n_bank)

out_file = trim(bank_prefix) // '.E.info'
call write_bank_info_file(out_file, bank, n_bank, ligand)

if (n_usc > 0) then
    out_file = trim(bank_prefix) // '.pdb'
    call open_write_pdb(pdb_unit, out_file)
    !
    do i_bank = 1, n_bank
        i_usc = 0
        protein = protein0
        do i_usc = 1, n_usc
            i_res = usc_list(i_usc)
            protein%residue(i_res)%R = bank(i_bank)%flex_sc(i_usc)%R
        end do
        do i_atm = 1, ligand%n_atm
            protein%ligand(ligand%lig_no)%R(:,i_atm) = bank(i_bank)%R(:,i_atm)
        end do
        !
        call write_pdb_model(pdb_unit, protein, i_bank)
    end do
    !
    call close_write_pdb(pdb_unit)
end if

end subroutine write_bank
!-------------------------------------------------------------------------------
subroutine write_bank_info_file(info_file, bank, n_bank, ligand)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: info_file
type(ligdock_bank_type), intent(in) :: bank(:)
integer, intent(in) :: n_bank
type(ligand_type), intent(in) :: ligand
!
character(len=20) :: atom_types(ligand%n_atm)
integer :: n_type, n_atm_s(ligand%n_atm), atm_id_s(ligand%n_atm, ligand%n_atm)
integer :: i_bank
real(dp) :: R_mat(n_bank), rmsd
!
integer :: f_unit, ioerror

call get_xlogP_atom_types(ref_lig(ligand%lig_type), atom_types)
call prep_RMSD_PM_calculation(atom_types, ligand%n_atm, n_type, n_atm_s, atm_id_s)

do i_bank = 1, n_bank
    call calc_RMSD_PM(bank(i_bank)%R, ref_lig(ligand%lig_type)%R, n_type, &
                      n_atm_s, atm_id_s, rmsd)
    R_mat(i_bank) = rmsd
end do

f_unit = 17
open(unit=f_unit, file=trim(info_file), action='write', status='replace', &
     iostat=ioerror)
if(ioerror /= 0) stop 'cannot create file'

340 format(A7,4X,A6,4X,A6,4X,A6,5X,A5,6X,A4,6X,A4,7X,A3,6X,A4)
350 format(2X,I3,2X,9(2X,F8.3))
write(f_unit, '(a)') '!----------------------------------------------------------------------------'
write(f_unit, 340) 'Bank No', 'Energy', 'l_RMSD', 'ATDK_E', 'INT_E', &
                   'DS_E', 'HM_E', 'PLP', 'PROT'
write(f_unit, '(a)') '!----------------------------------------------------------------------------'
do i_bank = 1, n_bank
   write(f_unit, 350) i_bank, bank(i_bank)%E_ML, R_mat(i_bank), bank(i_bank)%E(1:5), &
                      bank(i_bank)%E(6:7)
end do
write(f_unit, '(a)') '!----------------------------------------------------------------------------'
close(f_unit)

end subroutine write_bank_info_file
!-------------------------------------------------------------------------------
subroutine write_clust_info(prefix, cl_id, cl_size, n_cl)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: prefix
integer, intent(in) :: cl_id(:,:), cl_size(:), n_cl
character(len=len_fname) :: cl_file
integer :: i_cl, i_mem, i_unit, ioerror

cl_file = trim(prefix) // '.size.info'
i_unit = 33

open(unit=i_unit, file=trim(cl_file), action='write', status='replace', &
     iostat=ioerror)
if(ioerror /= 0) call terminate_with_error('ERROR: cannot create file')

do i_cl = 1, n_cl
    write(i_unit, '(A,I4)') 'Cluster No: ', i_cl
    write(i_unit, '(A,I4)') ' - Cluster Size: ', cl_size(i_cl)
    do i_mem = 1, cl_size(i_cl)
        write(i_unit, '(A,I4)') ' - Cluster member: ', cl_id(i_mem,i_cl)
    end do
end do
close(i_unit)

end subroutine write_clust_info
!-------------------------------------------------------------------------------
subroutine calc_colony_E(bank, n_bank, ligand, E_col)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: bank(:)
integer, intent(in) :: n_bank
type(ligand_type), intent(in) :: ligand
real(dp), intent(out) :: E_col(:)
!
integer :: i_bank, j_bank
!
character(len=20) :: atom_types(ligand%n_atm)
integer :: n_type, n_atm_s(ligand%n_atm), atm_id_s(ligand%n_atm, ligand%n_atm)
real(dp) :: R_mat(n_bank, n_bank), rmsd
!
real(dp) :: dist, Z, alpha, beta

call get_xlogP_atom_types(ref_lig(ligand%lig_type), atom_types)
call prep_RMSD_PM_calculation(atom_types, ligand%n_atm, n_type, n_atm_s, atm_id_s)

R_mat(:,:) = 0.0d0
do i_bank = 1, n_bank - 1
    do j_bank = i_bank + 1, n_bank
        call calc_RMSD_PM(bank(i_bank)%R, bank(j_bank)%R, n_type, &
                          n_atm_s, atm_id_s, rmsd)
        R_mat(i_bank,j_bank) = rmsd
        R_mat(j_bank,i_bank) = rmsd
    end do
end do

E_col(:) = 0.0d0
do i_bank = 1, n_bank
    Z = 0.0d0
    do j_bank = 1, n_bank
        dist = R_mat(j_bank, i_bank)
        alpha = exp(-(dist**3))
        beta = exp(-bank(j_bank)%E(0)/RT)
        Z = Z + alpha*beta
    end do
    E_col(i_bank) = -RT*log(Z)
end do

end subroutine calc_colony_E
!-------------------------------------------------------------------------------
subroutine write_colony_info(prefix, bank, n_bank, ligand)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: prefix
type(ligdock_bank_type), intent(in) :: bank(:)
integer, intent(in) :: n_bank
type(ligand_type), intent(in) :: ligand
!
real(dp) :: E_col(n_bank)
integer :: i_bank
integer :: key(n_bank)
character(len=len_fname) :: file_name
integer :: i_unit, ioerror

call calc_colony_E(bank, n_bank, ligand, E_col)
call sort2(n_bank, E_col(1:n_bank), key(1:n_bank))

i_unit = 19
file_name = trim(prefix) // '.info'
open(unit=i_unit, file=trim(file_name), status='replace', action='write', &
     iostat=ioerror)
if (ioerror /= 0) then
    write(log_msg,'(A,A)') 'Cannot make file: ', trim(file_name)
    call terminate_with_error(log_msg)
end if

360 format(A4,2X,I3,2X,A4,2X,I3,2X,A8,2X,F8.3)
do i_bank = 1, n_bank
    write(i_unit,360) 'Rank', i_bank, 'Bank', key(i_bank), 'Colony_E', &
                      E_col(i_bank)
end do
close(i_unit)

end subroutine write_colony_info
!-------------------------------------------------------------------------------
END MODULE LIGDOCK_CSA_LOG_OUT
!-------------------------------------------------------------------------------
