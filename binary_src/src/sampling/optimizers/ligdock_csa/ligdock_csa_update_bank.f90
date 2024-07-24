!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/sampling/optimizers/ligdock_csa/ligdock_csa_utils.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE LIGDOCK_CSA_UPDATE_BANK
!-------------------------------------------------------------------------------

use globals
use logger
!
use in_out_vars, only: infile_pdb, infile_pre_ML, infile_ML
use ligdock_csa_vars
use ligdock_csa_dist, only: calc_dist_btw_bank
use ligdock_csa_log_out, only: write_bank

implicit none
save
private

public :: sort_csa_bank
public :: update_bank
public :: update_min_max_bank_unit
public :: calc_ML_energy_ligand

CONTAINS
!===============================================================================
! Subroutines related to update bank
!===============================================================================
subroutine sort_csa_bank(bank, n_bank)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(inout) :: bank(:)
integer, intent(in) :: n_bank
integer :: i_bank, j_bank
integer :: min_loc
real :: min_E
type(ligdock_bank_type) :: start_bank, tmp_bank

do i_bank = 1, n_bank
    start_bank = bank(i_bank)
    min_E = bank(i_bank)%E_ML
    min_loc = i_bank
    do j_bank = i_bank, n_bank
        if(bank(j_bank)%E_ML < min_E) then
            min_loc = j_bank
            min_E = bank(j_bank)%E_ML
        end if
    end do
    tmp_bank = bank(min_loc)
    bank(min_loc) = start_bank
    bank(i_bank) = tmp_bank
end do

end subroutine sort_csa_bank
!-------------------------------------------------------------------------------
subroutine update_bank(bank, bank_new, n_bank, n_new_bank, n_lig_atm, &
                       i_Emax_bank_u, i_Emin_bank_u)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(inout) :: bank(:)
type(ligdock_bank_type), intent(in) :: bank_new(:)
integer, intent(in) :: n_bank, n_new_bank, n_lig_atm
integer, intent(inout) :: i_Emax_bank_u, i_Emin_bank_u
integer :: i_bank, j_bank, i_bank_min
real(dp) :: min_dist, dist

do i_bank = 1, n_new_bank
    if (bank_new(i_bank)%E_ML > bank(i_Emax_bank_u)%E_ML) cycle
    ! Find the closest conformation in bank
    min_dist = 99999.9d0
    do j_bank = 1, n_bank
        call calc_dist_btw_bank(bank_new(i_bank), bank(j_bank), n_lig_atm, dist)
        if (dist < min_dist) then
            min_dist = dist
            i_bank_min = j_bank
        end if
    end do
    call check_and_update_bank(min_dist, bank, bank_new, i_bank_min, i_bank,&
                               i_Emax_bank_u)
    call update_min_max_bank_unit(bank, n_bank, i_Emax_bank_u, i_Emin_bank_u)
end do

end subroutine update_bank
!-------------------------------------------------------------------------------
subroutine check_and_update_bank(min_dist, bank, bank_new, i_bank_min, i_new, &
                                 i_Emax_bank_u)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: min_dist
type(ligdock_bank_type), intent(inout) :: bank(:)
type(ligdock_bank_type), intent(in) :: bank_new(:)
integer, intent(in) :: i_bank_min, i_new, i_Emax_bank_u
logical :: is_used

320 format (6X,I4,1X,F10.3,A,I4,1X,F10.3,A)

if (min_dist < D_cut) then
    if (bank_new(i_new)%E_ML < bank(i_bank_min)%E_ML) then
        if ((bank(i_bank_min)%E_ML - bank_new(i_new)%E_ML) >= rep_conf_seed_cut) then
            write(log_msg,320) i_bank_min, bank(i_bank_min)%E_ML, " was replaced to ", &
                              i_new, bank_new(i_new)%E_ML, " in same group"
            call log_p(log_msg, level=30)
            bank(i_bank_min) = bank_new(i_new)
            bank(i_bank_min)%used = .false.
        else 
            write(log_msg,320) i_bank_min, bank(i_bank_min)%E_ML, " was replaced to ", &
                              i_new, bank_new(i_new)%E_ML, " in same group, no reset used"
            call log_p(log_msg, level=30)
            is_used = bank(i_bank_min)%used
            bank(i_bank_min) = bank_new(i_new)
            bank(i_bank_min)%used = is_used
        end if
    end if
else
    write(log_msg,320) i_Emax_bank_u, bank(i_Emax_bank_u)%E_ML, " was replaced to ", &
                      i_new, bank_new(i_new)%E_ML, " in new group"
    call log_p(log_msg, level=30)
    bank(i_Emax_bank_u) = bank_new(i_new)
    bank(i_Emax_bank_u)%used = .false.
end if

end subroutine check_and_update_bank
!-------------------------------------------------------------------------------
subroutine update_min_max_bank_unit(bank, n_bank, i_Emax_bank_u, i_Emin_bank_u)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: bank(:)
integer, intent(in) :: n_bank
integer, intent(out) :: i_Emax_bank_u, i_Emin_bank_u
integer :: i_bank

i_Emin_bank_u = 1
i_Emax_bank_u = 1

do i_bank = 2, n_bank
   if (bank(i_bank)%E_ML < bank(i_Emin_bank_u)%E_ML) then
      i_Emin_bank_u = i_bank
   end if
   if (bank(i_bank)%E_ML > bank(i_Emax_bank_u)%E_ML) then
      i_Emax_bank_u = i_bank
   end if
end do

end subroutine update_min_max_bank_unit
!-------------------------------------------------------------------------------
subroutine calc_ML_energy_ligand(bank, n_bank, protein, ligand, rerank)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(inout) :: bank(:)
integer, intent(in) :: n_bank
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
integer :: i_bank
character(len=len_fname), parameter :: prefix = 'calc_ml_energy'
character(len=10*len_fname) :: command
character(len=len_fname) :: result_fn
integer :: f_unit = 523, ioerror
real(dp) :: energy
logical :: rerank

call write_bank(prefix, bank, n_bank, protein, ligand)
if (rerank) then
    command = 'python ' // trim(infile_ML(2)) // &
            ' --infile_pre_ML ' // trim(infile_pre_ML) // &
            ' --mol2_prefix ' // trim(prefix) // &
            ' --load_model ' // trim(infile_ML(4))
else
    command = 'python ' // trim(infile_ML(1)) // &
            ' --infile_pre_ML ' // trim(infile_pre_ML) // &
            ' --mol2_prefix ' // trim(prefix) // &
            ' --load_model ' // trim(infile_ML(3))
end if
call execute_command_line(command)

! remove result file of write_bank
command = 'rm '// trim(prefix) // '.E.info'
call execute_command_line(command)

! read results
result_fn = trim(prefix) // '.mol2.th1'


open(f_unit, file = trim(result_fn))

i_bank = 0
do 
    read(f_unit, "(F10.3)", iostat = ioerror) energy
    if (ioerror < 0) exit
    !
    i_bank = i_bank + 1
    bank(i_bank)%E_ML = energy
end do

close(f_unit)

command = 'rm '// trim(result_fn) 
call execute_command_line(command)

end subroutine calc_ML_energy_ligand
!-------------------------------------------------------------------------------
END MODULE LIGDOCK_CSA_UPDATE_BANK
!-------------------------------------------------------------------------------
