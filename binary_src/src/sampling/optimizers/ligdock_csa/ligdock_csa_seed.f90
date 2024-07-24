!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/sampling/optimizers/ligdock_csa/ligdock_csa_seed.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE LIGDOCK_CSA_SEED
!-------------------------------------------------------------------------------

use globals
use logger
!
use ligdock_csa_vars

implicit none
save
private

public :: select_seeds
public :: done_seed_cycle
public :: bank_used_reset
public :: bank_all_used

CONTAINS
!===============================================================================
! Subroutines related to select seeds
!===============================================================================
subroutine select_seeds(bank, n_bank, Dij)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(inout) :: bank(:)
integer, intent(in) :: n_bank
real(dp), intent(in) :: Dij(:,:)
!
integer :: idx_unused(max_bank)
integer :: i_unused, i_bank, i_seed
logical :: select_used

call log_divider(level=30)
call log_p("  Select Seeds....", level=30)

! List of bank members (idx_unused) that were not used before
i_unused = 0
do i_bank = 1, n_bank
    if (.not. bank(i_bank)%used) then
        i_unused = i_unused + 1
        idx_unused(i_unused) = i_bank
    end if
end do

n_unused = i_unused
write(log_msg,"(6X,A,I4,A,I4)") "Unused seeds : ", n_unused, " / ", n_bank
call log_p(log_msg, level=30)

! Select n_seed seeds from the bank members

idx_seed(:) = 0

if (n_unused < n_seed) then
    ! select all the unused bank members as seed
    idx_seed(1:n_unused) = idx_unused(1:n_unused)
    ! select the remaining seeds out of the used bank members
    select_used = .true.
    do i_seed = n_unused + 1, n_seed
        idx_seed(i_seed) = select_a_seed(bank, n_bank, Dij, i_seed, select_used)
    end do
else if (n_unused == n_seed) then
    ! select all the unused bank members as seed
    idx_seed(1:n_unused) = idx_unused(1:n_unused)
else
    ! select from the unused bank members
    select_used = .false.
    do i_seed = 1, n_seed
        idx_seed(i_seed) = select_a_seed(bank, n_bank, Dij, i_seed, select_used)
    end do
end if

do i_seed = 1, n_seed
    bank(idx_seed(i_seed))%used = .true.
end do

end subroutine select_seeds
!-------------------------------------------------------------------------------
integer function select_a_seed(bank, n_bank, Dij, seed_idx, select_used)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: bank(:)
integer, intent(in) :: n_bank, seed_idx
real(dp), intent(in) :: Dij(:,:)
logical, intent(in) :: select_used
!
integer :: idx_list(max_bank)
integer :: i_list, n_list, i_seed
integer :: i_bank, j_bank
real(dp) :: min_dist_list(max_bank)
real(dp) :: dist, min_dist, sum_min_dist, ave_min_dist
real(dp) :: Emin

! get min distance from the already selected seeds
i_list = 0
sum_min_dist = 0.0d0
min_dist = 0.0d0

do i_bank = 1, n_bank
    if (bank(i_bank)%used .eqv. select_used) then
        do i_seed = 1, seed_idx - 1
            j_bank = idx_seed(i_seed)
            dist = Dij(i_bank,j_bank)
            if (i_seed == 1 .or. dist < min_dist) min_dist = dist
        end do
        i_list = i_list + 1
        min_dist_list(i_list) = min_dist
        idx_list(i_list) = i_bank
        sum_min_dist = sum_min_dist + min_dist
    end if
end do

n_list = i_list
ave_min_dist = sum_min_dist/n_list

! pick a seed from the list (min E with min distance greater than average)
Emin = e1max
do i_list = 1, n_list
    if (min_dist_list(i_list) >= ave_min_dist) then
        i_bank = idx_list(i_list)
        if (bank(i_bank)%E(0) < Emin) then
            Emin = bank(i_bank)%E(0)
            select_a_seed = i_bank
        end if
    end if
end do

end function select_a_seed
!===============================================================================
! subroutines related to check used seeds
!===============================================================================
logical function done_seed_cycle(bank, n_bank) 
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: bank(:)
integer, intent(in) :: n_bank
integer :: i
integer :: unused_seeds

unused_seeds = 0
i = 1
do i = 1, n_bank
   if (.not. bank(i)%used) unused_seeds = unused_seeds + 1
end do
if (unused_seeds <= cut_pos_seed_size) then
   done_seed_cycle = .true.
else
   done_seed_cycle = .false.
end if

end function done_seed_cycle
!-------------------------------------------------------------------------------
subroutine bank_used_reset(bank, n_bank) 
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(inout) :: bank(:)
integer, intent(in) :: n_bank
integer :: i

do i = 1, n_bank
   bank(i)%used = .false.
end do

end subroutine bank_used_reset
!-------------------------------------------------------------------------------
subroutine bank_all_used(bank, n_bank) 
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(inout) :: bank(:)
integer, intent(in) :: n_bank
integer :: i

do i = 1, n_bank
   bank(i)%used=.true.
end do

end subroutine bank_all_used
!-------------------------------------------------------------------------------
END MODULE LIGDOCK_CSA_SEED
!-------------------------------------------------------------------------------
