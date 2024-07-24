!-------------------------------------------------------------------------------
MODULE SYMM_RESTRAINT
!-------------------------------------------------------------------------------
use globals
use energy_vars, only: R

implicit none
save
private

real(dp), allocatable :: rsr_list(:)

public :: initialize_symm_rsr
public :: finalize_symm_rsr
public :: calc_symm_rsr_energy
CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_symm_rsr()
!-------------------------------------------------------------------------------
! Initialize hormonic restraint function for Z-axis symmetry.
!-------------------------------------------------------------------------------
integer :: i_res, i_atm

! allocate variables
allocate(rsr_list(tn%stdres))

do i_res = 1, tn%stdres
    i_atm = res_index(i_res)%Ca_id(1)
    rsr_list(i_res) = dot_product(R(1:2,i_atm), R(1:2,i_atm))
    rsr_list(i_res) = sqrt(rsr_list(i_res))
end do

end subroutine initialize_symm_rsr
!-------------------------------------------------------------------------------
subroutine finalize_symm_rsr()
!-------------------------------------------------------------------------------
deallocate(rsr_list)

end subroutine finalize_symm_rsr
!-------------------------------------------------------------------------------
subroutine calc_symm_rsr_energy(f_tmp, g_tmp, is_rsrres, calc_g)
!-------------------------------------------------------------------------------
! Harmonic restraint function for Z-axis symmetry.
!-------------------------------------------------------------------------------
logical, intent(in) :: is_rsrres(tn%stdres), calc_g
real(dp), intent(out) :: f_tmp, g_tmp(3,tn%stdatm)
integer :: i_res, i_atm
real(dp) :: rootxy

f_tmp = 0.0d0
g_tmp(:,:) = 0.0d0
do i_res = 1,tn%stdres
    if (.not. is_rsrres(i_res)) cycle

    i_atm = res_index(i_res)%Ca_id(1)
    rootxy = sqrt(R(1,i_atm)**2+R(2,i_atm)**2)
    f_tmp = f_tmp + (rootxy - rsr_list(i_res))**2
    if (calc_g) then
        g_tmp(1:2,i_atm) = (1.0d0-rsr_list(i_res)/rootxy)*2.0d0*R(1:2,i_atm)
    end if
end do

end subroutine calc_symm_rsr_energy
!-------------------------------------------------------------------------------
END MODULE SYMM_RESTRAINT
!-------------------------------------------------------------------------------
