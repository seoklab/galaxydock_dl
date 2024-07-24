!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/utils/FFT_utils.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE FFT_utils
!-------------------------------------------------------------------------------

use globals

implicit none
save
public

include 'fftw3.f'

CONTAINS
!-------------------------------------------------------------------------------
subroutine generate_FFT_plan_forward(grid, n_dim, plan)
!-------------------------------------------------------------------------------
complex, intent(in) :: grid(:,:,:)
integer, intent(in) :: n_dim(3)
integer(dp), intent(out) :: plan

call sfftw_plan_dft_3d(plan, n_dim(1), n_dim(2), n_dim(3), &
                       grid(:,:,:), grid(:,:,:), FFTW_FORWARD, FFTW_ESTIMATE)

end subroutine generate_FFT_plan_forward
!-------------------------------------------------------------------------------
subroutine generate_FFT_plan_backward(grid, n_dim, plan)
!-------------------------------------------------------------------------------
complex, intent(in) :: grid(:,:,:)
integer, intent(in) :: n_dim(3)
integer(dp), intent(out) :: plan

call sfftw_plan_dft_3d(plan, n_dim(1), n_dim(2), n_dim(3), &
                       grid(:,:,:), grid(:,:,:), FFTW_BACKWARD, FFTW_ESTIMATE)

end subroutine generate_FFT_plan_backward
!-------------------------------------------------------------------------------
subroutine execute_FFT(plan, grid, n_dim)
!-------------------------------------------------------------------------------
complex, intent(inout) :: grid(:,:,:)
integer, intent(in) :: n_dim(3)
integer(dp), intent(in) :: plan

call sfftw_execute_dft(plan, grid(:,:,:), grid(:,:,:))

end subroutine execute_FFT
!-------------------------------------------------------------------------------
subroutine destroy_FFT_plan(plan)
!-------------------------------------------------------------------------------
integer(dp), intent(inout) :: plan

call sfftw_destroy_plan(plan)

end subroutine destroy_FFT_plan
!-------------------------------------------------------------------------------
END MODULE FFT_utils
!-------------------------------------------------------------------------------
