!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/physics/cutoff_methods.f90
!
! Description: Subroutines of cutoff methods for nonbonded energy calculations
!
! Now it only contains switching function and shifting function but other types
! can be implemented
!
!-------------------------------------------------------------------------------
MODULE CUTOFF_METHODS
!-------------------------------------------------------------------------------

use globals
use energy_vars, only: Ron_sqr, Roff_sqr, R3on_off, Ron_off_3

implicit none
save
public

CONTAINS
!-------------------------------------------------------------------------------
subroutine nonbonded_switch_function(r_sqr, ftr, gtr, calc_g)
!-------------------------------------------------------------------------------
! TODO comment
! gtr must be used as gtr*r
!-------------------------------------------------------------------------------
real(dp), intent(in) :: r_sqr
logical, intent(in) :: calc_g
real(dp), intent(out) :: ftr, gtr
real(dp) :: dr1, dr2, tmp1

gtr = 0.0d0
if (r_sqr <= Ron_sqr) then
    ftr = 1.0d0
else
    dr1 = Roff_sqr - r_sqr
    dr2 = R3on_off + 2.0d0*r_sqr
    tmp1 = dr1*Ron_off_3
    ftr = dr1*tmp1*dr2
    if (calc_g) then
        gtr = 4.0d0*tmp1*(dr1-dr2)
    end if
end if

end subroutine nonbonded_switch_function
!-------------------------------------------------------------------------------
subroutine nonbonded_shift_function(r_sqr, Roff_sqr, ftr, gtr, calc_g)
!-------------------------------------------------------------------------------
! TODO comment
! gtr must be used as gtr*r
!-------------------------------------------------------------------------------
real(dp), intent(in) :: r_sqr, Roff_sqr
logical, intent(in) :: calc_g
real(dp), intent(out) :: ftr, gtr
real(dp) :: tmp

gtr = 0.0d0
if (r_sqr > Roff_sqr) then
    ftr = 0.0d0
else
    tmp = 1.0d0-r_sqr/Roff_sqr
    ftr = tmp*tmp
    if (calc_g) then
        gtr = -4.0d0*tmp/Roff_sqr
    endif
end if

end subroutine nonbonded_shift_function
!-------------------------------------------------------------------------------
END MODULE CUTOFF_METHODS
!-------------------------------------------------------------------------------
