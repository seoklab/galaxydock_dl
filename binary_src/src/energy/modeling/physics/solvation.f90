!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/physics/solvation.f90
!
! Description: Main solvation module to setup and call each solvation terms
! depending on the option
!
! TODO comment
!
!-------------------------------------------------------------------------------
MODULE SOLVATION
!-------------------------------------------------------------------------------

use globals
use energy_vars, only: use_solv, use_SA, use_SApp, solv_type, SA_type, appl_res, &
    SOLV_FACTS, SOLV_EEF1, SOLV_HABER, SOLV_HASEL
use facts, only: initialize_facts, calc_facts, finalize_facts
use eef1, only: initialize_eef1, calc_eef1, finalize_eef1
use surface_area, only: initialize_SA, calc_Hasel_SA, finalize_SA, refresh_SA_ref, &
                        SA_g_ref

implicit none
save
private

public :: initialize_solv
public :: finalize_solv
public :: calc_solv_energy

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_solv(molecule)
!-------------------------------------------------------------------------------
! Initialize solvation energy evaluation
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule

! If FACTS is activated (either polar or nonpolar)
if ((use_solv .and. solv_type == SOLV_FACTS) .or. (use_SA .and. SA_type == SOLV_FACTS)) then
    call initialize_facts(molecule)
end if

! If EEF1 is activated
if (use_solv .and. solv_type == SOLV_EEF1) then
    call initialize_eef1(molecule)
end if

if (use_SA .and. (SA_type == SOLV_HASEL .or. SA_type == SOLV_HABER)) then
    call initialize_SA(molecule)
end if


end subroutine initialize_solv
!-------------------------------------------------------------------------------
subroutine finalize_solv()
!-------------------------------------------------------------------------------

if ((use_solv .and. solv_type == SOLV_FACTS) .or. (use_SA .and. SA_type == SOLV_FACTS)) then
    call finalize_facts()
end if

if (use_solv .and. solv_type == SOLV_EEF1) then
    call finalize_eef1()
end if

if (use_SA .and. (SA_type == SOLV_HASEL .or. SA_type == SOLV_HABER)) then
    call finalize_SA()
end if

end subroutine finalize_solv
!-------------------------------------------------------------------------------
subroutine calc_solv_energy(f, g, appl_res, calc_g)
!-------------------------------------------------------------------------------
! TODO comment
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f(2), g(3,tn%atom,2)
logical, intent(in) :: appl_res(tn%residue), calc_g
real(dp) :: f_tmp, g_tmp(3,tn%atom)
real(dp) :: area(tn%stdatm)

f(:) = 0.0
g(:,:,:) = 0.0d0

if ((use_solv .and. solv_type == SOLV_FACTS) .or. (use_SA .and. SA_type == SOLV_FACTS)) then
    call calc_facts(f, g, appl_res, calc_g)
end if

if (use_solv .and. solv_type == SOLV_EEF1) then
    call calc_eef1(f_tmp, g_tmp, appl_res, calc_g)
    f(1) = f_tmp
    g(:,:,1) = g_tmp(:,:)
end if

if (use_SA .and. (SA_type == SOLV_HASEL .or. SA_type == SOLV_HABER)) then !for hasel or haber
    if (use_SApp) call refresh_SA_ref(calc_g)
    call calc_Hasel_SA(f_tmp, g_tmp, area, appl_res, calc_g)
    f(2) = f_tmp
    g(:,:,2) = -SA_g_ref(:,:) + g_tmp(:,:)
end if

end subroutine calc_solv_energy
!-------------------------------------------------------------------------------
END MODULE SOLVATION
