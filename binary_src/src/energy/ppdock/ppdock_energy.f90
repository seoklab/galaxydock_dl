!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ppdock/ppdock_energy.f90
!
! Description: Calculating ppDock originated energy terms.
!   1.   conserve_energy:
!   2-3. interev score:
!
!-------------------------------------------------------------------------------
MODULE PPDOCK_ENERGY
!-------------------------------------------------------------------------------
! global modules
use globals
use logger

! energy modules
use energy_vars

use conserve_score, only: initialize_conserve_energy, finalize_conserve_energy,&
                          calc_conserve_score
use interev, only: initialize_interev, finalize_interev, calc_interev_score
use symm_restraint, only: initialize_symm_rsr, finalize_symm_rsr, &
                          calc_symm_rsr_energy
implicit none
save
private

public :: initialize_ppdock_E
public :: finalize_ppdock_E
public :: ppdock_energy_and_gradient
public :: set_Escale_ppdock

CONTAINS
!===============================================================================
subroutine initialize_ppdock_E(protein)
!-------------------------------------------------------------------------------
! Initial setup for ppdock energy components
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein

if (use_conserve) then
    call log_p('- Energy: Conservation score activated.', me=me, level=10)    
    call initialize_conserve_energy(protein)
end if

if (use_interev) then
    call log_p('- Energy: InterEV score activated.', me=me, level=10)    
    call initialize_interev(protein)
end if

if (use_symm_rsr) then
    call log_p('- Energy: Symmertic restraint activated.', me=me, level=10)    
    call initialize_symm_rsr()
end if

end subroutine initialize_ppdock_E
!-------------------------------------------------------------------------------
subroutine finalize_ppdock_E()
!-------------------------------------------------------------------------------
if (use_conserve) call finalize_conserve_energy()
if (use_interev) call finalize_interev()
if (use_symm_rsr) call finalize_symm_rsr()

end subroutine finalize_ppdock_E
!-------------------------------------------------------------------------------
subroutine ppdock_energy_and_gradient(ff, g, calc_g, ierr)
!-------------------------------------------------------------------------------
! TODO: 
!-------------------------------------------------------------------------------
type(energy_type), intent(inout) :: ff
real(dp), intent(inout) :: g(:,:)
logical, intent(in) :: calc_g
integer, intent(inout) :: ierr
!
integer, parameter :: n_energy_output = 2
real(dp) :: f_tmp(n_energy_output), g_tmp(3, tn%atom, n_energy_output)

ff%ppdock(:) = 0.0d0
!Conservation Score
if (use_conserve .and. Escale%ppdock(1) /= 0.0d0) then
    call calc_conserve_score(f_tmp(1), g_tmp(:,:,1), appl_res, calc_g)
    ff%ppdock(1) = ff%ppdock(1) + Escale%ppdock(1)*f_tmp(1)
    g(:,:) = g(:,:) + Escale%ppdock(1)*g_tmp(:,:,1)
end if

! InterEV Score
if (use_interev .and. Escale%ppdock(2) /= 0.0d0) then
    call calc_interev_score(f_tmp(1:2), g_tmp(:,:,1:2), appl_respair, calc_g)
    ff%ppdock(2) = ff%ppdock(2) + Escale%ppdock(2)*f_tmp(1)
    ff%ppdock(3) = ff%ppdock(3) + Escale%ppdock(3)*f_tmp(2)
    g(:,:) = g(:,:) + Escale%ppdock(2)*g_tmp(:,:,1) + Escale%ppdock(3)*g_tmp(:,:,2)
end if

! Symmetric restraint Energy
! TODO: make build_appl_symm_rsr
if (use_symm_rsr .and. Escale%ppdock(4) /= 0.0d0) then
    call calc_symm_rsr_energy(f_tmp(1), g_tmp(:,:,1), appl_symm_rsr,calc_g)
    ff%ppdock(4) = ff%ppdock(4) + Escale%ppdock(4)*f_tmp(1)
    g(:,:) = g(:,:) + Escale%ppdock(4)*g_tmp(:,:,1) + Escale%ppdock(4)*g_tmp(:,:,1)
end if

ff%ppdock(0) = sum(ff%ppdock(1:n_E_component_ppdock))

end subroutine ppdock_energy_and_gradient
!-------------------------------------------------------------------------------
subroutine set_Escale_ppdock(Esch)
!-------------------------------------------------------------------------------
type(energy_type), intent(in) :: Esch

! Conservation Score
Escale%ppdock(1)    = Esch%ppdock(1)*conserve_w
! InterEVScore
Escale%ppdock(2:3)  = Esch%ppdock(2:3)*interev_w
! Symmetric restraint Energy
Escale%ppdock(4)  = Esch%ppdock(4)*symm_rsr_w

end subroutine set_Escale_ppdock
!-------------------------------------------------------------------------------
END MODULE PPDOCK_ENERGY
!-------------------------------------------------------------------------------
