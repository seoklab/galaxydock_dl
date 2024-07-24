!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/modeling_energy.f90
!
! Description: calculating modeling oriented energy terms.
!   1-4.   Bonded energy: bond(1), angle(2), torsion(3), and improper torsion(4)
!   5.     van der Waals energy
!   6.     Coulomb energy
!   7-8.   Solvation energy: polar(7) and apolar/SA(8)
!   9.     Rotamer preference score
!   10.    Ramachandran preference score
!   11.    Hydrogen bond energy
!   12-13. dDFIRE statistical potential: distance-dep.(12) and angle-dep.(13)
!   12-13. KGB energy: An alternative for dDFIRE with consideration of solv state
!   14-15. GOAP statistical potential: distance-dep.(14) and angle-dep.(15)
!   16-20. Restraints: CA-CA(16), N-O(17), BB-SC(18), SC-SC(19), and ligand(20)
!
!-------------------------------------------------------------------------------
MODULE MODELING_ENERGY
!-------------------------------------------------------------------------------
! global modules
use globals
use logger

! energy modules
use energy_vars
use cutoff_methods

! modeling/physics
use molecular_mechanics, only: initialize_MM, finalize_MM, calc_MM_energy
use solvation, only: initialize_solv, finalize_solv, calc_solv_energy

! modeling/knowledge
use ddfire,        only: initialize_ddfire, finalize_ddfire, calc_ddfire_score
use goap,          only: initialize_goap, finalize_goap, calc_goap_energy
use rotamer_score, only: calc_rotamer_score
use ramachandran_score,  only: initialize_rama_score, finalize_rama_score, calc_rama_score
use hbond_energy,  only: initialize_hbond, finalize_hbond, calc_hbond_energy
use knowledge_gb, only: initialize_kgb_energy, finalize_kgb_energy, calc_kgb_energy

! modeling/restraint
use restraint, only: initialize_rsr, finalize_rsr, calc_rsr_energy, calc_rsr_energy_by_meld
use distogram, only: initialize_distogram, finalize_distogram, calc_distogram_energy
use bb_torsion, only: initialize_bb_torsion, finalize_bb_torsion, calc_bb_torsion_energy
use ml_score

implicit none
save
private

public :: initialize_modeling_E
public :: finalize_modeling_E
public :: modeling_energy_and_gradient
public :: set_Escale_modeling

CONTAINS
!===============================================================================
subroutine initialize_modeling_E(protein)
!-------------------------------------------------------------------------------
! Initial setup for modeling energy components
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein

if (use_softsphere .or. use_LJ .or. use_vdw_scwrl4 .or. use_elec .or. use_bond) then
    call log_p('- Energy: MM energy activated.', me=me, level=10)    
    call initialize_MM(protein)
end if

if (use_solv .or. use_SA) then
    call log_p('- Energy: Solvation energy activated.', me=me, level=10)    
    call initialize_solv(protein)
end if

if (use_rotamer_score) then
    call log_p('- Energy: Rotamer score activated.', me=me, level=10)
end if

if (use_rama_score) then
    call log_p('- Energy: Ramachandran score activated.', me=me, level=10)
    call initialize_rama_score(protein)
end if

if (use_hbond) then
    call log_p('- Energy: Rosetta H-bond energy activated.', me=me, level=10)    
    call initialize_hbond(protein)
end if

if (use_ddfire) then
    call log_p('- Energy: dDFIRE energy activated.', me=me, level=10)    
    call initialize_ddfire(protein)
end if

if (use_kgb) then
    call log_p('- Energy: KGB energy activated.', me=me, level=10)    
    call initialize_kgb_energy(protein)
end if

if (use_goap) then
    call log_p('- Energy: GOAP energy activated.', me=me, level=10)    
    call initialize_goap(protein)
end if

if (use_rsr) then
    call log_p('- Energy: Restraint energy activated.', me=me, level=10)
    call initialize_rsr(protein)
end if

if (use_distogram) then
    call log_p('- Energy: Distogram energy activated.', me=me, level=10)
    call initialize_distogram(protein)
end if

if (use_bb_torsion) then
    call log_p('- Energy: bb_torsion energy activated.', me=me, level=10)
    call initialize_bb_torsion(protein)
end if

if (use_ml) then
    call log_p('- Energy: Ml energy activated.', me=me, level=10)
    call initialize_ml_score(protein)
end if

end subroutine initialize_modeling_E
!-------------------------------------------------------------------------------
subroutine finalize_modeling_E()
!-------------------------------------------------------------------------------

if (use_bond .or. use_softsphere .or. use_LJ .or. use_vdw_scwrl4 .or. use_elec) &
    call finalize_MM()
if (use_solv .or. use_SA) call finalize_solv()
if (use_hbond)  call finalize_hbond()
if (use_rama_score)   call finalize_rama_score()
if (use_ddfire) call finalize_ddfire()
if (use_kgb)    call finalize_kgb_energy()
if (use_goap)   call finalize_goap()
if (use_rsr)    call finalize_rsr()
if (use_distogram)    call finalize_distogram()
if (use_bb_torsion)   call finalize_bb_torsion()
if (use_ml)     call finalize_ml_score()

end subroutine finalize_modeling_E
!-------------------------------------------------------------------------------
subroutine modeling_energy_and_gradient(ff, g, calc_g, ierr)
!-------------------------------------------------------------------------------
! TODO: 
!-------------------------------------------------------------------------------
type(energy_type), intent(inout) :: ff
real(dp), intent(inout) :: g(:,:)
logical, intent(in) :: calc_g
integer, intent(inout) :: ierr
!
integer, parameter :: n_energy_output = 6
real(dp) :: f_tmp(n_energy_output), g_tmp(3, tn%atom, n_energy_output)

ff%modeling(:) = 0.0d0

!TODO: in original code, calc_MM_energy was done always (default on)
!should this be always calculated? !Escale option?
if ((use_bond .or. use_softsphere .or. use_LJ .or. use_vdw_scwrl4 .or. use_elec) .and. &
    (Escale%modeling(1) /= 0.0d0 .or. Escale%modeling(2) /= 0.0d0 .or. &
     Escale%modeling(3) /= 0.0d0 .or. Escale%modeling(4) /= 0.0d0 .or. &
     Escale%modeling(5) /= 0.0d0 .or. Escale%modeling(6) /= 0.0d0)) then
    call calc_MM_energy(f_tmp(1:6), g_tmp(:,:,1:6), appl_respair, skip_fixed, calc_g)
    ff%modeling(1) = ff%modeling(1) + Escale%modeling(1)*f_tmp(1)
    ff%modeling(2) = ff%modeling(2) + Escale%modeling(2)*f_tmp(2)
    ff%modeling(3) = ff%modeling(3) + Escale%modeling(3)*f_tmp(3)
    ff%modeling(4) = ff%modeling(4) + Escale%modeling(4)*f_tmp(4)
    ff%modeling(5) = ff%modeling(5) + Escale%modeling(5)*f_tmp(5)
    ff%modeling(6) = ff%modeling(6) + Escale%modeling(6)*f_tmp(6)
    g(:,:) = g(:,:) + Escale%modeling(1)*g_tmp(:,:,1) &
                    + Escale%modeling(2)*g_tmp(:,:,2) &
                    + Escale%modeling(3)*g_tmp(:,:,3) &
                    + Escale%modeling(4)*g_tmp(:,:,4) &
                    + Escale%modeling(5)*g_tmp(:,:,5) &
                    + Escale%modeling(6)*g_tmp(:,:,6)
end if

if ((use_solv .or. use_SA) .and. &
    (Escale%modeling(7) /= 0.0d0 .or. Escale%modeling(8) /= 0.0d0)) then
    call calc_solv_energy(f_tmp(1:2), g_tmp(:,:,1:2), appl_res, calc_g)
    ff%modeling(7) = ff%modeling(7) + Escale%modeling(7)*f_tmp(1)
    ff%modeling(8) = ff%modeling(8) + Escale%modeling(8)*f_tmp(2)
    g(:,:) = g(:,:) + Escale%modeling(7)*g_tmp(:,:,1) &
                    + Escale%modeling(8)*g_tmp(:,:,2)
end if

if (use_rotamer_score .and. Escale%modeling(9) /= 0.0d0) then
    call calc_rotamer_score(f_tmp(1), g_tmp(:,:,1), appl_res, calc_g)
    ff%modeling(9) = ff%modeling(9) + Escale%modeling(9)*f_tmp(1)
    g(:,:) = g(:,:) + Escale%modeling(9)*g_tmp(:,:,1)
end if

if (use_rama_score .and. Escale%modeling(10) /= 0.0d0) then
    call calc_rama_score(f_tmp(1), g_tmp(:,:,1), appl_res, calc_g)
    ff%modeling(10) = ff%modeling(10) + Escale%modeling(10)*f_tmp(1)
    g(:,:) = g(:,:) + Escale%modeling(10)*g_tmp(:,:,1)
end if

if (use_hbond .and. Escale%modeling(11) /= 0.0d0) then
    call calc_hbond_energy(f_tmp(1), g_tmp(:,:,1), appl_respair, calc_g)
    ff%modeling(11) = ff%modeling(11) + Escale%modeling(11)*f_tmp(1)
    g(:,:) = g(:,:) + Escale%modeling(11)*g_tmp(:,:,1)
end if

if (use_ddfire .and. Escale%modeling(12) /= 0.0d0) then
    call calc_ddfire_score(f_tmp(1:2), g_tmp(:,:,1:2), appl_respair, calc_g)
    ff%modeling(12) = ff%modeling(12) + Escale%modeling(12)*f_tmp(1)
    ff%modeling(13) = ff%modeling(13) + Escale%modeling(13)*f_tmp(2)
    g(:,:) = g(:,:) + Escale%modeling(12)*g_tmp(:,:,1) &
                    + Escale%modeling(13)*g_tmp(:,:,2)
end if

if (use_kgb .and. Escale%modeling(12) /= 0.0d0) then
    call calc_kgb_energy(f_tmp(1:2), g_tmp(:,:,1:2), appl_respair, calc_g)
    ff%modeling(12) = ff%modeling(12) + Escale%modeling(12)*f_tmp(1)
    ff%modeling(13) = ff%modeling(13) + Escale%modeling(13)*f_tmp(2)
    g(:,:) = g(:,:) + Escale%modeling(12)*g_tmp(:,:,1) &
                    + Escale%modeling(13)*g_tmp(:,:,2)
end if

if (use_goap .and. Escale%modeling(14) /= 0.0d0) then
    call calc_goap_energy(f_tmp(1:2), g_tmp(:,:,1:2), appl_respair, calc_g)
    ff%modeling(14) = ff%modeling(14) + Escale%modeling(14)*f_tmp(1)
    ff%modeling(15) = ff%modeling(15) + Escale%modeling(15)*f_tmp(2)
    g(:,:) = g(:,:) + Escale%modeling(14)*g_tmp(:,:,1) &
                    + Escale%modeling(15)*g_tmp(:,:,2)
end if

!TODO Escale option??
if (use_rsr .and. (Escale%modeling(16) /= 0.0d0 .or. &
    Escale%modeling(17) /= 0.0d0 .or. Escale%modeling(18) /= 0.0d0 .or. &
    Escale%modeling(19) /= 0.0d0 .or. Escale%modeling(20) /= 0.0d0)) then
    if (use_meld) then
        call calc_rsr_energy_by_meld(f_tmp(1:5), g_tmp(:,:,1:5), &
            max_n_features, max_n_grp, meld_active, calc_g)
    else
        call calc_rsr_energy(f_tmp(1:5), g_tmp(:,:,1:5), calc_g)
    end if
    ff%modeling(16) = ff%modeling(16) + Escale%modeling(16)*f_tmp(1)
    ff%modeling(17) = ff%modeling(17) + Escale%modeling(17)*f_tmp(2)
    ff%modeling(18) = ff%modeling(18) + Escale%modeling(18)*f_tmp(3)
    ff%modeling(19) = ff%modeling(19) + Escale%modeling(19)*f_tmp(4)
    ff%modeling(20) = ff%modeling(20) + Escale%modeling(20)*f_tmp(5)
    g(:,:) = g(:,:) + Escale%modeling(16)*g_tmp(:,:,1) &
                    + Escale%modeling(17)*g_tmp(:,:,2) &
                    + Escale%modeling(18)*g_tmp(:,:,3) &
                    + Escale%modeling(19)*g_tmp(:,:,4) &
                    + Escale%modeling(20)*g_tmp(:,:,5)
end if

if (use_distogram .and. Escale%modeling(21) /= 0.0d0) then
    call calc_distogram_energy(f_tmp(1), g_tmp(:,:,1), calc_g)
    ff%modeling(21) = ff%modeling(21) + Escale%modeling(21)*f_tmp(1)
    g(:,:) = g(:,:) + Escale%modeling(21)*g_tmp(:,:,1)
end if

if (use_bb_torsion .and. Escale%modeling(22) /= 0.0d0) then
    call calc_bb_torsion_energy(f_tmp(1), g_tmp(:,:,1), appl_res, calc_g)
    ff%modeling(22) = ff%modeling(22) + Escale%modeling(22)*f_tmp(1)
    g(:,:) = g(:,:) + Escale%modeling(22)*g_tmp(:,:,1)
end if

if (use_ml .and. Escale%modeling(23) /= 0.0d0) then
    call calc_ml_score(f_tmp(1), g_tmp(:,:,1), appl_res, calc_g)
    ff%modeling(23) = ff%modeling(23) + Escale%modeling(23)*f_tmp(1)
    g(:,:) = g(:,:) + Escale%modeling(23)*g_tmp(:,:,1)
end if

ff%modeling(0) = sum(ff%modeling(1:n_E_component_modeling))

end subroutine modeling_energy_and_gradient
!-------------------------------------------------------------------------------
subroutine set_Escale_modeling(Esch)
!-------------------------------------------------------------------------------
type(energy_type), intent(in) :: Esch

! Molecular mechanics (bond/angle/torsion/improper)
Escale%modeling(1:4)   = Esch%modeling(1:4)*bond_w
! Molecular mechanics (vdw)
Escale%modeling(5)     = Esch%modeling(5)*vdw_w
! Molecular mechanics (coulomb)
Escale%modeling(6)     = Esch%modeling(6)*elec_w
! Solvation
Escale%modeling(7)     = Esch%modeling(7)*solv_w
! Surface area
Escale%modeling(8)     = Esch%modeling(8)*SA_w
! Rotamer score
Escale%modeling(9)     = Esch%modeling(9)*rot_w
! Ramachandran
Escale%modeling(10)    = Esch%modeling(10)*rama_w
! Hbond energy
Escale%modeling(11)    = Esch%modeling(11)*hbond_w
! dDFIRE
if (use_ddfire) then
    Escale%modeling(12)    = Esch%modeling(12)*dfire_w
    Escale%modeling(13)    = Esch%modeling(13)*dfire_w*ddfire_add_scale
else if (use_kgb) then
    Escale%modeling(12)    = Esch%modeling(12)*kgb_w
    Escale%modeling(13)    = Esch%modeling(13)*kgb_w
end if
! GOAP
Escale%modeling(14:15) = Esch%modeling(14:15)*goap_w
! Restraint energy
Escale%modeling(16:20) = Esch%modeling(16:20)*rsr_w
! Distogram
Escale%modeling(21) = Esch%modeling(21)*distogram_w
! bb_torsion
Escale%modeling(22) = Esch%modeling(22)*bb_torsion_w
! Machine Learning
Escale%modeling(23) = Esch%modeling(23)*ml_w

end subroutine set_Escale_modeling
!-------------------------------------------------------------------------------
END MODULE MODELING_ENERGY
!-------------------------------------------------------------------------------
