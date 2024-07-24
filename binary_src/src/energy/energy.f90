!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/energy.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE ENERGY
!-------------------------------------------------------------------------------
! global modules
use globals
use logger
use symmetry, only: copy_symm_unit
use setup_molecule, only: apply_fix_type, setup_tot_num
!
! energy modules
use energy_vars
use energy_utils, only: protein_to_R, build_appl_res
use energy_setup
use pair_list, only: initialize_pair_list, finalize_pair_list, fill_pair_index, &
                     pairlist_update
!
use modeling_energy, only: initialize_modeling_E, finalize_modeling_E, &
                           modeling_energy_and_gradient, set_Escale_modeling
use ppdock_energy,   only: initialize_ppdock_E,   finalize_ppdock_E,   &
                           ppdock_energy_and_gradient,   set_Escale_ppdock
use ligdock_energy,  only: initialize_ligdock_E,  finalize_ligdock_E,  &
                           ligdock_energy_and_gradient,  set_Escale_ligdock

implicit none
save
private

public :: initialize_energy
public :: finalize_energy
public :: reinitialize_energy
!
public :: calc_GALAXY_E
public :: calc_energy_and_gradient
!
public :: set_Escale
public :: report_energy
public :: pdb_remark_energy

CONTAINS
!===============================================================================
! Initialize/finalize energy calculation
!   - allocate/deallocate related arrays.
!   - fill arrays related to energy calculation
!-------------------------------------------------------------------------------
subroutine initialize_energy(protein, ligand, reinitialize_in)
!-------------------------------------------------------------------------------
! Initialize arrays & variables related to energy calculation
! In this subroutine, arrays used to calculate energy are allocated, and
! values for variables depend only on the atom types and connectivity are assigned.
!
! Used variables:
! 1. Cutoff for energy calculation (in energy_vars)
!   - Ron, Roff, LRoff: It should be assigned at energy_input(default/user-define)
!   - Ron_sqr, Roff_sqr, LRoff_sqr, R3on_off, Ron_off_3 are calculated here.
! 2. tn: total number related. eg. total number of residues
!     (in globals module. see tot_num_type in globals module)
! 3. R array (in energy_vars)
! 4. index for R (in energy_vars)
!   - i_R, ii_R: index & inverse index of R
!   - i_L, ii_L: index & inverse index of R for ligand
! 5. res_index (in globals)
! 6. atm_parm (in energy_vars)
! 7. pair_list (i_P, pair_index_type in energy_vars)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(inout), optional :: ligand
logical, intent(in), optional :: reinitialize_in

! cutoff constants for switching function
Ron_sqr = Ron*Ron
Roff_sqr = Roff*Roff 
Ron_off_3 = 1.0d0/((Roff_sqr-Ron_sqr)**3)
R3on_off = Roff_sqr - 3.0d0*Ron_sqr
LRoff_sqr = LRoff*LRoff
if (present(reinitialize_in)) then
    reinitialize = reinitialize_in
else
    reinitialize = .false.
end if

! allocate variables related to energy calculation
allocate(R(3,tn%atom))
allocate(i_R(2, tn%atom))
allocate(ii_R(max_atm, tn%nonlig))
if (protein%n_lig > 0) then
    allocate(i_L(2, tn%ligatm))
    allocate(ii_L(max_lig_atom, protein%n_lig))
end if
!
allocate(atm_parm(tn%atom))
allocate(res_index(tn%residue))
!
allocate(appl_res(tn%residue))
allocate(appl_respair(tn%residue, tn%residue))
!
! Building energy-applicable residue info
call build_appl_res(fix_type, is_usc)

call initialize_pair_list()

! Get index & inverse index for R. (i_R, ii_R, i_L, ii_L)
call assign_index_for_R(protein)

! Fill res_index (see res_index_type in globals.f90)
call fill_res_index_array(protein)

! Fill pair_index, i_P (see pair_index_type in energy_vars.f90)
call fill_pair_index(protein)

! Fill residual res_index & atm_parm
call initialize_atm_parm(protein)

! Fill dof_id
call setup_dof(protein)

! Initialize coordinates and update pairlists
call protein_to_R(protein)
call pairlist_update()

! Initialize modeling energy
if (use_modeling_E) then
    call initialize_modeling_E(protein)
end if

! Initialize ppdock energy
if (use_ppdock_E) then
    call initialize_ppdock_E(protein)
end if

! Initialize ligdock energy
if (use_ligdock_E) then
    if (.not. present(ligand)) then
        ! TODO: modify error message
        call terminate_with_error('ERROR: Please specify ligand')
    end if
    call initialize_ligdock_E(protein, ligand, is_usc)
end if

call set_Escale()

end subroutine initialize_energy
!-------------------------------------------------------------------------------
subroutine reinitialize_energy(protein, ligand)
!-------------------------------------------------------------------------------
! Reinitializes energy-related arrays
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(inout), optional :: ligand

call finalize_energy(reinitialize_in=.true.)

call setup_tot_num(protein)
call apply_fix_type(protein)

if (present(ligand)) then
    call initialize_energy(protein, ligand=ligand, reinitialize_in=.true.)
else
    call initialize_energy(protein, reinitialize_in=.true.)
end if

end subroutine reinitialize_energy
!-------------------------------------------------------------------------------
subroutine finalize_energy(reinitialize_in)
!-------------------------------------------------------------------------------
! Deallocate arrays related to energy calculation
!-------------------------------------------------------------------------------
logical, intent(in), optional :: reinitialize_in

if (present(reinitialize_in)) then
    reinitialize = reinitialize_in
else
    reinitialize = .false.
end if

deallocate(R)
deallocate(i_R, ii_R)
if (allocated(i_L)) then
    deallocate(i_L, ii_L)
end if
!
deallocate(res_index, atm_parm)
deallocate(dof_id)
deallocate(appl_res, appl_respair)
!
call finalize_pair_list()
!
if (use_modeling_E) then
    call finalize_modeling_E()
end if

if (use_ppdock_E) then
    call finalize_ppdock_E()
end if

if (use_ligdock_E) then
    call finalize_ligdock_E()
end if

end subroutine finalize_energy
!===============================================================================
! Setup/Calculate GALAXY energy
!===============================================================================
subroutine calc_GALAXY_E(ff, g, calc_g, ierr)
!-------------------------------------------------------------------------------
! TODO : Fill this subroutine
!-------------------------------------------------------------------------------
type(energy_type), intent(out) :: ff
real(dp), intent(out) :: g(3,tn%atom)
logical, intent(in) :: calc_g
integer, intent(out) :: ierr

! Apply symmetry if it is defined
if (symmetric) call copy_symm_unit(R)

! Update pairlist
call pairlist_update()

! Calculate energy
call calc_energy_and_gradient(ff, g, calc_g, ierr)

end subroutine calc_GALAXY_E
!-------------------------------------------------------------------------------
subroutine calc_energy_and_gradient(ff, g, calc_g, ierr)
!-------------------------------------------------------------------------------
! TODO: 
!-------------------------------------------------------------------------------
type(energy_type), intent(out) :: ff
real(dp), intent(out) :: g(:,:)
logical, intent(in) :: calc_g
integer, intent(out) :: ierr

ierr = 0
g(:,:) = 0.0d0
ff%total = 0.0d0

if (use_modeling_E) then
    call modeling_energy_and_gradient(ff, g, calc_g, ierr)
    ff%total = ff%total + ff%modeling(0)
end if

if (use_ppdock_E) then
    call ppdock_energy_and_gradient(ff, g, calc_g, ierr)
    ff%total = ff%total + ff%ppdock(0)
end if

if (use_ligdock_E) then
    call ligdock_energy_and_gradient(ff, g, calc_g, ierr)
    ff%total = ff%total + ff%ligdock(0)
end if

end subroutine calc_energy_and_gradient
!-------------------------------------------------------------------------------
subroutine set_Escale(Esch_in)
!-------------------------------------------------------------------------------
! TODO:
!-------------------------------------------------------------------------------
type(energy_type), intent(in), optional :: Esch_in
type(energy_type) :: Esch

if (present(Esch_in)) then
    Esch = Esch_in
else
    Esch%modeling(:) = 1.0d0
    Esch%ppdock  (:) = 1.0d0
    Esch%ligdock (:) = 1.0d0
end if

if (use_modeling_E) then
    call set_Escale_modeling(Esch)
else
    Escale%modeling(:) = 0.0d0
end if

if (use_ppdock_E) then
    call set_Escale_ppdock(Esch)
else
    Escale%ppdock  (:) = 0.0d0
end if

if (use_ligdock_E) then
    call set_Escale_ligdock(Esch)
else
    Escale%ligdock (:) = 0.0d0
end if

end subroutine set_Escale
!===============================================================================
! Energy calculation related Utiliities
!===============================================================================
subroutine report_energy(ff, me, log_unit, level, prefix)
!-------------------------------------------------------------------------------
type(energy_type), intent(in) :: ff
integer, intent(in), optional :: me
integer, intent(in), optional :: log_unit
integer, intent(in), optional :: level
character(len=len_log), intent(in), optional :: prefix
character(len=len_log) :: prefix_out
integer :: energy_lv

if (present(me) .and. me /= king) return
!
if (present(level)) then
    energy_lv = level
else
    energy_lv = energy_print_level
end if
!
if (present(prefix)) then
    prefix_out = prefix
else
    prefix_out = '-'
end if
!
if (present(log_unit)) then
    write(log_unit, 200) trim(prefix_out), "ENERGY.0  ", ff%total
    if (energy_lv < 10) return
    !
    write(log_unit, 210) trim(prefix_out), "ENERGY.1  ", ff%modeling(0), ff%ppdock(0), ff%ligdock(0)
    if (energy_lv < 20) return
    !
    if (use_modeling_E) then
        write(log_unit, 220) trim(prefix_out), "ENERGY.2.0", ff%modeling(1:n_E_component_modeling)
    end if
    if (use_ppdock_E) then
        write(log_unit, 220) trim(prefix_out), "ENERGY.2.1", ff%ppdock(1:n_E_component_ppdock)
    end if
    if (use_ligdock_E) then
        write(log_unit, 220) trim(prefix_out), "ENERGY.2.2", ff%ligdock(1:n_E_component_ligdock)
    end if
    if (energy_lv < 30) return
else
    write(*, 200) trim(prefix_out), "ENERGY.0  ", ff%total
    if (energy_lv < 10) return
    !
    write(*, 210) trim(prefix_out), "ENERGY.1  ", ff%modeling(0), ff%ppdock(0), ff%ligdock(0)
    if (energy_lv < 20) return
    !
    if (use_modeling_E) then
        write(*, 220) trim(prefix_out), "ENERGY.2.0", ff%modeling(1:n_E_component_modeling)
    end if
    if (use_ppdock_E) then
        write(*, 220) trim(prefix_out), "ENERGY.2.1", ff%ppdock(1:n_E_component_ppdock)
    end if
    if (use_ligdock_E) then
        write(*, 220) trim(prefix_out), "ENERGY.2.2", ff%ligdock(1:n_E_component_ligdock)
    end if
    if (energy_lv < 30) return
end if

200 format (A,1x,A10,   1x,f12.3)
210 format (A,1x,A10, 3(1x,f12.3))
220 format (A,1x,A10,23(1x,f10.3))
221 format (A,1x,A10, 3(1x,f10.3))
222 format (A,1x,A10,15(1x,f10.3))

end subroutine report_energy
!-------------------------------------------------------------------------------
subroutine pdb_remark_energy(ff, remark, level)
!-------------------------------------------------------------------------------
type(energy_type), intent(in) :: ff
character(len=len_log), intent(out) :: remark(:)
integer, intent(in), optional :: level
integer :: energy_lv
!
if (present(level)) then
    energy_lv = level
else
    energy_lv = energy_print_level
end if

write(remark(1), 200) "ENERGY.0  ", ff%total
if (energy_lv < 10) return
!
write(remark(2), 210) "ENERGY.1  ", ff%modeling(0), ff%ppdock(0), ff%ligdock(0)
if (energy_lv < 20) return
!
if (use_modeling_E) then
    write(remark(3), 220) "ENERGY.2.0", ff%modeling(1:n_E_component_modeling)
else
    write(remark(3), '(A)') "ENERGY.2.0"
end if
if (use_ppdock_E) then
    write(remark(4), 220) "ENERGY.2.1", ff%ppdock(1:n_E_component_ppdock)
else
    write(remark(4), '(A)') "ENERGY.2.1"
end if
if (use_ligdock_E) then
    write(remark(5), 220) "ENERGY.2.2", ff%ligdock(1:n_E_component_ligdock)
else
    write(remark(5), '(A)') "ENERGY.2.2"
end if

200 format (A10,   1x,f12.3)
210 format (A10, 3(1x,f12.3))
220 format (A10,23(1x,f10.3))
221 format (A10, 3(1x,f10.3))
222 format (A10,15(1x,f10.3))

end subroutine pdb_remark_energy
!-------------------------------------------------------------------------------
END MODULE ENERGY
!-------------------------------------------------------------------------------
