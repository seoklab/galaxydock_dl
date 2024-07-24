!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ligdock/ligdock_energy.f90
!
! Description: calculating ligand docking oriented energy terms.
!
!-------------------------------------------------------------------------------
MODULE LIGDOCK_ENERGY
!-------------------------------------------------------------------------------
! global modules
use globals
use logger

! energy modules
use energy_vars
use energy_utils, only: protein_to_R
use pair_list, only: pairlist_update
use ligdock_E_vars
use ligdock_E_utils, only: read_aa_to_mol2
!
use restraint, only: construct_ligdock_rsr_grid, calc_ligdock_rsr_using_grid
use autodock3, only: initialize_atdk3_energy, finalize_atdk3_energy, &
                     construct_atdk3_grid, calc_atdk3_intrxn_E_using_grid, &
                     calc_atdk3_intern_E_using_table, calc_atdk3_inter_energy, &
                     calc_atdk3_intra_energy, check_clash_atdk3, &
                     calc_atdk3_prot_vdw, set_atdk3_hbond_flex, &
                     FFT_rec_grid_atdk3, FFT_lig_grid_atdk3
use autodock4, only: initialize_atdk4_energy, finalize_atdk4_energy, &
                     construct_atdk4_grid, calc_atdk4_intrxn_E_using_grid, &
                     calc_atdk4_intern_E_using_table, calc_atdk4_inter_energy, &
                     calc_atdk4_intra_energy, check_clash_atdk4, &
                     set_atdk4_hbond_flex, FFT_rec_grid_atdk4, FFT_lig_grid_atdk4
use drugscore, only: initialize_drugscore, finalize_drugscore, &
                     construct_drugscore_grid, calc_DS_intrxn_E_using_grid, &
                     calc_DrugScore, FFT_rec_grid_drugscore, FFT_lig_grid_drugscore
use Xscore, only: initialize_X_score, finalize_X_score, &
                  construct_HM_grid, calc_HM_score_using_grid, calc_HM_score, &
                  FFT_rec_grid_HM, FFT_lig_grid_HM
use PLP_tor, only: initialize_PLP_tor, finalize_PLP_tor, calc_PLP_tor
use ROTA, only: initialize_rota, finalize_ROTA, calc_ROTA_score
use hybrid_E, only: construct_prot_table, finalize_prot_table, &
                    calc_prot_E_using_table

implicit none
save
private

public :: initialize_ligdock_E
public :: finalize_ligdock_E
public :: construct_dock_grid
public :: ligdock_energy_using_grid
public :: ligdock_energy_and_gradient
public :: check_clash_within_ligand
public :: update_prot_prot_intrxn
public :: set_atdk_hbond_flex
public :: set_Escale_ligdock
public :: receptor_gridization
public :: ligand_gridization
public :: finalize_FFT_grid

CONTAINS
!===============================================================================
subroutine initialize_ligdock_E(protein, ligand, flex_res)
!-------------------------------------------------------------------------------
! Initial setup for ligdock energy components
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
logical, intent(in) :: flex_res(:)
type(residue_mol2_type) :: res_mol2_info(max_mol2_res)
integer :: num_mol2_info
!
! Read data file to assign mol2 type of protein atoms
if (use_atdk4 .or. use_drugscore) then
    infile_aa_to_mol2 = trim(data_dir) // 'X_score_res.prm'
    call read_aa_to_mol2(infile_aa_to_mol2, res_mol2_info, num_mol2_info)
end if

if (use_atdk3) then
    call log_p('- Energy: AutoDock3 energy activated.', me=me, level=10)
    call initialize_atdk3_energy(protein, ligand)
end if

if (use_atdk4) then
    call log_p('- Energy: AutoDock4 energy activated.', me=me, level=10)
    call initialize_atdk4_energy(protein, ligand, tn%atom, &
                                 res_mol2_info, num_mol2_info)
end if

if (use_drugscore) then
    call log_p('- Energy: DrugScore activated.', me=me, level=10)
    call initialize_drugscore(protein, ligand, tn%atom, flex_res, &
                              res_mol2_info, num_mol2_info)
end if

if (use_Xscore) then
    call log_p('- Energy: X-score activated.', me=me, level=10)
    call initialize_X_score(protein, ligand)
end if

if (use_PLP) then
    call log_p('- Energy: PLP_tor activated.', me=me, level=10)
    call initialize_PLP_tor(ligand)
end if

if (use_ROTA .and. n_usc > 0) then
    call log_P('- Energy: ROTA activated.', me=me, level=10)
    call initialize_rota(protein)
end if
!
!! TODO: consider restraint on protein-ligand intrxn.
!
end subroutine initialize_ligdock_E
!-------------------------------------------------------------------------------
subroutine finalize_ligdock_E()
!-------------------------------------------------------------------------------

if (use_atdk3) call finalize_atdk3_energy()
if (use_atdk4) call finalize_atdk4_energy()
if (use_drugscore) call finalize_drugscore()
if (use_Xscore) call finalize_X_Score()
if (use_ROTA .and. n_usc > 0) call finalize_rota()
if (use_ROTA .and. n_usc > 0) call finalize_prot_table()

end subroutine finalize_ligdock_E
!-------------------------------------------------------------------------------
subroutine construct_dock_grid(protein, ligand)
!-------------------------------------------------------------------------------
! Construct docking grid to calculate docking energy more faster.
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
integer :: ref_no

if (.not. use_input_cntr) then
    ref_no = ligand%lig_type
    dock_grid_info%grid_cntr(:) = ref_lig(ref_no)%R(:,ref_lig(ref_no)%cntr_atm)
end if
if (use_atdk3) then
    call log_p('- Construct autodock3 docking grid.', me=me, level=10)
    call construct_atdk3_grid(dock_grid_info, dock_grid, is_usc)
else if (use_atdk4) then
    call log_p('- Construct autodock4 docking grid.', me=me, level=10)
    call construct_atdk4_grid(dock_grid_info, dock_grid, is_usc)
end if

if (use_drugscore) then
    call log_p('- Construct DrugScore docking grid.', me=me, level=10)
    call construct_drugscore_grid(dock_grid_info, dock_grid, is_usc)
end if

if (use_Xscore) then
    call log_p('- Construct HM-Score docking grid.', me=me, level=10)
    call construct_HM_grid(dock_grid_info, dock_grid, is_usc)
end if

if (use_rsr) then
    call log_p('- Construct ligdock-rsr docking grid.', me=me, level=10)
    call construct_ligdock_rsr_grid(ligand, dock_grid_info, dock_grid, is_usc)
end if


if (use_ROTA .and. n_usc > 0) then
    call log_p('- Construct protein energy tables.', me=me, level=10)
    call construct_prot_table(protein, is_usc)
end if

call visualize_grid_boundary()

end subroutine construct_dock_grid
!-------------------------------------------------------------------------------
subroutine visualize_grid_boundary()
!-------------------------------------------------------------------------------
! To see grid boundary in viewer, write grid boundary to pdb format
!-------------------------------------------------------------------------------
character(len=len_fname) :: box_file
integer :: o_unit = 15
integer :: ioerror
real(dp), dimension(2,3) :: box_crd

box_crd(1,:) = dock_grid_info%grid_cntr(:) &
             - dock_grid_info%grid_width * dble((dock_grid_info%n_elem(:)-1)/2.0d0)
box_crd(2,:) = dock_grid_info%grid_cntr(:) &
             + dock_grid_info%grid_width * dble((dock_grid_info%n_elem(:)-1)/2.0d0)

box_file = 'box.pdb'

open(unit=o_unit, file=trim(box_file), action='write', status='replace', iostat=ioerror)
if(ioerror /= 0) stop 'cannot create box_file'

! print info
write(o_unit, '(a6, 4x, a15, 3f8.3)') 'REMARK', 'Center of box :', &
                                      dock_grid_info%grid_cntr(1:3)
! print crd
write(o_unit, '(a4, 6x, a1, 2x, a3, 1x, a3, 5x, a1, 4x, 3f8.3)') &
     'ATOM', '1', 'PT1', 'BOX', '1', box_crd(2, 1:3)
write(o_unit, '(a4, 6x, a1, 2x, a3, 1x, a3, 5x, a1, 4x, 3f8.3)') &
     'ATOM', '2', 'PT2', 'BOX', '1', box_crd(2, 1:2), box_crd(1,3)
write(o_unit, '(a4, 6x, a1, 2x, a3, 1x, a3, 5x, a1, 4x, 3f8.3)') &
     'ATOM', '3', 'PT3', 'BOX', '1', box_crd(2,1), box_crd(1,2), box_crd(2,3)
write(o_unit, '(a4, 6x, a1, 2x, a3, 1x, a3, 5x, a1, 4x, 3f8.3)') &
     'ATOM', '4', 'PT4', 'BOX', '1', box_crd(2,1), box_crd(1, 2:3)
write(o_unit, '(a4, 6x, a1, 2x, a3, 1x, a3, 5x, a1, 4x, 3f8.3)') &
     'ATOM', '5', 'PT5', 'BOX', '1', box_crd(1,1), box_crd(2, 2:3)
write(o_unit, '(a4, 6x, a1, 2x, a3, 1x, a3, 5x, a1, 4x, 3f8.3)') &
     'ATOM', '6', 'PT6', 'BOX', '1', box_crd(1,1), box_crd(2,2), box_crd(1,3)
write(o_unit, '(a4, 6x, a1, 2x, a3, 1x, a3, 5x, a1, 4x, 3f8.3)') &
     'ATOM', '7', 'PT7', 'BOX', '1', box_crd(1, 1:2), box_crd(2,3)
write(o_unit, '(a4, 6x, a1, 2x, a3, 1x, a3, 5x, a1, 4x, 3f8.3)') &
     'ATOM', '8', 'PT8', 'BOX', '1', box_crd(1, 1:3)
! print connet info
write(o_unit, '(a26)') 'CONECT    1    3    2    5'
write(o_unit, '(a26)') 'CONECT    2    4    6    1'
write(o_unit, '(a26)') 'CONECT    3    4    7    1'
write(o_unit, '(a26)') 'CONECT    4    8    3    2'
write(o_unit, '(a26)') 'CONECT    5    7    6    1'
write(o_unit, '(a26)') 'CONECT    6    8    2    5'
write(o_unit, '(a26)') 'CONECT    7    8    3    5'
write(o_unit, '(a26)') 'CONECT    8    4    7    6'

close(o_unit)
!-------------------------------------------------------------------------------
end subroutine visualize_grid_boundary
!-------------------------------------------------------------------------------
subroutine ligdock_energy_using_grid(protein, ligand, ff, g, calc_g)
!-------------------------------------------------------------------------------
! Calculste docking energy using docking grid.
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
type(energy_type), intent(inout) :: ff
real(dp), intent(inout) :: g(:,:)
logical, intent(in) :: calc_g
real(dp) :: f(4), g_tmp(3, ligand%n_atm), g_all(3, tn%atom, 6)
integer :: i, i_atm, atm_no

ff%ligdock(:) = 0.0d0
if (calc_g) then
    g = 0.0d0
end if

if (use_atdk3) then
    if (atdk_w > small_real) then
        call calc_atdk3_intrxn_E_using_grid(ligand, dock_grid_info, &
                                            dock_grid%atdk_grid, ff%ligdock(1),&
                                            g_tmp, calc_g)
        if (n_usc > 0) then
            call calc_atdk3_inter_energy(f(1:4), g_all(:,:,1), is_usc, .false.)
            ff%ligdock(1) = ff%ligdock(1) + sum(f(1:4))
            if (calc_g) then
                do i_atm = 1, ligand%n_atm
                    atm_no = ii_L(i_atm, ligand%lig_no)
                    g(:,i_atm) = g(:,i_atm) + atdk_w * g_all(:,atm_no,1)
                end do
            end if
        end if
        ff%ligdock(1) = atdk_w * ff%ligdock(1)
        if (calc_g) then
            g(:,:) = g(:,:) + atdk_w*g_tmp(:,:)
        end if
    end if
    if (atdk_int_w > small_real) then
        call calc_atdk3_intern_E_using_table(ligand, ff%ligdock(2), g_tmp, &
                                             calc_g)
        ff%ligdock(2) = atdk_int_w * ff%ligdock(2)
        if (calc_g) then
            g(:,:) = g(:,:) + atdk_int_w*g_tmp(:,:)
        end if
    end if
end if

if (use_atdk4) then
    if (atdk_w > small_real) then
        call calc_atdk4_intrxn_E_using_grid(ligand, dock_grid_info, &
                                            dock_grid%atdk_grid, ff%ligdock(1),&
                                            g_tmp, calc_g)
        if (n_usc > 0) then
            call calc_atdk4_inter_energy(f(1:4), g_all(:,:,1), is_usc, calc_g)
            ff%ligdock(1) = ff%ligdock(1) + sum(f(1:4))
            if (calc_g) then
                do i_atm = 1, ligand%n_atm
                    atm_no = ii_L(i_atm, ligand%lig_no)
                    g(:,i_atm) = g(:,i_atm) + atdk_w * g_all(:,atm_no,1)
                end do
            end if
        end if
        ff%ligdock(1) = atdk_w * ff%ligdock(1)
        if (calc_g) then
            g(:,:) = g(:,:) + atdk_w*g_tmp(:,:)
        end if
    end if
    if (atdk_int_w > small_real) then
        call calc_atdk4_intern_E_using_table(ligand, ff%ligdock(2), g_tmp, &
                                             calc_g)
        ff%ligdock(2) = atdk_int_w * ff%ligdock(2)
        if (calc_g) then
            g(:,:) = g(:,:) + atdk_int_w*g_tmp(:,:)
        end if
    end if
end if

if (use_drugscore) then
    call calc_DS_intrxn_E_using_grid(ligand, dock_grid_info, &
                                     dock_grid%drugscore_grid, ff%ligdock(3),&
                                     g_tmp, calc_g)
    if (n_usc > 0) then
        call calc_DrugScore(f(1), g, is_usc, calc_g)
        if (calc_g) then
            do i_atm = 1, ligand%n_atm
                atm_no = ii_L(i_atm, ligand%lig_no)
                g(:,i_atm) = g(:,i_atm) + drugscore_w * g_all(:,atm_no,1)
            end do
        end if
        ff%ligdock(3) = ff%ligdock(3) + f(1)
    end if
    ff%ligdock(3) = drugscore_w * ff%ligdock(3)
    if (calc_g) then
        g(:,:) = g(:,:) + drugscore_w*g_tmp(:,:)
    end if
end if

if (use_Xscore) then
    call calc_HM_score_using_grid(ligand, dock_grid_info, &
                                  dock_grid%HM_grid, ff%ligdock(4), g_tmp, &
                                  calc_g)
    if (n_usc > 0) then
        call calc_HM_score(f(1), g_all(:,:,1), is_usc, calc_g)
        ff%ligdock(4) = ff%ligdock(4) + f(1)
        if (calc_g) then
            do i_atm = 1, ligand%n_atm
                atm_no = ii_L(i_atm, ligand%lig_no)
                g(:,i_atm) = g(:,i_atm) + X_HM_w * g_all(:,atm_no,1)
            end do
        end if
    end if
    ff%ligdock(4) = X_HM_w * ff%ligdock(4)
    if (calc_g) then
        g(:,:) = g(:,:) + X_HM_w*g_tmp(:,:)
    end if
end if

if (use_rsr) then
    call calc_ligdock_rsr_using_grid(ligand, dock_grid_info, &
                                     dock_grid%rsr_grid, ff%ligdock(7), &
                                     g_tmp, calc_g)
    ff%ligdock(7) = rsr_w * ff%ligdock(7)
    if (calc_g) then
        g(:,:) = g(:,:) + rsr_w*g_tmp(:,:)
    end if
end if

if (use_PLP) then
    call calc_PLP_tor(ff%ligdock(5), g_tmp, calc_g)
    ff%ligdock(5) = PLP_w * ff%ligdock(5)
    if (calc_g) then
        g(:,:) = g(:,:) + PLP_w*g_tmp(:,:)
    end if
end if

ff%ligdock(0) = sum(ff%ligdock(1:n_E_component_ligdock))
ff%total = ff%ligdock(0)

end subroutine ligdock_energy_using_grid
!-------------------------------------------------------------------------------
subroutine update_prot_prot_intrxn(rot_idx, ligdock_E)
!-------------------------------------------------------------------------------
integer, intent(in) :: rot_idx(:)
real(dp), intent(inout) :: ligdock_E(0:n_E_component_ligdock)

call calc_prot_E_using_table(rot_idx, ligdock_E(6))
ligdock_E(0) = ligdock_E(0) + ligdock_E(6)

end subroutine update_prot_prot_intrxn
!-------------------------------------------------------------------------------
subroutine ligdock_energy_and_gradient(ff, g, calc_g, ierr)
!-------------------------------------------------------------------------------
! TODO: 
!-------------------------------------------------------------------------------
type(energy_type), intent(inout) :: ff
real(dp), intent(inout) :: g(:,:)
logical, intent(in) :: calc_g
integer, intent(inout) :: ierr
real(dp) :: f(4), g_tmp(3, tn%atom)

ff%ligdock(:) = 0.0d0

if (use_atdk3) then
    if (atdk_w > small_real) then
        call calc_atdk3_inter_energy(f, g, appl_res, .false.)
        ff%ligdock(1:4) = atdk_w * f(1:4)
    end if
    if (atdk_int_w > small_real) then
        call calc_atdk3_intra_energy(f, g, .false.)
        ff%ligdock(5:8) = atdk_int_w * f(1:4)
    end if
else if (use_atdk4) then
    if (atdk_w > small_real) then
        call calc_atdk4_inter_energy(f, g, appl_res, .false.)
        ff%ligdock(1:4) = atdk_w * f(1:4)
    end if
    if (atdk_int_w > small_real) then
        call calc_atdk4_intra_energy(f, g, .false.)
        ff%ligdock(5:8) = atdk_int_w * f(1:4)
    end if
end if

if (use_drugscore) then
    call calc_DrugScore(ff%ligdock(9), g, appl_res, .false.)
    ff%ligdock(9) = drugscore_w * ff%ligdock(9)
end if

if (use_Xscore) then
    call calc_HM_score(ff%ligdock(10), g, appl_res, .false.)
    ff%ligdock(10) = X_HM_w * ff%ligdock(10)
end if

if (use_PLP) then
    call calc_PLP_tor(ff%ligdock(11), g_tmp, calc_g)
    ff%ligdock(11) = PLP_w * ff%ligdock(11)
end if

if (use_ROTA .and. (n_usc > 0)) then
    call calc_ROTA_score(ff%ligdock(12), g, appl_respair, .false.)
    ff%ligdock(12) = ROTA_w * ff%ligdock(12)
    call calc_atdk3_prot_vdw(ff%ligdock(13), g, appl_respair, .false.)
    ff%ligdock(13) = prot_w * ff%ligdock(13)
end if

ff%ligdock(0) = sum(ff%ligdock(1:n_E_component_ligdock))

end subroutine ligdock_energy_and_gradient
!-------------------------------------------------------------------------------
subroutine check_clash_within_ligand(ligand, clash)
!-------------------------------------------------------------------------------
! Check clash within ligand with distance threshold
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
integer, intent(out) :: clash

if (use_atdk3) then
    call check_clash_atdk3(ligand, clash)
else
    call check_clash_atdk4(ligand, clash)
end if

end subroutine check_clash_within_ligand
!-------------------------------------------------------------------------------
subroutine set_atdk_hbond_flex(usc_list, n_usc)
!-------------------------------------------------------------------------------
integer, intent(in) :: usc_list(:), n_usc

if (use_atdk3) then
    call set_atdk3_hbond_flex(usc_list, n_usc)
else if (use_atdk4) then
    call set_atdk4_hbond_flex(usc_list, n_usc)
end if

end subroutine set_atdk_hbond_flex 
!-------------------------------------------------------------------------------
subroutine set_Escale_ligdock(Esch)
!-------------------------------------------------------------------------------
type(energy_type), intent(in) :: Esch

end subroutine set_Escale_ligdock
!-------------------------------------------------------------------------------
subroutine receptor_gridization(FFT_grid_info)
!----------------------------------------------------------------------
type(docking_grid_param_type), intent(in) :: FFT_grid_info

r_grid%n_atdk = 0
r_grid%n_drugscore = 0
r_grid%n_rsr = 0

if (use_atdk3) then
    call FFT_rec_grid_atdk3(dock_grid_info, dock_grid, FFT_grid_info, r_grid)
else if (use_atdk4) then
    call FFT_rec_grid_atdk4(dock_grid_info, dock_grid, FFT_grid_info, r_grid)
end if

if (use_drugscore) then
    call FFT_rec_grid_drugscore(dock_grid_info, dock_grid, &
                                FFT_grid_info, r_grid)
end if

if (use_Xscore) then
    call FFT_rec_grid_HM(dock_grid_info, dock_grid, FFT_grid_info, r_grid)
end if

!if (use_rsr) then
!    call FFT_rec_grid_rsr(dock_grid_info, dock_grid, r_grid, frag_size)
!end if

end subroutine receptor_gridization
!-------------------------------------------------------------------------------
subroutine ligand_gridization(FFT_grid_info, ligand, fragment, n_frag_atm,&
                              R_frag, ordering)
!-------------------------------------------------------------------------------
type(docking_grid_param_type), intent(in) :: FFT_grid_info
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: fragment(:)
integer, intent(in) :: n_frag_atm
real(dp), intent(in) :: R_frag(:,:)
logical, intent(in) :: ordering

if (use_atdk3) then
    call FFT_lig_grid_atdk3(FFT_grid_info, dock_grid_info, ligand, &
                            fragment, n_frag_atm, R_frag, l_grid, ordering)
else if (use_atdk4) then
    call FFT_lig_grid_atdk4(FFT_grid_info, dock_grid_info, ligand, &
                            fragment, n_frag_atm, R_frag, l_grid, ordering)
end if

if (use_drugscore) then
    call FFT_lig_grid_drugscore(FFT_grid_info, dock_grid_info, ligand,&
                                fragment, n_frag_atm, R_frag, l_grid, ordering)
end if

if (use_Xscore) then
    call FFT_lig_grid_HM(FFT_grid_info, dock_grid_info, ligand, &
                         fragment, n_frag_atm, R_frag, l_grid, ordering)
end if

!if (use_rsr) then
!    call FFT_lig_grid_rsr(dock_grid_info, dock_grid, l_grid, ordering)
!end if

end subroutine ligand_gridization
!-------------------------------------------------------------------------------
subroutine finalize_FFT_grid()
!-------------------------------------------------------------------------------
if (use_atdk3 .or. use_atdk4) then
    deallocate(r_grid%atdk_grid)
    deallocate(l_grid%atdk_grid)
end if

if (use_drugscore) then
    deallocate(r_grid%drugscore_grid)
    deallocate(l_grid%drugscore_grid)
end if

if (use_Xscore) then
    deallocate(r_grid%HM_grid)
    deallocate(l_grid%HM_grid)
end if

end subroutine finalize_FFT_grid
!-------------------------------------------------------------------------------
END MODULE LIGDOCK_ENERGY
!-------------------------------------------------------------------------------
