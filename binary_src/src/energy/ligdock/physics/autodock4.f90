!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ligdock/physics/AutoDock4.f90
!
! Description: AutoDock4 energy terms 
!   This module contains subroutines for AutoDock4 energy calculation.
!      1. parameter setup
!       - initialize_atdk4_energy: read data file, assign atom types, etc.
!       - finalize_atdk4_energy  : deallocate variables
!      2. Energy calculation using grid
!       - construct_atdk4_grid          : construct grid for ligand docking
!       - calc_atdk4_intrxn_E_using_grid: calculate energy using grid
!      3. Continuous energy calculation
!       - calc_atdk4_inter_energy
!       
!   TODO: Add gradient calculation
!-------------------------------------------------------------------------------
MODULE AUTODOCK4
!-------------------------------------------------------------------------------

use globals
use logger
use string, only: parse_string
use mathfunctions, only: v_norm
use sort, only: sort1
!
use energy_vars 
use ligdock_E_vars
use ligdock_E_utils, only: copy_to_grid, calc_E_using_grid
!
use autodock_vars
use autodock_utils, only: get_dist_idx, get_charges, calc_ramp, &
                          racc_rdon_H, racc_rdon_N, racc_rdon_O, &
                          HD_hbond_parameter, OA_SA_hbond_parameter, &
                          NA_hbond_parameter, get_lig_atom_type
use autodock4_setup, only: read_ad4_param_file, construct_atdk4_para, &
                           setup_atdk4_E
use autodock4_atom_typer, only: assign_atom_types

implicit none
save
private

logical, allocatable :: is_Hdon(:), is_Hacc_N(:), is_Hacc_O(:)

public :: initialize_atdk4_energy
public :: finalize_atdk4_energy
public :: construct_atdk4_grid
public :: calc_atdk4_intrxn_E_using_grid
public :: calc_atdk4_intern_E_using_table
public :: calc_atdk4_inter_energy
public :: calc_atdk4_intra_energy
public :: check_clash_atdk4
public :: set_atdk4_hbond_flex
public :: FFT_rec_grid_atdk4
public :: FFT_lig_grid_atdk4

CONTAINS
!===============================================================================
! Initialize/Finalize autodock4 energy
!===============================================================================
subroutine initialize_atdk4_energy(protein, ligand, tna, &
                                   res_mol2_info, num_mol2_info)
!-------------------------------------------------------------------------------
! setup AutoDock4 parameters and energy tables
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: tna
type(residue_mol2_type), intent(in) :: res_mol2_info(:)
integer, intent(in) :: num_mol2_info
type(ad4_param_type) :: ad4_param
integer :: atdk_types(tna)

! Read parameter file
call read_ad4_param_file(infile_atdk_prm, ad4_param)

allocate(is_Hdon(tna))
allocate(is_Hacc_O(tna))
allocate(is_Hacc_N(tna))

! Do atom typing - see autodock4_atom_typer.f90
call assign_atom_types(protein, ligand, ad4_param, atdk_types,&
                       is_Hdon, is_Hacc_N, is_Hacc_O, &
                       res_mol2_info, num_mol2_info)

! store needed parameters to atdk_para type - see autodock4_setup.f90
allocate(atdk_para_idx(tna))
allocate(is_ligand(tna))
target_lig_no = tn%nonlig + ligand%lig_no
call construct_atdk4_para(ad4_param, atdk_types, ligand%lig_no, &
                          is_Hdon, is_Hacc_N, is_Hacc_O, is_ligand)

! Get atomic charges - see autodock_utils.f90
allocate(q_s(tna))
call get_charges(protein, q_s)

! setup energy table and other parameters to calculate autodock4 energy
! - see autodock4_setup.f90
allocate(qiqj(ligand%n_atm*(ligand%n_atm-1)/2))
allocate(sum_SiVj(ligand%n_atm*(ligand%n_atm-1)/2))
call setup_atdk4_E(ligand, is_Hdon, is_Hacc_N, is_Hacc_O)

end subroutine initialize_atdk4_energy
!-------------------------------------------------------------------------------
subroutine finalize_atdk4_energy()
!-------------------------------------------------------------------------------
deallocate(is_Hdon, is_Hacc_O, is_Hacc_N)
deallocate(atdk_para_idx)
deallocate(q_s)
deallocate(atdk_hbond)
deallocate(T_vdw)
deallocate(T_diel)
deallocate(T_solv)
deallocate(int_E_table)
deallocate(qiqj)
deallocate(sum_SiVj)

end subroutine finalize_atdk4_energy
!===============================================================================
! Construct AutoDock4 Grid
!===============================================================================
subroutine construct_atdk4_grid(grid_info, dock_grid, flex_res)
!-------------------------------------------------------------------------------
! Construct autodock grid (dock_grid%atdk_grid)
! atdk_grid(:,:,1) for coulomb interaction
! atdk_grid(:,:,2) for part of desolvation
! atdk_grid(:,:,3:n_atm_types) for vdw/hbond/part of desolvation by atom types
!-------------------------------------------------------------------------------
type(docking_grid_param_type), intent(in) :: grid_info
type(docking_grid_type), intent(inout) :: dock_grid
logical, intent(in) :: flex_res(:)
integer :: n_hash
real(dp), allocatable :: tmp_grid_E(:)
real(dp) :: ref_pt(3), pt(3)
integer :: i_lig, i, x, y, z
logical :: lig_Hacc_O

n_hash = grid_info%n_elem(1) * grid_info%n_elem(2) * grid_info%n_elem(3) 
! allocate autodock grid
allocate(dock_grid%atdk_grid(4, n_hash, atdk_para%n_atm_types(2)+2))
allocate(tmp_grid_E(atdk_para%n_atm_types(2)+2))

lig_Hacc_O = .false.
do i = 1, atdk_para%n_atm_types(2)
    i_lig = atdk_para%lig_atom_types(i)
    if (atdk_para%is_Hacc_O(i_lig)) then
        lig_Hacc_O = .true.
        exit
    end if
end do
    
do i = 1, 3
    ref_pt(i) = grid_info%grid_cntr(i) &
              - (dble(grid_info%n_elem(i))-1.0d0)/2.0d0 * grid_info%grid_width
end do

do z = 1, grid_info%n_elem(3)
    pt(3) = (dble(z)-1.0)*grid_info%grid_width + ref_pt(3)
    do y = 1, grid_info%n_elem(2)
        pt(2) = (dble(y) - 1.0)*grid_info%grid_width + ref_pt(2)
        do x = 1, grid_info%n_elem(1)
            pt(1) = (dble(x) - 1.0)*grid_info%grid_width + ref_pt(1)
            call calc_atdk4_grid(pt, tmp_grid_E, flex_res, lig_Hacc_O)
            do i_lig = 1, atdk_para%n_atm_types(2) + 2
                call copy_to_grid(tmp_grid_E(i_lig), grid_info%n_elem(1:3),&
                                  dock_grid%atdk_grid(:,:,i_lig), x, y, z, n_hash)
            end do
        end do
    end do
end do

deallocate(tmp_grid_E)

end subroutine construct_atdk4_grid
!-------------------------------------------------------------------------------
subroutine calc_atdk4_grid(grid_pt, E_pt, flex_res, lig_Hacc_O)
!-------------------------------------------------------------------------------
! calculate autodock interaction energy on grid point
!-------------------------------------------------------------------------------
real(dp), intent(in) :: grid_pt(:)
real(dp), intent(out) :: E_pt(:)
logical, intent(in) :: flex_res(:), lig_Hacc_O
!
integer :: n_atm_types, i_atm_type, i_type
integer :: i_atm, res_no, atm_no, idx_r_l, idx_r_s
real(dp) :: dr(3), dist
! for H-bond energy calculation
real(dp) :: hbond_max(atdk_para%n_atm_types(2))
real(dp) :: hbond_min(atdk_para%n_atm_types(2))
logical :: hbond_flag(atdk_para%n_atm_types(2))
integer :: close_atm
real(dp) :: tmp_E

n_atm_types = atdk_para%n_atm_types(2)
E_pt(:) = 0.0d0

! for H-bond term
hbond_max(:) = -99999.99d0
hbond_min(:) = 99999.99d0
hbond_flag(:) = .false.
close_atm = 0
if (lig_Hacc_O) call find_closest_atm_idx(grid_pt, close_atm)

do res_no = 1, tn%residue
    if (res_no == target_lig_no) cycle
    if (res_no <= tn%stdres) then
        if (flex_res(res_no)) cycle
    end if

    do atm_no = 1, res_index(res_no)%n_atm
        if (res_no > tn%nonlig) then
            i_atm = ii_L(atm_no, res_no-tn%nonlig)
        else
            i_atm = ii_R(atm_no, res_no)
        end if
        
        dr(1:3) = R(1:3,i_atm) - grid_pt(1:3)
        dist = sqrt(dot_product(dr,dr))
        ! get idx_r(index for T_vdw, T_diel, T_solv)
        call get_dist_idx(dist, idx_r_l, idx_r_s)

        ! calc coulomb interaction energy
        E_pt(1) = E_pt(1) + T_diel(idx_r_l) * q_s(i_atm) / dist

        if (dist > NBCutoff) cycle  ! for desolvation/vdw/Hbond, energy cutoff is 8A.
        i_atm_type = atdk_para_idx(i_atm)  ! atom type index to access atdk_para
       
        ! calc part of desolvation (independent of ligand atom types)
        E_pt(2) = E_pt(2) + solpar_q * atdk_para%vol(i_atm_type) * T_solv(idx_r_s) 
        
        ! calc vdW / Hbond / part of desolvation which depend on ligand atom types
        ! get parameters related to hbond energy calculation
        do i_type = 1, n_atm_types
            call calc_vdw_hbond_solv(i_atm, i_atm_type, i_type, idx_r_s, dr,&
                                     q_s(i_atm), tmp_E, close_atm, hbond_min(i_type), &
                                     hbond_max(i_type), hbond_flag(i_type))
            E_pt(i_type + 2) = E_pt(i_type + 2) + tmp_E
        end do
    end do
end do

do i_type = 1, n_atm_types
    if (hbond_flag(i_type)) then
        E_pt(i_type + 2) = E_pt(i_type + 2) + hbond_min(i_type)
        E_pt(i_type + 2) = E_pt(i_type + 2) + hbond_max(i_type)
    end if
end do

end subroutine calc_atdk4_grid
!-------------------------------------------------------------------------------
subroutine calc_vdw_hbond_solv(i_atm, prot_type, i_lig, idx_r, dr, q, pt_E, &
                               close_atm, h_min, h_max, h_flag)
!-------------------------------------------------------------------------------
! calculate vdw/Hbond, and part of desolvation energy of given protein-ligand
! atom pair
!-------------------------------------------------------------------------------
integer, intent(in) :: i_atm, prot_type, i_lig
integer, intent(in) :: idx_r, close_atm
real(dp), intent(in) :: q, dr(3)
real(dp), intent(out) :: pt_E
real(dp), intent(inout) :: h_min, h_max
logical, intent(inout) :: h_flag
real(dp) :: vdw, rsph, racc, rdon, Hramp, tmp_E, SiVj, SjVi, norm_dr(3)
integer :: lig_type

lig_type = atdk_para%lig_atom_types(i_lig)
pt_E = 0.0d0
vdw = T_vdw(idx_r, i_lig, prot_type)
if (atdk_para%hbond_read(lig_type, prot_type)) then ! H-bond
    rsph = vdw / 100.0d0
    rsph = max(rsph, 0.0d0)
    rsph = min(rsph, 1.0d0)
    norm_dr = dr
    call v_norm(norm_dr)
    if (atdk_para%is_Hdon(prot_type)) then
        call racc_rdon_H(atdk_hbond(i_atm), norm_dr)
        racc = atdk_hbond(i_atm)%racc
        if (atdk_para%is_Hacc_O(lig_type)) then ! prot: 'HD' ; lig: 'OA'
            call calc_ramp(i_atm, close_atm, atdk_hbond(i_atm), &
                           atdk_hbond(close_atm), Hramp)
            pt_E = pt_E + vdw*Hramp*(racc + (1.0d0-racc)*rsph)
        else if (atdk_para%is_Hacc_N(lig_type)) then ! prot: 'HD' ; lig: 'NA'
            tmp_E = vdw*(racc + (1.0d0-racc)*rsph)
            h_min = min(h_min, tmp_E)
            h_max = max(h_max, tmp_E)
            h_flag = .true.
        end if
    else if (atdk_para%is_Hacc_O(prot_type)) then ! prot: 'OA' ; lig: 'HD'
        call racc_rdon_O(atdk_hbond(i_atm), norm_dr)
        rdon = atdk_hbond(i_atm)%rdon
        tmp_E = vdw*(rdon + (1.0d0-rdon)*rsph)
        h_min = min(h_min, tmp_E)
        h_max = max(h_max, tmp_E)
        h_flag = .true.
    else if (atdk_para%is_Hacc_N(prot_type)) then ! prot: 'OA' ; lig: 'HD'
        call racc_rdon_N(atdk_hbond(i_atm), norm_dr)
        rdon = atdk_hbond(i_atm)%rdon
        tmp_E = vdw*(rdon + (1.0d0-rdon)*rsph)
        h_min = min(h_min, tmp_E)
        h_max = max(h_max, tmp_E)
        h_flag = .true.
    end if
else ! vdw
    pt_E = pt_E + vdw
end if

SiVj = atdk_para%solpar(lig_type) * atdk_para%vol(prot_type)
SjVi = (atdk_para%solpar(prot_type) + solpar_q*abs(q)) * atdk_para%vol(lig_type)
pt_E = pt_E + (SiVj + SjVi) * T_solv(idx_r)

end subroutine calc_vdw_hbond_solv
!===============================================================================
! Calculate AutoDock4 energy
!===============================================================================
subroutine calc_atdk4_intrxn_E_using_grid(ligand, grid_info, atdk_grid, &
                                          atdk_E, g, calc_g)
!-------------------------------------------------------------------------------
! calculate autodock interaction energy between protein and ligand using grid.
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
type(docking_grid_param_type), intent(in) :: grid_info
real(dp), intent(in) :: atdk_grid(:,:,:)
real(dp), intent(out) :: atdk_E, g(:,:)
logical, intent(in) :: calc_g
integer :: i_atm, atm_idx, lig_type, i
real(dp) :: tmp_E, elec_E, solv_E, vdw_E
real(dp) :: g_elec(3,ligand%n_atm), g_solv(3,ligand%n_atm), g_vdw(3,ligand%n_atm)

atdk_E = 0.0
elec_E = 0.0
solv_E = 0.0
vdw_E = 0.0
if (calc_g) g = 0.0d0

! coulombic interaction
do i_atm = 1, ligand%n_atm
    atm_idx = ii_L(i_atm, ligand%lig_no)
    call calc_E_using_grid(R(1:3, atm_idx), grid_info, atdk_grid(:,:,1),&
                           tmp_E, g_elec(:,i_atm), calc_g)
    elec_E = elec_E + tmp_E*ref_lig(ligand%lig_type)%charge(i_atm)
end do

! part of desolvation
do i_atm = 1, ligand%n_atm
    atm_idx = ii_L(i_atm, ligand%lig_no)
    call calc_E_using_grid(R(1:3, atm_idx), grid_info, atdk_grid(:,:,2),&
                           tmp_E, g_solv(:,i_atm), calc_g)
    solv_E = solv_E + tmp_E*abs(ref_lig(ligand%lig_type)%charge(i_atm))
end do

! vdw/hbond and rest of desolvation
do i_atm = 1, ligand%n_atm
    atm_idx = ii_L(i_atm, ligand%lig_no)
    call get_lig_atom_type(atdk_para_idx(atm_idx), lig_type)
    call calc_E_using_grid(R(1:3, atm_idx), grid_info, atdk_grid(:,:,2+lig_type),&
                           tmp_E, g_vdw(:,i_atm), calc_g)
    vdw_E = vdw_E + tmp_E!*atdk_vdw_scale
end do

atdk_E = elec_E + solv_E + vdw_E
if (calc_g) then
    do i_atm = 1, ligand%n_atm
        g(:,i_atm) = ref_lig(ligand%lig_type)%charge(i_atm) * g_elec(:,i_atm) &
                   + abs(ref_lig(ligand%lig_type)%charge(i_atm))*g_solv(:,i_atm)&
                   + g_vdw(:,i_atm)
    end do
end if

end subroutine calc_atdk4_intrxn_E_using_grid
!-------------------------------------------------------------------------------
subroutine calc_atdk4_intern_E_using_table(ligand, int_E, g, calc_g)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
real(dp), intent(out) :: int_E, g(:,:)
logical, intent(in) :: calc_g
!
integer :: i_atm, j_atm, atm_no_1, atm_no_2, type_1, type_2
integer :: idx_r, idx1, idx2, i_pair
real(dp) :: w_scale, elec_E, solv_E, vdw_E
real(dp) :: q1, q2, SiVj, SjVi
real(dp) :: dr(3), dist_sqr, dist, delta_r, delta_E, dEdr
integer :: i, j 

real(dp) :: R_tmp(3)

int_E = 0.0d0
if (calc_g) then
    g = 0.0d0
end if

i_pair = 0
do i_atm = 1, ligand%n_atm - 1
    atm_no_1 = ii_L(i_atm, ligand%lig_no)
    call get_lig_atom_type(atdk_para_idx(atm_no_1), type_1)
    !
    do j_atm = i_atm + 1, ligand%n_atm
        if (nb_matrix(j_atm, i_atm) == 0) then
            cycle
        end if
        i_pair = i_pair + 1
        !
        w_scale = 1.0d0
        if (nb_matrix(j_atm, i_atm) == 4) w_scale = scale_1_4

        atm_no_2 = ii_L(j_atm, ligand%lig_no)
        dr(1:3) = R(1:3, atm_no_2) - R(1:3, atm_no_1)
        dist_sqr = dot_product(dr, dr)
        
        if (dist_sqr > NBcutoff2) cycle

        idx_r = int(dist_sqr*32.0d0)
        idx_r = min(idx_r, max_dist_long)
        idx_r = max(1, idx_r)
        
        call get_lig_atom_type(atdk_para_idx(atm_no_2), type_2)
        ! 
        elec_E = int_E_table(1, idx_r, type_2, type_1)
        solv_E = int_E_table(2, idx_r, type_2, type_1)
        vdw_E = int_E_table(3, idx_r, type_2, type_1)
        !
        solv_E = sum_SiVj(i_pair)*solv_E
        elec_E = qiqj(i_pair) * elec_E
        int_E = int_E + w_scale*(elec_E + vdw_E + solv_E)

        if (calc_g) then
            idx1 = max(1,idx_r-1)
            idx2 = min(idx_r+1, max_dist_long)
            ! 
            dist = sqrt(dist_sqr)
            dr(:) = dr(:)/dist
            delta_r = sqrt(idx2/32.0d0) - sqrt(idx1/32.0)
            !
            ! elec_E
            delta_E = int_E_table(1, idx2, type_2, type_1) &
                    - int_E_table(1, idx1, type_2, type_1)
            dEdr = qiqj(i_pair)*(delta_E/delta_r)
            g(:,j_atm) = g(:,j_atm) + w_scale*dEdr*dr(:)
            g(:,i_atm) = g(:,i_atm) - w_scale*dEdr*dr(:)
            !
            ! solv_E
            delta_E = int_E_table(2, idx2, type_2, type_1) &
                    - int_E_table(2, idx1, type_2, type_1)
            dEdr =  sum_SiVj(i_pair)*(delta_E/delta_r)
            g(:,j_atm) = g(:,j_atm) + w_scale*dEdr*dr(:)
            g(:,i_atm) = g(:,i_atm) - w_scale*dEdr*dr(:)
            !
            ! vdw_E
            delta_E = int_E_table(3, idx2, type_2, type_1) &
                    - int_E_table(3, idx1, type_2, type_1)
            dEdr = delta_E/delta_r
            g(:,j_atm) = g(:,j_atm) + w_scale*dEdr*dr(:)
            g(:,i_atm) = g(:,i_atm) - w_scale*dEdr*dr(:)
        end if
    end do
end do

end subroutine calc_atdk4_intern_E_using_table
!-------------------------------------------------------------------------------
subroutine check_clash_atdk4(ligand, clash)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
integer, intent(out) :: clash
integer :: i_atm, j_atm, atm_no_1, atm_no_2, type_1, type_2
real(dp) :: dr(3), dist, vdw_sum

clash = 0

do i_atm = 1, ligand%n_atm - 1
    atm_no_1 = ii_L(i_atm, ligand%lig_no)
    type_1 = atdk_para_idx(atm_no_1)

    do j_atm = i_atm + 1, ligand%n_atm
        if (nb_matrix(j_atm, i_atm) == 0) cycle
        if (nb_matrix(j_atm, i_atm) == 4) cycle

        atm_no_2 = ii_L(j_atm, ligand%lig_no)
        type_2 = atdk_para_idx(atm_no_2)

        dr(1:3) = R(1:3, atm_no_2) - R(1:3, atm_no_1)
        dist = sqrt(dot_product(dr, dr))

        vdw_sum = atdk_para%ljrij(type_2, type_1)
        if (dist < 0.5d0*vdw_sum) clash = clash + 1
    end do
end do

end subroutine check_clash_atdk4
!-------------------------------------------------------------------------------
subroutine calc_atdk4_inter_energy(f, g, appl_res, calc_g)
!-------------------------------------------------------------------------------
! calculate autodock4 interaction energy between protein and ligand in
! continuous way.
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f(4), g(3, tn%atom)
logical, intent(in) :: appl_res(:), calc_g
!
integer :: i_atm, j_atm, close_atm, lig_atm_no, atm_no, n_atm
integer :: res_no
integer :: l_atm_type, p_atm_type
real(dp) :: l_q, p_q
real(dp) :: dr(3), dist_sqr, dist
real(dp) :: Hramp, hb_factor, tmp_E

f(:) = 0.0d0
g(:,:) = 0.0d0
do lig_atm_no = 1, res_index(target_lig_no)%n_atm
    i_atm = ii_L(lig_atm_no, target_lig_no-tn%nonlig)
    
    ! atdk_para index of ligand
    l_atm_type = atdk_para_idx(i_atm)
    l_q = q_s(i_atm)

    ! check current ligand atom is 'OA'
    ! if it's 'OA' then find closest 'HD' atom of protein
    close_atm = -1
    if (atdk_para%is_Hacc_O(l_atm_type)) then
        call find_closest_atm_idx(R(1:3, i_atm), close_atm)
    end if

    do res_no = 1, tn%residue   ! for protein/cofactor
        if (.not. appl_res(res_no)) cycle
        if (res_no == target_lig_no) cycle
        
        n_atm = res_index(res_no)%n_atm
        do atm_no = 1, n_atm
            if (res_no > tn%nonlig) then
                j_atm = ii_L(atm_no, res_no-tn%nonlig)
            else
                j_atm = ii_R(atm_no, res_no)
            end if

            p_atm_type = atdk_para_idx(j_atm)
            p_q = q_s(j_atm)

            dr(1:3) = R(1:3, j_atm) - R(1:3, i_atm)
            dist_sqr = dot_product(dr, dr)
            dist = sqrt(dist_sqr)
            
            ! atdk coulomb interaction
            call calc_atdk4_elec(dist, p_q, l_q, tmp_E)
            f(1) = f(1) + tmp_E

            if (dist > NBcutoff) cycle

            ! atdk desolvation
            call calc_atdk4_solv(dist_sqr, p_atm_type, l_atm_type,&
                                 p_q, l_q, tmp_E)
            f(2) = f(2) + tmp_E

            ! atdk vdw/hbond
            if (.not. atdk_para%hbond_read(l_atm_type, p_atm_type)) then ! vdW
                call calc_atdk4_vdw(dist, p_atm_type, l_atm_type, tmp_E, &
                                    atdk_vdw_w)
                f(3) = f(3) + tmp_E
            else ! Hbond
                dr(1:3) = dr(1:3)/dist
                call set_atdk4_Hbond_param(dr, j_atm, p_atm_type, close_atm, Hramp, &
                                           hb_factor)
                call calc_atdk4_Hbond(dist, p_atm_type, l_atm_type, &
                                      Hramp, hb_factor, tmp_E, atdk_hbond_w)
                f(4) = f(4) + tmp_E
            end if
        end do
    end do
end do

f(1) = f(1) * atdk_elec_w
f(2) = f(2) * atdk_solv_w

end subroutine calc_atdk4_inter_energy
!-------------------------------------------------------------------------------
subroutine calc_atdk4_intra_energy(f, g, calc_g)
!-------------------------------------------------------------------------------
! calculate autodock4 interaction energy between protein and ligand in
! continuous way.
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f(4), g(3, tn%atom)
logical, intent(in) :: calc_g
!
integer :: i_atm, j_atm, atm_1, atm_2
integer :: atm_type_1, atm_type_2
real(dp) :: q_1, q_2
real(dp) :: dr(3), dist_sqr, dist
real(dp) :: Hramp, hb_factor, tmp_E, w_scale

f(:) = 0.0d0
g(:,:) = 0.0d0
do atm_1 = 1, res_index(target_lig_no)%n_atm - 1
    i_atm = ii_L(atm_1, target_lig_no-tn%nonlig)
    ! atdk_para index of ligand
    atm_type_1 = atdk_para_idx(i_atm)
    q_1 = q_s(i_atm)

    do atm_2 = atm_1 + 1, res_index(target_lig_no)%n_atm
        if (nb_matrix(atm_2, atm_1) == 0) cycle
        w_scale = 1.0d0
        if (nb_matrix(atm_2, atm_1) == 4) w_scale = scale_1_4 
        
        j_atm = ii_L(atm_2, target_lig_no-tn%nonlig)

        atm_type_2 = atdk_para_idx(j_atm)
        q_2 = q_s(j_atm)

        dr(1:3) = R(1:3, j_atm) - R(1:3, i_atm)
        dist_sqr = dot_product(dr, dr)
        dist = sqrt(dist_sqr)
        
        ! atdk coulomb interaction
        call calc_atdk4_elec(dist, q_2, q_1, tmp_E)
        f(1) = f(1) + w_scale * tmp_E
        
        if (dist > NBcutoff) cycle

        ! atdk desolvation
        call calc_atdk4_solv(dist_sqr, atm_type_2, atm_type_1, &
                             q_2, q_1, tmp_E)
        f(2) = f(2) + w_scale * tmp_E

        ! atdk vdw/hbond
        if (.not. atdk_para%hbond_read(atm_type_1, atm_type_2)) then ! vdW
            call calc_atdk4_vdw(dist, atm_type_2, atm_type_1, tmp_E, &
                                atdk_int_v_w)
            f(3) = f(3) + w_scale * tmp_E
        else ! Hbond
            Hramp = 1.0d0
            hb_factor = 1.0d0
            call calc_atdk4_Hbond(dist, atm_type_2, atm_type_1, &
                                  Hramp, hb_factor, tmp_E, atdk_int_h_w)
            f(4) = f(4) + w_scale * tmp_E
        end if
    end do
end do

f(1) = f(1) * atdk_int_e_w
f(2) = f(2) * atdk_int_s_w

end subroutine calc_atdk4_intra_energy
!-------------------------------------------------------------------------------
subroutine calc_atdk4_elec(dist, q1, q2, E)
!-------------------------------------------------------------------------------
! calculate autodock4 coulomb interaction energy
! Eqq = 332.0 * q1 * q2 / (eps(dist) * dist)
! eps(dist) = A + B/(1.0+rk*exp(lamB*dist))
!-------------------------------------------------------------------------------
real(dp), intent(in) :: dist, q1, q2
real(dp), intent(out) :: E
real(dp), parameter :: epsij0 = 78.4d0
real(dp), parameter :: lam = 0.003627d0
real(dp), parameter :: A = -8.5525d0
real(dp), parameter :: B = epsij0 - A
real(dp), parameter :: rk = 7.7839d0
real(dp), parameter :: lamB = -lam * B
real(dp) :: eps

E = 0.0d0
eps = A + B/(1.0d0 + rk*exp(lamB*dist))
E = 332.0d0*q1*q2 / (dist*eps)

end subroutine calc_atdk4_elec
!-------------------------------------------------------------------------------
subroutine calc_atdk4_solv(dist_sqr, p_type, l_type, q_p, q_l, E)
!-------------------------------------------------------------------------------
! calculate autodock4 desolvation energy
! Edesolv = (SiVj + SjVi) * exp(const*(dist**2))
! sol_const = -1.0/(2.0*(sol_sig**2)) in autodock_vars.f90
! Si = solpar(i_type) + solpar_q * abs(qi) | Vi = vol(i_type)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: dist_sqr, q_p, q_l
integer, intent(in) :: p_type, l_type
real(dp), intent(out) :: E
real(dp) :: SiVj, SjVi

E = 0.0d0

SiVj = (atdk_para%solpar(l_type) + solpar_q*abs(q_l)) * atdk_para%vol(p_type)
SjVi = (atdk_para%solpar(p_type) + solpar_q*abs(q_p)) * atdk_para%vol(l_type)
E = (SiVj + SjVi) * exp(dist_sqr*sol_const)

end subroutine calc_atdk4_solv
!-------------------------------------------------------------------------------
subroutine calc_atdk4_vdw(dist, p_type, l_type, E, vdw_w)
!-------------------------------------------------------------------------------
! calculate autodock4 van der Waals interaction energy
! It's a flat-bottom 12-6 LJ potential 
!-------------------------------------------------------------------------------
real(dp), intent(in) :: dist, vdw_w
integer, intent(in) :: p_type, l_type
real(dp), intent(out) :: E
real(dp) :: expA, expB, Rij, Eij, cA, cB, rA, rB
real(dp) :: Rij_12, Rij_6
real(dp) :: sm_dist

E = 0.0d0

expA = 12.0d0
expB = 6.0d0

Rij = atdk_para%ljrij(l_type, p_type)
Eij = atdk_para%ljeij(l_type, p_type) * vdw_w

Rij_6 = Rij**6
Rij_12 = Rij_6**2

cA = Eij/(expA-expB) * (Rij_12)*expB
cB = Eij/(expA-expB) * (Rij_6)*expA

sm_dist = dist
call smooth_dist(sm_dist, Rij)

rB = sm_dist**expB
rA = rB**2

E = cA/rA - cB/rB

end subroutine calc_atdk4_vdw
!-------------------------------------------------------------------------------
subroutine calc_atdk4_Hbond(dist, p_type, l_type, Hramp, hb_factor, E, hbond_w)
!-------------------------------------------------------------------------------
! calculate autodock4 Hbond interaction energy
! It's a flat-bottom 12-10 LJ potential with consideration of diretionality. 
!-------------------------------------------------------------------------------
real(dp), intent(in) :: dist, hbond_w
integer, intent(in) :: p_type, l_type
real(dp), intent(in) :: Hramp, hb_factor
real(dp), intent(out) :: E
!
real(dp) :: expA, expB, Rij, Eij, cA, cB, rA, rB
real(dp) :: sm_dist
!
real(dp) :: rsph

E = 0.0d0

expA = 12.0d0
expB = 10.0d0

Rij = atdk_para%hbrij(l_type, p_type)
Eij = atdk_para%hbeij(l_type, p_type) * hbond_w

cA = Eij/(expA-expB) * (Rij**expA)*expB
cB = Eij/(expA-expB) * (Rij**expB)*expA

sm_dist = dist
call smooth_dist(sm_dist, Rij)

rA = sm_dist**expA
rB = sm_dist**expB

E = cA/rA - cB/rB

rsph = E/100.0d0
rsph = max(rsph, 0.0)
rsph = min(rsph, 1.0)

E = Hramp * E * (hb_factor + (1.0d0 - hb_factor)*rsph)

end subroutine calc_atdk4_Hbond
!-------------------------------------------------------------------------------
subroutine set_atdk4_Hbond_param(dr, atm_idx, atm_type, close_atm, Hramp, hb_factor)
!-------------------------------------------------------------------------------
! to calculate autodock4 Hbond energy, get Hramp, racc/rdon value.
!-------------------------------------------------------------------------------
real(dp), intent(in) :: dr(3)
integer, intent(in) :: atm_idx, atm_type, close_atm
real(dp), intent(out) :: Hramp, hb_factor

Hramp = 1.0
hb_factor = 1.0

if (atdk_para%is_Hdon(atm_type)) then ! current protein atom is 'HD'
    call racc_rdon_H(atdk_hbond(atm_idx), dr)
    hb_factor = atdk_hbond(atm_idx)%racc
    if (close_atm /= -1) then
        call calc_ramp(atm_idx, close_atm, atdk_hbond(atm_idx), &
                       atdk_hbond(close_atm), Hramp)
    end if

else if (atdk_para%is_Hacc_N(atm_type)) then ! curr protein atom is 'NA'
    call racc_rdon_N(atdk_hbond(atm_idx), dr)
    hb_factor = atdk_hbond(atm_idx)%rdon

else if (atdk_para%is_Hacc_O(atm_type)) then ! curr protein atom is 'OA'
    call racc_rdon_O(atdk_hbond(atm_idx), dr)
    hb_factor = atdk_hbond(atm_idx)%rdon
end if

end subroutine set_atdk4_Hbond_param
!-------------------------------------------------------------------------------
subroutine smooth_dist(dist, Rij)
!-------------------------------------------------------------------------------
! make vdW/Hbond potential as flat-bottom potential
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: dist
real(dp), intent(in) :: Rij
real(dp), parameter :: smooth_width = 0.25

if (dist < Rij - smooth_width) then
    dist = dist + smooth_width
else if (dist > Rij + smooth_width) then
    dist = dist - smooth_width
else
    dist = Rij
end if

end subroutine smooth_dist
!===============================================================================
! Misc.
!===============================================================================
subroutine find_closest_atm_idx(pt, close_atm)
!-------------------------------------------------------------------------------
! find closest hydrogen donor atom of protein/cofactor to current point.
!-------------------------------------------------------------------------------
real(dp), intent(in) :: pt(3)
integer, intent(out) :: close_atm
integer :: i_atm, atm_type
real(dp) :: dr(3), dist_sqr, min_dist_sqr

min_dist_sqr = 99999.99d0
do i_atm = 1, tn%atom
    if (is_ligand(i_atm)) cycle
    atm_type = atdk_para_idx(i_atm)
    if (.not. atdk_para%is_Hdon(atm_type)) cycle
    dr(1:3) = pt(1:3) - R(1:3, i_atm)
    dist_sqr = dot_product(dr(1:3), dr(1:3))
    if (dist_sqr < min_dist_sqr) then
        min_dist_sqr = dist_sqr
        close_atm = i_atm
    end if
end do

end subroutine find_closest_atm_idx
!-------------------------------------------------------------------------------
subroutine set_atdk4_hbond_flex(usc_list, n_usc)
!-------------------------------------------------------------------------------
integer, intent(in) :: usc_list(:), n_usc
integer :: i_usc, i_res, i_atm, atm_no

do i_usc = 1, n_usc
    i_res = usc_list(i_usc)
    do i_atm = 1, res_index(i_res)%n_atm
        atm_no = ii_R(i_atm, i_res)
        if (is_Hdon(atm_no)) then
            call HD_hbond_parameter(atdk_hbond(atm_no), atm_no)
        else if (is_Hacc_N(atm_no)) then
            call NA_hbond_parameter(atdk_hbond(atm_no), atm_no)
        else if (is_Hacc_O(atm_no)) then
            call OA_SA_hbond_parameter(atdk_hbond(atm_no), atm_no)
        end if
    end do
end do

end subroutine set_atdk4_hbond_flex
!-------------------------------------------------------------------------------
subroutine FFT_rec_grid_atdk4(grid_info, dock_grid, FFT_grid_info, r_grid)
!-------------------------------------------------------------------------------
type(docking_grid_param_type), intent(in) :: grid_info, FFT_grid_info
type(docking_grid_type), intent(in) :: dock_grid
type(docking_FFT_grid_type), intent(inout) :: r_grid
!
integer :: n_type, i_type, j_type
integer :: x, y, z, hash, FFT_hash
real(dp) :: grid_v

if (mod(atdk_para%n_atm_types(2),2) == 0) then
    n_type = atdk_para%n_atm_types(2)/2 + 1
else
    n_type = atdk_para%n_atm_types(2)/2 + 2
end if
r_grid%n_atdk = n_type

FFT_hash = FFT_grid_info%n_elem(1)*FFT_grid_info%n_elem(2)*FFT_grid_info%n_elem(3)

allocate(r_grid%atdk_grid(FFT_hash, n_type))
r_grid%atdk_grid = cmplx(0.0,0.0)

do i_type = 1, atdk_para%n_atm_types(2)+2
    do z = 1, FFT_grid_info%n_elem(3)
        do y = 1, FFT_grid_info%n_elem(2)
            do x = 1, FFT_grid_info%n_elem(1)
                !
                FFT_hash = FFT_grid_info%n_elem(2)*FFT_grid_info%n_elem(1)*(z-1) &
                         + FFT_grid_info%n_elem(1)*(y-1) + x
                !
                if ((x > grid_info%n_elem(1)) .or. &
                    (y > grid_info%n_elem(2)) .or. &
                    (z > grid_info%n_elem(3))) then  ! outside of docking grid
                    if (mod(i_type,2) == 1) then
                        j_type = int(i_type/2) + 1
                        r_grid%atdk_grid(FFT_hash,j_type) = cmplx(max_energy, max_energy)
                    end if
                    cycle
                end if
                !
                hash = grid_info%n_elem(2)*grid_info%n_elem(1)*(z-1) &
                     + grid_info%n_elem(1)*(y-1) + x
                grid_v = dock_grid%atdk_grid(1,hash,i_type)
                !
                if (mod(i_type,2) == 1) then
                    j_type = int(i_type/2) + 1
                    r_grid%atdk_grid(FFT_hash,j_type) = &
                    r_grid%atdk_grid(FFT_hash,j_type) + cmplx(grid_v, 0.0)
                else
                    j_type = int(i_type/2)
                    r_grid%atdk_grid(FFT_hash,j_type) = &
                    r_grid%atdk_grid(FFT_hash,j_type) + cmplx(0.0,grid_v)
                end if
            end do
        end do
    end do
end do

end subroutine FFT_rec_grid_atdk4
!-------------------------------------------------------------------------------
subroutine FFT_lig_grid_atdk4(FFT_grid_info, grid_info, ligand,&
                              fragment, n_frag_atm, R_frag, l_grid, ordering)
!-------------------------------------------------------------------------------
type(docking_grid_param_type), intent(in) :: FFT_grid_info, grid_info
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: fragment(:), n_frag_atm
real(dp), intent(in) :: R_frag(:,:)
type(docking_FFT_grid_type), intent(inout) :: l_grid
logical, intent(in) :: ordering
integer :: n_type, i_type, j_type, lig_type
integer :: i_atm, atm_no
integer :: x, y, z, X0(3), FFT_hash
real(dp) :: dr(3), tcrd(3)
!
integer, allocatable, save :: order_by_type(:), type_s(:)

if (.not. allocated(l_grid%atdk_grid)) then
    if (mod(atdk_para%n_atm_types(2),2) == 0) then
        n_type = atdk_para%n_atm_types(2)/2 + 1
    else
        n_type = atdk_para%n_atm_types(2)/2 + 2
    end if
    l_grid%n_atdk = n_type
    
    FFT_hash = FFT_grid_info%n_elem(1)*FFT_grid_info%n_elem(2)*FFT_grid_info%n_elem(3)

    allocate(l_grid%atdk_grid(FFT_hash, n_type))
end if

if (ordering) then
    if (allocated(order_by_type)) then
        deallocate(type_s)
        deallocate(order_by_type)
    end if
    allocate(type_s(n_frag_atm))
    allocate(order_by_type(n_frag_atm))
    do i_atm = 1, n_frag_atm
        atm_no = ii_L(fragment(i_atm), ligand%lig_no)
        call get_lig_atom_type(atdk_para_idx(atm_no), i_type)
        type_s(i_atm) = i_type
    end do
    call sort1(n_frag_atm, type_s, order_by_type)
    return
end if

l_grid%atdk_grid = cmplx(0.0,0.0)

! first grid for coulomb interaction: Lp(i,j,k) = q or 0
! 2nd grid for desolvation: Lp(i,j,k) = |q| or 0
do i_atm = 1, n_frag_atm
    atm_no = fragment(i_atm)
    !
    dr(:) = R_frag(:,i_atm) - FFT_grid_info%grid_cntr(:)
    tcrd(:) = dr(:)/FFT_grid_info%grid_width
    tcrd(:) = tcrd(:) + (dble(grid_info%n_elem(:)-1)*0.5d0 + 1.0d0)
    X0(:) = nint(tcrd(:)) ! nearest grid point
    !
    x = X0(1)
    y = X0(2)
    z = X0(3)
    FFT_hash = FFT_grid_info%n_elem(2)*FFT_grid_info%n_elem(1)*(z-1) &
             + FFT_grid_info%n_elem(1)*(y-1) + x
    !
    l_grid%atdk_grid(FFT_hash,1) = l_grid%atdk_grid(FFT_hash,1) + &
                                cmplx(ref_lig(ligand%lig_type)%charge(atm_no),&
                                      abs(ref_lig(ligand%lig_type)%charge(atm_no)))
end do

! others for vdw/Hbond/desolv: Lp(i,j,k) = 1 or 0
do i_atm = 1, n_frag_atm
    atm_no = ii_L(fragment(i_atm), ligand%lig_no)
    call get_lig_atom_type(atdk_para_idx(atm_no), i_type)
    !
    dr(:) = R_frag(:,i_atm) - FFT_grid_info%grid_cntr(:)
    tcrd(:) = dr(:)/FFT_grid_info%grid_width
    tcrd(:) = tcrd(:) + (dble(grid_info%n_elem(:)-1)*0.5d0 + 1.0d0)
    X0(:) = nint(tcrd(:)) ! nearest grid point
    !
    x = X0(1)
    y = X0(2)
    z = X0(3)
    FFT_hash = FFT_grid_info%n_elem(2)*FFT_grid_info%n_elem(1)*(z-1) &
             + FFT_grid_info%n_elem(1)*(y-1) + x
    !
    if (mod(i_type,2) == 1) then
        j_type = int(i_type/2) + 2
        l_grid%atdk_grid(FFT_hash,j_type) = l_grid%atdk_grid(FFT_hash,j_type) + & 
                                         cmplx(1.0, 0.0)
    else
        j_type = int(i_type/2) + 1
        l_grid%atdk_grid(FFT_hash,j_type) = l_grid%atdk_grid(FFT_hash,j_type) + &
                                         cmplx(0.0, 1.0)
    end if
end do

end subroutine FFT_lig_grid_atdk4
!-------------------------------------------------------------------------------
END MODULE AUTODOCK4
