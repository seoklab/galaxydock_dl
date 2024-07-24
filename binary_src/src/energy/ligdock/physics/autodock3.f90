!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ligdock/physics/autodock3.f90
!
! Description: AutoDock3 energy terms 
!       
!-------------------------------------------------------------------------------
MODULE AUTODOCK3
!-------------------------------------------------------------------------------

use globals
!
use sort, only: sort1

use energy_vars
use ligdock_E_utils, only: copy_to_grid, calc_E_using_grid
!
use autodock_vars
use autodock_utils, only: get_dist_idx, get_charges, racc_rdon_H, racc_rdon_O, &
                          HD_hbond_parameter, OA_SA_hbond_parameter, &
                          get_lig_atom_type
!
use autodock3_setup, only: construct_atdk3_para, read_assign_atdk3_para, &
                           doping_prot_sol_par, setup_atdk3_E

implicit none
save
private

real(dp), allocatable :: atdk3_sol_V(:)
logical, allocatable :: is_Hdon(:), is_Hacc_O(:)

public :: initialize_atdk3_energy
public :: finalize_atdk3_energy
public :: construct_atdk3_grid
public :: calc_atdk3_intrxn_E_using_grid
public :: calc_atdk3_intern_E_using_table
public :: calc_atdk3_inter_energy
public :: calc_atdk3_intra_energy
public :: calc_atdk3_prot_vdw
public :: check_clash_atdk3
public :: set_atdk3_hbond_flex
public :: calc_atdk_rsr_weight
public :: FFT_rec_grid_atdk3
public :: FFT_lig_grid_atdk3

CONTAINS
!===============================================================================
! Initialize/Finalize autodock3 energy
!===============================================================================
subroutine initialize_atdk3_energy(protein, ligand)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand

allocate(atdk_para_idx(tn%atom))
allocate(is_ligand(tn%atom))
allocate(is_Hdon(tn%atom), is_Hacc_O(tn%atom))
! Check atom types in protein/cofactor/ligand first.
target_lig_no = ligand%lig_no + tn%nonlig
call construct_atdk3_para(protein, ligand, is_Hdon, is_Hacc_O)

! Read parameter file, and get parameters of used atom types
call read_assign_atdk3_para(infile_atdk_prm, atdk_para)
allocate(atdk3_sol_V(tn%atom))
call doping_prot_sol_par(infile_atdk3_sol_prm, protein, atdk3_sol_V)

! Get atomic charges - see autodock_utils.f90
allocate(q_s(tn%atom))
call get_charges(protein, q_s)

! setup energy table and other parameters to calculate autodock3 energy
call setup_atdk3_E(ligand, is_Hdon, is_Hacc_O)

end subroutine initialize_atdk3_energy
!-------------------------------------------------------------------------------
subroutine finalize_atdk3_energy()
!-------------------------------------------------------------------------------
deallocate(atdk_para_idx)
deallocate(is_ligand)
deallocate(is_Hdon, is_Hacc_O)
deallocate(atdk3_sol_V)
deallocate(q_s)
deallocate(atdk_hbond)
deallocate(T_vdw)
deallocate(T_diel)
deallocate(T_solv)
deallocate(int_E_table)

end subroutine finalize_atdk3_energy
!===============================================================================
! Construct AutoDock3 Grid
!===============================================================================
subroutine construct_atdk3_grid(grid_info, dock_grid, flex_res)
!-------------------------------------------------------------------------------
! Construct autodock grid (dock_grid%atdk_grid)
! atdk_grid(:,:,1) for coulomb interaction
! atdk_grid(:,:,2:n_atm_types+1) for vdw/hbond/part of desolvation by atom types
!-------------------------------------------------------------------------------
type(docking_grid_param_type), intent(in) :: grid_info
type(docking_grid_type), intent(inout) :: dock_grid
logical, intent(in) :: flex_res(:)
integer :: n_hash
real(dp), allocatable :: tmp_grid_E(:)
real(dp) :: ref_pt(3), grid_pt(3)
integer :: i_lig, i, x, y, z

n_hash = grid_info%n_elem(1) * grid_info%n_elem(2) * grid_info%n_elem(3) 
! allocate autodock grid
allocate(dock_grid%atdk_grid(4, n_hash, atdk_para%n_atm_types(2)+1))
allocate(tmp_grid_E(atdk_para%n_atm_types(2)+1))

tmp_grid_E(:) = 0.0d0

do i = 1, 3
    ref_pt(i) = grid_info%grid_cntr(i) &
              - (dble(grid_info%n_elem(i))-1.0d0)/2.0d0 * grid_info%grid_width
end do

do z = 1, grid_info%n_elem(3)
    do y = 1, grid_info%n_elem(2)
        do x = 1, grid_info%n_elem(1)
            grid_pt(:) = (dble((/x, y, z/))-1.0d0)*grid_info%grid_width + ref_pt(:)
            call calc_atdk3_grid(grid_pt(1:3), tmp_grid_E, flex_res)
            do i_lig = 1, atdk_para%n_atm_types(2) + 1
                call copy_to_grid(tmp_grid_E(i_lig), grid_info%n_elem(1:3),&
                                  dock_grid%atdk_grid(:,:,i_lig), x, y, z, n_hash)
            end do
        end do
    end do
end do
deallocate(tmp_grid_E)

end subroutine construct_atdk3_grid
!-------------------------------------------------------------------------------
subroutine calc_atdk3_grid(grid_pt, E_pt, flex_res)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: grid_pt(3)
real(dp), intent(inout) :: E_pt(:)
logical, intent(in) :: flex_res(:)
!
integer :: n_atm_types, i_lig, i_prot, lig_type
integer :: i_atm, res_no, atm_no
integer :: idx_r_l, idx_r_s
real(dp) :: dr(3), dist
real(dp) :: tmp_E, rsph, rcon

n_atm_types = atdk_para%n_atm_types(2)
E_pt(1) = 0.0d0
do i_lig = 1, n_atm_types
    lig_type = atdk_para%lig_atom_types(i_lig)
    E_pt(i_lig+1) = atdk_para%solcn(lig_type)
end do

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

        dr(1:3) = R(1:3, i_atm) - grid_pt(1:3)
        dist = sqrt(dot_product(dr,dr))
        call get_dist_idx(dist, idx_r_l, idx_r_s)

        E_pt(1) = E_pt(1) + T_diel(idx_r_l)*q_s(i_atm)/dist

        if (dist > NBCutoff) cycle
        i_prot = atdk_para_idx(i_atm)

        rcon = 1.0d0
        if (atdk_para%type_name(i_prot) == 'H ') then
            dr(1:3) = dr(1:3)/dist
            call racc_rdon_H(atdk_hbond(i_atm), dr)
            rcon = atdk_hbond(i_atm)%racc
        else if (atdk_para%type_name(i_prot) == 'O ') then
            dr(1:3) = dr(1:3)/dist
            call racc_rdon_O(atdk_hbond(i_atm), dr)
            rcon = atdk_hbond(i_atm)%rdon
        end if

        do lig_type = 1, n_atm_types
            i_lig = atdk_para%lig_atom_types(lig_type)
            
            tmp_E = T_vdw(idx_r_s, lig_type, i_prot)
            if (atdk_para%hbond_read(i_lig, i_prot)) then
                rsph = tmp_E/100.0d0
                rsph = max(0.0d0, rsph)
                rsph = min(1.0d0, rsph)
                tmp_E = tmp_E*(rcon + (1.0d0-rcon)*rsph)
            end if

            E_pt(lig_type+1) = E_pt(lig_type+1) + tmp_E

            tmp_E = T_solv(idx_r_s) * atdk_para%solpar(i_lig) * &
                    atdk3_sol_V(i_atm)
            E_pt(lig_type+1) = E_pt(lig_type+1) - tmp_E
        end do
    end do
end do

end subroutine calc_atdk3_grid
!===============================================================================
! Calculate AutoDock3 energy
!===============================================================================
subroutine calc_atdk3_intrxn_E_using_grid(ligand, grid_info, atdk_grid, &
                                          atdk_E, g, calc_g)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
type(docking_grid_param_type), intent(in) :: grid_info
real(dp), intent(in) :: atdk_grid(:,:,:)
real(dp), intent(out) :: atdk_E, g(:,:)
logical, intent(in) :: calc_g
integer :: i_atm, atm_idx, lig_type
real(dp) :: tmp_E, elec_E, vdw_E, g_elec(3,ligand%n_atm), g_vdw(3,ligand%n_atm)

atdk_E = 0.0
elec_E = 0.0
vdw_E = 0.0
if (calc_g) g(:,:) = 0.0d0

! coulombic interaction
do i_atm = 1, ligand%n_atm
    atm_idx = ii_L(i_atm, ligand%lig_no)
    call calc_E_using_grid(R(1:3, atm_idx), grid_info, atdk_grid(:,:,1),&
                           tmp_E, g_elec(:,i_atm), calc_g)
    elec_E = elec_E + tmp_E*ref_lig(ligand%lig_type)%charge(i_atm)
end do

! vdw/hbond and desolvation
do i_atm = 1, ligand%n_atm
    atm_idx = ii_L(i_atm, ligand%lig_no)
    call get_lig_atom_type(atdk_para_idx(atm_idx), lig_type)
    call calc_E_using_grid(R(1:3, atm_idx), grid_info, atdk_grid(:,:,1+lig_type),&
                           tmp_E, g_vdw(:,i_atm), calc_g)
    vdw_E = vdw_E + tmp_E
end do

if (calc_g) then
    do i_atm = 1, ligand%n_atm
        g(:,i_atm) = ref_lig(ligand%lig_type)%charge(i_atm)*g_elec(:,i_atm) &
                   + g_vdw(:,i_atm)
    end do
end if

atdk_E = elec_E + vdw_E

end subroutine calc_atdk3_intrxn_E_using_grid
!-------------------------------------------------------------------------------
subroutine calc_atdk3_intern_E_using_table(ligand, int_E, g, calc_g)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
real(dp), intent(out) :: int_E, g(:,:)
logical, intent(in) :: calc_g
!
integer :: i_atm, j_atm, atm_no_1, atm_no_2, type_1, type_2
integer :: idx_r, idx1, idx2
real(dp) :: vdw_E
real(dp) :: dr(3), dist_sqr, dist, delta_r, delta_E, dEdr

int_E = 0.0d0
vdw_E = 0.0d0
if (calc_g) g = 0.0d0

do i_atm = 1, ligand%n_atm - 1
    atm_no_1 = ii_L(i_atm, ligand%lig_no)
    call get_lig_atom_type(atdk_para_idx(atm_no_1), type_1)

    do j_atm = i_atm + 1, ligand%n_atm
        if (nb_matrix(j_atm, i_atm) == 0) cycle
        if (nb_matrix(j_atm, i_atm) == 4) cycle

        atm_no_2 = ii_L(j_atm, ligand%lig_no)
        call get_lig_atom_type(atdk_para_idx(atm_no_2), type_2)

        dr(1:3) = R(1:3, atm_no_2) - R(1:3, atm_no_1)
        dist_sqr = dot_product(dr, dr)
        if (dist_sqr > NBcutoff2) cycle

        idx_r = int(dist_sqr*32.0d0)
        idx_r = min(idx_r, max_dist_long)
        idx_r = max(1, idx_r)

        vdw_E = int_E_table(3, idx_r, type_2, type_1)
        int_E = int_E + vdw_E

        if (calc_g) then
            idx1 = max(1,idx_r-1)
            idx2 = min(idx_r+1, max_dist_long)
            delta_r = sqrt(idx2/32.0d0) - sqrt(idx1/32.0)
            delta_E = int_E_table(3, idx2, type_2, type_1) &
                    - int_E_table(3, idx1, type_2, type_1)
            dEdr = delta_E/delta_r
            !
            dist = sqrt(dist_sqr)
            dr(:) = dr(:)/dist
            !
            g(:,j_atm) = g(:,j_atm) + dEdr*dr(:)
            g(:,i_atm) = g(:,i_atm) - dEdr*dr(:)
        end if
    end do
end do

end subroutine calc_atdk3_intern_E_using_table
!-------------------------------------------------------------------------------
subroutine check_clash_atdk3(ligand, clash)
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

end subroutine check_clash_atdk3
!-------------------------------------------------------------------------------
subroutine calc_atdk3_inter_energy(f, g, appl_res, calc_g)
!-------------------------------------------------------------------------------
! calculate autodock3 interaction energy between protein and ligand in
! continuous way.
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f(4), g(3, tn%atom)
logical, intent(in) :: appl_res(:), calc_g
!
integer :: i_atm, j_atm, n_atm
integer :: res_no, atm_no, lig_atm_no
integer :: l_atm_type, p_atm_type
real(dp) :: l_q, p_q
real(dp) :: dr(3), dist_sqr, dist
real(dp) :: hb_factor, tmp_E

f(:) = 0.0d0
g(:,:) = 0.0d0

do lig_atm_no = 1, res_index(target_lig_no)%n_atm
    i_atm = ii_L(lig_atm_no, target_lig_no-tn%nonlig)
    
    ! atdk_para index of ligand
    l_atm_type = atdk_para_idx(i_atm)
    l_q = q_s(i_atm)

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
            call calc_atdk3_elec(dist, p_q, l_q, tmp_E)
            f(1) = f(1) + tmp_E

            if (dist > NBcutoff) cycle

            ! atdk desolvation
            call calc_atdk3_solv(dist_sqr, atdk_para%solpar(l_atm_type),&
                                 atdk3_sol_V(j_atm), tmp_E)
            f(2) = f(2) - tmp_E

            ! atdk vdw/hbond
            if (.not. atdk_para%hbond_read(l_atm_type, p_atm_type)) then ! vdW
                call calc_atdk3_vdw(dist, p_atm_type, l_atm_type, tmp_E)
                f(3) = f(3) + tmp_E
            else ! Hbond
                dr(1:3) = dr(1:3)/dist
                call set_atdk3_Hbond_param(dr, j_atm, p_atm_type, hb_factor)
                call calc_atdk3_Hbond(dist, p_atm_type, l_atm_type, &
                                      hb_factor, tmp_E)
                f(4) = f(4) + tmp_E
            end if
        end do
    end do
end do

f(1) = f(1) * atdk_elec_w
f(2) = f(2) * atdk_solv_w
f(3) = f(3) * atdk_vdw_w
f(4) = f(4) * atdk_hbond_w

end subroutine calc_atdk3_inter_energy
!-------------------------------------------------------------------------------
subroutine calc_atdk3_intra_energy(f, g, calc_g)
!-------------------------------------------------------------------------------
! calculate autodock3 interaction energy between protein and ligand in
! continuous way.
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f(4), g(3, tn%atom)
logical, intent(in) :: calc_g
!
integer :: i_atm, j_atm, atm_1, atm_2
integer :: atm_type_1, atm_type_2
real(dp) :: dr(3), dist_sqr, dist
real(dp) :: hb_factor, tmp_E

f(:) = 0.0d0
g(:,:) = 0.0d0
do atm_1 = 1, res_index(target_lig_no)%n_atm - 1
    i_atm = ii_L(atm_1, target_lig_no-tn%nonlig)
    
    ! atdk_para index of ligand
    atm_type_1 = atdk_para_idx(i_atm)

    do atm_2 = atm_1 + 1, res_index(target_lig_no)%n_atm
        
        j_atm = ii_L(atm_2, target_lig_no-tn%nonlig)
        if (nb_matrix(atm_2, atm_1) == 0) cycle
        if (nb_matrix(atm_2, atm_1) == 4) cycle

        atm_type_2 = atdk_para_idx(j_atm)

        dr(1:3) = R(1:3, j_atm) - R(1:3, i_atm)
        dist_sqr = dot_product(dr, dr)
        dist = sqrt(dist_sqr)
        
        if (dist > NBcutoff) cycle

        ! atdk vdw/hbond
        if (.not. atdk_para%hbond_read(atm_type_1, atm_type_2)) then ! vdW
            call calc_atdk3_vdw(dist, atm_type_2, atm_type_1, tmp_E)
            f(3) = f(3) + tmp_E
        else ! Hbond
            hb_factor = 1.0d0
            call calc_atdk3_Hbond(dist, atm_type_2, atm_type_1, &
                                  hb_factor, tmp_E)
            f(4) = f(4) + tmp_E
        end if
    end do
end do

f(1) = f(1) * atdk_int_e_w
f(2) = f(2) * atdk_int_s_w
f(3) = f(3) * atdk_int_v_w
f(4) = f(4) * atdk_int_h_w

end subroutine calc_atdk3_intra_energy
!-------------------------------------------------------------------------------
subroutine calc_atdk3_prot_vdw(f, g, appl_pair, calc_g)
!-------------------------------------------------------------------------------
! calculate autodock3 interaction energy between protein and ligand in
! continuous way.
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f, g(:,:)
logical, intent(in) :: appl_pair(:,:), calc_g
!
integer :: i_atm, j_atm, ia
integer :: atm_type_1, atm_type_2
integer :: n_pair, i_pair(max_neigh)
real(dp) :: dr(3), dist, tmp_E, hb_factor

f = 0.0d0
g(:,:) = 0.0d0

do i_atm = 1, tn%stdatm
    atm_type_1 = atdk_para_idx(i_atm)
    call get_atdk3_pair(i_atm, appl_pair, n_pair, i_pair)

    do ia = 1, n_pair
        j_atm = i_pair(ia)
        atm_type_2 = atdk_para_idx(j_atm)

        dr(1:3) = R(1:3, j_atm) - R(1:3, i_atm)
        dist = dot_product(dr, dr)

        if (dist > NBcutoff2) cycle
        dist = sqrt(dist)
        
        ! atdk vdw/hbond
        if (.not. atdk_para%hbond_read(atm_type_1, atm_type_2)) then
            call calc_atdk3_vdw(dist, atm_type_2, atm_type_1, tmp_E)
            f = f + tmp_E
        else ! Hbond
            dr(1:3) = dr(1:3)/dist
            call set_atdk3_Hbond_param(dr, j_atm, atm_type_2, hb_factor)
            call calc_atdk3_Hbond(dist, atm_type_2, atm_type_1, &
                                  hb_factor, tmp_E)
            f = f + tmp_E
        end if
    end do
end do

end subroutine calc_atdk3_prot_vdw
!-------------------------------------------------------------------------------
subroutine get_atdk3_pair(i_atm, appl_pair, n_pair, i_pair)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_atm
logical, intent(in) :: appl_pair(:,:)
integer, intent(out) :: n_pair
integer, intent(out) :: i_pair(:)

integer :: ia_start, ia
integer :: i_res, j_res, j_atm

n_pair = 0
i_res = i_R(1, i_atm)

ia_start = i_P(i_atm)%pair_end_index(2)+1

do ia = ia_start, i_P(i_atm)%n_pair
    j_atm = i_P(i_atm)%i_pair(ia)
    j_res = i_R(1, j_atm)
    if (j_res == i_res) cycle
    if (j_atm > tn%stdatm) cycle
    if (.not. appl_pair(i_res, j_res)) cycle
    n_pair = n_pair + 1
    i_pair(n_pair) = j_atm
end do

end subroutine get_atdk3_pair
!-------------------------------------------------------------------------------
subroutine calc_atdk3_elec(dist, q1, q2, E)
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

end subroutine calc_atdk3_elec
!-------------------------------------------------------------------------------
subroutine calc_atdk3_solv(dist_sqr, solpar, vol, E)
!-------------------------------------------------------------------------------
! calculate autodock4 desolvation energy
! Edesolv = (S_lig*V_prot) * exp(const*(dist**2))
! sol_const = -1.0/(2.0*(sol_sig**2)) in autodock_vars.f90
!-------------------------------------------------------------------------------
real(dp), intent(in) :: dist_sqr
real(dp), intent(in) :: solpar, vol
real(dp), intent(out) :: E

E = 0.0d0

E = solpar*vol * exp(dist_sqr*sol_const)

end subroutine calc_atdk3_solv
!-------------------------------------------------------------------------------
subroutine calc_atdk3_vdw(dist, p_type, l_type, E)
!-------------------------------------------------------------------------------
! calculate autodock4 van der Waals interaction energy
! It's a flat-bottom 12-6 LJ potential 
!-------------------------------------------------------------------------------
real(dp), intent(in) :: dist
integer, intent(in) :: p_type, l_type
real(dp), intent(out) :: E
real(dp) :: expA, expB, Rij, Eij, cA, cB, rA, rB
real(dp) :: Rij_12, Rij_6
real(dp) :: sm_dist

E = 0.0d0

expA = 12.0d0
expB = 6.0d0

Rij = atdk_para%ljrij(l_type, p_type)
Eij = atdk_para%ljeij(l_type, p_type)

Rij_6 = Rij**6
Rij_12 = Rij_6**2

cA = Eij/(expA-expB) * (Rij_12)*expB
cB = Eij/(expA-expB) * (Rij_6)*expA

sm_dist = dist
call smooth_dist(sm_dist, Rij)

rB = sm_dist**expB
rA = rB**2

E = cA/rA - cB/rB

end subroutine calc_atdk3_vdw
!-------------------------------------------------------------------------------
subroutine calc_atdk3_Hbond(dist, p_type, l_type, hb_factor, E)
!-------------------------------------------------------------------------------
! calculate autodock4 Hbond interaction energy
! It's a flat-bottom 12-10 LJ potential with consideration of diretionality. 
!-------------------------------------------------------------------------------
real(dp), intent(in) :: dist
integer, intent(in) :: p_type, l_type
real(dp), intent(in) :: hb_factor
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
Eij = atdk_para%hbeij(l_type, p_type)

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

E = E * (hb_factor + (1.0d0 - hb_factor)*rsph)

end subroutine calc_atdk3_Hbond
!-------------------------------------------------------------------------------
subroutine set_atdk3_Hbond_param(dr, atm_idx, atm_type, hb_factor)
!-------------------------------------------------------------------------------
! to calculate autodock3 Hbond energy, get racc/rdon value.
!-------------------------------------------------------------------------------
real(dp), intent(in) :: dr(3)
integer, intent(in) :: atm_idx, atm_type
real(dp), intent(out) :: hb_factor

hb_factor = 1.0

if (atdk_para%type_name(atm_type) == 'H ') then ! current protein atom is 'H '
    call racc_rdon_H(atdk_hbond(atm_idx), dr)
    hb_factor = atdk_hbond(atm_idx)%racc
else if (atdk_para%type_name(atm_type) == 'O ') then ! curr protein atom is 'O '
    call racc_rdon_O(atdk_hbond(atm_idx), dr)
    hb_factor = atdk_hbond(atm_idx)%rdon
end if

end subroutine set_atdk3_Hbond_param
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
!-------------------------------------------------------------------------------
subroutine set_atdk3_hbond_flex(usc_list, n_usc)
!-------------------------------------------------------------------------------
integer, intent(in) :: usc_list(:), n_usc
integer :: i_usc, i_res, i_atm, atm_no

do i_usc = 1, n_usc
    i_res = usc_list(i_usc)
    do i_atm = 1, res_index(i_res)%n_atm
        atm_no = ii_R(i_atm, i_res)
        if (is_Hdon(atm_no)) then
            call HD_hbond_parameter(atdk_hbond(atm_no), atm_no)
        else if (is_Hacc_O(atm_no)) then
            call OA_SA_hbond_parameter(atdk_hbond(atm_no), atm_no)
        end if
    end do
end do

end subroutine set_atdk3_hbond_flex
!-------------------------------------------------------------------------------
subroutine calc_atdk_rsr_weight(atom_pair, r0, weight)
!-------------------------------------------------------------------------------
integer, intent(in) :: atom_pair(2)
real(dp), intent(in) :: r0
real(dp), intent(inout) :: weight
integer :: i_prot, i_lig, lig_type, idx_r
    
i_prot = atdk_para_idx(atom_pair(1))
i_lig = atdk_para_idx(atom_pair(2))
call get_lig_atom_type(i_lig, lig_type)

idx_r = min(int(r0*divisor) + 1, max_dist_short)
weight = -1 * min(T_vdw(idx_r, lig_type, i_prot), 0.0d0)
if (atdk_para%hbond_read(i_lig, i_prot)) then
    weight = weight / atdk_para%hbeij(i_lig, i_prot)
else
    weight = weight / atdk_para%ljeij(i_lig, i_prot)
end if

end subroutine calc_atdk_rsr_weight
!-------------------------------------------------------------------------------
subroutine FFT_rec_grid_atdk3(grid_info, dock_grid, FFT_grid_info, r_grid)
!-------------------------------------------------------------------------------
type(docking_grid_param_type), intent(in) :: grid_info, FFT_grid_info
type(docking_grid_type), intent(in) :: dock_grid
type(docking_FFT_grid_type), intent(inout) :: r_grid
!
integer :: n_type, i_type, j_type
integer :: x, y, z, hash, FFT_hash
real(dp) :: grid_v

n_type = atdk_para%n_atm_types(2)/2 + 1
r_grid%n_atdk = n_type

FFT_hash = FFT_grid_info%n_elem(1)*FFT_grid_info%n_elem(2)*FFT_grid_info%n_elem(3)
allocate(r_grid%atdk_grid(FFT_hash, n_type))

r_grid%atdk_grid = cmplx(0.0,0.0)

do i_type = 1, atdk_para%n_atm_types(2)+1
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

end subroutine FFT_rec_grid_atdk3
!-------------------------------------------------------------------------------
subroutine FFT_lig_grid_atdk3(FFT_grid_info, grid_info, &
                              ligand, fragment, n_frag_atm, R_frag, l_grid, &
                              ordering)
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
    n_type = atdk_para%n_atm_types(2)/2 + 1
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
                                cmplx(ref_lig(ligand%lig_type)%charge(atm_no),0.0)
end do

! others for vdw/Hbond/desolv: Lp(i,j,k) = 1 or 0
do atm_no = 1, n_frag_atm
    i_atm = order_by_type(atm_no)
    i_type = type_s(atm_no)
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
        j_type = i_type/2 + 1
        l_grid%atdk_grid(FFT_hash,j_type) = l_grid%atdk_grid(FFT_hash,j_type) + &
                                         cmplx(0.0, 1.0)
    else
        j_type = i_type/2 + 1
        l_grid%atdk_grid(FFT_hash,j_type) = l_grid%atdk_grid(FFT_hash,j_type) + & 
                                         cmplx(1.0, 0.0)
    end if
end do

end subroutine FFT_lig_grid_atdk3
!-------------------------------------------------------------------------------
END MODULE AUTODOCK3
!-------------------------------------------------------------------------------
