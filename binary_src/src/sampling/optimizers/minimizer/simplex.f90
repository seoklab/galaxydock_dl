!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/sampling/optimizers/minimizer/simplex.f90
!
! Description:
!   This module contains subroutines for simplex algorithm used in ligand
!   docking.
!
!-------------------------------------------------------------------------------
MODULE SIMPLEX
!-------------------------------------------------------------------------------

use globals
use in_out_structure
use ran
use mathfunctions, only: bound_ang, v_norm
use geometry, only: pert_rotation, average_quaternion, inv_expmap_quat, &
                    expmap_quat
!
use energy_vars
use energy_utils, only: update_R_for_ligand
use ligdock_energy, only: ligdock_energy_using_grid
!
use transrot_operator, only: translate_ligand, rotate_ligand
use ligand_operator, only: construct_ligand

implicit none
save
private

integer, parameter :: max_gene = 6 + max_lig_br - 1 ! trans(3) + rot(4) &
                                                    ! tor(max_lig_br - 1)
!-------------------------------------------------------------------------------
type vertice_type
!-------------------------------------------------------------------------------
real(dp) :: gene(max_gene)
real(dp) :: E(0:n_E_component_ligdock)
!-------------------------------------------------------------------------------
end type vertice_type
!-------------------------------------------------------------------------------

real(dp), allocatable :: pert_range(:), gene_avg(:)
real(dp), allocatable :: pert_range_full(:)
type(vertice_type), allocatable :: vertice(:)

public :: initialize_simplex
public :: finalize_simplex
public :: do_simplex

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_simplex(gene_dim)
!-------------------------------------------------------------------------------
integer, intent(in) :: gene_dim

allocate(gene_avg(gene_dim))
allocate(pert_range(gene_dim))
allocate(pert_range_full(gene_dim))
allocate(vertice(gene_dim-2))

call set_pert_range(pert_range, pert_range_full, gene_dim)

end subroutine initialize_simplex
!-------------------------------------------------------------------------------
subroutine finalize_simplex()
!-------------------------------------------------------------------------------
deallocate(gene_avg)
deallocate(pert_range)
deallocate(pert_range_full)
deallocate(vertice)

end subroutine finalize_simplex
!-------------------------------------------------------------------------------
subroutine set_pert_range(pert_range, pert_range_full, gene_dim)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: pert_range(:), pert_range_full(:)
integer, intent(in) :: gene_dim
real(dp), parameter :: l_trans = 0.2d0
real(dp), parameter :: l_angle = 5.0d0
real(dp), parameter :: l_axis = 0.0d0
real(dp), parameter :: l_tor = 5.0d0
!real(dp), parameter :: l_trans = 2.0d0
!real(dp), parameter :: l_angle = 50.0d0
!real(dp), parameter :: l_axis = 0.0d0
!real(dp), parameter :: l_tor = 50.0d0

! perturbation range for translation
pert_range(1:3) = l_trans
! perturbation range for rotation
pert_range(4:6) = l_angle * deg2rad
! perturbation range for torsion rotation
pert_range(7:gene_dim) = l_tor * deg2rad

pert_range_full(1:3) = l_trans*10.0d0
! perturbation range for rotation
pert_range_full(4:6) = l_angle * deg2rad*12.0d0
! perturbation range for torsion rotation
pert_range_full(7:gene_dim) = l_tor * deg2rad*36.0d0

end subroutine set_pert_range
!-------------------------------------------------------------------------------
subroutine do_simplex(protein, ligand, gene, E, nft, max_iter_in, pert_full)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(in) :: ligand
real(dp), intent(inout) :: gene(:)
real(dp), intent(out) :: E(0:n_E_component_ligdock)
integer, intent(inout) :: nft
integer, intent(in), optional :: max_iter_in
logical, intent(in), optional :: pert_full
!
integer :: i_iter
!integer, parameter :: max_iter = 300
integer :: max_iter
real(dp), parameter :: conv_crit = 0.01d0 
real(dp) :: conv_check, g(3,ligand%n_atm)
!
integer :: gene_dim, n_vertice
integer :: max_loc, min_loc, max_2_loc
!
integer :: i
type(vertice_type) :: temp_unit_1, temp_unit_2
logical :: replace
type(energy_type) :: ff

if (present(max_iter_in)) then
    max_iter = max_iter_in
else
    max_iter = 300
end if

gene_dim = 6 + ligand%n_br-1
n_vertice = (gene_dim - 3) + 1

if (present(pert_full) .and. pert_full) then
    call generate_vertice(protein, ligand, vertice, n_vertice, gene, gene_dim,&
                          pert_range_full, nft)
else
    call generate_vertice(protein, ligand, vertice, n_vertice, gene, gene_dim,&
                          pert_range, nft)
end if

call find_Emax_Emin_vertice(vertice, n_vertice, max_loc, min_loc, max_2_loc)

do i_iter = 1, max_iter
    call get_gene_avg(vertice, n_vertice, gene_avg, gene_dim, max_loc, min_loc)
    call reflection(temp_unit_1%gene, gene_avg, vertice(max_loc)%gene, gene_dim)
    call construct_ligand(protein, ligand, temp_unit_1%gene)
    call update_R_for_ligand(protein%ligand(ligand%lig_no), ligand%lig_no)
    call ligdock_energy_using_grid(protein, ligand, ff, g, .false.)
    temp_unit_1%E(:) = ff%ligdock(:)
    nft = nft + 1

    temp_unit_2 = temp_unit_1

    if (temp_unit_2%E(0) < vertice(min_loc)%E(0)) then
        call expansion(temp_unit_1%gene, gene_avg, vertice(max_loc)%gene, &
                       gene_dim)
        call construct_ligand(protein, ligand, temp_unit_1%gene)
        call update_R_for_ligand(protein%ligand(ligand%lig_no), ligand%lig_no)
        call ligdock_energy_using_grid(protein, ligand, ff, g, .false.)
        temp_unit_1%E(:) = ff%ligdock(:)
        nft = nft + 1

        if (temp_unit_1%E(0) < temp_unit_2%E(0)) then
            vertice(max_loc) = temp_unit_1
        else
            vertice(max_loc) = temp_unit_2
        end if
    else if (temp_unit_2%E(0) > vertice(max_2_loc)%E(0)) then
        replace = .false.
        if (temp_unit_2%E(0) < vertice(max_loc)%E(0)) then
            vertice(max_loc) = temp_unit_2
            replace = .true.
        end if
        call contraction(temp_unit_1%gene, gene_avg, vertice(max_loc)%gene, &
                         gene_dim)
        call construct_ligand(protein, ligand, temp_unit_1%gene)
        call update_R_for_ligand(protein%ligand(ligand%lig_no), ligand%lig_no)
        call ligdock_energy_using_grid(protein, ligand, ff, g, .false.)
        temp_unit_1%E(:) = ff%ligdock(:)
        nft = nft+1

        if (temp_unit_1%E(0) < vertice(max_loc)%E(0)) then
            vertice(max_loc) = temp_unit_1
            replace = .true.
        end if

        if (.not. replace) then
            do i = 1, n_vertice
                if (i == min_loc) cycle
                call contraction(vertice(i)%gene, vertice(i)%gene, &
                                 vertice(min_loc)%gene, gene_dim)
                call construct_ligand(protein, ligand, vertice(i)%gene)
                call update_R_for_ligand(protein%ligand(ligand%lig_no), ligand%lig_no)
                call ligdock_energy_using_grid(protein, ligand, ff, g, .false.)
                vertice(i)%E(:) = ff%ligdock(:)
                nft = nft + 1
            end do
        end if
    else
        vertice(max_loc) = temp_unit_2
    end if
   
    call find_Emax_Emin_vertice(vertice, n_vertice, max_loc, min_loc, max_2_loc)

    conv_check = abs(vertice(max_loc)%E(0) - vertice(min_loc)%E(0))
    if (conv_check < conv_crit) exit
end do

gene = vertice(min_loc)%gene(1:gene_dim)
E(:) = vertice(min_loc)%E(:)
call construct_ligand(protein, ligand, gene)
call update_R_for_ligand(protein%ligand(ligand%lig_no), ligand%lig_no)

end subroutine do_simplex
!-------------------------------------------------------------------------------
subroutine generate_vertice(protein, ligand, vertice, n_vertice, &
                            gene, gene_dim, pert_range, nft)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
type(vertice_type), intent(out) :: vertice(:)
integer, intent(in) :: n_vertice, gene_dim
real(dp), intent(in) :: gene(:), pert_range(:)
integer, intent(inout) :: nft
!
type(molecule_type) :: tmp_prot
type(energy_type) :: ff
integer :: i
real(dp) :: g(3,ligand%n_atm)
real(dp) :: axis(3)
real(dp) :: q(4), p(3)

tmp_prot = protein

do i = 1, n_vertice
    vertice(i)%gene(1:gene_dim) = gene(1:gene_dim)
end do
!
do i = 1, n_vertice - 1
    if (i < 4) then
        vertice(i)%gene(i) = vertice(i)%gene(i) &
                           + (2.0d0*random()-1.0d0)*pert_range(i)
    else if (i==4) then
        p = gene(4:6)
        call expmap_quat(p, q)
        call pert_rotation(pert_range(4), q)
        call inv_expmap_quat(q,p)
        vertice(i)%gene(4:6) = p
    else
        vertice(i)%gene(i+2) = vertice(i)%gene(i+2) &
                           + (2.0d0*random()-1.0d0)*pert_range(i+2)
    end if
end do

do i = 1, n_vertice
    call construct_ligand(tmp_prot, ligand, vertice(i)%gene)
    call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
    call ligdock_energy_using_grid(protein, ligand, ff, g, .false.)
    vertice(i)%E(:) = ff%ligdock(:)
    nft = nft + 1
end do

end subroutine generate_vertice
!-------------------------------------------------------------------------------
subroutine find_Emax_Emin_vertice(vertice, n_vertice, max_loc, min_loc, max_2_loc)
!-------------------------------------------------------------------------------
type(vertice_type), intent(in) :: vertice(:)
integer, intent(in) :: n_vertice
integer, intent(out) :: max_loc, min_loc, max_2_loc
real :: max_E, min_E, max_2_E
integer :: i

max_E = vertice(1)%E(0)
max_2_E = vertice(1)%E(0)
min_E = vertice(1)%E(0)
max_loc = 1
max_2_loc = 1
min_loc = 1

do i = 2, n_vertice
    if (vertice(i)%E(0) > max_E) then
        max_2_E = vertice(max_loc)%E(0)
        max_2_loc = max_loc
        max_E = vertice(i)%E(0)
        max_loc = i
    else if (vertice(i)%E(0) > max_2_E) then
        max_2_E = vertice(i)%E(0)
        max_2_loc = i
    else if (vertice(i)%E(0) < min_E) then
        min_E = vertice(i)%E(0)
        min_loc = i
    end if
end do

end subroutine find_Emax_Emin_vertice
!-------------------------------------------------------------------------------
subroutine get_gene_avg(vertice, n_vertice, gene_avg, gene_dim, max_loc, min_loc)
!-------------------------------------------------------------------------------
type(vertice_type), intent(in) :: vertice(:)
real(dp), intent(out) :: gene_avg(:)
integer, intent(in) :: gene_dim, n_vertice, max_loc, min_loc
!
integer :: i, j, i_q
real(dp) :: tor_diff
real(dp) :: q_s(4, n_vertice-1), q_avg(4), q(4)
  
gene_avg(:) = 0.0d0

i_q = 0 
do i = 1, n_vertice
    if (i == max_loc) cycle
    ! for trans
    do j = 1, 3
        gene_avg(j) = gene_avg(j) + vertice(i)%gene(j)
    end do
    ! for rotation
    i_q = i_q + 1
    call expmap_quat(vertice(i)%gene(4:6), q)
    q_s(:,i_q) = q

    ! for torsion
    do j = 7, gene_dim
        tor_diff = (vertice(i)%gene(j) - vertice(min_loc)%gene(j))
        gene_avg(j) = gene_avg(j) + vertice(min_loc)%gene(j) + bound_ang(tor_diff)
    end do
end do

!average rest
gene_avg(:) = gene_avg(:)/dble(n_vertice-1)

!average rotation
call average_quaternion(q_s, n_vertice-1, q_avg)
call inv_expmap_quat(q_avg, gene_avg(4:6))

do j = 7, gene_dim
    gene_avg(j) = bound_ang(gene_avg(j))
end do

end subroutine get_gene_avg
!-------------------------------------------------------------------------------
subroutine reflection(tmp_gene, gene_avg, max_gene, gene_dim)
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: tmp_gene(:)
real(dp), intent(in) :: gene_avg(:), max_gene(:)
integer, intent(in) :: gene_dim
integer :: i
real(dp) :: q1(4), q2(4), q3(4), p(3)
real(dp) :: p1(3)
do i = 1, gene_dim
    if (i < 4 .or. i > 6) then ! not rotation
        tmp_gene(i) = 2.0d0*gene_avg(i) - 1.0d0*max_gene(i)
    end if
end do

call expmap_quat(gene_avg(4:6), q1)
call expmap_quat(max_gene(4:6), q2)
if (dot_product(q1, q2) < 0.0d0) q2 = -q2
do i=1,4
    q3(i) = 2.0d0*q1(i) - 1.0d0*q2(i)
end do
call v_norm(q3)
call inv_expmap_quat(q3, p)
tmp_gene(4:6) = p

do i = 7, gene_dim
    tmp_gene(i) = bound_ang(tmp_gene(i))
end do

end subroutine reflection
!-------------------------------------------------------------------------------
subroutine expansion(tmp_gene, gene_avg, max_gene, gene_dim)
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: tmp_gene(:)
real(dp), intent(in) :: gene_avg(:), max_gene(:)
integer, intent(in) :: gene_dim
integer :: i
real(dp) :: q1(4), q2(4), q3(4), p(3)
real(dp) :: p1(3)
do i = 1, gene_dim
    if (i < 4 .or. i > 6) then ! not rotation
       tmp_gene(i) = 3.0d0*gene_avg(i) - 2.0d0*max_gene(i)
    end if
end do

call expmap_quat(gene_avg(4:6), q1)
call expmap_quat(max_gene(4:6), q2)
if (dot_product(q1, q2) < 0.0d0) q2 = -q2
do i=1,4
    q3(i) = 3.0d0*q1(i) - 2.0d0*q2(i)
end do
call v_norm(q3)
call inv_expmap_quat(q3, p)
tmp_gene(4:6) = p

do i = 7, gene_dim
    tmp_gene(i) = bound_ang(tmp_gene(i))
end do

end subroutine expansion
!-------------------------------------------------------------------------------
subroutine contraction(tmp_gene, gene_avg, max_gene, gene_dim)
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: tmp_gene(:)
real(dp), intent(in) :: gene_avg(:), max_gene(:)
integer, intent(in) :: gene_dim
integer :: i
real(dp) :: tor_diff
real(dp) :: q1(4), q2(4), q3(4), p(3), p1(3)
do i = 1, gene_dim
    if (i < 7) then
        if (i < 4) then ! not rotation
            tmp_gene(i) = 0.5d0*(gene_avg(i)+max_gene(i))
        end if
        !tmp_gene(i) = 0.5d0*(gene_avg(i)+max_gene(i))
    else
        tor_diff = (max_gene(i) - gene_avg(i))
        tmp_gene(i) = 0.5d0*(2.0d0*gene_avg(i) + bound_ang(tor_diff))
        tmp_gene(i) = bound_ang(tmp_gene(i))
    end if
end do

call expmap_quat(gene_avg(4:6), q1)
call expmap_quat(max_gene(4:6), q2)
if (dot_product(q1, q2) < 0.0d0) q2 = -q2
do i=1,4
    q3(i) = 0.5d0*q1(i) + 0.5d0*q2(i)
end do
call v_norm(q3)
call inv_expmap_quat(q3, p)
tmp_gene(4:6) = p

end subroutine contraction
!-------------------------------------------------------------------------------
END MODULE SIMPLEX
!-------------------------------------------------------------------------------
