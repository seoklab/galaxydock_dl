!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/utils/rmsd.f90
!
! Description:
!  This module contains subroutines for RMSD calculation.
!  When molecules to be compared are not super-imposed, 
!  calc_rmsd can be used to find best-fitting matrix and according RMSD value.
!  calc_rmsd_simple can be used instead for calculation without super-imposition.
!
! Reference: Evangelos A. Coutsias, Chaok Seok*, and Ken A. Dill, Using 
!            Quaternions to Calculate RMSD, J. Comput. Chem. 25, 1849.
!
!
!-------------------------------------------------------------------------------
MODULE RMSD
!-------------------------------------------------------------------------------
use globals
use logger,   only: terminate_with_error
use geometry, only: rotation_matrix

implicit none
private

public :: calc_rmsd
public :: calc_rmsd_simple
public :: calc_rmsd_simple_protein
public :: calc_rmsd_het
public :: ls_rmsd
public :: prep_RMSD_PM_calculation
public :: mapping_type_for_RMSD_PM
public :: calc_RMSD_PM
public :: calc_RMSD_PM_2
public :: dstmev
public :: get_coord

CONTAINS
!-------------------------------------------------------------------------------
subroutine calc_rmsd(protein1, protein2, rmsd_type, option, U, error, res1_in, res2_in)
!-------------------------------------------------------------------------------
! This automatically superimpose two proteins and returns best-fitted
! RMSD and rotation matrix for superimposition
!-------------------------------------------------------------------------------
integer, intent(in) :: option
type(molecule_type), intent(in) :: protein1, protein2
character(len=2), intent(in) :: rmsd_type
integer, intent(in), optional :: res1_in, res2_in
real(dp), dimension(:,:), intent(out) :: U
real(dp), intent(out) :: error
real(dp), dimension(3,protein1%n_res*max_atm) :: coord1, coord2
integer :: num_atm1, num_atm2, res1, res2

if (present(res1_in)) then
    res1 = res1_in
else
    res1 = 1
end if

if (present(res2_in)) then
    res2 = res2_in
else
    res2 = protein1%n_res
end if

call get_coord(res1, res2, protein1%residue(:), coord1, num_atm1, rmsd_type) 
call get_coord(res1, res2, protein2%residue(:), coord2, num_atm2, rmsd_type) 

if (num_atm1 /= num_atm2) then
    call terminate_with_error('ERROR: RMSD Number of atoms in two proteins are not the same')
end if

call ls_rmsd(num_atm1, coord1, coord2, option, U, error)

end subroutine calc_rmsd
!-------------------------------------------------------------------------------
subroutine calc_rmsd_simple(crd1, crd2, n_atm, rmsd)
!-------------------------------------------------------------------------------
! Simply calculate RMSD WITHOUT superimposition (RMSD of given coordinates)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: rmsd
integer, intent(in) :: n_atm
real(dp), intent(in) :: crd1(3,n_atm), crd2(3,n_atm)
integer :: i_atm
real(dp) :: r(3), r_sqr

rmsd = 0.0d0
do i_atm = 1, n_atm
    r(:) = crd1(:,i_atm) - crd2(:,i_atm)
    r_sqr = dot_product(r(:),r(:))
    rmsd = rmsd + r_sqr
end do

rmsd = rmsd/n_atm
rmsd = rmsd**0.5

end subroutine calc_rmsd_simple
!-------------------------------------------------------------------------------
subroutine calc_rmsd_simple_protein(protein1, protein2, error, res1_in,res2_in,&
rmsd_type)
!-------------------------------------------------------------------------------
! Simply calculate RMSD WITHOUT superimposition (RMSD of given proteins)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein1, protein2
integer, intent(in), optional :: res1_in, res2_in
character(len=2), intent(in), optional :: rmsd_type
!res1_in : start residue
!res2_in : end residue
real(dp), intent(out) :: error
real(dp), dimension(3,protein1%n_res*max_atm) :: coord1, coord2
integer :: num_atm1, num_atm2, res1, res2, n_res
character(len=2) :: aa_type

if (present(res1_in)) then
    res1 = res1_in
else
    res1 = 1
end if

if (present(res2_in)) then
    res2 = res2_in
else
    res2 = protein1%n_res
end if

if (present(rmsd_type)) then
    aa_type = rmsd_type
else
    aa_type = 'aa'
end if
n_res = res2 - res1 + 1
call get_coord(res1, res2, protein1%residue(:), coord1, num_atm1, aa_type) 
call get_coord(res1, res2, protein2%residue(:), coord2, num_atm2, aa_type)

if (num_atm1 /= num_atm2) then
    call terminate_with_error('ERROR: RMSD Number of atoms in two proteins are not the same')
end if
call calc_rmsd_simple(coord1, coord2, num_atm1, error)

end subroutine calc_rmsd_simple_protein
!-------------------------------------------------------------------------------
subroutine calc_rmsd_het(hetmol1, hetmol2, rmsd_type, option, U, error)
!-------------------------------------------------------------------------------
! Similar to calc_rmsd, but in hetero molecule version
!-------------------------------------------------------------------------------
integer, intent(in) :: option
type(residue_type), intent(in) :: hetmol1, hetmol2
character(len=2), intent(in) :: rmsd_type
real(dp), dimension(:,:), intent(out) :: U
real(dp), intent(out) :: error
real(dp), dimension(3,hetmol1%n_atm) :: coord1, coord2
integer :: num_atm1, num_atm2

call get_coord_het(hetmol1, coord1, num_atm1, rmsd_type) 
call get_coord_het(hetmol2, coord2, num_atm2, rmsd_type) 

if (num_atm1 /= num_atm2) then
    call terminate_with_error('ERROR: RMSD. Number of atoms in two ligands are not the same')
end if

call ls_rmsd(num_atm1, coord1, coord2, option, U, error)

end subroutine calc_rmsd_het
!-------------------------------------------------------------------------------
subroutine ls_rmsd(n, coord1, coord2, option, U, error)
!-------------------------------------------------------------------------------
!  This subroutine calculates the least square rmsd of the two coordinate
!  sets coord1(3,n) and coord2(3,n) using a method based on quaternion.
!  If option=1, then the best-fit coord of coord1 is returned
!-------------------------------------------------------------------------------
integer, intent(in) :: n, option
real(dp), dimension(:,:), intent(inout) :: coord1, coord2
real(dp), dimension(:,:), intent(out) :: U
real(dp), intent(out) :: error
real(dp), dimension(3) :: x_center, y_center
integer :: i, j
real(dp), dimension(3,n) :: x, y
real(dp), dimension(n) :: xi, yi
real(dp) :: x_norm, y_norm, lambda
real(dp), dimension(3,3) :: Rmatrix
real(dp), dimension(4,4) :: S
real(dp), dimension(4) :: q

! make copies of the original coordinates
x(:,1:n) = coord1(:,1:n)  
y(:,1:n) = coord2(:,1:n)

! calculate the barycenters, centroidal coordinates, and the norms
x_norm = 0.0d0
y_norm = 0.0d0
do i = 1, 3
    xi = x(i,:)
    yi = y(i,:)
    x_center(i) = sum(xi)/dble(n)
    y_center(i) = sum(yi)/dble(n)
    xi = xi - x_center(i)
    yi = yi - y_center(i)
    x(i,:) = xi
    y(i,:) = yi
    x_norm = x_norm + dot_product(xi, xi)
    y_norm = y_norm + dot_product(yi, yi)
end do

! calculate the R matrix
do i = 1, 3
    do j = 1, 3
        Rmatrix(i,j) = dot_product(x(i,:),y(j,:))
    end do
end do

! lower triangular part of the S matrix
S(1, 1) = Rmatrix(1, 1) + Rmatrix(2, 2) + Rmatrix(3, 3)
S(2, 1) = Rmatrix(2, 3) - Rmatrix(3, 2)
S(1, 2) = S(2, 1)
S(2, 2) = Rmatrix(1, 1) - Rmatrix(2, 2) - Rmatrix(3, 3)
S(3, 1) = Rmatrix(3, 1) - Rmatrix(1, 3)
S(1, 3) = S(3, 1)
S(3, 2) = Rmatrix(1, 2) + Rmatrix(2, 1)
S(2, 3) = S(3, 2)
S(3, 3) =-Rmatrix(1, 1) + Rmatrix(2, 2) - Rmatrix(3, 3)
S(4, 1) = Rmatrix(1, 2) - Rmatrix(2, 1)
S(1, 4) = S(4, 1)
S(4, 2) = Rmatrix(1, 3) + Rmatrix(3, 1)
S(2, 4) = S(4, 2)
S(4, 3) = Rmatrix(2, 3) + Rmatrix(3, 2)
S(3, 4) = S(4, 3)
S(4, 4) =-Rmatrix(1, 1) - Rmatrix(2, 2) + Rmatrix(3, 3) 

! Calculate eigenvalues and eigenvectors. For explanation about arguments,
! see the end of this subroutine.
call dstmev(S, lambda, q)

if (option == 1) then
    call rotation_matrix(3, q, U)
    do i = 1, n 
        coord1(:,i) = matmul(U,x(:,i)) + y_center(:)
    end do
end if

! lambda is the largest eigenvalue
error = sqrt(max(0.0d0,(x_norm + y_norm - 2.0d0*lambda))/dble(n))

end subroutine ls_rmsd
!-------------------------------------------------------------------------------
subroutine prep_RMSD_PM_calculation(types, n_atm, n_type, n_atm_s, atm_id_s)
!-------------------------------------------------------------------------------
character(len=*), intent(in) :: types(:)
integer, intent(in) :: n_atm
integer, intent(out) :: n_type, n_atm_s(:), atm_id_s(:,:)
character(len=20) :: used_atm_types(n_atm)
integer :: i_atm, i_type
logical :: is_new

used_atm_types(:) = ''
n_atm_s(:) = 0
atm_id_s(:,:) = 0
n_type = 0

do i_atm = 1, n_atm
    if (types(i_atm)(1:2) == 'H.') cycle
    is_new = .true.
    do i_type = 1, n_type
        if (trim(used_atm_types(i_type)) == trim(types(i_atm))) then
            is_new = .false.
            n_atm_s(i_type) = n_atm_s(i_type) + 1
            atm_id_s(n_atm_s(i_type), i_type) = i_atm
            exit
        end if
    end do

    if (is_new) then
        n_type = n_type + 1
        used_atm_types(n_type) = types(i_atm)
        n_atm_s(n_type) = 1
        atm_id_s(n_atm_s(i_type), i_type) = i_atm
    end if
end do

end subroutine prep_RMSD_PM_calculation
!-------------------------------------------------------------------------------
subroutine mapping_type_for_RMSD_PM(type_1, type_2, n_type_1, n_type_2, &
                                    atm_id_s_1, atm_id_s_2, type_map)
!-------------------------------------------------------------------------------
character(len=*), intent(in) :: type_1(:), type_2(:)
integer, intent(in) :: n_type_1, n_type_2
integer, intent(in) :: atm_id_s_1(:,:), atm_id_s_2(:,:)
integer, intent(out) :: type_map(:)
integer :: i_type, j_type, idx1, idx2

type_map(:) = 0

!print*, n_type_1, n_type_2
do i_type = 1, n_type_1
!    print*, type_1(atm_id_s_1(1,i_type)), type_2(atm_id_s_2(1,i_type))
    do j_type = 1, n_type_2
        idx1 = atm_id_s_1(1,i_type)
        idx2 = atm_id_s_2(1,j_type)
        if (type_1(idx1) == type_2(idx2)) then
!            print*, type_1(idx1), type_2(idx2)
            type_map(i_type) = j_type
            exit
        end if
    end do
end do

end subroutine mapping_type_for_RMSD_PM
!-------------------------------------------------------------------------------
subroutine calc_RMSD_PM(crd1, crd2, n_type, n_atm_s, atm_id_s, rmsd)
!-------------------------------------------------------------------------------
! Calculate RMSD with consideration of symmetry (RMSD-PM)
! Instead of consideration of exact symmetry, it uses property(type)-match method.
! DOCK6 with hungarian algorithm
!-------------------------------------------------------------------------------
real(dp), intent(in) :: crd1(:,:), crd2(:,:)
integer, intent(in) :: n_type, n_atm_s(:), atm_id_s(:,:)
real(dp), intent(out) :: rmsd
integer :: n_heavy, i_type
real(dp) :: dr(3), d_sqr

rmsd = 0.0d0
n_heavy = 0

do i_type = 1, n_type
    n_heavy = n_heavy + n_atm_s(i_type)
    if (n_atm_s(i_type) == 1) then
        dr(1:3) = crd2(:,atm_id_s(1,i_type)) - crd1(:,atm_id_s(1,i_type))
        rmsd = rmsd + dot_product(dr,dr)
    else
        call hungarian_algorithm(crd1, crd2, n_atm_s(i_type), &
                                 atm_id_s(:,i_type), atm_id_s(:,i_type), d_sqr)
        rmsd = rmsd + d_sqr
    end if
end do

rmsd = sqrt(rmsd / n_heavy)

end subroutine calc_RMSD_PM
!-------------------------------------------------------------------------------
subroutine calc_RMSD_PM_2(crd1, crd2, n_type, type_mapping, &
                          n_atm_s_1, n_atm_s_2, atm_id_s_1, atm_id_s_2, rmsd)
!-------------------------------------------------------------------------------
! Calculate RMSD with consideration of symmetry (RMSD-PM)
! Instead of consideration of exact symmetry, it uses property(type)-match method.
! DOCK6 with hungarian algorithm
!-------------------------------------------------------------------------------
real(dp), intent(in) :: crd1(:,:), crd2(:,:)
integer, intent(in) :: n_type, type_mapping(:)
integer, intent(in) :: n_atm_s_1(:), n_atm_s_2(:)
integer, intent(in) :: atm_id_s_1(:,:), atm_id_s_2(:,:)
real(dp), intent(out) :: rmsd
integer :: n_heavy, i_type, j_type
real(dp) :: dr(3), d_sqr

rmsd = 0.0d0
n_heavy = 0

do i_type = 1, n_type
    n_heavy = n_heavy + n_atm_s_1(i_type)
    j_type = type_mapping(i_type)
    if (n_atm_s_1(i_type) /= n_atm_s_2(j_type)) then
        call terminate_with_error('ERROR: # of types are different in RMSD_PM')
    end if

    if (n_atm_s_1(i_type) == 1) then
        dr(1:3) = crd2(:,atm_id_s_2(1,j_type)) - crd1(:,atm_id_s_1(1,i_type))
        rmsd = rmsd + dot_product(dr,dr)
    else
        call hungarian_algorithm(crd1, crd2, n_atm_s_1(i_type), &
                                 atm_id_s_1(:,i_type), atm_id_s_2(:,j_type),&
                                 d_sqr)
        rmsd = rmsd + d_sqr
    end if
end do

rmsd = sqrt(rmsd / n_heavy)

end subroutine calc_RMSD_PM_2
!-------------------------------------------------------------------------------
subroutine hungarian_algorithm(crd1, crd2, n_atm, atm_id_1, atm_id_2, d_sqr)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: crd1(:,:), crd2(:,:)
integer, intent(in) :: n_atm, atm_id_1(:), atm_id_2(:)
real(dp), intent(out) :: d_sqr
!
real(dp) :: D_mat(n_atm, n_atm), M(n_atm, n_atm), dr(3)
integer :: matched(n_atm), matrix_case(n_atm, n_atm), n_assigned
integer :: i_atm, j_atm, idx1, idx2

call construct_dist_matrix(crd1, crd2, n_atm, atm_id_1, atm_id_2, D_mat)
M = D_mat
call subtract_row_elem(M, n_atm)
call subtract_col_elem(M, n_atm)

do
    call assign_matched_zeros(M, matched, n_atm, n_assigned)
    if (n_assigned == n_atm) exit
    call draw_line(M, matched, matrix_case, n_atm)
    call update_matrix(M, matrix_case, n_atm)
end do

d_sqr = 0.0d0
do i_atm = 1, n_atm
    j_atm = matched(i_atm)
    idx2 = atm_id_2(j_atm)
    idx1 = atm_id_1(i_atm)
    dr(1:3) = crd2(:,idx2) - crd1(:,idx1)
    d_sqr = d_sqr + dot_product(dr,dr)
end do

end subroutine hungarian_algorithm
!-------------------------------------------------------------------------------
subroutine construct_dist_matrix(crd1, crd2, n_atm, atm_id_1, atm_id_2, D_mat)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: crd1(:,:), crd2(:,:)
integer, intent(in) :: n_atm, atm_id_1(:), atm_id_2(:)
real(dp), intent(out) :: D_mat(:,:)
!
integer :: i_atm, j_atm, idx1, idx2
real(dp) :: dr(3)

do i_atm = 1, n_atm
    idx1 = atm_id_1(i_atm)
    do j_atm = 1, n_atm
        idx2 = atm_id_2(j_atm)
        dr(1:3) = crd2(:,idx2) - crd1(:,idx1)
        D_mat(i_atm, j_atm) = dot_product(dr,dr)
    end do
end do

end subroutine construct_dist_matrix
!-------------------------------------------------------------------------------
subroutine subtract_row_elem(M, n_atm)
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: M(:,:)
integer, intent(in) :: n_atm
real(dp) :: min_val
integer :: i_atm, j_atm

do i_atm = 1, n_atm
    ! Get minimum value in the row.
    min_val = minval(M(i_atm,1:n_atm))
    ! Subtract minimum value from all element in the row.
    do j_atm = 1, n_atm
        M(i_atm,j_atm) = M(i_atm,j_atm) - min_val
    end do
end do

end subroutine subtract_row_elem
!-------------------------------------------------------------------------------
subroutine subtract_col_elem(M, n_atm)
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: M(:,:)
integer, intent(in) :: n_atm
real(dp) :: min_val
integer :: i_atm, j_atm

do i_atm = 1, n_atm
    ! Get minimum value in the column.
    min_val = minval(M(1:n_atm, i_atm))
    ! Subtract minimum value from all element in the column.
    do j_atm = 1, n_atm
        M(j_atm,i_atm) = M(j_atm,i_atm) - min_val
    end do
end do

end subroutine subtract_col_elem
!-------------------------------------------------------------------------------
subroutine assign_matched_zeros(M, matched, n_atm, n_assigned)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: M(:,:)
integer, intent(out) :: matched(:), n_assigned
integer, intent(in) :: n_atm
!
logical :: row_assigned(n_atm), column_assigned(n_atm)
integer :: row_count(n_atm), column_count(n_atm)
integer :: i_atm, j_atm, min_value, row_or_column, row, column

matched(:) = -100
row_assigned(:) = .false.
column_assigned(:) = .false.

n_assigned = 0

do
    row_count(:) = 0
    column_count(:) = 0
    ! 1. count the number of 0s in each row and column
    do j_atm = 1, n_atm
        do i_atm = 1, n_atm
            if (M(i_atm,j_atm) == 0) then
                if (row_assigned(i_atm)) cycle
                if (column_assigned(j_atm)) cycle
                row_count(i_atm) = row_count(i_atm) + 1
                column_count(j_atm) = column_count(j_atm) + 1
            end if
        end do
    end do
    ! 2. Identify the row or column that has the least number of 0s
    min_value = n_atm + 1
    row_or_column = 0
    row = -1
    column = -1
    do i_atm = 1, n_atm
        if (row_count(i_atm) /= 0 .and. row_count(i_atm) < min_value) then
            min_value = row_count(i_atm)
            row = i_atm
            row_or_column = 1
            column = -1
        end if
        if (column_count(i_atm) /= 0 .and. column_count(i_atm) < min_value) then
            min_value = column_count(i_atm)
            column = i_atm
            row_or_column = 2
            row = -1
        end if
    end do
    ! 3. Make an assignment in the row or column with the least number of 0s
    if (row_or_column == 0) exit
    if (row_or_column == 1) then
        do i_atm = 1, n_atm
            if (M(row,i_atm) == 0 .and. (.not. column_assigned(i_atm))) then
                matched(row) = i_atm
                row_assigned(row) = .true.
                column_assigned(i_atm) = .true.
                exit
            end if
        end do
    else if (row_or_column == 2) then
        do i_atm = 1, n_atm
            if (M(i_atm,column) == 0 .and. (.not. row_assigned(i_atm))) then
                matched(i_atm) = column
                column_assigned(column) = .true.
                row_assigned(i_atm) = .true.
                exit
            end if
        end do
    end if
    ! 4. Increment the counter and repeat
    n_assigned = n_assigned + 1
    if (n_assigned == n_atm) exit
end do

end subroutine assign_matched_zeros
!-------------------------------------------------------------------------------
subroutine draw_line(M, matched, matrix_case, n_atm)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: M(:,:)
integer, intent(in) :: matched(:)
integer, intent(inout) :: matrix_case(:,:)
integer, intent(in) :: n_atm
logical :: row_assigned(n_atm), column_assigned(n_atm)
logical :: flag_update
integer :: i_atm, j_atm

row_assigned(:) = .false.
column_assigned(:) = .false.
matrix_case(:,:) = 0

! 1. Mark all rows that have no assignment
do i_atm = 1, n_atm
    if (matched(i_atm) == -100) then
        row_assigned(i_atm) = .true.
    end if
end do

do
    flag_update = .false.
    ! 2. Mark all columns that have 0s in those rows
    do i_atm = 1, n_atm
        if (.not. row_assigned(i_atm)) cycle
        do j_atm = 1, n_atm
            if (M(i_atm,j_atm) == 0.0) then
                column_assigned(j_atm) = .true.
            end if
        end do
    end do
    ! 3. Mark all rows having assignments in those columns
    do i_atm = 1, n_atm
        if (.not. column_assigned(i_atm)) cycle
        do j_atm = 1, n_atm
            if (matched(j_atm) == i_atm .and. (.not. row_assigned(j_atm))) then
                row_assigned(j_atm) = .true.
                flag_update = .true.
            end if
        end do
    end do
    if (.not. flag_update) exit
end do

! 4. Draw lines through marked columns and unmarked rows
do i_atm = 1, n_atm
    if (.not. row_assigned(i_atm)) then
        do j_atm = 1, n_atm
            matrix_case(i_atm, j_atm) = matrix_case(i_atm,j_atm) + 1
        end do
    end if
    if (column_assigned(i_atm)) then
        do j_atm = 1, n_atm
            matrix_case(j_atm, i_atm) = matrix_case(j_atm, i_atm) + 1
        end do
    end if
end do

end subroutine draw_line
!-------------------------------------------------------------------------------
subroutine update_matrix(M, matrix_case, n_atm)
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: M(:,:)
integer, intent(in) :: matrix_case(:,:)
integer, intent(in) :: n_atm
real(dp) :: min_value
integer :: i_atm, j_atm

min_value = 999999.9d0

! 1. find the minimum value in all cells of M
do j_atm = 1, n_atm
    do i_atm = 1, n_atm
        if (matrix_case(i_atm, j_atm) /= 0) cycle
        if (M(i_atm, j_atm) < min_value) then
            min_value = M(i_atm, j_atm)
        end if
    end do
end do
if (min_value == 0.0d0) stop

! 2. update all the cells of M
do i_atm = 1, n_atm
    do j_atm = 1, n_atm
        if (matrix_case(i_atm,j_atm) == 0) then
            M(i_atm,j_atm) = M(i_atm,j_atm) - min_value
        else if (matrix_case(i_atm,j_atm) == 1) then
            cycle
        else if (matrix_case(i_atm,j_atm) == 2) then
            M(i_atm, j_atm) = M(i_atm,j_atm) + min_value
        end if
    end do
end do

end subroutine update_matrix
!-------------------------------------------------------------------------------
subroutine get_coord(res1, res2, residue, coord, num_atm, atom_type)
!-------------------------------------------------------------------------------
! Extract coordinates according to the given atom_type
! This subroutine is quite useful, although not used frequently yet (2011 PHB)
! Someday this subroutine may replace inefficient parts of the code...
!-------------------------------------------------------------------------------
integer, intent(in) :: res1, res2
type(residue_type), intent(in) :: residue(:)
real(dp), dimension(:,:), intent(out) :: coord
integer, intent(out) :: num_atm
character(len=2), intent(in) :: atom_type
integer :: i_res, i_atm

num_atm = 0

if (atom_type == 'aa') then ! all atom
    do i_res = res1, res2
        do i_atm = 1, residue(i_res)%n_atm
            num_atm = num_atm + 1
            coord(:,num_atm) = residue(i_res)%R(:,i_atm)
        end do
    end do
else if (atom_type == 'bb') then ! all backbone atoms (N,CA,C,O)
    do i_res = res1, res2
        do i_atm = 1, 4
            num_atm = num_atm + 1
            coord(:,num_atm) = residue(i_res)%R(:,res_index(i_res)%bb_id(i_atm,2))
        end do
    end do
else if (atom_type == 'mc') then ! main-chain atoms (N,CA,C)
    do i_res = res1, res2
        do i_atm = 1, 3
            num_atm = num_atm + 1
            coord(:,num_atm) = residue(i_res)%R(:,res_index(i_res)%bb_id(i_atm,2))
        end do
    end do
else if (atom_type == 'cb') then ! backbone + CB
    do i_res = res1, res2
        do i_atm = 1, 4
            num_atm = num_atm + 1
            coord(:,num_atm) = residue(i_res)%R(:,res_index(i_res)%bb_id(i_atm,2))
        end do
        num_atm = num_atm + 1
        coord(:,num_atm) = residue(i_res)%R(:,res_index(i_res)%Cb_id(2))
    end do
else if (atom_type == 'ca') then ! CA
    do i_res = res1, res2
        num_atm = num_atm + 1
        coord(:,num_atm) = residue(i_res)%R(:,res_index(i_res)%Ca_id(2))
    end do
end if

end subroutine get_coord
!-------------------------------------------------------------------------------
subroutine get_coord_het(hetmol, coord, num_atm, atom_type)
!-------------------------------------------------------------------------------
! Similar to get_coord, but in hetero molecule version
! It should move to rmsd.f90
!-------------------------------------------------------------------------------
integer, intent(out) :: num_atm
character(len=2) :: atom_type
character(len=4) :: atom
type(residue_type), intent(in) :: hetmol
real(dp), dimension(:,:), intent(out) :: coord
integer :: i_ref_res, i_atm

num_atm = 0

i_ref_res = hetmol%res_type
do i_atm = 1, ref_res(i_ref_res)%n_atm
    atom = ref_res(i_ref_res)%atom_name(i_atm)
    if (atom_type == 'aa'  .or. &
       (atom_type == 'ha' .and. (atom(1:1) == 'C' .or. atom(1:1) == 'N' .or. &
                                 atom(1:1) == 'O' .or. atom(1:1) == 'S' .or. atom(1:1) == 'P'))) then 
        num_atm = num_atm + 1
        coord(:,num_atm) = hetmol%R(:,i_atm)
    end if
end do

end subroutine get_coord_het
!==============================================================================
! Subroutines for finding out best-alignment quaternion using SVD
!==============================================================================
subroutine dstmev(A, lambda, evec)
!-------------------------------------------------------------------------------
! a simple subroutine to compute the leading eigenvalue and eigenvector
! of a symmetric, traceless 4x4 matrix A by an inverse power iteration:
! (1) the matrix is converted to tridiagonal form by 3 Givens
! rotations;  V*A*V' = T
! (2) Gershgorin's theorem is used to estimate a lower
! bound for the leading negative eigenvalue:
! lambda_1 > g=min(T11-t12,-t21+T22-t23,-t32+T33-t34,-t43+T44)
!          =
! where tij=abs(Tij)
! (3) Form the positive definite matrix 
!     B = T-gI
! (4) Use svd (algorithm svdcmp from "Numerical Recipes")
!     to compute eigenvalues and eigenvectors for SPD matrix B
! (5) Shift spectrum back and keep leading singular vector
!     and largest eigenvalue.
! (6) Convert eigenvector to original matrix A, through 
!     multiplication by V'.  
!-------------------------------------------------------------------------------
real(dp), dimension(4,4) :: A, T, V, SV
integer :: i
integer, dimension(1) :: max_loc
real(dp), dimension(4) :: evec, SW
real(dp), dimension(8) :: rv1
real(dp):: lambda

!(I).Convert to tridiagonal form, keeping similarity transform
! (a product of 3 Givens rotations)
call givens4(A,T,V)

!(II)Estimate lower bound of smallest eigenvalue by Gershgorin's theorem
lambda = min(T(1,1)-abs(T(1,2)),-abs(T(2,1))+T(2,2)-abs(T(2,3)),&
        -abs(T(3,2))+T(3,3)-abs(T(3,4)),-abs(T(4,3))+T(4,4))

!(III). Form positive definite matrix     T <== lambda*I - T
do i = 1, 4
    T(i,i) = T(i,i) - lambda
end do

!(IV). Compute singular values/vectors of SPD matrix B
call svdcmp(4, T, 4, 4, SW, SV, rv1)

!(V). Shift spectrum back
max_loc = maxloc(SW) 
lambda = SW(max_loc(1)) + lambda

!(VI). Convert eigenvector to original coordinates: (V is transposed!)
evec = matmul(V,SV(:,max_loc(1)))

end subroutine dstmev
!-------------------------------------------------------------------------------
subroutine givens4(S,T,V)
!-------------------------------------------------------------------------------
! Called by 'dstmev'. Returns Givens rotations with a given 'S'.
! Don't ask me what this is :)
!-------------------------------------------------------------------------------
real(dp), dimension(4,4), intent(in)  :: S
real(dp), dimension(4,4), intent(out) :: T,V
real(dp) :: c1, c2, c3, s1, s2, s3, r1, r2, r3, c1c2, s1c2

!performs givens rotations to reduce symmetric 4x4 matrix to tridiagonal
T = S
V = 0.0d0

! Zero out entries T(4,1) and T(1,4)
! compute cos and sin of rotation angle in the 3-4 plane
r1 = pythag(T(3,1),T(4,1))
if(r1 .ne. 0.d0) then
    c1 = T(3,1)/r1; s1 = T(4,1)/r1
    V(3,3) = c1   ; V(3,4) = s1
    V(4,3) = -s1  ; V(4,4) = c1
    T(3,1) = r1   ; T(4,1) = 0.0d0
    T(3:4,2:4) = matmul(V(3:4,3:4),T(3:4,2:4))
    T(1:2,3:4) = transpose(T(3:4,1:2))
    T(3:4,3:4) = matmul(T(3:4,3:4),transpose(V(3:4,3:4)))
else
    c1 = 1.d0; s1 = 0.d0
end if

! Zero out entries T(3,1) and T(1,3)
! compute cos and sin of rotation angle in the 2-3 plane
r2 = pythag(T(3,1), T(2,1))
if(r2 .ne. 0.d0) then
    c2 = T(2,1)/r2; s2 = T(3,1)/r2
    V(2,2) = c2   ; V(2,3) = s2
    V(3,2) = -s2  ; V(3,3) = c2
    T(2,1) = r2   ; T(3,1) = 0.d0
    T(2:3,2:4) = matmul(V(2:3,2:3),T(2:3,2:4))
    T(1,2:3)   = T(2:3,1);  T(4,2:3) = T(2:3,4)
    T(2:3,2:3) = matmul(T(2:3,2:3),transpose(V(2:3,2:3)))
else
    c2 = 1.d0; s2 = 0.d0
end if

! Zero out entries T(4,2) and T(2,4)
! compute cos and sin of rotation angle in the 3-4 plane
r3 = pythag(T(4,2), T(3,2))
if(r3 /= 0.0d0) then
    c3 = T(3,2)/r3; s3 = T(4,2)/r3
    V(3,3) = c3   ; V(3,4) = s3
    V(4,3) = -s3   ; V(4,4) = c3
    T(3,2) = r3   ; T(4,2) = 0.0d0
    T(3:4,3:4) = matmul(V(3:4,3:4),T(3:4,3:4))
    T(1:2,3:4) = transpose(T(3:4,1:2))
    T(3:4,3:4) = matmul(T(3:4,3:4),transpose(V(3:4,3:4)))
else
    c3 = 1.0d0; s3 = 0.0d0
end if

! Compute net rotation matrix (accumulate similarity for evec. computation)
! To save transposing later, This is the transpose!
V(1,1) = 1.0d0  ; V(1,2:4) = 0.0d0 ; V(2:4,1) = 0.0d0
V(2,2) = c2     ; V(3,2)   = c1*s2         ; V(4,2) = s1*s2; c1c2 = c1*c2; s1c2=s1*c2
V(2,3) = -s2*c3 ; V(3,3)   = c1c2*c3-s1*s3 ; V(4,3) = s1c2*c3+c1*s3
V(2,4) =  s2*s3 ; V(3,4)   =-c1c2*s3-s1*c3 ; V(4,4) = -s1c2*s3+c1*c3
  
end subroutine givens4
!-------------------------------------------------------------------------------
subroutine svdcmp(mmax, a, m, n, w, v, rv1)
!-------------------------------------------------------------------------------
! TODO: Replace GOTO statements
! Algorithm svdcmp from 'Numerical Recipes' to compute eigenvalues and 
! eigenvectors for SPD matrix B
! Don't ask me what this is :)
!-------------------------------------------------------------------------------
integer :: mmax
integer :: m,n
integer :: i, its, j, jj, k, l, nm
real(dp) :: a(mmax,*), v(mmax,*), w(*), rv1(*)
real(dp) :: anorm, c, f, g, h, s, scale, x, y, z, one

g = 0.0d0
scale = 0.0d0
anorm = 0.0d0
do i = 1, n
    l = i + 1
    rv1(i) = scale*g
    g = 0.0d0
    s = 0.0d0
    scale = 0.0d0
    if (i <= m) then
        do k = i, m
            scale = scale + abs(a(k,i))
        end do
      
        if (scale /= 0.0d0) then
            do k = i, m
                a(k,i) = a(k,i)/scale
                s = s + a(k,i)*a(k,i)
            end do
            f = a(i,i)
            g = -sign(sqrt(s),f)
            h = f*g - s
            a(i,i) = f - g
            do j = l, n 
                s = 0.0d0
                do k = i, m
                    s = s + a(k,i)*a(k,j)
                end do
                f = s/h
                do k = i, m
                    a(k,j) = a(k,j) + f*a(k,i)
                end do
            end do
            do k = i, m
                a(k,i) = scale*a(k,i)
            end do
        end if
    end if
   
    w(i) = scale*g
    g = 0.0d0
    s = 0.0d0
    scale = 0.0d0
   
    if ((i.le.m).and.(i.ne.n))then
        do k = l,n
            scale = scale+abs(a(i,k))
        end do
        if (scale.ne.0.0d0)then
            do k = l,n
                a(i,k) = a(i,k)/scale
                s = s+a(i,k)*a(i,k)
            end do
            f = a(i,l)
            g = -sign(sqrt(s),f)
            h = f*g-s
            a(i,l) = f-g
            do k = l,n
                rv1(k) = a(i,k)/h
            end do
            do j = l,m
                s = 0.0d0
                do k = l,n
                    s = s+a(j,k)*a(i,k)
                end do
                do k = l,n
                    a(j,k) = a(j,k)+s*rv1(k)
                end do
            end do
            do k = l,n
                a(i,k) = scale*a(i,k)
            end do
        end if
    end if
    anorm = max(anorm,(abs(w(i))+abs(rv1(i))))
end do

do i = n, 1, -1
    if (i .lt. n) then
        if (g.ne.0.0d0) then
            do j = l,n
                v(j,i) = (a(i,j)/a(i,l))/g
            end do
            do j = l,n
                s = 0.0d0
                do k = l,n
                    s = s + a(i,k)*v(k,j)
                end do
                do k = l,n
                    v(k,j) = v(k,j) + s*v(k,i)
                end do
            end do
        end if
        do j = l,n
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
        end do
    end if
    v(i,i) = 1.0d0
    g = rv1(i)
    l = i
end do

do i = min(m,n), 1, -1
    l = i + 1
    g = w(i)
    do j = l,n
        a(i,j) = 0.0d0
    end do
    if (g /= 0.0d0) then
        g = 1.0d0/g
        do j = l, n
            s = 0.0d0
            do k = l,m
                s = s + a(k,i)*a(k,j)
            end do
            f = (s/a(i,i))*g
            do k = i, m
                a(k,j) = a(k,j) + f*a(k,i)
            end do
        end do
        do j = i,m
            a(j,i) = a(j,i)*g
        end do
    else
        do j = i,m
            a(j,i) = 0.0d0
        end do
    end if
    a(i,i) = a(i,i) + 1.0d0
end do

do k = n, 1, -1
    do its = 1, 30
        do l = k, 1, -1
            nm = l-1
            if ((abs(rv1(l)) + anorm) == anorm)  goto 2
            if ((abs(w(nm)) + anorm) == anorm)  goto 1 ! It's same as exit
        end do
1       c = 0.0d0
        s = 1.0d0
        do i = l, k
            f = s*rv1(i)
            rv1(i) = c*rv1(i)
            if ((abs(f) + anorm) == anorm) goto 2 ! It's same as exit
            g = w(i)
            h = pythag(f,g)
            w(i) = h
            h = 1.0d0/h
            c = g*h
            s = -f*h
            do j = 1, m     
                y = a(j,nm)
                z = a(j,i)
                a(j,nm) = y*c + z*s
                a(j,i) = -y*s + z*c
            end do
        end do
2       z = w(k)
        if (l == k) then
            if (z < 0.0d0) then
                w(k) = -z
                do j = 1, n
                    v(j,k) = -v(j,k)
                end do
            end if
            goto 3  ! It's same as exit
        end if
        if (its == 30) then
            call terminate_with_error('ERROR: no convergence in svdcmp')
        end if
      
        x = w(l)
        nm = k-1
        y = w(nm)
        g = rv1(nm)
        h = rv1(k)
        f = ((y-z)*(y+z)+(g-h)*(g+h)) / (2.0d0*h*y)
        one =  1.0
        g = pythag(f,one)
        f = ((x-z)*(x+z) + h*((y/(f + sign(g,f))) - h))/x
        c = 1.0d0
        s = 1.0d0
        do j = l, nm
            i = j + 1
            g = rv1(i)
            y = w(i)
            h = s*g
            g = c*g    
            z = pythag(f,h)
            rv1(j) = z
            c = f/z
            s = h/z
            f = x*c + g*s
            g = -x*s + g*c
            h = y*s
            y = y*c
            do jj = 1,n
                x = v(jj,j)
                z = v(jj,i)
                v(jj,j) = x*c + z*s
                v(jj,i) = -x*s + z*c
            end do
            z = pythag(f,h)
            w(j) = z
            if (z /= 0.0d0) then
                z = 1.0d0/z
                c = f*z
                s = h*z
            end if
            f = c*g + s*y
            x = -s*g + c*y
            do jj = 1,m
                y = a(jj,j)
                z = a(jj,i)
                a(jj,j) = y*c + z*s
                a(jj,i) = -y*s + z*c
            end do
        end do
        rv1(l) = 0.0d0
        rv1(k) = f       
        w(k) = x
    end do
3   continue
end do
  
end subroutine svdcmp
!-------------------------------------------------------------------------------
function pythag(a,b)
!-------------------------------------------------------------------------------
! Calculate c from given a and b based on Pythagorean theorem (c**2 = a**2 + b**2)
!-------------------------------------------------------------------------------
real(dp) :: a, b, pythag
real(dp) :: absa, absb

absa = abs(a)
absb = abs(b)
if (absa > absb) then
    pythag = absa*dsqrt(1.0d0+(absb/absa)**2)
else
    if (absb == 0.0d0) then
        pythag = 0.0d0
    else
        pythag = absb*dsqrt(1.0d0+(absa/absb)**2)
    end if
end if

end function pythag
!-------------------------------------------------------------------------------
END MODULE RMSD
!-------------------------------------------------------------------------------
