!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/pair_list.f90
!
! Description:
!  
!-------------------------------------------------------------------------------
MODULE PAIR_LIST
!-------------------------------------------------------------------------------

use globals
use in_out_utils, only: find_atom_idx
use energy_vars

implicit none
save
private

type position_int
    integer, dimension(3) :: r
end type position_int

type grid_hash_index
    integer :: first
    integer :: second
end type grid_hash_index

type position_index_pair
    integer :: i_atm
    real(dp) :: R(3)
end type position_index_pair

real(dp), dimension(:,:), allocatable :: R_sort
type(grid_hash_index), allocatable :: grid_particle_hash_index(:)

integer, parameter :: num_grid = 32, num_grid_half = num_grid/2
integer, parameter :: num_grid_cells = num_grid**3
integer :: cell_start(num_grid_cells), cell_end(num_grid_cells)
real(dp) :: cell_size(3)
real(dp) :: grid_origin(3)

type(position_int) :: init_grid_size
type(position_int) :: grid_size

public :: initialize_pair_list
public :: finalize_pair_list
public :: fill_pair_index
public :: pairlist_update
public :: pairlist_update_optSC
public :: get_pair_type
public :: position_index_pair

CONTAINS
!===============================================================================
! Initialize/finalize pair_list related variables
!===============================================================================
subroutine initialize_pair_list()
!-------------------------------------------------------------------------------
cell_size(:) = LRoff
grid_size%r(1:3) = num_grid

allocate(i_P(tn%atom))
allocate(R_sort(3,tn%atom))
allocate(grid_particle_hash_index(tn%atom))

end subroutine initialize_pair_list
!-------------------------------------------------------------------------------
subroutine finalize_pair_list()
!-------------------------------------------------------------------------------
deallocate(i_P)
deallocate(R_sort)
deallocate(grid_particle_hash_index)

end subroutine finalize_pair_list
!-------------------------------------------------------------------------------
subroutine get_connected_atom_info(protein, c_atm, cn)
!-------------------------------------------------------------------------------
! Get information about how atoms are connected.
! "cn(i_atm) = n_atm" means that the number of atoms connected to i_atm is n_atm.
! "c_atm(i_cn, i_atm) = j_atm" means i_cn th connected atom of i_atm is j_atm.
!  i_cn in above statement can have a value between 1 and cn(i_atm).
! This connected atom information is used to fill pair_index type by checking
! whether i_atm and j_atm are 1-2, 1-3, 1-4 connected or not.
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer, intent(out) :: c_atm(:,:)
integer, intent(out) :: cn(:)
integer :: atm_no, res_no, het_no, lig_no, ref_res_no
integer :: conn_no
integer :: i, i_atm, j_atm, i_bnd, i_link
integer :: start_bnd, min_idx, max_idx
integer :: indx(2), bnd_idx(2)
character(len=4) :: atom_name
character(len=6) :: error_mode
  
atm_no = 0
do res_no = 1, protein%n_res + protein%n_het
    if (res_no <= protein%n_res) then
        ref_res_no = protein%residue(res_no)%res_type
    else
        het_no = res_no - protein%n_res
        ref_res_no = protein%hetmol(het_no)%res_type
    end if

    do i_atm = 1, ref_res(ref_res_no)%n_atm
        atm_no = atm_no + 1
        conn_no = 0 ! initialize no of connected atoms
        if (res_no == 1 .or. res_no > protein%n_res) then !first residue or hetmol
            start_bnd = 1
        else
            start_bnd = -2
        end if
        do i_bnd = start_bnd, ref_res_eng(ref_res_no)%n_bnd_E
            bnd_idx(1:2) = ref_res_eng(ref_res_no)%atm_in_bnd_E(1:2,i_bnd)
            indx(1:2) = res_index(res_no)%loc_no(bnd_idx(1:2))
            do i = 1, 2
                if (indx(i) == atm_no) then
                conn_no = conn_no + 1
                c_atm(conn_no,atm_no) = indx(3-i) ! bonded atm
                end if
            end do
        end do
        cn(atm_no) = conn_no ! no of bonded atm
    end do
end do

! To consider disulfide bond
do i_link = 1, protein%n_link
    if (protein%link(i_link)%link_name /= 'DISU') cycle
    indx(:) = 0
    atom_name = 'SG'
    do i = 1, 2
        res_no = protein%link(i_link)%link_res_no(i)
        call find_atom_idx(res_no, protein%residue(res_no)%res_type, atom_name, &
                           atm_no, error_mode)
        indx(i) = ii_R(atm_no,res_no)
    end do
    min_idx = minval(indx(1:2))
    max_idx = maxval(indx(1:2))
    cn(min_idx) = cn(min_idx) + 1
    c_atm(cn(min_idx),min_idx) = max_idx
    cn(max_idx) = cn(max_idx) + 1
    c_atm(cn(max_idx),max_idx) = min_idx
end do

do lig_no = 1, protein%n_lig
    ref_res_no = protein%ligand(lig_no)%lig_type
    do atm_no = 1, ref_lig(ref_res_no)%n_atm
        i_atm = ii_L(atm_no, lig_no)
        do i_bnd = 1, ref_lig(ref_res_no)%n_bnd
            if (ref_lig(ref_res_no)%bnd(2,i_bnd) == atm_no) then
                j_atm = ii_L(ref_lig(ref_res_no)%bnd(3,i_bnd), lig_no)
                cn(i_atm) = cn(i_atm) + 1
                c_atm(cn(i_atm),i_atm) = j_atm
                cn(j_atm) = cn(j_atm) + 1
                c_atm(cn(j_atm),j_atm) = i_atm
            else if (ref_lig(ref_res_no)%bnd(3,i_bnd) == atm_no) then
                j_atm = ii_L(ref_lig(ref_res_no)%bnd(2,i_bnd), lig_no)
                cn(i_atm) = cn(i_atm) + 1
                c_atm(cn(i_atm),i_atm) = j_atm
                cn(j_atm) = cn(j_atm) + 1
                c_atm(cn(j_atm),j_atm) = i_atm
            end if
        end do
    end do
end do

end subroutine get_connected_atom_info
!-------------------------------------------------------------------------------
subroutine fill_pair_index(protein)
!-------------------------------------------------------------------------------
! Fill n_bnd, i_bnd information in pair_index type
! i_P(i_atm)%n_bnd(1): # of bonds which 1-2 connected.
! i_P(i_atm)%n_bnd(2): # of bonds which 1-3 connected.
! i_P(i_atm)%n_bnd(3): # of bonds which 1-4 connected.
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer :: c_atm(12,tn%atom), cn(tn%atom)
integer :: i_atm, j_atm, conn_12, conn_13, conn_14
integer :: i_cn, j_cn, k_cn
logical :: excluded, nb14

! Get coonected atom information to generate pair index.
cn(:) = 0
call get_connected_atom_info(protein, c_atm, cn)

! initialize
do i_atm = 1, tn%atom
    i_P(i_atm)%n_bnd(:) = 0
    i_P(i_atm)%i_bnd(:,:) = 0
end do

do i_atm = 1, tn%atom - 1
    do j_atm = i_atm + 1, tn%atom
        excluded = .false.
        nb14 = .false.
        ! Start with connected atom to j_atm
        do i_cn = 1, cn(j_atm)
            ! conn_12 is atom index for i_cn th connected atom to j_atm
            ! j_atm - conn_12
            conn_12 = c_atm(i_cn, j_atm)
            if (conn_12 == i_atm) then ! i_atm is connected to j_atm.
                                       ! It's 1-2 pair
                excluded = .true.
                i_P(i_atm)%n_bnd(1) = i_P(i_atm)%n_bnd(1) + 1
                i_P(i_atm)%i_bnd(i_P(i_atm)%n_bnd(1),1) = j_atm
                i_P(j_atm)%n_bnd(1) = i_P(j_atm)%n_bnd(1) + 1
                i_P(j_atm)%i_bnd(i_P(j_atm)%n_bnd(1),1) = i_atm
                exit
            end if
            ! search connected atom to conn_12
            ! j_atm - conn_12 - conn_13
            do j_cn = 1, cn(conn_12)
                conn_13 = c_atm(j_cn, conn_12)
                if (conn_13 == j_atm) then
                    cycle
                end if
                if (conn_13 == i_atm) then ! i_atm & j_atm : 1-3 pair!
                    excluded = .true.
                    i_P(i_atm)%n_bnd(2) = i_P(i_atm)%n_bnd(2) + 1
                    i_P(i_atm)%i_bnd(i_P(i_atm)%n_bnd(2),2) = j_atm
                    i_P(j_atm)%n_bnd(2) = i_P(j_atm)%n_bnd(2) + 1
                    i_P(j_atm)%i_bnd(i_P(j_atm)%n_bnd(2),2) = i_atm
                    exit
                end if
                ! search connected atom to conn_13
                ! j_atm - conn_12 - conn_13 - conn_14
                do k_cn = 1, cn(conn_13)
                    conn_14 = c_atm(k_cn, conn_13)
                    if (conn_14 == i_atm) then ! i_atm & j_atm : 1-4 pair
                        nb14 = .true.
                    end if
                end do
            end do
            if (excluded) exit
        end do
        if (.not. excluded) then
            if (nb14) then ! i_atm & j_atm : 1-4 pair
                i_P(i_atm)%n_bnd(3) = i_P(i_atm)%n_bnd(3) + 1
                i_P(i_atm)%i_bnd(i_P(i_atm)%n_bnd(3),3) = j_atm
                i_P(j_atm)%n_bnd(3) = i_P(j_atm)%n_bnd(3) + 1
                i_P(j_atm)%i_bnd(i_P(j_atm)%n_bnd(3),3) = i_atm
            end if
        end if
    end do 
end do 

end subroutine fill_pair_index
!===============================================================================
! Update pair_list
!===============================================================================
subroutine pairlist_update()
!-------------------------------------------------------------------------------
! Update pair_list, i_P, using current structure (R array) using cell lists.
! Please refer Wikipedia about cell lists
!-------------------------------------------------------------------------------
grid_origin(1) = sum(R(1,:))/tn%atom
grid_origin(2) = sum(R(2,:))/tn%atom
grid_origin(3) = sum(R(3,:))/tn%atom
grid_origin(:) = grid_origin(:) - cell_size(:)*num_grid_half

cell_start(:) = -1
cell_end(:) = -1

call update_grid(R, R_sort, grid_particle_hash_index, grid_size, grid_origin,&
                 cell_size, cell_start, cell_end)
call update_neighbors(tn%atom, R_sort, grid_size, grid_origin, cell_size,&
                      cell_start, cell_end, grid_particle_hash_index)

end subroutine pairlist_update
!-------------------------------------------------------------------------------
subroutine pairlist_update_optSC(tna_i, tna_j, R_i, R_j, use_self)
!-------------------------------------------------------------------------------
! TODO: Comments
!-------------------------------------------------------------------------------
integer, intent(in) :: tna_i, tna_j
type(position_index_pair), intent(in) :: R_i(tna_i), R_j(tna_j)
logical, intent(in) :: use_self
integer :: i, j, ii, np
real(dp) :: dR(3), Ri(3), dsq

do i = 1, tn%atom
    i_P(i)%pair_end_index(1:3) = 0
    i_P(i)%n_pair = 0
    i_P(i)%n_Lpair = 0
end do

if (use_self) then
    do i = 1, tna_i-1
        ii = R_i(i)%i_atm
        np = i_P(ii)%n_pair
        Ri(1:3) = R_i(i)%R(1:3)

        do j = i+1, tna_i
            dR(1:3) = R_i(j)%R(1:3) - Ri(1:3)
            dsq = dot_product(dR,dR)
            if (dsq < Roff_sqr) then
                np = np + 1
                i_P(ii)%d(np) = sqrt(dsq)
                i_P(ii)%i_pair(np) = R_i(j)%i_atm
            end if
        end do
        i_P(ii)%n_pair = np
    end do
end if

do i = 1, tna_i
    ii = R_i(i)%i_atm
    np = i_P(ii)%n_pair
    Ri(1:3) = R_i(i)%R(1:3)

    do j = 1, tna_j
        dR(1:3) = R_j(j)%R(1:3) - Ri(1:3)
        dsq = dot_product(dR,dR)
        if (dsq < Roff_sqr) then
            np = np + 1
            i_P(ii)%d(np) = sqrt(dsq)
            i_P(ii)%i_pair(np) = R_j(j)%i_atm
        end if
    end do

    i_P(ii)%n_pair = np
    call sort_by_pair_type(i_P(ii))
end do

end subroutine pairlist_update_optSC
!-------------------------------------------------------------------------------
subroutine update_grid(pos, sorted_pos, grid_particle_hash_index, grid_size,&
                       origin, cell_size, cell_start, cell_end)
!-------------------------------------------------------------------------------
! TODO: Comments
!-------------------------------------------------------------------------------
real(dp), intent(in) :: pos(:,:), origin(3), cell_size(3)
real(dp), intent(inout) :: sorted_pos(:,:)
type(grid_hash_index), intent(inout) :: grid_particle_hash_index(:)
type(position_int), intent(in) :: grid_size
integer, intent(inout) :: cell_start(:), cell_end(:)

call calc_hash(grid_particle_hash_index, pos, grid_size, origin, cell_size)
call sort_hash(tn%atom, grid_particle_hash_index)
call reorder_data_find_cell_start(cell_start, cell_end, sorted_pos, &
                                  grid_particle_hash_index, pos)

end subroutine update_grid
!-------------------------------------------------------------------------------
subroutine update_neighbors(n_atm, old_pos, grid_size, origin, cell_size,&
                           cell_start, cell_end, grid_particle_hash_index)
!-------------------------------------------------------------------------------
! TODO: Comments
!-------------------------------------------------------------------------------
real(dp), intent(in) :: old_pos(:,:), origin(3), cell_size(3)
type(position_int), intent(in) :: grid_size
integer, intent(in) :: n_atm, cell_start(:), cell_end(:)
type(grid_hash_index), intent(in) :: grid_particle_hash_index(:)

real(dp) :: pos(3)
type(position_int) :: grid_pos, neighbor_pos
integer :: i, ii, x, y, z, grid_hash
integer :: ip1, ip2

do i = 1, n_atm
    call calc_grid_pos(old_pos(1:3,i), origin, cell_size, grid_pos)
    call calc_grid_hash(grid_pos, grid_size, grid_hash)
end do

do i = 1, n_atm
    pos = old_pos(1:3,i)
    ii = grid_particle_hash_index(i)%second

    ip1 = 0
    ip2 = 0
    call calc_grid_pos(pos, origin, cell_size, grid_pos)

    do z = -1, 1
        do y = -1, 1
            do x = -1, 1
                neighbor_pos%r = grid_pos%r + (/x,y,z/)
                call calc_grid_hash(neighbor_pos, grid_size, grid_hash)
                call check_cell(grid_hash, i, ii, pos, old_pos, cell_start, &
                                cell_end, ip1, ip2, grid_particle_hash_index)
            end do
        end do
    end do
    i_P(ii)%n_pair = ip1
    i_P(ii)%n_Lpair = ip2
    call sort_by_pair_type(i_P(ii))
end do

end subroutine update_neighbors
!-------------------------------------------------------------------------------
subroutine calc_hash(grid_particle_hash_index, pos, grid_size, origin, cell_size)
!-------------------------------------------------------------------------------
! TODO: Comments
!-------------------------------------------------------------------------------
type(grid_hash_index), intent(inout) :: grid_particle_hash_index(:)
real(dp), intent(in) :: pos(:,:), origin(3), cell_size(3)
type(position_int), intent(in) :: grid_size

type(position_int) :: grid_pos
integer :: i, hash

do i = 1, tn%atom
    call calc_grid_pos(pos(1:3,i), origin, cell_size, grid_pos)
    call calc_grid_hash(grid_pos, grid_size, hash)
    grid_particle_hash_index(i)%first = hash
    grid_particle_hash_index(i)%second = i
end do

end subroutine calc_hash
!-------------------------------------------------------------------------------
subroutine calc_grid_pos(p, origin, cell_size, grid_pos)
!-------------------------------------------------------------------------------
! TODO: Comments
!-------------------------------------------------------------------------------
real(dp), intent(in) :: p(3), origin(3), cell_size(3)
type(position_int), intent(out) :: grid_pos

grid_pos%r = (p(1:3) - origin(1:3)) / cell_size(1:3) + 1

end subroutine calc_grid_pos
!-------------------------------------------------------------------------------
subroutine calc_grid_hash(grid_pos, grid_size, hash)
!-------------------------------------------------------------------------------
! TODO: Comments
!-------------------------------------------------------------------------------
type(position_int), intent(in) :: grid_pos, grid_size
integer, intent(out) :: hash
type(position_int) :: pos

! changed to use standard function: "iand" is same in & operator in C
!pos%r = (grid_pos%r-1) .and. (grid_size%r -1)
pos%r = iand(grid_pos%r-1, grid_size%r -1)
hash = pos%r(3) * grid_size%r(2)* grid_size%r(1) + &
       pos%r(2) * grid_size%r(1) + &
       pos%r(1) + 1

end subroutine calc_grid_hash
!-------------------------------------------------------------------------------
subroutine sort_hash(n,hash)
!-------------------------------------------------------------------------------
! TODO: Comments
!-------------------------------------------------------------------------------
integer, intent(in) :: n
type(grid_hash_index), intent(inout) :: hash(n)

integer :: i, j, k, index
type(grid_hash_index) :: tmp

k = n/2 + 1
index = n
do while (n > 1)
    if (k > 1) then
        k = k - 1
        tmp = hash(k)
    else
        tmp = hash(index)
        hash(index) = hash(1)
        index = index - 1
        if (index <= 1) then
            hash(1) = tmp
            return
        end if
    end if

    i = k
    j = k + k
   
    do while (j <= index)
        if (j < index) then
            if (cmp_hash(hash(j), hash(j+1))) j = j+1
        end if
        if (cmp_hash(tmp, hash(j))) then
            hash(i) = hash(j)
            i = j
            j = j + j
        else
            j = index + 1
        end if
    end do
    hash(i) = tmp
end do

end subroutine sort_hash
!-------------------------------------------------------------------------------
logical function cmp_hash(X,Y)
!-------------------------------------------------------------------------------
! TODO: Comments
!-------------------------------------------------------------------------------
type(grid_hash_index), intent(in) :: X, Y

if (X%first < Y%first) then
    cmp_hash = .true.
else if (X%first == Y%first .and. X%second < Y%second) then
    cmp_hash = .true.
else
    cmp_hash = .false.
end if

end function cmp_hash
!-------------------------------------------------------------------------------
subroutine reorder_data_find_cell_start(cell_start, cell_end, sorted_pos, &
                                        grid_particle_hash_index, old_pos)
!-------------------------------------------------------------------------------
! TODO: Comments
!-------------------------------------------------------------------------------
integer, intent(inout) :: cell_start(:), cell_end(:)
real(dp), intent(inout) :: sorted_pos(:,:)
type(grid_hash_index), intent(in) :: grid_particle_hash_index(:)
real(dp), intent(in) :: old_pos(:,:)

integer :: i, hash, sorted_index

do i = 1, tn%atom
    hash = grid_particle_hash_index(i)%first

    if (i == 1) then
        cell_start(hash) = i
    else if (hash /= grid_particle_hash_index(i-1)%first) then
        cell_start(hash) = i
        if (i > 1) cell_end(grid_particle_hash_index(i-1)%first) = i-1
    end if
    if (i == tn%atom) cell_end(hash) = i

    sorted_index = grid_particle_hash_index(i)%second
    sorted_pos(1:3,i) = old_pos(1:3,sorted_index)
end do

end subroutine reorder_data_find_cell_start
!-------------------------------------------------------------------------------
subroutine check_cell(grid_hash, i_atm, ii, pos1, old_pos, cell_start, &
                      cell_end, ip1, ip2, grid_particle_hash_index)
!-------------------------------------------------------------------------------
! TODO: Comments
!-------------------------------------------------------------------------------
integer, intent(in) :: grid_hash, i_atm, ii, cell_start(:), cell_end(:)
real(dp), intent(in) :: pos1(3), old_pos(:,:)
integer, intent(inout) :: ip1, ip2
type(grid_hash_index), intent(in) :: grid_particle_hash_index(:)
integer :: start_idx, end_idx, j
real(dp) :: dR(3), dsq

start_idx = cell_start(grid_hash)
if (start_idx == -1) return

end_idx = min(cell_end(grid_hash), i_atm-1)

do j = start_idx, end_idx
    dR(1:3) = pos1(1:3) - old_pos(1:3,j)
    dsq = dot_product(dR,dR)

    if (dsq < LRoff_sqr) then
        if (dsq < Roff_sqr) then
            ip1 = ip1 + 1
            i_P(ii)%d(ip1) = sqrt(dsq)
            i_P(ii)%i_pair(ip1) = grid_particle_hash_index(j)%second
        else
            ip2 = ip2+1
            i_P(ii)%i_Lpair(ip2) = grid_particle_hash_index(j)%second
        end if
    end if
end do

end subroutine check_cell
!-------------------------------------------------------------------------------
subroutine sort_by_pair_type(pair)
!-------------------------------------------------------------------------------
! This subroutine sorts pair information by pair type in order of
! 1-2 -> 1-3 -> 1-4 -> other pairs.
! It replaces function of pair_type array.
!-------------------------------------------------------------------------------
type(pair_index_type), intent(inout) :: pair
integer :: i, j, k, l, n_bnd
integer :: i_pair
real(dp) :: d

l = 0
do i = 1, 3
    n_bnd = pair%n_bnd(i)

    do j = 1+l, pair%n_pair
        do k = 1, n_bnd
            if (pair%i_pair(j) == pair%i_bnd(k,i)) then
                l = l + 1

                i_pair = pair%i_pair(l)
                pair%i_pair(l) = pair%i_pair(j)
                pair%i_pair(j) = i_pair

                d = pair%d(l)
                pair%d(l) = pair%d(j)
                pair%d(j) = d

                exit
            end if
        end do
    end do
    pair%pair_end_index(i) = l
end do

end subroutine sort_by_pair_type
!===============================================================================
! Pair_type related
!===============================================================================
subroutine get_pair_type(i_atm, j_atm, pair_type)
!-------------------------------------------------------------------------------
! TODO: Comments
!-------------------------------------------------------------------------------
integer, intent(in)  :: i_atm, j_atm
integer, intent(out) :: pair_type
integer :: j

if (i_atm == j_atm) then
    pair_type = 0
    return
end if

do pair_type = 1, 3
    do j = 1, i_P(i_atm)%n_bnd(pair_type)
        if (i_P(i_atm)%i_bnd(j,pair_type) == j_atm) return
    end do
end do
pair_type = 4

end subroutine get_pair_type
!-------------------------------------------------------------------------------
END MODULE PAIR_LIST
!-------------------------------------------------------------------------------
