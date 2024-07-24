!-------------------------------------------------------------------------------
! Copyright (C) 2021, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/knowledge/ml_energy.f90
!
! Description: Machine-learning-based scoring function
!              for monomer structure prediction
!-------------------------------------------------------------------------------
MODULE ML_SCORE
!-------------------------------------------------------------------------------
use globals
use logger, only: log_p, terminate_with_error
use sort, only: sort2
use convert_res_name, only: convert_to_stdres
!
use energy_vars, only: R, reinitialize

implicit none
private

!-------------------------------------------------------------------------------
! PARAMETERS
!-------------------------------------------------------------------------------
type weight_table
!-------------------------------------------------------------------------------
real, allocatable :: w(:,:)
real, allocatable :: wt(:,:)
real, allocatable :: b(:)

end type weight_table
!-------------------------------------------------------------------------------
type amino_acid_table
!-------------------------------------------------------------------------------
real :: f(8)
character(len=4) :: res_name

end type amino_acid_table
!-------------------------------------------------------------------------------

integer, parameter :: N_NEIGHBOR = 40
integer, parameter :: n_layer = 11
type(weight_table) :: ml_param(n_layer)
type(amino_acid_table) :: ml_aa_table(21)
real, allocatable :: amino_param(:,:) ! dim = 8, n_res
integer :: idx_NA
logical, allocatable :: is_glycine(:)
integer, allocatable :: bb_id(:,:) ! dim = 5, n_res

public :: initialize_ml_score
public :: finalize_ml_score
public :: calc_ml_score

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_ml_score(protein)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
character(len=len_fname) :: infile_ml_param

infile_ml_param = trim(data_dir) // 'score_ml.data'
!
if (.not. reinitialize) then
    write(log_msg, "(A)") "  Reading ML parameter table"
    call log_p(log_msg, me=me, level=40)
    call read_ml_param(infile_ml_param)
end if
!
call setup_protein_for_ml(protein)

end subroutine initialize_ml_score
!-------------------------------------------------------------------------------
subroutine finalize_ml_score()
!-------------------------------------------------------------------------------
integer :: i

deallocate(amino_param)
deallocate(is_glycine)
deallocate(bb_id)

if (.not. reinitialize) then
    do i = 1, n_layer
        deallocate(ml_param(i)%w)
        deallocate(ml_param(i)%b)
    end do
end if
    
end subroutine finalize_ml_score
!-------------------------------------------------------------------------------
subroutine calc_ml_score(f, g, appl_res, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f, g(3,tn%atom)
logical, intent(in) :: appl_res(tn%residue), calc_g

integer :: i_res, j_res, n_res, i_neigh, i_atm, j_atm, n_atm, m_atm, i
real(dp) :: dist_sq_ca(tn%residue, tn%residue)
real(dp) :: dist, dr(3), dr_s(3, 5, 5, N_NEIGHBOR), dfdr(3)
real :: f_tmp, dfdd(25, N_NEIGHBOR)
integer :: neigh_idx(tn%residue), n_neigh
real :: input(8+8+26, N_NEIGHBOR)
real :: R_bb(3,5,tn%residue)
integer, parameter :: idx_ca = 1

f = 0.0d0
g = 0.0d0

n_res = tn%residue

! backbone coordinates
do i_res = 1, n_res
    if (.not. appl_res(i_res)) cycle
    if (is_glycine(i_res)) then
        n_atm = 4
    else
        n_atm = 5
    end if
    do i_atm = 1, n_atm
        R_bb(:,i_atm,i_res) = R(:,bb_id(i_atm,i_res))
    end do
end do

! calculate CA dist
dist_sq_ca = 1.0d10
do i_res = 2, n_res
    if (.not. appl_res(i_res)) cycle
    do j_res = 1, i_res - 1
        if (.not. appl_res(j_res)) cycle
        dr = R_bb(:,idx_ca,i_res) - R_bb(:,idx_ca,j_res)
        dist = dot_product(dr, dr)
        dist_sq_ca(i_res, j_res) = dist
        dist_sq_ca(j_res, i_res) = dist
    end do
end do

! energy calculation
do i_res = 1, n_res
    if (.not. appl_res(i_res)) cycle
    if (is_glycine(i_res)) then
        n_atm = 4
    else
        n_atm = 5
    end if
    input = 0.0

    do i = 1, N_NEIGHBOR
        input(1:8,i) = amino_param(1:8, i_res)
    end do
    !
    ! sort dist and get neighbor list
    call sort2(n_res, dist_sq_ca(1:n_res, i_res), neigh_idx(1:n_res))
    !
    n_neigh = 0
    do i_neigh = 1, N_NEIGHBOR
        if (i_neigh >= n_res .or. dist_sq_ca(i_neigh,i_res) > 1.0d5) exit
        n_neigh = i_neigh
        !
        j_res = neigh_idx(i_neigh)
        input(9:16, i_neigh) = amino_param(1:8, j_res)
        if (is_glycine(j_res)) then
            m_atm = 4
        else
            m_atm = 5
        end if
        !
        do j_atm = 1, m_atm
            do i_atm = 1, n_atm
                dr = R_bb(:,j_atm,j_res) - R_bb(:,i_atm,i_res)
                dist = sqrt(dot_product(dr, dr))
                input(11+j_atm*5+i_atm, i_neigh) = dist
                if (calc_g) &
                    dr_s(:,i_atm,j_atm,i_neigh) = dr / dist
            end do
        end do
        !
        input(42, i_neigh) = j_res - i_res
    end do
    !
    call pass_layer(f_tmp, dfdd, input, n_neigh, calc_g)
    f = f + f_tmp
    if (.not. calc_g) cycle
    !
    do i_neigh = 1, n_neigh
        j_res = neigh_idx(i_neigh)
        if (is_glycine(j_res)) then
            m_atm = 4
        else
            m_atm = 5
        end if
        !
        do j_atm = 1, m_atm
            do i_atm = 1, n_atm
                dfdr = dr_s(:,i_atm,j_atm,i_neigh) * dfdd((j_atm-1)*5+i_atm,i_neigh)
                g(:,bb_id(i_atm,i_res)) = g(:,bb_id(i_atm,i_res)) - dfdr
                g(:,bb_id(j_atm,j_res)) = g(:,bb_id(j_atm,j_res)) + dfdr
            end do
        end do
    end do    
end do
        
end subroutine calc_ml_score
!-------------------------------------------------------------------------------
subroutine pass_layer(fout, gout, input, n_neigh, calc_g)
! Architecture: input(42, n_neigh) -> fc2_0 -> (n1, n_neigh)
!                                  -> fc2_1 -> fc2_2 -> fc2_3 -> (n1, n_neigh)
!                                  -> fc2_4 -> (n12, n_neigh)
!                                  -> max axis=2 -> (n12)
!                                  -> fc3_0 -> (n2)
!                                  -> fc2_1 -> fc2_2 -> fc2_3 -> fc2_4 -> (n2)
!                                  -> fc4_0 -> (1)
!   activation: ReLU
!-------------------------------------------------------------------------------
real, intent(in) :: input(42, N_NEIGHBOR)
integer, intent(in) :: n_neigh
real, intent(out) :: fout, gout(25, N_NEIGHBOR)
logical, intent(in) :: calc_g

integer, parameter :: n1 = 128
real :: f2(n1, n_neigh)
integer, parameter :: n12 = 256
real :: f24(n12, n_neigh)
integer, parameter :: n2 = 128
real :: f30(n12)
real :: f3(n2)

integer :: i, j, k

! variables for gradient calculation
logical :: l20(n1, n_neigh)
logical :: l21(n1, n_neigh)
logical :: l22(n1, n_neigh)
logical :: l23(n1, n_neigh)
logical :: l24(n12)
integer :: max_sel(n12)
logical :: l30(n2)
logical :: l31(n2)
logical :: l32(n2)
logical :: l33(n2)
logical :: l34(n2)
real :: tmp(n1, 2)
real :: tmp2(n2, 2)
real :: g3(n12)
logical :: flag_pass

! fc2_0
f2 = matmul(ml_param(1)%w, input(:,1:n_neigh))
do i = 1, n_neigh
    f2(:,i) = f2(:,i) + ml_param(1)%b
end do
f2 = max(f2, 0.0)
if (calc_g) l20 = (f2 > 0.0)
! fc2_1
f2 = matmul(ml_param(2)%w, f2)
do i = 1, n_neigh
    f2(:,i) = f2(:,i) + ml_param(2)%b
end do
f2 = max(f2, 0.0)
if (calc_g) l21 = (f2 > 0.0)
! fc2_2
f2 = matmul(ml_param(3)%w, f2)
do i = 1, n_neigh
    f2(:,i) = f2(:,i) + ml_param(3)%b
end do
f2 = max(f2, 0.0)
if (calc_g) l22 = (f2 > 0.0)
! fc2_3
f2 = matmul(ml_param(4)%w, f2)
do i = 1, n_neigh
    f2(:,i) = f2(:,i) + ml_param(4)%b
end do
f2 = max(f2, 0.0)
if (calc_g) l23 = (f2 > 0.0)
! fc2_4
f24 = matmul(ml_param(5)%w, f2)
!
f30 = max(maxval(f24, dim=2) + ml_param(5)%b, 0.0)
if (calc_g) then
    max_sel = maxloc(f24, dim=2)
    l24 = (f30 > 0.0)
end if
! fc3_0
f3 = matmul(ml_param(6)%w, f30) + ml_param(6)%b
f3 = max(f3, 0.0)
if (calc_g) l30 = (f3 > 0.0)
! fc3_1
f3 = matmul(ml_param(7)%w, f3) + ml_param(7)%b
f3 = max(f3, 0.0)
if (calc_g) l31 = (f3 > 0.0)
! fc3_2
f3 = matmul(ml_param(8)%w, f3) + ml_param(8)%b
f3 = max(f3, 0.0)
if (calc_g) l32 = (f3 > 0.0)
! fc3_3
f3 = matmul(ml_param(9)%w, f3) + ml_param(9)%b
f3 = max(f3, 0.0)
if (calc_g) l33 = (f3 > 0.0)
! fc3_4
f3 = matmul(ml_param(10)%w, f3) + ml_param(10)%b
f3 = max(f3, 0.0)
if (calc_g) l34 = (f3 > 0.0)
! fc4_0
f3(1:1) = matmul(ml_param(11)%w, f3) !+ ml_param(11)%b
fout = f3(1)

if (.not. calc_g) return

! gradient calculation
tmp2 = 0.0
! g34 * g40
do k = 1, n2
    if (l34(k)) &
        tmp2(:,1) = tmp2(:,1) + ml_param(10)%wt(:,k) * ml_param(11)%w(1,k)
end do
! g33 * now
do k = 1, n2
    if (l33(k)) &
        tmp2(:,2) = tmp2(:,2) + ml_param(9)%wt(:,k) * tmp2(k,1)
end do
! g32 * now
tmp2(:,1) = 0.0
do k = 1, n2
    if (l32(k)) &
        tmp2(:,1) = tmp2(:,1) + ml_param(8)%wt(:,k) * tmp2(k,2)
end do
! g31 * now
tmp2(:,2) = 0.0
do k = 1, n2
    if (l31(k)) &
        tmp2(:,2) = tmp2(:,2) + ml_param(7)%wt(:,k) * tmp2(k,1)
end do
! g30 * now
g3 = 0.0
do k = 1, n2
    if (l30(k)) &
        g3 = g3 + ml_param(6)%wt(:,k) * tmp2(k,2)
end do

gout = 0.0
do i = 1, n_neigh
    tmp = 0.0
    ! g24 * g3
    flag_pass = .true.
    do j = 1, n12
        if (l24(j) .and. i == max_sel(j)) then
            tmp(:,2) = tmp(:,2) + ml_param(5)%wt(:,j) * g3(j)
            flag_pass = .false.
        end if
    end do
    if (flag_pass) cycle
    ! g23 * now
    do k = 1, n1
        if (l23(k,i)) &
            tmp(:,1) = tmp(:,1) + ml_param(4)%wt(:,k) * tmp(k,2)
    end do
    tmp(:,2) = 0.0
    ! g22 * now
    do k = 1, n1
        if (l22(k,i)) &
            tmp(:,2) = tmp(:,2) + ml_param(3)%wt(:,k) * tmp(k,1)
    end do
    ! g21 * now
    tmp(:,1) = 0.0
    do k = 1, n1
        if (l21(k,i)) &
            tmp(:,1) = tmp(:,1) + ml_param(2)%wt(:,k) * tmp(k,2)
    end do
    ! g20 * now
    do k = 1, n1
        if (l20(k,i)) &
            gout(:,i) = gout(:,i) + ml_param(1)%wt(17:41,k) * tmp(k,1)
    end do
end do

end subroutine pass_layer
!-------------------------------------------------------------------------------
subroutine read_ml_param(infile_ml_param)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: infile_ml_param

character(len=len_fname) :: line
integer :: f_unit, io
integer :: i_res, i_layer, mode, i1, i2, n1, n2

f_unit = 42
open(f_unit, file=trim(infile_ml_param), iostat=io)
if (io < 0) then
    call terminate_with_error("Error: No such file for infile_ml_param")
end if

i_res = 0
i_layer = 0
do
    read(f_unit, "(A120)", iostat=io) line
    if (io < 0) exit

    if (line(1:1) == '@') then
        if (line(3:9) == 'RESIDUE') then
            mode = 0
            i_res = i_res + 1
            read(line(11:14), '(A)') ml_aa_table(i_res)%res_name
            if (line(11:12) == 'NO') then
                idx_NA = i_res
            end if
            i2 = 0
        else if (line(3:7) == 'LAYER') then
            if (line(9:9) == 'W') then ! weight
                mode = 1
                i_layer = i_layer + 1
                read(line(22:24),'(I3)') n1
                read(line(26:28),'(I3)') n2
                allocate(ml_param(i_layer)%w(n1,n2))
                allocate(ml_param(i_layer)%wt(n2,n1))
                allocate(ml_param(i_layer)%b(n1))                
                i1 = 1
                i2 = 0
            else ! bias
                mode = 2
                i2 = 0
            end if
        end if
    else
        i2 = i2 + 1
        if (mode == 0) then
            read(line, '(F10.7)') ml_aa_table(i_res)%f(i2)
        else if (mode == 1) then
            if (i2 > n2) then
                i2 = 1
                i1 = i1 + 1
            end if
            read(line, '(F10.7)') ml_param(i_layer)%wt(i2,i1)
        else
            read(line, '(F10.7)') ml_param(i_layer)%b(i2)
        end if
    end if
end do

close(f_unit)

do i_layer = 1, n_layer
    ml_param(i_layer)%w = transpose(ml_param(i_layer)%wt)
end do

end subroutine read_ml_param
!-------------------------------------------------------------------------------
subroutine setup_protein_for_ml(protein)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein

integer :: i_ref, i_res, i
character(len=4) :: resName

allocate(amino_param(8, protein%n_res))
allocate(is_glycine(protein%n_res))
allocate(bb_id(5, protein%n_res))
is_glycine = .false.

do i_res = 1, protein%n_res
    i_ref = protein%residue(i_res)%res_type
    resName = ref_res(i_ref)%res_name
    call convert_to_stdres(resName)
    do i = 1, 21
        if (resName == ml_aa_table(i)%res_name) then
            amino_param(:, i_res) = ml_aa_table(i)%f(:)
            exit
        end if
    end do
    !
    bb_id(1, i_res) = res_index(i_res)%bb_id(2,1) ! CA
    bb_id(2, i_res) = res_index(i_res)%bb_id(1,1) ! N
    bb_id(3, i_res) = res_index(i_res)%bb_id(3,1) ! C
    bb_id(4, i_res) = res_index(i_res)%bb_id(4,1) ! O
    if (resName == 'GLY') then
        is_glycine(i_res) = .true.
    else
        bb_id(5, i_res) = res_index(i_res)%Cb_id(1) ! CB
    end if
end do

end subroutine setup_protein_for_ml
!-------------------------------------------------------------------------------
END MODULE ML_SCORE
!-------------------------------------------------------------------------------
