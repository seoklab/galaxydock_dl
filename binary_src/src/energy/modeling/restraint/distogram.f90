!-------------------------------------------------------------------------------
! Copyright (C) 2020, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/restraint/distogram.f90
!
! Description: Reading and calculating (predicted) distogram potential
!
!-------------------------------------------------------------------------------
MODULE DISTOGRAM

use globals
use in_out_utils, only: find_atom_idx
use logger, only: terminate_with_error, log_p
use string, only: parse_longstring
use mathfunctions, only: cubic_spline
!
use energy_vars, only: energy_type, R, ii_R, distogram_file, &
                       use_distogram, distogram_pdb

implicit none
save
private

! Parameters for setting up arrays
integer, parameter :: max_rsrcolumn = 100   ! Maximum # of parameters of a restraint
integer, parameter :: max_rsrlinelen = 500 ! Maximum length of a line in restraint file
integer, parameter :: n_bin = 16
real(dp), parameter :: first_dist = 2.5d0
real(dp), parameter :: dist_interval = 1.0d0
real(dp), parameter :: last_dist = first_dist + (n_bin-1)*dist_interval

integer :: n_dist         ! # of restraints

integer, allocatable :: dist_atm(:,:)
real(dp), allocatable :: dist_table(:,:)
real(dp), allocatable :: distpp_table(:,:)

public :: initialize_distogram
public :: finalize_distogram
public :: calc_distogram_energy

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_distogram(molecule)
!-------------------------------------------------------------------------------
! Initialize distogram function
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule

! For the case using only segment fixing
if (trim(distogram_file) == '') return

!Read distogram variables
call read_distogram_input(distogram_file)

if (n_dist > 0) then
    ! Match restraint index with program index
    call distid_to_topolid(molecule)
    ! Precalculate spline related variable
    call build_spline()
end if

if (n_dist == 0) then
    use_distogram = .false.
    call log_p('- Energy: Distogram energy deactivated, n_dist=0.', me=me, level=10)
end if

end subroutine initialize_distogram
!-------------------------------------------------------------------------------
subroutine finalize_distogram()
!-------------------------------------------------------------------------------
if (allocated(dist_atm)) then
    deallocate(dist_atm)
    deallocate(dist_table, distpp_table)
end if

end subroutine finalize_distogram
!-------------------------------------------------------------------------------
subroutine read_distogram_input(file_name)
!-------------------------------------------------------------------------------
! Read the distogram input file
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: file_name
integer :: f_unit, ioerror, num_word, openstat
integer :: i_rsr, i
character(len=len_fname) :: word(max_rsrcolumn)
character(len=max_rsrlinelen) :: line

!1. First get # dist for allocation
f_unit = 38
open(f_unit, file = trim(file_name), iostat=openstat)
if (openstat > 0) then
    call terminate_with_error('Terminate with error: No distogram file found. Check {distogram_file}.')
end if
  
n_dist = 0
do
    read(f_unit,"(A500)", iostat = ioerror) line
    if (ioerror < 0) exit
    n_dist = n_dist + 1
end do

if (n_dist == 0) then
    call terminate_with_error('No distogram rsr was found from the file. Please check.')
else
    write(log_msg,"(I10,A)") n_dist, ' Restraints found from {distogram_file}.'
    call log_p(log_msg, me=me, level=10)
end if

!2. Allocate array sizes with given input file
allocate(dist_atm(2, n_dist))
allocate(dist_table(n_bin, n_dist))
allocate(distpp_table(n_bin, n_dist))
close(f_unit)

!3. Then read input file again and assign parameters
open(f_unit, file = trim(file_name))

i_rsr = 0
do
    read(f_unit,"(A500)", iostat = ioerror) line
    if (ioerror < 0) exit

    call parse_longstring(line, num_word, word, max_rsrlinelen)
    if (num_word < 2) exit

    i_rsr = i_rsr + 1
    ! Atom index
    do i = 1, 2
        read(word(8+i),"(I6)") dist_atm(i, i_rsr)
    end do
    ! Restraint parameter
    do i = 1, n_bin
        read(word(i+16),"(F10.4)") dist_table(i, i_rsr)
    end do
end do

close(f_unit)

end subroutine read_distogram_input
!-------------------------------------------------------------------------------
subroutine distid_to_topolid(molecule)
!-------------------------------------------------------------------------------
! Map atom index of rsrs into atom index of the program
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule
integer :: max_rsr_atm, max_rsr_atm_prev
integer :: f_unit, ioerror, openstat
integer :: i_atm, i_rsr, res_no, atm_no, ref_res_no, resatm_no
integer :: res_no_pdb, res_no_prv
character(len=6) :: atom_error_mode
character(len=4) :: atmtype
character(len=max_rsrlinelen) :: line
! Map of (atmno,resno) <-> rsr atom index
integer, allocatable :: map_atomid(:,:), map(:,:)

max_rsr_atm = tn%atom
allocate(map_atomid(2,max_rsr_atm))

! Parse pdb file
i_atm = 0
res_no_prv = -9999
res_no = 0
f_unit = 25
atom_error_mode = 'ignore'

open(f_unit, file = trim(distogram_pdb), iostat = openstat)
if (openstat > 0) then
    call terminate_with_error('Terminate with error: No distogram_pdb file found. Check {distogram_pdb}.')
end if

do
    read(f_unit,"(A500)", iostat = ioerror) line
    if (ioerror < 0) exit
    if (line(1:4) /= 'ATOM') then
        cycle
    end if
     
    i_atm = i_atm + 1

    if (line(13:13) == ' ') then
        read(line(14:16),"(A4)") atmtype
    else !If atmtype is written like "1HH1"
        read(line(13:16),"(A4)") atmtype
        atmtype = trim(atmtype(2:4))//atmtype(1:1)
    end if
    read(line(23:26),"(I4)") res_no_pdb
    if (res_no_pdb /= res_no_prv) then
        res_no = res_no + 1
    end if
    res_no_prv = res_no_pdb

    if (i_atm > max_rsr_atm) then
        ! reallocate map_atomid by 2*(max_rsr_atm+tn%atom)
        max_rsr_atm_prev = max_rsr_atm
        max_rsr_atm = max_rsr_atm + tn%atom
        allocate(map(2,max_rsr_atm))
        map(:,1:max_rsr_atm_prev) = map_atomid(:,1:max_rsr_atm_prev)
        !
        call move_alloc(map, map_atomid)
    end if

    ref_res_no = molecule%residue(res_no)%res_type
    call find_atom_idx(res_no, ref_res_no, atmtype, resatm_no, atom_error_mode)
    map_atomid(1,i_atm) = res_no
    map_atomid(2,i_atm) = resatm_no
end do

close(f_unit)

!Reset atom indices
do i_rsr = 1, n_dist
    do i_atm = 1, 2
        res_no = map_atomid(1, dist_atm(i_atm, i_rsr))
        atm_no = map_atomid(2, dist_atm(i_atm, i_rsr))
        dist_atm(i_atm, i_rsr) = ii_R(atm_no,res_no)
    end do
end do
    
deallocate(map_atomid)

end subroutine distid_to_topolid
!-------------------------------------------------------------------------------
subroutine build_spline()
!-------------------------------------------------------------------------------
! Prepare spline array spl_ypp (second derivative)
!-------------------------------------------------------------------------------
integer :: i, k, i_rsr, n
real(dp) :: c(n_bin), c_pp(n_bin)
real(dp) :: p, u(n_bin), c_p0, c_pn
real(dp), parameter :: sig = 0.5d0

do i_rsr = 1, n_dist
    n = n_bin

    ! Build second derivative recursively
    c = dist_table(:,i_rsr)
    c_p0 = 0.0d0
    c_pn = 0.0d0

    c_pp(:) = 0.0d0
    ! with a specified first derivative c_p0,
    c_pp(1) = -0.5d0
    u(1) = 3.0d0/dist_interval * ((c(2)-c(1))/dist_interval-c_p0)

    do i = 2, n-1
        p = sig*c_pp(i-1) + 2.0d0
        c_pp(i) = (sig - 1.0d0)/p
        u(i) = (6.0d0*((c(i+1) - c(i))/dist_interval - (c(i)-c(i-1))/dist_interval) / &
             (2.0d0*dist_interval) - sig*u(i-1)) / p
    end do

    ! with a specified first derivative c_pn,
    c_pp(n) = 0.5d0
    u(n) = 3.0d0/dist_interval * (c_pn-(c(n)-c(n-1))/dist_interval)
    do k = n-1, 1, -1
        c_pp(k) = c_pp(k)*c_pp(k+1) + u(k)
    end do

    ! Save obtained spl_ypp values
    distpp_table(:,i_rsr) = c_pp
end do

end subroutine build_spline
!-------------------------------------------------------------------------------
subroutine calc_distogram_energy(ff, g, calc_g)
!-------------------------------------------------------------------------------
! Main subroutine for distogram potential
!-------------------------------------------------------------------------------
real(dp), intent(out) :: ff
real(dp), intent(out) :: g(3,tn%atom)
logical, intent(in) :: calc_g
integer :: i_rsr, i_atm, j_atm, i_bin
real(dp) :: dist, dr(3), dcdf
real(dp) :: e_val
! Spline variable
real(dp) :: argv(2)
real(dp) :: spl_a

ff = 0.0d0
g = 0.0d0

do i_rsr = 1, n_dist
    i_atm = dist_atm(1, i_rsr)
    j_atm = dist_atm(2, i_rsr)
    dr = R(:,j_atm) - R(:,i_atm)
    dist = sqrt(dot_product(dr, dr))

    ! Check the location of feature in spline table
    i_bin = int((dist - first_dist)/dist_interval) + 1
    if (i_bin < 1) then
        e_val = dist_table(1, i_rsr)
    else if (i_bin >= n_bin) then
        e_val = dist_table(n_bin, i_rsr)
    else
        spl_a = (first_dist - dist)/dist_interval + i_bin
        argv(:) = cubic_spline(dist_interval, dist_table(i_bin:i_bin+1, i_rsr), &
                               distpp_table(i_bin:i_bin+1, i_rsr), spl_a, calc_g)
        e_val = argv(1)
        !
        if (calc_g) then
            dcdf = argv(2)
            dr = dr / dist
            g(:,i_atm) = g(:,i_atm) - dcdf*dr
            g(:,j_atm) = g(:,j_atm) + dcdf*dr
        end if
    end if
    ff = ff + e_val
end do

end subroutine calc_distogram_energy
!-------------------------------------------------------------------------------
END MODULE DISTOGRAM
