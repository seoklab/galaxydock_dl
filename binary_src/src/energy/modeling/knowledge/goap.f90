!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/knowledge/goap.f90
!
! Description: GOAP statistical potential
!
! Reference : Hongyi Zhou and Jeffrey Skolnick, GOAP: A Generalized
! Orientation-Dependent, All-Atom Statistical Potential for Protein
! Structure Prediction, Biophysical Journal, 2011, 101, 2043-2052.
!
!-------------------------------------------------------------------------------
MODULE GOAP
!-------------------------------------------------------------------------------

use globals
use logger,           only: log_p
use string,           only: parse_string
use geometry,         only: calc_tor_and_grad
use in_out_utils,     only: find_atom_idx
use mathfunctions,    only: cubic_spline, cross, v_norm
use convert_res_name, only: convert_to_stdres
!
use energy_vars,      only: R, i_R, ii_R, i_P, max_neigh2

implicit none
save
private

!===============================================================================
! VARIABLES for GOAP score
!===============================================================================
logical :: use_goap_discrete    ! To calculate GOAP with discrete version or not
                                ! The original GOAP uses discrete version,
                                !  but it can be continuous version with splines.
                                ! However, it is not implemented due to speed issue

!===============================================================================
! PARAMETERS for GOAP score
!===============================================================================
integer, parameter :: tot_n_goap_atom_type = 167   ! Total number of GOAP atom types
integer, parameter :: tot_n_rbin = 20              ! Total number of bins for radial components
integer, parameter :: tot_n_abin = 12              ! Total number of bins for angular components
integer, parameter :: tot_n_atype = 5              ! Total number of angular terms
integer, parameter :: sequence_separation = 7
integer, parameter :: rmap(60) = &                 ! mapping r -> ir (rbin index)
                        (/ 1, 1, 1, 1, 1,  1, 1, 1, 2, 2, &
                           3, 3, 4, 4, 5,  5, 6, 6, 7, 7, &
                           8, 8, 9, 9,10, 10,11,11,12,12, &
                          13,13,13,14,14, 14,14,15,15,15, &
                          15,16,16,16,16, 17,17,17,17,18, &
                          18,18,18,19,19, 19,19,20,20,20 /)
integer, parameter :: rmap_discrete(30) = &        ! mapping r -> ir
                        (/ 1, 1, 1, 1, 2,  3, 4, 5, 6, 7, &
                           8, 9,10,11,12, 13,14,14,15,15, &
                          16,16,17,17,18, 18,19,19,20,20 /)
real(dp), parameter :: r_radial(tot_n_rbin) = &    ! x for radial score binning
                        (/  0.00d0, 2.25d0, 2.75d0, 3.25d0, 3.75d0, &
                            4.25d0, 4.75d0, 5.25d0, 5.75d0, 6.25d0, &
                            6.75d0, 7.25d0, 7.75d0, 8.50d0, 9.50d0, &
                           10.50d0,11.50d0,12.50d0,13.50d0,14.50d0 /)
real(dp), parameter :: dx_radial(tot_n_rbin-1) = & ! dx for radial score binning
                        (/  2.25d0, 0.50d0, 0.50d0, 0.50d0, 0.50d0, &
                            0.50d0, 0.50d0, 0.50d0, 0.50d0, 0.50d0, &
                            0.50d0, 0.50d0, 0.75d0, 1.00d0, 1.00d0, &
                            1.00d0, 1.00d0, 1.00d0, 1.00d0         /)
real(dp), parameter :: dx_angular = 30.0d0         ! dx for angular score binning

!===============================================================================
! DERIVED TYPES for GOAP score
!===============================================================================
type goap_atom_type
!-------------------------------------------------------------------------------
! A type to save GOAP atom type information
!-------------------------------------------------------------------------------
character(len=6) :: atom
integer :: i2, i3
logical :: two_bonds
!-------------------------------------------------------------------------------
end type goap_atom_type
!-------------------------------------------------------------------------------
type radial_score_type
!-------------------------------------------------------------------------------
! A type to save GOAP radial score data
!-------------------------------------------------------------------------------
real :: f
real :: fpp                   ! d2f/dx2
!-------------------------------------------------------------------------------
end type radial_score_type
!-------------------------------------------------------------------------------
type angular_score_type
!-------------------------------------------------------------------------------
! A type to save GOAP angular score data
!-------------------------------------------------------------------------------
real :: f
real :: fpp(3)                ! df/dx, df/dy, d2f/dxdy
!-------------------------------------------------------------------------------
end type angular_score_type
!-------------------------------------------------------------------------------
type goap_index_type
!-------------------------------------------------------------------------------
! A type to save Protein-GOAP information for evaluations
!-------------------------------------------------------------------------------
logical :: disabled
integer :: ig
integer :: i2, i3
logical :: two_bonds
!-------------------------------------------------------------------------------
end type goap_index_type
!-------------------------------------------------------------------------------

type(radial_score_type),  allocatable :: radial_table(:,:,:)       ! i_dist, ib, ia
type(angular_score_type), allocatable :: angular_table(:,:,:,:,:)  ! i_ang, i_dist, ang_type, ib, ia

type(goap_index_type), allocatable :: i_G(:)

!-------------------------------------------------------------------------------

public :: initialize_goap
public :: finalize_goap
!
public :: calc_goap_energy

CONTAINS
!===============================================================================
subroutine initialize_goap(protein)
!-------------------------------------------------------------------------------
! Initialization Processes
!  1. Reading GOAP data files
!   - GOAP atom definitions
!   - GOAP score files (radial, angular)
!  2. Setting up proteins
!   - Assignment of GOAP atom indices
!   - Setting up GOAP geometry indices
!  3. Pre-calculations of splines
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein

character(len=len_fname) :: infile_goap_atom_list
character(len=len_fname) :: infile_goap_radial_score
character(len=len_fname) :: infile_goap_angular_score
type(goap_atom_type) :: goap_atom_list(tot_n_goap_atom_type)

use_goap_discrete = .true.

infile_goap_atom_list     = trim(data_dir) // 'goap_atomtype.list'
infile_goap_radial_score  = trim(data_dir) // 'goap_radial.score'
infile_goap_angular_score = trim(data_dir) // 'goap_angular.score'

allocate(i_G(tn%stdatm))
allocate(radial_table(tot_n_rbin, tot_n_goap_atom_type, tot_n_goap_atom_type))
allocate(angular_table(0:tot_n_abin, tot_n_rbin, 5, tot_n_goap_atom_type, tot_n_goap_atom_type))

write(log_msg, "(A)") "  Reading GOAP Atom List"
call log_p(log_msg, me=me, level=40)
call read_goap_atom_list(infile_goap_atom_list, goap_atom_list)

write(log_msg, "(A)") "  Reading GOAP Radial Table"
call log_p(log_msg, me=me, level=40)
call read_radial_table(infile_goap_radial_score)

write(log_msg, "(A)") "  Reading GOAP Angular Table"
call log_p(log_msg, me=me, level=40)
call read_angular_table(infile_goap_angular_score)

write(log_msg, "(A)") "  Setting up the Protein for GOAP Evaluation"
call log_p(log_msg, me=me, level=40)
call setup_protein_for_goap(protein, goap_atom_list)

if (.not. use_goap_discrete) then
    write(log_msg, "(A)") "  Interpolation of GOAP data files"
    call log_p(log_msg, me=me, level=40)
    call setup_radial_spline()
    call setup_angular_spline()
end if

end subroutine initialize_goap
!-------------------------------------------------------------------------------
subroutine finalize_goap()
!-------------------------------------------------------------------------------
! Finalization Processes
!-------------------------------------------------------------------------------
deallocate(i_G)
deallocate(radial_table)
deallocate(angular_table)

end subroutine finalize_goap
!-------------------------------------------------------------------------------
subroutine calc_goap_energy(ff, gg, appl_respair, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: ff(2), gg(3,tn%stdatm,2)
logical, intent(in) :: appl_respair(tn%residue,tn%residue), calc_g

integer :: i_atm, j_atm, ia, ig, jg, ir
integer :: n_pair, i_pair(max_neigh2)
logical :: use_ang(max_neigh2)
real(dp) :: Ri(3), Rj(3), Rij(3), ar
real :: dij
real(dp) :: fr(2), fr_pp(2), ff_r(2), dff_r(3)
real(dp) :: Xn(3, tn%stdatm), Xp(3, tn%stdatm)

ff(:) = 0.0d0
gg(:,:,:) = 0.0d0
 
if (use_goap_discrete) then
    call calc_goap_energy_discrete(ff, appl_respair)
    return
end if

call update_dplanes(Xn, Xp)

do i_atm = 1, tn%stdatm
    if (i_G(i_atm)%disabled) cycle
    ig = i_G(i_atm)%ig
    Ri(1:3) = R(1:3,i_atm)
    call get_goap_pair(i_atm, appl_respair, n_pair, i_pair, use_ang)

    do ia = 1, n_pair
        j_atm = i_pair(ia)
        jg = i_G(j_atm)%ig
        Rj(1:3) = R(1:3, j_atm)
        Rij = Rj - Ri
        dij = sqrt(dot_product(Rij,Rij))

        ! Calculating GOAP radial term
        ir = rmap(int(dij*4.0d0))
        ar = 1.0-(dij-r_radial(ir))/dx_radial(ir)
        fr(1:2)    = radial_table(ir:ir+1, jg, ig)%f
        fr_pp(1:2) = radial_table(ir:ir+1, jg, ig)%fpp
        ff_r = cubic_spline(dx_radial(ir), fr, fr_pp, ar, calc_g)

        ff(1) = ff(1) + ff_r(1)
        if (calc_g) then
            dff_r(1:3) = ff(2)*Rij/dij
            gg(:,i_atm,1) = gg(:,i_atm,1) - dff_r
            gg(:,j_atm,1) = gg(:,j_atm,1) + dff_r
        end if
    end do
end do

end subroutine calc_goap_energy
!-------------------------------------------------------------------------------
subroutine calc_goap_energy_discrete(ff, appl_respair)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: ff(2)
logical, intent(in) :: appl_respair(tn%residue, tn%residue)

integer :: i_atm, j_atm, ia, ig, jg, ir
integer :: n_pair, i_pair(max_neigh2)
logical :: use_ang(max_neigh2)
real(dp) :: Ri(3), Rj(3), Rij(3)
real :: dij, ff_r
integer :: it(5)
real(dp) :: Xn(3, tn%stdatm), Xd(3, tn%stdatm)
real(dp) :: theta(5), ff_t(5)
real(dp) :: Rd(3,4), dRddr(3,4)

ff(:) = 0.0d0

call update_dplanes(Xn, Xd)

do i_atm = 1, tn%stdatm
    if (i_G(i_atm)%disabled) cycle
    ig = i_G(i_atm)%ig
    Ri(1:3) = R(1:3,i_atm)
    call get_goap_pair(i_atm, appl_respair, n_pair, i_pair, use_ang)
    if (i_R(1,i_atm) == 1) cycle

    do ia = 1, n_pair
        j_atm = i_pair(ia)
        if (i_R(1,j_atm) == 1) cycle
        jg = i_G(j_atm)%ig
        Rj(1:3) = R(1:3, j_atm)
        Rij = Rj - Ri
        dij = sqrt(dot_product(Rij,Rij))
        Rij = Rij/dij

        ! Calculating GOAP radial term
        ir = rmap_discrete(int(dij*2.0d0)+1)
        ff_r = radial_table(ir, jg,ig)%f
        ff(1) = ff(1) + ff_r

        ! Calculating GOAP angular term
        if (.not. use_ang(ia)) cycle
        call calc_goap_ang(Xn(:,i_atm), Rij, theta(1))
        call calc_goap_phi(Xn(:,i_atm), Xd(:,i_atm), Rij, theta(2))
        call calc_goap_ang(Xn(:,j_atm),-Rij, theta(3))
        call calc_goap_phi(Xn(:,j_atm), Xd(:,j_atm),-Rij, theta(4))
        !
        Rd(1:3,1) = Ri + Xn(:,j_atm)
        Rd(1:3,2) = Ri
        Rd(1:3,3) = Rj
        Rd(1:3,4) = Rj + Xn(:,i_atm)
        call calc_tor_and_grad(Rd, theta(5), dRddr, .false.)
        theta(5) = -theta(5)
        !
        it(1) = max(min(int((theta(1)+1.001)*6.0)+1, 12), 1)
        it(2) = max(min(int((theta(2)*rad2deg+180.001)/dx_angular)+1, 12), 1)
        it(3) = max(min(int((theta(3)+1.001)*6.0)+1, 12), 1)
        it(4) = max(min(int((theta(4)*rad2deg+180.001)/dx_angular)+1, 12), 1)
        it(5) = max(min(int((theta(5)*rad2deg+180.001)/dx_angular)+1, 12), 1)
        !
        ff_t(1) = angular_table(it(1), ir, 1, jg, ig)%f
        ff_t(2) = angular_table(it(2), ir, 2, jg, ig)%f
        ff_t(3) = angular_table(it(3), ir, 3, jg, ig)%f
        ff_t(4) = angular_table(it(4), ir, 4, jg, ig)%f
        ff_t(5) = angular_table(it(5), ir, 5, jg, ig)%f
        !
        ff(2) = ff(2) + sum(ff_t)
    end do
end do

end subroutine calc_goap_energy_discrete
!-------------------------------------------------------------------------------
subroutine calc_pairwise_goap_energy(ff, appl_respair)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: ff(tn%residue,tn%residue)
logical, intent(in) :: appl_respair(tn%residue,tn%residue)

ff(:,:) = 0.0d0
 
if (use_goap_discrete) then
    call calc_pairwise_goap_energy_discrete(ff, appl_respair)
    return
end if

end subroutine calc_pairwise_goap_energy
!-------------------------------------------------------------------------------
subroutine calc_pairwise_goap_energy_discrete(ff, appl_respair)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: ff(tn%residue,tn%residue)
logical, intent(in) :: appl_respair(tn%residue,tn%residue)

integer :: i_res, j_res
integer :: i_atm, j_atm, ia, ig, jg, ir
integer :: n_pair, i_pair(max_neigh2)
logical :: use_ang(max_neigh2)
real(dp) :: Ri(3), Rj(3), Rij(3)
real :: dij, ff_r
integer :: it(5)
real(dp) :: Xn(3, tn%stdatm), Xd(3, tn%stdatm)
real(dp) :: theta(5), ff_t(5)
real(dp) :: Rd(3,4), dRddr(3,4)

ff(:,:) = 0.0d0

call update_dplanes(Xn, Xd)

do i_atm = 1, tn%stdatm
    if (i_G(i_atm)%disabled) cycle
    ig = i_G(i_atm)%ig
    i_res = i_R(1,i_atm)
    !if (i_res == 1) cycle  ! Modified from the original GOAP

    Ri(1:3) = R(1:3,i_atm)
    call get_goap_pair(i_atm, appl_respair, n_pair, i_pair, use_ang)

    do ia = 1, n_pair
        j_atm = i_pair(ia)
        j_res = i_R(1,j_atm)
        !if (j_res == 1) cycle  ! Modified from the origianl GOAP
        jg = i_G(j_atm)%ig
        Rj(1:3) = R(1:3, j_atm)
        Rij = Rj - Ri
        dij = sqrt(dot_product(Rij,Rij))
        Rij = Rij/dij

        ! Calculating GOAP radial term
        ir = rmap_discrete(int(dij*2.0d0)+1)
        ff_r = radial_table(ir, jg,ig)%f
        ff(i_res, j_res) = ff(i_res, j_res) + ff_r
        ff(j_res, i_res) = ff(j_res, i_res) + ff_r

        ! Calculating GOAP angular term
        if (.not. use_ang(ia)) cycle
        call calc_goap_ang(Xn(:,i_atm), Rij, theta(1))
        call calc_goap_phi(Xn(:,i_atm), Xd(:,i_atm), Rij, theta(2))
        call calc_goap_ang(Xn(:,j_atm),-Rij, theta(3))
        call calc_goap_phi(Xn(:,j_atm), Xd(:,j_atm),-Rij, theta(4))
        !
        Rd(1:3,1) = Ri + Xn(:,j_atm)
        Rd(1:3,2) = Ri
        Rd(1:3,3) = Rj
        Rd(1:3,4) = Rj + Xn(:,i_atm)
        call calc_tor_and_grad(Rd, theta(5), dRddr, .false.)
        theta(5) = -theta(5)
        !
        it(1) = max(min(int((theta(1)+1.001)*6.0)+1, 12), 1)
        it(2) = max(min(int((theta(2)*rad2deg+180.001)/dx_angular)+1, 12), 1)
        it(3) = max(min(int((theta(3)+1.001)*6.0)+1, 12), 1)
        it(4) = max(min(int((theta(4)*rad2deg+180.001)/dx_angular)+1, 12), 1)
        it(5) = max(min(int((theta(5)*rad2deg+180.001)/dx_angular)+1, 12), 1)
        !
        ff_t(1) = angular_table(it(1), ir, 1, jg, ig)%f
        ff_t(2) = angular_table(it(2), ir, 2, jg, ig)%f
        ff_t(3) = angular_table(it(3), ir, 3, jg, ig)%f
        ff_t(4) = angular_table(it(4), ir, 4, jg, ig)%f
        ff_t(5) = angular_table(it(5), ir, 5, jg, ig)%f
        !
        ff(i_res, j_res) = ff(i_res, j_res) + sum(ff_t)
        ff(j_res, i_res) = ff(j_res, i_res) + sum(ff_t)
        !
    end do
end do

end subroutine calc_pairwise_goap_energy_discrete
!-------------------------------------------------------------------------------
subroutine get_goap_pair(i_atm, appl_respair, n_pair, i_pair, use_ang)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_atm
logical, intent(in) :: appl_respair(tn%residue,tn%residue)
integer, intent(out) :: n_pair
integer, intent(out) :: i_pair(max_neigh2)
logical, intent(out) :: use_ang(max_neigh2)

integer :: ia_start
integer :: i_res, j_res, j_atm, ia

n_pair = 0
i_res = i_R(1,i_atm)

! The original GOAP score definition
!if (use_goap_discrete) then
!    ia_start = 1
!else
!    ia_start = i_P(i_atm)%pair_end_index(2)+1
!end if
ia_start = i_P(i_atm)%pair_end_index(2)+1

do ia = ia_start, i_P(i_atm)%n_pair
    j_atm = i_P(i_atm)%i_pair(ia)
    j_res = i_R(1,j_atm)
    if (j_res /= i_res .and. (.not. i_G(j_atm)%disabled) .and. appl_respair(i_res,j_res)) then
        n_pair = n_pair + 1
        i_pair(n_pair) = j_atm
        if (abs(j_res-i_res) >= sequence_separation) then
            use_ang(n_pair) = .true.
        else
            use_ang(n_pair) = .false.
        end if
        if (n_pair == max_neigh2) return
    end if
end do

do ia = 1, i_P(i_atm)%n_Lpair
    j_atm = i_P(i_atm)%i_Lpair(ia)
    j_res = i_R(1,j_atm)
    if (j_res /= i_res .and. (.not. i_G(j_atm)%disabled) .and. appl_respair(i_res,j_res)) then
        n_pair = n_pair + 1
        i_pair(n_pair) = j_atm
        if (abs(j_res-i_res) >= sequence_separation) then
            use_ang(n_pair) = .true.
        else
            use_ang(n_pair) = .false.
        end if
        if (n_pair == max_neigh2) return
    end if
end do

end subroutine get_goap_pair
!-------------------------------------------------------------------------------
subroutine read_goap_atom_list(infile_goap_atom_list, goap_atom_list)
!-------------------------------------------------------------------------------
! Reading GOAP atom list 
!  In the file,
!   RES ATM BND ANG TWO_BONDS
!   RES: Residue Name
!   ATM: Atom Name
!   BND: Bond defining GOAP atom index
!   ANG: Angle defining GOAP atom index
!   TWO_BONDS: Whether it is a tip atom or not
!   ie, 
!    CYS N    -1   2 T
!    CYS CA    1   3 T
!    CYS C     2   4 T
!    CYS O     3   2 F
!    ...
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: infile_goap_atom_list
type(goap_atom_type), intent(inout) :: goap_atom_list(tot_n_goap_atom_type)

character(len=len_fname) :: line, word(10)
character(len=6) :: atom
integer :: f_unit, io, n_word, ia

f_unit = 41
open(f_unit, file=trim(infile_goap_atom_list))

ia = 0
do
    read(f_unit, '(A120)', iostat=io) line
    if (io < 0) exit

    call parse_string(line, n_word, word)

    write(atom, "(A3,A3)") word(1)(1:3), word(2)(1:3)
    ia = ia + 1
    goap_atom_list(ia)%atom = atom
    read(word(3), "(I3)") goap_atom_list(ia)%i2
    read(word(4), "(I3)") goap_atom_list(ia)%i3
    if (trim(word(5)) == 'T') then
        goap_atom_list(ia)%two_bonds = .true.
    else
        goap_atom_list(ia)%two_bonds = .false.
    end if
end do

close(f_unit)

end subroutine read_goap_atom_list
!-------------------------------------------------------------------------------
subroutine read_radial_table(infile_goap_radial_score)
!-------------------------------------------------------------------------------
! Reading GOAP radial score data file
!  The file contains f(r=i_r*0.25; iG_1,iG_2) = ff
!  In the file,
!   RES_1 ATM_1 iG_1 RES_2 ATM_2 iG_2 i_r  ff
!  ie,
!   CYS N     1  CYS N     1   1    1.55106779
!   CYS N     1  CYS N     1   2    1.55661016
!   CYS N     1  CYS N     1   3    1.69516949
!   ...
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: infile_goap_radial_score

character(len=len_fname) :: line
integer :: f_unit, io
integer :: ia, ib, ir
real(dp) :: ff

f_unit = 42
open(f_unit, file=trim(infile_goap_radial_score))

do
    read(f_unit, "(A120)", iostat=io) line
    if (io < 0) exit

    read(line( 9:11), "(I3)") ia
    read(line(22:24), "(I3)") ib
    read(line(27:28), "(I2)") ir
    read(line(31:42), "(f12.8)") ff
    radial_table(ir,ib,ia)%f = ff
end do

do ia = 1, tot_n_goap_atom_type-1
    do ib = ia+1, tot_n_goap_atom_type
        radial_table(:,ia,ib) = radial_table(:,ib,ia)
    end do
end do

close(f_unit)

end subroutine read_radial_table
!-------------------------------------------------------------------------------
subroutine read_angular_table(infile_goap_angular_score)
!-------------------------------------------------------------------------------
! Reading GOAP angular score data file
!  The file contains f(r=i_r*0.25, j_a; iG_1,iG_2,iT) = ff
!   ,where iT is torsion angle type.
!  In the file,
!   RES_1 ATM_1 iG_1 RES_2 ATM_2 iG_2 i_r iT ff(j_a=1..12)
!  ie,   
!   CYS N     1  CYS N     1   1  1   0.011  0.014 -0.013 -0.010  0.011  0.011  0.004  0.010  0.010  0.004  0.009  0.003
!   CYS N     1  CYS N     1   2  1   0.011  0.014 -0.013 -0.010  0.010  0.011  0.004  0.010  0.010  0.004  0.009  0.003
!   CYS N     1  CYS N     1   3  1   0.008  0.009 -0.009 -0.007  0.007  0.007  0.003  0.007  0.007  0.003  0.006  0.002
!   ...
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: infile_goap_angular_score

character(len=len_fname) :: line
integer :: f_unit, io
integer :: ia, ib, ir, it, i
real(dp) :: ff(tot_n_abin)

f_unit = 43
open(f_unit, file=trim(infile_goap_angular_score))

do
    read(f_unit, "(A120)", iostat=io) line
    if (io < 0) exit

    read(line( 9:11), "(I3)") ia
    read(line(22:24), "(I3)") ib
    read(line(27:28), "(I2)") ir
    read(line(31:31), "(I1)") it
    read(line(33:116),"(12f7.3)") ff
    angular_table(0,ir,it,ib,ia)%f = ff(tot_n_abin)
    do i = 1, tot_n_abin
        angular_table(i,ir,it,ib,ia)%f = ff(i)
    end do
end do

do ia = 1, tot_n_goap_atom_type-1
    do ib = ia+1, tot_n_goap_atom_type
        angular_table(:,:,1,ia,ib) = angular_table(:,:,3,ib,ia)
        angular_table(:,:,2,ia,ib) = angular_table(:,:,4,ib,ia)
        angular_table(:,:,3,ia,ib) = angular_table(:,:,1,ib,ia)
        angular_table(:,:,4,ia,ib) = angular_table(:,:,2,ib,ia)
        angular_table(:,:,5,ia,ib) = angular_table(:,:,5,ib,ia)
    end do
end do

close(f_unit)

end subroutine read_angular_table
!-------------------------------------------------------------------------------
subroutine find_goap_atom_index(goap_atom_list, atomName, i_atm)
!-------------------------------------------------------------------------------
type(goap_atom_type), intent(in) :: goap_atom_list(tot_n_goap_atom_type)
character(len=6), intent(in) :: atomName
integer, intent(out) :: i_atm
integer :: ia

i_atm = -1
do ia = 1, tot_n_goap_atom_type
    if (atomName == goap_atom_list(ia)%atom) then
        i_atm = ia
        return
    end if
end do

end subroutine find_goap_atom_index
!-------------------------------------------------------------------------------
subroutine setup_protein_for_goap(protein, goap_atom_list)
!-------------------------------------------------------------------------------
! Initialization of the protein to evaluate GOAP
!  1. Assigning GOAP atom type indices for every protein atoms
!  2. Finding connected atoms defining dihedral angles
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(goap_atom_type), intent(in) :: goap_atom_list(tot_n_goap_atom_type)

integer :: i_res, i_ref, i_atm, j_atm, ig
integer :: i_res_prev, i_ref_prev
integer :: ia, ia_s(max_atm)
integer :: i2, i3
character(len=4) :: resName, atmName
character(len=6) :: goapName
character(len=6), parameter :: error_mode = "ignore"

! Initialization
do i_atm = 1, tn%stdatm
    i_G(i_atm)%disabled = .false.
    i_G(i_atm)%ig = 0
    i_G(i_atm)%i2 = 0
    i_G(i_atm)%i3 = 0
    i_G(i_atm)%two_bonds = .false.
end do

do i_res = 1, protein%n_res
    i_ref = protein%residue(i_res)%res_type
    resName = ref_res(i_ref)%res_name
    call convert_to_stdres(resName)
    !
    ia_s(:) = -1    ! To save atom numbers for a residue, i_res
    ! Assign GOAP atom indices 
    do i_atm = 1, ref_res(i_ref)%n_atm
        ia = ii_R(i_atm, i_res)
        ia_s(i_atm) = ia
        atmName = ref_res(i_ref)%atom_name(i_atm)
        if (atmName(1:3) == 'OXT') then
            atmName = 'O'
        else if (atmName(1:2) == 'OT') then
            atmName = 'O'
        else if (resName(1:3) == 'ILE' .and. atmName(1:3) == 'CD1') then
            atmName = 'CD'
        end if

        write(goapName, "(A3,A3)") resName(1:3), atmName(1:3)
        call find_goap_atom_index(goap_atom_list, goapName, ig)
        if (ig == -1) then
            i_G(ia)%disabled = .true.
            cycle
        end if

        i_G(ia)%ig = ig
        i_G(ia)%two_bonds = goap_atom_list(ig)%two_bonds
    end do

    ! Assign bond/angle defining atom indices
    do i_atm = 1, ref_res(i_ref)%n_atm
        ia = ia_s(i_atm)
        ig = i_G(ia)%ig
        if (ig <= 0) cycle
        !
        i2 = -1
        if (goap_atom_list(ig)%i2 /= -1) then
            do j_atm = 1, ref_res(i_ref)%n_atm
                if (i_G(ia_s(j_atm))%ig == goap_atom_list(ig)%i2) then
                        i2 = ia_s(j_atm)
                        exit
                end if
            end do
        else if (protein%residue(i_res)%ter_type /= 'N') then
            i_res_prev = i_res-1
            i_ref_prev = protein%residue(i_res_prev)%res_type

            call find_atom_idx(i_res_prev, i_ref_prev, "C   ", j_atm, error_mode)
            i2 = ii_R(j_atm,i_res_prev)
        end if
        !
        i3 = -1
        do j_atm = 1, ref_res(i_ref)%n_atm
            if (i_G(ia_s(j_atm))%ig == goap_atom_list(ig)%i3) then
                i3 = ia_s(j_atm)
            end if
        end do
        !
        i_G(ia)%i2 = i2
        i_G(ia)%i3 = i3
        if (i2 == -1 .or. i3 == -1) then
            i_G(ia)%disabled = .true.
        end if
    end do
end do

end subroutine setup_protein_for_goap
!-------------------------------------------------------------------------------
subroutine setup_radial_spline()
!-------------------------------------------------------------------------------
integer :: ig, jg
real(dp) :: f(tot_n_rbin), fpp(tot_n_rbin)

! Evaluation TEST
logical :: calc_g = .true.

do ig = 1, tot_n_goap_atom_type
    do jg = 1, tot_n_goap_atom_type
        f(1:tot_n_rbin) = radial_table(:,jg,ig)%f
        call get_spline_derivatives(tot_n_rbin, f, r_radial, fpp)
        radial_table(1:tot_n_rbin,jg,ig)%fpp = fpp(1:tot_n_rbin)
    end do
end do

end subroutine setup_radial_spline
!-------------------------------------------------------------------------------
subroutine setup_angular_spline()
!-------------------------------------------------------------------------------
integer :: ig, jg, it, ir, ia
real(dp) :: f(0:tot_n_abin, tot_n_rbin), fxy(0:tot_n_rbin, tot_n_rbin)
real(dp) :: fx(0:tot_n_abin, tot_n_rbin), fy(0:tot_n_abin, tot_n_rbin)

do ig = 1, tot_n_goap_atom_type
    do jg = 1, tot_n_goap_atom_type
        do it = 1, tot_n_atype
            f(0:tot_n_abin,1:tot_n_rbin) = &
                        angular_table(0:tot_n_abin,1:tot_n_rbin,it,jg,ig)%f
            !
            do ir = 1, tot_n_rbin
                fx(0,ir)              = f(1,ir) - f(tot_n_abin-1,ir)
                fx(1:tot_n_abin-1,ir) = f(2:tot_n_abin,ir) - f(0:tot_n_abin-2,ir)
                fx(tot_n_abin,ir)     = f(0,ir)
                fx(0:tot_n_abin,ir)   = fx(0:tot_n_abin,ir) / (2.0d0*dx_angular)
            end do
            !
            do ia = 0, tot_n_abin
                call get_spline_derivatives(tot_n_rbin,f(ia,1:tot_n_rbin), &
                                                        r_radial, fy(ia,1:tot_n_rbin))
                call get_spline_derivatives(tot_n_rbin,fx(ia,1:tot_n_rbin), &
                                                        r_radial, fxy(ia,1:tot_n_rbin))
            end do
            !
            angular_table(:,:,it,jg,ig)%fpp(1) =  fx(:,:)
            angular_table(:,:,it,jg,ig)%fpp(2) =  fy(:,:)
            angular_table(:,:,it,jg,ig)%fpp(3) = fxy(:,:)
            !
        end do
    end do
end do
end subroutine setup_angular_spline
!-------------------------------------------------------------------------------
subroutine get_spline_derivatives(n_bin, f, x, fpp)
!-------------------------------------------------------------------------------
integer, intent(in) :: n_bin
real(dp), intent(in) :: f(n_bin), x(n_bin)
real(dp), intent(out) :: fpp(n_bin)

real(dp) :: u(n_bin), p, sig
integer :: i

fpp(1:n_bin) = 0.0d0
u(1) = 0.0d0

do i = 2, n_bin-1
    sig = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
    p = sig*fpp(i-1) + 2.0d0
    fpp(i) = (sig-1.0d0) / p
    u(i) = (6.0d0*((f(i+1) - f(i))/(x(i+1)-x(i)) - (f(i)-f(i-1))/(x(i)-x(i-1))) / &
            (x(i+1)-x(i-1)) - sig*u(i-1)) / p
end do

do i = n_bin-1, 1, -1
    fpp(i) = fpp(i)*fpp(i+1) + u(i)
end do

end subroutine get_spline_derivatives
!-------------------------------------------------------------------------------
subroutine fill_inter_neigh(ig,jg,it,ir,ia, f)
!-------------------------------------------------------------------------------
integer, intent(in) :: ig, jg, it, ir, ia
real(dp), intent(out) :: f(4,4)

f(1,1) = angular_table(ia,  ir,  it,jg,ig)%f
f(1,2) = angular_table(ia+1,ir,  it,jg,ig)%f
f(1,3) = angular_table(ia,  ir+1,it,jg,ig)%f
f(1,4) = angular_table(ia+1,ir+1,it,jg,ig)%f

f(2:4,1) = angular_table(ia,  ir,  it,jg,ig)%fpp(1:3)
f(2:4,2) = angular_table(ia+1,ir,  it,jg,ig)%fpp(1:3)
f(2:4,3) = angular_table(ia,  ir+1,it,jg,ig)%fpp(1:3)
f(2:4,4) = angular_table(ia+1,ir+1,it,jg,ig)%fpp(1:3)

end subroutine fill_inter_neigh
!-------------------------------------------------------------------------------
subroutine update_dplanes(Xn, Xd)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: Xn(3,tn%stdatm), Xd(3, tn%stdatm)
integer :: i_atm

do i_atm = 1, tn%stdatm
    if (i_G(i_atm)%disabled) cycle
    call get_dplane(i_atm, xn(1:3,i_atm), xd(1:3,i_atm))
end do

end subroutine update_dplanes
!-------------------------------------------------------------------------------
subroutine get_dplane(i_atm, xn,xd)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_atm
real(dp), intent(out) :: xn(3), xd(3)
real(dp) :: v12(3), v13(3), vtmp(3)

v12(1:3) = R(1:3, i_G(i_atm)%i2) - R(1:3, i_atm)
v13(1:3) = R(1:3, i_G(i_atm)%i3) - R(1:3, i_atm)
if (i_G(i_atm)%two_bonds) then
    xn = v12 + v13
else
    xn = v12
end if

call cross(xn, v13, vtmp)
call cross(vtmp, xn, xd)

call v_norm(xd)
call v_norm(xn)

end subroutine get_dplane
!-------------------------------------------------------------------------------
subroutine calc_goap_ang(v1,v2, cs)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: v1(3), v2(3)
real(dp), intent(out) :: cs

cs = dot_product(v1,v2)

end subroutine calc_goap_ang
!-------------------------------------------------------------------------------
subroutine calc_goap_phi(v1,v2,v3, phi)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: v1(3), v2(3), v3(3)
real(dp), intent(out) :: phi
real(dp) :: cos1, cos2, cos3, v4(3), v5(3)

cos1 = dot_product(v1,v3)

v4 = v3 - cos1*v1
call v_norm(v4)
cos2 = min(1.0d0, max(-1.0d0, dot_product(v4,v2)))

call cross(v1,v2,v5)
cos3 = dot_product(v5,v4)

phi = acos(cos2)
if (cos3 < 0.0d0) phi = -phi

end subroutine calc_goap_phi
!-------------------------------------------------------------------------------
END MODULE GOAP
!-------------------------------------------------------------------------------
