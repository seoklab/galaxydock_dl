!-------------------------------------------------------------------------------
! Copyright (C) 2020, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/restraint/bb_torsion.f90
!
! Description: Subroutines related to predicted ramachandran preference
!
!-------------------------------------------------------------------------------
MODULE BB_TORSION

use globals
use in_out_utils, only: find_atom_idx
use logger, only: terminate_with_error, log_p
use string, only: parse_longstring
use mathfunctions, only: cross, bicubic_interpolation, multiple_binormal, bound_ang
use ramachandran,  only: fill_inter_neigh
!
use energy_vars, only: energy_type, R, ii_R, bb_torsion_file, &
                       use_bb_torsion
use ramachandran_score, only: calc_phipsi_using_R

implicit none
save
private

! Parameters for setting up arrays
integer, parameter :: max_rsrcolumn = 20   ! Maximum # of parameters of a restraint
integer, parameter :: max_rsrlinelen = 200 ! Maximum length of a line in restraint file
integer, parameter :: n_bin = 18
real(dp), parameter :: first_ang = -180.d0 * deg2rad
real(dp), parameter :: ang_interval = 20.d0 * deg2rad
!real(dp), parameter :: last_ang = first_ang + (n_bin-1)*ang_interval

! General info of all the restraints is stored in type rsr,
! while parameters are stored separately depending on the restraint types
real(dp), allocatable :: bbtor_table(:,:,:,:) ! 4, 0:18, 0:18, n_res

public :: initialize_bb_torsion
public :: finalize_bb_torsion
public :: calc_bb_torsion_energy
public :: bbtor_table

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_bb_torsion(molecule)
!-------------------------------------------------------------------------------
! Initialize backbone phi/psi torsion potential
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule

! For the case using only segment fixing
if (trim(bb_torsion_file) == '') return

!Read bb_torsion variables
call read_bb_torsion_input(molecule, bb_torsion_file)

end subroutine initialize_bb_torsion
!-------------------------------------------------------------------------------
subroutine finalize_bb_torsion()
!-------------------------------------------------------------------------------
if (allocated(bbtor_table)) then
    deallocate(bbtor_table)
end if

end subroutine finalize_bb_torsion
!-------------------------------------------------------------------------------
subroutine read_bb_torsion_input(molecule, file_name)
!-------------------------------------------------------------------------------
! Read the phi/psi preference input file
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule
character(len=len_fname), intent(in) :: file_name
integer :: f_unit, ioerror, num_word, openstat
integer :: i_res, i_psi, i_phi
character(len=len_fname) :: word(max_rsrcolumn)
character(len=max_rsrlinelen) :: line

f_unit = 39
open(f_unit, file = trim(file_name), iostat=openstat)
if (openstat > 0) then
    call terminate_with_error('Terminate with error: No bb_torsion file found. Check {bb_torsion_file}.')
end if  

allocate(bbtor_table(4, 0:n_bin, 0:n_bin, molecule%n_res))
bbtor_table = 0.0d0

i_res = 0
i_psi = -1
do
    read(f_unit,"(A200)", iostat = ioerror) line
    if (ioerror < 0) exit

    call parse_longstring(line, num_word, word, max_rsrlinelen)
    if (num_word < 2) exit
    
    if (word(1) == 'RESIDUE') then
        read(word(2),"(I)") i_res
        i_psi = 0
    else
        ! potential parameter
        do i_phi = 0, 17
            read(word(i_phi+1),"(F)") bbtor_table(1, i_phi, i_psi, i_res)
        end do
        i_psi = i_psi + 1
    end if
end do
if (i_psi == -1) then
    call terminate_with_error('Terminate with error: No bb_torsion file found. Check {bb_torsion_file}.')
end if  

close(f_unit)

! table parsing
bbtor_table(1,:,18,:) = bbtor_table(1,:,0,:)
bbtor_table(1,18,:,:) = bbtor_table(1,0,:,:)

call get_derivative_bbtor_table(bbtor_table, molecule%n_res)

end subroutine read_bb_torsion_input
!-------------------------------------------------------------------------------
subroutine get_derivative_bbtor_table(bbtor_table, n_res)
!-------------------------------------------------------------------------------
! get derivatives using finite difference of neighbors
!-------------------------------------------------------------------------------
integer, intent(in) :: n_res
real(dp), intent(inout) :: bbtor_table(4, 0:18, 0:18, n_res)

real(dp) :: tmp(2,0:17,0:17,n_res)

! dphi
bbtor_table(2,1:17,:,:) = (bbtor_table(1,2:18,:,:) - bbtor_table(1,0:16,:,:)) &
                         / (2.d0 * ang_interval)
bbtor_table(2,0,:,:) = (bbtor_table(1,1,:,:) - bbtor_table(1,17,:,:)) &
                      / (2.d0 * ang_interval)
bbtor_table(2,18,:,:) = bbtor_table(2,0,:,:)

! dpsi
bbtor_table(3,:,1:17,:) = (bbtor_table(1,:,2:18,:) - bbtor_table(1,:,0:16,:)) &
                         / (2.d0 * ang_interval)
bbtor_table(3,:,0,:) = (bbtor_table(1,:,1,:) - bbtor_table(1,:,17,:)) &
                      / (2.d0 * ang_interval)
bbtor_table(3,:,18,:) = bbtor_table(3,:,0,:)

! dphi -> dpsi
tmp(1,1:17,:,:) = (bbtor_table(2,2:18,:,:) - bbtor_table(2,0:16,:,:)) &
                   / (2.d0 * ang_interval)
tmp(2,0,:,:) = (bbtor_table(2,1,:,:) - bbtor_table(2,17,:,:)) &
                / (2.d0 * ang_interval)

! dpsi -> dphi
tmp(2,:,1:17,:) = (bbtor_table(3,:,2:18,:) - bbtor_table(3,:,0:16,:)) &
                   / (2.d0 * ang_interval)
tmp(2,:,0,:) = (bbtor_table(3,:,1,:) - bbtor_table(3,:,17,:)) &
                / (2.d0 * ang_interval)

bbtor_table(4,0:17,0:17,:) = sum(tmp, dim=1) / 2.0d0
bbtor_table(4,18,:,:) = bbtor_table(4,0,:,:)
bbtor_table(4,:,18,:) = bbtor_table(4,:,0,:)

end subroutine get_derivative_bbtor_table
!-------------------------------------------------------------------------------
subroutine calc_bb_torsion_energy(f, g, use_res, calc_g)
!-------------------------------------------------------------------------------
! Main subroutine for backbone phi/psi preference potential
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f
real(dp), intent(out) :: g(3,tn%atom)
logical, intent(in) :: use_res(tn%stdres), calc_g

integer :: i_res, i_atm, i_ref, i_phi, i_psi, atmno
real(dp) :: phipsi(2), dphipsi(3,4,2)
real(dp) :: ff(4,4), x(2), arg(3)

f = 0.0d0
g(:,:) = 0.0d0

do i_res = 2, tn%stdres-1
    if (.not. use_res(i_res)) cycle
    !
    i_ref = res_index(i_res)%ref_res_no
    if (ref_res(i_ref)%ter_type == 'N' .or. &
        ref_res(i_ref)%ter_type == 'C') cycle

    call calc_phipsi_using_R(i_res, phipsi, dphipsi, calc_g)
    i_phi = max(0, min(n_bin-1, int((phipsi(1)-first_ang)/ang_interval)))
    i_psi = max(0, min(n_bin-1, int((phipsi(2)-first_ang)/ang_interval)))
    !
    ff(:,1) = bbtor_table(:,i_phi,  i_psi,  i_res)
    ff(:,2) = bbtor_table(:,i_phi+1,i_psi,  i_res)
    ff(:,4) = bbtor_table(:,i_phi  ,i_psi+1,i_res)
    ff(:,3) = bbtor_table(:,i_phi+1,i_psi+1,i_res)

    x(1) = phipsi(1) - first_ang - dble(i_phi)*ang_interval
    x(2) = phipsi(2) - first_ang - dble(i_psi)*ang_interval
    !
    arg(:) = bicubic_interpolation(x(1), x(2), ang_interval, ang_interval, ff)
    !
    f = f + arg(1)
    if (calc_g) then
        do i_atm = 1, 4
            atmno = res_index(i_res)%phi_atmid(i_atm,1)
            g(:,atmno) = g(:,atmno) + arg(2)*dphipsi(:,i_atm,1)
            atmno = res_index(i_res)%psi_atmid(i_atm,1)
            g(:,atmno) = g(:,atmno) + arg(3)*dphipsi(:,i_atm,2)
        end do
    end if
end do

end subroutine calc_bb_torsion_energy
!-------------------------------------------------------------------------------
END MODULE BB_TORSION
