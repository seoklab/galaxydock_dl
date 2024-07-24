!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/utils/ramachandran.f90
!
! Description:
!-------------------------------------------------------------------------------
MODULE RAMACHANDRAN
!-------------------------------------------------------------------------------
use globals
use logger
use mathfunctions, only: cross, bicubic_interpolation, multiple_binormal, bound_ang
use string,        only: parse_string, parse_longstring
use geometry,      only: calc_tor_and_grad
use ran, only: random

implicit none
save
private

logical :: use_rama
!
character(len=len_fname) :: rama_mode
character(len=len_fname) :: infile_rama_score(3) = ''
character(len=len_fname) :: infile_psipred, infile_psrm

logical :: psrm_mode
integer, parameter :: max_templ = 50
integer, parameter :: max_psrmlinelen = 3000

type psrm_type
  ! Residue type: 1:Standard; 2:GLY; 3:PRO 
  integer :: n_rama                           ! Number of rama maps for a residue
  real(dp) :: w(max_templ)                    ! Weight
  real(dp) :: x(2,max_templ)                  ! phi/psi mean
  real(dp) :: s(2,max_templ)                  ! phi/psi stdev
  real(dp) :: rho(max_templ)                  ! Corr. btw phi and psi
  real(dp) :: max_p
end type psrm_type

integer, parameter :: n_rama_table = 81
real(dp), parameter :: delta_angle = 10.0d0
real(dp), parameter :: delta_angle_rad = delta_angle*deg2rad

character(len=4), allocatable :: ndrama_table_name(:)
real(dp), allocatable :: protein_rama(:,:,:,:)

logical :: use_psipred
real(dp), allocatable :: protein_ss_prob(:,:)

!-------------------------------------------------------------------------------
public :: use_rama
public :: rama_mode
!
public :: initialize_ramachandran
public :: finalize_ramachandran
!
public :: protein_rama
public :: protein_ss_prob
!
public :: fill_inter_neigh
!
public :: calc_single_residue_rama
public :: calc_fragment_rama
public :: phipsi_consistency_to_psipred
public :: ss_consistency_to_psipred
public :: random_sample_rama
!
public :: infile_psipred, infile_psrm
public :: delta_angle_rad

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_ramachandran(protein)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein

logical :: use_table   = .false.
logical :: use_ndrama  = .false.

type(psrm_type), allocatable :: psrm(:)
real(dp), allocatable :: rama_table(:,:,:,:)

allocate(protein_rama(4, 0:36, 0:36, protein%n_res))

use_psipred = .false.
if (trim(rama_mode) == 'psipred') then
    write(log_msg, "(A)") '  Ramachandran score mode is set to "psipred".'
    use_psipred = .true.
else if (trim(rama_mode) == 'hybrid') then
    write(log_msg, "(A)") '  Ramachandran score mode is set to "hybrid".'
    use_table   = .true.
    use_psipred = .true.
    !
    infile_rama_score(1) = trim(data_dir) // '/rama_standard_18_mod.prob'
    infile_rama_score(2) = trim(data_dir) // '/rama_standard_GLY_mod.prob'
    infile_rama_score(3) = trim(data_dir) // '/rama_standard_PRO_mod.prob'
else if (trim(rama_mode) == 'standard') then
    write(log_msg, "(A)") '  Ramachandran score mode is set to "standard".'
    use_table   = .true.
    !
    infile_rama_score(1) = trim(data_dir) // '/rama_standard_18_mod.prob'
    infile_rama_score(2) = trim(data_dir) // '/rama_standard_GLY_mod.prob'
    infile_rama_score(3) = trim(data_dir) // '/rama_standard_PRO_mod.prob'
else if (trim(rama_mode) == 'correct') then
    write(log_msg, "(A)") '  Ramachandran score mode is set to "correct".'
    use_table   = .true.
    !
    infile_rama_score(1) = trim(data_dir) // '/rama_correct_18.prob'
    infile_rama_score(2) = trim(data_dir) // '/rama_correct_GLY.prob'
    infile_rama_score(3) = trim(data_dir) // '/rama_correct_PRO.prob'
else if (trim(rama_mode) == 'neighbor') then
    use_table   = .true.
    use_ndrama  = .true.
    if (trim(infile_psipred) == '') then
        use_psipred = .false.
        write(log_msg, "(A)") '  Ramachandran score mode is set to "neighbor-dependent".'
    else
        use_psipred = .true.
        write(log_msg, "(A)") '  Ramachandran score mode is set to "neighbor,SS-dependent".'
    end if
    !
    infile_rama_score(1) = trim(data_dir) // '/rama_soluble_01Apr14.prob'
    infile_rama_score(2) = trim(data_dir) // '/rama_soluble_01Apr14.list'
else
    write(log_msg, "(A)") "Error: Wrong rama_mode; rama_mode can be one of"
    call log_p(log_msg, me=me, level=00)
    write(log_msg, "(A)") "       neighbor, standard, correct, hybrid, psipred"
    call terminate_with_error(log_msg)
end if
call log_p(log_msg, me=me, level=30)

if (trim(infile_psrm) /= '') then
    psrm_mode = .true.
else
    psrm_mode = .false.
end if

if (use_table) then
    allocate(rama_table(4, 0:36, 0:36, n_rama_table))
    if (use_ndrama) then
        allocate(ndrama_table_name(n_rama_table))
        call read_ndrama_score(infile_rama_score(1), rama_table)
    else
        call read_rama_score(infile_rama_score(1), rama_table(:,:,:,1))
        call read_rama_score(infile_rama_score(2), rama_table(:,:,:,2))
        call read_rama_score(infile_rama_score(3), rama_table(:,:,:,3))
    end if
end if

allocate(protein_ss_prob(3,protein%n_res))
if (use_psipred) then
    call read_psipred(infile_psipred, protein%n_res, protein_ss_prob)
end if

if (use_ndrama) then
    if (.not. use_psipred) call read_ndrama_ss_list(infile_rama_score(2), protein, protein_ss_prob)
    call setup_protein_ndrama(protein, protein_ss_prob, rama_table)
else
    if (use_psipred) then
        call setup_protein_psipred_score(protein, protein_ss_prob)
    end if
    if (use_table) then
        call setup_protein_rama(protein, rama_table(:,:,:,1:3), use_psipred)
    end if
end if

if (psrm_mode) then
    write(log_msg, "(A)") "  PSRM mode is turned on."
    call log_p(log_msg, me=me, level=30)

    allocate(psrm(protein%n_res))
    call read_psrm_score(protein, psrm)
    call update_psrm_rama_score(protein, psrm)
    deallocate(psrm)
end if

if (use_table)   deallocate(rama_table)
if (use_ndrama)  deallocate(ndrama_table_name)
if (.not. use_psipred) deallocate(protein_ss_prob)

end subroutine initialize_ramachandran
!-------------------------------------------------------------------------------
subroutine finalize_ramachandran()
!-------------------------------------------------------------------------------
if (use_rama)    deallocate(protein_rama)
if (use_psipred) deallocate(protein_ss_prob)

end subroutine finalize_ramachandran
!-------------------------------------------------------------------------------
subroutine read_rama_score(infile_rama, rama_table)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: infile_rama
real(dp), intent(out) :: rama_table(4,0:36,0:36)

character(len=len_fname) :: line, word(6)
integer, parameter :: f_unit = 35
integer :: io, n_word, i_phi, i_psi
real(dp) :: f, g(3), phi, psi

open(f_unit, file=trim(infile_rama), iostat=io)
if (io /= 0) then
    call terminate_with_error("ERROR: No file found for ramachandran score")
end if

do 
    read(f_unit, "(A120)", iostat=io) line
    if (io < 0) exit

    call parse_string(line, n_word, word)
    read(word(1), "(F10.1)") phi
    read(word(2), "(F10.1)") psi
    read(word(3), "(F10.4)") f
    read(word(4), "(F10.4)") g(1)
    read(word(5), "(F10.4)") g(2)
    read(word(6), "(F10.4)") g(3)

    i_phi = int((phi+180.0d0)/delta_angle)
    i_psi = int((psi+180.0d0)/delta_angle)
    g(1:3) = g(1:3) / deg2rad

    rama_table(1, i_phi, i_psi) = f
    rama_table(2:4, i_phi, i_psi) = g(1:3)

end do

close(f_unit)

end subroutine read_rama_score
!-------------------------------------------------------------------------------
subroutine read_ndrama_score(infile_rama, ndrama_table)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: infile_rama
real(dp), intent(out) :: ndrama_table(4,0:36,0:36,81)

integer, parameter :: f_unit = 35
integer :: io, n_word
character(len=len_fname) :: line, word(6)
integer :: i_phi, i_psi, i_table
real(dp) :: f, g(3), phi, psi

open(f_unit, file=trim(infile_rama), iostat=io)
if (io /= 0) then
    call terminate_with_error("ERROR: No file found for neighbor,SS-dependent ramachandran score")
end if

i_table = 0
do
    read(f_unit, '(A120)', iostat = io) line
    if (io /= 0) exit

    if (line(1:3) == 'RES') then
        i_table = i_table + 1
        call parse_string(line, n_word, word)
        write(ndrama_table_name(i_table),"(A3,A1)") word(2)(1:3), word(4)(1:1)
    else if (line(1:3) == 'TER') then
        cycle
    else
        call parse_string(line, n_word, word)
        read(word(1), "(F6.1)") phi
        read(word(2), "(F6.1)") psi
        read(word(3), "(F12.7)") f
        read(word(4), "(F12.7)") g(1)
        read(word(5), "(F12.7)") g(2)
        read(word(6), "(F12.7)") g(3)

        i_phi = int((phi+180.0d0)/delta_angle)
        i_psi = int((psi+180.0d0)/delta_angle)
        g(1:3) = g(1:3) / deg2rad

        ndrama_table(1,   i_phi, i_psi, i_table) = f
        ndrama_table(2:4, i_phi, i_psi, i_table) = g(1:3)

    end if
end do

close(f_unit)

end subroutine read_ndrama_score
!-------------------------------------------------------------------------------
subroutine read_psrm_score(protein, psrm)
!-------------------------------------------------------------------------------
! Read Position-specific Ramachandran Map (PSRM)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(psrm_type), intent(out) :: psrm(:)
integer :: f_unit=40, io, n_word
integer :: i_res, i_rama
real(dp) :: w, phi, psi, max_p_i
character(len=max_psrmlinelen) :: line
character(len=len_fname) :: word(100)
character(len=1) :: ss_type

do i_res = 2, protein%n_res-1
    if (protein%residue(i_res)%ter_type == 'N' .or. &
        protein%residue(i_res)%ter_type == 'C')  cycle
    psrm(i_res)%n_rama = 0
end do

open(f_unit, file=trim(infile_psrm), iostat=io)
if (io /= 0) then
    call terminate_with_error('Terminate with error: No psrm file found. Check {infile_psrm}.')
end if

do
    read(f_unit, "(A3000)", iostat=io) line
    if (io < 0) exit
    if (line(1:1) == '#' .or. line(1:1) == '!') cycle

    call parse_longstring(line, n_word, word, max_psrmlinelen)
    if (n_word < 5) cycle
    read(word(1), "(I10)") i_res
    read(word(2), '(I10)') psrm(i_res)%n_rama
    psrm(i_res)%max_p = 0.0d0

    do i_rama = 1, psrm(i_res)%n_rama
        read(word((i_rama-1)*4 + 3),'(A1)') ss_type
        read(word((i_rama-1)*4 + 4),'(F10.5)') w
        read(word((i_rama-1)*4 + 5),'(F10.5)') phi
        read(word((i_rama-1)*4 + 6),'(F10.5)') psi

        psrm(i_res)%w(i_rama) = w
        psrm(i_res)%x(1,i_rama) = phi*deg2rad
        psrm(i_res)%x(2,i_rama) = psi*deg2rad
        psrm(i_res)%rho(i_rama) = 0.0d0

        if (trim(ss_type) == 'H') then
            psrm(i_res)%s(1:2,i_rama) = (/15.0d0, 15.0d0/)*deg2rad
        else if (trim(ss_type) == 'E') then                        
            psrm(i_res)%s(1:2,i_rama) = (/15.0d0, 25.0d0/)*deg2rad
        else if (trim(ss_type) == 'L') then                        
            psrm(i_res)%s(1:2,i_rama) = (/25.0d0, 25.0d0/)*deg2rad
        else if (trim(ss_type) == 'G') then                        
            psrm(i_res)%s(1:2,i_rama) = (/60.0d0, 60.0d0/)*deg2rad
        else if (trim(ss_type) == 'P') then                        
            psrm(i_res)%s(1:2,i_rama) = (/15.0d0, 15.0d0/)*deg2rad
        else                                                       
            psrm(i_res)%s(1:2,i_rama) = (/60.0d0, 60.0d0/)*deg2rad
        end if

        max_p_i = two_pi * psrm(i_res)%s(1,i_rama) * psrm(i_res)%s(2,i_rama)
        max_p_i = max_p_i * sqrt(1.0d0-psrm(i_res)%rho(i_rama)**2)
        psrm(i_res)%max_p = psrm(i_res)%max_p + psrm(i_res)%w(i_rama) / max_p_i

    end do
end do

close(f_unit)

end subroutine read_psrm_score
!-------------------------------------------------------------------------------
subroutine read_psipred(file_name, n_res, ss_prob)
!-------------------------------------------------------------------------------
! read psipred
!-------------------------------------------------------------------------------
  character(len=len_fname), intent(in) :: file_name
  integer, intent(in) :: n_res
  real(dp), intent(out) :: ss_prob(:,:)
  integer :: f_unit, ioerror, num_word, openstat, resno
  character(len=len_fname) :: word(10)
  character(len=len_fname) :: line

  f_unit = 40
  open(f_unit, file = trim(file_name), iostat=openstat)

  if (openstat > 0) then
     call terminate_with_error('Terminate with error: No psipred file found. Check {psipred_file}.')
  end if

  do
     read(f_unit,"(A120)", iostat = ioerror) line
     if (ioerror < 0) exit

     call parse_string(line, num_word, word)

     if (num_word /= 6) cycle

     read(word(1),"(I10)") resno
     if (resno > n_res) then
         close(f_unit)
         write(log_msg,'(A)') "Terminate with error: The number of sequence exceeds tnr."
         call terminate_with_error(log_msg)
     end if

     read(word(4),"(F10.4)") ss_prob(1,resno)
     read(word(5),"(F10.4)") ss_prob(2,resno)
     read(word(6),"(F10.4)") ss_prob(3,resno)
  end do

  close(f_unit)

end subroutine read_psipred
!-------------------------------------------------------------------------------
subroutine update_psrm_rama_score(protein, psrm)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(psrm_type), intent(in) :: psrm(:)
integer :: i_res, i_phi, i_psi, n
real(dp) :: x(2), f(4), f0(4), f1(4)

do i_res = 2, protein%n_res-1
    if (protein%residue(i_res)%ter_type == 'N' .or. &
        protein%residue(i_res)%ter_type == 'C')  cycle

    n = psrm(i_res)%n_rama
    if (n == 0) cycle

    do i_psi = 0, 36
        x(2) = dble(i_psi-18)*delta_angle_rad
        do i_phi = 0, 36
            x(1) = dble(i_phi-18)*delta_angle_rad
            f0 = protein_rama(1:4,i_phi,i_psi,i_res)
            !
            f1 = multiple_binormal_ext(n, x(1:2), psrm(i_res)%x(1:2,1:n), &
                                       psrm(i_res)%w(1:n), psrm(i_res)%s(1:2,1:n),&
                                       psrm(i_res)%rho(1:n), .false., .true., .true.)
            f1 = f1 / psrm(i_res)%max_p
            !
            f(1) = f0(1)*f1(1)
            f(2) = f0(1)*f1(2) + f0(2)*f1(1)
            f(3) = f0(1)*f1(3) + f0(3)*f1(1)
            f(4) = f0(1)*f1(4) + f0(2)*f1(3) + f0(3)*f1(2) + f0(4)*f1(1)
            !
            protein_rama(1:4,i_phi,i_psi,i_res) = f(1:4)
        end do
    end do
end do

end subroutine update_psrm_rama_score
!-------------------------------------------------------------------------------
subroutine read_ndrama_ss_list(infile_rama_score, protein, protein_ss_prob)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: infile_rama_score
type(molecule_type), intent(in) :: protein
real(dp), intent(out) :: protein_ss_prob(:,:)
integer :: io, f_unit = 77, n_word
integer :: i_res, i, i_type
character(len=len_fname) :: line, word(4)
character(len=3) :: resType_s(27), resType
real(dp) :: prob_s(3,27)

open(f_unit, file=trim(infile_rama_score), iostat=io)
if (io /= 0) then
    call terminate_with_error("Terminate with error: No NDRAMA SS preference file.")
end if

i = 0
do 
    read(f_unit, "(A120)", iostat=io) line
    if (io /= 0) exit

    i = i + 1
    call parse_string(line, n_word, word)
    resType_s(i) = word(1)(1:3)
    read(word(2),"(f7.5)") prob_s(1,i)
    read(word(3),"(f7.5)") prob_s(2,i)
    read(word(4),"(f7.5)") prob_s(3,i)
end do

close(f_unit)

protein_ss_prob(:,:) = 0.0d0

do i_res = 2, protein%n_res-1
    if (protein%residue(i_res)%ter_type == 'N' .or. &
        protein%residue(i_res)%ter_type == 'C')  cycle
    do i = 1, 3
        if (protein%residue(i_res-2+i)%pdb_res_name == 'PRO') then
            resType(i:i) = 'P'
        else if (protein%residue(i_res-2+i)%pdb_res_name == 'GLY') then
            resType(i:i) = 'G'
        else
            resType(i:i) = 'X'
        end if
    end do

    do i_type = 1, 27
        if (resType == resType_s(i_type)) exit
    end do

    protein_ss_prob(1:3,i_res) = prob_s(1:3,i_type)
end do

end subroutine read_ndrama_ss_list
!-------------------------------------------------------------------------------
subroutine setup_protein_psipred_score(protein, protein_ss_prob)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
real(dp), intent(in) :: protein_ss_prob(:,:)
integer :: i_res, i_phi, i_psi, n
real(dp) :: x(2), f(4)
real(dp) :: w(5), x0(2,5), sigma(2,5), rho(5), probs(2), confidence
!
real(dp), parameter :: x_h(1:2) = (/ -65.0d0,  -41.0d0 /)*deg2rad
real(dp), parameter :: x_e(1:2) = (/-130.0d0,  135.0d0 /)*deg2rad
real(dp), parameter :: x_i(1:2) = (/  65.0d0,   41.0d0 /)*deg2rad
real(dp), parameter :: x_l(1:2) = (/ 130.0d0,  180.0d0 /)*deg2rad
real(dp), parameter :: x_p(1:2) = (/ -65.0d0,  140.0d0 /)*deg2rad
real(dp), parameter :: x_t(1:2) = (/ -80.0d0,   20.0d0 /)*deg2rad
!
real(dp), parameter :: s_h(1:2) = (/   7.5d0,    7.5d0 /)*deg2rad
real(dp), parameter :: s_e(1:2) = (/  25.0d0,   20.0d0 /)*deg2rad
real(dp), parameter :: s_i(1:2) = (/  10.0d0,   10.0d0 /)*deg2rad
real(dp), parameter :: s_l(1:2) = (/  15.0d0,   15.0d0 /)*deg2rad
real(dp), parameter :: s_p(1:2) = (/  10.0d0,   10.0d0 /)*deg2rad
real(dp), parameter :: s_t(1:2) = (/  10.0d0,   10.0d0 /)*deg2rad
!
real(dp), parameter :: rho_h = -0.3d0
real(dp), parameter :: rho_e = 0.0d0
!
real(dp), parameter :: c_xhh  = 0.8d0
real(dp), parameter :: c_xhi  = 0.2d0
real(dp), parameter :: c_xch  = 0.0d0
real(dp), parameter :: c_xce  = 0.5d0 
real(dp), parameter :: c_xct  = 0.3d0 
real(dp), parameter :: c_xci  = 0.2d0 
!
real(dp), parameter :: c_ghh = 1.0d0/3.0d0
real(dp), parameter :: c_gee = 1.0d0/3.0d0
real(dp), parameter :: c_gch = 0.0d0
real(dp), parameter :: c_gce = 0.0d0
real(dp), parameter :: c_gcl = 1.0d0/3.0d0
real(dp), parameter :: c_gci = 1.0d0/3.0d0
real(dp), parameter :: c_gct = 1.0d0/3.0d0
!
real(dp), parameter :: c_phh = 1.0d0/3.0d0
real(dp), parameter :: c_pch = 0.0d0
real(dp), parameter :: c_pcp = 1.0d0

do i_res = 2, protein%n_res-1
    if (protein%residue(i_res)%ter_type == 'N' .or. &
        protein%residue(i_res)%ter_type == 'C')  cycle

    w(:) = 0.0d0
    x0(:,:) = 0.0d0
    sigma(:,:) = 0.0d0

    if (protein%residue(i_res)%pdb_res_name == 'PRO') then
        n = 2
        !
        w(1) = c_pch*protein_ss_prob(1,i_res) + c_phh*protein_ss_prob(2,i_res)
        w(2) = c_pcp*protein_ss_prob(1,i_res)
        !
        rho(1) = rho_h
        rho(2) = 0.0d0
        !
        x0(:,1) = x_h(:)
        x0(:,2) = x_p(:)
        sigma(:,1) = s_h(:)
        sigma(:,2) = s_p(:)

    else if (protein%residue(i_res)%pdb_res_name == 'GLY') then
        n = 5
        !
        w(1) = c_gch*protein_ss_prob(1,i_res) + c_ghh*protein_ss_prob(2,i_res)
        w(2) = c_gce*protein_ss_prob(1,i_res) + c_gee*protein_ss_prob(3,i_res)
        w(3) = c_gcl*protein_ss_prob(1,i_res)
        w(4) = c_gci*protein_ss_prob(1,i_res)
        w(5) = c_gct*protein_ss_prob(1,i_res)
        !
        x0(:,1) = x_h(:)
        x0(:,2) = x_e(:)
        x0(:,3) = x_l(:)
        x0(:,4) = x_i(:)
        x0(:,5) = x_t(:)
        !
        sigma(:,1) = s_h(:)
        sigma(:,2) = s_e(:)
        sigma(:,3) = s_l(:)
        sigma(:,4) = s_i(:)
        sigma(:,5) = s_t(:)
        !
        rho(1) = rho_h
        rho(2) = rho_e
        rho(3:5) = 0.0d0

    else
        n = 4
        !
        w(1) = c_xch*protein_ss_prob(1,i_res) + c_xhh*protein_ss_prob(2,i_res)
        w(2) = c_xce*protein_ss_prob(1,i_res) +       protein_ss_prob(3,i_res)
        w(3) = c_xci*protein_ss_prob(1,i_res) + c_xhi*protein_ss_prob(2,i_res)
        w(4) = c_xct*protein_ss_prob(1,i_res)
        !
        x0(:,1) = x_h(:)
        x0(:,2) = x_e(:)
        x0(:,3) = x_i(:)
        x0(:,4) = x_t(:)
        !
        sigma(:,1) = s_h(:)
        sigma(:,2) = s_e(:)
        sigma(:,3) = s_i(:)
        sigma(:,4) = s_t(:)
        !
        rho(1) = rho_h
        rho(2) = rho_e
        rho(3:4) = 0.0d0

    end if

    probs = (/protein_ss_prob(2,i_res), protein_ss_prob(1,i_res)*0.5d0+protein_ss_prob(3,i_res)/)
    confidence = maxval(probs) - minval(probs)

    do i_psi = 0, 36
        x(2) = dble(i_psi-18)*delta_angle_rad
        do i_phi = 0, 36
            x(1) = dble(i_phi-18)*delta_angle_rad
            f(:) = multiple_binormal_ext(n, x(1:2), x0(1:2,1:n), w(1:n), &
                sigma(1:2,1:n), rho(1:n), .true., .true., .true.)*confidence
            protein_rama(1:4,i_phi,i_psi,i_res) = f(1:4)
        end do
    end do
end do

end subroutine setup_protein_psipred_score
!-------------------------------------------------------------------------------
subroutine setup_protein_rama(protein, rama_table, use_psipred)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
real(dp), intent(in) :: rama_table(4,0:36,0:36,3)
logical, intent(in) :: use_psipred
integer :: i_res, i_table

do i_res = 2, protein%n_res-1
    if (protein%residue(i_res)%ter_type == 'N' .or. &
        protein%residue(i_res)%ter_type == 'C')  cycle

    if (protein%residue(i_res)%pdb_res_name == 'PRO') then
        i_table = 3
    elseif (protein%residue(i_res)%pdb_res_name == 'GLY') then
        i_table = 2
    else
        i_table = 1
    end if
    if (.not. use_psipred) then
        protein_rama(:,:,:,i_res) = rama_table(:,:,:,i_table)
    else
        protein_rama(:,:,:,i_res) = protein_rama(:,:,:,i_res) + &
                                          0.05d0*rama_table(:,:,:,i_table)
    end if
end do

end subroutine setup_protein_rama
!-------------------------------------------------------------------------------
subroutine setup_protein_ndrama(protein, protein_ss_prob, rama_table)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
real(dp), intent(in) :: rama_table(:,:,:,:)
real(dp), intent(in) :: protein_ss_prob(:,:)
integer :: i_res, i, i_table
character(len=4) :: table_name
character(len=3) :: aa_name

do i_res = 2, protein%n_res-1
    if (protein%residue(i_res)%ter_type == 'N' .or. &
        protein%residue(i_res)%ter_type == 'C')  cycle
    do i = 1, 3
        if (protein%residue(i_res-2+i)%pdb_res_name == 'PRO') then
            aa_name(i:i) = 'P'
        else if (protein%residue(i_res-2+i)%pdb_res_name == 'GLY') then
            aa_name(i:i) = 'G'
        else
            aa_name(i:i) = 'X'
        end if
    end do
    write(log_msg, "(A,I4,1x,A3,A,A3)") "  Ramachandran map: ", i_res, &
        protein%residue(i_res)%pdb_res_name, " --> ", aa_name
    call log_p(log_msg, me=me, level=40)

    protein_rama(:,:,:,i_res) = 0.0d0

    write(table_name,'(A3,A1)') aa_name(1:3), 'C'
    call find_table_index(table_name, i_table)
    protein_rama(:,:,:,i_res) = protein_rama(:,:,:,i_res) + & 
                                      protein_ss_prob(1,i_res) * rama_table(:,:,:,i_table)

    write(table_name,'(A3,A1)') aa_name(1:3), 'H'
    call find_table_index(table_name, i_table)
    protein_rama(:,:,:,i_res) = protein_rama(:,:,:,i_res) + & 
                                      protein_ss_prob(2,i_res) * rama_table(:,:,:,i_table)

    write(table_name,'(A3,A1)') aa_name(1:3), 'E'
    call find_table_index(table_name, i_table)
    protein_rama(:,:,:,i_res) = protein_rama(:,:,:,i_res) + & 
                                      protein_ss_prob(3,i_res) * rama_table(:,:,:,i_table)

end do

end subroutine setup_protein_ndrama
!-------------------------------------------------------------------------------
subroutine find_table_index(table_name, i_table)
!-------------------------------------------------------------------------------
character(len=4), intent(in) :: table_name
integer, intent(out) :: i_table

do i_table = 1, 81
    if (ndrama_table_name(i_table) == table_name) return 
end do

end subroutine find_table_index
!-------------------------------------------------------------------------------
subroutine fill_inter_neigh(i_res, i_phi, i_psi, f)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_res, i_phi, i_psi
real(dp), intent(out) :: f(:,:)

f(:,1) = protein_rama(:,i_phi,  i_psi,  i_res)
f(:,2) = protein_rama(:,i_phi+1,i_psi,  i_res)
f(:,4) = protein_rama(:,i_phi  ,i_psi+1,i_res)
f(:,3) = protein_rama(:,i_phi+1,i_psi+1,i_res)

end subroutine fill_inter_neigh
!-------------------------------------------------------------------------------
subroutine calc_single_residue_rama(resNo, residue, f)
!-------------------------------------------------------------------------------
! resNo(2) shoule be resNo and resNo+1
!-------------------------------------------------------------------------------
integer, intent(in) :: resNo(2)
type(residue_type), intent(in) :: residue(2)
real(dp), intent(out) :: f

integer :: i_ref(2), i_phi, i_psi
real(dp) :: phipsi(2), ff(4,4), x(2), arg(3)

i_ref(1) = residue(1)%res_type
i_ref(2) = residue(2)%res_type

!Phi angle for residue(1)
phipsi(1) = bound_ang(residue(1)%t_ang(ref_res(i_ref(1))%iphi))
!Psi angle for residue(1)
phipsi(2) = bound_ang(residue(2)%t_ang(ref_res(i_ref(2))%ipsi))

i_phi = min(35, int((phipsi(1)+pi)/delta_angle_rad))
i_psi = min(35, int((phipsi(2)+pi)/delta_angle_rad))

x(1) = phipsi(1) + pi - dble(i_phi)*delta_angle_rad
x(2) = phipsi(2) + pi - dble(i_psi)*delta_angle_rad

call fill_inter_neigh(resNo(1), i_phi, i_psi, ff(:,:))
arg(:) = bicubic_interpolation(x(1), x(2), delta_angle_rad, delta_angle_rad, ff)
f = arg(1)

call fill_inter_neigh(resNo(2), i_phi, i_psi, ff(:,:))
arg(:) = bicubic_interpolation(x(1), x(2), delta_angle_rad, delta_angle_rad, ff)

f = (f + arg(1))/2.0d0

end subroutine calc_single_residue_rama
!-------------------------------------------------------------------------------
subroutine calc_fragment_rama(f, n_res, resNo, phipsi)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f
integer, intent(in) :: n_res, resNo(n_res)
real(dp), intent(in) :: phipsi(2,n_res)

integer :: i, i_res, i_ref, i_phi, i_psi
real(dp) :: ff(4,4), x(2), arg(3)
logical, parameter :: calc_g = .false.

f = 0.0d0
do i = 1, n_res
    i_res = resNo(i)
    i_ref = res_index(i_res)%ref_res_no
    if (ref_res(i_ref)%ter_type == 'N' .or. &
        ref_res(i_ref)%ter_type == 'C') cycle
    !
    i_phi = min(35, int((phipsi(1,i)+pi)/delta_angle_rad))
    i_psi = min(35, int((phipsi(2,i)+pi)/delta_angle_rad))
    !
    call fill_inter_neigh(i_res, i_phi, i_psi, ff(:,:))
    x(1) = phipsi(1,i) + pi - dble(i_phi)*delta_angle_rad
    x(2) = phipsi(2,i) + pi - dble(i_psi)*delta_angle_rad
    !
    arg(:) = bicubic_interpolation(x(1), x(2), delta_angle_rad, delta_angle_rad, ff)
    !
    f = f + arg(1)
end do

end subroutine calc_fragment_rama
!-------------------------------------------------------------------------------
subroutine random_sample_rama(i_res, rama_ang, rama_in)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_res
real(dp), intent(in), optional :: rama_in(2)
real(dp), intent(out) :: rama_ang(2)
real(dp) :: p_rama(36,36), p_norm, p_random, p_sum
real(dp) :: phipsi(2), sigma(2,1), w(1), rho(1), f_in(3)
integer :: i_phi, i_psi, i_max(2)

do i_phi = 1, 36
    do i_psi = 1, 36
        p_rama(i_phi,i_psi) = exp(-protein_rama(1, i_phi, i_psi, i_res))
    end do
end do
if (present(rama_in)) then
    sigma(1:2,1) = 9.0d0*delta_angle_rad
    w(1) = 1.0d0
    rho(1) = 0.0d0
    do i_phi = 1, 36
        phipsi(1) = dble(i_phi-18)*delta_angle_rad
        do i_psi = 1, 36
            phipsi(2) = dble(i_psi-18)*delta_angle_rad
            f_in = multiple_binormal(1, phipsi, rama_in, w, sigma, rho, .false., .false.)
            p_rama(i_phi,i_psi) = p_rama(i_phi, i_psi) * f_in(1)
        end do
    end do
end if
p_norm = sum(p_rama)
p_rama = p_rama / p_norm
!
p_random = random()
p_sum = 0.0d0
do i_phi = 1, 36
    do i_psi = 1, 36
        p_sum = p_sum + p_rama(i_phi,i_psi)
        if (p_sum > p_random) then
            rama_ang(1) = (i_phi+0.5)*delta_angle_rad-pi
            rama_ang(2) = (i_psi+0.5)*delta_angle_rad-pi
            return
        end if
    end do
end do

i_max = maxloc(p_rama)
rama_ang(1) = (i_max(1)+0.5)*delta_angle_rad-pi
rama_ang(2) = (i_max(2)+0.5)*delta_angle_rad-pi

end subroutine random_sample_rama
!-------------------------------------------------------------------------------
subroutine phipsi_consistency_to_psipred(resno, phi, psi, score, extend_correction_in)
!-------------------------------------------------------------------------------
integer, intent(in) :: resno
real(dp), intent(in) :: phi, psi
real(dp), intent(out) :: score
logical, intent(in), optional :: extend_correction_in
logical :: extend_correction
real(dp) :: w(1), rho(1), x(2), x0(2,1), sigma(2,1)
real(dp) :: f(3), score_h, score_e, score_c
real(dp) :: prob(3)

real(dp), parameter :: x_h(1:2) = (/ -65.0d0,  -41.0d0 /)*deg2rad
real(dp), parameter :: s_h(1:2) = (/   7.5d0,    7.5d0 /)*deg2rad
real(dp), parameter :: x_e(1:2) = (/-130.0d0,  135.0d0 /)*deg2rad
real(dp), parameter :: s_e(1:2) = (/  25.0d0,   20.0d0 /)*deg2rad

if (present(extend_correction_in)) then
    extend_correction = extend_correction_in
else
    extend_correction = .true.
end if

score = 0.0d0
rho(:) = 0.0d0
w(:) = 1.0d0
x(1:2) = (/phi, psi/)

if (extend_correction) then
    prob(1) = max(0.0d0,protein_ss_prob(1,resno) - 1.0d0/3.0d0)
    prob(2) = max(0.0d0,protein_ss_prob(2,resno) - 1.0d0/3.0d0)
    prob(3) = protein_ss_prob(3,resno)
else
    prob(1) = max(0.0d0,protein_ss_prob(1,resno) - 1.0d0/3.0d0)
    prob(2) = max(0.0d0,protein_ss_prob(2,resno) - 1.0d0/3.0d0)
    prob(3) = max(0.0d0,protein_ss_prob(3,resno) - 1.0d0/3.0d0)
end if

! Coil
score_c = prob(1)/3.0d0

! Helix
sigma(:,1) = s_h(1:2)*1.2d0
x0(:,1) = x_h(1:2)
f(:) = multiple_binormal(1, x(1:2), x0(1:2,1), w(1:1), sigma(1:2,1:1), rho(1:1), &
                        .false., .false.)

score_h = min(1.0d0,max(prob(2)-prob(3),0.0d0)*f(1))

! Extended
sigma(:,1) = s_e(1:2)*1.2d0
x0(1:2,1) = x_e(1:2)

f(:) = multiple_binormal(1, x(1:2), x0(1:2,1), w(1:1), sigma(1:2,1:1), rho(1:1), &
                         .false., .false.)
score_e = min(1.0d0,max(prob(3)-prob(2),0.0d0)*f(1))

! Sum
score = score_h + score_e + score_c

end subroutine phipsi_consistency_to_psipred
!-------------------------------------------------------------------------------
subroutine ss_consistency_to_psipred(resno, sec, score, extend_correction_in, SS_enhance_in)
!-------------------------------------------------------------------------------
integer, intent(in) :: resno
character(len=1), intent(in) :: sec
real(dp), intent(out) :: score
logical, intent(in), optional :: extend_correction_in, SS_enhance_in
logical :: extend_correction, SS_enhance
real(dp) :: prob(3)

if (present(extend_correction_in)) then
    extend_correction = extend_correction_in
else
    extend_correction = .true.
end if

if (present(SS_enhance_in)) then
    SS_enhance = SS_enhance_in
else
    SS_enhance = .false.
end if

score = 0.0d0
prob(:) = protein_ss_prob(:,resno)
if (extend_correction) then
    prob(3) = min(1.0d0,prob(3) + 0.5d0*prob(1))
end if

if (SS_enhance) then
    prob(1) = min(0.1d0,prob(1))!0.5d0*prob(1)
end if

if (sec == '_' .or. sec == 'T' .or. sec == 'S') then
    score = prob(1)
else if (sec == 'G' .or. sec == 'H' .or. sec == 'I') then
    score = prob(2)
else if (sec == 'E' .or. sec == 'B') then
    score = prob(3)+0.5d0*prob(1)
end if

end subroutine ss_consistency_to_psipred
!-------------------------------------------------------------------------------
function multiple_binormal_ext(n, x, x0, w, sigma, rho, get_log, calc_g, calc_h)
!-------------------------------------------------------------------------------
real(dp) :: multiple_binormal_ext(4)
integer, intent(in) :: n
real(dp), intent(in) :: x(2), x0(2,n), w(n)
real(dp), intent(in) :: sigma(2,n), rho(n)
logical, intent(in) :: get_log, calc_g, calc_h
integer :: i_m
real(dp) :: x1diff, x2diff, isig1, isig2
real(dp) :: sin1, sin2, cos1, cos2, ss1, ss2
real(dp) :: arg, tmp1, tmp2, tmp3
real(dp) :: g1_tmp, g2_tmp

multiple_binormal_ext(:) = 0.0d0

do i_m = 1, n
    x1diff = x(1) - x0(1,i_m)
    x2diff = x(2) - x0(2,i_m)
    isig1 = 1.0d0/sigma(1,i_m)
    isig2 = 1.0d0/sigma(2,i_m)

    sin1 = sin(x1diff)
    sin2 = sin(x2diff)
    cos1 = cos(x1diff)
    cos2 = cos(x2diff)
    tmp1 = 1.0d0 - rho(i_m)**2
    ss1 = sin1*isig1
    ss2 = sin2*isig2
    tmp2 = (1.0d0 - cos1)*isig1*isig1 - rho(i_m)*ss1*ss2 + &
           (1.0d0 - cos2)*isig2*isig2

    arg = w(i_m)*exp(-tmp2/tmp1)*isig1*isig2 / (two_pi*sqrt(tmp1))

    multiple_binormal_ext(1) = multiple_binormal_ext(1) + arg
    if (calc_g) then
        g1_tmp = ss1 - rho(i_m)*cos1*ss2
        g2_tmp = ss2 - rho(i_m)*cos2*ss1

        tmp3 = arg/tmp1
        multiple_binormal_ext(2) = multiple_binormal_ext(2) + tmp3*g1_tmp*isig1
        multiple_binormal_ext(3) = multiple_binormal_ext(3) + tmp3*g2_tmp*isig2

        if (calc_h) then
            multiple_binormal_ext(4) = multiple_binormal_ext(4) &
                + (arg*isig1*isig2/tmp1)*(rho(i_m)*cos1*cos2 + g1_tmp*g2_tmp/tmp1)
        end if
    end if
end do

if (get_log) then
    if (calc_g) then
        multiple_binormal_ext(2:3) = multiple_binormal_ext(2:3)/multiple_binormal_ext(1)
        if (calc_h) then
            multiple_binormal_ext(4) = -multiple_binormal_ext(4)/multiple_binormal_ext(1) &
                                       +multiple_binormal_ext(2)*multiple_binormal_ext(3)
        end if
    end if
    multiple_binormal_ext(1) = -log(multiple_binormal_ext(1))
else
    if (calc_g) multiple_binormal_ext(2:3) = -1.0d0*multiple_binormal_ext(2:3)
end if

end function multiple_binormal_ext
!-------------------------------------------------------------------------------
END MODULE RAMACHANDRAN
!-------------------------------------------------------------------------------
