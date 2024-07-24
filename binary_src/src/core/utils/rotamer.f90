!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/utils/rotamer.f90
!! Description:
!  - This module roles as 'hub' for rotamer modelings. 
!  - This module contains rotamer variables, as well as utilities for rotamer
!    operations. 
!  - Subroutines in this module can be classified as:
!    - Setup subroutines
!    - Utilities
!  - Initializing rotamer variables can be simply achieved by calling 
!    'setup_rotamer' (others are sub-components for setup_rotamer)
!  - Various utilities can facilitate your programming.
!       
!-------------------------------------------------------------------------------
MODULE ROTAMER
!-------------------------------------------------------------------------------
use globals
use logger,           only: log_p, terminate_with_error, log_thick_divider
use geometry,         only: internal2cartesian, cartesian2internal, &
                            internal2cartesian_sc, rotate_torsion_angle
use string,           only: parse_string
use ran,              only: random
use mathfunctions,    only: v_norm, cross, bound_ang
use convert_res_name, only: convert_to_stdres
!
use in_out_vars,      only: infile_rotamer_lib, infile_chi_def

implicit none
save
private

!===============================================================================
! PARAMETERS/VARIABLES
!===============================================================================
integer, parameter :: num_rot_res = 18  ! 20 - GLY,ALA
real(dp) :: max_rotamer_prob            ! rotamer cummulative probability cutoff
integer :: max_rotamer_state            ! maximum number of rotamer states
integer :: n_rot_res                    ! no. rotatable residues in the protein
  
!===============================================================================
! DERIVED TYPES
!===============================================================================
type rotamer_record_type
!-------------------------------------------------------------------------------
! TODO: Comment
!-------------------------------------------------------------------------------
integer  :: i_phi, i_psi
real(dp) :: prob
real(dp) :: chi(4)
real(dp) :: sigma(4)
!-------------------------------------------------------------------------------
end type rotamer_record_type
!-------------------------------------------------------------------------------
type rotamer_index_type
!-------------------------------------------------------------------------------
! TODO: Comment
!-------------------------------------------------------------------------------
character(len=4) :: res_name
character(len=4) :: chi_atoms(4,4)
integer :: start, end, n_chi, n_rottot
integer :: n_rot(-18:18, -18:18)
integer :: bbdep_id(2,-18:18,-18:18)
!-------------------------------------------------------------------------------
end type rotamer_index_type
!-------------------------------------------------------------------------------
type(rotamer_index_type),  dimension(:), allocatable :: rotamer_index  ! Index about i-type residue
type(rotamer_record_type), dimension(:), allocatable :: rotamer_record ! Full rotamer library record
!
public :: rotamer_index_type
public :: rotamer_record_type
!
public :: max_rotamer_state
public :: max_rotamer_prob
public :: n_rot_res
!
public :: initialize_rotamer
public :: finalize_rotamer
!
public :: get_rotamer_index
!
public :: place_rotamer
public :: rotate_chi_angles
public :: perturb_rotamer_state
public :: identify_rotamer_state
public :: identify_rotamer_state_single
!
public :: is_symmetric_chi
!
public :: rotamer_index
public :: rotamer_record

CONTAINS
!===============================================================================
subroutine initialize_rotamer()
!-------------------------------------------------------------------------------
! Initialize rotamer related arrays by reading protein info and rotamer library
!-------------------------------------------------------------------------------
call log_p('- Rotamer setup', me=me, level=20)
!
call read_rotamer_library()
call read_chi_definition()
!
call setup_ref_res_rotamer_index()
call setup_ref_res_chi_angle()
call setup_ref_res_chi_dep()

end subroutine initialize_rotamer
!-------------------------------------------------------------------------------
subroutine finalize_rotamer()
!-------------------------------------------------------------------------------
if (allocated(rotamer_record)) then
    deallocate(rotamer_record)
    deallocate(rotamer_index)
end if

end subroutine finalize_rotamer
!-------------------------------------------------------------------------------
subroutine read_rotamer_library()
!-------------------------------------------------------------------------------
! Read the Dunbrack's backbone-dependent rotamer library
!-------------------------------------------------------------------------------
integer :: n_rotamer

! Fill up rotamer_index array
call read_rotamer_library_index(n_rotamer)

! Fill up rotamer_record array
call read_rotamer_library_record(n_rotamer)

end subroutine read_rotamer_library
!-------------------------------------------------------------------------------
subroutine read_rotamer_library_index(n_rotamer)
!-------------------------------------------------------------------------------
! TODO: Comment
!-------------------------------------------------------------------------------
integer, intent(out) :: n_rotamer
integer, parameter :: f_unit = 137
integer :: ioerr
integer :: i_aa, i_phi, i_psi
integer :: phi, psi, phi_prev, psi_prev
integer :: n_rotamer_state
real(dp) :: cumm_prob, prob
character(len=4) :: res_name, res_name_prev

allocate(rotamer_index(num_rot_res))

100 format(A4, T6,I4, T11,I4, T34,F8.6)

call log_p('  Reading rotamer library file  '//trim(infile_rotamer_lib), me=me, level=20)

i_aa = 0
n_rotamer = 0
!
res_name_prev = ''
phi_prev = 9999
psi_prev = 9999

! Count the number of available rotamer states
open(f_unit, file=trim(infile_rotamer_lib), iostat=ioerr)

do 
    read(f_unit, 100, iostat=ioerr) res_name, phi, psi, prob
    if (ioerr < 0) exit

    i_phi = phi/10
    i_psi = psi/10

    if (trim(res_name) /= trim(res_name_prev)) then
        res_name_prev = res_name
        phi_prev = phi
        psi_prev = psi
        !
        i_aa = i_aa + 1
        rotamer_index(i_aa)%res_name = res_name
        rotamer_index(i_aa)%start = n_rotamer + 1
        rotamer_index(i_aa)%bbdep_id(1, i_phi, i_psi) = n_rotamer + 1
        !
        cumm_prob = 0.0
        n_rotamer_state = 1
    else if (phi /= phi_prev .or. psi /= psi_prev) then
        phi_prev = phi
        psi_prev = psi
        !
        rotamer_index(i_aa)%bbdep_id(1, i_phi, i_psi) = n_rotamer + 1
        !
        cumm_prob = 0.0
        n_rotamer_state = 1
    end if

    if (n_rotamer_state > max_rotamer_state) cycle
    if (cumm_prob       > max_rotamer_prob ) cycle

    n_rotamer  = n_rotamer + 1
    cumm_prob = cumm_prob + prob
    n_rotamer_state = n_rotamer_state + 1
    !
    rotamer_index(i_aa)%bbdep_id(2, i_phi, i_psi) = n_rotamer
end do

close(f_unit)

do i_aa = 1, num_rot_res-1
    rotamer_index(i_aa)%end = rotamer_index(i_aa+1)%start - 1
end do
rotamer_index(num_rot_res)%end = n_rotamer

do i_aa = 1, num_rot_res
    rotamer_index(i_aa)%n_rottot = rotamer_index(i_aa)%end - rotamer_index(i_aa)%start
    rotamer_index(i_aa)%bbdep_id(:, 18,-18: 18) = rotamer_index(i_aa)%bbdep_id(:,-18,-18: 18)
    rotamer_index(i_aa)%bbdep_id(:,-18: 18, 18) = rotamer_index(i_aa)%bbdep_id(:,-18: 18,-18)
    do i_phi = -18, 18
        do i_psi = -18, 18
            rotamer_index(i_aa)%n_rot(i_phi,i_psi) = &
                rotamer_index(i_aa)%bbdep_id(2,i_phi,i_psi) - rotamer_index(i_aa)%bbdep_id(1,i_phi,i_psi) + 1
        end do
    end do
end do

end subroutine read_rotamer_library_index
!-------------------------------------------------------------------------------
subroutine read_rotamer_library_record(n_rotamer)
!-------------------------------------------------------------------------------
! TODO: Comment
!-------------------------------------------------------------------------------
integer, intent(in) :: n_rotamer
integer, parameter :: f_unit = 137
integer :: ioerr
integer :: i_phi, i_psi, i_rot
integer :: phi, psi, phi_prev, psi_prev
integer :: n_rotamer_state
real(dp) :: cumm_prob, prob, chi(4), sigma(4)
character(len=4) :: res_name, res_name_prev

allocate(rotamer_record(n_rotamer))

110 format(A4, T6,I4, T11,I4, T34,F8.6, &
           T44,F6.1, T52,F6.1, T60,F6.1, T68,F6.1,&
           T78,F6.1, T86,F6.1, T94,F6.1,T102,F6.1)     

write(log_msg,"(A,I10)") '  Total number of rotamers in rot lib file  ', n_rotamer
call log_p(log_msg, me=me, level=20)

! Re-read the rotamer library file to fill up "real" data
i_rot = 0
res_name_prev = ''
phi_prev = 9999
psi_prev = 9999

open(f_unit, file=trim(infile_rotamer_lib), iostat=ioerr)

do 
    read(f_unit, 110, iostat=ioerr) res_name, phi, psi, prob, chi(1:4), sigma(1:4)
    if (ioerr < 0) exit

    if (trim(res_name) /= trim(res_name_prev) .or. phi /= phi_prev .or. psi /= psi_prev) then
        res_name_prev = res_name
        phi_prev = phi
        psi_prev = psi
        !
        cumm_prob = 0.0
        n_rotamer_state = 1
    end if

    if (n_rotamer_state > max_rotamer_state) cycle
    if (cumm_prob       > max_rotamer_prob ) cycle

    cumm_prob = cumm_prob + prob
    n_rotamer_state = n_rotamer_state + 1

    i_phi = phi/10
    i_psi = psi/10
    !
    i_rot = i_rot + 1
    rotamer_record(i_rot)%i_phi  = i_phi
    rotamer_record(i_rot)%i_psi  = i_psi
    rotamer_record(i_rot)%prob  = prob
    rotamer_record(i_rot)%chi   = chi*deg2rad
    rotamer_record(i_rot)%sigma = sigma*deg2rad
end do

close(f_unit)

end subroutine read_rotamer_library_record
!-------------------------------------------------------------------------------
subroutine read_chi_definition()
!-------------------------------------------------------------------------------
! Read definition of chi angles in terms of consisting atoms
!-------------------------------------------------------------------------------
integer :: ioerror, i_aa, i_res, f_unit = 20
character(len=len_fname) :: line, word(4)
character(1) :: chr
character(4) :: res_name
integer :: n_chi, n_word, k

call log_p('  Reading infile_chi_def file  '//trim(infile_chi_def), me=me, level=20)

open(f_unit, file = trim(infile_chi_def))

do
    read (f_unit, "(A120)", iostat =  ioerror) line
    if (ioerror < 0 .or. line(1:3) == 'END') exit

    chr = line(1:1)
    if (chr == " " .or. chr ==  "" .or. chr == "#") cycle
        
    call parse_string(line, n_word, word)

    if (n_word == 1) then
        if (line(1:4) == 'DONE') then
            rotamer_index(i_res)%n_chi = n_chi
        else
            res_name = line(1:4)
            do i_aa = 1, num_rot_res
                if (rotamer_index(i_aa)%res_name == res_name) then
                    i_res = i_aa
                    exit
                end if
            end do
            n_chi = 0
        end if
    else if (n_word == 4) then
        n_chi = n_chi + 1
        do k = 1, 4
            rotamer_index(i_res)%chi_atoms(k,n_chi) = word(k)(1:4)
        end do
    else
        write(log_msg,"(A,A)") 'Error: in reading infile_chi_def file  ', trim(infile_chi_def)
        call terminate_with_error(log_msg)
    end if
end do
  
close(f_unit)

end subroutine read_chi_definition
!-------------------------------------------------------------------------------
subroutine setup_ref_res_rotamer_index()
!-------------------------------------------------------------------------------
! TODO: Comment
!-------------------------------------------------------------------------------
integer :: i_ref, i_res, i_atm, i_rot_res
logical :: N_found, CA_found, C_found, rot_res_found
character(len=4) :: res_name, atom_name

do i_ref = 1, num_ref_res
    res_name = ref_res(i_ref)%res_name

    if (len(trim(res_name)) == 3) then
        N_found = .false.
        CA_found = .false.
        C_found = .false.
        do i_atm = 1, ref_res(i_ref)%n_atm
            atom_name = ref_res(i_ref)%atom_name(i_atm)
            if (atom_name == 'N') then
                N_found = .true.
            else if (atom_name == 'CA') then
                CA_found = .true.
            else if (atom_name == 'C') then
                C_found = .true.
            end if
        end do

        if (N_found .and. CA_found .and. C_found) then
            ! now we made sure that this ref_res is non-terminal amino acid
            ! find corresponding rotamer
            rot_res_found = .false.
            call convert_to_stdres(res_name)
            do i_res = 1, num_rot_res
                if (res_name == rotamer_index(i_res)%res_name) then
                    i_rot_res = i_res
                    rot_res_found = .true.
                    exit
                end if
            end do
            if (rot_res_found) then
                ref_res(i_ref)%i_rot_res = i_rot_res
            else
                if (res_name /= 'ALA' .and. res_name /= 'GLY') then
                    write(log_msg,"(A,A5)") 'Error. Residue in rot lib for ref_res has not been found.', res_name
                    call terminate_with_error(log_msg)
                end if
                ref_res(i_ref)%i_rot_res = 0
            end if
        else
            ref_res(i_ref)%i_rot_res = 0
        end if
    else
        ref_res(i_ref)%i_rot_res = 0
    end if
end do

end subroutine setup_ref_res_rotamer_index
!-------------------------------------------------------------------------------
subroutine setup_ref_res_chi_angle()
!-------------------------------------------------------------------------------
! Fill in arrays in 'ref_res' by information from rotamer library
!-------------------------------------------------------------------------------
integer :: i_ref, i_rot_res
integer :: i_chi, n_chi, i_ang, i_chi_ang
character(len=4) :: atom4, chi_atom4
logical :: chi_atom_found

! chi angle indexing
do i_ref = 1, num_ref_res
    i_rot_res = ref_res(i_ref)%i_rot_res
    if (i_rot_res <= 0) cycle
    n_chi = rotamer_index(i_rot_res)%n_chi
    ref_res(i_ref)%n_chi = n_chi

    do i_chi = 1, n_chi
        chi_atom4 = rotamer_index(i_rot_res)%chi_atoms(4,i_chi)
        chi_atom_found = .false.

        ! t_ang index for each chi (check only the last atom in the t_ang)
        do i_ang = 1, ref_res(i_ref)%n_atm
            atom4 = ref_res(i_ref)%atom_name(ref_res(i_ref)%atm_in_t_ang(4,i_ang))
            if (chi_atom4 == atom4) then
                ! i_chi_ang is the index of the first atom of the chi angle
                i_chi_ang = i_ang
                chi_atom_found = .true.
                exit
            end if
        end do

        if (.not. chi_atom_found) then
            write(log_msg,"(A,I4,A5,A5)") 'Error. chi_atom not found', i_ref, &
                  ref_res(i_ref)%res_name, chi_atom4
            call terminate_with_error(log_msg)
        end if

        ref_res(i_ref)%t_ang_for_chi(i_chi) = i_chi_ang
    end do
end do

end subroutine setup_ref_res_chi_angle
!-------------------------------------------------------------------------------
subroutine setup_ref_res_chi_dep()
!-------------------------------------------------------------------------------
logical :: found
integer :: i_ref, i_chi, i_atm
integer :: atm_prv, b_no, n_chi, chi_vertex(4)
type(ref_res_type) :: ref

do i_ref = 1, num_ref_res
    ref = ref_res(i_ref)
    ref_res(i_ref)%dep_on_chi(:,:) = .false.
    n_chi = ref%n_chi

    if (n_chi == 0) cycle

    ! Upper vertex point
    chi_vertex(:) = 0
    do i_chi = 1, n_chi
        chi_vertex(i_chi) = ref%atm_in_t_ang(3,ref%t_ang_for_chi(i_chi))
    end do

    do i_atm = 1, ref%n_atm
        if (ref%atom_name(i_atm) == 'N   ' .or. ref%atom_name(i_atm) == 'CA  ' .or. &
            ref%atom_name(i_atm) == 'C   ' .or. ref%atom_name(i_atm) == 'O   ' .or. &
            ref%atom_name(i_atm) == 'CB  ' .or. ref%atom_name(i_atm) == 'H   ') cycle

        ! Backtrace from the sidechain end using bond connection
        found = .false.
        atm_prv = i_atm
        do 
            b_no = atm_prv
            if (b_no < 1) cycle

            if (ref%atom_name(ref%atm_in_bnd(1,b_no)) == 'N   ') exit

            do i_chi = 1, n_chi
                if (ref%atm_in_bnd(1,b_no) == chi_vertex(i_chi)) then
                    ref_res(i_ref)%dep_on_chi(1:i_chi,i_atm) = .true.
                    found = .true.
                    exit
                end if
            end do
           
            atm_prv = ref%atm_in_bnd(1,atm_prv)
            if (found) exit
        end do
    end do
end do

end subroutine setup_ref_res_chi_dep
!-------------------------------------------------------------------------------
subroutine get_rotamer_index(residue, rot_index_start_bb, rot_index_end_bb)
!-------------------------------------------------------------------------------
! This returns rotamer index (rot_index_start_bb & rot_index_end_bb)
!  based on current backbone phi/psi angle.
! The residue of interest together with the following residue should be given as 
!  input in order to calculate the psi angle.
!-------------------------------------------------------------------------------
type(residue_type), intent(in) :: residue(2)
integer, intent(out) :: rot_index_start_bb, rot_index_end_bb
integer :: i_rot_res, i_ref_res, i_ref_next_res, i_phi, i_psi
real(dp) :: phi, psi

i_ref_res = residue(1)%res_type
i_ref_next_res = residue(2)%res_type
i_rot_res = ref_res(i_ref_res)%i_rot_res

rot_index_start_bb = 0
rot_index_end_bb = 0
if (i_rot_res == 0) return
  
phi = bound_ang(residue(1)%t_ang(ref_res(i_ref_res)%iphi))
psi = bound_ang(residue(2)%t_ang(ref_res(i_ref_next_res)%ipsi))
  
i_phi = anint(rad2deg*phi/10.0d0)
i_psi = anint(rad2deg*psi/10.0d0)

rot_index_start_bb = rotamer_index(i_rot_res)%bbdep_id(1,i_phi,i_psi)
rot_index_end_bb   = rotamer_index(i_rot_res)%bbdep_id(2,i_phi,i_psi)

end subroutine get_rotamer_index
!-------------------------------------------------------------------------------
subroutine identify_rotamer_state(protein, rotamer_idx, rotamer_idx_full)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer, intent(out) :: rotamer_idx(protein%n_res), rotamer_idx_full(protein%n_res)
integer :: i_res

rotamer_idx(:) = 0
rotamer_idx_full(:) = 0

do i_res = 1, protein%n_res-1
    call identify_rotamer_state_single(protein, i_res, rotamer_idx(i_res), &
                                       rotamer_idx_full(i_res))
end do

end subroutine identify_rotamer_state
!-------------------------------------------------------------------------------
subroutine identify_rotamer_state_single(protein, i_res, rotamer_idx, &
                                         rotamer_idx_full)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer, intent(in) :: i_res
integer, intent(out) :: rotamer_idx, rotamer_idx_full
integer :: i_start_rot, i_end_rot, i_ref_res, i_rot, i_chi
integer :: chi_idx(4)
real(dp) :: min_ang_diff, ang_diff, tmp_ang_diff
real(dp) :: chi_ang, rotamer_ang

rotamer_idx = 0
rotamer_idx_full = 0

call get_rotamer_index(protein%residue(i_res:i_res+1), i_start_rot, i_end_rot)
if (i_start_rot == 0) return

i_ref_res = protein%residue(i_res)%res_type
chi_idx(1:4) = ref_res(i_ref_res)%t_ang_for_chi(1:4)
min_ang_diff = 999999999.9d0

do i_rot = i_start_rot, i_end_rot
    ! calculate angle difference
    ang_diff = 0.0d0

    do i_chi = 1, 4
        if (chi_idx(i_chi) == 0) exit ! No chi angle definition
       
        chi_ang     = rad2deg * protein%residue(i_res)%t_ang(chi_idx(i_chi))
        rotamer_ang = rad2deg * rotamer_record(i_rot)%chi(i_chi)
       
        tmp_ang_diff = abs(180.0d0 - mod(abs(chi_ang-rotamer_ang)+180.0d0, 360.0d0))
        ang_diff = ang_diff + tmp_ang_diff
    end do       

    ! Minimimum angle difference is the rotamer state
    if (ang_diff < min_ang_diff) then
        min_ang_diff     = ang_diff
        rotamer_idx      = i_rot - i_start_rot + 1
        rotamer_idx_full = i_rot
    end if
end do

end subroutine identify_rotamer_state_single
!-------------------------------------------------------------------------------
subroutine place_rotamer(protein, i_res, rotamer_idx) 
!-------------------------------------------------------------------------------
! This subroutine places the sidechains based on the rotamer number.
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
integer, intent(in) :: i_res
integer, intent(in) :: rotamer_idx

if (rotamer_idx <= 0) return

! change internal coordinates
call rotate_chi_angles(protein, i_res, rotamer_record(rotamer_idx)%chi(1:4))

! place cartesian coordinate
call internal2cartesian_sc(protein%residue(i_res))
!call internal2cartesian(i_res, i_res, protein%residue(1:protein%n_res))

end subroutine place_rotamer
!-------------------------------------------------------------------------------
subroutine perturb_rotamer_state(protein, i_res, select_mode)
!-------------------------------------------------------------------------------
! Read protein and randomly perturb sidechain at 'i_res' into another rotamer
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
integer, intent(in) :: i_res
character(15), intent(in) :: select_mode
integer :: i_rot_start, i_rot_end, i_rot_selected

call get_rotamer_index(protein%residue(i_res:i_res+1), i_rot_start, i_rot_end)
if (i_rot_start == 0) return

call select_rotamer(i_rot_start, i_rot_end, select_mode, i_rot_selected)
if (i_rot_selected == 0) return

! reflect the change in chi angles to geometry
call rotate_chi_angles(protein, i_res, rotamer_record(i_rot_selected)%chi(1:4))

end subroutine perturb_rotamer_state
!-------------------------------------------------------------------------------
subroutine select_rotamer(i_rot_start, i_rot_end, select_mode, i_rot_selected)
!-------------------------------------------------------------------------------
! Select backbone-dependent rotamer to return chi(4) and rotamer index:
! i_rot_selected  => full rotamer index
!-------------------------------------------------------------------------------
integer, intent(in) :: i_rot_start, i_rot_end
integer, intent(out) :: i_rot_selected
character(len=15), intent(in) :: select_mode
integer :: num_rot, i_rot
real(dp) :: cumm_prob, prev_prob, ran_prob

ran_prob = random()

if (select_mode == "random") then
    num_rot = i_rot_end - i_rot_start + 1
    i_rot_selected = i_rot_start + int(ran_prob*num_rot)
     
else if (select_mode == "probablistic") then
    cumm_prob = 0.0d0
    i_rot_selected = 0

    do i_rot = i_rot_start, i_rot_end
        prev_prob = cumm_prob
        cumm_prob = prev_prob + rotamer_record(i_rot)%prob
        if (ran_prob >= prev_prob .and. ran_prob < cumm_prob) then
            i_rot_selected = i_rot
            exit
        end if
    end do
     
    if (i_rot_selected == 0) then
        call log_p('Warning. No selected rotamer.',me=me, level=20)
        return
    end if

else
    write(log_msg,"(A,A)") 'Unknown select_mode in select_rotamer: ', select_mode
    call terminate_with_error(log_msg)

end if

end subroutine select_rotamer
!-------------------------------------------------------------------------------
subroutine rotate_chi_angles(protein, i_res, chi, preserve_dep_in)
!-------------------------------------------------------------------------------
! Reconstruct sidechain of the residue by reflecting the change in chi angles
! Replace depending angles using "ideal topology"
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
integer, intent(in) :: i_res
real(dp), intent(in) :: chi(4)
logical, intent(in), optional :: preserve_dep_in
type(ref_res_type) :: ref
integer :: i_chi, i_ang
logical :: preserve_dep

preserve_dep = .false.
if (present(preserve_dep_in)) preserve_dep = .true.

ref = ref_res(protein%residue(i_res)%res_type)
do i_chi = 1, ref%n_chi
    i_ang = ref%t_ang_for_chi(i_chi)
    ! Actually, the second argument for rotate_torsion_angle should be 
    !  the previous residue. However, it is only for -O atom, and it CANNOT
    !  be occured while rotating chi angles. So, to prevent indexing error
    !  it is replaced with the current residue
    !call rotate_torsion_angle(protein%residue(i_res), protein%residue(i_res-1), &
    call rotate_torsion_angle(protein%residue(i_res), protein%residue(i_res), &
        i_ang, chi(i_chi), preserve_dep_in=preserve_dep)
end do

end subroutine rotate_chi_angles
!-------------------------------------------------------------------------------
function is_symmetric_chi(i_ref, i_chi)
!-------------------------------------------------------------------------------
! Check if i-th chi in given amino acid type is symmetric, 
!  following Dunbrack BBdep library definition
!-------------------------------------------------------------------------------
logical :: is_symmetric_chi
integer, intent(in) :: i_ref, i_chi
  
if ((trim(ref_res(i_ref)%res_name) == 'ASP' .and. i_chi == 2) .or. &
    (trim(ref_res(i_ref)%res_name) == 'PHE' .and. i_chi == 2) .or. &
    (trim(ref_res(i_ref)%res_name) == 'TYR' .and. i_chi == 2) .or. &
    (trim(ref_res(i_ref)%res_name) == 'PTR' .and. i_chi == 2) .or. &
    (trim(ref_res(i_ref)%res_name) == 'GLU' .and. i_chi == 3)) then
    is_symmetric_chi = .true.
else
    is_symmetric_chi = .false.
end if

end function is_symmetric_chi
!-------------------------------------------------------------------------------
END MODULE ROTAMER
!-------------------------------------------------------------------------------
