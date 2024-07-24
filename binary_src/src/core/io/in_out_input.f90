!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/io/in_out_input.f90
!
! Description:
!   This module contains subroutines for read input files, set default
!   parameters/libraries, and check consistency of options in input files.
!-------------------------------------------------------------------------------
MODULE IN_OUT_INPUT
!-------------------------------------------------------------------------------
use globals
use in_out_vars
use logger,  only: log_unit, open_logfile, log_p, terminate_with_error
use string,  only: parse_string
use ran,     only: seed, seed_given
!
use ramachandran, only: use_rama, rama_mode, infile_psipred, infile_psrm
use rotamer, only: max_rotamer_state, max_rotamer_prob
use fragment, only: use_fragment, max_frag_lib, max_frag_len, &
                    n_frag_lib,  frag_lib
use symmetry, only: n_symm, symm_resrange

implicit none
private

integer :: max_usc_res
logical, allocatable :: is_usc_prev(:)

public :: read_global_input
public :: read_pdblist
public :: read_initial_qa
public :: set_default_global
public :: set_default_library
public :: reinitialize_usc

CONTAINS
!-------------------------------------------------------------------------------
subroutine read_global_input(file_name)
!-------------------------------------------------------------------------------
! Read global options from input file
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: file_name
integer :: f_unit = 20, ioerror, num_word
character(len=len_fname) :: line, word(25), keyword
integer :: i_topo, i_topoe, i_range, i_word
integer :: i_res

! default
call set_default_global()

! initialize
i_topo = 0
i_topoe = 0
i_range = 0

n_frag_lib = 0

open(f_unit, file = trim(file_name))

infile_user_input = file_name

do 
    read(f_unit, "(A250)", iostat = ioerror) line
    if (ioerror < 0) exit
    if (line(1:1) == '!' .or. line(1:1) == '#') cycle
    !
    call parse_string(line, num_word, word)
    keyword = word(1)

    ! Mode selection
    if (keyword == 'molecule_type') then
        mol_type = word(2)
    else if (keyword == 'task_mode') then
        task_mode = word(2)
    else if (keyword == 'top_type') then
        top_type = word(2)
        if (trim(top_type) == 'polarh' .or. trim(top_type) == 'eef1' .or. &
            trim(top_type) == 'coarse' .or. trim(top_type) == 'allh_ch22' .or. &
            trim(top_type) == 'cbeta') then
            force_field_type = 'CHARMM'
        end if
    else if (keyword == 'read_het') then
        if (word(2) == 'yes') then
            read_het = .true.
        else
            read_het = .false.
        end if

    ! File I/O
    else if (keyword == 'data_directory') then
        write(data_dir,'(A,A)') trim(word(2)),'/'
    else if (keyword == 'infile_pdb') then
        infile_pdb = word(2)
        multiple_task = .false.
    else if (keyword == 'infile_pre_ML') then
        infile_pre_ML = word(2)
    else if (keyword == 'infile_pdblist') then
        infile_pdblist = word(2)
        multiple_task = .true.
    else if (keyword == 'read_multiple_models') then
        if (word(2) == 'yes') then
            multiple_models = .true.
        else
            multiple_models = .false.
        end if
    else if (keyword == 'infile_ligand') then 
        read_het = .true.
        infile_ligand(1) = word(2)
        if (num_word == 2) then
            infile_ligand(2) = 'LIG'
        else
            infile_ligand(2) = word(3)
        end if
        n_mol2_top = n_mol2_top + 1
        infile_mol2_topo(1:2,n_mol2_top) = infile_ligand(1:2)
    else if (keyword == 'infile_mol2_topo') then
        n_mol2_top = n_mol2_top + 1
        infile_mol2_topo(1,n_mol2_top) = word(2)
        infile_mol2_topo(2,n_mol2_top) = word(3)
    else if (keyword == 'outfile_prefix') then
        outfile_prefix = word(2)
    else if (keyword == 'single_output') then
        if (word(2) == 'yes') then
            single_output = .true.
        else
            single_output = .false.
        end if
    else if (keyword == 'logfile') then
        logfile = word(2)
        log_unit = 98
        call open_logfile(logfile, log_unit)
    else if (keyword == 'print_level') then
        read(word(2), "(I10)") print_level_global

    else if (keyword == 'force_field_type') then
        force_field_type = word(2)
    else if (keyword == 'infile_topology') then
        i_topo = i_topo + 1
        infile_topo(i_topo) = trim(data_dir)//word(2)
        topo_file_mol_type(i_topo) = word(3)
        topo_file_ter_type(i_topo) = word(4)
    else if (keyword == 'infile_parameter') then
        infile_parameter = trim(data_dir)//word(2)
    else if (keyword == 'infile_topo_eng') then
        i_topoe = i_topoe + 1
        infile_topo_eng(i_topoe) = trim(data_dir)//word(2)
        topo_eng_file_mol_type(i_topoe) = word(3)
    else if (keyword == 'infile_topology_het') then
        infile_topo(4) = trim(data_dir)//word(2)
        topo_file_mol_type(4) = 'hetmol'
        topo_file_ter_type(4) = ''
    else if (keyword == 'infile_topo_eng_het') then
        infile_topo_eng(2) = trim(data_dir)//word(2)
        topo_eng_file_mol_type(2) = 'hetmol'
    else if (keyword == 'infile_rotamer_lib') then
        infile_rotamer_lib = trim(data_dir)//word(2)
    else if (keyword == 'infile_chi_def') then
        infile_chi_def = trim(data_dir)//word(2)
    else if (keyword == 'max_rotamer_state') then
        read(word(2), '(I4)') max_rotamer_state
        if (max_rotamer_state > 81 .or. max_rotamer_state <= 0) then
            max_rotamer_state = 81
            call log_p("  WARNING: maximum number of rotamer states can be [1, 81].")
            call log_p("           max_rotamer_state is set to 81 (default).")
        end if
    else if (keyword == 'max_rotamer_prob') then
        read(word(2), '(F12.5)') max_rotamer_prob
        if (max_rotamer_prob < 0.0 .or. max_rotamer_prob > 1.0) then
            max_rotamer_prob = 1.1
            call log_p("  WARNING: rotamer cummulative probability cutoff can be [0.0, 1.0].")
            call log_p("           max_rotamer_prob is set to 1.0 (default).")
        end if
    else if (keyword == 'allow_prot_state') then
        if (word(2) == 'yes') then
            allow_prot_state = .true.
        else if (word(2) == 'no') then
            allow_prot_state = .false.
        end if
    else if (keyword == 'infile_frag' .or. keyword == 'frag_lib_file') then
        use_fragment = .true.
        !
        n_frag_lib = n_frag_lib + 1
        if (n_frag_lib > max_frag_lib) then
            call terminate_with_error("ERROR: n_frag_lib > max_frag_lib.")
        end if
        !
        frag_lib(n_frag_lib)%frag_lib_file = word(2)
        if (num_word > 2) then
            read(word(3),"(I2)") frag_lib(n_frag_lib)%n_len
        else
            frag_lib(n_frag_lib)%n_len = 9
        end if
        !
        if (n_frag_lib == 1 .or. frag_lib(n_frag_lib)%n_len > max_frag_len) then
             max_frag_len = frag_lib(n_frag_lib)%n_len
        end if
        frag_lib(n_frag_lib)%pos_shift = int(frag_lib(n_frag_lib)%n_len/2)
    else if (keyword == 'rama_mode') then
        use_rama = .true.
        rama_mode = word(2)
    else if (keyword == 'infile_psipred') then
        use_rama = .true.
        infile_psipred = word(2)
    else if (keyword == 'infile_psrm') then
        use_rama = .true.
        infile_psrm = word(2)

    ! Fixing atoms
    else if (keyword == 'fix_type') then
        if (word(2) == 'all') then
            fix_type = FIX_ALL
        else if (word(2) == 'backbone') then
            fix_type = FIX_BB
        else
            fix_type = FIX_NONE
        end if
    else if (keyword == 'fix_atom_file') then
        fix_atom_file = word(2)
    else if (keyword == 'ULR') then
        n_ulr = n_ulr + 1
        ULR(n_ulr)%type = word(2)
        ULR(n_ulr)%subtype = word(3)
        read(word(4),"(F15.5)") ULR(n_ulr)%importance
        if (word(2) == 'L') then ! loop modeling type
            read(word(5),"(I10)") ULR(n_ulr)%resrange(1)
            read(word(6),"(I10)") ULR(n_ulr)%resrange(2)
        else if (word(2) == 'S') then ! segment modeling type
            read(word(5),"(I10)") ULR(n_ulr)%resrange(1)
            read(word(6),"(I10)") ULR(n_ulr)%segrange(1)
            read(word(7),"(I10)") ULR(n_ulr)%segrange(2)
            read(word(8),"(I10)") ULR(n_ulr)%resrange(2)
        end if
    else if (keyword == 'USC') then
        do i_word = 2, num_word
            read(word(i_word),"(I4)") i_res
            if (i_res > max_usc_res) then
                allocate(is_usc_prev(i_res))
                is_usc_prev(:) = .false.
                is_usc_prev(1:max_usc_res) = is_usc(1:max_usc_res)
                call move_alloc(is_usc_prev, is_usc)
                max_usc_res = i_res
            end if
            is_usc(i_res) = .true.
            n_usc = n_usc + 1
        end do
    else if (keyword == 'infile_qa') then
        infile_qa = word(2)
    else if (keyword == 'infile_ppDock_opr') then
        infile_ppDock_opr = word(2)
    else if (keyword == 'infile_native') then
        evaluate_TMscore = .true.
        infile_native = word(2)
    
    ! Define receptor & ligand protein
    else if (keyword == 'n_res_receptor') then
        read(word(2), "(I10)") tn%recres

    ! Etc.
    else if (keyword == 'random') then
        read(word(2),"(I10)") seed
        seed_given = .true.
    ! Symmetry
    else if (keyword == 'symm_unit') then
        n_symm = n_symm + 1
        read(word(2),"(I10)") symm_resrange(1,n_symm)
        read(word(3),"(I10)") symm_resrange(2,n_symm)
    else if (keyword == 'allow_broken_bond') then
        if (word(2) == 'yes') then
            allow_broken_bond = .true.
        else
            allow_broken_bond = .false.
        end if
    else if (keyword == 'check_disulfide_pdb') then
        if (word(2) == 'yes') then
            check_disulfide_pdb = .true.
        else
            check_disulfide_pdb = .false.
        end if
    end if
end do

close(f_unit)

call set_default_library()
call set_usc_information()
call set_default_ML()
call check_input_consistency()

end subroutine read_global_input
!-------------------------------------------------------------------------------
subroutine set_default_global()
!-------------------------------------------------------------------------------
! Set default parameters before reading input
!-------------------------------------------------------------------------------
! 1. General Option
infile_topo(:) = ''
data_dir = ''
mol_type = 'protein'
top_type = ''
force_field_type = 'AMBER'
print_level_global = 0

! 2. I/O filename
outfile_prefix = 'out'
infile_pdb = ''
infile_pdblist = ''
infile_ligand(:) = ''
infile_mol2_topo(:,:) = ''
multiple_task = .false.
multiple_models = .true.
log_unit = 6 ! print at screen
logfile = ''
single_output = .true.
infile_qa = ''
infile_native = ''
evaluate_TMscore = .false.
use_fragment = .false.

! 3. Library
infile_rotamer_lib = ''
infile_chi_def = ''
max_rotamer_state = 81
max_rotamer_prob  = 1.001
!
use_rama = .false.
rama_mode = 'neighbor'
infile_psipred = ''
infile_psrm = ''

! 6. Hetmol
read_het = .false.
n_mol2_top = 0

! 7. Atom fix type
allocate(is_usc(max_res))
fix_type = FIX_NONE    ! all bb/sc free
n_ulr = 0              ! no ulrs
n_usc = 0              ! no uscs
is_usc(:) = .false.    ! no uscs
fix_atom_file = ''     ! no user-defined fix atom file
are_defined = .false.  ! no selected residues
max_usc_res = max_res  ! initialize max_usc_res

! 8. Etc.
seed_given = .false.
allow_prot_state = .false.
n_symm = 0
tn%recres = 0
allow_broken_bond = .false.
check_disulfide_pdb = .false.

end subroutine set_default_global
!-------------------------------------------------------------------------------
subroutine set_default_library()
!-------------------------------------------------------------------------------
! Set library file according to given topology type
!-------------------------------------------------------------------------------
! Set default parameter (based on topology type selected)
if (trim(top_type) == 'heavy') then
    infile_topo(1) = trim(data_dir)//'topology_amber94mod.in'
    infile_topo(2) = trim(data_dir)//'topology_amber94_ntmod.in'
    infile_topo(3) = trim(data_dir)//'topology_amber94_ctmod.in'
    topo_file_mol_type(1:3) = 'protein'
    topo_file_ter_type(2) = 'N'
    topo_file_ter_type(3) = 'C'
    infile_topo_eng(1) = trim(data_dir)//'topology_eng_amber94mod.in'
    topo_eng_file_mol_type(1) = 'protein'
    infile_parameter = trim(data_dir)//'parameters_amber96mod.in'

else if (trim(top_type) == 'polarh') then
    infile_topo(1) = trim(data_dir)//'topology_charmm22ph.in'
    infile_topo(2) = trim(data_dir)//'topology_charmm22ph_nt.in'
    infile_topo(3) = trim(data_dir)//'topology_charmm22ph_ct.in'
    topo_file_mol_type(1:3) = 'protein'
    topo_file_ter_type(2) = 'N'
    topo_file_ter_type(3) = 'C'
    infile_topo_eng(1) = trim(data_dir)//'topology_eng_charmm22ph.in'
    topo_eng_file_mol_type(1) = 'protein'
    infile_parameter = trim(data_dir)//'parameters_charmm22ph.in'
    infile_parameter_ligand = trim(data_dir)//'parameters_charmm22cgenff.in'

else if (trim(top_type) == 'cbeta') then
    infile_topo(1) = trim(data_dir)//'topology_charmm22cb.in'
    infile_topo(2) = trim(data_dir)//'topology_charmm22cb_nt.in'
    infile_topo(3) = trim(data_dir)//'topology_charmm22cb_ct.in'
    topo_file_mol_type(1:3) = 'protein'
    topo_file_ter_type(2) = 'N'
    topo_file_ter_type(3) = 'C'
    infile_topo_eng(1) = trim(data_dir)//'topology_eng_charmm22cb.in'
    topo_eng_file_mol_type(1) = 'protein'
    infile_parameter = trim(data_dir)//'parameters_charmm22cb.in'
    
else if (trim(top_type) == 'eef1') then
    infile_topo(1) = trim(data_dir)//'topology_eef1_22.in'
    infile_topo(2) = trim(data_dir)//'topology_eef1_22_nt.in'
    infile_topo(3) = trim(data_dir)//'topology_eef1_22_ct.in'
    topo_file_mol_type(1:3) = 'protein'
    topo_file_ter_type(2) = 'N'
    topo_file_ter_type(3) = 'C'
    infile_topo_eng(1) = trim(data_dir)//'topology_eng_charmm22ph.in'
    topo_eng_file_mol_type(1) = 'protein'
    infile_parameter = trim(data_dir)//'parameters_charmm22ph.in'

else if (trim(top_type) == 'allh') then
    infile_topo(1) = trim(data_dir)//'topology_amber94.in'
    infile_topo(2) = trim(data_dir)//'topology_amber94_nt.in'
    infile_topo(3) = trim(data_dir)//'topology_amber94_ct.in'
    topo_file_mol_type(1:3) = 'protein'
    topo_file_ter_type(2) = 'N'
    topo_file_ter_type(3) = 'C'
    infile_topo_eng(1) = trim(data_dir)//'topology_eng_amber94.in'
    topo_eng_file_mol_type(1) = 'protein'
    infile_parameter = trim(data_dir)//'parameters_amber96.in'
    allow_prot_state = .true.

else if (trim(top_type) == 'allh_ch22') then
    infile_topo(1) = trim(data_dir)//'topology_charmm22.in'
    infile_topo(2) = trim(data_dir)//'topology_charmm22_nt.in'
    infile_topo(3) = trim(data_dir)//'topology_charmm22_ct.in'
    topo_file_mol_type(1:3) = 'protein'
    topo_file_ter_type(2) = 'N'
    topo_file_ter_type(3) = 'C'
    infile_topo_eng(1) = trim(data_dir)//'topology_eng_charmm22.in'
    topo_eng_file_mol_type(1) = 'protein'
    infile_parameter = trim(data_dir)//'parameters_charmm22.in'
    allow_prot_state = .true.
    infile_parameter_ligand = trim(data_dir)//'parameters_charmm22cgenff.in'

else if (trim(top_type) == 'coarse') then
    infile_topo(1) = trim(data_dir)//'topology_martini.in'
    infile_topo(2) = trim(data_dir)//'topology_martini_nt.in'
    infile_topo(3) = trim(data_dir)//'topology_martini_ct.in'
    topo_file_mol_type(1:3) = 'protein'
    topo_file_ter_type(2) = 'N'
    topo_file_ter_type(3) = 'C'
    infile_topo_eng(1) = trim(data_dir)//'topology_eng_martini.in'
    topo_eng_file_mol_type(1) = 'protein'
    infile_parameter = trim(data_dir)//'parameters_martini.in'

end if

n_topo_file = 3
n_topo_eng_file = 1

! Top/Parm for hetero molecules
if (read_het) then
    if (trim(force_field_type) == 'AMBER') then
        infile_topo(4) = trim(data_dir)//'topology_amber_het.in'
        infile_topo_eng(2) = trim(data_dir)//'topology_amber_eng_het.in'
    else if (trim(force_field_type) == 'CHARMM') then
        infile_topo(4) = trim(data_dir)//'topology_charmm22_het.in'
        infile_topo_eng(2) = trim(data_dir)//'topology_charmm22_eng_het.in'
    end if

    if (infile_topo(4) /= '') then
        topo_file_mol_type(4) = 'hetmol'
        topo_eng_file_mol_type(2) = 'hetmol'
        topo_file_ter_type(4) = ''
        n_topo_file = n_topo_file + 1
        n_topo_eng_file = n_topo_eng_file + 1
    end if
end if

if (trim(infile_rotamer_lib) == '') then
    infile_rotamer_lib = trim(data_dir)//'bbdep02.May.sortlib'
    !infile_rotamer_lib = trim(data_dir)//'bbdep10.lib'
end if
if (trim(infile_chi_def) == '') then
    infile_chi_def = trim(data_dir)//'def_chi.in'
end if

end subroutine set_default_library
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine set_default_ML()
!-------------------------------------------------------------------------------
! Set python scripts and models
!-------------------------------------------------------------------------------
infile_ML(1) = trim(data_dir)//'ML_inference.py'
infile_ML(2) = trim(data_dir)//'ML_rerank.py'
infile_ML(3) = trim(data_dir)//'sampling_model.pt'
infile_ML(4) = trim(data_dir)//'rerank_model.pt'

end subroutine set_default_ML
!-------------------------------------------------------------------------------
subroutine set_usc_information()
!-------------------------------------------------------------------------------
integer :: i_ulr, i_res

if (trim(fix_atom_file) /= '') then
    if (n_usc > 0) then
        call terminate_with_error('ERROR: USC are both in input and fix_atom_file.')
    else
        call read_fix_atom_file()
    end if
end if

! set ULR residues as USC
do i_ulr = 1, n_ulr
    if (ULR(i_ulr)%resrange(2) > max_usc_res) then
        allocate(is_usc_prev(ULR(i_ulr)%resrange(2)))
        is_usc_prev(:) = .false.
        is_usc_prev(1:max_usc_res) = is_usc(1:max_usc_res)
        call move_alloc(is_usc_prev, is_usc)
        max_usc_res = ULR(i_ulr)%resrange(2)
    end if
    do i_res = ULR(i_ulr)%resrange(1), ULR(i_ulr)%resrange(2)
        if (is_usc(i_res)) cycle
        is_usc(i_res) = .true.
        n_usc = n_usc + 1
    end do
end do

end subroutine set_usc_information
!-------------------------------------------------------------------------------
subroutine reinitialize_usc(n_res)
!-------------------------------------------------------------------------------
integer, intent(in) :: n_res

if (max_usc_res == n_res) return

allocate(is_usc_prev(n_res))
is_usc_prev(:) = .false.
is_usc_prev(1:min(n_res, max_usc_res)) = is_usc(1:min(n_res,max_usc_res))
call move_alloc(is_usc_prev, is_usc)

end subroutine reinitialize_usc
!-------------------------------------------------------------------------------
subroutine check_input_consistency()
!-------------------------------------------------------------------------------
! Check input consistency about global parameters in GALAXY
! This subroutine provides user the reason for termination 
!-------------------------------------------------------------------------------
logical :: error_found

error_found = .false.

call log_p('Check global input consistency...', me=me, level=10)

if (trim(data_dir) == '') then
    call log_p('  ERROR: {data_directory} not defined.',me=me)
    error_found = .true.
end if

if (trim(top_type) /= 'heavy' .and. trim(top_type) /= 'allh' .and. &
    trim(top_type) /= 'polarh' .and. trim(top_type) /= 'eef1' .and. &
    trim(top_type) /= 'coarse' .and. trim(top_type) /= 'allh_ch22' .and. &
    trim(top_type) /= 'cbeta') then
    call log_p('  ERROR: top_type should be [heavy/allh/allh_ch22/polarh/eef1/coarse/cbeta].',me=me)
    error_found = .true.
end if

if (fix_type == FIX_BB .and. n_ulr > 0) then
    call log_p('  ERROR: When fix_type is set to backbone, only USC can be declared.', me=me)
    error_found = .true.
end if

if (error_found) then
    call terminate_with_error('Terminate with error: Fatal error during reading input.')
end if

call log_p('Done.', me=me, level=10)

end subroutine check_input_consistency
!-------------------------------------------------------------------------------
subroutine read_pdblist(infile_pdblist, fnamelist, nfile)
!-------------------------------------------------------------------------------
! Read multiple number of PDB files
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: infile_pdblist
character(len=len_fname), intent(out) :: fnamelist(max_pdbfile)
integer, intent(out) :: nfile
integer :: f_unit, ioerror, openstat
character(len=len_fname) :: line

f_unit = 30
open(f_unit, file = trim(infile_pdblist), iostat=openstat)
if (openstat > 0) then
    call log_p('Error: Input pdblist file not found.')
    stop
end if

nfile = 0
do
    read(f_unit, "(A120)", iostat = ioerror) line
    if (ioerror < 0) exit

    nfile = nfile + 1
    read(line,'(A120)') fnamelist(nfile)
end do

close(f_unit)

end subroutine read_pdblist
!-------------------------------------------------------------------------------
subroutine read_fix_atom_file()
!-------------------------------------------------------------------------------
! Reading user-assigned atom fixation from a file
!-------------------------------------------------------------------------------
integer :: f_unit, ioerror, num_word, openstat
character(len=len_fname) :: word(2), line
integer :: i_res

! Ignore USC info in input file.
is_usc(:) = .false.
n_usc = 0
write(log_msg,"(A)") " Reading fix_atom_file."
call log_p(log_msg, level=10,me=me)

f_unit = 41
open(f_unit, file = trim(fix_atom_file), iostat=openstat)
if (openstat > 0) then
    call terminate_with_error('Terminate with error: No fix atom file found. Check {fix_atom_file}.')
end if

do
    read(f_unit,"(A120)", iostat = ioerror) line
    if (ioerror < 0) exit
    if (line(1:1) == '!' .or. line(1:1) == '#') cycle
    !
    call parse_string(line, num_word, word)

    if (num_word /= 2) cycle

    if (trim(word(2)) == 'USC') then
        read(word(1),"(I4)") i_res
        if (i_res > max_usc_res) then
            allocate(is_usc_prev(i_res))
            is_usc_prev(:) = .false.
            is_usc_prev(1:max_usc_res) = is_usc(1:max_usc_res)
            call move_alloc(is_usc_prev, is_usc)
            max_usc_res = i_res
        end if
        is_usc(i_res) = .true.
        n_usc = n_usc + 1
    end if
end do

close(f_unit)

end subroutine read_fix_atom_file
!-------------------------------------------------------------------------------
subroutine read_initial_qa(infile_qa)
!-------------------------------------------------------------------------------
character(len=len_fname) :: infile_qa
integer, parameter :: f_unit = 28
integer :: ioerror, num_word
character(len=len_fname) :: line, word(25)
integer :: res_no
real(dp) :: qa_score

allocate(init_qa(tn%residue))
init_qa(:) = 0.0d0

open(f_unit, file = trim(infile_qa), iostat=ioerror)
if (ioerror /= 0) then
    call terminate_with_error("Error: cannot open the QA file.")
end if

do
    read(f_unit, "(A120)", iostat = ioerror) line
    if (ioerror < 0) exit
    if (line(1:1) == '!' .or. line(1:1) == '#') cycle
    !
    call parse_string(line, num_word, word)
    !
    read(word(1), '(I5)') res_no
    read(word(2), '(F12.7)') qa_score
    !
    if (res_no < 1 .or. res_no > tn%residue) then
        call terminate_with_error("ERROR: residue number in the QA file is not exists.")
    end if
    init_qa(res_no) = max(small_real, min(1.0d0, qa_score))
end do

init_qa(:) = init_qa(:) / sum(init_qa)

end subroutine read_initial_qa
!!-------------------------------------------------------------------------------
END MODULE IN_OUT_INPUT
!-------------------------------------------------------------------------------
