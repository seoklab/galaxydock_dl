!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/io/in_out_structure.f90
!
! Description:
!   This module contains subroutines related to read and write pdb file.
!-------------------------------------------------------------------------------
MODULE IN_OUT_STRUCTURE
!-------------------------------------------------------------------------------
use globals
use in_out_vars
use in_out_utils
use allocate_molecule
use logger, only: log_p, terminate_with_error, log_divider
use string, only: parse_string
use symmetry, only: symm_U, symm_T, n_symm, symm_resrange, max_symm, &
                    n_chain_symm_unit, symm_chain, write_symm_info

implicit none
private

! To read multiple models in a PDB file
integer, save :: n_model_in_buffer
integer, save :: model_no_in_buffer(max_pdbfile)
character(len=len_fname), save :: infile_pdb_in_buffer(2)
!
real(dp), allocatable, save :: residue_buffer(:,:,:,:)
real(dp), allocatable, save :: hetmol_buffer(:,:,:,:)
real(dp), allocatable, save :: ligand_buffer(:,:,:,:)
character(len=len_fname), allocatable, save :: remark350_buffer(:,:)

public :: read_pdb
public :: read_sequence_from_pdb
!
public :: write_pdb
public :: write_pdb_model
!
public :: open_write_pdb
public :: close_write_pdb
!
public :: finalize_structure
public :: n_models_in_pdb
public :: read_ppDock_operations

CONTAINS
!-------------------------------------------------------------------------------
subroutine read_sequence_from_pdb(pdb_file, protein)
!-------------------------------------------------------------------------------
! Get sequence from a pdb file
! Adjust residue names (considering terminal residues, protonation states, etc)
! before reading coordinates.
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: pdb_file
type(molecule_type), intent(inout) :: protein
type(molecule_type), save :: protein_buffer
character(len=4) :: res_name, atom_name, het_name
character(len=5) :: res_prev, res_curr, het_prev, het_curr
character(len=1) :: chain_prev, chain_curr
character(len=len_fname) :: line, word(20)
integer :: res_no, f_unit = 10, ioerror, het_no, lig_no, chain_no, openstat, chain_noprv
integer :: link_no, res_curr_no, het_curr_no
integer :: i, i_res, j_res, i_ss, num_word, i_mol2
logical :: is_ligand
logical :: new_res_no_found, HD1_read, HE2_read, new_ss_bond
logical :: HD_read, HE_read
real(dp) :: dr(3), dr2, r_n(3), r_c(3), r_c_prev(3)
real(dp), allocatable :: r_s(:,:), r_s_prev(:,:)
character(len=1) :: chain
integer :: res_prev_no
integer :: cis_res_list(10)
integer :: cis_res_num, i_cis_res
integer :: curr_model, n_model
logical :: read_model
character(len=4) :: ori_res_name

if ((trim(pdb_file) == trim(infile_pdb_in_buffer(1)))) then
    protein = protein_buffer
    return
end if

call log_p('  Reading sequence from pdb file: '//trim(pdb_file), me=me, level=30)

open(f_unit, file = trim(pdb_file), iostat = openstat)

if (openstat > 0) then
    write(log_msg,"(A,A)") 'Terminate with error: No pdbfile found with filename: ', trim(pdb_file)
    call log_p(log_msg,me=me)
    call terminate_with_error('Please check {infile_pdb} or {infile_pdblist} or {csa_first_bank_fname}.')
end if

res_no = 0
link_no = 0
het_no = 0
lig_no = 0
chain_no = 0
res_prev = '-1000' ! something unlikely
chain_prev = '!' ! something unlikely
het_prev = '-1000'
res_name = ''
protein%residue(:)%ter_type = ''
protein%residue(:)%res_added = '' ! character right after the residue number (column 27)
cis_res_num = 0
cis_res_list(:) = 0
symmetric = .false. 
use_remark350 = .false.
n_h_lines = 0 ! initialize number of header lines

n_model = 0
if (multiple_models) then
    read_model = .true.
else
    read_model = .false.
end if
read_model = .true.

allocate(r_s(3,max_res))

do  
    read(f_unit, "(A100)", iostat=ioerror) line
    if (ioerror < 0) exit ! end of file

    if (multiple_models .and. line(1:5) == 'MODEL') then
        call parse_string(line, num_word, word)
        read(word(2), '(I5)') curr_model
        n_model = n_model + 1
        model_no_in_buffer(n_model) = curr_model
        cycle
    else if (line(1:6) == 'ENDMDL') then
        if (read_model) then
            read_model = .false.
            cycle
        else if (.not. multiple_models) then
            exit
        end if
    end if
    if (.not. read_model) cycle

    if (line(1:4) == 'ATOM') then
        read(line(23:26),"(I4)") res_curr_no
        res_curr = line(23:27)
        chain_curr = line(22:22)
        
        if ((chain_prev /= chain_curr) .or. (res_curr /= res_prev)) then !new residue
            res_prev_no = res_no
            res_no = res_no + 1
          
            if ((chain_prev /= chain_curr)) then
                chain_no = chain_no + 1
                chain_prev = chain_curr
            end if
         
            if (res_no > protein%n_res) then
                allocate(r_s_prev(3,protein%n_res+max_res))
                r_s_prev(:,1:protein%n_res) = r_s(:,1:protein%n_res)
                call move_alloc(r_s_prev, r_s)
                !
                call reallocate_residue_in_molecule(protein, protein%n_res+max_res)
            end if
            protein%residue(res_no)%pdb_res_no = res_curr_no
            protein%residue(res_no)%res_added = res_curr(5:5)
            res_prev = res_curr
            chain_prev = chain_curr
            res_name = line(18:20)
            protein%residue(res_no)%pdb_res_name = res_name
            protein%residue(res_no)%res_name = res_name
            protein%residue(res_no)%chain = line(22:22)
            protein%residue(res_no)%i_chain = chain_no
            protein%residue(res_no)%is_broken = .false.
            ! initialize protonation state of HIS
            if (res_name == 'HIS') then
                HD1_read = .false.
                HE2_read = .false.
            else if (res_name == 'ASP') then
                HD_read = .false.
            else if (res_name == 'GLU') then
                HE_read = .false.
            end if
            ! when allowing different protonated states
            ! XX_read will change into true if XX atom exists in the residue
            if (allow_prot_state) then
                if (res_name == 'HID' .or. res_name == 'HIE' .or. res_name == 'HIP') then
                    HD1_read = .false. 
                    HE2_read = .false.
                else if (res_name == 'ASH') then
                    HD_read = .false.
                else if (res_name == 'GLH') then
                    HE_read = .false.
                end if
            end if
        end if !end of residue intialization procedure
      
        ! save S coord for CYS for SS bond check
        if (res_name == 'CYS') then
            atom_name = line(13:16)
            if (atom_name == ' SG ') then 
                read(line(31:), "(3F8.3)") r_s(1:3,res_no)
            end if
        end if
      
        ! save N coord for C-N bond check
        atom_name = line(13:16)
        ! 1HH1 => HH11
        if ((atom_name(1:1) == '1' .or. atom_name(1:1) == '2') .and. atom_name(2:2) == 'H') then
            atom_name = atom_name(2:4)//atom_name(1:1)
        end if
      
        ! setup N and C terminal and find broken residue connection
        if ( (atom_name == ' N  ') .or. &
              ((atom_name == ' CH3') .and. (res_name == 'ACE'))  ) then
            read(line(31:), "(3F8.3)") r_n(1:3)
            if (res_no > 1) then
                r_c_prev(:) = r_c(:) ! assume that N is read earlier than C
                dr(:) = r_n(:) - r_c_prev(:)
                dr2 = dot_product(dr,dr)
               
                if (protein%residue(res_no)%chain /= protein%residue(res_prev_no)%chain) then
                    protein%residue(res_no)%ter_type = 'N'
                    protein%residue(res_prev_no)%ter_type = 'C'
                else if (dr2 > cn2_cut) then
                    write(log_msg, "(A,2(I5,A2,1X,A4),A,F10.3)") '    There seems to be break in chain.', &
                          protein%residue(res_prev_no)%pdb_res_no, &
                          protein%residue(res_prev_no)%chain, &
                          protein%residue(res_prev_no)%pdb_res_name, &
                          protein%residue(res_no)%pdb_res_no, &
                          protein%residue(res_no)%chain, &
                          protein%residue(res_no)%pdb_res_name, &
                          ' C-N bond length is', sqrt(dr2)
                    call log_p(log_msg,me=me,level=20)
                    protein%residue(res_prev_no)%is_broken = .true.
                    protein%residue(res_no)%is_broken = .true.
                end if
            end if
        else if (atom_name == ' C') then
            read(line(31:), "(3F8.3)") r_c(1:3)
        end if
      
        ! check protonation state of HIS
        if (allow_prot_state) then
            if (res_name=='HIS' .or. res_name=='HID' .or. res_name=='HIE' .or. res_name=='HIP') then
                atom_name = line(13:16)
                if (atom_name == ' HD1') HD1_read = .true.
                if (atom_name == ' HE2') HE2_read = .true.
                ! update protonation state of HIS
                if (HD1_read .and. .not.HE2_read) then
                    protein%residue(res_no)%res_name = 'HID'
                else if (.not.HD1_read .and. HE2_read) then
                    protein%residue(res_no)%res_name = 'HIE'
                else if (HD1_read .and. HE2_read) then
                    protein%residue(res_no)%res_name = 'HIP'
                else
                    protein%residue(res_no)%res_name = 'HID'
                end if
            else if (res_name == 'ASP') then
                atom_name = line(13:16)
                if (atom_name(1:3) == ' HD') HD_read = .true.
                ! update protonation state of ASP
                if (HD_read) then
                    protein%residue(res_no)%res_name = 'ASH'
                end if
            else if (res_name == 'GLU') then
                atom_name = line(13:16)
                if (atom_name(1:3) == ' HE') HE_read = .true.
                ! update protonation state of GLU
                if (HE_read) then
                    protein%residue(res_no)%res_name = 'GLH'
                end if
            end if
        end if

    else if (line(1:6) == 'SSBOND') then
        link_no = link_no + 1
        protein%link(link_no)%link_name = 'DISU'
        call parse_string(line, num_word, word)
        if (num_word == 3) then
            protein%link(link_no)%chain(1) = ' '
            protein%link(link_no)%chain(2) = ' '
            read(word(2),"(I4)") protein%link(link_no)%link_res_no(1)
            read(word(3),"(I4)") protein%link(link_no)%link_res_no(2)
        else if (num_word == 5) then
            protein%link(link_no)%chain(1) = word(2)
            read(word(3),"(I4)") protein%link(link_no)%link_res_no(1)
            protein%link(link_no)%chain(2) = word(4)
            read(word(5),"(I4)") protein%link(link_no)%link_res_no(2)
        end if
        call find_res_idx(protein%link(link_no)%link_name, protein%link(link_no)%link_type)

    else if (line(1:7) == 'CISBOND') then
        call parse_string(line, num_word, word)
        cis_res_num = num_word-1
        do i = 2, num_word
            read(word(i),*) cis_res_list(i-1)
        end do
        n_h_lines = n_h_lines + 1
        if (n_h_lines > max_h_lines) then
            write(log_msg,"(A,I5)") 'Error. n_h_lines > max_h_lines', n_h_lines
            call terminate_with_error(log_msg)
        end if
        pdb_h_lines(n_h_lines) = line

    else if (line(1:8) == 'SYMMETRY') then
        symmetric = .true.
        n_h_lines = n_h_lines + 1
        if (n_h_lines > max_h_lines) then
            write(log_msg,"(A,I5)") 'Error. n_h_lines > max_h_lines', n_h_lines
            call terminate_with_error(log_msg)
        end if
        pdb_h_lines(n_h_lines) = line
          
    else if (line(1:3) == 'TER') then
        protein%residue(res_no)%ter_type = 'C'

    else if (line(1:6) == 'HETATM' .and. read_het) then
        read(line(23:26),*) het_curr_no
        het_curr = line(23:27)
        
        if (het_curr /= het_prev) then
            het_name = adjustl(line(18:20))
            het_prev = het_curr
            is_ligand = .false.
            do i_mol2 = 1, n_mol2_top
                if (het_name == trim(infile_mol2_topo(2,i_mol2))) then
                    is_ligand = .true.
                end if
            end do
         
            if (is_ligand) then
                lig_no = lig_no + 1
                if (lig_no > protein%n_lig) then
                    call reallocate_ligand_in_molecule(protein, protein%n_lig+max_lig)
                end if
                !
                protein%ligand(lig_no)%pdb_res_no = het_curr_no
                protein%ligand(lig_no)%res_added = het_curr(5:5)
                protein%ligand(lig_no)%pdb_res_name = het_name
                protein%ligand(lig_no)%res_name = het_name
                protein%ligand(lig_no)%chain = line(22:22)
            else
                het_no = het_no + 1
                if (het_no > protein%n_het) then
                    call reallocate_hetmol_in_molecule(protein, protein%n_het+max_het)
                end if
                !
                protein%hetmol(het_no)%pdb_res_no = het_curr_no
                protein%hetmol(het_no)%res_added = het_curr(5:5)
                protein%hetmol(het_no)%pdb_res_name = het_name
                protein%hetmol(het_no)%res_name = het_name
                protein%hetmol(het_no)%chain = line(22:22)
            end if
        end if

    else if (line(1:3) == 'END') then
        exit
       
    else if (line(1:6) == 'HEADER' .or. line(1:6) == 'COMPND' .or. line(1:6) == 'SOURCE' .or. &
             line(1:6) == 'REMARK' .or. line(1:6) == 'SEQRES' .or. line(1:5) == 'HELIX' .or. &
             line(1:5) == 'SHEET'  .or. line(1:4) == 'TURN') then
        if (line(1:10) == 'REMARK 350') then
            symmetric = .true.
            use_remark350 = .true.
        else
            n_h_lines = n_h_lines + 1
            if (n_h_lines > max_h_lines) then
                write(log_msg,"(A,I5)") 'Error. n_h_lines > max_h_lines', n_h_lines
                call terminate_with_error(log_msg)
            end if
            !
            pdb_h_lines(n_h_lines) = line
        end if
    end if
end do

close(f_unit)

if (res_no == 0 .and. het_no == 0 .and. lig_no == 0) then
    write(log_msg,"(A,3(1x,I5))") 'Error. Could not find ATOM or HETATM lines', res_no, het_no, lig_no
    call terminate_with_error(log_msg)
end if

protein%n_chain = chain_no
protein%n_link = link_no
!
protein%n_res = res_no
protein%n_het = het_no
protein%n_lig = lig_no
!
call reallocate_residue_in_molecule(protein, protein%n_res)
call reallocate_hetmol_in_molecule(protein,  protein%n_het)
call reallocate_ligand_in_molecule(protein,  protein%n_lig)

! modify the residue numbers in links
! (because res_no starts from 1, but it may not be the case in the input pdb.)
do link_no = 1, protein%n_link
    do i = 1, 2
        res_no = protein%link(link_no)%link_res_no(i)
        chain = protein%link(link_no)%chain(i)
        new_res_no_found = .false.
        do i_res = 1, protein%n_res
            if (res_no == protein%residue(i_res)%pdb_res_no .and. &
                 chain == protein%residue(i_res)%chain) then
                protein%link(link_no)%link_res_no(i) = i_res
                new_res_no_found = .true.
                exit
            end if
        end do
        if (.not. new_res_no_found) then
            write(log_msg,*) 'Error. Could not find the residue in SSBOND', res_no 
            call log_p(log_msg, me=me, level=30)
        end if
    end do
end do

! check for additional SS bond from SG coordinates
if (check_disulfide_pdb) then
    link_no = protein%n_link
    do i_res = 1, protein%n_res - 1
        res_name = protein%residue(i_res)%res_name
        if (res_name == 'CYS') then
            do j_res = i_res + 1, protein%n_res
                res_name = protein%residue(j_res)%res_name
                if (res_name == 'CYS') then
                    dr(:) = r_s(:,i_res) - r_s(:,j_res)
                    dr2 = dot_product(dr,dr)
                    if (dr2 <= ss2_cut) then
                        write(log_msg,"(A,2I5,F10.3)") 'Disulfide bond detected between', &
                              protein%residue(i_res)%pdb_res_no, &
                              protein%residue(j_res)%pdb_res_no, sqrt(dr2)
                        call log_p(log_msg, me=me, level=30)
                        new_ss_bond = .true.
                        do i = 1, protein%n_link
                            if (i_res == protein%link(i)%link_res_no(1) .or. &
                                i_res == protein%link(i)%link_res_no(2) .or. &
                                j_res == protein%link(i)%link_res_no(1) .or. &
                                j_res == protein%link(i)%link_res_no(2)) then
                                new_ss_bond = .false.
                            end if
                        end do
                        if (new_ss_bond) then
                            link_no = link_no + 1
                            protein%link(link_no)%link_name = 'DISU'
                            call find_res_idx('DISU', protein%link(link_no)%link_type)
                            protein%link(link_no)%link_res_no(1:2) = (/ i_res, j_res /)
                            protein%link(link_no)%chain(1) = protein%residue(i_res)%chain
                            protein%link(link_no)%chain(2) = protein%residue(j_res)%chain
                        end if
                    end if
                end if
            end do
        end if
    end do
    protein%n_link = link_no
end if
deallocate(r_s)

!to allow cis peptide bond (currently works for proline)
if (cis_res_num /= 0 .and. (top_type == 'eef1' .or. top_type == 'polarh')) then
    do i_cis_res = 1, cis_res_num
        do i_res = 1, protein%n_res
            if (i_res == cis_res_list(i_cis_res)) then
                ori_res_name = protein%residue(i_res)%res_name
                protein%residue(i_res)%res_name = 'C'//ori_res_name(1:2)//' '
                exit
            end if
        end do
    end do
end if

i_ss = 0
do link_no = 1, protein%n_link
    if (protein%link(link_no)%link_name == 'DISU') then
        i_ss = i_ss + 1
      
        ! modify the residue names in SSBOND
        ! This is only activated in all hydrogen topology.
        if (top_type == 'allh') then
            do i = 1, 2
                protein%link(link_no)%link_res_name(i) = 'CYX'
                res_no = protein%link(link_no)%link_res_no(i)
                if (protein%residue(res_no)%res_name /= "CYS") then
                    write(log_msg,"(A,A5,A,2I5)") 'ERR : SSBOND data is abnormal : ', &
                          protein%residue(res_no)%res_name, & 
                          " cannot have SSBOND!", res_no, protein%residue(res_no)%pdb_res_no
                    call log_p(log_msg, me=me, level=30)
                else
                    protein%residue(res_no)%res_name = 'CYX'
                end if
            end do
        end if
    end if
end do

write(log_msg,"(A,I4)") 'Number of disulfide bonds = ', i_ss
call log_p(log_msg,me=me,level=40)
call log_divider(me=me,level=40)

! chain index
chain_noprv = 0
do i_res = 1, protein%n_res
    chain_no = protein%residue(i_res)%i_chain
    if (chain_no > chain_noprv) then
        protein%chain_res(1,chain_no) = i_res
    end if
    protein%chain_res(2,chain_no) = i_res
    chain_noprv = chain_no
end do

protein_buffer = protein
infile_pdb_in_buffer(1) = pdb_file
if (n_model == 0) then
    n_model_in_buffer = 1
else
    n_model_in_buffer = n_model
end if

end subroutine read_sequence_from_pdb
!-------------------------------------------------------------------------------
subroutine read_pdb(pdb_file, molecule, model_no_in)
!-------------------------------------------------------------------------------
! The amino acid sequence is assumed to be known before calling read_pdb.
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: pdb_file
type(molecule_type), intent(inout) :: molecule
integer, intent(in), optional :: model_no_in
type(molecule_type), save :: protein_buffer
integer :: model_no
character(len=4) :: atom_name, pdb_atom_name, het_name
character(len=5) :: res_prev, res_curr, het_prev, het_curr
character(len=1) :: chain_prev, chain_curr
character(len=len_fname) :: line, word(25)
integer :: f_unit = 10, ioerror, openstat, num_word
integer :: res_no, atm_no, het_no, lig_no
integer :: i_ref_res, i_mol2, idx, i_symm
character(len=1) :: str
character(len=6) :: atom_error_mode
logical :: is_ligand
integer :: i_model, i_remark350
logical :: read_remark350

if (present(model_no_in)) then
    model_no = model_no_in
else
    model_no = 1
end if

if (allocated(residue_buffer)) then
    if (trim(pdb_file) == trim(infile_pdb_in_buffer(2))) then  ! use R in buffer
        molecule = protein_buffer
        call fill_coord_using_buffer(molecule, model_no)
        if (use_remark350) call set_symm_using_buffer(model_no)
        return
    else
        deallocate(residue_buffer)
        deallocate(hetmol_buffer)
        deallocate(ligand_buffer)
        if (use_remark350) deallocate(remark350_buffer)
    end if
end if
!
if (.not. allocated(residue_buffer)) then
    allocate(residue_buffer(3, max_atm,      molecule%n_res, n_model_in_buffer))
    allocate(hetmol_buffer (3, max_atm,      molecule%n_het, n_model_in_buffer))
    allocate(ligand_buffer (3, max_lig_atom, molecule%n_lig, n_model_in_buffer))
    if (use_remark350) allocate(remark350_buffer(10*max_symm, n_model_in_buffer))
end if

atom_error_mode = 'ignore'

call log_p('  Reading pdb file: '//trim(pdb_file),me=me,level=10)

open(f_unit, file=pdb_file, iostat=openstat)

if (openstat > 0) then
    write(log_msg,"(A,A)") 'Terminate with error: No pdbfile found, ', trim(pdb_file)
    call terminate_with_error(log_msg)
end if

! read ATOM lines only.
! ignore headers, disulfide bonds, and other information for now.
! read a single chain only

res_no = 0
res_prev = '-1000'     ! start with something unlikely.
het_no = 0
het_prev = '-1000'     ! start with something unlikely.
lig_no = 0
chain_prev= ''
i_model = 0
i_remark350 = 0
read_remark350 = .true.

do    
    read(f_unit, "(A120)", iostat=ioerror) line
    if (ioerror < 0) exit ! end of file

    if (multiple_models .and. line(1:5) == 'MODEL') then
        i_model = i_model + 1
        res_no = 0
        res_prev = '-1000'     ! start with something unlikely.
        het_no = 0
        het_prev = '-1000'     ! start with something unlikely.
        lig_no = 0
        chain_prev= ''
        i_remark350 = 0
        cycle
    else if (line(1:6) == 'ENDMDL') then
        if (multiple_models) then
            cycle
        else
            exit
        end if
    end if

    if (line(1:4) == 'ATOM') then
        ! res_no
        res_curr = line(23:27)
        chain_curr = line(22:22)
        if (i_model == 0) then
            i_model = 1
        end if
        if ((chain_curr /= chain_prev) .or. (res_curr /= res_prev)) then
            res_no = res_no+1
            res_prev = res_curr
            chain_prev = chain_curr
        end if
      
        str = line(17:17)
        if (str == ' ' .or. str == '1' .or. str == '2' .or. str == '3' .or. str == 'A') then
            ! skip lines something like BLYS (alternative conformation) for now.
            ! atom_no
            atom_name = line(13:16)
            pdb_atom_name = atom_name
            if (atom_name(1:1) /= 'H') then
                atom_name = trim(atom_name(2:4))//atom_name(1:1)
            end if
          
            call find_atom_idx(res_no, molecule%residue(res_no)%res_type, &
                               atom_name, atm_no, atom_error_mode)
            if (atm_no/=-100) then
                molecule%residue(res_no)%pdb_atom_name(atm_no) = pdb_atom_name
                molecule%residue(res_no)%atm_read(atm_no) = .true.
         
                ! coordinate
                read(line(31:), "(3f8.3)") molecule%residue(res_no)%R(1:3,atm_no)
                residue_buffer(1:3, atm_no, res_no, i_model) = molecule%residue(res_no)%R(1:3,atm_no)
          
                molecule%residue(res_no)%chain = line(22:22)
                molecule%residue(res_no)%code = line(75:78)
            end if
        end if
   
    ! Read Hetero atom only if read_het is turned on
    else if (line(1:6) == 'HETATM' .and. read_het) then
        ! res_no
        het_name = adjustl(line(18:20))
        het_curr = line(23:27)
        if (het_curr /= het_prev) then
            ! start of lines for a new residue
            is_ligand = .false.
            do i_mol2 = 1, n_mol2_top
                if (het_name == trim(infile_mol2_topo(2,i_mol2))) is_ligand = .true.
            end do
            
            if (is_ligand) then
                lig_no = lig_no + 1
            else
                het_no = het_no + 1
            end if
            het_prev = het_curr
        end if
        
        atom_name = line(13:16)
        pdb_atom_name = atom_name
        if (atom_name(1:1) == ' ') then
            atom_name = trim(atom_name(2:4))//atom_name(1:1)
        end if
      
        if (is_ligand) then
            i_ref_res = molecule%ligand(lig_no)%lig_type
            call find_ligand_atom_idx(lig_no, i_ref_res, atom_name, atm_no, atom_error_mode)
            if (atm_no == -100) cycle
            molecule%ligand(lig_no)%pdb_atom_name(atm_no) = pdb_atom_name
            molecule%ligand(lig_no)%atm_read(atm_no) = .true.
            
            ! coordinate
            read(line(31:), "(3F8.3)") molecule%ligand(lig_no)%R(1:3,atm_no)
            ligand_buffer(1:3, atm_no, lig_no, i_model) = molecule%ligand(lig_no)%R(1:3,atm_no)
            molecule%ligand(lig_no)%chain = line(22:22)
            molecule%ligand(lig_no)%code = line(75:78)
        else
            i_ref_res = molecule%hetmol(het_no)%res_type
            if (ref_res(i_ref_res)%n_atm /= molecule%hetmol(het_no)%n_atm) then
                ref_res(i_ref_res)%filled = .true.
            end if

            ! atom_no
            call find_atom_idx(het_no, molecule%hetmol(het_no)%res_type, atom_name, atm_no, atom_error_mode)
            molecule%hetmol(het_no)%pdb_atom_name(atm_no) = pdb_atom_name
            molecule%hetmol(het_no)%atm_read(atm_no) = .true.
         
            ! coordinate
            read(line(31:), "(3F8.3)") molecule%hetmol(het_no)%R(1:3,atm_no)
            molecule%hetmol(het_no)%chain = line(22:22)
            molecule%hetmol(het_no)%code = line(75:78)
            !
            hetmol_buffer(1:3, atm_no, het_no, i_model) = molecule%hetmol(het_no)%R(1:3,atm_no)
        end if
    else if (use_remark350 .and. line(1:10) == 'REMARK 350') then
        if (i_model == 0) then
            i_model = 1
        end if
        i_remark350 = i_remark350 + 1
        remark350_buffer(i_remark350, i_model) = line
        !
        call parse_string(line, num_word, word)
        if (word(3)(1:12) == 'BIOMOLECULE:') then
            read(word(4), '(I2)') idx
            if (idx /= 1) then
                read_remark350 = .false.
                write(log_msg, '(A,I2,A)') 'Warning: BIOMOLECULE: ', idx, ' is ignored.'
                call log_p(log_msg, me=me, level=20)
            end if
            cycle
        end if
        if (.not. read_remark350) cycle
        if (word(3)(1:5) == 'APPLY') then
            n_chain_symm_unit = 0
            do idx = 8, num_word
                n_chain_symm_unit = n_chain_symm_unit + 1
                symm_chain(n_chain_symm_unit) = word(idx)(1:1)
            end do
        else if (word(3)(1:5) == 'BIOMT') then
            read(word(3)(6:6),'(I1)') idx
            read(word(4),'(I2)') i_symm
            read(word(5),'(F10.6)') symm_U(idx,1,i_symm)
            read(word(6),'(F10.6)') symm_U(idx,2,i_symm)
            read(word(7),'(F10.6)') symm_U(idx,3,i_symm)
            read(word(8),'(F15.5)') symm_T(idx,i_symm)
            n_symm = i_symm
        end if

    else if (line(1:3) == 'END') then
        exit
   
    end if
end do

ref_res(1:num_ref_res)%filled = .true.

close(f_unit)

infile_pdb_in_buffer(2) = pdb_file
protein_buffer = molecule
if (multiple_models) then
    call fill_coord_using_buffer(molecule, model_no)
end if

end subroutine read_pdb
!-------------------------------------------------------------------------------
subroutine fill_coord_using_buffer(protein, model_no)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
integer, intent(in) :: model_no
integer :: i_model, i, i_res, i_atm

if (n_model_in_buffer == 1) then
    i_model = 1
else
    i_model = -1
    do i = 1, n_model_in_buffer
        if (model_no_in_buffer(i) == model_no) then
            i_model = i
            exit
        end if
    end do
    if (i_model == -1) then
        call terminate_with_error("  ERROR: There is no matched model_no")
    end if
end if

do i_res = 1, protein%n_res
    do i_atm = 1, protein%residue(i_res)%n_atm
        protein%residue(i_res)%R(1:3,i_atm) = residue_buffer(1:3,i_atm,i_res,i_model)
    end do
end do

do i_res = 1, protein%n_het
    do i_atm = 1, protein%hetmol(i_res)%n_atm
        protein%hetmol(i_res)%R(1:3,i_atm)  = hetmol_buffer(1:3,i_atm,i_res,i_model)
    end do
end do

do i_res = 1, protein%n_lig
    do i_atm = 1, protein%ligand(i_res)%n_atm
        protein%ligand(i_res)%R(1:3,i_atm)  = ligand_buffer(1:3,i_atm,i_res,i_model)
    end do
end do

end subroutine fill_coord_using_buffer
!-------------------------------------------------------------------------------
subroutine set_symm_using_buffer(model_no)
!-------------------------------------------------------------------------------
integer, intent(in) :: model_no
integer :: i_model, i, i_line, n_word, idx, i_symm
character(len=len_fname) :: line, word(25)
logical :: read_remark350

read_remark350 = .true.

if (n_model_in_buffer == 1) then
    i_model = 1
else
    i_model = -1
    do i = 1, n_model_in_buffer
        if (model_no_in_buffer(i) == model_no) then
            i_model = i
            exit
        end if
    end do
    if (i_model == -1) then
        call terminate_with_error("  ERROR: There is no matched model_no")
    end if
end if

do i_line = 1, 10*max_symm
    line = remark350_buffer(i_line, i_model)
    call parse_string(line, n_word, word)
    if (n_word < 3) cycle
    if (word(3)(1:12) == 'BIOMOLECULE:') then
        read(word(4), '(I2)') idx
        if (idx /= 1) then
            read_remark350 = .false.
            write(log_msg, '(A,I2,A)') 'Warning: BIOMOLECULE: ', idx, ' is ignored.'
            call log_p(log_msg, me=me, level=20)
        end if
        cycle
    end if
    if (.not. read_remark350) cycle
    if (word(3)(1:5) == 'APPLY') then
        n_chain_symm_unit = 0
        do idx = 8, n_word
            n_chain_symm_unit = n_chain_symm_unit + 1
            symm_chain(n_chain_symm_unit) = word(idx)(1:1)
        end do
    else if (word(3)(1:5) == 'BIOMT') then
        read(word(3)(6:6),'(I1)') idx
        read(word(4),'(I2)') i_symm
        read(word(5),'(F10.6)') symm_U(idx,1,i_symm)
        read(word(6),'(F10.6)') symm_U(idx,2,i_symm)
        read(word(7),'(F10.6)') symm_U(idx,3,i_symm)
        read(word(8),'(F15.5)') symm_T(idx,i_symm)
        n_symm = i_symm
    end if
end do

end subroutine set_symm_using_buffer            
!-------------------------------------------------------------------------------
subroutine read_ppDock_operations(infile_ppDock_opr, n_ppDock_oprs, res_i, res_j, Ts, Rs)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: infile_ppDock_opr
integer, intent(out) :: n_ppDock_oprs, res_i, res_j
real(dp), intent(out), allocatable :: Ts(:,:), Rs(:,:,:)

integer, parameter :: f_unit=176
integer :: io, n_word, k
character(len=len_fname) :: line, word(20)

open(f_unit, file=trim(infile_ppDock_opr), iostat=io)

read(f_unit, '(A120)', iostat=io) line
call parse_string(line, n_word, word)
read(word(2), '(I10)') n_ppDock_oprs

read(f_unit, '(A120)', iostat=io) line
call parse_string(line, n_word, word)
read(word(2), '(I10)') res_i
read(word(3), '(I10)') res_j

allocate(Ts(3,   n_ppDock_oprs))
allocate(Rs(3,3, n_ppDock_oprs))

k = 0
do
    read(f_unit, '(A120)', iostat=io) line
    if (io < 0) exit
    !
    call parse_string(line, n_word, word)
    !
    k = k + 1
    !
    read(word(1),'(F8.3)') Ts(1,k)
    read(word(2),'(F8.3)') Ts(2,k)
    read(word(3),'(F8.3)') Ts(3,k)
    !
    read(word( 4),'(F10.6)') Rs(1,1,k)
    read(word( 5),'(F10.6)') Rs(1,2,k)
    read(word( 6),'(F10.6)') Rs(1,3,k)
    read(word( 7),'(F10.6)') Rs(2,1,k)
    read(word( 8),'(F10.6)') Rs(2,2,k)
    read(word( 9),'(F10.6)') Rs(2,3,k)
    read(word(10),'(F10.6)') Rs(3,1,k)
    read(word(11),'(F10.6)') Rs(3,2,k)
    read(word(12),'(F10.6)') Rs(3,3,k)
end do

close(f_unit)

end subroutine read_ppDock_operations
!-------------------------------------------------------------------------------
subroutine finalize_structure()
!-------------------------------------------------------------------------------
if (allocated(residue_buffer)) then
    deallocate(residue_buffer)
    deallocate(hetmol_buffer)
    deallocate(ligand_buffer)
end if
if (allocated(remark350_buffer)) deallocate(remark350_buffer)

end subroutine finalize_structure
!-------------------------------------------------------------------------------
subroutine n_models_in_pdb(infile_pdb, n_model, model_no)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: infile_pdb
integer, intent(out) :: n_model, model_no(max_pdbfile)

integer,parameter :: f_unit=10
integer :: ioerr, n_word
character(len=len_fname) :: line, word(10)

open(f_unit, file=trim(infile_pdb), iostat=ioerr)
if (ioerr /= 0) then
    write(log_msg,"(A,A)") 'Terminate with error: No pdbfile found, ', trim(infile_pdb)
    call terminate_with_error(log_msg)
end if

n_model = 0
model_no(:) = 1
do
    read(f_unit, '(A120)', iostat=ioerr) line
    if (ioerr < 0) exit
    !
    if (line(1:5) == 'MODEL') then
        n_model = n_model + 1
        call parse_string(line, n_word, word)
        read(word(2), "(I5)") model_no(n_model)
    end if
end do

if (n_model == 0) n_model = 1

close(f_unit)

end subroutine n_models_in_pdb
!-------------------------------------------------------------------------------
subroutine write_protein_to_pdb(pdb_unit, molecule)
!-------------------------------------------------------------------------------
integer, intent(in) :: pdb_unit
type(molecule_type), intent(in) :: molecule
!
integer :: i_res, i_ref_res, i_atm, k, res_no, i_het, i_lig
character(len=4) :: res_name, atom_name
character(len=1) :: prev_chain

k = 0
prev_chain = ''
do i_res = 1, molecule%n_res 
    i_ref_res = molecule%residue(i_res)%res_type
    res_name = molecule%residue(i_res)%pdb_res_name
    res_no = molecule%residue(i_res)%pdb_res_no
    
    if (symmetric) then
        if (i_res > symm_resrange(2,1)) exit
    end if

    if (i_res > 1 .and. molecule%residue(i_res)%chain /= prev_chain) then
        write(pdb_unit, '(A3)') 'TER'
    end if
    prev_chain = molecule%residue(i_res)%chain

    do i_atm = 1, ref_res(i_ref_res)%n_atm
        k = k + 1
        !
        if (molecule%residue(i_res)%atm_read(i_atm)) then
            atom_name = molecule%residue(i_res)%pdb_atom_name(i_atm)
        else
            atom_name = ref_res(i_ref_res)%atom_name(i_atm)
            atom_name = atom_name(4:4) // atom_name(1:3)
        end if
        write(pdb_unit, 72) 'ATOM  ', k, atom_name, res_name, & 
            molecule%residue(i_res)%chain, res_no, & 
            molecule%residue(i_res)%res_added, &
            molecule%residue(i_res)%R(:,i_atm)
    end do
end do

write(pdb_unit, '(A3)') 'TER'

do i_het = 1, molecule%n_het
    i_ref_res = molecule%hetmol(i_het)%res_type
    res_name = molecule%hetmol(i_het)%pdb_res_name
    res_no = molecule%hetmol(i_het)%pdb_res_no
   
    do i_atm = 1, ref_res(i_ref_res)%n_atm
        k = k + 1
        if (molecule%hetmol(i_het)%atm_read(i_atm)) then
            atom_name = molecule%hetmol(i_het)%pdb_atom_name(i_atm)
        else
            atom_name = ref_res(i_ref_res)%atom_name(i_atm)
            atom_name = atom_name(4:4)//atom_name(1:3)
        end if
        write(pdb_unit, 72) 'HETATM', k, atom_name, res_name, & 
            molecule%hetmol(i_het)%chain, res_no, & 
            molecule%hetmol(i_het)%res_added, molecule%hetmol(i_het)%R(:,i_atm)
    end do
end do

do i_lig = 1, molecule%n_lig
    i_ref_res = molecule%ligand(i_lig)%lig_type
    res_name  = molecule%ligand(i_lig)%pdb_res_name
    res_no    = molecule%ligand(i_lig)%pdb_res_no
   
    do i_atm = 1, ref_lig(i_ref_res)%n_atm
        k = k + 1
        atom_name = ref_lig(i_ref_res)%atom_name(i_atm)
        atom_name = atom_name(4:4)//atom_name(1:3)
        write(pdb_unit, 72) 'HETATM', k, atom_name, res_name, & 
            molecule%ligand(i_lig)%chain, res_no, & 
            molecule%ligand(i_lig)%res_added, molecule%ligand(i_lig)%R(:,i_atm)
    end do
end do

if (molecule%n_lig + molecule%n_het > 0) then
    write(pdb_unit,"(A)") "TER" 
end if

72  format(a6,i5,1x,a4,1x,a3,1x,a1,i4,a1,3x,3f8.3)

!-------------------------------------------------------------------------------
end subroutine write_protein_to_pdb
!-------------------------------------------------------------------------------
subroutine write_pdb(pdb_unit, molecule, n_remark, remark)
!-------------------------------------------------------------------------------
! Write PDB file for the molecule given
!-------------------------------------------------------------------------------
integer, intent(in) :: pdb_unit
type(molecule_type), intent(in) :: molecule
integer, intent(in), optional :: n_remark
character(len=len_log), optional :: remark(:)
integer :: i

if (present(n_remark)) then
    do i = 1, n_remark
        write(pdb_unit, "(A,1x,A)") "REMARK", trim(remark(i))
    end do
end if

if (symmetric) then ! write chainA only. Add REMARK 350
    call write_symm_info(pdb_unit, molecule)
end if

do i = 1, n_h_lines
    write(pdb_unit,"(a)") trim(pdb_h_lines(i))
end do

! write ss bnd
do i = 1, molecule%n_link
    if (molecule%link(i)%link_name == 'DISU') then
        write(pdb_unit,69) molecule%link(i)%chain(1), &
             molecule%residue(molecule%link(i)%link_res_no(1))%pdb_res_no, &
             molecule%link(i)%chain(2), &
             molecule%residue(molecule%link(i)%link_res_no(2))%pdb_res_no
    end if
end do
69  format('SSBOND',1X,A1,2X,I3,4X,1X,A1,2X,I3)

call write_protein_to_pdb(pdb_unit, molecule)

end subroutine write_pdb
!-------------------------------------------------------------------------------
subroutine write_pdb_model(pdb_unit, molecule, i_conf, n_remark, remark)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule
integer, intent(in) :: pdb_unit, i_conf
integer, intent(in), optional :: n_remark
character(len=len_log), intent(in), optional :: remark(:)
integer :: i

if (i_conf == 1) then
    do i = 1, n_h_lines
        write(pdb_unit,"(a)") trim(pdb_h_lines(i))
    end do

    ! write ss bnd
    do i = 1, molecule%n_link
        if (molecule%link(i)%link_name == 'DISU') then
            write(pdb_unit,69) molecule%link(i)%chain(1), &
                 molecule%residue(molecule%link(i)%link_res_no(1))%pdb_res_no, &
                 molecule%link(i)%chain(2), &
                 molecule%residue(molecule%link(i)%link_res_no(2))%pdb_res_no
        end if
    end do
    69  format('SSBOND',1X,A1,2X,I3,4X,1X,A1,2X,I3)
end if

write(pdb_unit, '(A,1x,I4)') "MODEL", i_conf

if (present(n_remark)) then
    do i = 1, n_remark
        write(pdb_unit, "(A,1x,A)") "REMARK", trim(remark(i))
    end do
end if

if (symmetric) then ! write chainA only. Add REMARK 350
    call write_symm_info(pdb_unit, molecule)
end if

call write_protein_to_pdb(pdb_unit, molecule)
write(pdb_unit, '(A)') "ENDMDL"

end subroutine write_pdb_model
!-------------------------------------------------------------------------------
subroutine open_write_pdb(pdb_unit, pdb_file)
!-------------------------------------------------------------------------------
integer, intent(in) :: pdb_unit
character(len=len_fname), intent(in) :: pdb_file

call log_p('Writing pdb file: '//trim(pdb_file), level=40)
open(pdb_unit, file=trim(pdb_file))

end subroutine open_write_pdb
!-------------------------------------------------------------------------------
subroutine close_write_pdb(pdb_unit)
!-------------------------------------------------------------------------------
integer, intent(in) :: pdb_unit

write(pdb_unit, '(A)') 'END'
close(pdb_unit)

end subroutine close_write_pdb
!-------------------------------------------------------------------------------
subroutine write_pdbm(molecule)
!-------------------------------------------------------------------------------
! TODO: Have to write down!
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule

end subroutine write_pdbm
!-------------------------------------------------------------------------------
END MODULE IN_OUT_STRUCTURE
!-------------------------------------------------------------------------------
