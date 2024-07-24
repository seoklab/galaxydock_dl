!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/io/setup_molecule.f90
!
! Description:
!  This module contains subroutines fundamental for initializing molecule_type.
!
!-------------------------------------------------------------------------------
MODULE SETUP_MOLECULE
!-------------------------------------------------------------------------------
use globals
use allocate_molecule
use logger,           only: log_p, terminate_with_error, log_divider
!
use in_out_vars
use in_out_utils,     only: find_res_idx, find_het_idx, find_ligand_idx, &
                            find_atom_idx, find_bnd_prm, find_ang_prm
use in_out_input,     only: reinitialize_usc
use in_out_structure, only: read_sequence_from_pdb, read_pdb
use in_out_ligand,    only: build_tree, get_ligand_center
!
use sort,             only: sort2
use rmsd,             only: calc_rmsd
use symmetry,         only: setup_symmetry, n_symm, n_chain_symm_unit, &
                            set_USC_for_symm, set_ULR_for_symm, &
                            apply_fix_type_symm
use mathfunctions,    only: cross, bound_ang, v_norm
use convert_res_name, only: find_long_res_name
use geometry,         only: calc_tor_angle, internal2cartesian, cartesian2internal, &
                            internal2cartesian_het, cartesian2internal_het, &
                            quaternion, q_product, rotation_matrix

implicit none
private

public :: initialize_molecule
!
public :: place_missing_atoms
public :: place_missing_atoms_het
public :: place_missing_atoms_lig
public :: c2i_i2c_reconstruction_het
public :: place_an_atom
public :: mutate_aa
public :: setup_tot_num
public :: apply_fix_type
  
CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_molecule(molecule, ligand, model_no)
!-------------------------------------------------------------------------------
! Setup a molecule_type based on given coordinate and reference topology
! 1. read sequence and pdb file(if inmol_filetype is pdb)
! 2. setup link which indice for 3 residues in the previous residue that is
!    linked to current residue
! 3. construct structure from canonical geometry
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
type(ligand_type), intent(inout), optional :: ligand
integer, intent(in), optional :: model_no
integer :: res_no, n_unplaced, het_no, lig_no

call allocate_molecule_type(molecule)

!get sequence & pdb file (if exists)
if (mol_type == 'protein') then
    call read_sequence_from_pdb(infile_pdb, molecule)
    call setup_residue_indices(molecule)
end if
if (present(model_no)) then
    call read_pdb(infile_pdb, molecule, model_no_in=model_no)
else
    call read_pdb(infile_pdb, molecule)
end if

! Set link_type variable inside molecule
! Find indices for the three atoms in the previous residue
! that are directly linked to the current residue.
call setup_link_btw_res(molecule)

call log_p('Finished with setup_link_btw_res.', me=me, level=40)

do res_no = 1, molecule%n_res
    ! set all the atoms free (for the purpose of E calc) for now.
    molecule%residue(res_no)%atom_fixed(:) = .false.
    ! missing atoms are not placed yet.
    molecule%residue(res_no)%atm_placed(:) = .false.
end do

do het_no = 1, molecule%n_het
    ! set all the atoms free (for the purpose of E calc) for now.
    molecule%hetmol(het_no)%atom_fixed(:) = .false.
    ! missing atoms are not placed yet.
    molecule%hetmol(het_no)%atm_placed(:) = .false.
end do

do lig_no = 1, molecule%n_lig
    ! set all the atoms free (for the purpose of E calc) for now.
    molecule%ligand(lig_no)%atom_fixed(:) = .false.
    ! missing atoms are not placed yet.
    molecule%ligand(lig_no)%atm_placed(:) = .false.
end do

call place_missing_atoms(molecule, 1, molecule%n_res, n_unplaced)
if (n_unplaced > 0) then
    call c2i_i2c_reconstruction(molecule)
    do res_no = 1, molecule%n_res
        molecule%residue(res_no)%atm_placed(:) = .true.
    end do
end if
if (molecule%n_het > 0) then
    call place_missing_atoms_het(molecule, 1, molecule%n_het, n_unplaced)
    if (n_unplaced > 0) then
        call c2i_i2c_reconstruction_het(molecule)
        do het_no = 1, molecule%n_het
            molecule%hetmol(het_no)%atm_placed(:) = .true.
        end do
    end if
end if
if (molecule%n_lig > 0) then
    call place_missing_atoms_lig(molecule, 1, molecule%n_lig, n_unplaced)
    if (n_unplaced > 0) then
        ! TODO: have to write _lig version (?)
        call c2i_i2c_reconstruction_het(molecule)
        do het_no = 1, molecule%n_het
            molecule%hetmol(het_no)%atm_placed(:) = .true.
        end do
    end if
end if

! setup ligand for ligand docking
if (present(ligand) .and. (infile_ligand(1) /= '')) then
    call setup_docking_ligand(molecule, ligand, infile_ligand(2)(1:3))
end if

call setup_tot_num(molecule)

!symmetric: defined at globals. sets at in_out_structure
if (n_symm > 1) symmetric = .true.
if (symmetric) call reinitialize_molecule_w_symm(molecule)
if (n_symm == 0) n_symm = 1

call apply_fix_type(molecule)

end subroutine initialize_molecule
!-------------------------------------------------------------------------------
subroutine reinitialize_molecule_w_symm(molecule)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule

call setup_symmetry(molecule)
if (use_remark350) then
    call setup_link_btw_res(molecule)
    call setup_tot_num(molecule)
end if

call cartesian2internal(1, molecule%n_res, molecule%residue(1:molecule%n_res))
call internal2cartesian(1, molecule%n_res, molecule%residue(1:molecule%n_res))

!do i_res = 1, molecule%n_res
!    print*, 'Residue: ', i_res
!    do i_atm = -1, molecule%residue(i_res)%n_atm
!        write(*,'(A,3F10.5)'), '  - bond vector: ', molecule%residue(i_res)%b(:, i_atm)
!    end do
!    do i_atm = 0, molecule%residue(i_res)%n_atm
!        write(*,'(A,4F10.5)'), '  - quat: ', molecule%residue(i_res)%quat(:, i_atm)
!    end do
!    do i_atm = 1, molecule%residue(i_res)%n_atm
!        write(*,'(A,3F20.5)'), '  - len & ang: ', &
!            molecule%residue(i_res)%b_len(i_atm), &
!            molecule%residue(i_res)%b_ang(i_atm), &
!            molecule%residue(i_res)%t_ang(i_atm)
!    end do
!    print*, ' - link residue: ', molecule%residue(i_res)%link_res_no(:)
!    print*, ' - link atom: ', molecule%residue(i_res)%link_atm_no(:)
!end do
!stop

end subroutine reinitialize_molecule_w_symm
!-------------------------------------------------------------------------------
subroutine setup_docking_ligand(molecule, ligand, lig_name)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
type(ligand_type), intent(inout) :: ligand
character(len=3), intent(in) :: lig_name
integer :: i_lig, i_atm, ref_lig_no, n_atm

ref_lig_no = -1
do i_lig = 1, n_mol2_top
    if (ref_lig(i_lig)%lig_name == lig_name) then
        ref_lig_no = i_lig
        ligand%lig_type = ref_lig_no
        ligand%n_atm = ref_lig(ref_lig_no)%n_atm
        exit
    end if
end do
if (ref_lig_no == -1) then
    write(log_msg, '(A,A)') 'ERROR: Please check ligand name: ', lig_name
    call terminate_with_error(log_msg)
end if

call get_ligand_center(ref_lig(ref_lig_no))
call build_tree(ligand)

! Reallocate ligand size in molecule
i_lig = molecule%n_lig + 1
call reallocate_ligand_in_molecule(molecule, molecule%n_lig+1)
call allocate_ligand_residue_type(molecule%ligand(i_lig), ref_lig(ref_lig_no)%n_atm)

! copy ligand to molecule%ligand
ligand%lig_no = i_lig
molecule%n_lig = i_lig
molecule%ligand(i_lig)%lig_type = ref_lig_no
molecule%ligand(i_lig)%pdb_res_no = i_lig
molecule%ligand(i_lig)%res_added = ' '
molecule%ligand(i_lig)%pdb_res_name = lig_name
molecule%ligand(i_lig)%res_name = ref_lig(ref_lig_no)%lig_name
molecule%ligand(i_lig)%chain = ' '
molecule%ligand(i_lig)%n_atm = ref_lig(ref_lig_no)%n_atm

molecule%ligand(i_lig)%R(1:3,0) = 0.0d0
do i_atm = 1, ref_lig(ref_lig_no)%n_atm
    molecule%ligand(i_lig)%pdb_atom_name(i_atm) = &
                            ref_lig(ref_lig_no)%atom_name(i_atm)
    molecule%ligand(i_lig)%atm_read(i_atm) = .true.
    molecule%ligand(i_lig)%R(1:3, i_atm) = ref_lig(ref_lig_no)%R(1:3, i_atm)
end do
n_atm = ref_lig(ref_lig_no)%n_atm
molecule%ligand(i_lig)%b_len(1:2*n_atm) = ref_lig(ref_lig_no)%b_len0(1:2*n_atm)
molecule%ligand(i_lig)%b_ang(1:2*n_atm) = ref_lig(ref_lig_no)%b_ang0(1:2*n_atm)
molecule%ligand(i_lig)%t_ang(1:4*n_atm) = ref_lig(ref_lig_no)%d_ang0(1:4*n_atm)

write(log_msg, '(A,A)') 'Setup docking ligand is done. ', lig_name
call log_p(log_msg, level=10)

end subroutine setup_docking_ligand
!-------------------------------------------------------------------------------
subroutine setup_residue_indices(protein)
!-------------------------------------------------------------------------------
! Set res_type index of the protein by searching ref_res
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
character(len=4) :: res_name
integer :: res_no, het_no

! Mark terminal residues
! N terminal
res_no = 1
res_name = protein%residue(res_no)%res_name
if (res_name == 'ACE') then
    protein%n_capped = .true.
    protein%residue(res_no)%ter_type = 'N'
else
    protein%n_capped = .false.
    protein%residue(res_no)%ter_type = 'N'
end if

! C terminal
res_no = protein%n_res
res_name = protein%residue(res_no)%res_name
protein%c_capped = .false.
if (res_name == 'NME' .or. res_name == 'NHE') then
    protein%c_capped = .true.
    protein%residue(res_no)%ter_type = 'C'
else
    protein%c_capped = .true.
    protein%residue(res_no)%ter_type = 'C'
end if

! sequence breaks
do res_no = 1, protein%n_res
    if (protein%residue(res_no)%ter_type == 'N') then
        protein%residue(res_no)%res_name = 'N'//protein%residue(res_no)%res_name
    else if (protein%residue(res_no)%ter_type == 'C') then
        protein%residue(res_no)%res_name = 'C'//protein%residue(res_no)%res_name
    end if
    ! topology not available for NGLH and CGLH
    ! change to GLU for now...
    if (protein%residue(res_no)%res_name == 'NGLH') protein%residue(res_no)%res_name = 'NGLU'
    if (protein%residue(res_no)%res_name == 'CGLH') protein%residue(res_no)%res_name = 'CGLU'
end do

! assume non-cyclic for now.
protein%cyclic = .false.

! find residue index
do res_no = 1, protein%n_res
    call find_res_idx(protein%residue(res_no)%res_name, protein%residue(res_no)%res_type)
    call allocate_residue_type(protein%residue(res_no),&
                               ref_res(protein%residue(res_no)%res_type)%n_atm)
end do

do res_no = 1, protein%n_het
    call find_het_idx(protein%hetmol(res_no))
    call allocate_residue_type(protein%hetmol(res_no),&
                               ref_res(protein%hetmol(res_no)%res_type)%n_atm)
end do

do res_no = 1, protein%n_lig
    call find_ligand_idx(protein%ligand(res_no))
    call allocate_ligand_residue_type(protein%ligand(res_no),&
                                      ref_lig(protein%ligand(res_no)%lig_type)%n_atm)
end do

! report sequences
if (print_level_global >= 40 .and. me==king) then
    call log_divider(me=me, level=40)
    call log_p('Amino acid sequence:', me=me, level=40)
    do res_no = 1, protein%n_res
        write(log_msg,"(I6,2X,A4)") res_no, ref_res(protein%residue(res_no)%res_type)%res_name
        call log_p(log_msg, me=me, level=40)
    end do
    if (protein%n_het > 0) then
        call log_p('Hetero-molecules:', me=me, level=40)
        do het_no = 1, protein%n_het
            write(log_msg,"(I6,2X,A4)") het_no, ref_res(protein%hetmol(het_no)%res_type)%res_name
            call log_p(log_msg, me=me, level=40)
        end do
    end if
    if (protein%n_lig > 0) then
        call log_p('Ligand-molecules:', me=me, level=40)
        do het_no = 1, protein%n_lig
            write(log_msg,"(I6,2X,A4)") het_no, ref_lig(protein%ligand(het_no)%lig_type)%lig_name
            call log_p(log_msg, me=me, level=40)
        end do
    end if
    call log_divider(me=me, level=40)
end if

end subroutine setup_residue_indices
!-------------------------------------------------------------------------------
subroutine mutate_aa(molecule, nres, reslist, newaa)
!-------------------------------------------------------------------------------
! Mutate [nres] residues into newly assigned amino acids in [newaa]
! Caution: This is only for 'molecule'. For energy evaluation, corresponding
! arrays for example R, i_R, ... MUST BE REFRESHED afterwards!!
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer, intent(in) :: nres, reslist(nres)
character(len=1), intent(in) :: newaa(nres)
integer, parameter :: nbb = 6
character(len=4), parameter :: bbatom(nbb) = (/ 'N   ', 'CA  ', 'C   ', 'CB  ', 'O   ', 'OXT ' /)
logical :: atom_built(max_atm)
integer :: i_atm, i_res, resno, atmno1, atmno2
integer :: i_ref, i_ref_old
type(residue_type) :: residue_new, residue_save
character(len=1) :: aa
character(len=4) :: laa, atom_name
character(len=6) :: error_mode

error_mode = 'ignore'
call cartesian2internal(1, molecule%n_res, molecule%residue(1:molecule%n_res))

do i_res = 1, nres
    aa = newaa(i_res)
    resno = reslist(i_res)
    call find_long_res_name(aa, laa)
    if (molecule%residue(resno)%ter_type == 'N') laa = 'N'//laa(1:3)
    if (molecule%residue(resno)%ter_type == 'C') laa = 'C'//laa(1:3)
   
    ! Set new and old ref_res indices
    residue_save = molecule%residue(resno)
    residue_new = residue_save
   
    !! PHB 11/11/25: misc. arrays in residue_type like link_atm_no not copied.
    !! Would it be problematic someday? (otherwise you can just gen pdb and read again...)
    i_ref_old = residue_save%res_type
    call find_res_idx(laa, i_ref)
    residue_new%res_type = i_ref
    residue_new%res_name = ref_res(i_ref)%res_name
    residue_new%pdb_res_name = ref_res(i_ref)%res_name
   
    ! Leave backbone coordinate as it is
    residue_new%b_len(:) = 0.0d0
    residue_new%b_ang(:) = 0.0d0
    residue_new%t_ang(:) = 0.0d0
    atom_built(:) = .false.
   
    do i_atm = 1, nbb
        call find_atom_idx(resno, i_ref_old, bbatom(i_atm), atmno1, error_mode)
        call find_atom_idx(resno, i_ref, bbatom(i_atm), atmno2, error_mode)
        if (atmno1 < 0 .or. atmno2 < 0) cycle
        residue_new%b_len(atmno2) = residue_save%b_len(atmno1)
        residue_new%b_ang(atmno2) = residue_save%b_ang(atmno1)
        residue_new%t_ang(atmno2) = residue_save%t_ang(atmno1)
        residue_new%quat(:,atmno2) = residue_save%quat(:,atmno1)
        atom_name = ref_res(i_ref)%atom_name(atmno2)
        atom_name = atom_name(4:4) // atom_name(1:3)
        residue_new%pdb_atom_name(atmno2) = atom_name
        
        atom_built(atmno2) = .true.
    end do
   
    ! Copy default internal coordinates for remaining sidechain atoms
    do i_atm = 1, ref_res(i_ref)%n_atm
        if (.not. atom_built(i_atm)) then
            residue_new%b_len(i_atm) = ref_res(i_ref)%b_len0(i_atm)
            residue_new%b_ang(i_atm) = ref_res(i_ref)%b_ang0(i_atm)
            residue_new%t_ang(i_atm) = ref_res(i_ref)%t_ang0(i_atm)
            atom_name = ref_res(i_ref)%atom_name(i_atm)
            atom_name = atom_name(4:4) // atom_name(1:3)
            residue_new%pdb_atom_name(i_atm) = atom_name
            atom_built(i_atm) = .true.
        end if
    end do
   
    residue_new%atm_read(:) = .true.
    residue_new%atm_placed(:) = .true.
    ! Finally copy the new residue generated
    molecule%residue(resno) = residue_new
end do

! Refresh links
call setup_link_btw_res(molecule)

! Fully reconstruct
call internal2cartesian(1, molecule%n_res, molecule%residue(1:molecule%n_res))
call cartesian2internal(1, molecule%n_res, molecule%residue(1:molecule%n_res))

end subroutine mutate_aa
!-------------------------------------------------------------------------------
subroutine c2i_i2c_reconstruction(molecule)
!-------------------------------------------------------------------------------
! make internal coordinate from cartesian coordinate,
! then make cartesian coordinate from internal coordinate.
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
type(molecule_type) :: molecule_initial
character(len=2) :: rmsd_type
real(dp) :: U(3,3), error
integer :: rmsd_option

molecule_initial = molecule !for reconstruction RMSD test

!make internal coordinate from cartesian coordinate
!cartesian2internal(start_residue no, end_residue no,residue type data)
call cartesian2internal(1, molecule%n_res, molecule%residue(1:molecule%n_res))
call log_p('Finished with cartesian2internal.', me=me, level=20)
 
!make cartesian coordinate from internal coordinate
call internal2cartesian(1, molecule%n_res, molecule%residue(1:molecule%n_res))
call log_p('Finished with internal2cartesian.', me=me, level=20)

!! report reconstruction test result
rmsd_type = 'ha'
rmsd_option = 0
call calc_rmsd(molecule_initial, molecule, rmsd_type, rmsd_option, U, error)
call log_divider(me=me)
write(log_msg,"(A,F9.3)") 'All-atom RMSD after reconstruction =', error
call log_p(log_msg, me=me)
call log_divider(me=me)

end subroutine c2i_i2c_reconstruction
!-------------------------------------------------------------------------------
subroutine c2i_i2c_reconstruction_het(molecule)
!-------------------------------------------------------------------------------
! make internal coordinate from cartesian coordinate,
! then make cartesian coordinate from internal coordinate. (for hetero molecule)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
!type(molecule_type) :: molecule_initial
integer :: i_het

!molecule_initial = molecule

do i_het = 1, molecule%n_het
    !make internal coordinate from cartesian coordinate
    call cartesian2internal_het(molecule%hetmol(i_het))
    !make cartesian coordinate from internal coordinate
    call internal2cartesian_het(molecule%hetmol(i_het))
end do
call log_p('  Finished with i2c and c2i for hetero molecules.', me=me, level=20)

end subroutine c2i_i2c_reconstruction_het
!-------------------------------------------------------------------------------
subroutine place_missing_atoms(molecule, i_res1, i_res2, n_unplaced)
!-------------------------------------------------------------------------------
! Place missing atoms in the molecule using topology information
! defined on top_eng file
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer, intent(in) :: i_res1, i_res2
integer, intent(out) :: n_unplaced
integer :: i_add, res_no, atm_no, ref_res_no, i_add_heavy, ierr
character(len=4) :: atom
logical :: placed

i_add = 0
i_add_heavy = 0
n_unplaced = 0

call log_p('- Add missing atoms activated.',me=me,level=30)

! place heavy atoms first
do res_no = i_res1, i_res2
    ref_res_no = molecule%residue(res_no)%res_type
    do atm_no = 1, ref_res(ref_res_no)%n_atm
        if (.not. molecule%residue(res_no)%atm_read(atm_no) .and. &
            .not. molecule%residue(res_no)%atm_placed(atm_no)) then
            atom = ref_res(ref_res_no)%atom_name(atm_no)
           
            if (atom(1:1) /= 'H') then
                call place_an_atom(molecule, res_no, atm_no, placed, ierr)
                if (ierr == 1) then
                    write(log_msg,"(A,I4,1X,A4)") 'ERROR (PLACE_MISSING_ATOM): Parent Bond/Angle not found for residue ', res_no, atom
                    call terminate_with_error(log_msg)
                else if (ierr == 2) then
                    write(log_msg,"(A,I4,1X,A4)") 'ERROR (PLACE_MISSING_ATOM): Hybridization not defined for ', res_no, atom
                    call terminate_with_error(log_msg)
                else if (ierr == 3) then
                    write(log_msg,"(A,I4,1X,A4)") 'ERROR (PLACE_MISSING_ATOM): Connected atom not defined for ', res_no, atom
                    call terminate_with_error(log_msg)
                end if
               
                if (placed) then
                    write(log_msg,"(A,I4,1X,A4,1X,A4,1X,A)") '  Placed non-hydrogen atom ', & 
                          res_no, ref_res(ref_res_no)%res_name, atom
                    call log_p(log_msg, me=me, level=30)
                    i_add = i_add + 1
                    i_add_heavy = i_add_heavy + 1
                    molecule%residue(res_no)%atm_placed(atm_no) = .true.
                end if
            end if
        end if
    end do
end do

! place hydrogen atom
do res_no = i_res1, i_res2
    ref_res_no = molecule%residue(res_no)%res_type
    do atm_no = 1, ref_res(ref_res_no)%n_atm
        if (.not.molecule%residue(res_no)%atm_read(atm_no) .and. &
            .not.molecule%residue(res_no)%atm_placed(atm_no)) then
            atom = ref_res(ref_res_no)%atom_name(atm_no)
            if (atom(1:1) == 'H') then
                call place_an_atom(molecule, res_no, atm_no, placed, ierr)
              
                if (ierr == 1) then
                    write(log_msg,"(A,I4,1X,A4)") 'ERROR (PLACE_MISSING_ATOM): Parent Bond/Angle not found for residue ', res_no, atom
                    call terminate_with_error(log_msg)
                else if (ierr == 2) then
                    write(log_msg,"(A,I4,1X,A4)") 'ERROR (PLACE_MISSING_ATOM): Hybridization not defined for ', res_no, atom
                    call terminate_with_error(log_msg)
                else if (ierr == 3) then
                    write(log_msg,"(A,I4,1X,A4)") 'ERROR (PLACE_MISSING_ATOM): Connected atom not defined for ', res_no, atom
                    call terminate_with_error(log_msg)
                end if
              
                if (placed) then
                    i_add = i_add + 1
                    molecule%residue(res_no)%atm_placed(atm_no) = .true.
                else
                    n_unplaced = n_unplaced + 1
                end if
            end if
        end if
    end do
end do

write(log_msg,"(A,I5)") '  Number of added heavy atoms = ', i_add_heavy
call log_p(log_msg, me=me, level=30)
write(log_msg,"(A,I5)") '  Total number of added atoms = ', i_add
call log_p(log_msg, me=me, level=30)
write(log_msg,"(A,I5)") '  n_unplaced = ', n_unplaced
call log_p(log_msg, me=me, level=30)
  
end subroutine place_missing_atoms
!-------------------------------------------------------------------------------
subroutine place_missing_atoms_het(molecule, i_het1, i_het2, n_unplaced)
!-------------------------------------------------------------------------------
! Similar to place_missing_atoms, hetero molecule version
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer, intent(in) :: i_het1, i_het2
integer, intent(out) :: n_unplaced
integer :: i_add, atm_no, ref_res_no, i_add_heavy, het_no
character(len=4) :: atom
logical :: placed

i_add = 0
i_add_heavy = 0
n_unplaced = 0

  ! place heavy atoms first
do het_no = i_het1, i_het2
    ref_res_no = molecule%hetmol(het_no)%res_type
    do atm_no = 1, ref_res(ref_res_no)%n_atm
        if (.not. molecule%hetmol(het_no)%atm_read(atm_no) .and. &
            .not. molecule%hetmol(het_no)%atm_placed(atm_no)) then
            atom = ref_res(ref_res_no)%atom_name(atm_no)
          
            if (atom(1:1) /= 'H') then
                call place_an_atom_het(molecule, het_no, atm_no, placed)
                if (placed) then
                    write(log_msg,"(A,I4,1X,A4,1X,A4,1X,A)") '  Placed non-hydrogen atom ', & 
                          het_no, ref_res(ref_res_no)%res_name, atom
                    call log_p(log_msg, me=me, level=30)
                    i_add = i_add + 1
                    i_add_heavy = i_add_heavy + 1
                    molecule%hetmol(het_no)%atm_placed(atm_no) = .true.
                end if
            end if
        end if
    end do
end do

do het_no = i_het1, i_het2
    ref_res_no = molecule%hetmol(het_no)%res_type
    do atm_no = 1, ref_res(ref_res_no)%n_atm
        if (.not. molecule%hetmol(het_no)%atm_read(atm_no) .and. &
            .not. molecule%hetmol(het_no)%atm_placed(atm_no)) then
            atom = ref_res(ref_res_no)%atom_name(atm_no)
            call place_an_atom_het(molecule, het_no, atm_no, placed)
            if (atom(1:1) == 'H') then
                write(log_msg,"(A,I4,1X,A4,1X,A4,1X,A)") '  Placed hydrogen atom ', & 
                      het_no, ref_res(ref_res_no)%res_name, atom
                call log_p(log_msg, me=me, level=30)
              
                if (placed) then
                    i_add = i_add + 1
                    molecule%hetmol(het_no)%atm_placed(atm_no) = .true.
                else
                    n_unplaced = n_unplaced + 1
                end if
            end if
        end if
    end do
end do

write(log_msg,"(A,I5)") '  Number of added heavy atoms = ', i_add_heavy
call log_p(log_msg, me=me, level=30)
write(log_msg,"(A,I5)") '  Total number of added atoms = ', i_add
call log_p(log_msg, me=me, level=30)
write(log_msg,"(A,I5)") '  n_unplaced = ', n_unplaced
call log_p(log_msg, me=me, level=30)

end subroutine place_missing_atoms_het
!-------------------------------------------------------------------------------
subroutine place_missing_atoms_lig(molecule, i_lig1, i_lig2, n_unplaced)
!-------------------------------------------------------------------------------
! Similar to place_missing_atoms, ligand molecule version
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer, intent(in) :: i_lig1, i_lig2
integer, intent(out) :: n_unplaced
integer :: i_add, atm_no, ref_res_no, i_add_heavy, lig_no
character(len=4) :: atom
logical :: placed

i_add = 0
i_add_heavy = 0
n_unplaced = 0

! place heavy atoms first
do lig_no = i_lig1, i_lig2
    ref_res_no = molecule%ligand(lig_no)%lig_type
    do atm_no = 1, ref_lig(ref_res_no)%n_atm
        if (molecule%ligand(lig_no)%atm_read(atm_no)) then
            molecule%ligand(lig_no)%atm_placed(atm_no) = .true.
        end if
    end do
    do atm_no = 1, ref_lig(ref_res_no)%n_atm
        if (.not. molecule%ligand(lig_no)%atm_read(atm_no) .and. &
            .not. molecule%ligand(lig_no)%atm_placed(atm_no)) then
            atom = ref_lig(ref_res_no)%atom_name(atm_no)
           
            if (atom(1:1) /= 'H') then
                call place_an_atom_lig(molecule, lig_no, atm_no, placed)
                if (placed) then
                    write(log_msg,"(A,I4,1X,A4,1X,A4,1X,A)") '  Placed non-hydrogen atom ', & 
                          lig_no, ref_lig(ref_res_no)%lig_name, atom
                    call log_p(log_msg, me=me, level=30)
                    i_add = i_add + 1
                    i_add_heavy = i_add_heavy + 1
                    molecule%ligand(lig_no)%atm_placed(atm_no) = .true.
                end if
            end if
        end if
    end do
end do

do lig_no = i_lig1, i_lig2
    ref_res_no = molecule%ligand(lig_no)%lig_type
    do atm_no = 1, ref_lig(ref_res_no)%n_atm
        if (.not. molecule%ligand(lig_no)%atm_read(atm_no) .and. &
            .not. molecule%ligand(lig_no)%atm_placed(atm_no)) then
            atom = ref_lig(ref_res_no)%atom_name(atm_no)
            call place_an_atom_lig(molecule, lig_no, atm_no, placed)
            if (atom(1:1) == 'H') then
                write(log_msg,"(A,I4,1X,A4,1X,A4,1X,A)") '  Placed hydrogen atom ', & 
                      lig_no, ref_lig(ref_res_no)%lig_name, atom
                call log_p(log_msg, me=me, level=30)
              
                if (placed) then
                    i_add = i_add + 1
                    molecule%ligand(lig_no)%atm_placed(atm_no) = .true.
                else
                    n_unplaced = n_unplaced + 1
                end if
            end if
        end if
    end do
end do

write(log_msg,"(A,I5)") '  Number of added heavy atoms = ', i_add_heavy
call log_p(log_msg, me=me, level=30)
write(log_msg,"(A,I5)") '  Total number of added atoms = ', i_add
call log_p(log_msg, me=me, level=30)
write(log_msg,"(A,I5)") '  n_unplaced = ', n_unplaced
call log_p(log_msg, me=me, level=30)

end subroutine place_missing_atoms_lig
!-------------------------------------------------------------------------------
subroutine place_an_atom(molecule, res_no, atom_no, placed, ierr)
!-------------------------------------------------------------------------------
! Place a missing atom (at atom_no & res_no) on molecule,
! according to hybridization state, connectivity, and so on.
!-------------------------------------------------------------------------------
integer, intent(in) :: res_no, atom_no
type(molecule_type), intent(inout) :: molecule
integer, intent(out) :: ierr
logical, intent(out) :: placed
type(ref_res_type) :: ref
integer :: i, j, k, i_bnd, i_ang, i_res, atom_idx(4), res_idx(4)
integer :: n_others, atm_in_bnd(2), idx, others(2,2), i1, i2, i_atm
integer :: ref_res_no, i_ref
logical :: NterH, bond_found, angle_found, NterHA, NterHD, NterCB, ter5H
character(len=4) :: atom(4), place_type, atom_other
character(len=6) :: atom_cls(4)
real(dp) :: r_atm(3,3), b_len, r(3,3), angle, angle_others(2)
real(dp) :: bnd_ang, tor_ang, d_ang, x(3), dr(3), z(3), p(4), q(4), s(4), U(3,3)
character(len=1) :: ter_type
character(len=6) :: atom_error_mode = ''
logical :: status

ref_res_no = molecule%residue(res_no)%res_type
ref = ref_res(ref_res_no)
atom(4) = ref%atom_name(atom_no)
ter_type = molecule%residue(res_no)%ter_type
ierr = 0

!! get three atoms (A1-A2-A3-atom4) preceeding atom(4)
! check if Nterminal H or HA
NterH = .false.
NterHA = .false.
NterHD = .false.
NterCB = .false.
ter5H = .false.
if (ter_type == 'N') then
    if (atom(4) == 'H1' .or. atom(4) == 'H2' .or. atom(4) == 'H3') then
        NterH = .true.
    end if
    if (atom(4)(1:2) == 'HA') then
        NterHA = .true.
    else if (atom(4)(1:2) == 'CB') then
        NterCB = .true.
    end if
    if (ref%res_name == 'NPRO' .and. atom(4)(1:2) == 'HD') then
        NterHD = .true.
    end if
else if (ter_type == '5') then
    if (atom(4) == 'H5T') then
        ter5H = .true.
    end if
end if

! get A1, A2, A3 and place_type
if (NterH) then
    ! if Nterminal H, return C, CA, N
    atom(1) = 'C'
    atom(2) = 'CA'
    atom(3) = 'N'
    place_type = 'sp3'
else if (NterHA) then
    ! if Nterminal HA, return +N, C, CA
    atom(1) = '+N'
    atom(2) = 'C'
    atom(3) = 'CA'
    place_type = 'sp3'
else if (NterHD) then
    ! if Nterminal HD of Proline
    atom(1) = 'CA'
    atom(2) = 'N'
    atom(3) = 'CD'
    place_type = 'sp3'
else if (NterCB) then
    ! if Nterminal CB, return +N, C, CA
    atom(1) = '+N'
    atom(2) = 'C'
    atom(3) = 'CA'
    place_type = 'sp3'
else if (ref%res_name == 'NACE' .and. (atom(4) == 'HH31' .or. atom(4) == 'HH32' .or. atom(4) == 'HH33')) then
    atom(1) = '+N'
    atom(2) = 'C'
    atom(3) = 'CH3'
    place_type = 'sp3'
else if (ref%res_name == 'PRO' .and. atom(4) == 'CB  ') then
    atom(1) = '-C'
    atom(2) = 'N'
    atom(3) = 'CA'
    place_type = 'sp3'
else if (ter5H) then
    ! if 5' H, return C4', C5', O5'
    atom(1) = "C4'"
    atom(2) = "C5'"
    atom(3) = "O5'"
    place_type = 'sp3'
else
    ! get bond idx, i_bnd
    bond_found = .false.
    do i = 1, ref%n_atm
        if (ref%atm_in_bnd(2,i) == atom_no) then
            bond_found = .true.
            i_bnd = i
            idx = ref%atm_in_bnd(1,i)
            atom(3) = ref%atom_name(idx)
            place_type = eng_para%coord_type(ref_res_eng(ref_res_no)%atm_cls(idx))
            exit
        end if
    end do
   
    if (.not. bond_found) then
        ierr = 1
        return
        !write(log_msg,"(A,I5,A5,A5)") 'Error. Parent bond not found for ', &
        !     res_no, ref%res_name, atom(4)
        !call terminate_with_error(log_msg)
    end if
    if (place_type /= 'sp2' .and. place_type /= 'sp3') then
        ierr = 2
        return
        !write(log_msg,"(A,A5,A5)") 'Error. Placement type could not be determined for ', &
        !     atom(4), place_type
        !call terminate_with_error(log_msg)
    end if
   
    ! get atoms in parent dih angle of this atom
    ! get parent angle of this bond
    angle_found = .false.
    do i = 1, ref%n_atm
        if (ref%bnd_in_t_ang(3,i) == i_bnd) then
            angle_found = .true.
            i_ang = i
            exit
        end if
    end do
    if (.not. angle_found) then
        !write(log_msg,"(A,A5)") 'Error. First parent angle not found for ', atom
        !call terminate_with_error(log_msg)
        ierr = 1
        return
    end if
    atom(1:2) = ref%atom_name(ref%atm_in_t_ang(1:2,i_ang))
end if

do i = 1, 4
    if (atom(i)(1:1) == '-') then
        if (molecule%cyclic .and. res_no == 1) then
            i_res = molecule%n_res
        else
            i_res = res_no - 1
        end if
        atom(i) = atom(i)(2:)
        if (ref_res(molecule%residue(i_res)%res_type)%res_name == 'NACE' .and. atom(i) == 'CA') then
            atom(i) = 'CH3'
        end if
    else if (atom(i)(1:1) == '+') then
        if (molecule%cyclic .and. res_no == molecule%n_res) then
            i_res = 1
        else
            i_res = res_no + 1
        end if
        atom(i) = atom(i)(2:)
    else
        i_res = res_no
    end if
   
    ! get atom cls
    i_ref = molecule%residue(i_res)%res_type
    call find_atom_idx(i_res, i_ref, atom(i), idx, atom_error_mode)
    atom_cls(i) = eng_para%atom_cls(ref_res_eng(i_ref)%atm_cls(idx))
    res_idx(i) = i_res
    atom_idx(i) = idx
   
    ! get coord of the bonded atoms
    if (i < 4) then
        if (molecule%residue(i_res)%atm_read(idx) .or. molecule%residue(i_res)%atm_placed(idx)) then
            r_atm(:,i) = molecule%residue(i_res)%R(:,idx)
        else
            call log_p('Placing this way failed. Will place in another way.', me=me, level=20)
            placed = .false.
            return
        end if
    end if
end do

! get b_len for placement
call find_bnd_prm(atom_cls(3:4), i_bnd, status)
b_len = eng_para%bnd_para(2,i_bnd)

! get bnd_ang for placement
call find_ang_prm(atom_cls(2:4), 2, i_ang, status)
bnd_ang = eng_para%ang_para(2,1,i_ang)

if (.not. status) then
    write(log_msg, "(A,I4,A5)") 'ERROR: failed in adding', res_no, atom(4)
    call terminate_with_error(log_msg)
end if

! get atoms connected to atom(3), except for atom(2) and atom(4)
n_others = 0
do i = max(1, res_no - 1), min(molecule%n_res, res_no + 1)
    i_ref = molecule%residue(i)%res_type
    do j = 1, ref_res_eng(i_ref)%n_bnd_E
        atm_in_bnd(1:2) = ref_res_eng(i_ref)%atm_in_bnd_E(1:2,j)
        if (atm_in_bnd(1) == atom_idx(3) .and. i == res_idx(3)) then
            k = 2
        else if (atm_in_bnd(2) == atom_idx(3) .and. i == res_idx(3)) then
            k = 1
        else
            k = 0
        end if
      
        if (k > 0) then
            idx = atm_in_bnd(k)
            if (idx <= 0 .and. i == 1) then
                k = 0
            end if
        end if
      
        if (k > 0) then
            if (idx == 0 .and. i > 1) then ! -C
                if (molecule%cyclic .and. res_no == 1) then
                    i_res = molecule%n_res
                else
                    i_res = i - 1
                end if
                atom_other = 'C'
            else if (idx == -1) then ! -CA
                if (molecule%cyclic .and. res_no == 1) then
                    i_res = molecule%n_res
                else
                    i_res = i - 1
                end if
                if (ref_res(molecule%residue(i_res)%res_type)%res_name == 'NACE') then
                    atom_other = 'CH3'
                else
                    atom_other = 'CA'
                end if
            else if (idx == ref_res(i_ref)%n_atm + 1) then ! +N
                if (molecule%cyclic .and. res_no == molecule%n_res) then
                    i_res = 1
                else
                    i_res = i + 1
                end if
                atom_other = 'N'
            else
                i_res = i
                i_atm = idx
                atom_other = ''
            end if
            if (atom_other /= '') then
                call find_atom_idx(i_res, molecule%residue(i_res)%res_type, atom_other, i_atm, atom_error_mode)
            end if
            if ((i_atm /= atom_idx(2) .or. i_res /= res_idx(2)) .and. & 
                (i_atm /= atom_idx(4) .or. i_res /= res_idx(4)) .and. &
                (i_atm /= atom_idx(1) .or. i_res /= res_idx(1))) then
                if (molecule%residue(i_res)%atm_read(i_atm) .or. molecule%residue(i_res)%atm_placed(i_atm)) then
                    if (n_others == 1 .and. others(1,1) == i_res .and. others(2,1) == i_atm) cycle ! for backbone oxygen
                    !
                    n_others = n_others + 1
                    if (n_others > 2) then
                        ierr = 3
                        return
                    !   call terminate_with_error('Error. N_others > 2 on placing atom.')
                    end if
                    others(1:2,n_others) = (/ i_res, i_atm /)
                    ! torsion angle of other bnds connected to the same atom(3)
                    r(:,1) = r_atm(:,2) - r_atm(:,1)
                    r(:,2) = r_atm(:,3) - r_atm(:,2)
                    r(:,3) = molecule%residue(i_res)%R(:,i_atm) - r_atm(:,3)
                    call calc_tor_angle(r(:,1), r(:,2), r(:,3), angle)
                    ! place angle in [0,two_pi]
                    if (angle < 0.0d0) then
                        angle = angle + two_pi
                    end if
                    angle_others(n_others) = angle
                end if
            end if
        end if
    end do
end do

if (n_others == 0) then
!     tor_ang = pi ! trans
    tor_ang = ref%t_ang0(atom_no) ! modeled by topology file
else if (n_others == 1) then
    if (place_type == 'sp3') then
        if (atom(4) == 'CB  ' .and. (.not. NterCB)) then
            tor_ang = angle_others(1) + deg2rad*240.0d0
        else
            tor_ang = angle_others(1) + deg2rad*120.0d0
        end if
    else 
        tor_ang = angle_others(1) + pi
    end if
else if (n_others == 2) then
     if (place_type == 'sp2') then
         ierr = 3
         return
         !call terminate_with_error('Error. Inconsistency found on placing atom.')
     end if
   
    ! find smaller angle
    if (angle_others(1) < angle_others(2)) then
        i1 = 1; i2 = 2
    else
        i1 = 2; i2 = 1
    end if
    d_ang = angle_others(i2) - angle_others(i1) ! always positive
    if (d_ang < pi) then
        tor_ang = angle_others(i2) + 0.5d0*(two_pi - d_ang)
    else
        tor_ang = angle_others(i2) - 0.5d0*d_ang
    end if
end if

! define x, y, and z axes
x(:) = r_atm(:,3) - r_atm(:,2)
x(:) = x(:)/sqrt(dot_product(x,x))
dr(:) = r_atm(:,1) - r_atm(:,2)
call cross(x, dr, z)
z(:) = z(:)/sqrt(dot_product(z,z))

! first rotate about z by (pi - bnd_ang)
call quaternion(z, (pi-bnd_ang), p)
! next rotate about x by tor_ang
tor_ang = bound_ang(tor_ang)
call quaternion(x, tor_ang, q)
call q_product(3, q, p, s)
call rotation_matrix(3, s, U)

molecule%residue(res_idx(4))%R(:,atom_idx(4)) = r_atm(:,3) + b_len*matmul(U, x)
placed = .true.

end subroutine place_an_atom
!-------------------------------------------------------------------------------
subroutine place_an_atom_het(molecule, het_no, atom_no, placed)
!-------------------------------------------------------------------------------
! Similar to place_an_atom_het, hetero molecule version
!-------------------------------------------------------------------------------
integer, intent(in) :: het_no, atom_no
type(molecule_type), intent(inout) :: molecule
logical, intent(out) :: placed
type(ref_res_type) :: ref
integer :: i, j, k, i_bnd, i_ang, atom_idx(4)
integer :: n_others, atm_in_bnd(2), idx, others(2), i1, i2, i_atm
integer :: ref_res_no
logical :: bond_found, angle_found, atom_found, atom1_found, atom2_found
character(len=4) :: atom(4), place_type, atom_con1, atom_con2, atom_con3
character(len=6) :: atom_cls(4)
real(dp) :: r_atm(3,3), b_len, r(3,3), angle, angle_others(2)
real(dp) :: bnd_ang, tor_ang, d_ang, x(3), dr(3), z(3), p(4), q(4), s(4), U(3,3)
character(len=6) :: atom_error_mode = ''
integer :: i_b1, i_b2, i_b3, atm1(2), atm2(2), atm3(2), atm_con1, atm_con2, atm_con3
logical :: status

ref_res_no = molecule%hetmol(het_no)%res_type
ref = ref_res(ref_res_no)
atom(4) = ref%atom_name(atom_no)

! Need atom(1:3) and place_type  
! get bond idx, i_bnd
bond_found = .false.
do i = 1, ref%n_atm
    if (ref%atm_in_bnd(2,i) == atom_no .and. ref%atm_in_bnd(1,i) > 0) then
        bond_found = .true.
        i_bnd = i
        idx = ref%atm_in_bnd(1,i)
        atom(3) = ref%atom_name(idx)
        place_type = eng_para%coord_type(ref_res_eng(ref_res_no)%atm_cls(idx))
        exit
    end if
end do

if (bond_found) then
    if (place_type /= 'sp2' .and. place_type /= 'sp3') then
        write(log_msg,"(A,A5,A4)") 'Error. Placement type could not be determined for ', &
              atom(4), place_type
        call terminate_with_error(log_msg)
    end if
    
    ! get atoms in parent dih angle of this atom
    ! get parent angle of this bond
    angle_found = .false.
    do i = 1, ref%n_atm
        if (ref%bnd_in_t_ang(3,i) == i_bnd .and. ref%atm_in_bnd(1,ref%bnd_in_t_ang(1,i)) > 0) then
            angle_found = .true.
            i_ang = i
            exit
        end if
    end do
   
    if (angle_found) then
        atom(1:2) = ref%atom_name(ref%atm_in_t_ang(1:2,i_ang))
    end if
end if

write(log_msg,"(A,A5,L2,A5,L2,A5,A5)") '  Placing', atom(4), &
      bond_found, atom(3), angle_found, atom(2), atom(1)
call log_p(log_msg, level=40)

if (.not.bond_found .or. .not.angle_found) then
    !! try to find alternative bond
    ! if found, angle_found = .true.
    atom_found = .false.
    atom1_found = .false.
    atom2_found = .false.
   
    do i_b1 = 1, ref%n_atm
        atm1(1:2) = ref%atm_in_bnd(1:2, i_b1)
        atm_con1 = 0
        if (atm1(1) == atom_no) then
            atm_con1 = atm1(2)
        else if (atm1(2) == atom_no) then
            atm_con1 = atm1(1)
        end if
      
        ! first connected atom
        if (atm_con1 > 0) then
            atom_con1 = ref%atom_name(atm_con1)
            ! require heavy atom
            if (atom_con1(1:1) /= 'H') then
                atom1_found = .true.
                atom(3) = ref%atom_name(atm_con1)
                do i_b2 = 1, ref%n_atm
                    atm2(1:2) = ref%atm_in_bnd(1:2, i_b2)
                    atm_con2 = 0
                    if (atm2(1) == atm_con1) then
                        atm_con2 = atm2(2)
                    else if (atm2(2) == atm_con1) then
                        atm_con2 = atm2(1)
                    end if
                 
                    ! second connected atom
                    if (atm_con2 > 0 .and. atm_con2 /= atom_no) then
                        atom_con2 = ref%atom_name(atm_con2)
                        ! require heavy atom
                        if (atom_con2(1:1) /= 'H') then
                            atom2_found = .true.
                            atom(2) = ref%atom_name(atm_con2)
                            do i_b3 = 1, ref%n_atm
                                atm3(1:2) = ref%atm_in_bnd(1:2, i_b3)
                                atm_con3 = 0
                                if (atm3(1) == atm_con2) then
                                    atm_con3 = atm3(2)
                                else if (atm3(2) == atm_con2) then
                                    atm_con3 = atm3(1)
                                end if
                              
                                if (atm_con3 > 0 .and. atm_con3 /= atm_con1) then
                                    atom_con3 = ref%atom_name(atm_con3)
                                    ! require heavy atom
                                    if (atom_con3(1:1) /= 'H') then
                                        atom(1) = ref%atom_name(atm_con3)
                                        atom_found = .true.
                                        goto 111
                                    end if
                                end if
                                
                            end do
                        end if
                        
                    end if
                end do
          
            end if
        end if
    end do
   
111 continue

else
    atom_found = .true.
end if

if (.not.atom_found) then
    atom(1) = 'CA'
    if (.not.atom2_found) then
        atom(2) = 'CA'
    end if
    if (.not.atom1_found) then
        atom(3) = 'CA'
    end if
   
    ! find atom_cls, atom_idx, r_atm
    placed = .false.
    return
   
else
    do i = 1, 4
        ! get atom cls
        call find_atom_idx(het_no, ref_res_no, atom(i), idx, atom_error_mode)
        atom_cls(i) = eng_para%atom_cls(ref_res_eng(ref_res_no)%atm_cls(idx))
        atom_idx(i) = idx
       
        ! get coord of the bonded atoms
        if (i < 4) then
            if (molecule%hetmol(het_no)%atm_read(idx) .or. molecule%hetmol(het_no)%atm_placed(idx)) then
                r_atm(:,i) = molecule%hetmol(het_no)%R(:,idx)
            else
                placed = .false.
                return
            end if
        end if
    end do
end if

! get b_len for placement
call find_bnd_prm(atom_cls(3:4), i_bnd, status)
b_len = eng_para%bnd_para(2,i_bnd)

! get bnd_ang for placement
call find_ang_prm(atom_cls(2:4), 2, i_ang, status)
bnd_ang = eng_para%ang_para(2,1,i_ang)

if (.not. status) then
    write(log_msg, "(A,I4,A5)") 'ERROR: failed in adding', het_no, atom(4)
    call terminate_with_error(log_msg)
end if

! get atoms connected to atom(3), except for atom(2) and atom(4)
n_others = 0
do j = 1, ref_res_eng(ref_res_no)%n_bnd_E
    atm_in_bnd(1:2) = ref_res_eng(ref_res_no)%atm_in_bnd_E(1:2,j)
    if (atm_in_bnd(1) == atom_idx(3)) then
        k = 2
    else if (atm_in_bnd(2) == atom_idx(3)) then
        k = 1
    else
        k = 0
    end if
   
    if (k > 0) then
        idx = atm_in_bnd(k)
        i_atm = idx
        if ((i_atm /= atom_idx(2)) .and. (i_atm /= atom_idx(4)) .and. (i_atm /= atom_idx(1)) ) then
            if (molecule%hetmol(het_no)%atm_read(i_atm) .or. molecule%hetmol(het_no)%atm_placed(i_atm)) then
                n_others = n_others + 1
                if (n_others > 2) then
                    call terminate_with_error('Error. N_others > 2 on placing atom.')
                end if
                others(n_others) = i_atm
                ! torsion angle of other bnds connected to the same atom(3)
                r(:,1) = r_atm(:,2) - r_atm(:,1)
                r(:,2) = r_atm(:,3) - r_atm(:,2)
                r(:,3) = molecule%hetmol(het_no)%R(:,i_atm) - r_atm(:,3)
                call calc_tor_angle(r(:,1), r(:,2), r(:,3), angle)
                ! place angle in [0,two_pi]
                if (angle < 0.0d0) then
                    angle = angle + two_pi
                end if
                angle_others(n_others) = angle
            end if
        end if
    end if
end do

if (n_others == 0) then
    tor_ang = pi ! trans
else if (n_others == 1) then
    if (place_type == 'sp3') then
        tor_ang = angle_others(1) + deg2rad*120.0d0
    else 
        tor_ang = angle_others(1) + pi
    end if
else if (n_others == 2) then
    if (place_type == 'sp2') then
        call terminate_with_error('Error. Inconsistency found on placing atom.')
    end if
    ! find smaller angle
    if (angle_others(1) < angle_others(2)) then
        i1 = 1; i2 = 2
    else
        i1 = 2; i2 = 1
    end if
    d_ang = angle_others(i2) - angle_others(i1) ! always positive
    if (d_ang < pi) then
        tor_ang = angle_others(i2) + 0.5d0*(two_pi - d_ang)
    else
        tor_ang = angle_others(i2) - 0.5d0*d_ang
    end if
end if

! define x, y, and z axes
x(:) = r_atm(:,3) - r_atm(:,2)
x(:) = x(:)/sqrt(dot_product(x,x))
dr(:) = r_atm(:,1) - r_atm(:,2)
call cross(x, dr, z)
z(:) = z(:)/sqrt(dot_product(z,z))

! first rotate about z by (pi - bnd_ang)
call quaternion(z, (pi-bnd_ang), p)
! next rotate about x by tor_ang
call quaternion(x, tor_ang, q)
call q_product(3, q, p, s)
call rotation_matrix(3, s, U)

molecule%hetmol(het_no)%R(:,atom_idx(4)) = r_atm(:,3) + b_len*matmul(U, x)
placed = .true.

end subroutine place_an_atom_het
!-------------------------------------------------------------------------------
subroutine place_an_atom_lig(molecule, lig_no, atom_no, placed)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer, intent(in) :: lig_no, atom_no
logical, intent(out) :: placed

type(ref_lig_type) :: ref
integer, parameter :: max_candidate = 10
integer :: ref_res_no, i, ia, n_tor
integer :: atom(4, max_candidate)
character(len=4) :: atom_name(4, max_candidate)

integer :: n_placed
logical :: info_found
real(dp) :: Rs(3, max_candidate)
real(dp) :: b_len, b_ang, t_ang, R_atm(3,3), dR(3,3)
real(dp) :: p(4), q(4), s(4), U(3,3)

placed = .false.
ref_res_no = molecule%ligand(lig_no)%lig_type
ref = ref_lig(ref_res_no)

n_tor = 0
do ia = 1, ref%n_dih
    if (ref%dih(2,ia) == atom_no) then
        if (.not. molecule%ligand(lig_no)%atm_placed(ref%dih(3,ia))) cycle
        if (.not. molecule%ligand(lig_no)%atm_placed(ref%dih(4,ia))) cycle
        if (.not. molecule%ligand(lig_no)%atm_placed(ref%dih(5,ia))) cycle
        n_tor = n_tor + 1
        atom(4,n_tor) = ref%dih(2,ia)
        atom(3,n_tor) = ref%dih(3,ia)
        atom(2,n_tor) = ref%dih(4,ia)
        atom(1,n_tor) = ref%dih(5,ia)
        atom_name(4,n_tor) = ref%atom_name(ref%dih(2,ia))
        atom_name(3,n_tor) = ref%atom_name(ref%dih(3,ia))
        atom_name(2,n_tor) = ref%atom_name(ref%dih(4,ia))
        atom_name(1,n_tor) = ref%atom_name(ref%dih(5,ia))
    else if (ref%dih(5,ia) == atom_no) then
        if (.not. molecule%ligand(lig_no)%atm_placed(ref%dih(2,ia))) cycle
        if (.not. molecule%ligand(lig_no)%atm_placed(ref%dih(3,ia))) cycle
        if (.not. molecule%ligand(lig_no)%atm_placed(ref%dih(4,ia))) cycle
        n_tor = n_tor + 1
        atom(4,n_tor) = ref%dih(5,ia)
        atom(3,n_tor) = ref%dih(4,ia)
        atom(2,n_tor) = ref%dih(3,ia)
        atom(1,n_tor) = ref%dih(2,ia)
        atom_name(4,n_tor) = ref%atom_name(ref%dih(5,ia))
        atom_name(3,n_tor) = ref%atom_name(ref%dih(4,ia))
        atom_name(2,n_tor) = ref%atom_name(ref%dih(3,ia))
        atom_name(1,n_tor) = ref%atom_name(ref%dih(2,ia))
    end if
    if (n_tor >= max_candidate) exit
end do

n_placed = 0
Rs(:,:) = 0.0d0
do ia = 1, n_tor
    R_atm(:,1) = molecule%ligand(lig_no)%R(:, atom(1,ia))
    R_atm(:,2) = molecule%ligand(lig_no)%R(:, atom(2,ia))
    R_atm(:,3) = molecule%ligand(lig_no)%R(:, atom(3,ia))
   
    info_found = .false.
    do i = 1, ref%n_bnd
        if ((ref%bnd(2,i) == atom(3,ia) .and. ref%bnd(3,i) == atom(4,ia)) .or. &
            (ref%bnd(2,i) == atom(4,ia) .and. ref%bnd(3,i) == atom(3,ia))) then
             b_len = ref%b_len0(i)
             info_found = .true.
             exit
        end if
    end do
    if (.not. info_found) cycle
   
    info_found = .false.
    do i = 1, ref%n_ang
        if ((ref%ang(2,i) == atom(2,ia) .and. ref%ang(3,i) == atom(3,ia) .and. &
             ref%ang(4,i) == atom(4,ia)) .or. &
            (ref%ang(2,i) == atom(4,ia) .and. ref%ang(3,i) == atom(3,ia) .and. &
             ref%ang(4,i) == atom(2,ia))) then
             b_ang = ref%b_ang0(i)
             info_found = .true.
         end if
    end do
    if (.not. info_found) cycle
   
    info_found = .false.
    do i = 1, ref%n_dih
        if ((ref%dih(2,i) == atom(1,ia) .and. ref%dih(3,i) == atom(2,ia) .and. &
             ref%dih(4,i) == atom(3,ia) .and. ref%dih(5,i) == atom(4,ia)) .or. &
            (ref%dih(2,i) == atom(4,ia) .and. ref%dih(3,i) == atom(3,ia) .and. &
             ref%dih(4,i) == atom(2,ia) .and. ref%dih(5,i) == atom(1,ia))) then
             t_ang = ref%d_ang0(i)
             info_found = .true.
         end if
    end do
    if (.not. info_found) cycle
   
    dR(:,1) = R_atm(:,3) - R_atm(:,2)
    call v_norm(dR(:,1))
    dR(:,2) = R_atm(:,1) - R_atm(:,2)
    call v_norm(dR(:,2))
    call cross(dR(:,1),dR(:,2),dR(:,3))
    call v_norm(dR(:,3))
   
    call quaternion(dR(:,3), (pi-b_ang), p)
    call quaternion(dR(:,1), t_ang, q)
    call q_product(3, q, p, s)
   
    call rotation_matrix(3, s, U)
   
    n_placed = n_placed + 1
    Rs(:,n_placed) = R_atm(:,3) + b_len*matmul(U, dR(:,1))
end do

placed = .true.
molecule%ligand(lig_no)%atm_placed(atom_no) = .true.
molecule%ligand(lig_no)%R(1, atom_no) = sum(Rs(1,:))/float(n_placed)
molecule%ligand(lig_no)%R(2, atom_no) = sum(Rs(2,:))/float(n_placed)
molecule%ligand(lig_no)%R(3, atom_no) = sum(Rs(3,:))/float(n_placed)

end subroutine place_an_atom_lig
!-------------------------------------------------------------------------------
subroutine setup_link_btw_res(molecule)
!-------------------------------------------------------------------------------
! Set link_type variable inside molecule
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer :: i_res, i_res_prev, i_ref_res_prev, k_atm, i_atm, atm_idx, link_res_no
character(len=4) :: atom, atom_ref, prev_res_name

do i_res = 2, molecule%n_res
    ! Find indices for the three atoms in the previous residue
    ! that are directly linked to the current residue.
    ! Here we consider the case of regular peptide bonds only,
    ! and the three atoms are  (-N, -CA, -C)
    i_res_prev = i_res - 1
    i_ref_res_prev = molecule%residue(i_res_prev)%res_type
    prev_res_name = ref_res(i_ref_res_prev)%res_name
   
    do k_atm = 1, 3
        link_res_no = i_res - 1
        if (mol_type == 'protein') then
            if (k_atm == 1) then
                if (prev_res_name == 'NACE' .and. trim(top_type) == 'allh_ch22') then
                    atom = 'HH31'
                else if (prev_res_name == 'NACE' .and. trim(top_type) == 'polarh') then
                    atom = 'CH3' 
                else if (prev_res_name == 'CNHE') then
                    atom = 'CA'
                    link_res_no = i_res - 2
                else
                    atom = 'N'
                end if
            else if (k_atm == 2) then
                if (prev_res_name == 'NACE') then
                    atom = 'CH3'
                else if (prev_res_name == 'CNHE') then
                    atom = 'C'
                    link_res_no = i_res - 2
                else
                    atom = 'CA'
                end if
            else
                if (prev_res_name == 'CNHE') then
                    atom = 'N'
                else
                    atom = 'C'
                end if
            end if
        end if
      
        i_ref_res_prev = molecule%residue(link_res_no)%res_type
        do i_atm = 1, ref_res(i_ref_res_prev)%n_atm
            atom_ref = ref_res(i_ref_res_prev)%atom_name(i_atm)
            if (atom == atom_ref) then
                atm_idx = i_atm
                exit
            end if
        end do
      
        molecule%residue(i_res)%link_atm_no(k_atm) = atm_idx 
        molecule%residue(i_res)%link_res_no(k_atm) = link_res_no
    end do
end do

end subroutine setup_link_btw_res
!-------------------------------------------------------------------------------
subroutine apply_fix_type(molecule)
!-------------------------------------------------------------------------------
! Using input ULR/USC information & fix_type, update ULR/USC information for
! sampling and fix atoms.
! This subroutine should be called before [initialize_energy] to activate the
! degree of freedom variables.
! 
! e.g.
!
!   fix_type   | ULR(input) | USC(input) |  final USC |  atom_fixed  
!  ------------------------------------------------------------------
!     all      |     x      |     x      |      x     | all stdres 
!              |     o      |     o      |  as input  | all stdres except ULR/USC
!              |     x      |     o      |  as input  | all stdres except USC
!  ------------------------------------------------------------------
!     none     |     x      |     x      | all stdres | none
!              |     o      |     o      |  as input  | none
!              |     x      |     o      |  as input  | none
!  ------------------------------------------------------------------
!    backbone  |     x      |     x      | all stdres | backbone of all stdres
!              |     x      |     o      |  as input  | backbone of all stdres
!
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer :: i_ulr, i_res, i_atm, ref_res_no
logical, allocatable, save :: is_ulr(:)    ! To store ULR residues.
                                           ! By using 'save', it retains the value when
                                           ! execution of subroutine is ended.

if (.not. are_defined) then
    call reinitialize_usc(tn%stdres)
    allocate(is_ulr(tn%stdres))
    is_ulr(:) = .false.
    are_defined  = .true.
    ! If USC info. is not given as an input, set those based on fix_type.
    if (n_usc == 0) then
        if ((fix_type == FIX_NONE) .or. (fix_type == FIX_BB)) then
            is_usc(:) = .true.
            n_usc = molecule%n_res
        end if
    end if
    !
    if (symmetric) then
        ! is_usc, ULR is only set to the first symmetric unit.
        call set_USC_for_symm(is_usc, n_usc)
        call set_ULR_for_symm(ULR)
    end if
    ! setup is_ulr to assign atom_fixed easily.
    if (n_ulr > 0) then
        do i_ulr = 1, n_ulr
            is_ulr(ULR(i_ulr)%resrange(1):ULR(i_ulr)%resrange(2)) = .true.
        end do
    end if

    ! For log messages
    if (fix_type == FIX_ALL) then
        call log_p('- Fix type "all" activated: All atoms except ULR/USC are fixed.',\
                   me=me, level=20)
    else if (fix_type == FIX_NONE) then
        call log_p('- Fix type "none" activated: All atoms are free to move.',\
                   me=me, level=20)
    else if (fix_type == FIX_BB) then
        call log_p('- Fix type "backbone" activated: All backbone atoms except ULR are fixed.',\
                   me=me, level=20)
    end if
end if

! Fix atoms based on ULR/USC info and fix_type.
! For all residues, atom_fixed has already initialized to false.
if (fix_type == FIX_ALL) then ! fix all atoms except ULR/USC
    do i_res = 1, molecule%n_res
        if (is_ulr(i_res)) cycle 
        if (is_usc(i_res)) then  ! For USC, fix backbone only.
            ref_res_no = molecule%residue(i_res)%res_type
            do i_atm = 1, molecule%residue(i_res)%n_atm
                if (ref_res(ref_res_no)%is_sc_atom(i_atm)) cycle
                molecule%residue(i_res)%atom_fixed(i_atm) = .true.
            end do
        else                     
            molecule%residue(i_res)%atom_fixed(:) = .true.
        end if
    end do

else if (fix_type == FIX_BB) then ! fix all backbone except ULR
    do i_res = 1, molecule%n_res
        ref_res_no = molecule%residue(i_res)%res_type
        do i_atm = 1, molecule%residue(i_res)%n_atm
            if (ref_res(ref_res_no)%is_sc_atom(i_atm)) cycle
            molecule%residue(i_res)%atom_fixed(i_atm) = .true.
        end do
    end do
end if

call apply_fix_type_symm(molecule)

end subroutine apply_fix_type
!-------------------------------------------------------------------------------
subroutine setup_tot_num(molecule)
!-------------------------------------------------------------------------------
! Fill tot_num_type, tn (except dof information).
! Total number of residues, standard residue(amino acid + nucleic acid), hetero
! molecule, and ligand are counted.
! Total number of atoms, atoms in standard residues, atoms in hetmol, and atoms
! in ligands are counted.
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule
integer :: i_res, ref_res_no

! count the number of atoms
tn%stdatm = 0
tn%hetatm = 0
tn%ligatm = 0
do i_res = 1, molecule%n_res
    ref_res_no = molecule%residue(i_res)%res_type
    tn%stdatm = tn%stdatm + ref_res(ref_res_no)%n_atm
end do
do i_res = 1, molecule%n_het
    ref_res_no = molecule%hetmol(i_res)%res_type
    tn%hetatm = tn%hetatm + ref_res(ref_res_no)%n_atm
end do
tn%nonligatm = tn%stdatm + tn%hetatm
do i_res = 1, molecule%n_lig
    ref_res_no = molecule%ligand(i_res)%lig_type
    tn%ligatm = tn%ligatm + ref_lig(ref_res_no)%n_atm
end do
tn%atom = tn%nonligatm + tn%ligatm

tn%recatm = 0
do i_res = 1, tn%recres ! count the num of atoms in receptor protein
    ref_res_no = molecule%residue(i_res)%res_type
    tn%recatm = tn%recatm + ref_res(ref_res_no)%n_atm
enddo

! count the number of residues
tn%residue = molecule%n_res + molecule%n_het + molecule%n_lig
tn%stdres = molecule%n_res
tn%hetmol = molecule%n_het
tn%nonlig = molecule%n_res + molecule%n_het
tn%ligand = molecule%n_lig

end subroutine setup_tot_num
!-------------------------------------------------------------------------------
END MODULE SETUP_MOLECULE
!-------------------------------------------------------------------------------
