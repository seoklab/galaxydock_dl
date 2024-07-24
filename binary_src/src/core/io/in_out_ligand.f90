!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/io/in_out_ligand.f90
!
! Description:
!   This module contains subroutines to read/write mol2 file, and setup geometry
!   of ligand.
!-------------------------------------------------------------------------------
MODULE IN_OUT_LIGAND
!-------------------------------------------------------------------------------
use globals
use in_out_vars
use logger, only: log_p, terminate_with_error
use string, only: parse_string

implicit none
private

integer, parameter :: max_atoms_in_frag  = 50
integer, parameter :: max_arms_in_frag   = 10
integer, parameter :: max_connected_atom = 6
!-------------------------------------------------------------------------------
! derived type to describe fragments in ligand.
!-------------------------------------------------------------------------------
type lig_frag_type
!-------------------------------------------------------------------------------
integer, dimension(3, max_arms_in_frag) :: arm     ! Index for rotatable bond connecting
                                                   ! between two fragments 
                                                   ! #1 - bond index
                                                   ! #2 - one end of arm that is in the
                                                   ! fragment
                                                   ! #3 - the other end of arm that is in
                                                   ! other fragment
integer, dimension(max_atoms_in_frag) :: atm_list  ! index of atoms in fragment
logical, dimension(max_arms_in_frag)  :: arm_used  ! flag to define branches 
integer :: n_atm, n_arm
!-------------------------------------------------------------------------------
end type lig_frag_type
!-------------------------------------------------------------------------------

public :: read_mol2
public :: check_n_lig
public :: read_multiple_mol2
public :: merge_nonpolar_hydrogens
public :: write_mol2
public :: write_multiple_mol2
public :: build_tree
public :: find_lig_bnd_idx
public :: get_ligand_center

CONTAINS
!===============================================================================
! Subroutines related to read/write mol2 file.
!===============================================================================
subroutine read_mol2(mol2_file, ligand)
!-------------------------------------------------------------------------------
! This subroutine reads coordinates, atom types, charges, and bond information
! from mol2 file. Further treatment such as rigid fragment detection would be 
! covered in other subroutines.
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(inout) :: ligand
character(len=len_fname), intent(in) :: mol2_file
! related to file I/O
integer :: i_unit, ioerror
! related to read lines
integer :: i_line, n_word
character(len=len_fname) :: word(25), line
! related to information in mol2 file
logical :: molecule_info, atm_info, bnd_info

! Open file
i_unit = 31
open(file=trim(mol2_file), unit=i_unit, action='read', status='old', iostat=ioerror)
if (ioerror /= 0) then
    call terminate_with_error("Error: MOL2 file does not exist!")
end if

! Initialize variables
molecule_info = .false.
atm_info = .false.
bnd_info = .false.
i_line = 0

! Read file
do
    read(i_unit, '(a80)', iostat=ioerror) line
    if (ioerror /= 0) exit ! End of file
    if (line(1:1) == '#') cycle ! Skip unnecessary information

    ! Read molecule information
    if (line(1:17) == '@<TRIPOS>MOLECULE') then
        if (bnd_info) exit
        molecule_info = .true.
    end if

    if (molecule_info) then
        i_line = i_line + 1
    end if
    
    if (molecule_info .and. i_line == 3) then ! information aboun No. of atoms & bonds
        call parse_string(line, n_word, word)
        read(word(1), '(i3)') ligand%n_atm
        if (ligand%n_atm > max_lig_atom) then
            write(log_msg,'(A,I3,A)') 'Ligand atom number exceeds ',&
                                      max_lig_atom, '. Terminate.'
            call terminate_with_error(log_msg)
        end if
        read(word(2), '(i3)') ligand%n_bnd
        molecule_info = .false.
        i_line = 0
    end if

    ! Read atom information
    if (line(1:13) == '@<TRIPOS>ATOM') then
        atm_info = .true.
    end if
    
    if (atm_info .and. line(1:1) /= '@') then
        call parse_string(line, n_word, word)
        i_line = i_line + 1
        if (i_line == 1) then
            read(word(8)(1:3), '(a3)') ligand%lig_name
        end if
        read(word(2), '(a4)')   ligand%atom_name(i_line)
        read(word(3), '(f8.3)') ligand%R(1,i_line)
        read(word(4), '(f8.3)') ligand%R(2,i_line)
        read(word(5), '(f8.3)') ligand%R(3,i_line)
        read(word(6), '(a6)')   ligand%mol2_type(i_line)
        read(word(9), '(f8.3)') ligand%charge(i_line)
    end if

    ! Read bond information
    if (line(1:13) == '@<TRIPOS>BOND') then
        bnd_info = .true.
        atm_info = .false.
        i_line = 0
    end if
    
    if (bnd_info .and. line(1:1) /= '@') then
        call parse_string(line, n_word, word)
        i_line = i_line + 1
        read(word(1), '(i3)') ligand%bnd(1,i_line)
        read(word(2), '(i3)') ligand%bnd(2,i_line)
        read(word(3), '(i3)') ligand%bnd(3,i_line)
        if(n_word == 4) then
            ligand%bnd_type(i_line) = word(4)
        end if
    end if

    if (line(1:21) == '@<TRIPOS>SUBSTRUCTURE') then
        bnd_info = .false.
        exit
    end if
end do

close(i_unit)

end subroutine read_mol2
!-------------------------------------------------------------------------------
subroutine read_multiple_mol2(mol2_file, ligand_s)
!-------------------------------------------------------------------------------
! Read mol2 file which contains more than one ligand conformation.
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: mol2_file
type(ref_lig_type) :: ligand_s(:)
! related file I/O
integer :: i_unit, ioerr
! related to read line
integer :: i_line, n_word
character(len=len_fname) :: line, word(25)
! related to informatino in mol2 file
integer :: i_lig
logical :: mol_info, atm_info, bnd_info

i_unit = 31
open(unit=i_unit, file=trim(mol2_file), action='read', status='old', iostat=ioerr)
if (ioerr/=0) stop

! initialize_variables
mol_info = .false.
atm_info = .false.
bnd_info = .false.
i_lig = 0
i_line = 0

! Read file
do
    read(i_unit, '(a120)', iostat=ioerr) line
    if (ioerr/=0) exit ! end of file
    if (line(1:1) == '!' .or. line(1:1) == '#') cycle

    ! Read molecule info.
    if (line(1:17) == '@<TRIPOS>MOLECULE') then
        mol_info = .true.
        atm_info = .false.
        bnd_info = .false.
        i_lig = i_lig + 1
        i_line = 0
    end if

    if (mol_info) i_line = i_line + 1
    
    if (mol_info .and. i_line==3) then
        call parse_string(line, n_word, word)
        read(word(1),'(i3)') ligand_s(i_lig)%n_atm
        if (ligand_s(i_lig)%n_atm > max_lig_atom) then
            write(log_msg,'(A,I3,A)') 'Ligand atom number exceeds ',&
                                      max_lig_atom, '. Terminate.'
            call terminate_with_error(log_msg)
        end if
        read(word(2), '(i3)') ligand_s(i_lig)%n_bnd
        mol_info = .false.
        i_line = 0
    end if

    ! Read atom information
    if (line(1:13) == '@<TRIPOS>ATOM') then
        atm_info = .true.
        mol_info = .false.
        bnd_info = .false.
        i_line = 0
    end if

    if (atm_info .and. line(1:1) /= '@') then
        call parse_String(line, n_word, word)
        i_line = i_line + 1
        if (i_line == 1) read(word(8)(1:3), '(a3)') ligand_s(i_lig)%lig_name
        read(word(2), '(a4)')   ligand_s(i_lig)%atom_name(i_line)
        read(word(3), '(f8.3)') ligand_s(i_lig)%R(1,i_line)
        read(word(4), '(f8.3)') ligand_s(i_lig)%R(2,i_line)
        read(word(5), '(f8.3)') ligand_s(i_lig)%R(3,i_line)
        read(word(6), '(a6)')   ligand_s(i_lig)%mol2_type(i_line)
        read(word(9), '(f8.3)') ligand_s(i_lig)%charge(i_line)
    end if

    ! Read bond information
    if (line(1:13) == '@<TRIPOS>BOND') then
        atm_info = .false.
        mol_info = .false.
        bnd_info = .true.
        i_line = 0
    end if

    if (bnd_info .and. line(1:1) /= '@') then
        call parse_string(line, n_word, word)
        i_line = i_line + 1
        read(word(1), '(i3)') ligand_s(i_lig)%bnd(1,i_line)
        read(word(2), '(i3)') ligand_s(i_lig)%bnd(2,i_line)
        read(word(3), '(i3)') ligand_s(i_lig)%bnd(3,i_line)
        if (n_word == 4) then
            ligand_s(i_lig)%bnd_type(i_line) = word(4)
        end if
    end if
     
    if(line(1:21) == '@<TRIPOS>SUBSTRUCTURE') then
        bnd_info = .false.
    end if
end do

close(i_unit)
  
end subroutine read_multiple_mol2
!------------------------------------------------------------------------------
subroutine check_n_lig(mol2_file, n_lig)
!------------------------------------------------------------------------------
! Count the number of ligand in mol2_file.
!------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: mol2_file
integer, intent(out) :: n_lig
character(len=len_fname) :: line
integer :: i_unit, ioerr

i_unit = 49
open(unit=i_unit,file=trim(mol2_file),action='read',status='old',iostat=ioerr)
if (ioerr/=0) stop

n_lig = 0
do
    read(i_unit,'(a120)',iostat=ioerr) line
    if (ioerr /=0) exit
    if (line(1:17) == '@<TRIPOS>MOLECULE') then
        n_lig = n_lig + 1
    end if
end do

close(i_unit)

end subroutine check_n_lig
!------------------------------------------------------------------------------
subroutine write_mol2(mol2_file, ligand_crd, ref)
!-------------------------------------------------------------------------------
! Write coordinates, atom types, charges, and bond information into mol2 file.
!-------------------------------------------------------------------------------
real(dp), intent(in) :: ligand_crd(:,:)
type(ref_lig_type), intent(in) :: ref
character(len=len_fname), intent(in) :: mol2_file
integer :: i_unit, i_atm, i_bnd

! Open file
call log_p('Writing mol2 file: '//trim(mol2_file), level=40)

i_unit = 32
open(i_unit, file = mol2_file, action = 'write', status = 'replace')

! write molecule info
write(i_unit, '(a)') '@<TRIPOS>MOLECULE'
write(i_unit, '(a)') 'query molecule'
write(i_unit, '(I5,1X,I5,1X,I5,1X,I5,1X,I5)') ref%n_atm, ref%n_bnd, 1,0,0
write(i_unit, '(a)') 'SMALL'
write(i_unit, '(a)') 'NO_CHARGE'
write(i_unit, '(a)') ''

! write atom info
write(i_unit, '(a)') '@<TRIPOS>ATOM'
do i_atm = 1, ref%n_atm
    write(i_unit,'(I7,1X,A4,4X,3F10.4,1X,A6,1X,I4,1X,A3,1X,F10.4)') &
          i_atm, ref%atom_name(i_atm), ligand_crd(:,i_atm),&
          ref%mol2_type(i_atm), 1, 'MOL', ref%charge(i_atm)
end do

! write bond info
write(i_unit, '(a)') '@<TRIPOS>BOND'
do i_bnd = 1, ref%n_bnd
    write(i_unit,'(I6,2(1X,I4),1X,A2)') ref%bnd(:,i_bnd),ref%bnd_type(i_bnd)
end do

!write substructure info
write(i_unit, '(a)') '@<TRIPOS>SUBSTRUCTURE'
write(i_unit, '(a)') '1 MOL         1'

close(i_unit)

end subroutine write_mol2
!------------------------------------------------------------------------------
subroutine write_multiple_mol2(mol2_file, ligand_crd_s, ref, n_lig)
!-------------------------------------------------------------------------------
! Write coordinates, atom types, charges, and bond information into mol2 file.
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: mol2_file
real(dp), intent(in) :: ligand_crd_s(:,:,:)
type(ref_lig_type), intent(in) :: ref
integer, intent(in) :: n_lig
integer :: i_unit, i_atm, i_bnd, i_lig

! Open file
call log_p('Writing mol2 file: '//trim(mol2_file), level=40)

i_unit = 32
open(i_unit, file = mol2_file, action = 'write', status = 'replace')

do i_lig = 1, n_lig
    ! write molecule info
    write(i_unit, '(a)') '@<TRIPOS>MOLECULE'
    write(i_unit, '(a)') 'query molecule'
    write(i_unit, '(I5,1X,I5,1X,I5,1X,I5,1X,I5)') ref%n_atm, ref%n_bnd, 1,0,0
    write(i_unit, '(a)') 'SMALL'
    write(i_unit, '(a)') 'NO_CHARGE'
    write(i_unit, '(a)') ''

    ! write atom info
    write(i_unit, '(a)') '@<TRIPOS>ATOM'
    do i_atm = 1, ref%n_atm
        write(i_unit,'(I7,1X,A4,4X,3F10.4,1X,A6,1X,I4,1X,A3,1X,F10.4)') &
              i_atm, ref%atom_name(i_atm), ligand_crd_s(:,i_atm,i_lig),&
              ref%mol2_type(i_atm), 1, 'MOL', ref%charge(i_atm)
    end do

    ! write bond info
    write(i_unit, '(a)') '@<TRIPOS>BOND'
    do i_bnd = 1, ref%n_bnd
        write(i_unit,'(I6,2(1X,I4),1X,A2)') ref%bnd(:,i_bnd),ref%bnd_type(i_bnd)
    end do

    !write substructure info
    write(i_unit, '(a)') '@<TRIPOS>SUBSTRUCTURE'
    write(i_unit, '(a)') '1 MOL         1'
end do

close(i_unit)

end subroutine write_multiple_mol2
!===============================================================================
! Subroutines related to merge nonpolar hydrogens
!===============================================================================
subroutine mark_nonpolar_hydrogens(ligand, n_nphs, nph)
!-------------------------------------------------------------------------------
! Find nonpolar hydrogens attached to Carbon
! n_nphs: No. of nonpolar hydrogens
! nph: flags to represent polarity of hydrogens.
!      if i_th atom in ligand is nonpolar hydrogen, then nph(i) = .true.
!      Otherwise, nph(i) = .false.
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(inout) :: ligand
integer, intent(out) :: n_nphs
logical, intent(out) :: nph(ligand%n_atm)
integer :: i_atm, i_bnd
integer :: stem_atm, partner_atm
character(len=6) :: mol2_type
logical :: nonpolar

n_nphs = 0
do i_atm = 1, ligand%n_atm
    mol2_type = ligand%mol2_type(i_atm)
    nonpolar = .false.
    
    ! If atom is hydrogen, then check whether connected atom is carbon or not.
    if(trim(mol2_type) == 'H') then
        stem_atm = i_atm
        do i_bnd = 1, ligand%n_bnd
            call find_atom_connected_by_bond(ligand, stem_atm, partner_atm)
            call check_hydrogen_polarity(ligand, partner_atm, nonpolar)
        end do
    end if
       
    if(nonpolar) then
        nph(i_atm) = .true.
        n_nphs = n_nphs + 1
    else if(.not. nonpolar) then
        nph(i_atm) = .false.
    end if
end do

end subroutine mark_nonpolar_hydrogens
!-------------------------------------------------------------------------------
subroutine check_hydrogen_polarity(ligand, partner_atm, nonpolar)
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(in) :: ligand
integer, intent(in) :: partner_atm
logical, intent(inout) :: nonpolar
character(len=6) :: mol2_type

mol2_type = ligand%mol2_type(partner_atm)

if(mol2_type(1:1) == 'C') then
    nonpolar = .true.
else
    nonpolar = .false.
end if

end subroutine check_hydrogen_polarity
!-------------------------------------------------------------------------------
subroutine merge_nonpolar_hydrogens(ligand)
!-------------------------------------------------------------------------------
! Merge nonpolar hydrogens to its connected atom(carbon).
! Charges of nonpolar hydrogens are added to connected atom, and atom/bond
! information of ligand is rearranged.
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(inout) :: ligand
integer :: i_atm, i_bnd, j_atm, j_bnd
integer :: stem_atm, partner_atm
! related to mark nonpolar hydrogens
logical :: nph(ligand%n_atm) 
integer :: n_nphs
! related to merge nonpolar hydrogens
type(ref_lig_type) :: tmp_ligand
integer :: org_atm_idx(ligand%n_atm)

! Find nonpolar hydrogens in ligand
call mark_nonpolar_hydrogens(ligand, n_nphs, nph)

do i_atm = 1, ligand%n_atm
    ! If atom is nonpolar hydrogen, then add its charge to connected atom.
    if(nph(i_atm)) then
        stem_atm = i_atm
        call find_atom_connected_by_bond(ligand, stem_atm, partner_atm)
        ligand%charge(partner_atm) = ligand%charge(partner_atm) + ligand%charge(stem_atm)
    end if
end do

! Merge atom record
tmp_ligand = ligand
j_atm = 0
tmp_ligand%n_atm = ligand%n_atm - n_nphs

do i_atm = 1, ligand%n_atm
    if(.not. nph(i_atm)) then
        j_atm = j_atm + 1
        tmp_ligand%atom_name(j_atm) = ligand%atom_name(i_atm)
        tmp_ligand%R(1:3,j_atm)     = ligand%R(1:3,i_atm)
        tmp_ligand%mol2_type(j_atm) = ligand%mol2_type(i_atm)
        tmp_ligand%charge(j_atm)    = ligand%charge(i_atm)
        org_atm_idx(j_atm)          = i_atm
    end if
end do

! Update bond information
j_bnd = 0
do i_bnd = 1, ligand%n_bnd
    stem_atm = ligand%bnd(2,i_bnd)
    partner_atm = ligand%bnd(3,i_bnd)
    if(nph(stem_atm) .or. nph(partner_atm)) cycle
    j_bnd = j_bnd + 1
    tmp_ligand%bnd(1,j_bnd) = j_bnd
    tmp_ligand%bnd_type(j_bnd) = ligand%bnd_type(i_bnd)
    
    do i_atm = 1, tmp_ligand%n_atm
        if(org_atm_idx(i_atm) == stem_atm) then
            tmp_ligand%bnd(2,j_bnd) = i_atm
            exit
        end if
    end do
   
    do i_atm = 1, tmp_ligand%n_atm
        if(org_atm_idx(i_atm) == partner_atm) then
            tmp_ligand%bnd(3,j_bnd) = i_atm
            exit
        end if
    end do
end do

tmp_ligand%n_bnd = j_bnd
ligand = tmp_ligand

end subroutine merge_nonpolar_hydrogens
!-------------------------------------------------------------------------------
subroutine get_ligand_center(ligand)
!-------------------------------------------------------------------------------
! Find geometrically centered atom of reference ligand
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(inout) :: ligand
real(dp) :: d, min_d, dr(3), center(3)
integer :: i_cnt, ia

! Find geometrical center
center(:) = 0.0d0
do ia = 1, ligand%n_atm
    center(:) = center(:) + ligand%R(:,ia)
end do
center(:) = center(:) / dble(ligand%n_atm)

min_d = 1000.d0
do ia = 1, ligand%n_atm
    dr(:) = ligand%R(:,ia) - center(:)
    d = sqrt(dot_product(dr,dr))
    if (d < min_d) then
        min_d = d
        i_cnt = ia
    end if
end do
ligand%cntr_atm = i_cnt
print*, "ligand center atom: ", ligand%cntr_atm

end subroutine get_ligand_center
!===============================================================================
! Subroutines related to setup geometry of ligand
! e.g. detect ring structure/rotatable bond
!      build tree to reconstruct ligand when torsion angles are perturbed
!      etc.
!===============================================================================
! Subroutines to detect rotatable bond in ligand
!===============================================================================
subroutine detect_rotatable_bond(ligand, rotatable)
!-------------------------------------------------------------------------------
! Detect rotatable bond of ligand from bond information
! If the bond type is '2', '3', 'am', or 'ar' then it is designated to
! non-rotatable
! If the bond type is '1' (ie. single bond), then check atoms connected by the
! bond to determine the bond is rotatable or not.
! Finally, detect ring structures in ligand to set bonds in rings as
! non-rotatable.
!-------------------------------------------------------------------------------
type(ligand_type), intent(inout) :: ligand  
logical, intent(inout) :: rotatable(:)      ! flags for rotatable bond
integer :: ref_lig_no
integer :: i_bnd        ! bond index
integer :: atm1, atm2   ! atoms connected by bond (atm1-atm2)
character(len=2) :: bond_type

rotatable(:) = .true.

ref_lig_no = ligand%lig_type
! First, check bond type & atoms connected by bond.
do i_bnd = 1, ref_lig(ref_lig_no)%n_bnd
    bond_type = ref_lig(ref_lig_no)%bnd_type(i_bnd)
    ! unknown bond type => ERROR
    if (bond_type == 'du' .or. bond_type == 'un' .or. bond_type == 'nc') then
        call terminate_with_error('Wrong bond type. Please check your mol2 file.')
    ! double/triple bond => this bond is not rotatable.
    else if (bond_type == 'am' .or. bond_type == 'ar' .or. &
             bond_type == '2 ' .or. bond_type == '3 ') then
        rotatable(i_bnd) = .false.
    ! single bond => Check atoms connected by the bond.
    else if(bond_type == '1 ') then
        atm1 = ref_lig(ref_lig_no)%bnd(2,i_bnd)
        atm2 = ref_lig(ref_lig_no)%bnd(3,i_bnd)
        call check_connected_atom(ref_lig_no, atm1, atm2, rotatable(i_bnd))
    end if
end do

! Second, detect ring structures in ligand.
! Bonds in ring structure are treated as non-rotatable, because of absence of
! ring-closure algorithm.
call detect_cyclic_bonds(ligand, rotatable)

end subroutine detect_rotatable_bond
!-------------------------------------------------------------------------------
subroutine check_connected_atom(ref_lig_no, atm1, atm2, rotatable_flag)
!-------------------------------------------------------------------------------
! If one of atoms in bond has only one bond, this bond is set to
! non-rotatable because there's no conformation change when rotate torsion
! around this.
!-------------------------------------------------------------------------------
integer, intent(in) :: ref_lig_no, atm1, atm2
logical, intent(inout) :: rotatable_flag
integer :: n_bnd

call count_n_bond_of_atom(ref_lig_no, atm1, n_bnd)
if (n_bnd == 1) rotatable_flag = .false.

call count_n_bond_of_atom(ref_lig_no, atm2, n_bnd)
if (n_bnd == 1) rotatable_flag = .false.

end subroutine check_connected_atom
!-------------------------------------------------------------------------------
subroutine count_n_bond_of_atom(ref_lig_no, atm_no, n_bnd)
!-------------------------------------------------------------------------------
integer, intent(in) :: ref_lig_no, atm_no
integer, intent(out) :: n_bnd
integer :: i_bnd

n_bnd = 0
do i_bnd = 1, ref_lig(ref_lig_no)%n_bnd
    if (ref_lig(ref_lig_no)%bnd(2,i_bnd) == atm_no .or. &
        ref_lig(ref_lig_no)%bnd(3,i_bnd) == atm_no) then
        n_bnd = n_bnd + 1
    end if
end do

end subroutine count_n_bond_of_atom
!-------------------------------------------------------------------------------
subroutine detect_cyclic_bonds(ligand, rotatable)
!-------------------------------------------------------------------------------
! Detect ring sturctures in ligand using graph algorithm and depth first search. 
! If there's ring structure, bonds in ring structure are set to non-rotatable.
! If you want to know how to detect cycles using DFS, please google "Depth first
! search cycle detection".
!-------------------------------------------------------------------------------
type(ligand_type), intent(inout) :: ligand
logical, intent(inout) :: rotatable(:)   ! flag for rotatable bond
integer :: i_atm, i_cn
integer :: atm1, atm2 ! atom index consisting the bond. (atm1-atm2)
integer :: idx1, idx2 ! cn_list index, e.g. (cn_list(idx1,atm2) = atm1)
! Connection list
integer :: cn_list(max_connected_atom, max_lig_atom) ! Connection list (See make_connection_list) 
integer :: n_cn(max_lig_atom)       ! No. of connected atom to each atom
! For Depth First Search
integer :: depth                    ! TODO: Comment by mink
integer :: stack(max_lig_atom)      ! TODO: Comment by mink
integer :: visited(max_connected_atom,max_lig_atom)  ! TODO: Comment by mink
! For ring information
integer :: n_ring, min_ring          ! No. of ring in ligand, minimum ring size
integer :: ring_atm(max_lig_atom)    ! atom index consisting ring structure
logical :: is_ring                   ! flag for ring structure
logical :: is_ring_atm(max_lig_atom) ! flag for atoms consisting ring or not.

! Initialize variables
cn_list(:,:) = 0
n_cn(:) = 0
is_ring_atm(:) = .false.
n_ring = 0

! Make connection list
! .false. means hydrogens are not included in connection list.
call make_connection_list(ligand, cn_list, n_cn, .false.)

do i_atm = 1, ligand%n_atm
    ! initialize cycle detection using Depth First Search
    atm1 = i_atm
    depth = 1
    stack(:) = 0
    visited(:,:) = 0
    stack(depth) = i_atm
    is_ring = .false.
    min_ring = ligand%n_atm 
    ring_atm(:) = 0
    ring_atm(depth) = i_atm

    ! Hydrogen or atoms which connected less than two atoms couldn't make ring
    ! structure.
    if (ref_lig(ligand%lig_type)%mol2_type(i_atm) == 'H ') cycle
    if (n_cn(i_atm) < 2) then
        cycle
    else if (n_cn(i_atm) > max_connected_atom) then
        write(log_msg,'(A,I4,A,I3,A)') "ERROR : Atom, ", i_atm, &
                                ", has more than ", max_connected_atom," bonds. Terminate."
        call terminate_with_error(log_msg)
    end if
  
    ! Depth first search to find ring structure
    do i_cn = 1, n_cn(i_atm)
        atm2 = cn_list(i_cn,i_atm)
        idx2 = i_cn
        !
        visited(idx2,atm1) = visited(idx2,atm1) + 1
        call find_cn_list_idx(atm1, atm2, cn_list, n_cn, idx1)
        visited(idx1,atm2) = visited(idx1,atm2) + 1

        call check_cycle(atm1, atm2, visited, cn_list, n_cn, &
                         depth, stack, min_ring, ring_atm, is_ring)
      
        if (is_ring) then
            call set_cyclic_bonds_rigid(ligand%lig_type, rotatable, &
                                        min_ring, ring_atm, is_ring_atm)
            call update_ring_information(ligand, ring_atm, min_ring, n_ring)
            depth = 1
            stack(:) = 0
            visited(:,:) = 0
            stack(depth) = i_atm
            is_ring = .false.
            min_ring = ligand%n_atm 
            ring_atm(:) = 0
            ring_atm(depth) = i_atm
            cycle
        end if
      
        visited(idx2, atm1) = visited(idx2, atm1) - 1
        call find_cn_list_idx(atm1, atm2, cn_list, n_cn, idx1)
        visited(idx1, atm2) = visited(idx1, atm2) - 1
    end do
end do
 
end subroutine detect_cyclic_bonds
!-------------------------------------------------------------------------------
recursive subroutine check_cycle(origin, v, visited, cn_list, n_cn,&
                                 depth, stack, min_ring, ring_atm, is_ring)
!-------------------------------------------------------------------------------
! Find the smallest ring in ligand using DFS
! Note. "RECURSIVE" subroutine
!-------------------------------------------------------------------------------
integer, intent(inout) :: origin, depth, min_ring
integer, intent(in) :: v
integer, intent(inout) :: visited(:,:), cn_list(:,:)
integer, intent(inout) :: n_cn(:), stack(:), ring_atm(:)
logical, intent(inout) :: is_ring
integer :: i_cn, cn_idx, new_w, new_v, new_depth

new_depth = depth + 1
stack(new_depth) = v

new_w = v
do i_cn = 1, n_cn(v)
    new_v = cn_list(i_cn,v)
    call find_cn_list_idx(new_w, new_v, cn_list, n_cn, cn_idx)
    if (visited(i_cn, new_w) /= 0 .or. &
        visited(cn_idx, new_v) /= 0) cycle
    if (new_v == origin) then
        is_ring = .true.
        if (new_depth <= min_ring) then
            min_ring = new_depth
            ring_atm(:) = 0
            ring_atm(1:min_ring) = stack(1:min_ring)
            return
        end if
    end if

    visited(i_cn, new_w) = visited(i_cn, new_w) + 1
    call find_cn_list_idx(new_w, new_v, cn_list, n_cn, cn_idx)
    visited(cn_idx, new_v) = visited(cn_idx, new_v) + 1
   
    call check_cycle(origin, new_v, visited, cn_list, n_cn, &
                     new_depth, stack, min_ring, ring_atm, is_ring)
   
    visited(i_cn, new_w) = visited(i_cn, new_w) - 1
    call find_cn_list_idx(new_w, new_v, cn_list, n_cn, cn_idx)
    visited(cn_idx, new_v) = visited(cn_idx, new_v) - 1
end do
  
end subroutine check_cycle
!-------------------------------------------------------------------------------
subroutine set_cyclic_bonds_rigid(ref_lig_no, rotatable, &
                                  n_ring_atm, ring_atm_list, is_ring_atm)
!-------------------------------------------------------------------------------
! set bonds connected atoms in ring sturcture as non-rotable.
!-------------------------------------------------------------------------------
integer, intent(in) :: ref_lig_no
logical, intent(inout) :: rotatable(:)
! related to ring structure info.
integer, intent(in) :: n_ring_atm
integer, intent(in) :: ring_atm_list(:)
logical, intent(inout) :: is_ring_atm(:)
integer :: i_atm, j_atm
integer :: atm_1, atm_2, bnd_idx

bnd_idx = 0
do i_atm = 1, n_ring_atm - 1
    atm_1 = ring_atm_list(i_atm)
    is_ring_atm(atm_1) = .true.
    do j_atm = i_atm + 1, n_ring_atm
        atm_2 = ring_atm_list(j_atm)
        call find_lig_bnd_idx(ref_lig_no, atm_1, atm_2, bnd_idx)
        if (bnd_idx == 0) cycle
        rotatable(bnd_idx) = .false.
    end do
end do
is_ring_atm(ring_atm_list(n_ring_atm)) = .true.
atm_1 = ring_atm_list(n_ring_atm)

end subroutine set_cyclic_bonds_rigid
!-------------------------------------------------------------------------------
subroutine update_ring_information(ligand, ring_atoms, n_ring_atm, n_ring)
!-------------------------------------------------------------------------------
! Update ring information(ligand_ring_type in ligand_type).
!-------------------------------------------------------------------------------
type(ligand_type), intent(inout) :: ligand
integer, intent(in) :: ring_atoms(:) ! atom index consisting ring structure
integer, intent(in) :: n_ring_atm    ! number of atoms in ring.
integer, intent(inout) :: n_ring     ! number of rings
character(len=4) :: hybrid
integer :: i_atm, atm_idx
integer :: i_ring, j_atm, n_exist
logical :: aromatic, is_exist

! check the atoms in detected ring is already set to ring member or not.
n_exist = 0
do i_atm = 1, n_ring_atm
    is_exist = .false.
    do i_ring = 1, n_ring
        do j_atm = 1, ligand%rings(i_ring)%n_member
            if(ring_atoms(i_atm) == ligand%rings(i_ring)%member(j_atm)) then
                is_exist = .true.
                exit
            end if
        end do
    end do
    if (is_exist) then
        n_exist = n_exist + 1
    end if
end do

! If all atoms in current ring structure has already set to ring member, just
! return.
if (n_exist == n_ring_atm) return

! New ring structure
n_ring = n_ring + 1
ligand%rings(n_ring)%n_member = n_ring_atm
ligand%rings(n_ring)%member(1:n_ring_atm) = ring_atoms(1:n_ring_atm)

! check aromaticity
aromatic = .true.
do i_atm = 1, n_ring_atm
    atm_idx = ring_atoms(i_atm)
    hybrid = ref_lig(ligand%lig_type)%mol2_type(atm_idx)(3:6)
    if (trim(hybrid) == '2' .or. trim(hybrid) == 'ar' &
        .or. trim(hybrid) == 'am' .or. trim(hybrid) == 'pl3') then
        cycle
    else if(trim(ref_lig(ligand%lig_type)%mol2_type(atm_idx)) == 'O.3') then
        cycle
    else
        aromatic = .false.
    end if
end do

if (aromatic) then
    ligand%rings(n_ring)%aromatic = .true.
else if (.not. aromatic) then
    ligand%rings(n_ring)%aromatic = .false.
end if
ligand%n_ring = n_ring

end subroutine update_ring_information
!===============================================================================
! Subroutines to generate a tree which is used to construct ligand structure 
! when rotate torsions around rotatable bond.
!===============================================================================
subroutine build_tree(ligand)
!-------------------------------------------------------------------------------
! Build tree used to set order for ligand conformation construction based on
! torsion angles.
!-------------------------------------------------------------------------------
type(ligand_type), intent(inout) :: ligand
logical, allocatable :: rotatable(:)  ! flags for rotatable bond
integer :: cn_list(max_connected_atom, max_lig_atom)             ! connection list
integer :: n_cn(max_lig_atom)                   ! number of connected atom
type(lig_frag_type), allocatable :: fragments(:) ! information about rigid fragment
integer, allocatable :: frag_cn_list(:,:), n_frag_cn(:)
integer :: root_idx, i_bnd, br_idx
integer :: i_atm, i_br

allocate(rotatable(ref_lig(ligand%lig_type)%n_bnd))
call detect_rotatable_bond(ligand, rotatable)

! Gen number of fragments(root and branches) in ligand
call get_number_of_branches(ligand, rotatable)

! initialize variables
allocate(fragments(ligand%n_br))
allocate(frag_cn_list(max_arms_in_frag, ligand%n_br))
allocate(n_frag_cn(ligand%n_br))
cn_list(:,:) = 0
n_cn(:) = 0

! Generate connection list include H atoms (what .true. means)
call make_connection_list(ligand, cn_list, n_cn, .true.)
! From rotatable bond information, do fragmentation of ligand
call detect_rigid_fragment(ligand, rotatable, cn_list, n_cn, fragments) 

! Among fragments in ligand, define root which is a starting point for
! constructing ligand conformation from torsion angles
call make_frag_connection_list(ligand, fragments, frag_cn_list, n_frag_cn)
call define_root(fragments, ligand%n_br, frag_cn_list, n_frag_cn, root_idx)

! set root fragments as a 1st branch
br_idx = 1
ligand%n_atm_br(br_idx) = fragments(root_idx)%n_atm
ligand%atm_in_br(1:fragments(root_idx)%n_atm,br_idx) &
     = fragments(root_idx)%atm_list(1:fragments(root_idx)%n_atm)
ligand%bridge(1, br_idx) = 0
ligand%bridge(2, br_idx) = 1

! set branches by walking along rotatable bonds
ligand%n_core_br = fragments(root_idx)%n_arm
do i_bnd = 1, fragments(root_idx)%n_arm
    ligand%core_bridge(i_bnd) = br_idx + 1
    call walk_along_rotatable_bond(ligand, fragments, ligand%n_br, root_idx, br_idx)
end do

call find_idx_for_rotatable_tors(ligand)
call calc_len(ligand)

do i_br = 1, ligand%n_br
    if (ligand%bridge(1, i_br) == 0 .and. ligand%bridge(2, i_br) == 0) then
        call terminate_with_error(&
        "Error: the input ligand MOL2 file has some error; maybe some bond info are missing.")
    end if
end do

ref_lig(ligand%lig_type)%cntr_atm = ligand%atm_in_br(1,1)
print*, "ligand center atom: ", ref_lig(ligand%lig_type)%cntr_atm

deallocate(rotatable)
deallocate(fragments)
deallocate(frag_cn_list)
deallocate(n_frag_cn)

end subroutine build_tree
!-------------------------------------------------------------------------------
subroutine get_number_of_branches(ligand, rotatable)
!-------------------------------------------------------------------------------
! Gen number of branches in ligand.
! No. of branches = No. of rotatable bonds + 1
!-------------------------------------------------------------------------------
type(ligand_type), intent(inout) :: ligand
logical, dimension(:), intent(in) :: rotatable
integer :: i_bnd, i_tor

i_tor = 0
do i_bnd = 1, ref_lig(ligand%lig_type)%n_bnd
    if (.not. rotatable(i_bnd)) cycle
    i_tor = i_tor + 1
end do

ligand%n_br = i_tor + 1 
if (ligand%n_br > max_lig_br) then
    write(log_msg,'(A,I3,A)') "ERROR : the number of branches in ligand exceeds ,",&
                                 max_lig_br, ". Terminate."
    call terminate_with_error(log_msg)
end if
write(log_msg, '(A,I4)') "The number of rotatable torsion angles in target ligand: ", i_tor
call log_p(log_msg, level=30)

end subroutine get_number_of_branches
!-------------------------------------------------------------------------------
subroutine detect_rigid_fragment(ligand, rotatable, cn_list, n_cn, fragments)
!-------------------------------------------------------------------------------
! Detect rigid fragments in ligand by walking along non-rotatable bond
! There are # of rotatable bonds + 1 fragments.
! TODO : Comments by mink
!-------------------------------------------------------------------------------
type(ligand_type), intent(inout) :: ligand
integer, intent(in) :: cn_list(:,:)
integer, intent(in) :: n_cn(:)
logical, intent(in) :: rotatable(:) 
type(lig_frag_type), intent(out) :: fragments(:)
integer :: atm_idx, partner_idx, inv_idx, bnd_idx, i_frag
integer :: start_atm, partner_atm, prev_atm
integer :: n_atm
integer, dimension(max_lig_atom) :: atm_list
integer, dimension(max_connected_atom, max_lig_atom) :: visited
logical :: reverse, new_atm

visited(:,:) = 0

do i_frag = 1, ligand%n_br
    atm_list(:) = 0
    call find_starting_atom(ligand%lig_type, visited, start_atm, &
                            cn_list, n_cn, .true., partner_atm, partner_idx)
   
    if (start_atm == 0) then
        call pick_left_atom(ligand%n_atm, fragments, start_atm, partner_atm, &
                            partner_idx, visited, i_frag-1, n_cn, cn_list)
    end if

    atm_idx = 1
    n_atm = 1
    atm_list(atm_idx) = start_atm
    fragments(i_frag)%atm_list(n_atm) = start_atm
    reverse = .false.
    prev_atm = 0

    ! Walk along rigid bond
    do
        if (atm_idx /= 1 .or. reverse) then
            call find_partner(start_atm, cn_list, n_cn, visited, partner_atm, &
                              partner_idx, prev_atm)
        end if
        call find_lig_bnd_idx(ligand%lig_type, start_atm, partner_atm, bnd_idx)
      
        if (partner_atm == 0) then ! Dead end: move back
            reverse = .true.
            atm_idx = atm_idx - 1
            if (atm_idx == 0) exit
            start_atm = atm_list(atm_idx)
            if (atm_idx == 1) then
                prev_atm = 0
            else
                prev_atm = atm_list(atm_idx - 1)
            end if
        else if(rotatable(bnd_idx) .and. partner_atm /= 0) then ! rotatable bond: move back
            visited(partner_idx, start_atm) = visited(partner_idx, start_atm) + 1
            call find_cn_list_idx(start_atm, partner_atm, cn_list, n_cn, inv_idx)
            visited(inv_idx, partner_atm) = visited(inv_idx, partner_atm) + 1
            reverse = .true.
            start_atm = atm_list(atm_idx)
            if (atm_idx == 1) then
                prev_atm = 0
            else
                prev_atm = atm_list(atm_idx - 1)
            end if
        else ! move forward
            reverse = .false.
            visited(partner_idx, start_atm) = visited(partner_idx, start_atm) + 1
            call find_cn_list_idx(start_atm, partner_atm, cn_list, n_cn, inv_idx)
            visited(inv_idx, partner_atm) = visited(inv_idx, partner_atm) + 1
            atm_idx = atm_idx + 1
            call detect_new_member(fragments(i_frag)%atm_list(1:n_atm),&
                                   n_atm, partner_atm, new_atm)
            if (new_atm) then
                n_atm = n_atm + 1
                fragments(i_frag)%atm_list(n_atm) = partner_atm
            end if
            atm_list(atm_idx) = partner_atm
            start_atm = partner_atm
            prev_atm = atm_list(atm_idx - 1)
        end if
    end do
    fragments(i_frag)%n_atm = n_atm
end do

do i_frag = 1, ligand%n_br
    call set_arms_of_fragment(ligand%lig_type, rotatable, fragments(i_frag))
end do

! assign fragment info in ligand
do i_frag = 1, ligand%n_br
    do atm_idx = 1, fragments(i_frag)%n_atm
        ligand%piece(fragments(i_frag)%atm_list(atm_idx)) = i_frag
    end do
end do

end subroutine detect_rigid_fragment
!-------------------------------------------------------------------------------
subroutine set_arms_of_fragment(ref_lig_no, rotatable, fragment)
!-------------------------------------------------------------------------------
! find arms(=rotatable bonds) of given fragment
!-------------------------------------------------------------------------------
integer, intent(in) :: ref_lig_no
logical, intent(in) :: rotatable(:)
type(lig_frag_type), intent(inout) :: fragment
integer :: i_bnd, i_atm
integer :: atm_idx, n_flex

n_flex = 0
do i_atm = 1, fragment%n_atm
    atm_idx = fragment%atm_list(i_atm)
    do i_bnd = 1, ref_lig(ref_lig_no)%n_bnd
        if (.not. rotatable(i_bnd)) cycle
      
        if (ref_lig(ref_lig_no)%bnd(2,i_bnd) == atm_idx) then
            n_flex = n_flex + 1
            fragment%arm(1,n_flex) = i_bnd
            fragment%arm(2,n_flex) = atm_idx
            fragment%arm(3,n_flex) = ref_lig(ref_lig_no)%bnd(3,i_bnd)
            fragment%arm_used(n_flex) = .false.
        else if(ref_lig(ref_lig_no)%bnd(3,i_bnd) == atm_idx) then
            n_flex = n_flex + 1
            fragment%arm(1,n_flex) = i_bnd
            fragment%arm(2,n_flex) = atm_idx
            fragment%arm(3,n_flex) = ref_lig(ref_lig_no)%bnd(2,i_bnd)
            fragment%arm_used(n_flex) = .false.
        end if
    end do
end do

fragment%n_arm = n_flex

end subroutine set_arms_of_fragment
!-------------------------------------------------------------------------------
subroutine define_root(fragments, n_frag, frag_cn_list, n_frag_cn, root_idx)
!-------------------------------------------------------------------------------
! Define root as a center fragment
!-------------------------------------------------------------------------------
type(lig_frag_type), dimension(:) :: fragments
integer, intent(in) :: n_frag, frag_cn_list(:,:), n_frag_cn(:)
integer, intent(out) :: root_idx
!
integer :: max_dist, dist, stack(n_frag), path(n_frag)
integer :: far_1, far_2, tmp_idx

max_dist = 0
dist = 1
stack(dist) = 1
call find_farthest_frag(fragments, n_frag, 1, frag_cn_list, n_frag_cn, far_1, &
                        stack, dist, path, max_dist)
!
max_dist = 0
dist = 1
stack(dist) = far_1
call find_farthest_frag(fragments, n_frag, far_1, frag_cn_list, n_frag_cn, far_2, &
                        stack, dist, path, max_dist)
!
root_idx = path(int((max_dist+1)/2))
if (mod(max_dist+1,2) == 1) then
    tmp_idx = path(int((max_dist+1)/2)+1)
    if (fragments(root_idx)%n_atm < fragments(tmp_idx)%n_atm) then
        root_idx = tmp_idx
    end if
end if

end subroutine define_root
!-------------------------------------------------------------------------------
recursive subroutine find_farthest_frag(fragments, n_frag, start_frag, &
                                        frag_cn_list, n_frag_cn, farthest, &
                                        stack, dist, path, max_dist)
!-------------------------------------------------------------------------------
type(lig_frag_type), intent(in) :: fragments(:)
integer, intent(in) :: n_frag, start_frag, frag_cn_list(:,:), n_frag_cn(:)
integer, intent(inout) :: farthest, stack(:), dist, path(:), max_dist
!
integer :: i_frag, frag_idx
logical :: visited

do i_frag = 1, n_frag_cn(start_frag)
    frag_idx = frag_cn_list(i_frag, start_frag)
    call check_already_visited(frag_idx, stack, dist, visited)
    if (visited) cycle
    dist = dist + 1
    stack(dist) = frag_idx
    call find_farthest_frag(fragments, n_frag, frag_idx, frag_cn_list, &
                            n_frag_cn, farthest, stack, dist, path, max_dist)
    dist = dist - 1
end do

if (dist > max_dist) then
    max_dist = dist
    path(:) = stack(:)
    farthest = start_frag
end if

end subroutine find_farthest_frag
!-------------------------------------------------------------------------------
subroutine check_already_visited(frag_idx, stack, dist, visited)
!-------------------------------------------------------------------------------
integer, intent(in) :: frag_idx
integer, intent(in) :: stack(:), dist
logical, intent(out) :: visited
integer :: i_frag

visited = .false.
do i_frag = 1, dist
    if (stack(i_frag) == frag_idx) then
        visited = .true.
        exit
    end if
end do

end subroutine check_already_visited
!-------------------------------------------------------------------------------
subroutine walk_along_rotatable_bond(ligand, fragments, n_frag, root_idx, br_idx)
!-------------------------------------------------------------------------------
! walk along rotatable bond and find another branch that connected to root?
!-------------------------------------------------------------------------------
type(ligand_type), intent(inout) :: ligand
type(lig_frag_type), dimension(:), intent(inout) :: fragments
integer, intent(inout) :: br_idx
integer, intent(in) :: root_idx, n_frag
integer, allocatable, dimension(:) :: frag_list, walk_path, prev_frag, history
integer :: frag_idx, m_frag, i_frag, j_frag, k_frag
integer :: start_frag, partner_frag
integer :: stem_1, stem_2, n_atm_br, prev_atm_sum
logical :: reverse, new_frag, connected

allocate(frag_list(n_frag))
allocate(prev_frag(n_frag))
allocate(history(n_frag))
allocate(walk_path(n_frag*2))
reverse = .false.

frag_list(:) = 0
walk_path(:) = 0
frag_idx = 0
m_frag = 0
start_frag = root_idx

! build fragment connection tree starting from root_fragment
! frag_list(m_frag) th fragment is connected to prev_frag(m_frag) th fragment
do
    call identify_next_frag(fragments, start_frag, partner_frag, n_frag)
    if (partner_frag /= 0) then ! connected fragment is found!!
        reverse = .false.
        frag_idx = frag_idx + 1
        walk_path(frag_idx) = partner_frag
     
        ! check partner_frag is new one or not.
        call detect_new_member(frag_list(1:m_frag), m_frag, partner_frag, new_frag)
        if (new_frag) then
            ! if partner_frag is new one, save the fragment index
            m_frag = m_frag + 1
            frag_list(m_frag) = partner_frag 
            prev_frag(m_frag) = start_frag
        end if
        start_frag = partner_frag
    else if(partner_frag == 0) then ! there's no connected fragment. Move back!
        reverse = .true.
        frag_idx = frag_idx - 1
        if (frag_idx == 0) exit ! dead end. 
        start_frag = walk_path(frag_idx) ! return previous fragment index
    end if
end do

do i_frag = 1, m_frag
    call find_stem(fragments(prev_frag(i_frag)), fragments(frag_list(i_frag)), stem_1, stem_2)
    br_idx = br_idx + 1
    ligand%bridge(1,br_idx) = stem_1 ! bridge atom in previous connected fragment  
    ligand%bridge(2,br_idx) = stem_2 ! bridge atom in current fragment

    k_frag = 0
    n_atm_br = 0
    prev_atm_sum = 0
    history(:) = 0

    ! specify all atoms in each branches.
    ! When you rotate torsion angle with index i, 
    ! the atoms in ligand%atm_in_br(:,i+1) are rotated together.
    do j_frag = i_frag, m_frag
        connected = .false.
        if (j_frag /= i_frag) then
            call check_connected(fragments, history, k_frag, frag_list(j_frag), connected)
        else
            connected = .true.
        end if

        if (connected) then
            prev_atm_sum = n_atm_br
            n_atm_br = n_atm_br + fragments(frag_list(j_frag))%n_atm
            ligand%atm_in_br(prev_atm_sum + 1 : n_atm_br, br_idx) &
               = fragments(frag_list(j_frag))%atm_list(1 : fragments(frag_list(j_frag))%n_atm)
            k_frag = k_frag + 1
            history(k_frag) = frag_list(j_frag)
        end if
    end do
    ligand%n_atm_br(br_idx) = n_atm_br
end do
   
deallocate(history)
deallocate(frag_list)
deallocate(prev_frag)
deallocate(walk_path)

end subroutine walk_along_rotatable_bond
!-------------------------------------------------------------------------------
subroutine identify_next_frag(fragments, start_frag, partner_frag, n_frag)
!-------------------------------------------------------------------------------
! Identify undetected partner fragment connected to start fragment.
!-------------------------------------------------------------------------------
type(lig_frag_type), dimension(:), intent(inout) :: fragments
integer, intent(in) :: start_frag, n_frag
integer, intent(out) :: partner_frag
integer :: i_frag, i_bnd, bnd_idx, inv_idx
integer :: stem_1, stem_2
integer :: tmp_stem_1, tmp_stem_2
logical :: partner_found

partner_frag = 0
stem_1 = 0
stem_2 = 0

! find bond index that connect start_frag to unidentified partner_frag
do i_bnd = 1, fragments(start_frag)%n_arm
    ! partner_frag connected to start_frag by i_bnd has been aleady identified.
    if (fragments(start_frag)%arm_used(i_bnd)) cycle
    stem_1 = fragments(start_frag)%arm(2,i_bnd)
    stem_2 = fragments(start_frag)%arm(3,i_bnd)
    bnd_idx = i_bnd
    exit
end do

if (stem_1 == 0 .and. stem_2 == 0) then
    ! There's no unidentified partner_frag. Exit subroutine
    partner_frag = 0
    return
end if

partner_found = .false.
do i_frag = 1, n_frag
    if(i_frag == start_frag) cycle

    do i_bnd = 1, fragments(i_frag)%n_arm
        if(fragments(i_frag)%arm_used(i_bnd)) cycle
      
        tmp_stem_1 = fragments(i_frag)%arm(3,i_bnd)
        tmp_stem_2 = fragments(i_frag)%arm(2,i_bnd)

        ! check i_frag is connected to start_frag using their arm info.
        if(stem_1 == tmp_stem_1 .and. stem_2 == tmp_stem_2) then
            partner_frag = i_frag
            inv_idx = i_bnd
            partner_found = .true.
            exit
        end if
    end do
    if(partner_found) exit
end do
if (.not. partner_found) return

! set arm of fragments as a used one
fragments(start_frag)%arm_used(bnd_idx) = .true.
fragments(partner_frag)%arm_used(inv_idx) = .true.
    
end subroutine identify_next_frag
!-------------------------------------------------------------------------------
subroutine make_connection_list(ligand, cn_list, n_cn, H_include)
!-------------------------------------------------------------------------------
! Make connection list of each ligand atom.
! n_cn(i_atm) = j : the number of connected atom to i_atm is j
! cn_list(i_idx, i_atm) = j_atm : i_idx th connected atom to i_atm is j_atm
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
logical, intent(in) :: H_include
integer, dimension(:,:), intent(inout) :: cn_list
integer, dimension(:), intent(inout) :: n_cn
integer :: ref_lig_no
integer :: i_atm, i_bnd
integer :: i_cn, p_atm_idx

ref_lig_no = ligand%lig_type

do i_atm = 1, ligand%n_atm
    i_cn = 0
    if ((.not. H_include) .and. (ref_lig(ref_lig_no)%mol2_type(i_atm) == 'H ')) cycle
    
    do i_bnd = 1, ref_lig(ref_lig_no)%n_bnd
        if (ref_lig(ref_lig_no)%bnd(2,i_bnd) == i_atm) then
            p_atm_idx = ref_lig(ref_lig_no)%bnd(3,i_bnd)
            if (.not. H_include .and. ref_lig(ref_lig_no)%mol2_type(p_atm_idx) == 'H ') cycle
            i_cn = i_cn + 1
            cn_list(i_cn, i_atm) = p_atm_idx
        else if(ref_lig(ref_lig_no)%bnd(3,i_bnd) == i_atm) then
            p_atm_idx = ref_lig(ref_lig_no)%bnd(2,i_bnd)
            if (.not. H_include .and. ref_lig(ref_lig_no)%mol2_type(p_atm_idx) == 'H ') cycle
            i_cn = i_cn + 1
            cn_list(i_cn, i_atm) = p_atm_idx
        end if
    end do
    n_cn(i_atm) = i_cn
end do

end subroutine make_connection_list
!-------------------------------------------------------------------------------
subroutine make_frag_connection_list(ligand, fragments, frag_cn_list, n_frag_cn)
!-------------------------------------------------------------------------------
! make connection list between fragments
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
type(lig_frag_type), intent(in) :: fragments(:)
integer, intent(inout) :: frag_cn_list(:,:)
integer, intent(inout) :: n_frag_cn(:)
!
integer :: i_frag, j_frag
logical :: connected

n_frag_cn(:) = 0
do i_frag = 1, ligand%n_br-1
    do j_frag = i_frag + 1, ligand%n_br
        call check_connectivity_of_frag(fragments(i_frag), fragments(j_frag), connected)
        if (connected) then
            n_frag_cn(i_frag) = n_frag_cn(i_frag) + 1
            frag_cn_list(n_frag_cn(i_frag),i_frag) = j_frag
            n_frag_cn(j_frag) = n_frag_cn(j_frag) + 1
            frag_cn_list(n_frag_cn(j_frag),j_frag) = i_frag
        end if
    end do
end do

end subroutine make_frag_connection_list
!-------------------------------------------------------------------------------
subroutine check_connectivity_of_frag(frag1, frag2, connected)
!-------------------------------------------------------------------------------
type(lig_frag_type), intent(in) :: frag1, frag2
logical, intent(out) :: connected
!
integer :: stem_1, stem_2

connected = .false.
call find_stem(frag1, frag2, stem_1, stem_2)
if (stem_1 /= 0 .and. stem_2 /= 0) then
    connected = .true.
end if

end subroutine check_connectivity_of_frag
!-------------------------------------------------------------------------------
subroutine calc_len(ligand)
!-------------------------------------------------------------------------------
type(ligand_type), intent(inout) :: ligand
integer :: i_tor, tor_init, tor_end
integer :: i_atm, atm_idx
real(dp) :: dr(3), length

do i_tor = 1, ligand%n_br-1
    tor_init = ligand%bridge(1, i_tor+1)
    tor_end = ligand%bridge(2, i_tor+1)

    do i_atm = 1, ligand%n_atm_br(i_tor+1)
        atm_idx = ligand%atm_in_br(i_atm, i_tor+1)
        dr(1:3) = ref_lig(ligand%lig_type)%R(1:3,atm_idx) &
                - ref_lig(ligand%lig_type)%R(1:3,tor_init)
        length = dot_product(dr,dr)
        ligand%vec_len(i_atm, i_tor+1) = length
    end do
end do

end subroutine calc_len
!===============================================================================
! Related to find index
!===============================================================================
subroutine find_idx_for_rotatable_tors(ligand)
!-------------------------------------------------------------------------------
type(ligand_type), intent(inout) :: ligand
integer :: i_tor, j_tor
integer :: stem_1, stem_2, cand_1, cand_2

do i_tor = 1, ligand%n_br - 1
    stem_1 = ligand%bridge(1, i_tor+1) 
    stem_2 = ligand%bridge(2, i_tor+1)

    do j_tor = 1, ref_lig(ligand%lig_type)%n_dih
        cand_1 = ref_lig(ligand%lig_type)%dih(3,j_tor)
        cand_2 = ref_lig(ligand%lig_type)%dih(4,j_tor)
        if (((stem_1 == cand_1) .and. (stem_2 == cand_2)) &
            .or. ((stem_2 == cand_1) .and. (stem_1 == cand_2))) then
            ligand%i_rot_tor(i_tor) = j_tor
            exit
        end if
    end do
end do

end subroutine find_idx_for_rotatable_tors
!-------------------------------------------------------------------------------
subroutine find_cn_list_idx(atm1, atm2, cn_list, n_cn, cn_idx)
!-------------------------------------------------------------------------------
! Find index of cn_list, cn_idx, which means cn_list(cn_idx, atm2) = atm1
!-------------------------------------------------------------------------------
integer, intent(in) :: atm1, atm2
integer, dimension(:), intent(in) :: n_cn
integer, dimension(:,:), intent(in) :: cn_list
integer, intent(out) :: cn_idx
integer :: i_cn

cn_idx = 0
do i_cn = 1, n_cn(atm2)
    if (cn_list(i_cn, atm2) == atm1) then
        cn_idx = i_cn
        return
    end if
end do

end subroutine find_cn_list_idx
!-------------------------------------------------------------------------------
subroutine find_starting_atom(ref_lig_no, visited, start_atm, cn_list, n_cn, &
                              H_include, partner_atm, partner_idx)
!-------------------------------------------------------------------------------
! TODO : Comments by mink
!-------------------------------------------------------------------------------
integer, intent(in) :: ref_lig_no
integer, dimension(:), intent(in) :: n_cn(:)
integer, intent(in) :: cn_list(:,:)
integer, intent(inout) :: visited(:,:)
logical, intent(in) :: H_include
integer, intent(out) :: start_atm, partner_atm, partner_idx
integer :: i_atm, j_atm, i_idx

start_atm = 0
partner_atm = 0

do i_atm = 1, ref_lig(ref_lig_no)%n_atm
    if(.not. H_include) then
        if(trim(ref_lig(ref_lig_no)%mol2_type(i_atm)) == 'H') then
            cycle
        end if
    end if
   
    do i_idx = 1, n_cn(i_atm)
        j_atm = cn_list(i_idx, i_atm)
        if(visited(i_idx, i_atm) == 0) then
            start_atm = i_atm
            partner_atm = j_atm
            partner_idx = i_idx
            return
        end if
    end do
end do

end subroutine find_starting_atom
!-------------------------------------------------------------------------------
subroutine pick_left_atom(n_lig_atm, fragments, start_atm, partner_atm, partner_idx, &
                          visited, n_frag, n_cn, cn_list)
!-------------------------------------------------------------------------------
! TODO : Comments by mink
!-------------------------------------------------------------------------------
integer, intent(in) :: n_lig_atm
type(lig_frag_type), intent(in) :: fragments(:)
integer, intent(in) :: n_cn(:)
integer, intent(in) :: cn_list(:,:)
integer, intent(in) :: n_frag
integer, intent(in) :: visited(:,:)
integer, intent(out) :: start_atm, partner_atm, partner_idx
logical, allocatable :: in_frag(:)
integer :: i_frag, i_atm, atm_idx

allocate(in_frag(n_lig_atm))

in_frag(:) = .false.
do i_frag = 1, n_frag
    do i_atm = 1, fragments(i_frag)%n_atm
        atm_idx = fragments(i_frag)%atm_list(i_atm)
        in_frag(atm_idx) = .true.
    end do
end do

do i_atm = 1, n_lig_atm
    if (.not. in_frag(i_atm)) then
        start_atm = i_atm
        exit
    end if
end do

deallocate(in_frag)

do i_atm = 1, n_cn(start_atm)
    atm_idx = cn_list(i_atm, start_atm)
    if (visited(i_atm, start_atm) == 1) then
        partner_atm = atm_idx
        partner_idx = i_atm
        return
    end if
end do

end subroutine pick_left_atom
!-------------------------------------------------------------------------------
subroutine find_partner(start_atm, cn_list, n_cn, visited, &
                        partner_atm, partner_idx, prev_atm)
!-------------------------------------------------------------------------------
! TODO : Comments by mink
!-------------------------------------------------------------------------------
integer, intent(in) :: cn_list(:,:)
integer, intent(in) :: n_cn(:)
integer, intent(in) :: start_atm, prev_atm
integer, intent(inout) :: visited(:,:)
integer, intent(inout) :: partner_atm, partner_idx
integer :: i_idx, tmp_atm, min_visited

partner_atm = 0
partner_idx = 0
min_visited = 2
do i_idx = 1, n_cn(start_atm)
    tmp_atm = cn_list(i_idx, start_atm)
    if (tmp_atm == prev_atm) cycle
    if (visited(i_idx, start_atm) < min_visited) then
        min_visited = visited(i_idx, start_atm)
        partner_atm = tmp_atm
        partner_idx = i_idx
    end if
end do

end subroutine find_partner
!-------------------------------------------------------------------------------
subroutine find_lig_bnd_idx(ref_lig_no, start_atm, partner_atm, bnd_idx)
!-------------------------------------------------------------------------------
! TODO : Comments by mink
!-------------------------------------------------------------------------------
integer, intent(in) :: ref_lig_no
integer, intent(in) :: start_atm, partner_atm
integer, intent(out) :: bnd_idx
integer :: i_bnd, atm1, atm2

bnd_idx = 0 ! if bnd_idx = 0, then start_atm and partner_atm don't make a bond
do i_bnd = 1, ref_lig(ref_lig_no)%n_bnd
    atm1 = ref_lig(ref_lig_no)%bnd(2,i_bnd)
    atm2 = ref_lig(ref_lig_no)%bnd(3,i_bnd)

    if ((atm1 == start_atm .and. atm2 == partner_atm) .or. &
        (atm1 == partner_atm .and. atm2 == start_atm)) then
        bnd_idx = i_bnd
        return
    end if
end do

end subroutine find_lig_bnd_idx
!-------------------------------------------------------------------------------
subroutine find_stem(prev_frag, curr_frag, stem_1, stem_2)
!-------------------------------------------------------------------------------
! find atom idx participated in connecting prev_frag to curr_frag
! stem_1 in prev_frag, stem_2 in curr_frag
!-------------------------------------------------------------------------------
type(lig_frag_type), intent(in) :: prev_frag, curr_frag
integer, intent(out) :: stem_1, stem_2
integer :: i_bnd, j_bnd
integer :: prev_tmp_1, prev_tmp_2, curr_tmp_1, curr_tmp_2

stem_1 = 0
stem_2 = 0
do i_bnd = 1, prev_frag%n_arm
   prev_tmp_1 = prev_frag%arm(2,i_bnd)
   prev_tmp_2 = prev_frag%arm(3,i_bnd)
   do j_bnd = 1, curr_frag%n_arm
      curr_tmp_1 = curr_frag%arm(3,j_bnd)
      curr_tmp_2 = curr_frag%arm(2,j_bnd)
      if(prev_tmp_1 == curr_tmp_1 .and. prev_tmp_2 == curr_tmp_2) then
         stem_1 = prev_tmp_1
         stem_2 = prev_tmp_2
         return
      end if
   end do
end do

end subroutine find_stem
!-------------------------------------------------------------------------------
subroutine find_atom_connected_by_bond(ligand, stem_atom, connected_atom)
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(in) :: ligand
integer, intent(in) :: stem_atom
integer, intent(out) :: connected_atom
integer :: i_bnd

do i_bnd = 1, ligand%n_bnd
    if(ligand%bnd(2,i_bnd) == stem_atom) then
        connected_atom = ligand%bnd(3,i_bnd)
        exit
    end if
    if(ligand%bnd(3,i_bnd) == stem_atom) then
        connected_atom = ligand%bnd(2,i_bnd)
        exit
    end if
end do

end subroutine find_atom_connected_by_bond
!-------------------------------------------------------------------------------
subroutine detect_new_member(member_list, n_mem, query_mem, is_new)
!-------------------------------------------------------------------------------
! check query_mem is in member_list or not.
! if query_mem is aleady a member of member_list, is_new = .false.
!-------------------------------------------------------------------------------
integer, dimension(:), intent(in) :: member_list
integer, intent(in) :: n_mem, query_mem
logical, intent(out) :: is_new
integer :: i_mem

is_new = .true.
do i_mem = 1, n_mem
    if(member_list(i_mem) == query_mem) then
        is_new = .false.
    end if
end do

end subroutine detect_new_member
!-------------------------------------------------------------------------------
subroutine check_connected(fragments, history, n_frag, frag_idx, connected)
!-------------------------------------------------------------------------------
! check connectivity between fragments(frag_idx) and fragments(idx in history)
!-------------------------------------------------------------------------------
type(lig_frag_type), dimension(:), intent(in) :: fragments
integer, dimension(:), intent(in) :: history
integer, intent(in) :: n_frag, frag_idx
logical, intent(inout) :: connected
integer :: i_frag, stem_1, stem_2

do i_frag = 1, n_frag
    call find_stem(fragments(history(i_frag)), fragments(frag_idx), stem_1, stem_2)
    if(stem_1 /= 0 .or. stem_2 /= 0) then
        connected = .true.
    end if
end do

end subroutine check_connected
!-------------------------------------------------------------------------------
END MODULE IN_OUT_LIGAND
!-------------------------------------------------------------------------------
