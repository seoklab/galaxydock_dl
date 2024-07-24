!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ligdock/physics/autodock4_atom_typer.f90
!
! Description: subroutines are related to autodock4 atom typing.
!-------------------------------------------------------------------------------
MODULE AUTODOCK4_ATOM_TYPER
!-------------------------------------------------------------------------------

use globals
use logger
use convert_res_name, only: convert_to_stdres
!
use energy_vars
use ligdock_E_vars
use ligdock_E_utils, only: find_mol2_res_idx, find_mol2_atm_idx
!
use autodock_vars

implicit none
save
private

public :: assign_atom_types

CONTAINS
!-------------------------------------------------------------------------------
subroutine assign_atom_types(protein, ligand, ad4_param, atdk_types, &
                             is_Hdon, is_Hacc_N, is_Hacc_O, &
                             res_mol2_info, num_mol2_info)
!-------------------------------------------------------------------------------
! Do atom typing of protein/cofactor and target ligand
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
type(ad4_param_type), intent(in) :: ad4_param
integer, intent(out) :: atdk_types(:)
logical, intent(inout) :: is_Hdon(:), is_Hacc_N(:), is_Hacc_O(:)
type(residue_mol2_type), intent(in) :: res_mol2_info(:)
integer, intent(in) :: num_mol2_info

! initialize
atdk_types(:) = -1
is_Hdon(:) = .false.
is_Hacc_N(:) = .false.
is_Hacc_O(:) = .false.

! atom typing for protein/small cofactor
call assign_prot_atom_types(protein, ad4_param, atdk_types, &
                            is_Hdon, is_Hacc_N, is_Hacc_O, &
                            res_mol2_info, num_mol2_info)

! atom typing for protein/large cofactor
call assign_lig_atom_types(protein, ligand, ad4_param, atdk_types, &
                           is_Hdon, is_Hacc_N, is_Hacc_O)

end subroutine assign_atom_types
!-------------------------------------------------------------------------------
subroutine assign_prot_atom_types(protein, ad4_param, atdk_types, &
                                  is_Hdon, is_Hacc_N, is_Hacc_O, &
                                  res_mol2_info, num_mol2_info)
!-------------------------------------------------------------------------------
! Assign AutoDock4 atom type to protein/hetmol
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ad4_param_type), intent(in) :: ad4_param
integer, intent(inout) :: atdk_types(:)
logical, intent(inout) :: is_Hdon(:), is_Hacc_N(:), is_Hacc_O(:)
type(residue_mol2_type), intent(in) :: res_mol2_info(:)
integer, intent(in) :: num_mol2_info
integer :: atm_no, res_no, ref_res_no, i_atm, idx
integer :: mol2_res_idx, mol2_atm_idx
character(len=4) :: res_name, atm_name
character(len=2) :: atm_type
character(len=6) :: mol2_type

do i_atm = 1, tn%nonligatm
    res_no = i_R(1, i_atm)
    atm_no = i_R(2, i_atm)
    
    if (res_no > protein%n_res) then
        ref_res_no = protein%hetmol(res_no - protein%n_res)%res_type
    else
        ref_res_no = protein%residue(res_no)%res_type
    end if

    res_name = ref_res(ref_res_no)%res_name 
    if ((.not. res_name == 'PTR') .and. &
        (.not. res_name == 'SEP') .and. &
        (.not. res_name == 'TPO')) then
        call convert_to_stdres(res_name)
    end if
    atm_name = ref_res(ref_res_no)%atom_name(atm_no)

    ! assign mol2 type of atom
    call find_mol2_res_idx(res_name, mol2_res_idx, res_mol2_info, &
                           num_mol2_info)
    call find_mol2_atm_idx(atm_name, mol2_res_idx, mol2_atm_idx, &
                           res_mol2_info)
    mol2_type = res_mol2_info(mol2_res_idx)%mol2_type(mol2_atm_idx)

    if (mol2_type(1:1) == 'H') then
        ! check hydrogen atom is 'H ', or 'HD'.
        call define_H_type_prot(ref_res_no, i_atm, atm_type)
        if (atm_type == 'HD') then
            is_Hdon(i_atm) = .true.
        end if
    else if (mol2_type(1:2) == 'N.') then
        ! check nitrogen atom is 'N ', or 'NA'.
        if (res_name(1:3) == 'PRO') then
            atm_type = 'N '
        else
            call define_N_type_prot(mol2_type, i_atm, atm_type)
        end if
        if (atm_type == 'NA') then
            is_Hacc_N(i_atm) = .true.
        end if
    else if (mol2_type(1:2) == 'O.') then
        atm_type = 'OA'
        is_Hacc_O(i_atm) = .true.
    else if (mol2_type(1:2) == 'S.') then
        atm_type = 'SA'
        is_Hacc_O(i_atm) = .true.
    else if (mol2_type(1:2) == 'C.') then
        ! check carbon atom is 'C ', or 'A '.
        call define_C_type_prot(ref_res(ref_res_no)%is_sc_atom(atm_no),&
                                res_name, atm_name, mol2_type, atm_type)
    else
        call define_misc_type_prot(atm_name, atm_type)
    end if
    call find_atom_type_idx(ad4_param%atm_type(:), atm_type, idx)
    atdk_types(i_atm) = idx
end do

end subroutine assign_prot_atom_types
!-------------------------------------------------------------------------------
subroutine define_H_type_prot(ref_res_no, atm_idx, atm_type)
!-------------------------------------------------------------------------------
! Assign autodock4 atom type of hydrogen atom in protein/cofactor.
! If hydrogen is connected to carbon, atom type of this atom is 'H '.
! Otherwise, atom type of hydrogein atom is 'HD' (H-bond donor).
!-------------------------------------------------------------------------------
integer, intent(in) :: ref_res_no, atm_idx
character(len=2), intent(out) :: atm_type
!
integer :: partner_idx, partner_no
character(len=4) :: atm_name

partner_idx = i_P(atm_idx)%i_bnd(1,1)
partner_no = i_R(2, partner_idx)

atm_name = ref_res(ref_res_no)%atom_name(partner_no)

if (atm_name(1:1) == 'C') then
    atm_type = 'H '
else
    atm_type = 'HD'  ! H-bond donor
end if


end subroutine define_H_type_prot
!-------------------------------------------------------------------------------
subroutine define_N_type_prot(mol2_type, atm_idx, atm_type)
!-------------------------------------------------------------------------------
! Assign autodock4 atom type of nitrogen atom in protein/cofactor.
! If mol2 type of this nitrogen is N.3, N.2, or N.1, atom type is 'NA'.
! If mol2 type is N.ar, N.am, or N.pl3 and the number of bond around this atom is 2,
! atom type is 'NA'. (otherwise, atom type is 'N ')
!-------------------------------------------------------------------------------
character(len=6), intent(in) :: mol2_type
integer, intent(in) :: atm_idx
character(len=2), intent(out) :: atm_type

if (mol2_type(1:3) == 'N.3' .or. mol2_type(1:3) == 'N.2' &
    .or. mol2_type(1:3) == 'N.1') then
    atm_type = 'NA'
else if (mol2_type(1:4) == 'N.ar' .or. mol2_type(1:4) == 'N.am' &
         .or. mol2_type(1:5) == 'N.pl3') then
    if (i_P(atm_idx)%n_bnd(1) == 2) then
        atm_type = 'NA'     ! H-bond acceptor
    else
        atm_type = 'N '
    end if
else
    atm_type = 'N '
end if

end subroutine define_N_type_prot
!-------------------------------------------------------------------------------
subroutine define_C_type_prot(is_sc_atom, res_name, atm_name, mol2_type, atm_type)
!-------------------------------------------------------------------------------
! Assign autodock4 atom type of carbon atom in protein/cofactor
! If carbon is aromatic carbon, assign 'A '.
! otherwise, atom type is 'C '.
!-------------------------------------------------------------------------------
logical, intent(in) :: is_sc_atom
character(len=4), intent(in) :: res_name, atm_name
character(len=6), intent(in) :: mol2_type
character(len=2), intent(out) :: atm_type

if (res_name(1:3) == 'PHE' .or. res_name(1:3) == 'TYR' .or. &
    res_name(1:3) == 'TRP' .or. res_name(1:3) == 'HIS' .or. &
    res_name(1:3) == 'HID' .or. res_name(1:3) == 'HIE' .or. &
    res_name(1:3) == 'HIP') then
    if (is_sc_atom .and. atm_name(1:2) /= 'CB') then
        atm_type = 'A '
    else
        atm_type = 'C '
    end if
else if (mol2_type(1:4) == 'C.ar') then
    atm_type = 'A '
else
    atm_type = 'C '
end if

end subroutine define_C_type_prot
!-------------------------------------------------------------------------------
subroutine define_misc_type_prot(atm_name, atm_type)
!-------------------------------------------------------------------------------
! If atom is not H, N, O, S, and C, assign atom type based on their name.
!-------------------------------------------------------------------------------
character(len=4), intent(in) :: atm_name
character(len=2), intent(out) :: atm_type

if (atm_name(1:1) == 'P' .or. atm_name(1:1) == 'I' .or. &
    atm_name(1:1) == 'F') then
    atm_type = atm_name(1:1) // ' '
else
    atm_type = atm_name(1:2)
end if

end subroutine define_misc_type_prot
!-------------------------------------------------------------------------------
subroutine assign_lig_atom_types(protein, ligand, ad4_param, atdk_types, &
                                 is_Hdon, is_Hacc_N, is_Hacc_O)
!-------------------------------------------------------------------------------
! Assign AutoDock4 atom type to ligand/ligand_cofactor
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
type(ad4_param_type), intent(in) :: ad4_param
integer, intent(inout) :: atdk_types(:)
logical, intent(inout) :: is_Hdon(:), is_Hacc_N(:), is_Hacc_O(:)
integer :: i_atm, lig_no, atm_no, ref_lig_no, idx
character(len=6) :: mol2_type
character(len=2) :: atm_type
logical :: is_target_lig

do i_atm = tn%nonligatm+1, tn%atom ! for ligand/large cofactor atoms
    is_target_lig = .false.

    lig_no = i_L(1, i_atm - tn%nonligatm)
    atm_no = i_L(2, i_atm - tn%nonligatm)
    ref_lig_no = protein%ligand(lig_no)%lig_type

    if (lig_no == ligand%lig_no) is_target_lig = .true.

    mol2_type = ref_lig(ref_lig_no)%mol2_type(atm_no)
    if(trim(mol2_type) == 'H') then
        call define_H_type(ref_lig(ref_lig_no), atm_no, atm_type)
        if (atm_type == 'HD') then
            is_Hdon(i_atm) = .true.
        end if
    else if(mol2_type(1:2) == 'C.') then
        if (is_target_lig) then
            call define_C_type(mol2_type, atm_type, ligand, atm_no)
        else
            call define_C_type(mol2_type, atm_type)
        end if
    else if(mol2_type(1:2) == 'N.') then
        call define_N_type(ref_lig(ref_lig_no), atm_no, mol2_type, atm_type)
        if (atm_type == 'NA') then
            is_Hacc_N(i_atm) = .true.
        end if
    else if(mol2_type(1:2) == 'O.') then
        call define_O_type(ref_lig(ref_lig_no), atm_no, atm_type)
        if (atm_type == 'OA') then
            is_Hacc_O(i_atm) = .true.
        end if
    else if(mol2_type(1:2) == 'S.') then
        call define_S_type(ref_lig(ref_lig_no), atm_no, atm_type)
        if (atm_type == 'SA') then
            is_Hacc_O(i_atm) = .true.
        end if
    else
        if(mol2_type(2:2) == '.') then
            atm_type = mol2_type(1:1) // ' '
        else
            atm_type = mol2_type(1:2)
        end if
    end if
    call find_atom_type_idx(ad4_param%atm_type(:), atm_type, idx)
    atdk_types(i_atm) = idx
end do

end subroutine assign_lig_atom_types
!-------------------------------------------------------------------------------
subroutine define_H_type(ref, atm_idx, atm_type)
!-------------------------------------------------------------------------------
! Determine atom type of H
! H: Non H-bonding hydorgen
! HD: H-honding hydrogen (as a donor) w/ directionality
! HS: H-bonding hydrogen (as a donor) w/o directionality (not used)
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(in) :: ref
integer, intent(in) :: atm_idx
character(len=2), intent(out) :: atm_type
integer :: i_bnd, partner_idx

do i_bnd = 1, ref%n_bnd
    if (ref%bnd(2,i_bnd) == atm_idx) then
        partner_idx = ref%bnd(3,i_bnd)
        exit
    else if (ref%bnd(3,i_bnd) == atm_idx) then
        partner_idx = ref%bnd(2,i_bnd)
        exit
    end if
end do

if(ref%mol2_type(partner_idx)(1:2) == 'C.') then
    atm_type = 'H '
else
    atm_type = 'HD'
end if

end subroutine define_H_type
!-------------------------------------------------------------------------------
subroutine define_C_type(mol2_type, atm_type, ligand, atm_idx)
!-------------------------------------------------------------------------------
! Determine atom type of C
! A: aromatic carbon
! C: non-aromatic carbon
! If ring information is provided (= present(ligand)), use ring information.
!-------------------------------------------------------------------------------
character(len=6), intent(in) :: mol2_type
character(len=2), intent(out) :: atm_type
type(ligand_type), intent(in), optional :: ligand
integer, intent(in), optional :: atm_idx
integer :: i_ring, i_atm
logical :: aromatic

if (.not. present(ligand)) then
    if (mol2_type == 'C.ar') then
        atm_type = 'A '
    else
        atm_type = 'C '
    end if
else
    aromatic = .false.
    do i_ring = 1, ligand%n_ring
        if (.not. ligand%rings(i_ring)%aromatic) cycle
        do i_atm = 1, ligand%rings(i_ring)%n_member
            if (atm_idx == ligand%rings(i_ring)%member(i_atm)) then
                aromatic = .true.
                exit
            end if
        end do
    end do
    
    if (aromatic) then
        atm_type = 'A '
    else
        atm_type = 'C '
    end if
end if

end subroutine define_C_type
!-------------------------------------------------------------------------------
subroutine define_N_type(ref, atm_idx, mol2_type, atm_type)
!-------------------------------------------------------------------------------
! Determine atom type of N
! N: Non H-bonding nitrogen
! NA: H-bonding nitrogen w/ directionality
! NS: H-bonding nitrogen w/o directionality (not used)
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(in) :: ref
integer, intent(in) :: atm_idx
character(len=6), intent(in) :: mol2_type
character(len=2), intent(out) :: atm_type
integer :: n_bond

if (trim(mol2_type) == 'N.3' .or. trim(mol2_type) == 'N.2' &
    .or. trim(mol2_type) == 'N.1') then
    atm_type = 'NA'
else if(trim(mol2_type) == 'N.ar') then
    call count_bond_number(ref, atm_idx, n_bond)
    if(n_bond == 2) then
        atm_type = 'NA'
    else
        atm_type = 'N '
    end if
else if(trim(mol2_type) == 'N.am' .or. trim(mol2_type) == 'N.pl3') then
    call count_bond_number(ref, atm_idx, n_bond)
    if(n_bond == 2) then
        atm_type = 'NA'
    else
        atm_type = 'N '
    end if
else
    atm_type = 'N '
end if

end subroutine define_N_type
!-------------------------------------------------------------------------------
subroutine define_O_type(ref, atm_idx, atm_type)
!-------------------------------------------------------------------------------
! Determine atom type of O
! OA: H-bonding oxygen w/ directionality
! OS: H-bonding oxygen w/o directionality (not used)
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(in) :: ref
integer, intent(in) :: atm_idx
character(len=2), intent(out) :: atm_type

atm_type = 'OA'

end subroutine define_O_type
!-------------------------------------------------------------------------------
subroutine define_S_type(ref, atm_idx, atm_type)
!-------------------------------------------------------------------------------
! Determine atom type of S
! S: Non-hydrogen bonding sulfur
! SA: H-bonding sulfur
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(in) :: ref
integer, intent(in) :: atm_idx
character(len=2), intent(out) :: atm_type
integer :: n_bond, n_oxygen
logical :: H_accept

H_accept = .true.
call count_bond_number(ref, atm_idx, n_bond)

if (n_bond == 4) then
    call count_free_O(ref, atm_idx, n_oxygen)
    if(n_oxygen > 0) then
        H_accept = .false.
    end if
else if(n_bond == 3) then
    call count_free_O(ref, atm_idx, n_oxygen)
    if(n_oxygen > 0) then
        H_accept = .false.
    end if
end if

if(H_accept) then
    atm_type = 'SA'
else if(.not. H_accept) then
    atm_type = 'S '
end if

end subroutine define_S_type
!-------------------------------------------------------------------------------
subroutine count_free_O(ref, atm_idx, n_oxygen)
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(in) :: ref
integer, intent(in) :: atm_idx
integer, intent(out) :: n_oxygen
integer :: i_bnd, partner_idx, n_heavy

n_oxygen = 0
do i_bnd = 1, ref%n_bnd
    partner_idx = 0
    if (ref%bnd(2,i_bnd) == atm_idx) then
        partner_idx = ref%bnd(3,i_bnd)
    else if (ref%bnd(3,i_bnd) == atm_idx) then
        partner_idx = ref%bnd(2,i_bnd)
    end if

    if(partner_idx /= 0) then
        if(ref%mol2_type(partner_idx)(1:2) == 'O.') then
            call count_heavy_atm(ref, partner_idx, n_heavy)
            if(n_heavy == 1) then
                n_oxygen = n_oxygen + 1
            end if
        end if
    end if
end do

end subroutine count_free_O
!-------------------------------------------------------------------------------
subroutine count_heavy_atm(ref, atm_idx, n_heavy)
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(in) :: ref
integer, intent(in) :: atm_idx
integer, intent(out) :: n_heavy
integer :: n_bond, n_hydrogen, i_bnd, partner_idx

n_bond = 0
n_hydrogen = 0
do i_bnd = 1, ref%n_bnd
    if (ref%bnd(2,i_bnd) == atm_idx) then
        partner_idx = ref%bnd(3, i_bnd)
        n_bond = n_bond + 1
    else if (ref%bnd(3,i_bnd) == atm_idx) then
        partner_idx = ref%bnd(2, i_bnd)
        n_bond = n_bond + 1
    else
        cycle
    end if

    if(trim(ref%mol2_type(partner_idx)) == 'H') then
        n_hydrogen = n_hydrogen + 1
    end if
end do

n_heavy = n_bond - n_hydrogen

end subroutine count_heavy_atm
!-------------------------------------------------------------------------------
subroutine count_bond_number(ref, atm_idx, n_bond)
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(in) :: ref
integer, intent(in) :: atm_idx
integer, intent(out) :: n_bond
integer :: i_bnd

n_bond = 0
do i_bnd = 1, ref%n_bnd
    if (ref%bnd(2,i_bnd) == atm_idx) then
        n_bond = n_bond + 1
    else if (ref%bnd(3,i_bnd) == atm_idx) then
        n_bond = n_bond + 1
    end if
end do

end subroutine count_bond_number
!-------------------------------------------------------------------------------
subroutine find_atom_type_idx(atm_type_list, atm_type, idx)
!-------------------------------------------------------------------------------
character(len=2), intent(in) :: atm_type_list(:), atm_type
integer, intent(out) :: idx
integer :: i

idx = -1
do i = 1, max_atdk_atm_type
    if (atm_type == atm_type_list(i) .or. atm_type_list(i)(1:1) == 'X') then
        idx = i
        return
    end if
end do

write(log_msg,'(A,A)') 'ERROR: There is unknown AutoDock atom type, ', atm_type
call terminate_with_error(log_msg)

end subroutine find_atom_type_idx
!-------------------------------------------------------------------------------
END MODULE AUTODOCK4_ATOM_TYPER
!-------------------------------------------------------------------------------
