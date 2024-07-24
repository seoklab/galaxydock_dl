!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ligdock/physics/AutoDock3_setup.f90
!
! Description: Subroutines related setup autodock3 parameters
!
!-------------------------------------------------------------------------------
MODULE AUTODOCK3_SETUP
!-------------------------------------------------------------------------------

use globals
use logger
use string, only: parse_string 
use convert_res_name, only: convert_to_stdres
!
use energy_vars
use ligdock_E_utils, only: check_is_new_type_string, check_is_new_type
! 
use autodock_vars
use autodock_utils, only: initialize_atdk_para, set_nb_matrix, set_energy_table,&
                          HD_hbond_parameter, OA_SA_hbond_parameter

implicit none
save
private

public :: construct_atdk3_para
public :: read_assign_atdk3_para
public :: doping_prot_sol_par
public :: setup_atdk3_E

CONTAINS
!-------------------------------------------------------------------------------
subroutine construct_atdk3_para(protein, ligand, is_Hdon, is_Hacc_O)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
logical, intent(out) :: is_Hdon(:), is_Hacc_O(:)
integer :: n_tot_type, n_lig_type
integer :: i_atm, res_no, atm_no, ref_res_no
integer :: type_idx
character(len=2) :: atom
logical :: new_type

call initialize_atdk_para(atdk_para)

n_tot_type = 0
n_lig_type = 0
is_ligand(:) = .false.

do i_atm = 1, tn%nonligatm ! for protein/small cofactor
    res_no = i_R(1, i_atm)
    atm_no = i_R(2, i_atm)
    if (res_no > tn%stdres) then
        ref_res_no = protein%hetmol(res_no-tn%stdres)%res_type
        call check_2letter_name(ref_res(ref_res_no)%atom_name(atm_no),&
                                protein%hetmol(res_no-tn%stdres)%pdb_res_name(1:2), atom)
    else
        ref_res_no = protein%residue(res_no)%res_type
        call check_2letter_name(ref_res(ref_res_no)%atom_name(atm_no),&
                                protein%residue(res_no)%pdb_res_name(1:2), atom)
    end if

    call check_is_new_type_string(atdk_para%type_name(:), n_tot_type, atom, &
                                  new_type, type_idx)

    if (new_type) then
        n_tot_type = n_tot_type + 1
        atdk_para%type_idx(n_tot_type) = n_tot_type
        atdk_para%type_name(n_tot_type) = atom
        atdk_para_idx(i_atm) = n_tot_type
    else
        atdk_para_idx(i_atm) = type_idx
    end if
end do

do i_atm = tn%nonligatm+1, tn%atom ! for ligand/large cofactor
    res_no = i_L(1, i_atm-tn%nonligatm)
    atm_no = i_L(2, i_atm-tn%nonligatm)
    ref_res_no = protein%ligand(res_no)%lig_type
    if (res_no == ligand%lig_no) then ! for target ligand
        is_ligand(i_atm) = .true.
        call check_2letter_name(ref_lig(ref_res_no)%mol2_type(atm_no)(1:2), &
                                ref_lig(ref_res_no)%mol2_type(atm_no)(1:2), atom)
        if (atom == 'C ') then
            call check_aromatic(ligand, atm_no, atom)
        end if
        call check_is_new_type_string(atdk_para%type_name(:), n_tot_type, &
                                      atom, new_type, type_idx)
        if (new_type) then
            n_tot_type = n_tot_type + 1
            atdk_para%type_idx(n_tot_type) = n_tot_type
            atdk_para%type_name(n_tot_type) = atom
            atdk_para_idx(i_atm) = n_tot_type
        else
            atdk_para_idx(i_atm) = type_idx
        end if

        call check_is_new_type(atdk_para%lig_atom_types(:), n_lig_type, &
                               atdk_para_idx(i_atm), new_type, type_idx)
        if (new_type) then
            n_lig_type = n_lig_type + 1
            atdk_para%lig_atom_types(n_lig_type) = atdk_para_idx(i_atm)
        end if
                
    else  ! for large cofactor
        call check_2letter_name(ref_lig(ref_res_no)%mol2_type(atm_no)(1:2), &
                                ref_lig(ref_res_no)%mol2_type(atm_no)(1:2), atom)
        if (ref_lig(ref_res_no)%mol2_type(atm_no)(1:4) == 'C.ar') then
            atom = 'A '
        end if
        call check_is_new_type_string(atdk_para%type_name(:), n_tot_type, &
                                      atom, new_type, type_idx)
        if (new_type) then
            n_tot_type = n_tot_type + 1
            atdk_para%type_idx(n_tot_type) = n_tot_type
            atdk_para%type_name(n_tot_type) = atom
            atdk_para_idx(i_atm) = n_tot_type
        else
            atdk_para_idx(i_atm) = type_idx
        end if
    end if
end do

atdk_para%n_atm_types(1) = n_tot_type
atdk_para%n_atm_types(2) = n_lig_type

call assign_Hdon_Hacc(atdk_para_idx, atdk_para, is_Hdon, is_Hacc_O)

end subroutine construct_atdk3_para
!-------------------------------------------------------------------------------
subroutine read_assign_atdk3_para(file_name, para)
!-------------------------------------------------------------------------------
! Read autodock3 energy parameter file, and store parameters in atdk_para.
! For parameters used in solvation, 
!   - solvation parameters of ligand : store them to atdk_para
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: file_name
type(atdk_eng_para_type), intent(inout) :: para
integer :: f_unit, ioerr
character(len=3) :: comp1
character(len=2) :: comp2
character(len=2) :: atm1, atm2
real(dp) :: value
logical :: yn1, yn2
integer :: i_prot, i_lig

f_unit = 111

open(unit=f_unit, file=trim(file_name), position='rewind', status='old',&
     action='read', iostat=ioerr)
if (ioerr /=0 ) then
    call terminate_with_error("ERROR: Can not read autodock3 energy parameter file.")
end if

102 format(A3,T6,A2,T10,A2,T13,A2,T17,F11.8)

do
    read(f_unit, 102, iostat=ioerr) comp1, comp2, atm1, atm2, value
    if (ioerr/=0) exit
    if (comp1 == "EOF") exit

    call check_is_new_type_string(para%type_name(:), para%n_atm_types(1),&
                                  atm1, yn1, i_prot)
    call check_is_new_type_string(para%type_name(:), para%n_atm_types(1),&
                                  atm2, yn2, i_lig)

    if ((.not. yn2) .and. comp1 == 'Sol') then
        if (comp2 == 'pa') then
            para%solpar(i_lig) = value/1000.0d0
        else if (comp2 == 'cn') then
            para%solcn(i_lig) = value
        end if
    else if ((.not. yn1) .and. (.not. yn2)) then
        if (comp1 == 'Rij') then
            if (comp2 == 'lj') then
                para%ljrij(i_lig, i_prot) = value
            else if (comp2 == 'hb') then
                para%hbrij(i_lig, i_prot) = value
                para%hbond_read(i_lig, i_prot) = .true.
            end if
        else if (comp1 == 'Eij') then
            if (comp2 == 'lj') then
                para%ljeij(i_lig, i_prot) = value
            else if (comp2 == 'hb') then
                para%hbeij(i_lig, i_prot) = value
                para%hbond_read(i_lig, i_prot) = .true.
            end if
        end if
    end if
end do

close(f_unit)

end subroutine read_assign_atdk3_para
!-------------------------------------------------------------------------------
subroutine doping_prot_sol_par(file_name, protein, sol_V)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: file_name
type(molecule_type), intent(in) :: protein
real(dp), intent(out) :: sol_V(:)
!
integer :: f_unit, ioerr
integer :: num_word
character(len=len_fname) :: line, word(20)
!
integer :: i_line, n_line
real(dp) :: vol(200)
character(len=3) :: p_res_name(200), p_atm_name(200)
integer :: i_atm, res_no, atm_no, ref_res_no
character(len=4) :: res_name, atom_name
logical :: is_assigned

sol_V(:) = 0.0d0

f_unit = 35
open(unit=f_unit, file=trim(file_name), status='old', action='read', &
     iostat = ioerr)
if (ioerr /=0) then
    call terminate_with_error('ERROR: Can not find solvation parameter file.')
end if

i_line = 0
do
    read(f_unit, '(a25)', iostat=ioerr)line
    if (ioerr /= 0) exit
    call parse_string(line, num_word, word)
    i_line = i_line + 1
    p_atm_name(i_line) = word(1)
    p_res_name(i_line) = word(2)
    read(word(3), '(f8.3)') vol(i_line)
end do

close(f_unit)

n_line = i_line

do i_atm = 1, tn%nonligatm
    res_no = i_R(1, i_atm)
    atm_no = i_R(2, i_atm)
    if (res_no > protein%n_res) then
        !debug
        if (.not. protein%hetmol(res_no-protein%n_res)%atm_read(atm_no)) cycle
        res_name = protein%hetmol(res_no-protein%n_res)%res_name
        ref_res_no = protein%hetmol(res_no-protein%n_res)%res_type
    else
        !debug
        if (.not. protein%residue(res_no)%atm_read(atm_no)) cycle
        res_name = protein%residue(res_no)%res_name
        ref_res_no = protein%residue(res_no)%res_type
    end if
    call convert_to_stdres(res_name)

    atom_name = ref_res(ref_res_no)%atom_name(atm_no)
    if (trim(atom_name) == 'C' .and. trim(res_name) /= 'CXL') then
        sol_V(i_atm) = 9.82d0
    else if (trim(atom_name) == 'O') then
        sol_V(i_atm) = 8.17d0
    else if (trim(atom_name) == 'N' .and. trim(res_name) /= 'AMN') then
        sol_V(i_atm) = 9.00d0
    else if (trim(atom_name) == 'CA' .and. trim(res_name) /= 'CA') then
        sol_V(i_atm) = 9.40d0
    else
        is_assigned = .false.
        do i_line = 1, n_line
            if (trim(atom_name) == trim(p_atm_name(i_line)) .and. &
                trim(res_name) == trim(p_res_name(i_line))) then
                sol_V(i_atm) = vol(i_line)
                is_assigned = .true.
                exit
            end if
        end do
        if (.not. is_assigned) then
            if (atom_name(1:1) == 'C') then
                sol_V(i_atm) = 12.77d0
            end if
        end if
    end if
end do 
    
do i_atm = tn%nonligatm+1, tn%atom
    if (is_ligand(i_atm)) cycle
    res_no = i_L(1, i_atm - tn%nonligatm)
    atm_no = i_L(2, i_atm - tn%nonligatm)
    ref_res_no = protein%ligand(res_no)%lig_type
    res_name = ref_lig(ref_res_no)%lig_name
    atom_name = ref_lig(ref_res_no)%atom_name(atm_no)
    
    if (trim(atom_name) == 'C' .and. trim(res_name) /= 'CXL') then
        sol_V(i_atm) = 9.82d0
    else if (trim(atom_name) == 'O') then
        sol_V(i_atm) = 8.17d0
    else if (trim(atom_name) == 'N' .and. trim(res_name) /= 'AMN') then
        sol_V(i_atm) = 9.00d0
    else if (trim(atom_name) == 'CA' .and. trim(res_name) /= 'CA') then
        sol_V(i_atm) = 9.40d0
    else
        is_assigned = .false.
        do i_line = 1, n_line
            if (trim(atom_name) == trim(p_atm_name(i_line)) .and. &
                trim(res_name) == trim(p_res_name(i_line))) then
                sol_V(i_atm) = vol(i_line)
                is_assigned = .true.
                exit
            end if
        end do
        if (.not. is_assigned) then
            if (atom_name(1:1) == 'C') then
                sol_V(i_atm) = 12.77d0
            end if
        end if
    end if
end do 

end subroutine doping_prot_sol_par
!-------------------------------------------------------------------------------
subroutine assign_Hdon_Hacc(atdk_para_idx, atdk_para, is_Hdon, is_Hacc_O)
!-------------------------------------------------------------------------------
integer, intent(in) :: atdk_para_idx(:)
type(atdk_eng_para_type), intent(in) :: atdk_para
logical, intent(out) :: is_Hdon(:), is_Hacc_O(:)
integer :: i_atm, stem_atm
character(len=2) :: atm_type, stem_type

is_Hdon(:) = .false.
is_Hacc_O(:) = .false.

do i_atm = 1, tn%atom
    atm_type = atdk_para%type_name(atdk_para_idx(i_atm))

    if (trim(atm_type) == 'O') then
        is_Hacc_O(i_atm) = .true.
    else if (trim(atm_type) == 'H') then
        stem_atm = i_P(i_atm)%i_bnd(1,1)
        stem_type = atdk_para%type_name(atdk_para_idx(stem_atm))
        if (trim(stem_type) == 'N' .or. trim(stem_type) == 'O') then
            is_Hdon(i_atm) = .true.
        end if
    end if
end do

end subroutine assign_Hdon_Hacc
!-------------------------------------------------------------------------------
subroutine check_2letter_name(org_atom, org_res, mod_atom)
!-------------------------------------------------------------------------------
character(len=2), intent(in) :: org_atom
character(len=2), intent(in) :: org_res
character(len=2), intent(out) :: mod_atom

if (org_atom == 'CL' .or. org_atom == 'Cl') then
    mod_atom = 'c '
else if (org_atom == 'BR' .or. org_atom == 'Br') then
    mod_atom = 'b '
else if (org_res == 'FE' .or. org_res == 'Fe') then
    mod_atom = 'f '
else if (org_res == 'MG' .or. org_res == 'Mg') then
    mod_atom = 'M '
else if (org_res == 'CA' .or. org_res == 'Ca') then
    mod_atom = 'L '
else
    mod_atom = org_atom(1:1)
end if
end subroutine check_2letter_name
!-------------------------------------------------------------------------------
subroutine check_aromatic(ligand, atm_no, atom)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: atm_no
character(len=2), intent(inout) :: atom
integer :: i_ring, i_atm

do i_ring = 1, ligand%n_ring
    if (.not. ligand%rings(i_ring)%aromatic) cycle
    do i_atm = 1, ligand%rings(i_ring)%n_member
        if (atm_no == ligand%rings(i_ring)%member(i_atm)) then
            atom = 'A '
            return
        end if
    end do
end do

end subroutine check_aromatic
!-------------------------------------------------------------------------------
subroutine setup_atdk3_E(ligand, is_Hdon, is_Hacc_O)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
logical, intent(in) :: is_Hdon(:), is_Hacc_O(:)

call set_nb_matrix(ligand)
call set_energy_table(atdk_para%n_atm_types(1), atdk_para%n_atm_types(2))

allocate(atdk_hbond(tn%atom))
call initialize_atdk3_hbond(atdk_hbond, is_Hdon, is_Hacc_O)

end subroutine setup_atdk3_E
!-------------------------------------------------------------------------------
subroutine initialize_atdk3_hbond(atdk_hbond, is_Hdon, is_Hacc_O)
!-------------------------------------------------------------------------------
type(atdk_hbond_type), intent(inout) :: atdk_hbond(:)
logical, intent(in) :: is_Hdon(:), is_Hacc_O(:)
integer :: i_atm

! initialize
do i_atm = 1, tn%atom
    atdk_hbond(i_atm)%rexp = 1
    atdk_hbond(i_atm)%rvec(:,:) = 0.0d0
    atdk_hbond(i_atm)%racc = 1.0d0
    atdk_hbond(i_atm)%rdon = 1.0d0
end do

do i_atm = 1, tn%atom
    if (is_Hdon(i_atm)) then
        call HD_hbond_parameter(atdk_hbond(i_atm), i_atm)
    else if (is_Hacc_O(i_atm)) then
        call OA_SA_hbond_parameter(atdk_hbond(i_atm), i_atm)
    end if
end do

end subroutine initialize_atdk3_hbond
!-------------------------------------------------------------------------------
END MODULE AUTODOCK3_SETUP
!-------------------------------------------------------------------------------
