!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/energy_setup.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE ENERGY_SETUP
!-------------------------------------------------------------------------------
use globals
use logger
use in_out_utils, only: find_atom_idx
!
use energy_vars

implicit none
private

public :: assign_index_for_R
public :: fill_res_index_array
public :: initialize_atm_parm
public :: setup_dof

CONTAINS
!===============================================================================
subroutine assign_index_for_R(protein)
!-------------------------------------------------------------------------------
! Fill i_R, ii_R, i_L, ii_L.
! 1. i_R: index for R to molecule type
!   R(:,idx) is cartesian coord of i_R(2,idx)th atom in i_R(1,idx)th residue.
!   eg. If i_R(1,idx) <= protein%n_res,
!       R(:,idx) = protein%residue(i_R(1,idx))%R(:,i_R(2,idx))
! 2. ii_R: index for molecule type to R
!   idx th atom in res_no th residue has coord R(:,ii_R(idx,res_no)).
!   eg. If res_no <= protein%n_res,
!       R(:, ii_R(idx,res_no)) = protein%residue(res_no)%R(:,idx)
! 3. i_L: index for R to ligand 
!   If idx > tn%nonligatm,
!       R(:,idx)  = protein%residue(i_L(1,idx))%R(:,i_L(2,idx))
! 4. ii_L : index for ligand to R
!   R(:,ii_L(idx, lig_no)) = protein%ligand(lig_no)%R(:,idx)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer :: atm_no, res_no, het_no, lig_no, ref_res_no, i

atm_no = 0
do res_no = 1, protein%n_res + protein%n_het
    if (res_no <= protein%n_res) then
        ref_res_no = protein%residue(res_no)%res_type
    else
        het_no = res_no - protein%n_res
        ref_res_no = protein%hetmol(het_no)%res_type
    end if

    do i = 1, ref_res(ref_res_no)%n_atm
        atm_no = atm_no + 1
        i_R(1:2, atm_no) = (/ res_no, i /)
        ii_R(i,res_no) = atm_no
    end do
end do

do lig_no = 1, protein%n_lig
    ref_res_no = protein%ligand(lig_no)%lig_type
    do i = 1, ref_lig(ref_res_no)%n_atm
        atm_no = atm_no + 1
        res_no = lig_no + protein%n_res + protein%n_het
        i_R(1:2,atm_no) = (/ res_no, i /)
        i_L(1:2,atm_no - tn%nonligatm) = (/ lig_no, i /)
        ii_L(i,lig_no) = atm_no
    end do
end do

end subroutine assign_index_for_R
!-------------------------------------------------------------------------------
subroutine fill_res_index_array(protein)
!-------------------------------------------------------------------------------
! Fill loc_no information in res_index array
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein

call fill_loc_no_in_res_index(protein)
call fill_res_info_in_res_index(protein) 
!
if (top_type /= 'coarse') then
    call fill_res_index_array_rotamer(protein)
end if

end subroutine fill_res_index_array
!-------------------------------------------------------------------------------
subroutine fill_loc_no_in_res_index(protein)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer :: atm_no, res_no, het_no, ref_res_no, i_res_no, i_ref_res_no
integer :: i

atm_no = 0
do res_no = 1, protein%n_res + protein%n_het
    if (res_no <= protein%n_res) then
        ref_res_no = protein%residue(res_no)%res_type
    else
        het_no = res_no - protein%n_res
        ref_res_no = protein%hetmol(het_no)%res_type
    end if
    do i = 1, ref_res(ref_res_no)%n_atm
        atm_no = atm_no + 1
        res_index(res_no)%loc_no(i) = atm_no ! atom position
    end do
end do

do res_no = 1, protein%n_res
    ref_res_no = protein%residue(res_no)%res_type
    if (res_no > 1 .and. res_no < protein%n_res) then
        i_ref_res_no = protein%residue(res_no-1)%res_type
        call fill_loc_no_w_prev_res(res_no, res_no-1, i_ref_res_no)
        i_ref_res_no = protein%residue(res_no+1)%res_type
        call fill_loc_no_w_next_res(res_no, res_no+1, ref_res_no, i_ref_res_no)

    else if (res_no == 1) then
        i_ref_res_no = protein%residue(res_no+1)%res_type
        call fill_loc_no_w_next_res(res_no, res_no+1, ref_res_no, i_ref_res_no)
        if (protein%cyclic) then
            i_res_no = protein%n_res
            i_ref_res_no = protein%residue(i_res_no)%res_type
            call fill_loc_no_w_prev_res(res_no, i_res_no, i_ref_res_no)
        end if

    else if (res_no == protein%n_res) then
        i_ref_res_no = protein%residue(res_no-1)%res_type
        call fill_loc_no_w_prev_res(res_no, res_no-1, i_ref_res_no)
        if (protein%cyclic) then
            i_ref_res_no = protein%residue(1)%res_type
            call fill_loc_no_w_next_res(res_no, 1, ref_res_no, i_ref_res_no)
        end if
    end if
end do

end subroutine fill_loc_no_in_res_index
!-------------------------------------------------------------------------------
subroutine fill_res_index_array_rotamer(protein)
!-------------------------------------------------------------------------------
! Fill chi_atom indices in res_index
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer :: i_ref, i_rot_res
integer :: i_res, i_ang, i_chi, i_atm, j_atm
character(len=4) :: chi_atoms(4)

do i_res = 1, protein%n_res
    i_ref = protein%residue(i_res)%res_type
    i_rot_res = ref_res(i_ref)%i_rot_res
    !
    res_index(i_res)%n_chi = ref_res(i_ref)%n_chi
    res_index(i_res)%chi_atmid(:,:,:) = 0

    if (res_index(i_res)%n_chi == 0) cycle

    do i_chi = 1, ref_res(i_ref)%n_chi
        i_ang = ref_res(i_ref)%t_ang_for_chi(i_chi)
        chi_atoms = ref_res(i_ref)%atom_name(ref_res(i_ref)%atm_in_t_ang(1:4,i_ang))
        !
        do i_atm = 1, 4
            do j_atm = 1, ref_res(i_ref)%n_atm
                if (trim(ref_res(i_ref)%atom_name(j_atm)) == &
                                        trim(chi_atoms(i_atm))) exit
            end do
            res_index(i_res)%chi_atmid(i_atm,i_chi,1) = ii_R(j_atm,i_res)
            res_index(i_res)%chi_atmid(i_atm,i_chi,2) = j_atm
        end do
    end do
end do

end subroutine fill_res_index_array_rotamer
!-------------------------------------------------------------------------------
subroutine fill_loc_no_w_prev_res(res_no, prev_res_no, prev_ref_res_no)
!-------------------------------------------------------------------------------
! Fill loc_no in res_index with information from previous residue
! loc_no(0) = atom index of 'C' in previous residue
! loc_no(-1) = atom index of 'CA' in previous residue
! loc_no(-3) = atom index of 'O' in previous residue
!-------------------------------------------------------------------------------
integer, intent(in) :: res_no, prev_res_no, prev_ref_res_no

res_index(res_no)%loc_no(0) &
    = res_index(prev_res_no)%loc_no(ref_res(prev_ref_res_no)%i_atm_prev(3)) !-C
res_index(res_no)%loc_no(-1) &
    = res_index(prev_res_no)%loc_no(ref_res(prev_ref_res_no)%i_atm_prev(2)) !-CA
if (top_type /= 'coarse') then
    res_index(res_no)%loc_no(-3) &
        = res_index(prev_res_no)%loc_no(ref_res(prev_ref_res_no)%i_atm_o) !-O
end if

end subroutine fill_loc_no_w_prev_res
!-------------------------------------------------------------------------------
subroutine fill_loc_no_w_next_res(res_no, next_res_no, ref_res_no, next_ref_res_no)
!-------------------------------------------------------------------------------
! Fill loc_no in res_index with information from next residue
! loc_no(n_atm + 1) = atom index of 'N' in next residue
!-------------------------------------------------------------------------------
integer, intent(in) :: res_no, next_res_no, ref_res_no, next_ref_res_no

res_index(res_no)%loc_no(ref_res(ref_res_no)%n_atm+1) &
    = res_index(next_res_no)%loc_no(ref_res(next_ref_res_no)%i_atm_prev(1)) !+N

end subroutine fill_loc_no_w_next_res
!-------------------------------------------------------------------------------
subroutine fill_res_info_in_res_index(protein)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer :: i_res, i_atm, res_no, atm_no, ref_res_no

res_index(:)%n_atm = 0
do i_atm = 1, tn%stdatm
    res_no = i_R(1,i_atm)
    atm_no = i_R(2,i_atm)
    ref_res_no = protein%residue(res_no)%res_type
    ! set basic info.
    res_index(res_no)%ref_res_no = ref_res_no
    res_index(res_no)%n_atm = res_index(res_no)%n_atm + 1
    ! set Ca_id, Cb_id, bb_id
    if (trim(ref_res(ref_res_no)%atom_name(atm_no)) == 'CA') then
        res_index(res_no)%Ca_id(1) = i_atm
        res_index(res_no)%Ca_id(2) = atm_no
        res_index(res_no)%bb_id(2,1) = i_atm
        res_index(res_no)%bb_id(2,2) = atm_no
        if (protein%residue(res_no)%res_name(1:3) == 'GLY' .or. &
            protein%residue(res_no)%res_name(1:4) == 'NGLY' .or. &
            protein%residue(res_no)%res_name(1:4) == 'CGLY') then
            res_index(res_no)%Cb_id(1) = i_atm
            res_index(res_no)%Cb_id(2) = atm_no
        end if
    else if (trim(ref_res(ref_res_no)%atom_name(atm_no)) == 'N') then
        res_index(res_no)%bb_id(1,1) = i_atm
        res_index(res_no)%bb_id(1,2) = atm_no
    else if (trim(ref_res(ref_res_no)%atom_name(atm_no)) == 'C') then
        res_index(res_no)%bb_id(3,1) = i_atm
        res_index(res_no)%bb_id(3,2) = atm_no
    else if (trim(ref_res(ref_res_no)%atom_name(atm_no)) == 'O') then
        res_index(res_no)%bb_id(4,1) = i_atm
        res_index(res_no)%bb_id(4,2) = atm_no
    else if (trim(ref_res(ref_res_no)%atom_name(atm_no)) == 'H') then
        res_index(res_no)%bb_id(5,1) = i_atm
        res_index(res_no)%bb_id(5,2) = atm_no
    else if (trim(ref_res(ref_res_no)%atom_name(atm_no)) == 'CB') then
        res_index(res_no)%Cb_id(1) = i_atm
        res_index(res_no)%Cb_id(2) = atm_no
    end if
end do

if (protein%n_het > 0) then
    do i_atm = tn%stdatm + 1, tn%nonligatm
        res_no = i_R(1,i_atm)
        res_index(res_no)%n_atm = res_index(res_no)%n_atm + 1
    end do
end if

if (protein%n_lig > 0) then
    do i_atm = tn%nonligatm + 1, tn%atom
        res_no = i_R(1,i_atm)
        res_index(res_no)%n_atm = res_index(res_no)%n_atm + 1
    end do
end if

do i_res = 1, tn%stdres
    if (protein%residue(i_res)%ter_type /= 'N') then
        res_index(i_res)%phi_atmid(1,1) = res_index(i_res-1)%bb_id(3,1)
        res_index(i_res)%phi_atmid(2,1) = res_index(i_res)%bb_id(1,1)
        res_index(i_res)%phi_atmid(3,1) = res_index(i_res)%bb_id(2,1)
        res_index(i_res)%phi_atmid(4,1) = res_index(i_res)%bb_id(3,1)
        !
        res_index(i_res)%phi_atmid(1,2) = res_index(i_res-1)%bb_id(3,2)
        res_index(i_res)%phi_atmid(2,2) = res_index(i_res)%bb_id(1,2)
        res_index(i_res)%phi_atmid(3,2) = res_index(i_res)%bb_id(2,2)
        res_index(i_res)%phi_atmid(4,2) = res_index(i_res)%bb_id(3,2)
    else
        res_index(i_res)%phi_atmid(1:4,1) = 0
        res_index(i_res)%phi_atmid(1:4,2) = 0
    end if

    if (protein%residue(i_res)%ter_type /= 'C') then
        res_index(i_res)%psi_atmid(1,1) = res_index(i_res)%bb_id(1,1)
        res_index(i_res)%psi_atmid(2,1) = res_index(i_res)%bb_id(2,1)
        res_index(i_res)%psi_atmid(3,1) = res_index(i_res)%bb_id(3,1)
        res_index(i_res)%psi_atmid(4,1) = res_index(i_res+1)%bb_id(1,1)
        !
        res_index(i_res)%psi_atmid(1,2) = res_index(i_res)%bb_id(1,2)
        res_index(i_res)%psi_atmid(2,2) = res_index(i_res)%bb_id(2,2)
        res_index(i_res)%psi_atmid(3,2) = res_index(i_res)%bb_id(3,2)
        res_index(i_res)%psi_atmid(4,2) = res_index(i_res+1)%bb_id(1,2)
    else
        res_index(i_res)%psi_atmid(1:4,1) = 0
        res_index(i_res)%psi_atmid(1:4,2) = 0
    end if
end do

end subroutine fill_res_info_in_res_index
!-------------------------------------------------------------------------------
subroutine initialize_atm_parm(protein)
!-------------------------------------------------------------------------------
! initialize atm_parm
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer :: i_atm, atm_no, res_no, ref_res_no
character(len=4) :: atm_name

atm_parm(:)%is_H = .false.

do i_atm = 1, tn%stdatm
    res_no = i_R(1,i_atm)
    atm_no = i_R(2,i_atm)
    ref_res_no = protein%residue(res_no)%res_type
    !
    atm_name = trim(ref_res(ref_res_no)%atom_name(atm_no))
    if (atm_name(1:1) == 'H') then
        atm_parm(i_atm)%is_H = .true.
    end if
    !
    atm_parm(i_atm)%fixed = protein%residue(res_no)%atom_fixed(atm_no)
    atm_parm(i_atm)%mass = 0
    atm_parm(i_atm)%excluded = .false.
end do

if (protein%n_het > 0) then
    do i_atm = tn%stdatm + 1, tn%nonligatm
        res_no = i_R(1,i_atm) - protein%n_res
        atm_no = i_R(2,i_atm)
        ref_res_no = protein%hetmol(res_no)%res_type
        !
        !TODO: check if always hydrogen atom startswith 'H'
        atm_name = trim(ref_res(ref_res_no)%atom_name(atm_no))
        if (atm_name(1:1) == 'H') then
            atm_parm(i_atm)%is_H = .true.
        end if
        !
        atm_parm(i_atm)%fixed = protein%hetmol(res_no)%atom_fixed(atm_no)
        atm_parm(i_atm)%mass = 0
    end do
end if
if (protein%n_lig > 0) then
    do i_atm = 1, tn%ligatm
        res_no = i_L(1,i_atm)
        atm_no = i_L(2,i_atm)
        ref_res_no = protein%ligand(res_no)%lig_type
        !
        !TODO: check if always hydrogen atom startswith 'H'
        atm_name = trim(ref_lig(ref_res_no)%atom_name(atm_no))
        if (atm_name(1:1) == 'H') then
            atm_parm(i_atm+tn%nonligatm)%is_H = .true.
        end if
        !
        atm_parm(i_atm+tn%nonligatm)%fixed = protein%ligand(res_no)%atom_fixed(atm_no)
        atm_parm(i_atm+tn%nonligatm)%mass = 0
    end do
end if

end subroutine initialize_atm_parm
!-------------------------------------------------------------------------------
subroutine setup_dof(molecule)
!-------------------------------------------------------------------------------
! Fill dof_id array & count total number of dof(tn%dof). 
! dof_id contains non-fixed atom indices.
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule
integer :: dof_no, atm_no
integer :: i_res, i_atm

allocate(dof_id(tn%atom))
dof_id(:) = 0
dof_no = 0
atm_no = 0

do i_res = 1, molecule%n_res
    do i_atm = 1, ref_res(molecule%residue(i_res)%res_type)%n_atm
        atm_no = atm_no + 1
        if (.not. molecule%residue(i_res)%atom_fixed(i_atm)) then
            dof_no = dof_no + 1
            dof_id(dof_no) = atm_no
        end if
    end do
end do

do i_res = 1, molecule%n_het
    do i_atm = 1, ref_res(molecule%hetmol(i_res)%res_type)%n_atm
        atm_no = atm_no + 1
        if (.not. molecule%hetmol(i_res)%atom_fixed(i_atm)) then
            dof_no = dof_no + 1
            dof_id(dof_no) = atm_no
        end if
    end do
end do

do i_res = 1, molecule%n_lig
    do i_atm = 1, ref_lig(molecule%ligand(i_res)%lig_type)%n_atm
        atm_no = atm_no + 1
        if (.not. molecule%ligand(i_res)%atom_fixed(i_atm)) then
            dof_no = dof_no + 1
            dof_id(dof_no) = atm_no
        end if
    end do
end do

tn%dof = dof_no*3

end subroutine setup_dof
!-------------------------------------------------------------------------------
END MODULE ENERGY_SETUP
!-------------------------------------------------------------------------------
