!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/io/in_out_utils.f90
!
! Description:
!  This module contains subroutines can be shared only between in_out modules.
!
!-------------------------------------------------------------------------------
MODULE IN_OUT_UTILS
!-------------------------------------------------------------------------------
use globals
use in_out_vars
use logger, only: log_p, terminate_with_error

implicit none
public

CONTAINS
!-------------------------------------------------------------------------------
subroutine find_res_idx(resname, ref_res_no)
!-------------------------------------------------------------------------------
! find res_idx corresponding to resname
!-------------------------------------------------------------------------------
character(len=4), intent(in) :: resname
integer, intent(out) :: ref_res_no
integer :: i

do i = 1, num_ref_res
    if (resname == ref_res(i)%res_name) then
        ref_res_no = i
        return
    end if
end do

call terminate_with_error('Error: unknown residue name -'//resname//'-')

end subroutine find_res_idx
!-------------------------------------------------------------------------------
subroutine find_het_idx(hetmol)
!-------------------------------------------------------------------------------
! find res_idx corresponding to resname
!-------------------------------------------------------------------------------
type(residue_type), intent(inout) :: hetmol
integer :: i, ref_res_no

ref_res_no = 0
do i = 1, num_ref_res
    if (hetmol%res_name == ref_res(i)%res_name) then
        hetmol%res_type = i
        return
    end if
end do

call terminate_with_error('Error: unknown residue name -'//hetmol%res_name//'-')

end subroutine find_het_idx
!-------------------------------------------------------------------------------
subroutine find_ligand_idx(ligand)
!-------------------------------------------------------------------------------
! find lig_type/res_type corresponding to ligand name
!-------------------------------------------------------------------------------
type(ligand_residue_type), intent(inout) :: ligand
integer :: i, ref_res_no

ref_res_no = 0
do i = 1, n_mol2_top
    if (ligand%res_name == ref_lig(i)%lig_name) then
        ligand%lig_type = i
        ligand%res_type = i+num_ref_res
        return
    end if
end do

! Assign unknown
if (ref_res_no == 0) then
    write(log_msg, "(A,A4,A7,1X,A4)") '- Unknown ligand found', ligand%res_name
    call terminate_with_error(log_msg)
end if

end subroutine find_ligand_idx
!-------------------------------------------------------------------------------
subroutine find_atom_idx(resno, i_ref_res, atom_name, atom_no, error_mode)
!-------------------------------------------------------------------------------
! Find local atomic index (atom_no) corresponding to atom_name from i_ref_res
!-------------------------------------------------------------------------------
integer, intent(in) :: resno, i_ref_res
character(len=6), intent(in) :: error_mode
character(len=4), intent(in) :: atom_name
character(len=4) :: atom, res_name
integer, intent(out) :: atom_no
integer :: i_atm, n_tail, n
character(len=len_fname) :: error_line
logical :: alternative_atom_found

! initialize
atom = atom_name

do
    do i_atm = -3, ref_res(i_ref_res)%n_atm + 1 
        if (atom == ref_res(i_ref_res)%atom_name(i_atm)) then
            atom_no = i_atm
            return
        end if
    end do
    alternative_atom_found = .false.

    n_tail = len(trim(atom))
    res_name = ref_res(i_ref_res)%res_name
    if ((res_name == 'ILE' .or. res_name(2:) == 'ILE') .and. atom == 'CD') then
        atom = 'CD1'
        alternative_atom_found = .true.
    else if ((res_name == 'GLH' .or. res_name(2:) == 'GLH') .and. atom == 'HE1') then
        atom = 'HE2'
        alternative_atom_found = .true.
    else if ((res_name == 'ASH' .or. res_name(2:) == 'ASH') .and. atom == 'HD1') then
        atom = 'HD2'
        alternative_atom_found = .true.
    else if (res_name == 'HOH' .and. atom == 'O') then
        atom = 'OW'
        alternative_atom_found = .true.
    else if (res_name == 'HOH' .and. atom == 'H1') then
        atom = 'HW1'
        alternative_atom_found = .true.
    else if (res_name == 'HOH' .and. atom == 'H2') then
        atom = 'HW2'
        alternative_atom_found = .true.
    else if (atom(1:1) == 'H' .and. atom(n_tail:n_tail) == '1') then
        atom = atom(1:n_tail-1)//'3'
        alternative_atom_found = .true.
    else if (ref_res(i_ref_res)%ter_type == 'N' .and. atom == 'H') then
        atom = 'H1'
        alternative_atom_found = .true.
    end if

    ! In nucleic acids, * for ' may be used for atom names.
    n = len(trim(atom))
    if (atom(n:n) == '*') then
        atom = atom(1:n-1)//"'"
        alternative_atom_found = .true.
    else if (res_name(1:2) == 'DT' .and. atom == 'C5M') then
        atom = 'C7'
        alternative_atom_found = .true.
    end if

    if (.not. alternative_atom_found) exit
end do

if (i_ref_res < 1) then
    write(error_line, "(A,I4)") "WARNING: unknown residue type, ", resno
    call log_p(error_line, me=me, level=0)
end if

if (error_mode /= 'ignore') then
    write(error_line, "(A,A4,A,A4,I4)") "ERROR: unknown atom name ",&
          atom_name, " in ", ref_res(i_ref_res)%res_name, resno
    call log_p(error_line, me=me, level=0)
    call terminate_with_error(error_line)
else
    write(error_line, "(A,A4,A,A4,I4)") "WARNING: unknown atom name ",&
          atom_name, " in ", ref_res(i_ref_res)%res_name, resno
    call log_p(error_line, me=me, level=30)
    atom_no = -100
end if

end subroutine find_atom_idx
!-------------------------------------------------------------------------------
subroutine find_ligand_atom_idx(resNo, i_ref, atom_name, atom_no, error_mode)
!-------------------------------------------------------------------------------
integer, intent(in) :: resNo, i_ref
character(len=6), intent(in) :: error_mode
character(len=4), intent(in) :: atom_name
integer, intent(out) :: atom_no
character(len=len_fname) :: error_line
integer :: i_atm

atom_no = 0
do i_atm = 1, ref_lig(i_ref)%n_atm
    if (atom_name == trim(ref_lig(i_ref)%atom_name(i_atm))) then
        atom_no = i_atm
        return
    end if
end do

if (atom_no == 0) then
    atom_no = -100
    if (error_mode /= 'ignore') then
        write(error_line, "(A,A4,A,A4,I4)") "ERROR: unknown ligand atom name ",&
              atom_name, " in ", ref_lig(i_ref)%lig_name, resNo
        call terminate_with_error(error_line)
    else
        write(error_line, "(A,A4,A,A4,I4)") "WARNING: unknown ligand atom name ",&
              atom_name, " in ", ref_lig(i_ref)%lig_name, resNo
        call log_p(error_line, me=me, level=0)
    end if
end if

end subroutine find_ligand_atom_idx
!-------------------------------------------------------------------------------
subroutine find_atom_cls(eng_para, atom, index)
!-------------------------------------------------------------------------------
! find index for atm als in eng_para
!-------------------------------------------------------------------------------
type(eng_para_type), intent(in) :: eng_para
character(len=6), intent(in) :: atom
integer, intent(out) :: index
integer :: i

do i = 1, eng_para%n_atom_cls
    if (atom == eng_para%atom_cls(i)) then
        index = i
        return
    end if
end do

write(log_msg,"(A,1X,A4)") 'Error in finding atom cls no for ', atom
call terminate_with_error(log_msg)

end subroutine find_atom_cls
!-------------------------------------------------------------------------------
subroutine find_bnd_prm(atm, idx, status)
!-------------------------------------------------------------------------------
! find index for bnd parameter in eng_para structure
!-------------------------------------------------------------------------------
character(len=6), dimension(:), intent(in) :: atm
integer, intent(out) :: idx
logical, intent(out) :: status
character(len=6), dimension(2) :: atom
integer :: i

status = .true.
do i = 1, eng_para%n_bnd
    atom(1:2) = eng_para%atm_in_bnd(1:2,i)
    if ((atm(1) == atom(1) .and. atm(2) == atom(2)) .or. &
        (atm(2) == atom(1) .and. atm(1) == atom(2))) then
        idx = i
        return
    end if
end do

do i = 1, eng_para%n_bnd
    atom(1:2) = eng_para%atm_in_bnd(1:2,i)
    if (atom(2) == 'X' .and. (atom(1) == atm(1) .or. atom(1) == atm(2))) then
        idx = i
        return
    end if
end do

status = .false.
write(log_msg, "(A,2A6)") 'ERROR: failed in getting bnd parm for ', atm(1:2)
call log_p(log_msg, me=me)

end subroutine find_bnd_prm
!-------------------------------------------------------------------------------
subroutine find_ang_prm(atm, ang_type, idx, status)
!-------------------------------------------------------------------------------
! find index for dih ang parameter in eng_para structure
!-------------------------------------------------------------------------------
character(len=6), dimension(:), intent(in) :: atm
integer, intent(in) :: ang_type
integer, intent(out) :: idx
logical, intent(out) :: status
character(len=6), dimension(4) :: atom
integer :: i, imax

status = .true.
if (ang_type == 1) then ! dih ang
    imax = 4 
    idx = 0
    do i = 1, eng_para%n_ang
        if (ang_type == eng_para%ang_type(i)) then
            atom(1:imax) = eng_para%atm_in_ang(1:imax,i)
            if ( ((atom(1) == 'X'    .and. atom(4) == 'X')      .and. &
                 ((atom(2) == atm(2) .and. atom(3) == atm(3))   .or.  &
                  (atom(2) == atm(3) .and. atom(3) == atm(2)))) .or.  &
              
                  (atom(1) == atm(1) .and. atom(2) == atm(2)    .and. &
                   atom(3) == atm(3) .and. atom(4) == atm(4))   .or.  &
                                                                
                  (atom(1) == atm(4) .and. atom(2) == atm(3)    .and. &
                   atom(3) == atm(2) .and. atom(4) == atm(1))   .or.  &
                    
                 ((atom(1) == 'X' .and. atom(3) == 'X' .and. atom(4) == 'X') .and. &
                  (atom(2) == atm(2) .or. atom(2) == atm(3)))    ) then
                idx = i
                if (force_field_type == 'CHARMM') then
                    return ! take the prm that appears first
                end if
            end if
        end if
    end do
    if (idx > 0) return

else if (ang_type == 3) then ! imp ang
    imax = 4
    idx = 0
    do i = 1, eng_para%n_ang
        if (ang_type == eng_para%ang_type(i)) then
            atom(1:imax) = eng_para%atm_in_ang(1:imax,i)
            if (force_field_type == 'CHARMM') then
                if (((atom(2) == 'X'    .and. atom(3) == 'X')      .and. &
                    ((atom(1) == atm(1) .and. atom(4) == atm(4))   .or.  &
                     (atom(1) == atm(4) .and. atom(4) == atm(1)))) .or.  &
                 
                     (atom(1) == atm(1) .and. atom(2) == atm(2)    .and. &
                      atom(3) == atm(3) .and. atom(4) == atm(4))   .or.  &
                    
                     (atom(1) == atm(4) .and. atom(2) == atm(3)    .and. &
                      atom(3) == atm(2) .and. atom(4) == atm(1))   .or.  &
                 
                    ((atom(2) == 'X' .and. atom(3) == 'X' .and. atom(4) == 'X') .and. &
                     (atom(1) == atm(1) .or. atom(1) == atm(4)))    ) then
                    idx = i
                    return
                end if
            else if (force_field_type == 'AMBER') then
                if (((atom(1) == 'X'    .and. atom(2) == 'X')      .and. &
                    ((atom(3) == atm(3) .and. atom(4) == atm(4))   .or.  &
                     (atom(1) == atm(4) .and. atom(2) == atm(3)))) .or.  &
                    
                     (atom(1) == 'X'     .and. &
                    ((atom(2) == atm(2) .and. atom(3) == atm(3) .and. atom(4) == atm(4))   .or. &
                     (atom(1) == atm(4) .and. atom(2) == atm(3) .and. atom(3) == atm(2)))) .or. &
                    
                     (atom(1) == atm(1) .and. atom(2) == atm(2)  .and. &
                      atom(3) == atm(3) .and. atom(4) == atm(4)) .or. &
                    
                     (atom(1) == atm(4) .and. atom(2) == atm(3)  .and. &
                      atom(3) == atm(2) .and. atom(4) == atm(1))  ) then
                    idx = i
                end if
            end if
        end if
    end do
    if (idx > 0) return

else if (ang_type == 2) then
    imax = 3
    do i = 1, eng_para%n_ang ! get bnd angles
        if (ang_type == eng_para%ang_type(i)) then
            atom(1:imax) = eng_para%atm_in_ang(1:imax,i)
            if ( atom(2) == atm(2) .and. &
                ((atom(1) == atm(1) .and. atom(3) == atm(3)) .or. &
                 (atom(1) == atm(3) .and. atom(3) == atm(1)) .or. &
                  atom(1) == 'X'    .and. atom(3) == 'X'   ) ) then
                idx = i
                return
            end if
        end if
    end do
end if

status = .false.
write(log_msg,"(A,I4,4(1x,A6))") 'ERROR: failed in getting ang parm for ', ang_type, atm(1:imax)
call log_p(log_msg, me=me)

end subroutine find_ang_prm
!-------------------------------------------------------------------------------
END MODULE IN_OUT_UTILS
!-------------------------------------------------------------------------------
