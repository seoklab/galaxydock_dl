!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/utils/convert_res_name.f90
!
! Description:
!  This module contains residue name related subroutines
!   - Conversions between 3-letter code and 1-letter code
!   - Conversion from non-standard residue name to standard residue name
!
!-------------------------------------------------------------------------------
MODULE CONVERT_RES_NAME
!-------------------------------------------------------------------------------
use logger, only: terminate_with_error

implicit none
public

CONTAINS
!-------------------------------------------------------------------------------
subroutine find_long_res_name(short_name, long_name)
!-------------------------------------------------------------------------------
! Convert one-letter amino acid code into three-letter code
!-------------------------------------------------------------------------------
character(len=1), intent(in) :: short_name
character(len=4), intent(out) :: long_name
integer, parameter :: n_res = 20
character(len=1), parameter :: res_name(n_res) = &
                        (/ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', &
                           'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' /)
character(len=4), parameter :: residue_name(n_res) = &
    (/ 'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', &
       'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR' /)
integer :: i_res

do i_res = 1, n_res
    if (short_name == res_name(i_res)) then
        long_name = residue_name(i_res)
        return
    end if
end do

call terminate_with_error('Error: unknown short residue name '//short_name)

end subroutine find_long_res_name
!-------------------------------------------------------------------------------
subroutine find_short_res_name(short_name, long_name)
!-------------------------------------------------------------------------------
! Convert three-letter amino acid code into one-letter code
!-------------------------------------------------------------------------------
character(len=1), intent(out) :: short_name
character(len=4), intent(in) :: long_name
integer, parameter :: n_res = 22
character(len=4), parameter :: residue_name(n_res) = &
    (/ 'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'HIP', 'HID', 'ILE', 'LYS', &
       'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR' /)
character(len=1), parameter :: res_name(n_res) = &
    (/ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'H', 'H', 'I', 'K', &
       'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' /)
integer :: i_res

do i_res = 1, n_res
    if (long_name == residue_name(i_res)) then
        short_name = res_name(i_res)
        return
    end if
end do

call terminate_with_error('Error: unknown short name of residue : '//long_name)

end subroutine find_short_res_name
!-------------------------------------------------------------------------------
subroutine convert_to_stdres(resname)
!-------------------------------------------------------------------------------
! Return residue name as one of 20 standard amino acids
!-------------------------------------------------------------------------------
character(len=4), intent(inout) :: resname

if (resname(4:4) /= ' ' .and. (resname(1:1) == 'N' .or. (resname(1:1) == 'C'))) then
    resname = resname(2:4)
end if

if (resname == 'HID' .or. resname == 'HIE' .or. resname == 'HIP') then
    resname = 'HIS'
else if (resname == 'CYM' .or. resname == 'CYX') then
    resname = 'CYS'
else if (resname == 'LYN') then
    resname = 'LYS'
else if (resname == 'ASH') then
    resname = 'ASP'
else if (resname == 'GLH') then
    resname = 'GLU'
else if (resname == 'CPR') then
    resname = 'PRO'
else if (resname == 'PTR') then
    resname = 'TYR'
else if (resname == 'SEP') then
    resname = 'SER'
else if (resname == 'TPO') then
    resname = 'THR'
end if
  
end subroutine convert_to_stdres
!-------------------------------------------------------------------------------
END MODULE CONVERT_RES_NAME
!-------------------------------------------------------------------------------
