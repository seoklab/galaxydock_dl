!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ligdock/ligdock_E_utils.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE LIGDOCK_E_UTILS
!-------------------------------------------------------------------------------

use globals
use logger
use string
use mathfunctions
use energy_vars
use ligdock_E_vars

implicit none
save
public

CONTAINS
!-------------------------------------------------------------------------------
subroutine read_aa_to_mol2(file_name, residue_mol2_info, num_mol2_info)
!-------------------------------------------------------------------------------
! Read data file to map mol2 type of atoms in standard amino acids.
! data_file: infile_aa_to_mol2 (Default: $(data_dir)/X_score_res.prm)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: file_name
type(residue_mol2_type), intent(out) :: residue_mol2_info(:)
integer, intent(out) :: num_mol2_info
!
integer :: f_unit, ioerror
integer :: num_word
character(len=len_fname) :: line, word(25)
! 
integer :: res_idx, atm_idx
character(len=4) :: res_name

f_unit = 28
open(unit=f_unit, file=trim(file_name), action='read', iostat=ioerror)
if (ioerror /= 0) then
    call terminate_with_error("Cannot find data file to assign mol2 type of protein")
end if
  
res_idx = 0
atm_idx = 0

do
    read(f_unit, '(a80)', iostat=ioerror) line
    if(ioerror /= 0) exit
   
    if (line(1:1) == "!" .or. line(1:1) == "#") cycle

    call parse_string(line, num_word, word)
    if(word(1) == 'RESI') then 
        res_name = word(2)
        res_idx = res_idx + 1
        atm_idx = 0
        residue_mol2_info(res_idx)%res_name = res_name
    end if

    if(word(1) == 'ATOM') then
        atm_idx = atm_idx + 1
        residue_mol2_info(res_idx)%n_atm = atm_idx
        residue_mol2_info(res_idx)%atom_name(atm_idx) = word(2)
        residue_mol2_info(res_idx)%mol2_type(atm_idx) = word(3)
    end if
end do
num_mol2_info = res_idx

close(f_unit)

end subroutine read_aa_to_mol2
!-------------------------------------------------------------------------------
subroutine find_mol2_res_idx(res_name, res_type, residue_mol2_info, num_mol2_info)
!-------------------------------------------------------------------------------
character(len=4), intent(inout) :: res_name
integer, intent(out) :: res_type
type(residue_mol2_type), intent(in) :: residue_mol2_info(:)
integer, intent(in) :: num_mol2_info
integer :: i_res

res_type = 0

if (trim(res_name) == 'HID' .or. trim(res_name) == 'HIE'&
    .or. trim(res_name) == 'HIP') then
    res_name = 'HIS'
else if (trim(res_name) == 'CYM' .or. trim(res_name) == 'CYX') then
    res_name = 'CYS'
else if (trim(res_name) == 'CPR') then
    res_name = 'PRO'
else if (trim(res_name) == 'ZN' .or. trim(res_name) == 'LI') then
    res_name = 'HET'
else if (trim(res_name) == 'NA' .or. trim(res_name) == 'K') then
    res_name = 'HET'
else if (trim(res_name) == 'CA' .or. trim(res_name) == 'MG') then
    res_name = 'HET'
else if (trim(res_name) == 'AL' .or. trim(res_name) == 'MN') then
    res_name = 'HET'
else if (trim(res_name) == 'FE' .or. trim(res_name) == 'NI') then
    res_name = 'HET'
else if (trim(res_name) == 'CD' .or. trim(res_name) == 'CO') then
    res_name = 'HET'
else if (trim(res_name) == 'CU' .or. trim(res_name) == 'HG') then
    res_name = 'HET'
else if (trim(res_name) == 'F' .or. trim(res_name) == 'CL') then
    res_name = 'HET'
else if (trim(res_name) == 'BR' .or. trim(res_name) == 'I') then
    res_name = 'HET'
end if

do i_res = 1, num_mol2_info
    if (trim(res_name) == trim(residue_mol2_info(i_res)%res_name)) then
        res_type = i_res
        return
    end if
end do

write(log_msg,'(A,A)') 'ERROR: Failed to map mol2 type of this residue, ', res_name
call terminate_with_error(log_msg)

end subroutine find_mol2_res_idx
!-------------------------------------------------------------------------------
subroutine find_mol2_atm_idx(atom_name, res_type, atm_type, &
                             residue_mol2_info)
!-------------------------------------------------------------------------------
character(len=4), intent(inout) :: atom_name
integer, intent(in) :: res_type
integer, intent(out) :: atm_type
type(residue_mol2_type), intent(in) :: residue_mol2_info(:)
character(len=4) :: res_name
integer :: i_atm

atm_type = 0
res_name = residue_mol2_info(res_type)%res_name

if (trim(res_name) == 'ACE' .or. trim(res_name) == 'NME') then
   if(trim(atom_name) == 'CH3') then
      atom_name = 'CA'
   else if(trim(atom_name) == 'HH31') then
      atom_name = 'HA1'
   else if(trim(atom_name) == 'HH32') then
      atom_name = 'HA2'
   else if(trim(atom_name) == 'HH33') then
      atom_name = 'HA3'
   end if
! It's for considering different atom name between amber_topo and
! residue_mol2_info.
else if(trim(res_name) /= 'ALA' .and. trim(res_name) /= 'VAL' &
        .and. trim(res_name) /= 'THR') then
   if(trim(atom_name) == 'HA3') then
      atom_name = 'HA1'
   else if(trim(atom_name) == 'HB3') then
      atom_name = 'HB1'
   else if(trim(atom_name) == 'HG3') then
      atom_name = 'HG1'
   else if(trim(atom_name) == 'HD3') then
      atom_name = 'HD1'
   else if(trim(res_name) == 'LYS' .and. trim(atom_name) == 'HE3') then
      atom_name = 'HE1'
   else if(trim(atom_name) == 'HG13') then
      atom_name = 'HG11'
   end if
end if

do i_atm = 1, residue_mol2_info(res_type)%n_atm 
   if (trim(atom_name) == trim(residue_mol2_info(res_type)%atom_name(i_atm))) then
      atm_type = i_atm
      return
   end if
end do

write(log_msg,'(A,A,A,A)') 'ERROR: Failed to map mol2 type of this atom, ', atom_name, &
                           ' in this residue, ', res_name
call terminate_with_error(log_msg)

end subroutine find_mol2_atm_idx
!-------------------------------------------------------------------------------
subroutine check_is_new_type(type_list, n_type, atom_type, is_new, type_idx)
!-------------------------------------------------------------------------------
! check atom_type is in type_list or not.
! If atom_type is in type_list, then return index of atom_type in type_list.
! Otherwise, return is_new=.true.
!-------------------------------------------------------------------------------
integer, intent(in) :: type_list(:), atom_type
integer, intent(in) :: n_type
logical, intent(out) :: is_new
integer, intent(out) :: type_idx
integer :: i_type

is_new = .true.
type_idx = -1
do i_type = 1, n_type
    if (atom_type == type_list(i_type)) then
        is_new = .false.
        type_idx = i_type
        return
    end if
end do

end subroutine check_is_new_type
!-------------------------------------------------------------------------------
subroutine check_is_new_type_string(type_list, n_type, atom_type, is_new, type_idx)
!-------------------------------------------------------------------------------
! string version of check_is_new_type
!-------------------------------------------------------------------------------
character(len=2), intent(in) :: type_list(:), atom_type
integer, intent(in) :: n_type
logical, intent(out) :: is_new
integer, intent(out) :: type_idx
integer :: i_type

is_new = .true.
type_idx = -1
do i_type = 1, n_type
    if (atom_type == type_list(i_type)) then
        is_new = .false.
        type_idx = i_type
        return
    end if
end do

end subroutine check_is_new_type_string
!-------------------------------------------------------------------------------
subroutine copy_to_grid(E, n_elem, grid, x, y, z, n_hash)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: E
integer, intent(in) :: n_elem(3)
real(dp), intent(inout) :: grid(:,:)
integer, intent(in) :: x, y, z, n_hash
integer :: hash_org, hash, idx, diff
integer :: i_y, i_z

hash_org = n_elem(2)*n_elem(1)*(z-1) + n_elem(1)*(y-1) + x
do i_z = z-1, z
    do i_y = y-1, y
        hash = n_elem(2)*n_elem(1)*(i_z-1) + n_elem(1)*(i_y-1) + x
        if (hash < 1) cycle
        if (hash > n_hash) cycle
        diff = hash_org - hash
        if (diff == 0) then
            idx = 1
        else if (diff == n_elem(1)) then
            idx = 2
        else if (diff == n_elem(2)*n_elem(1)) then
            idx = 3
        else
            idx = 4
        end if
        grid(idx, hash) = E
    end do
end do

end subroutine copy_to_grid
!-------------------------------------------------------------------------------
subroutine calc_E_using_grid(coord, grid_info, dock_grid, E, g, calc_g)
!!-------------------------------------------------------------------------------
real(dp), intent(in) :: coord(3)
type(docking_grid_param_type), intent(in) :: grid_info
real(dp), intent(in) :: dock_grid(:,:)
real(dp), intent(out) :: E
real(dp) :: dr(3), tcrd(3)
real(dp) :: grid_value(8), Xtmp(3)
integer :: i_crd, X0(3), hash
real(dp), intent(out) :: g(3)
logical, intent(in) :: calc_g

dr(:) = coord(:) - grid_info%grid_cntr(:)
tcrd(:) = dr(:)/grid_info%grid_width
tcrd(:) = tcrd(:) + (dble(grid_info%n_elem(:)-1)*0.5d0 + 1.0d0)

!TODO: using IFDEF (dble/real)
do i_crd = 1, 3
   if (tcrd(i_crd) < 1.0d0 .or. tcrd(i_crd) > dble(grid_info%n_elem(i_crd))) then
      E = dot_product(dr,dr)*500.0
      g(:) = 1000.0d0*dr(:)
      return
   end if
end do

X0(:) = int(tcrd(:))
hash = grid_info%n_elem(2)*grid_info%n_elem(1)*(X0(3)-1) &
     + grid_info%n_elem(1)*(X0(2)-1) + X0(1)
grid_value(1:4) = dock_grid(1:4,hash)
grid_value(5:8) = dock_grid(1:4,hash+1)

!TODO: using IFDEF (dble/real)
Xtmp = dble(X0)
E = 0.0d0
!call trilinterp(tcrd, Xtmp, grid_value, 1.0d0, E)
call trilinterp(tcrd, Xtmp, grid_value, 1.0d0, E, g, calc_g)
if (calc_g) then
    g(:) = g(:) * 1.0d0 / grid_info%grid_width
end if

end subroutine calc_E_using_grid
!-------------------------------------------------------------------------------
END MODULE LIGDOCK_E_UTILS
!-------------------------------------------------------------------------------
