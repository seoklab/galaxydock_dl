!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/vars/allocate_molecule.f90
!
! Description:
!  This file deals molecule_type-related array sizes. It allocates/deallocates
!   array related to molecule_type or simple_molecule_type. 
!-------------------------------------------------------------------------------
MODULE ALLOCATE_MOLECULE
!-------------------------------------------------------------------------------
use globals

implicit none
public

CONTAINS
!===============================================================================
! MOLECULE_TYPE-related
!===============================================================================
subroutine allocate_molecule_type(molecule, n_res, n_het, n_lig)
!-------------------------------------------------------------------------------
! Allocate molecule_type with default size
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer, intent(in), optional :: n_res, n_het, n_lig

if (allocated(molecule%residue)) return

if (.not. present(n_res)) then
    allocate(molecule%residue(max_res))
    molecule%n_res = max_res
else
    allocate(molecule%residue(n_res))
    molecule%n_res = n_res
end if
!
if (.not. present(n_het)) then
    allocate(molecule%hetmol(max_het))
    molecule%n_het = max_het
else
    allocate(molecule%hetmol(n_het))
    molecule%n_het = n_het
end if
!
if (.not. present(n_lig)) then
    allocate(molecule%ligand(max_lig))
    molecule%n_lig = max_lig
else
    allocate(molecule%ligand(n_lig))
    molecule%n_lig = n_lig
end if

end subroutine allocate_molecule_type
!-------------------------------------------------------------------------------
subroutine deallocate_molecule_type(molecule)
!-------------------------------------------------------------------------------
! Deallocate the molecule_type variable
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule

deallocate(molecule%residue)
deallocate(molecule%hetmol)
deallocate(molecule%ligand)
!
molecule%n_res = 0
molecule%n_het = 0
molecule%n_lig = 0

end subroutine deallocate_molecule_type
!-------------------------------------------------------------------------------
subroutine allocate_residue_type(residue, n_atm)
!-------------------------------------------------------------------------------
! Allocate residue_type with given n_atm
!-------------------------------------------------------------------------------
type(residue_type), intent(inout) :: residue
integer, intent(in) :: n_atm

! Initializes atm_read for reading "read_pdb@in_out_structure.f90"
!  only if the variable residue is already allocated
if (allocated(residue%R)) then
    residue%atm_read(:) = .false.
    return
end if

residue%n_atm = n_atm
!
allocate(residue%pdb_atom_name(n_atm))
!
! Atom geometry information 
allocate(residue%R(3, 0:n_atm))
allocate(residue%b(3,-1:n_atm))
allocate(residue%quat(4, 0:n_atm))
allocate(residue%b_len(n_atm))
allocate(residue%b_ang(n_atm))
allocate(residue%t_ang(n_atm))
!
! Atom misc information
allocate(residue%atm_read(n_atm))
allocate(residue%atm_placed(n_atm))
allocate(residue%bnd_read(-1:n_atm))
allocate(residue%atom_fixed(n_atm))
!allocate(residue%atm_in_bnd_E(2,n_atm))
!
residue%atm_read(:) = .false.
residue%atm_placed(:) = .false.
residue%chain = ''
residue%code = ''
residue%res_added = ''

end subroutine allocate_residue_type
!-------------------------------------------------------------------------------
subroutine reallocate_residue_type(residue, n_atm)
!-------------------------------------------------------------------------------
! Changes the size of residue varaible to n_atm
! This subroutine MUST be tested before use it; it can be used for mutation, ...
!-------------------------------------------------------------------------------
type(residue_type), intent(inout) :: residue
integer, intent(in) :: n_atm
type(residue_type) :: residue_prev

! Allocate a temporary variable with n_atm
call allocate_residue_type(residue_prev, n_atm)
! Copy the whole info into the temporary variable
residue_prev = residue
!
residue%n_atm = n_atm
call move_alloc(residue_prev%pdb_atom_name, residue%pdb_atom_name)
!
call move_alloc(residue_prev%R,     residue%R)
call move_alloc(residue_prev%b,     residue%b)
call move_alloc(residue_prev%quat,  residue%quat)
call move_alloc(residue_prev%b_len, residue%b_len)
call move_alloc(residue_prev%b_ang, residue%b_ang)
call move_alloc(residue_prev%t_ang, residue%t_ang)
!
call move_alloc(residue_prev%atm_read,     residue%atm_read)
call move_alloc(residue_prev%atm_placed,   residue%atm_placed)
call move_alloc(residue_prev%bnd_read,     residue%bnd_read)
call move_alloc(residue_prev%atom_fixed,   residue%atom_fixed)
!call move_alloc(residue_prev%atm_in_bnd_E, residue%atm_in_bnd_E)

end subroutine reallocate_residue_type
!-------------------------------------------------------------------------------
subroutine deallocate_residue_type(residue)
!-------------------------------------------------------------------------------
type(residue_type), intent(inout) :: residue

residue%n_atm = 0
!
deallocate(residue%pdb_atom_name)
!
! Atom geometry information 
deallocate(residue%R)
deallocate(residue%b)
deallocate(residue%quat)
deallocate(residue%b_len)
deallocate(residue%b_ang)
deallocate(residue%t_ang)
!
! Atom misc information
deallocate(residue%atm_read)
deallocate(residue%atm_placed)
deallocate(residue%bnd_read)
deallocate(residue%atom_fixed)
!deallocate(residue%atm_in_bnd_E)

end subroutine deallocate_residue_type
!-------------------------------------------------------------------------------
subroutine allocate_ligand_residue_type(residue, n_atm)
!-------------------------------------------------------------------------------
type(ligand_residue_type), intent(inout) :: residue
integer, intent(in) :: n_atm

if (allocated(residue%R)) then
    residue%atm_read(:) = .false.
    return
end if

residue%n_atm = n_atm
!
allocate(residue%pdb_atom_name(n_atm))
!
! Atom geometry information 
allocate(residue%R(3, 0:n_atm))
allocate(residue%b(3,-1:n_atm))
allocate(residue%quat(4, 0:n_atm))
allocate(residue%b_len(n_atm*2))
allocate(residue%b_ang(n_atm*2))
allocate(residue%t_ang(n_atm*4))
!
! Atom misc information
allocate(residue%atm_read(n_atm))
allocate(residue%atm_placed(n_atm))
allocate(residue%bnd_read(-1:n_atm))
allocate(residue%atom_fixed(n_atm))
allocate(residue%atm_in_bnd_E(2,n_atm))
!
residue%atm_placed(:) = .false.
residue%chain = ''
residue%code = ''

end subroutine allocate_ligand_residue_type
!-------------------------------------------------------------------------------
subroutine reallocate_ligand_residue_type(residue, n_atm)
!-------------------------------------------------------------------------------
type(ligand_residue_type), intent(inout) :: residue
integer, intent(in) :: n_atm
type(ligand_residue_type) :: residue_prev

call allocate_ligand_residue_type(residue_prev, n_atm)
residue_prev = residue
!
residue%n_atm = n_atm
call move_alloc(residue_prev%pdb_atom_name, residue%pdb_atom_name)
!
call move_alloc(residue_prev%R,     residue%R)
call move_alloc(residue_prev%b,     residue%b)
call move_alloc(residue_prev%quat,  residue%quat)
call move_alloc(residue_prev%b_len, residue%b_len)
call move_alloc(residue_prev%b_ang, residue%b_ang)
call move_alloc(residue_prev%t_ang, residue%t_ang)
!
call move_alloc(residue_prev%atm_read,     residue%atm_read)
call move_alloc(residue_prev%atm_placed,   residue%atm_placed)
call move_alloc(residue_prev%bnd_read,     residue%bnd_read)
call move_alloc(residue_prev%atom_fixed,   residue%atom_fixed)
call move_alloc(residue_prev%atm_in_bnd_E, residue%atm_in_bnd_E)

end subroutine reallocate_ligand_residue_type
!-------------------------------------------------------------------------------
subroutine deallocate_ligand_residue_type(residue)
!-------------------------------------------------------------------------------
type(ligand_residue_type), intent(inout) :: residue

residue%n_atm = 0
!
deallocate(residue%pdb_atom_name)
!
! Atom geometry information 
deallocate(residue%R)
deallocate(residue%b)
deallocate(residue%quat)
deallocate(residue%b_len)
deallocate(residue%b_ang)
deallocate(residue%t_ang)
!
! Atom misc information
deallocate(residue%atm_read)
deallocate(residue%atm_placed)
deallocate(residue%bnd_read)
deallocate(residue%atom_fixed)
deallocate(residue%atm_in_bnd_E)

end subroutine deallocate_ligand_residue_type
!-------------------------------------------------------------------------------
subroutine reallocate_residue_in_molecule(molecule, n_res)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer, intent(in) :: n_res
type(residue_type), allocatable :: residue(:)

allocate(residue(n_res))
residue(1:min(n_res,molecule%n_res)) = molecule%residue(1:min(n_res,molecule%n_res))
!
call move_alloc(residue, molecule%residue)
molecule%n_res = n_res

end subroutine reallocate_residue_in_molecule
!-------------------------------------------------------------------------------
subroutine reallocate_hetmol_in_molecule(molecule, n_het)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer, intent(in) :: n_het
type(residue_type), allocatable :: hetmol(:)

allocate(hetmol(n_het))
hetmol(1:min(n_het,molecule%n_het)) = molecule%hetmol(1:min(n_het,molecule%n_het))
!
call move_alloc(hetmol, molecule%hetmol)
molecule%n_het = n_het

end subroutine reallocate_hetmol_in_molecule
!-------------------------------------------------------------------------------
subroutine reallocate_ligand_in_molecule(molecule, n_lig)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer, intent(in) :: n_lig
type(ligand_residue_type), allocatable :: ligand(:)

allocate(ligand(n_lig))
ligand(1:min(n_lig,molecule%n_lig)) = molecule%ligand(1:min(n_lig,molecule%n_lig))
!
call move_alloc(ligand, molecule%ligand)
molecule%n_lig = n_lig

end subroutine reallocate_ligand_in_molecule
!===============================================================================
! SIMPLE_MOLECULE_TYPE-related
!===============================================================================
subroutine allocate_simple_molecule_type(molecule, n_res, n_het, n_lig)
!-------------------------------------------------------------------------------
! Allocate simple_molecule_type with given size
!-------------------------------------------------------------------------------
type(simple_molecule_type), intent(inout) :: molecule
integer, intent(in), optional :: n_res, n_het, n_lig

if (allocated(molecule%residue)) return

if (.not. present(n_res)) then
    allocate(molecule%residue(max_res))
else
    allocate(molecule%residue(n_res))
end if
!
if (.not. present(n_het)) then
    allocate(molecule%hetmol(max_het))
else
    allocate(molecule%hetmol(n_het))
end if
!
if (.not. present(n_lig)) then
    allocate(molecule%ligand(max_lig))
else
    allocate(molecule%ligand(n_lig))
end if

end subroutine allocate_simple_molecule_type
!-------------------------------------------------------------------------------
subroutine deallocate_simple_molecule_type(molecule)
!-------------------------------------------------------------------------------
! Deallocate the simple_molecule_type variable
!-------------------------------------------------------------------------------
type(simple_molecule_type), intent(inout) :: molecule

deallocate(molecule%residue)
deallocate(molecule%hetmol)
deallocate(molecule%ligand)

end subroutine deallocate_simple_molecule_type
!-------------------------------------------------------------------------------
subroutine allocate_simple_residue_type(residue, n_atm)
!-------------------------------------------------------------------------------
! Allocate residue_type with given n_atm
!-------------------------------------------------------------------------------
type(simple_residue_type), intent(inout) :: residue
integer, intent(in) :: n_atm

residue%n_atm = n_atm
!
! Atom geometry information 
allocate(residue%R(3, 0:n_atm))
allocate(residue%b_len(n_atm))
allocate(residue%b_ang(n_atm))
allocate(residue%t_ang(n_atm))

end subroutine allocate_simple_residue_type
!-------------------------------------------------------------------------------
subroutine reallocate_simple_residue_type(residue, n_atm)
!-------------------------------------------------------------------------------
! Changes the size of residue varaible to n_atm
!-------------------------------------------------------------------------------
type(simple_residue_type), intent(inout) :: residue
integer, intent(in) :: n_atm
type(simple_residue_type) :: residue_prev

! Allocate a temporary variable with n_atm
call allocate_simple_residue_type(residue_prev, n_atm)
! Copy the whole info into the temporary variable
residue_prev = residue
!
residue%n_atm = n_atm
call move_alloc(residue_prev%R,     residue%R)
call move_alloc(residue_prev%b_len, residue%b_len)
call move_alloc(residue_prev%b_ang, residue%b_ang)
call move_alloc(residue_prev%t_ang, residue%t_ang)

end subroutine reallocate_simple_residue_type
!-------------------------------------------------------------------------------
subroutine deallocate_simple_residue_type(residue)
!-------------------------------------------------------------------------------
type(simple_residue_type), intent(inout) :: residue

residue%n_atm = 0
!
! Atom geometry information 
deallocate(residue%R)
deallocate(residue%b_len)
deallocate(residue%b_ang)
deallocate(residue%t_ang)

end subroutine deallocate_simple_residue_type
!-------------------------------------------------------------------------------
subroutine allocate_simple_ligand_residue_type(residue, n_atm)
!-------------------------------------------------------------------------------
type(simple_ligand_residue_type), intent(inout) :: residue
integer, intent(in) :: n_atm

residue%n_atm = n_atm
!
! Atom geometry information 
allocate(residue%R(3, 0:n_atm))
allocate(residue%b_len(n_atm*2))
allocate(residue%b_ang(n_atm*2))
allocate(residue%t_ang(n_atm*3))

end subroutine allocate_simple_ligand_residue_type
!-------------------------------------------------------------------------------
subroutine reallocate_simple_ligand_residue_type(residue, n_atm)
!-------------------------------------------------------------------------------
type(simple_ligand_residue_type), intent(inout) :: residue
integer, intent(in) :: n_atm
type(simple_ligand_residue_type) :: residue_prev

call allocate_simple_ligand_residue_type(residue_prev, n_atm)
residue_prev = residue
!
residue%n_atm = n_atm
call move_alloc(residue_prev%R,     residue%R)
call move_alloc(residue_prev%b_len, residue%b_len)
call move_alloc(residue_prev%b_ang, residue%b_ang)
call move_alloc(residue_prev%t_ang, residue%t_ang)

end subroutine reallocate_simple_ligand_residue_type
!-------------------------------------------------------------------------------
subroutine deallocate_simple_ligand_residue_type(residue)
!-------------------------------------------------------------------------------
type(simple_ligand_residue_type), intent(inout) :: residue

residue%n_atm = 0
!
! Atom geometry information 
deallocate(residue%R)
deallocate(residue%b_len)
deallocate(residue%b_ang)
deallocate(residue%t_ang)

end subroutine deallocate_simple_ligand_residue_type
!-------------------------------------------------------------------------------
subroutine update_molecule_coordinates(prot_1, prot_2)
!-------------------------------------------------------------------------------
! Copy coord info in the prot_1 to prot_2
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: prot_1
type(molecule_type), intent(inout) :: prot_2
!
integer :: i_res, n

do i_res = 1, prot_1%n_res
    n = prot_1%residue(i_res)%n_atm
    !
    prot_2%residue(i_res)%R(1:3, 0:n) = prot_1%residue(i_res)%R(1:3, 0:n)   
    prot_2%residue(i_res)%b(1:3,-1:n) = prot_1%residue(i_res)%b(1:3,-1:n)   
    prot_2%residue(i_res)%quat(1:4, 0:n) = prot_1%residue(i_res)%quat(1:4, 0:n)
    prot_2%residue(i_res)%b_len(1:n) = prot_1%residue(i_res)%b_len(1:n)    
    prot_2%residue(i_res)%b_ang(1:n) = prot_1%residue(i_res)%b_ang(1:n)    
    prot_2%residue(i_res)%t_ang(1:n) = prot_1%residue(i_res)%t_ang(1:n)    
end do

do i_res = 1, prot_1%n_het
    n = prot_1%hetmol(i_res)%n_atm
    !
    prot_2%hetmol(i_res)%R(1:3, 0:n) = prot_1%hetmol(i_res)%R(1:3, 0:n)   
    prot_2%hetmol(i_res)%b(1:3,-1:n) = prot_1%hetmol(i_res)%b(1:3,-1:n)   
    prot_2%hetmol(i_res)%quat(1:4, 0:n) = prot_1%hetmol(i_res)%quat(1:4, 0:n)
    prot_2%hetmol(i_res)%b_len(1:n) = prot_1%hetmol(i_res)%b_len(1:n)    
    prot_2%hetmol(i_res)%b_ang(1:n) = prot_1%hetmol(i_res)%b_ang(1:n)    
    prot_2%hetmol(i_res)%t_ang(1:n) = prot_1%hetmol(i_res)%t_ang(1:n)    
end do

do i_res = 1, prot_1%n_lig
    n = prot_1%ligand(i_res)%n_atm
    !
    prot_2%ligand(i_res)%R(1:3, 0:n) = prot_1%ligand(i_res)%R(1:3, 0:n)   
    prot_2%ligand(i_res)%b(1:3,-1:n) = prot_1%ligand(i_res)%b(1:3,-1:n)   
    prot_2%ligand(i_res)%quat(1:4, 0:n) = prot_1%ligand(i_res)%quat(1:4, 0:n)
    prot_2%ligand(i_res)%b_len(1:n) = prot_1%ligand(i_res)%b_len(1:n)    
    prot_2%ligand(i_res)%b_ang(1:n) = prot_1%ligand(i_res)%b_ang(1:n)    
    prot_2%ligand(i_res)%t_ang(1:n) = prot_1%ligand(i_res)%t_ang(1:n)    
end do

end subroutine update_molecule_coordinates
!-------------------------------------------------------------------------------
subroutine update_residue_coordinates(res_1, res_2)
!-------------------------------------------------------------------------------
type(residue_type), intent(in) :: res_1
type(residue_type), intent(inout) :: res_2
!
integer :: n

n = res_1%n_atm
!
res_2%R(1:3, 0:n) = res_1%R(1:3, 0:n)   
res_2%b(1:3,-1:n) = res_1%b(1:3,-1:n)   
res_2%quat(1:4, 0:n) = res_1%quat(1:4, 0:n)
res_2%b_len(1:n) = res_1%b_len(1:n)    
res_2%b_ang(1:n) = res_1%b_ang(1:n)    
res_2%t_ang(1:n) = res_1%t_ang(1:n)    

end subroutine update_residue_coordinates
!-------------------------------------------------------------------------------
END MODULE ALLOCATE_MOLECULE
!-------------------------------------------------------------------------------
