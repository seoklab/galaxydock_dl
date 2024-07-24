!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/utils/symmetry.f90
!
! Description: Subroutines used for symmetric modeling.
!
!-------------------------------------------------------------------------------
MODULE SYMMETRY
!-------------------------------------------------------------------------------
use globals
use allocate_molecule 
use logger,       only: terminate_with_error, log_p
use rmsd,         only: ls_rmsd
use in_out_utils, only: find_atom_idx
use geometry,     only: cartesian2internal

implicit none
save
private

integer, parameter :: max_symm = 12     ! Max no. symmetry
!
integer :: n_symm                      ! # of symmetric unit in molecule
integer  :: symm_resrange(2,max_symm)  ! residue range for each unit 
real(dp) :: symm_matrix(3,3,max_symm)  ! Rotation matrix for each unit
real(dp) :: symm_center(3,0:max_symm)  ! Center for each unit
!
! Variables used to set symmetry with REMARK 350 
real(dp) :: symm_U(3,3,max_symm)      
real(dp) :: symm_T(3,max_symm)
!
integer :: n_chain_symm_unit
character(len=1) :: symm_chain(max_symm)

public :: max_symm
public :: n_symm
public :: symm_resrange
!
public :: symm_U
public :: symm_T
public :: n_chain_symm_unit
public :: symm_chain
!
public :: setup_symmetry
!
public :: set_USC_for_symm
public :: set_ULR_for_symm
public :: apply_fix_type_symm
public :: build_appl_res_symm
!
public :: copy_symm_unit
public :: copy_symm_unit_molecule
!
public :: write_symm_info
!
public :: find_residue_in_first_symm_unit

CONTAINS
!===============================================================================
! Subroutines to set symmetry matrix
!===============================================================================
subroutine setup_symmetry(protein)
!-------------------------------------------------------------------------------
! Initialize symmetric operation
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
integer :: n_res

if (use_remark350) then
    if (protein%n_chain /= n_chain_symm_unit) then
        write(log_msg, '(A)') 'ERROR: Failed to set up symmetry. Chains in REMARK 350 and chains in PDB are different.'
        call terminate_with_error(log_msg)
    end if
    call setup_symmetry_using_REMARK(protein)
    n_res = protein%n_res
    call reallocate_residue_in_molecule(protein, n_res*n_symm)
    call copy_symm_unit_res_info(protein)
    protein%n_chain = n_symm * n_chain_symm_unit
else
    call setup_symmetry_using_PDB(protein)
end if

end subroutine setup_symmetry
!-------------------------------------------------------------------------------
subroutine setup_symmetry_using_REMARK(protein)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer :: i_symm, i, i_res, n_res
real(dp), allocatable :: ref_crd(:,:), crd_n(:,:)
logical :: is_symm_chain

write(log_msg, '(A,I3,A)') "Setting symmetry with ", n_symm,&
                           " symm_units using REMARK 350."
call log_p(log_msg, level=10, me=me)

! Check chainID in input PDB is matched to that in REMARK 350
do i_res = 1, protein%n_res
    if (protein%residue(i_res)%ter_type /= 'N') cycle
    do i = 1, n_chain_symm_unit
        if (protein%residue(i_res)%chain == symm_chain(i)) then
            is_symm_chain = .true.
            exit
        end if
    end do
    if (.not. is_symm_chain) then
        call terminate_with_error('ERROR: Failed to set up symmetry. Please check chain ID in REMARK 350 and ATOM line')
    end if
end do

! setup symm_resrange
do i_symm = 1, n_symm
    symm_resrange(1,i_symm) = protein%n_res*(i_symm-1) + 1
    symm_resrange(2,i_symm) = protein%n_res*(i_symm-1) + protein%n_res
end do

symm_matrix(:,:,:) = 0.0d0
symm_matrix(:,:,:) = symm_U(:,:,:)
symm_center(:,0:n_symm) = 0.0d0

n_res = protein%n_res
allocate(ref_crd(3,n_res), crd_n(3,n_res))

call copy_CA_atm(protein, ref_crd, symm_resrange(1,1), symm_resrange(2,1))
do i = 1, 3
    symm_center(i,1) = sum(ref_crd(i,:))/dble(n_res)
end do

do i_symm = 2, n_symm
    do i = 1, protein%n_res
        crd_n(:,i) = matmul(symm_U(:,:,i_symm), ref_crd(:,i)) + symm_T(:,i_symm)
    end do
    do i = 1, 3
        symm_center(i, i_symm) = sum(crd_n(i,:))/dble(n_res)
    end do
end do
deallocate(ref_crd, crd_n)

symm_center(1,0) = sum(symm_center(1,1:n_symm)) / dble(n_symm)
symm_center(2,0) = sum(symm_center(2,1:n_symm)) / dble(n_symm)
symm_center(3,0) = sum(symm_center(3,1:n_symm)) / dble(n_symm)

end subroutine setup_symmetry_using_REMARK
!-------------------------------------------------------------------------------
subroutine setup_symmetry_using_PDB(protein)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer :: i, i_res, i_symm, n_res
real(dp), allocatable :: crd_1(:,:), crd_n(:,:), ref_crd(:,:)
real(dp) :: rmsd
character(len=1) :: prev_chain

if (n_symm == 0) call setup_symm_resrange(protein)

write(log_msg, '(A,I3,A)') "Setting symmetry with ", n_symm, " symm_units."
call log_p(log_msg, level=10, me=me)

symm_matrix(:,:,:) = 0.0d0
symm_matrix(1,:,1) = (/1.0d0, 0.0d0, 0.0d0/)
symm_matrix(2,:,1) = (/0.0d0, 1.0d0, 0.0d0/)
symm_matrix(3,:,1) = (/0.0d0, 0.0d0, 1.0d0/)
!
symm_center(:,0:n_symm) = 0.0d0

!Check if has identical residues within symmetric units
n_res = symm_resrange(2,1) - symm_resrange(1,1) + 1 ! #of resiudes in symm. unit.
do i_symm = 2, n_symm
    if (symm_resrange(2,i_symm) - symm_resrange(1,i_symm) + 1 /= n_res) then
        call terminate_with_error('ERROR: Failed to set up symmetry.')
    end if
end do

! Fill symm_chain which used to write output PDB files.
prev_chain = ''
n_chain_symm_unit = 0
do i_res = symm_resrange(1,1), symm_resrange(2,1)
    if (protein%residue(i_res)%chain /= prev_chain) then
        n_chain_symm_unit = n_chain_symm_unit + 1
        symm_chain(n_chain_symm_unit) = protein%residue(i_res)%chain
    end if
    prev_chain = protein%residue(i_res)%chain
end do

! Get center for each unit, coordinate system, rotation matrix
allocate(crd_1(3,n_res), crd_n(3,n_res), ref_crd(3,n_res))
! copy 1st symmetric unit CA atm to ref_crd
call copy_CA_atm(protein, ref_crd, symm_resrange(1,1), symm_resrange(2,1))
do i = 1, 3
    symm_center(i,1) = sum(ref_crd(i,:))/dble(n_res)
end do

do i_symm = 2, n_symm ! Iter over symmetric units
    ! copy n-th symmetric unit CA atm to crd_n
    crd_1(:,:) = ref_crd(:,:)
    call copy_CA_atm(protein, crd_n, symm_resrange(1,i_symm), symm_resrange(2,i_symm))
    do i = 1, 3
        symm_center(i,i_symm) = sum(crd_n(i,:))/dble(n_res)
    end do
    ! Get best-fitting rotation matrix with SVD
    call ls_rmsd(n_res, crd_1, crd_n, 1, symm_matrix(:,:,i_symm), rmsd)
end do
deallocate(crd_1, crd_n, ref_crd)

if (mod(tn%stdatm, n_symm)/= 0) then
    call terminate_with_error('ERROR: Failed to set up symmetry.')
end if
  
symm_center(1,0) = sum(symm_center(1,1:n_symm)) / dble(n_symm)
symm_center(2,0) = sum(symm_center(2,1:n_symm)) / dble(n_symm)
symm_center(3,0) = sum(symm_center(3,1:n_symm)) / dble(n_symm)

end subroutine setup_symmetry_using_PDB
!-------------------------------------------------------------------------------
subroutine setup_symm_resrange(protein)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer :: i_res, res_no_N, res_no_C

res_no_N = 9999
res_no_C = -9999
do i_res = 1, protein%n_res
    if (protein%residue(i_res)%ter_type == 'N') then
        res_no_N = i_res
    elseif (protein%residue(i_res)%ter_type == 'C') then
        res_no_C = i_res
        if (res_no_C < res_no_N) call terminate_with_error("Error. Setup symmetry units.")
        n_symm = n_symm + 1
        symm_resrange(1,n_symm) = res_no_N
        symm_resrange(2,n_symm) = res_no_C
    end if
end do

end subroutine setup_symm_resrange
!-------------------------------------------------------------------------------
subroutine copy_CA_atm(protein,crd,start_res,end_res)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer, intent(in) :: start_res,end_res
real(dp),intent(out) :: crd(:,:)
integer :: res_no, i_res, ca_atm
character(len=6) :: error_mode !for find_atom_idx. if the value is NOT 'ignore',
!                              !it will prints error when atom/residue type is
!                              !unknown
res_no = 0
crd(:,:) = 0.0d0

do i_res = start_res, end_res
    res_no = res_no + 1
    call find_atom_idx(i_res, protein%residue(i_res)%res_type, 'CA  ', ca_atm, error_mode) 
    crd(:,res_no) = protein%residue(i_res)%R(:,ca_atm)
end do

end subroutine copy_CA_atm
!-------------------------------------------------------------------------------
subroutine copy_symm_unit_res_info(molecule)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer :: i_symm, i_res, res_no, i_atm, n_atm

do i_symm = 2, n_symm
    do i_res = symm_resrange(1,1), symm_resrange(2,1)
        res_no = i_res + symm_resrange(2,i_symm-1)
        n_atm = molecule%residue(i_res)%n_atm
        !
        call reallocate_residue_type(molecule%residue(res_no), n_atm)
        !
        molecule%residue(res_no) = molecule%residue(i_res)
        do i_atm = 1, n_atm
            molecule%residue(res_no)%R(:,i_atm) &
                = matmul(symm_U(:,:,i_symm),molecule%residue(i_res)%R(:,i_atm)) &
                + symm_T(:,i_symm)
        end do
    end do
end do

end subroutine copy_symm_unit_res_info
!===============================================================================
! Subroutines to apply symmetric operations
!===============================================================================
subroutine center_symm(n_res, R)
!-------------------------------------------------------------------------------
real(dp),intent(in) :: R(:,:)
integer, intent(in) :: n_res
real(dp) :: U(3,3), vec(3)
integer :: i_symm, i_res

! Calculate symmetry center
symm_center(:,1) = 0.0d0
do i_res = 1, n_res
    symm_center(:,1) = symm_center(:,1) + R(:,res_index(i_res)%Ca_id(1))
end do
symm_center(:,1) = symm_center(:,1) / dble(n_res)

vec(:) = symm_center(:,1) - symm_center(:,0) 
do i_symm = 2, n_symm
    U(:,:) = symm_matrix(:,:,i_symm)
    symm_center(:,i_symm) = matmul(U,vec(:)) + symm_center(:,0)
end do

end subroutine center_symm
!-------------------------------------------------------------------------------
subroutine center_symm_molecule(protein, n_res)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer, intent(in) :: n_res
real(dp) :: ref_crd(3,n_res), U(3,3), vec(3)
integer :: i, i_symm 

! Calculate symmetry center
symm_center(:,1) = 0.0d0
call copy_CA_atm(protein, ref_crd, symm_resrange(1,1), symm_resrange(2,1))
do i = 1, 3
    symm_center(i,1) = sum(ref_crd(i,:))/dble(n_res)
end do

vec(:) = symm_center(:,1) - symm_center(:,0) 
do i_symm = 2, n_symm
    U(:,:) = symm_matrix(:,:,i_symm)
    symm_center(:,i_symm) = matmul(U,vec(:)) + symm_center(:,0)
end do

end subroutine center_symm_molecule
!-------------------------------------------------------------------------------
subroutine copy_symm_unit(R)
!-------------------------------------------------------------------------------
! Copy structure of first symmetric unit onto other symmetric units
! This needs to be called whenever structure is changed at reference unit
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: R(:,:)
  
real(dp) :: U(3,3), atmcrd(3)
integer :: i_atm, i_symm, i_res
integer :: n_res, symm_atm, atmno

n_res = symm_resrange(2,1) - symm_resrange(1,1) + 1
if (fix_type == FIX_NONE) call center_symm(n_res, R)
!call center_symm(n_res,R)

symm_atm = 0 ! # of atom in one symmetric unit 
do i_res = 1, n_res
    symm_atm = symm_atm + res_index(i_res)%n_atm
end do

! Then copy reference coordinate into other symmetric units
do i_symm = 2, n_symm
    U(:,:) = symm_matrix(:,:,i_symm)
    do i_atm = 1, symm_atm
        atmno = symm_atm*(i_symm-1)+i_atm
        atmcrd(:) = R(:,i_atm) - symm_center(:,1)
        R(:,atmno) = matmul(U,atmcrd(:)) + symm_center(:,i_symm)
    end do
end do
  
end subroutine copy_symm_unit
!-------------------------------------------------------------------------------
subroutine copy_symm_unit_molecule(protein)
!-------------------------------------------------------------------------------
! Copy structure of first symmetric unit onto other symmetric units
! This needs to be called whenever structure is changed at reference unit
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
real(dp) :: U(3,3), atmcrd(3)
integer :: i_atm, i_symm, i_res, res_no, n_res

n_res = symm_resrange(2,1) - symm_resrange(1,1) + 1
if (fix_type == FIX_NONE) call center_symm_molecule(protein, n_res)

! Then copy reference coordinate into other symmetric units
do i_symm = 2, n_symm
    U(:,:) = symm_matrix(:,:,i_symm)
    do i_res = symm_resrange(1,1), symm_resrange(2,1)
        res_no = i_res + symm_resrange(2, i_symm-1)
        do i_atm = 1, protein%residue(i_res)%n_atm
            atmcrd(:) = protein%residue(i_res)%R(:,i_atm) - symm_center(:,1)
            protein%residue(res_no)%R(:,i_atm) = &
                matmul(U,atmcrd(:)) + symm_center(:,i_symm)
        end do
    end do
end do
  
end subroutine copy_symm_unit_molecule
!===============================================================================
! Subroutines to set USC, ULR, appl_res symmetric.
!===============================================================================
subroutine set_USC_for_symm(is_usc, n_usc)
!-------------------------------------------------------------------------------
! Check current is_usc has any USC in other symmetric unit.
! If yes, set corresponding side-chain in first symmetric unit set to USC.
! After that, all side-chains in other symmetric units are turned off.
!-------------------------------------------------------------------------------
logical, intent(inout) :: is_usc(:)
integer, intent(inout) :: n_usc
integer :: i_symm, i_res, res_no

do i_symm = 2, n_symm
    do i_res = symm_resrange(1,1), symm_resrange(2,1)
        res_no = i_res + symm_resrange(2,i_symm-1)
        if (is_usc(res_no)) then 
            is_usc(i_res) = is_usc(res_no)
        end if
    end do
end do

is_usc(symm_resrange(1,2):symm_resrange(2,n_symm)) = .false.

n_usc = 0
do i_res = symm_resrange(1,1), symm_resrange(2,1)
    if (is_usc(i_res)) n_usc = n_usc + 1
end do

end subroutine set_USC_for_symm
!-------------------------------------------------------------------------------
subroutine set_ULR_for_symm(ULR)
!-------------------------------------------------------------------------------
! Check current ULR has any ULR in other symmetric unit.
! If yes, set corresponding residues in first symmetric unit set to ULR.
!-------------------------------------------------------------------------------
type(ulr_assign_type), intent(inout) :: ULR(:)
integer :: i_ulr

do i_ulr = 1, n_ulr
    call find_residue_in_first_symm_unit(ULR(i_ulr)%resrange(1))
    call find_residue_in_first_symm_unit(ULR(i_ulr)%resrange(2))
end do

end subroutine set_ULR_for_symm
!-------------------------------------------------------------------------------
subroutine apply_fix_type_symm(molecule)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer :: i_symm, i_res, res_no

do i_symm = 2, n_symm
    do i_res = symm_resrange(1,1), symm_resrange(2,1)
        res_no = i_res + symm_resrange(2,i_symm-1)
        molecule%residue(res_no)%atom_fixed(:) &
            = molecule%residue(i_res)%atom_fixed(:)
    end do
end do

end subroutine apply_fix_type_symm
!-------------------------------------------------------------------------------
subroutine build_appl_res_symm(appl_res)
!-------------------------------------------------------------------------------
logical, intent(inout) :: appl_res(:)
integer :: i_symm

do i_symm = 2, n_symm
    appl_res(symm_resrange(1,i_symm):symm_resrange(2,i_symm)) &
        = appl_res(symm_resrange(1,1):symm_resrange(2,1))
end do

end subroutine build_appl_res_symm
!===============================================================================
! Subroutines to write REMARK 350
!===============================================================================
subroutine write_symm_info(pdb_unit, molecule)
!-------------------------------------------------------------------------------
integer, intent(in) :: pdb_unit
type(molecule_type), intent(in) :: molecule
integer :: i_symm, i, n_res
real(dp) :: center(3,n_symm), rotated_A(3), rmsd
real(dp) :: matrix(3,3,n_symm)  ! Rotation matrix for each unit
real(dp), allocatable :: ref_crd(:,:), crd_1(:,:), crd_n(:,:)
integer :: i_chain
character(len=3) :: chain

n_res = symm_resrange(2,1) - symm_resrange(1,1) + 1

! Get center for each unit, coordinate system, rotation matrix
allocate(crd_1(3,n_res), crd_n(3,n_res), ref_crd(3,n_res))
! copy 1st symmetric unit CA atm to ref_crd
call copy_CA_atm(molecule, ref_crd, symm_resrange(1,1), symm_resrange(2,1))
do i = 1, 3
    center(i,1) = sum(ref_crd(i,:))/dble(n_res)
end do

matrix(1,:,1) = (/1.0d0, 0.0d0, 0.0d0/)
matrix(2,:,1) = (/0.0d0, 1.0d0, 0.0d0/)
matrix(3,:,1) = (/0.0d0, 0.0d0, 1.0d0/)
do i_symm = 2, n_symm ! Iter over symmetric units
    ! copy n-th symmetric unit CA atm to crd_n
    crd_1(:,:) = ref_crd(:,:)
    call copy_CA_atm(molecule, crd_n, symm_resrange(1,i_symm), symm_resrange(2,i_symm))
    do i = 1, 3
        center(i,i_symm) = sum(crd_n(i,:))/dble(n_res)
    end do
    ! Get best-fitting rotation matrix with SVD
    call ls_rmsd(n_res, crd_1, crd_n, 1, matrix(:,:,i_symm), rmsd)
end do

79 format(A,I1,I4,3F10.6,F15.5)
write(pdb_unit, '(A)') "REMARK 350 BIOMOLECULE: 1"
log_msg = "REMARK 350 APPLY THE FOLLOWING TO CHAINS:"
do i_chain = 1, n_chain_symm_unit - 1
    chain = ' '//symm_chain(i_chain)//','
    log_msg = trim(log_msg) // chain
end do
chain = ' '//symm_chain(n_chain_symm_unit)
log_msg = trim(log_msg) // chain
write(pdb_unit, '(A)') trim(log_msg)

do i_symm = 1, n_symm
    rotated_A(:) = matmul(matrix(:,:,i_symm), center(:,1)) 
    do i = 1, 3
        write(pdb_unit, 79) "REMARK 350   BIOMT", i, i_symm, &
                            matrix(i,:,i_symm), &
                            center(i,i_symm) - rotated_A(i)
    end do
end do
deallocate(ref_crd, crd_1, crd_n)

end subroutine write_symm_info
!===============================================================================
! Misc.
!===============================================================================
subroutine find_residue_in_first_symm_unit(res_no)
!-------------------------------------------------------------------------------
integer, intent(inout) :: res_no
integer :: i_symm, i_res

if (res_no <= symm_resrange(2,1)) return

do i_symm = 2, n_symm
    i_res = res_no - symm_resrange(2,i_symm-1)
    if (i_res <= symm_resrange(2,1)) then
        res_no = i_res
        exit
    end if
end do

end subroutine find_residue_in_first_symm_unit
!-------------------------------------------------------------------------------
END MODULE SYMMETRY
!-------------------------------------------------------------------------------
