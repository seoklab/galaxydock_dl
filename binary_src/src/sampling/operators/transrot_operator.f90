!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/sampling/operators/transrot_operator.f90
!
! Description:
!   This module constains subroutines related to translate/rotate protein,
!   hetmol, and ligand. Subroutine to rotate ligand torsion angle is also
!   included. 
!
!-------------------------------------------------------------------------------
MODULE TRANSROT_OPERATOR
!-------------------------------------------------------------------------------
use globals
use in_out_utils, only: find_atom_idx
use mathfunctions, only: v_norm
use geometry, only: quaternion, rotation, rotation_matrix

implicit none
save
public

CONTAINS
!-------------------------------------------------------------------------------
subroutine translate_protein(molecule, start_res, end_res, trans_v)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer, intent(in) :: start_res, end_res
real(dp), intent(in) :: trans_v(3)
integer :: i_res, i_atm

do i_res = start_res, end_res
    do i_atm = 1, molecule%residue(i_res)%n_atm
        molecule%residue(i_res)%R(:,i_atm) = molecule%residue(i_res)%R(:,i_atm) &
                                           + trans_v(:)
    end do
end do

end subroutine translate_protein
!-------------------------------------------------------------------------------
subroutine translate_hetmol(molecule, start_res, end_res, trans_v)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer, intent(in) :: start_res, end_res
real(dp), intent(in) :: trans_v(3)
integer :: i_res, i_atm

do i_res = start_res, end_res
    do i_atm = 1, molecule%hetmol(i_res)%n_atm
        molecule%hetmol(i_res)%R(:,i_atm) = molecule%hetmol(i_res)%R(:,i_atm) &
                                           + trans_v(:)
    end do
end do

end subroutine translate_hetmol
!-------------------------------------------------------------------------------
subroutine translate_ligand(molecule, start_res, end_res, trans_v)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer, intent(in) :: start_res, end_res
real(dp), intent(in) :: trans_v(3)
integer :: i_res, i_atm

do i_res = start_res, end_res
    do i_atm = 1, molecule%ligand(i_res)%n_atm
        molecule%ligand(i_res)%R(:,i_atm) = molecule%ligand(i_res)%R(:,i_atm) &
                                           + trans_v(:)
    end do
end do

end subroutine translate_ligand
!-------------------------------------------------------------------------------
subroutine rotate_protein(molecule, start_res, end_res, U, input_center)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer, intent(in) :: start_res, end_res
real(dp), intent(in) :: U(3,3)
real(dp), intent(in), optional :: input_center(3)
real(dp) :: tmp_R(3), center(3)
integer :: i_res, i_atm

if (present(input_center)) then
    center(:) = input_center(:)
else
    center(:) = 0.0d0
    do i_res = start_res, end_res
        i_atm = res_index(i_res)%Ca_id(2)
        center(:) = center(:) + molecule%residue(i_res)%R(:,i_atm)
    end do
    center(:) = center(:) / dble(end_res - start_res + 1)
endif

do i_res = start_res, end_res
    do i_atm = 1, molecule%residue(i_res)%n_atm
        tmp_R(:) = molecule%residue(i_res)%R(:,i_atm) - center(:)
        call rotation(tmp_R(:), U(:,:))
        molecule%residue(i_res)%R(:,i_atm) = tmp_R(:) + center(:)
    end do
end do

end subroutine rotate_protein
!-------------------------------------------------------------------------------
subroutine rotate_hetmol(molecule, start_res, end_res, U)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer, intent(in) :: start_res, end_res
real(dp), intent(in) :: U(3,3)
real(dp) :: tmp_R(3), center(3)
integer :: i_res, i_atm, tot_atm

center(:) = 0.0d0

tot_atm = 0
do i_res = start_res, end_res
    tot_atm = tot_atm + molecule%hetmol(i_res)%n_atm
    do i_atm = 1, molecule%hetmol(i_res)%n_atm
        center(:) = center(:) + molecule%hetmol(i_res)%R(:,i_atm)
    end do
end do
center(:) = center(:) / dble(tot_atm)

do i_res = start_res, end_res
    do i_atm = 1, molecule%hetmol(i_res)%n_atm
        tmp_R(:) = molecule%hetmol(i_res)%R(:,i_atm) - center(:)
        call rotation(tmp_R(:), U(:,:))
        molecule%hetmol(i_res)%R(:,i_atm) = tmp_R(:) + center(:)
    end do
end do

end subroutine rotate_hetmol
!-------------------------------------------------------------------------------
subroutine rotate_ligand(molecule, start_res, end_res, U)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
integer, intent(in) :: start_res, end_res
real(dp), intent(in) :: U(3,3)
real(dp) :: tmp_R(3), center(3)
integer :: i_res, i_atm, tot_atm

do i_res = start_res, end_res
    center(:) = molecule%ligand(i_res)%R(:,ref_lig(molecule%ligand(i_res)%lig_type)%cntr_atm)
    do i_atm = 1, molecule%ligand(i_res)%n_atm
        tmp_R(:) = molecule%ligand(i_res)%R(:,i_atm) - center(:)
        call rotation(tmp_R(:), U(:,:))
        molecule%ligand(i_res)%R(:,i_atm) = tmp_R(:) + center(:)
    end do
end do

end subroutine rotate_ligand
!-------------------------------------------------------------------------------
subroutine transrot_crd(crd, n_crd, trans_v, center, U)
!-------------------------------------------------------------------------------
! Rotate and translate crd about U and trans_v
! x'(3) = U(3,3)*(x(3)-center(3)) + trans_v(3) + center(3)
!-------------------------------------------------------------------------------
integer, intent(in) :: n_crd
real(dp), intent(inout) :: crd(3,n_crd)
real(dp), intent(in) :: trans_v(3)
real(dp), intent(in) :: center(3), U(3,3)
integer :: i_crd
real(dp) :: argv(3)
 
do i_crd = 1, n_crd
    argv(:) = crd(:,i_crd) - center(:)
    call rotation(argv(:), U(:,:))
    crd(:,i_crd) = argv(:) + center(:)
    crd(:,i_crd) = crd(:,i_crd) + trans_v(:)
end do

end subroutine transrot_crd
!-------------------------------------------------------------------------------
END MODULE TRANSROT_OPERATOR
!-------------------------------------------------------------------------------
