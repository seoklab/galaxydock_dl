!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/sampling/operators/ligand_operator.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE LIGAND_OPERATOR
!-------------------------------------------------------------------------------

use globals
use mathfunctions, only: v_norm
use geometry, only: quaternion, rotation_matrix, rotation, expmap_quat, reparam_expmap
!
use transrot_operator, only: translate_ligand, rotate_ligand

public :: construct_ligand
public :: get_initial_ligand_gene

CONTAINS
!-------------------------------------------------------------------------------
subroutine construct_ligand(protein, ligand, gene)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(in) :: ligand
real(dp), intent(in) :: gene(:)
!
real(dp) :: p(3), q(4)
real(dp) :: U(3,3)
integer :: i_tor, i_atm

! initialize ligand coordinates
do i_atm = 1, protein%ligand(ligand%lig_no)%n_atm
    protein%ligand(ligand%lig_no)%R(:,i_atm) = &
            ref_lig(ligand%lig_type)%R(:,i_atm)
end do

! rotate rotatable torsions in ligand
do i_tor = 1, ligand%n_br - 1
    call rotate_ligand_tor(protein%ligand(ligand%lig_no), ligand, &
                           i_tor, gene(6+i_tor))
end do

! rotate ligand
!angle = gene(4)
!axis(1:3) = gene(5:7)
!call v_norm(axis)
!call quaternion(axis, angle, p)
p=gene(4:6)
call reparam_expmap(p)
call expmap_quat(p, q)
call rotation_matrix(3, q, U)
call rotate_ligand(protein, ligand%lig_no, ligand%lig_no, U)

! translate ligand by gene(1:3)
call translate_ligand(protein, ligand%lig_no, ligand%lig_no, gene(1:3))

end subroutine construct_ligand
!-------------------------------------------------------------------------------
subroutine get_initial_ligand_gene(ligand_residue, ligand, gene)
!-------------------------------------------------------------------------------
type(ligand_residue_type), intent(in) :: ligand_residue
type(ligand_type), intent(in) :: ligand
real(dp), intent(out) :: gene(:)
integer :: center_idx, i_tor, tor_idx

! gene for translation
center_idx = ref_lig(ligand%lig_type)%cntr_atm
gene(1:3) = ligand_residue%R(:,center_idx) - ref_lig(ligand%lig_type)%R(:,center_idx)

! gene for rotation
gene(4:6) = (/0.0d0, 0.0d0, 0.0d0/)

! gene for torsions
do i_tor = 1, ligand%n_br-1
    tor_idx = ligand%i_rot_tor(i_tor)
    gene(i_tor+6) = ligand_residue%t_ang(tor_idx)
end do

end subroutine get_initial_ligand_gene
!-------------------------------------------------------------------------------
subroutine rotate_ligand_tor(ligand_residue, ligand, i_tor, tor_ang)
!-------------------------------------------------------------------------------
type(ligand_residue_type), intent(inout) :: ligand_residue
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: i_tor
real(dp), intent(in) :: tor_ang
!
integer :: tor_idx, tor_init, tor_end
integer :: i_atm, atm_idx
real(dp) :: theta, q(4), U(3,3)
real(dp) :: R_init(3), axis(3), tmp_R(3)
real(dp) :: length

tor_idx = ligand%i_rot_tor(i_tor)
ligand_residue%t_ang(tor_idx) = tor_ang

theta = tor_ang - ref_lig(ligand%lig_type)%d_ang0(tor_idx)

tor_init = ligand%bridge(1, i_tor+1)
tor_end = ligand%bridge(2, i_tor+1)

R_init(1:3) = ligand_residue%R(1:3, tor_init)

axis(1:3) = ligand_residue%R(1:3, tor_end) - R_init(1:3)

call v_norm(axis(1:3))
call quaternion(axis(1:3), theta, q)
call rotation_matrix(3, q, U)

do i_atm = 1, ligand%n_atm_br(i_tor+1)
    atm_idx = ligand%atm_in_br(i_atm, i_tor+1)

    tmp_R(1:3) = ligand_residue%R(1:3, atm_idx) - R_init(1:3)
    call rotation(tmp_R, U)

    length = dot_product(tmp_R, tmp_R)
    
    tmp_R(1:3) = tmp_R(1:3)/length * ligand%vec_len(i_atm, i_tor+1)

    tmp_R(1:3) = tmp_R(1:3) + R_init(1:3)
    ligand_residue%R(1:3, atm_idx) = tmp_R(1:3)
end do

end subroutine rotate_ligand_tor
!-------------------------------------------------------------------------------
END MODULE LIGAND_OPERATOR
!-------------------------------------------------------------------------------
