!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/energy_utils.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE ENERGY_UTILS
!-------------------------------------------------------------------------------
use globals
use energy_vars
!
use symmetry, only: build_appl_res_symm

implicit none
save
public

CONTAINS
!-------------------------------------------------------------------------------
subroutine R_to_protein(protein)
!-------------------------------------------------------------------------------
! Copy current Cartesian coordinate into molecule type
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
integer :: res_no, atm_no
integer :: i_atm

do i_atm = 1, tn%stdatm
    res_no = i_R(1,i_atm)
    atm_no = i_R(2,i_atm)
    protein%residue(res_no)%R(:,atm_no) = R(:,i_atm)
end do   

do i_atm = tn%stdatm + 1, tn%nonligatm
    res_no = i_R(1,i_atm)
    atm_no = i_R(2,i_atm)
    protein%hetmol(res_no-protein%n_res)%R(:,atm_no) = R(:,i_atm)
end do

do i_atm = 1, tn%ligatm
    res_no = i_L(1,i_atm)
    atm_no = i_L(2,i_atm)
    protein%ligand(res_no)%R(:,atm_no) = R(:,i_atm+tn%nonligatm)
end do

end subroutine R_to_protein
!-------------------------------------------------------------------------------
subroutine protein_to_R(protein)
!-------------------------------------------------------------------------------
! Copy current molecule type info into Cartesian coordinate
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer :: res_no, atm_no
integer :: i_atm

do i_atm = 1, tn%stdatm
    res_no = i_R(1,i_atm)
    atm_no = i_R(2,i_atm)
    R(:,i_atm) = protein%residue(res_no)%R(:,atm_no) 
end do

do i_atm = tn%stdatm + 1, tn%nonligatm
    res_no = i_R(1,i_atm)
    atm_no = i_R(2,i_atm)
    R(:,i_atm) = protein%hetmol(res_no-protein%n_res)%R(:,atm_no) 
end do

do i_atm = 1, tn%ligatm
    res_no = i_L(1,i_atm)
    atm_no = i_L(2,i_atm)
    R(:,i_atm+tn%nonligatm) = protein%ligand(res_no)%R(:,atm_no)
end do

end subroutine protein_to_R
!-------------------------------------------------------------------------------
subroutine update_R_for_sc(residue, res_no)
!-------------------------------------------------------------------------------
type(residue_type), intent(in) :: residue
integer, intent(in) :: res_no
integer :: i_atm, atm_no, i_ref_res

i_ref_res = residue%res_type

do atm_no = 1, ref_res(i_ref_res)%n_atm
    if (ref_res(i_ref_res)%is_sc_atom(atm_no)) then
        i_atm = ii_R(atm_no,res_no)
        R(:,i_atm) = residue%R(:,atm_no) 
    end if
end do

end subroutine update_R_for_sc
!-------------------------------------------------------------------------------
subroutine update_R_for_ligand(ligand_residue, lig_no)
!-------------------------------------------------------------------------------
type(ligand_residue_type), intent(in) :: ligand_residue
integer, intent(in) :: lig_no
integer :: i_atm, atm_no

do atm_no = 1, ligand_residue%n_atm
    i_atm = ii_L(atm_no, lig_no)
    R(:,i_atm) = ligand_residue%R(:,atm_no)
end do

end subroutine update_R_for_ligand
!-------------------------------------------------------------------------------
subroutine build_appl_res(fix_type, is_usc)
!-------------------------------------------------------------------------------
! Setup appl_res and appl_res_pair.
!
! e.g.
!
!   fix_type   | ULR(input) | USC(input) |  appl_res
!  ------------------------------------------------------------------
!     all      |     x      |     x      |  hetmol/ligand
!              |     o      |     o      |  ULR/USC + hetmol/ligand
!              |     x      |     o      |  USC + hetmol/ligand
!  ------------------------------------------------------------------
!     none     |     x      |     x      |  ALL 
!              |     o      |     o      |  ALL
!              |     x      |     o      |  ALL
!  ------------------------------------------------------------------
!    backbone  |     x      |     x      |  ALL
!              |     x      |     o      |  ALL
!
!-------------------------------------------------------------------------------
integer, intent(in) :: fix_type
logical, intent(in) :: is_usc(:)
integer :: i_res

if (energy_on_full) then
    appl_res(:) = .true.
    appl_respair(:,:) = .true.
else
    appl_res(:) = .false.
    appl_respair(:,:) = .false.
    if (fix_type == FIX_ALL) then
        do i_res = 1, tn%stdres
            if (is_usc(i_res)) then
                appl_res(i_res) = .true.
            end if
        end do
        ! non-standard residues are always free to move.
        do i_res = tn%stdres + 1, tn%residue
            appl_res(i_res) = .true.
        end do
        if (symmetric) call build_appl_res_symm(appl_res)
    else ! fix_type 'none' or 'backbone'
        appl_res(:) = .true.
    end if

    do i_res = 1, tn%residue
        if (appl_res(i_res)) then
            appl_respair(i_res,:) = .true.
            appl_respair(:,i_res) = .true.
        end if
    end do
end if

end subroutine build_appl_res
!-------------------------------------------------------------------------------
END MODULE ENERGY_UTILS
!-------------------------------------------------------------------------------
