!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ligdock/ligdock_E_vars.f90
!
! Description:
!   varialbes/derived types used in ligand docking energy calculation.
!
!-------------------------------------------------------------------------------
MODULE LIGDOCK_E_VARS
!-------------------------------------------------------------------------------

use globals
use energy_vars

implicit none
save
public

! for mapping mol2 type to amino acids
integer, parameter :: max_mol2_res = 40
character(len=len_fname) :: infile_aa_to_mol2  ! name of data file

!-------------------------------------------------------------------------------
type residue_mol2_type
!-------------------------------------------------------------------------------
! Doping mol2 type to protein
!-------------------------------------------------------------------------------
character(len=4) :: res_name
character(len=4), dimension(max_atm) :: atom_name
character(len=6), dimension(max_atm) :: mol2_type
integer :: n_atm
!-------------------------------------------------------------------------------
end type residue_mol2_type
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
END MODULE LIGDOCK_E_VARS
!-------------------------------------------------------------------------------
