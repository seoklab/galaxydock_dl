!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/apps/main_ligdock.f90
!
! description: do ligand docking using CSA
!
! references: LigDock-CSA / GalaxyDock / GalaxyDock2
!
!-------------------------------------------------------------------------------
PROGRAM MAIN_LIGDOCK
!-------------------------------------------------------------------------------
use globals
!
use logger, only: log_thick_divider, my_timer, finish_timer
use setup_molecule, only: initialize_molecule
use geometry, only: cartesian2internal, internal2cartesian
!
use energy, only: initialize_energy
use ligdock_energy, only: construct_dock_grid
!
use apps_utils, only: initialize_apps, finalize_apps
use ligdock_input, only: read_ligdock_input
use ligand_docking, only: ligdock_runner

implicit none

character(len=len_fname) :: infile_input
type(molecule_type) :: protein
type(ligand_type) :: ligand
logical :: include_init_conf

call initialize_apps(infile_input)
call read_ligdock_input(infile_input, include_init_conf)
!
call initialize_molecule(protein, ligand)
call initialize_energy(protein, ligand)
!
call cartesian2internal(1, protein%n_res, protein%residue(1:protein%n_res))
call internal2cartesian(1, protein%n_res, protein%residue(1:protein%n_res))
!
call construct_dock_grid(protein, ligand)
!
call log_thick_divider()
!
call ligdock_runner(protein, ligand, include_init_conf)
!
call finalize_apps()

!-------------------------------------------------------------------------------
END PROGRAM MAIN_LIGDOCK
!-------------------------------------------------------------------------------
