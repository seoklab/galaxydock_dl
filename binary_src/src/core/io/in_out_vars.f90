!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/io/in_out_vars.f90
!
! Description:
! 
!-------------------------------------------------------------------------------
MODULE IN_OUT_VARS
!-------------------------------------------------------------------------------
use globals

implicit none
public

!===============================================================================
! LOCAL VARIABLES
!===============================================================================
! Input Files
!-------------------------------------------------------------------------------
character(len=len_fname) :: infile_pdb        ! input pdbfile name
character(len=len_fname) :: infile_pre_ML     ! input for ML
character(len=len_fname) :: infile_pdblist    ! list file name of input pdbs
character(len=len_fname) :: infile_seq        ! input seqfile name
character(len=len_fname) :: infile_ligand(2)  ! input ligand file
                                              ! 1 => file name ; 2 => lig name
character(len=len_fname) :: outfile_prefix    ! prefix for output pdb files
character(len=len_fname) :: logfile           ! log file name
character(len=len_fname) :: fix_atom_file     ! user-assigned fixing atom file name 
logical :: single_output                      ! number of output PDB files
!
logical :: evaluate_TMscore                   ! Evaluate TM-score with native
character(len=len_fname) :: infile_native     ! Native PDB file for quality eval
character(len=len_fname) :: infile_qa         ! input QA file
character(len=len_fname) :: infile_ppDock_opr ! input ppDock rigid-body oprs
real(dp), allocatable :: init_qa(:)           ! QA-score
!-------------------------------------------------------------------------------
!integer :: n_gap                         ! # of gaps
!character(len=len_fname) :: gap_type     ! Info for sequence gaps
!
integer :: n_h_lines                     ! # line of pdb HEADER
integer, parameter :: max_h_lines = 5000 ! max no of PDB header lines
character(len=len_fname) :: pdb_h_lines(max_h_lines) ! pdb HEADER
!
integer :: n_mol2_top                    ! # of different ligands
logical :: read_het                      ! whether to read hetero
logical :: multiple_task                 ! Input pdb given as single file or as a list
logical :: multiple_models               ! Input PDB has multiple models
!
logical :: allow_prot_state              ! Whether to allow diff protonation states
logical :: check_disulfide_pdb           ! Check disulfide bond by distance criteria automatically
!
real(dp), parameter :: ss2_cut = 6.25d0  ! disulfide, 2.5**2 (SSBOND length criterion)
real(dp), parameter :: cn2_cut = 2.56d0  ! sequence cut, 1.6**2 (peptide bond length criterion)
!
!-------------------------------------------------------------------------------
! Reference topology/energy-related parameters
!-------------------------------------------------------------------------------
character(len=len_fname) :: infile_parameter        ! Parameter file
!
integer, parameter :: max_topo_file = 10 ! max no of input topology files
!
integer :: n_topo_file !number of topology files.
!                      !topology file: topology_XXX.in (or XXX_nt.in /
!                      !                                   XXX_ct.in)
!
character(len=len_fname) :: infile_topo(max_topo_file)  ! array of topology file 
!
!if the value is ' ', the topology file is for average residue.
!if the value is 'N', the topology file is for N-terminal residue.
!if the value is 'C', the topology file is for C-terminal residue.
character(len=1)  :: topo_file_ter_type(max_topo_file)
!
!if the value is 'protein', the topology file is for protein.
!if the value is 'hetmol', the topology file is for hetero molecules.
character(len=10) :: topo_file_mol_type(max_topo_file)  
!
character(len=len_fname) :: infile_mol2_topo(2,max_lig) ! mol2 files for lig topology
!
integer :: n_topo_eng_file ! number of topology_eng file
!
!if the value is 'protein', the topology_eng file is for protein.
!if the value is 'hetmol', the topology_eng file is for hetero molecules.
character(len=10) :: topo_eng_file_mol_type(max_topo_file)
!
character(len=len_fname) :: infile_topo_eng(max_topo_file)! array of
!                                                         ! topology_eng files
!
character(len=len_fname) :: infile_parameter_ligand ! Parameter file for ligands
!
character(len=len_fname) :: infile_rotamer_lib  ! Rotamer library file
character(len=len_fname) :: infile_chi_def      ! Chi-angle definition file
!-------------------------------------------------------------------------------
! python scripts and models for ML
!-------------------------------------------------------------------------------
character(len=2*len_fname) :: infile_ML(4)
!-------------------------------------------------------------------------------
END MODULE IN_OUT_VARS
!-------------------------------------------------------------------------------
