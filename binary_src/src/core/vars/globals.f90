!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/vars/globals.f90
!
! Description:
!  This file defines global parameters and derived types.
!-------------------------------------------------------------------------------
MODULE GLOBALS
!-------------------------------------------------------------------------------

implicit none
save
public

!===============================================================================
! PARAMETERS
!===============================================================================
! MPI parameters
!-------------------------------------------------------------------------------
integer, parameter :: king = 0                   ! MPI king processor index
integer :: n_proc = 1                            ! # of processor
integer :: me                                    ! MPI processor index
integer, parameter :: mpi_tag_type  = 10         ! MPI tag for default
integer, parameter :: mpi_type_termination = 999 ! Send this for termination for any MPI programs

!-------------------------------------------------------------------------------
! I/O parameters
!-------------------------------------------------------------------------------
integer :: print_level_global                 ! print level adjustor
integer, parameter :: len_fname = 250         ! max length of a file name
integer, parameter :: len_log = 300 !TODO temporary for lcsa
character(len=len_fname) :: data_dir          ! GALAXY data directory
!
character(len=len_fname) :: infile_user_input ! input file name
character(len=len_fname) :: task_mode         ! type of task to be performed

!-------------------------------------------------------------------------------
! Constants
!-------------------------------------------------------------------------------
#ifdef SP
    integer, parameter :: dp = kind(1.0)         
#else
    integer, parameter :: dp = kind(1.0d0)      
#endif
real(dp), parameter :: small_real = 1.0d-10      ! To determine it is zero or not.
!
real(dp), parameter :: pi = 3.141592653589793238462643383279502884197d0
real(dp), parameter :: pi_half = 0.5d0*pi
real(dp), parameter :: two_pi = 2.0d0*pi
!
real(dp), parameter :: deg2rad = pi/180.0d0      ! Degree -> Radian
real(dp), parameter :: rad2deg = 180.0d0/pi      ! Radian -> Degree
!
real(dp), parameter :: E_conv = 418.4d0          ! Conversion from kcal to g*Ang**2/ps**2
real(dp), parameter :: gasconst = 1.9872065d-3   ! Ideal gas constant (R) in kcal/mole/K
real(dp), parameter :: boltzmann = 0.83143435d0  ! Boltzmann constant (kB) in g*Ang**2/ps**2/K/mole
real(dp), parameter :: electron2kcal = 18.2223d0 ! for consistency in charge with AMBER
real(dp), parameter :: RT = 0.5900991d0          ! kT at 300K

!-------------------------------------------------------------------------------
! Molecule-related parameters
!-------------------------------------------------------------------------------
character(len=len_fname) :: mol_type     ! Types for the molecule
!
integer, parameter :: max_chain = 30     ! max no of chain in a protein
!
integer, parameter :: max_pdbfile = 10000 ! max no of PDBs in [infile_pdblist]
integer, parameter :: max_res = 1000     ! max no of residues in a protein
integer, parameter :: max_het = 500      ! max no of hetero molecules in a protein
integer, parameter :: max_lig = 10       ! max no of ligand molecules in a protein
integer, parameter :: len_res_name = 4   ! max length of a residue name
!
integer, parameter :: max_atm = 30       ! max no of atoms in a residue
integer, parameter :: max_lig_atom = 200 ! max no of atoms in a ligand
integer, parameter :: len_atom_name = 4  ! max length of an atomic name
!
integer, parameter :: max_link = 100     ! max no of links (see link_type) 
!
integer, parameter :: max_lig_br = 71    ! max no of branches in ligand
!
logical :: symmetric                     ! Whether to allow symmetric modeling
logical :: use_remark350                 ! Whether REMARK 350 for symmetry
!
logical :: allow_broken_bond             ! allow missing residues

!-------------------------------------------------------------------------------
! Reference topology/energy-related parameters
!-------------------------------------------------------------------------------
integer :: num_ref_res                   ! # of reference residues
integer, parameter :: max_ref_res = 100  ! max no of reference residues
integer, parameter :: max_ang_E = 500    ! max no of ang in a res for energy calc
!
character(len=len_fname) :: top_type     ! Topology type to be used
character(len=10) :: force_field_type    ! Force field type (AMBER/CHARMM)
!
integer, parameter :: max_ang_prm = 1200  ! max no of ang in parameters
integer, parameter :: max_atm_q_typ = 1500! max no of atm typ with different charges
!integer, parameter :: max_atm_cls = 100   ! max no of extended atm class
!TEMP
integer, parameter :: max_atm_cls = 200   ! max no of extended atm class
integer, parameter :: max_bnd_prm = 500   ! max no of bnd in parameters

!-------------------------------------------------------------------------------
! Fixing atom parameters (related to fix type and ulr,ulr_sc assignment)
!-------------------------------------------------------------------------------
integer :: fix_type                           ! Atom fix type (FIX_ALL/NONE/BB)
integer, parameter :: FIX_ALL  = 1            ! Fix all atoms except ULR/USC
integer, parameter :: FIX_NONE = 2            ! Make free to move all atoms
integer, parameter :: FIX_BB   = 3            ! Fix all backbone atoms
logical :: are_defined                         
!
integer, parameter :: max_n_ulr = 100         ! Max no. of ulrs allowed
integer :: n_ulr                              ! No. of ulrs
integer :: n_usc                              ! No. of unreliable sidechains 
logical, allocatable :: is_usc(:)       ! Whether the sc is assigned to move

!-------------------------------------------------------------------------------
! Log parameters
!-------------------------------------------------------------------------------
character(len=len_log) :: log_msg        ! Temporary variable for log print out
!
character(len=60)  :: divider = &
   '------------------------------------------------------------'
character(len=100) :: longdivider = &
   '----------------------------------------------------------------------------------------------------'
character(len=60)  :: thick_divider = &
   '============================================================'

!===============================================================================
! DERIVED TYPES
!===============================================================================
type residue_type
!-------------------------------------------------------------------------------
! A derived type for dealing "Residues" like residue, hetmol in molecule_type.
! This type can have UPTO 30 atoms limited by max_atm.
!-------------------------------------------------------------------------------
! Residue information with reference residues
integer :: res_type                    ! index for ref_res
character(len=4) :: res_name           ! res_name from ref_res
!
! Residue information in the PDB file.
integer :: pdb_res_no                  ! res_no in input pdb
character(len=4) :: pdb_res_name       ! res_name in input_pdb 
character(len=1) :: ter_type           ! terminal type (N/C/ )
character(len=1) :: res_added          ! inserted residue index
character(len=4) :: code               ! PDB code (at 75:78)
character(len=4), allocatable :: pdb_atom_name(:) ! atom name in input_pdb 
character(len=1) :: chain              ! chain index where the residue belongs to 
integer :: i_chain                     ! chain index by integer
integer :: n_atm                       ! # of atoms in residue
!
! Atom geometry information 
real(dp), allocatable :: R(:,:)      ! (x,y,z)
real(dp), allocatable :: b(:,:)      ! bond vector
real(dp), allocatable :: quat(:,:)   ! quaternions
real(dp), allocatable :: b_len(:)      ! bond length
real(dp), allocatable :: b_ang(:)       ! bond angles
real(dp), allocatable :: t_ang(:)       ! torsion angles
!
! Atom misc information
integer, dimension(3) :: link_atm_no, link_res_no  ! 3 previous atoms linked to the current residue
logical, allocatable :: atm_read(:)     ! whether an atom exists in PDB file
logical, allocatable :: atm_placed(:)   ! whether an atom is placed or not.
logical, allocatable :: bnd_read(:)  ! whether two atoms exist in PDB file
logical, allocatable :: atom_fixed(:)   ! fixed atoms during the simulations
!integer, allocatable :: atm_in_bnd_E(:,:)! indices for atom pairs for bond energy
logical :: is_broken                    ! broken peptide-bond

end type residue_type
!-------------------------------------------------------------------------------
type ligand_residue_type
!-------------------------------------------------------------------------------
! A derived type for dealing "Ligands" like ligand in molecule_type.
! This type can have UPTO 100 atoms limited by max_lig_atm.
!-------------------------------------------------------------------------------
! Ligand information with reference ligands
integer :: lig_type                     ! index for ref_lig
integer :: res_type                     ! index for ref_res, rarely used
character(len=4) :: res_name            ! res_name from ref_res
!
! Ligand information in the PDB file.
integer :: pdb_res_no                   ! res_no in input pdb
character(len=4) :: pdb_res_name        ! res_name in input_pdb 
character(len=1) :: ter_type            ! terminal type
character(len=1) :: res_added           ! inserted residue index
character(len=4) :: code                ! PDB code (at 75:78)
character(len=4), allocatable :: pdb_atom_name(:) ! atom name in input_pdb 
character(len=3) :: sec_str             ! secondary structure
character(len=1) :: chain               ! chain index where the residue belongs to 
integer :: i_chain                      ! chain index by integer
integer :: n_atm                        ! # of atoms in residue
!
! Atom geometry information
real(dp), allocatable :: R(:,:)      ! (x,y,z)
real(dp), allocatable :: b(:,:)      ! bond vector
real(dp), allocatable :: quat(:,:)   ! quaternions
real(dp), allocatable :: b_len(:)      ! bond length
real(dp), allocatable :: b_ang(:)       ! bond angles
real(dp), allocatable :: t_ang(:)       ! torsion angles
!
! Atom misc information
integer, dimension(3) :: link_atm_no, link_res_no  ! 3 previous atoms linked to the current residue
logical, allocatable :: atm_read(:)     ! whether an atom exists in PDB file
logical, allocatable :: atm_placed(:)   ! whether an atom is placed or not.
logical, allocatable :: bnd_read(:)  ! whether two atoms exist in PDB file
logical, allocatable :: atom_fixed(:)   ! fixed atoms during the simulations
integer, allocatable :: atm_in_bnd_E(:,:)! indices for atom pairs for bond energy

end type ligand_residue_type
!-------------------------------------------------------------------------------
type link_type
!-------------------------------------------------------------------------------
! A derived type for dealing "Links" like disulfide bonds.
! (in fact link seems to dealing disulfide bonds only.)
! since each link has 2 residue, link_res_name, chain, link_res_no 
! has 2 elements.
!-------------------------------------------------------------------------------
integer :: link_type ! ref_res number of link_name
character(len=4) :: link_name ! type of link. Ex: disulfide bond -> 'DISU'
character(len=4) :: link_res_name(2) ! residue name
character(len=1) :: chain(2)  !default : ' '
integer :: link_res_no(2) !residue number
   
end type link_type
!-------------------------------------------------------------------------------
type molecule_type
!-------------------------------------------------------------------------------
! A derived type for dealing "Molecules" like proteins
!-------------------------------------------------------------------------------
! Molecule type = { protein }
character(len=len_fname) :: mol_type
!
! Molecule information
! Chains
integer :: n_chain                      ! # of chains
integer :: chain_res(2,max_chain)       ! indices for chain breaks
!
! Standard residues -> 20 amino acids, 5 nucleic acids
integer :: n_res                        ! # of residues
type(residue_type), allocatable :: residue(:)
!
! Hetero molecules -> non-standard residues having less than 30 atoms
!                     ie, water, metals, small ions
integer :: n_het                        ! # of hetero molecules
type(residue_type), allocatable :: hetmol(:)
!
! Ligand molecules -> non-standard residues having less than 150 atoms
!                     ie, binding ligands, cofactors, ...
integer :: n_lig                        ! # of ligand molecules
type(ligand_residue_type), allocatable :: ligand(:)
!
! Misc
integer :: n_link                       ! # of links
type(link_type), dimension(max_link) :: link
!
logical :: n_capped        ! N-term. capping by ACE
logical :: c_capped        ! C-term. capping by NME or NHE
logical :: cyclic          ! Whether it is cyclic peptide or not.

end type molecule_type
!-------------------------------------------------------------------------------
type simple_residue_type
!-------------------------------------------------------------------------------
! A derived type for dealing "Residues" like residue, hetmol in molecule_type
! with minial information. This type can have UPTO 30 atoms limited by max_atm.
!-------------------------------------------------------------------------------
integer :: n_atm
real(dp), allocatable :: R(:,:)      ! (x,y,z)
real(dp), allocatable :: b_len(:)      ! bond length
real(dp), allocatable :: b_ang(:)      ! bond angles
real(dp), allocatable :: t_ang(:)      ! dihedral angles

end type simple_residue_type
!-------------------------------------------------------------------------------
type simple_ligand_residue_type
!-------------------------------------------------------------------------------
! A derived type for dealing "Ligands" like ligand in molecule_type with minimal
! information. This type can have UPTO 100 atoms limited by max_lig_atm.
!-------------------------------------------------------------------------------
integer :: n_atm
real(dp), allocatable :: R(:,:)      ! (x,y,z)
real(dp), allocatable :: b_len(:)      ! bond length
real(dp), allocatable :: b_ang(:)      ! bond angles
real(dp), allocatable :: t_ang(:)      ! dihedral angles

end type simple_ligand_residue_type
!-------------------------------------------------------------------------------
type simple_molecule_type
!-------------------------------------------------------------------------------
! A derived type for dealing "Molecules" like proteins with minimal information
!-------------------------------------------------------------------------------
type(simple_residue_type), allocatable :: residue(:)
type(simple_residue_type), allocatable :: hetmol(:)
type(simple_ligand_residue_type), allocatable :: ligand(:)

end type simple_molecule_type
!-------------------------------------------------------------------------------
type ref_res_type
!-------------------------------------------------------------------------------
! A derived type for dealing "reference" residue geometry information
!-------------------------------------------------------------------------------
character(len=4) :: res_name                       ! 3 or 4 letter residue type
character(len=1) :: ter_type                       ! terminal type={N/C/None}
!
integer :: n_atm, n_gr                             ! no. atoms / groups in ref_res
! TODO: remove group, and find related components and remove them also
integer, dimension(max_atm) :: group               ! group number (for CHARMM)
character(len=4), dimension(-3:max_atm+1) :: atom_name  ! reference atom name
!
integer, dimension(2, -1:max_atm) :: atm_in_bnd    ! atom indices in bonds
real(dp), dimension(-1:max_atm)   :: b_len0        ! reference bond lengths
!
integer, dimension(3, 0:max_atm) :: atm_in_b_ang   ! atom indices in b_ang
integer, dimension(2, 0:max_atm) :: bnd_in_b_ang   ! bond indices in b_ang
real(dp), dimension(0:max_atm)   :: b_ang0         ! reference bond angles
!
integer, dimension(4, max_atm) :: atm_in_t_ang     ! atom indices in t_ang
integer, dimension(3, max_atm) :: bnd_in_t_ang     ! bond indices in t_ang
real(dp), dimension(max_atm)   :: t_ang0           ! reference torsion angles
!
integer :: i_atm_prev(3), i_atm_o                  ! atom indices for -N,-CA,-C,-O
integer :: i_bnd_prev(2)                           ! bond indices for -CA..N, -C..-CA
integer :: iphi, ipsi, iomg                        ! angle indices for phi,psi,omg
!integer :: iphi_curr, ipsi_prev, iomg_curr         ! angle indices for phi,psi,omg
integer :: i_rot_res                               ! residue index for rotamer
!
logical, dimension(max_atm) :: is_sc_atom          ! whether it is SC or not
integer :: n_chi, t_ang_for_chi(4)                 ! t_ang index corresponding to chi
logical, dimension(4,max_atm) :: dep_on_chi        ! whether it is changed by chi_s
!
integer :: n_ang_dep_on_t_ang(max_atm)             ! # of t_ang_s changed by t_ang
integer :: ang_dep_on_t_ang(4, max_atm)            ! index of t_ang_s changed by t_ang
real(dp) :: d_ang_by_t_ang(4, max_atm)             ! t_ang_diff_s
!
logical :: filled

end type ref_res_type 
!-------------------------------------------------------------------------------
type ref_res_eng_type
!-------------------------------------------------------------------------------
! A derived type for dealing "reference" residue energy information
!-------------------------------------------------------------------------------
integer :: n_bnd_E, n_ang_E     ! # of bonds, angles in a residue to be evaluated
!
integer, dimension(max_atm) :: qq_idx !charge type index for eng_para%charge 
!                                     !eng_para%charge( qq_idx(i_atm) ) 
!                                     ! = atomic charge for i_atm
!
integer, dimension(max_atm) :: atm_cls! atom type index for atom_cls
character(len=6), dimension(-2:max_atm+1) :: atom_cls! atom type for parameter assignment
!
integer, dimension(max_atm) :: bnd_prm !bond parameter index
!                                      ! eng_para%n_bnd(1:2, bnd_prm(i) ) 
!                                      ! = atm_in_bnd_E(1:2, i )
!
integer, dimension(2,-2:max_atm) :: atm_in_bnd_E ! atom index pair in bond
character(len=6), dimension(2,-2:max_atm) :: atom_in_bnd_E! atom name pair in bond
!                                                         ! Ex) (/'-CA ','-C '/)
!
integer, dimension(max_ang_E) :: ang_prm ! angle parameter index
integer, dimension(max_ang_E) :: ang_E_type !type of angle
!                                           !1 - dihedral
!                                           !2 - angle
!                                           !3 - improper torsion
!
integer, dimension(max_ang_E) :: n_atm_in_ang   ! # of atom in angle (3 or 4)
integer, dimension(4,max_ang_E) :: atm_in_ang_E ! atom index in angle 
character(len=6), dimension(4,max_ang_E) :: atom_in_ang_E ! atom name in angle

end type ref_res_eng_type
!-------------------------------------------------------------------------------
type eng_para_type
!-------------------------------------------------------------------------------
! A derived type for dealing energy parameters
!-------------------------------------------------------------------------------
integer :: n_atom_cls, n_q_typ ! # of atm class, charge type in a residue
integer :: n_bnd, n_ang        ! # of bonds, angles in a residue
!
character(len=6), dimension(max_atm_cls) :: atom_cls     ! atom type
character(len=4), dimension(max_atm_cls) :: coord_type   ! crd type
character(len=6), dimension(2,max_bnd_prm) :: atm_in_bnd ! atm in bnd
character(len=6), dimension(4,max_ang_prm) :: atm_in_ang ! atm in ang
!
real(dp), dimension(max_atm_cls) :: mass         ! atomic mass
real(dp), dimension(max_atm_q_typ) :: charge     ! atomic charge
logical, dimension(0:max_atm_q_typ) :: non0q     ! whether it is neutral or not
!                                                
real(dp), dimension(2,max_bnd_prm) :: bnd_para   ! bond parameter
!
! ang_type = {1: Proper dih.,  2: Bond ang, 3: Improper dih.}
integer, dimension(max_ang_prm) :: ang_type
integer, dimension(max_ang_prm) :: n_ang_fold    ! angular periodicity
real(dp), dimension(3,5,max_ang_prm) :: ang_para ! angle parameter
!
real(dp), dimension(4,max_atm_cls) :: LJ_para    ! Lennard-Jones parameter
real(dp) :: E14fac, V14fac                       ! prefactor for Coulomb / LJ potential
real(dp), dimension(5,max_atm_cls,max_atm_cls) :: aux_para ! auxiliary for fast evaluations
real(dp), dimension(max_atm_cls,max_atm_cls) :: vdwsum  ! sum of vdw radii

end type eng_para_type
!-------------------------------------------------------------------------------
type res_index_type
!-------------------------------------------------------------------------------
! A derived type for dealing various "Indices" for molecules
!  - Indices for Ca_id, Cb_id, bb_id, phi/psi/chi_atmid
!    1 -> indices for R
!    2 -> indices for residue_type
!-------------------------------------------------------------------------------
integer :: ref_res_no         ! ref_res index
!
integer :: n_atm              ! no. atoms 
!
integer :: Ca_id(2)           ! C-alpha index
integer :: Cb_id(2)           ! C-beta index
integer :: bb_id(5,2)         ! backbone index; N,CA,C,O,H
integer :: loc_no(-3:max_atm) ! Local atomic index
!
logical :: is_het             ! non-standard amino acid
logical :: ignore             ! ignore for energy evaluation
!
integer :: n_chi
integer :: phi_atmid(4,2)     ! atomic index constituting phi-angle 
integer :: psi_atmid(4,2)     ! atomic index constituting psi-angle 
integer :: chi_atmid(4,4,2)
   
end type res_index_type
!-------------------------------------------------------------------------------
type tot_num_type
!-------------------------------------------------------------------------------
! TODO: Comment
!-------------------------------------------------------------------------------
integer :: residue
integer :: stdres
integer :: hetmol
integer :: nonlig
integer :: ligand
integer :: recres

integer :: atom
integer :: stdatm
integer :: hetatm
integer :: nonligatm
integer :: ligatm
integer :: recatm

integer :: dof
!-------------------------------------------------------------------------------
end type tot_num_type
!-------------------------------------------------------------------------------
type ulr_assign_type
!-------------------------------------------------------------------------------
! ULR type for assigning ulr by reading input. Simple type of ULR containing 
! variables which are used only for setting up molecule and energy. 
! For modeling, use type 'ulr_model_type'
!-------------------------------------------------------------------------------
integer :: resrange(2)                    ! init / fin resno in protein
integer :: segrange(2)                    ! init / fin segment resno in protein
character(len=2) :: type                  ! whether ulr is loop type or segment type
character(len=2) :: subtype               ! whether ulr is loop or terminus
real(dp) :: importance                    ! relative sampling importance
!-------------------------------------------------------------------------------
end type ulr_assign_type
!-------------------------------------------------------------------------------
type ref_lig_type
!-------------------------------------------------------------------------------
character(len=3) :: lig_name   ! 3-code name of ligand
!
! Atom information in the mol2 file.
integer :: n_atm                                        ! No. of atoms
character(len=4), dimension(max_lig_atom) :: atom_name  
character(len=6), dimension(max_lig_atom) :: mol2_type    
real(dp), dimension(3,max_lig_atom) :: R                ! Atom coordinates
real(dp), dimension(max_lig_atom) :: charge             ! Atomic charge
integer :: cntr_atm                                     ! center of ref_lig
!
! For bond information from mol2 file
integer :: n_bnd
! Bond type index: read this information if it exists in mol2 file
! 1 - single bond / 2 - double bond / 3- triple bond
! am - amide / ar - aromatic
! du - dummy / un - unknown / nc - not connected
character(len=2), dimension(max_lig_atom*2) :: bnd_type
integer, dimension(3,max_lig_atom*2) :: bnd  ! idx, i_atm, j_atm
real(dp), dimension(max_lig_atom*2) :: b_len0
!
! For angle information derived from mol2
integer :: n_ang
integer, dimension(4,max_lig_atom*2) :: ang  ! idx, i_atm, j_atm, k_atm
real(dp), dimension(max_lig_atom*2) :: b_ang0
!
! For dihedral information derived from mol2
integer :: n_dih
integer, dimension(5,max_lig_atom*2) :: dih  ! idx, i_atm, j_atm, k_atm, l_atm
real(dp), dimension(max_lig_atom*2) :: d_ang0
!
! Atom class for FF-based energy evaluations
character(len=6) :: atom_type(max_lig_atom)
integer :: eng_para(max_lig_atom*2, 3)
integer :: atm_cls(max_lig_atom)             ! Index for eng_para
integer :: q_typ_no(max_lig_atom)

end type ref_lig_type
!-------------------------------------------------------------------------------
type ligand_ring_type
!-------------------------------------------------------------------------------
! A derived type for dealing ring structure in ligand.
! It is used to detect rotatable bond and assign aromatic atom type in Docking.
!-------------------------------------------------------------------------------
integer, dimension(30) :: member ! index of atoms which are component of ring
integer :: n_member              ! No. of atoms in ring
logical :: aromatic              ! aromatic character of ring

end type ligand_ring_type
!-------------------------------------------------------------------------------
type ligand_type 
!-------------------------------------------------------------------------------
! A derived type for dealing ligand information used in docking
!-------------------------------------------------------------------------------
integer :: lig_type      ! index to find ref_lig
integer :: lig_no        ! index to find in molecule type
integer :: n_atm         ! # of atoms in ligand molecule
!
! Info. about ring within ligand (eg. ring structure like benzene, cyclohexane, etc)
integer :: n_ring                               ! No. of ring structures
type(ligand_ring_type), dimension(10) :: rings  ! Info. about rings
!
! Info. about branches which are move togather when rotate bond.
integer :: n_br
integer, dimension(max_lig_atom, max_lig_br) :: atm_in_br ! Atom index 
!                                                         ! in same branch
integer, dimension(2,max_lig_br) :: bridge    ! atoms consisting bond which
                                              ! connects different branches
integer, dimension(max_lig_br) :: n_atm_br    ! number of atom in the branch
integer :: n_core_br                          
integer, dimension(max_lig_br) :: core_bridge ! branch idx attached to root
integer, dimension(max_lig_br-1) :: i_rot_tor ! index for rotatable torsion bond
                                              ! tor_ang =
                                              !         t_ang(i_rot_tor(i_tor))
integer, dimension(max_lig_atom) :: piece     ! fragmentation info.
                                              ! atom i is in fragment piece(i)
!
! length between br_init atom and other atom in same branch  
real(dp), dimension(max_lig_atom,max_lig_br) :: vec_len

end type ligand_type
!-------------------------------------------------------------------------------
type(ref_res_type), dimension(max_ref_res) :: ref_res
type(ref_res_eng_type), dimension(max_ref_res) :: ref_res_eng
type(ref_lig_type), dimension(max_lig) :: ref_lig
!
type(eng_para_type) :: eng_para
type(res_index_type), allocatable :: res_index(:)
type(tot_num_type) :: tn
type(ulr_assign_type) :: ULR(max_n_ulr)
!
!-------------------------------------------------------------------------------
END MODULE GLOBALS
!-------------------------------------------------------------------------------
