!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/energy_weight.f90
!
! Description:
!-------------------------------------------------------------------------------
MODULE ENERGY_WEIGHT
!-------------------------------------------------------------------------------
use globals
use logger, only: log_p
use ramachandran, only: rama_mode
use energy_vars

implicit none
save
private

public :: set_default_energy_param
public :: setup_typical_weight
public :: report_weight

CONTAINS
!-------------------------------------------------------------------------------
subroutine set_default_energy_param()
!-------------------------------------------------------------------------------
! Set default variables for modeling before reading input file
! This values will be used unless parameters are specified on input file
!-------------------------------------------------------------------------------

energy_print_level = 0

! Pairlist
Ron = 8.0d0
Roff = 10.0d0
LRoff = 14.5d0

report_pairwise = .false.
energy_on_full = .false.

!
use_modeling_E = .false.
use_ppdock_E = .false.
use_ligdock_E = .false.

! Molecular Mechanics
use_bond = .true. !There is no part for reading in 'bond_w' so this should be default
use_LJ = .true. !LJ_type should be allocated for solv/sa even when vdw is off
use_softsphere = .false.
use_vdw_scwrl4 = .false.
use_elec = .false.
rattle_type = 0
!
bond_w = 1.0d0
vdw_w = 1.0d0
elec_w = 0.0d0
!
LJ_rep_w = 1.0d0
LJ_att_w = 1.0d0
!
rdie = .false.
die_const_pro = 1.0d0
die_const_slv = 78.5d0
!
soften_short = .true. ! TODO: conflict!
vdw_type = 'LJ'
vdw_r_scale = 1.0d0
trclash_scale = 0.8d0
scclash_scale = 0.7d0
vdw_soft_k = 200.0d0 ! It's needed when softsphere is on
!
truncate_by_shift = .true.
skip_fixed = .true.
pair_type_cut = 3

! solvation
use_solv = .false.
use_SA = .false.
use_factsmem = .false.
use_SApp = .false.
!
solv_w = 0.0d0
SA_w = 0.0d0
!
r_probe = 1.4d0
solv_type = SOLV_FACTS
SA_type = SOLV_FACTS  ! TODO: conflict!
ASP_type = 0
update_reff = .true.

! dDFIRE/DFIRE
use_ddfire = .false.
ddfire_add_scale = 1.0d0
!
dfire_w = 0.0d0
!
if (top_type == 'coarse') then
    dfire_score_file = trim(data_dir)//'score_ddfire_cg3.list'
    dfire_atype_file = trim(data_dir)//'ddfire_atomtype_cg.list'
else if (top_type == 'cbeta') then
    dfire_score_file = trim(data_dir)//'score_ddfire_cb.list'
    dfire_atype_file = trim(data_dir)//'ddfire_atomtype_cb.list'
else
    dfire_score_file = trim(data_dir)//'score_ddfire.list'
    dfire_atype_file = trim(data_dir)//'ddfire_atomtype.list'
end if

! Knowledge_GB
use_kgb = .false.
kgb_w = 0.0d0
wkgb_w = 0.1d0
kgb_water_file = ''
use_water_kgb = .false.
use_water_kgb = .false.

! GOAP
use_goap = .false.
goap_w  = 0.0d0
use_goap_discrete = .false.

! rsr from Machine Learning
use_distogram = .false.
use_bb_torsion = .false.

! Machine Learning score
use_ml = .false.

! Restraints
use_rsr = .false.
use_meld = .false.
!
rsr_w = 1.0d0
!
rsr_pdb = ''
rsr_file = ''
segrsr_scale = 100.0d0
!
use_update_rsr = .false.
update_NO_rsr = .false.
update_lig_rsr = .false.
update_rsr_dcut = 10.0d0
update_rsr_sig = 1.4142d0
update_lig_rsr_sig = 1.4142d0
infile_update_rsr_pdblist = ''

! Hydrogen bond
use_hbond = .false.
hbond_w =  0.0d0
use_coHbond = .false.
use_SAHbond = .false.

! Ramachandran score
use_rama_score = .false.
rama_w = 0.0d0

! Rotamer score
use_rotamer_score = .false.
rot_w = 0.0d0

! Conservation Score
use_conserve = .false.
conserve_w = 0.0d0
blastfile_defined = .false.
blastfile = ''

! InterEV score
use_interev = .false.
interev_w = 0.0d0

! Autodock score
use_atdk3 = .false.
use_atdk4 = .false.

! Drugscore
use_drugscore = .false.

! Xscore
use_Xscore = .false.

! PLP_tor
use_PLP = .false.

! dSAS
use_dSAS = .false.

! ROTA score
use_ROTA = .false.

! Docking Grid
use_input_cntr = .false.
dock_grid_info%n_elem(:) = 61
dock_grid_info%grid_width = 0.375

! softening option for ligand docking
soften_dock_E = .false.

end subroutine set_default_energy_param
!-------------------------------------------------------------------------------
subroutine setup_typical_weight(weight_type)
!-------------------------------------------------------------------------------
! Set pre-defined energy parameters
!-------------------------------------------------------------------------------
character(len=len_fname) :: weight_type
! PSloop & PSterminal are removed.

! for PS2loop protocol
if (trim(weight_type) == 'PS2loop') then
    call log_p('- Energy parameter is set to PS2loop.',me=me,level=20)
    call set_PS2loop()

! for PS1 TBM version
else if (trim(weight_type) == 'PS1tbm') then
    call log_p('- Energy parameter is set to PS1tbm.',me=me,level=20)
    call set_PS1tbm()

! for FM (not used)
else if (trim(weight_type) == 'Perseus') then
    call log_p('- Energy parameter is set to Perseus.',me=me,level=20)
    call set_Perseus()

else if (trim(weight_type) == 'Cassiopeia') then
    call log_p('- Energy parameter is set to Cassiopeia.',me=me,level=20)
    call set_Cassiopeia()

else if (trim(weight_type) == 'Prot2016') then
    call log_p('- Energy parameter is set to Prot2016.',me=me,level=20)
    call set_Prot2016()

else if (trim(weight_type) == 'C_allH') then
    call log_p('- Energy parameter is set to C_allH.',me=me,level=20)
    call set_C_allH()

else if (trim(weight_type) == 'Cassio7TM') then
    call log_p('- Energy parameter is set for membrane proteins.',me=me,level=20)
    call set_Cassio7TM()

else if (trim(weight_type) == 'PPDock') then
    call log_p('- Energy parameter is set for Protein-protein docking.',me=me,level=20)
    call set_PPDOCK()

else if (trim(weight_type) == 'PepDock') then
    call log_p('- Energy parameter is set for Protein-peptide docking.',me=me,level=20)
    call set_PEPDOCK()

else if (trim(weight_type) == 'GalaxySite') then
    call log_p('- Energy parameter is set to GalaxySite.',me=me,level=20)
    call set_GalaxySite()

else if (trim(weight_type) == 'GalaxyDock') then
    call log_p('- Energy parameter is set to GalaxyDock.',me=me,level=20)
    call set_GalaxyDock()

else if (trim(weight_type) == 'GalaxyDock2') then
    call log_p('- Energy parameter is set to GalaxyDock2.',me=me,level=20)
    call set_GalaxyDock2()

else if (trim(weight_type) == 'AutoDock3') then
    call log_p('- Energy parameter is set to AutoDock3.',me=me,level=20)
    call set_AutoDock3()

else if (trim(weight_type) == 'AutoDock4') then
    call log_p('- Energy parameter is set to AutoDock4.',me=me,level=20)
    call set_AutoDock4()

else if (trim(weight_type) == 'DrugScore') then
    call log_p('- Energy parameter is set to DrugScore.',me=me,level=20)
    call set_DrugScore()

else
    call log_p('WARNING: Currently available {weight_type} is', me=me, level=20)
    call log_p(' - PS2loop', me=me, level=20)
    call log_p(' - PS1tbm', me=me, level=20)
    call log_p(' - Perseus', me=me, level=20)
    call log_p(' - Cassiopeia', me=me, level=20)
    call log_p(' - Prot2016', me=me, level=20)
    call log_p(' - C_allH', me=me, level=20)
    call log_p(' - Cassio7TM', me=me, level=20)
    call log_p(' - PPDock', me=me, level=20)
    call log_p(' - PepDock', me=me, level=20)
    call log_p(' - GalaxyDock', me=me, level=20)
    call log_p(' - GalaxyDock2', me=me, level=20)
    call log_p(' - AutoDock3', me=me, level=20)
    call log_p(' - AutoDock4', me=me, level=20)
    call log_p(' - DrugScore', me=me, level=20)
end if

end subroutine setup_typical_weight
!-------------------------------------------------------------------------------
subroutine set_PS2loop()
!-------------------------------------------------------------------------------
use_modeling_E = .true.

! molecular mechanics
use_bond = .true.
use_LJ = .true.
use_elec = .true.
!
bond_w = 1.0d0
vdw_w = 1.0d0
elec_w = 0.2d0   !weight changed from 0.2d0-0.16d0
LJ_att_w = 1.0d0    ! default
LJ_rep_w = 1.0d0    ! default
!
rdie = .false.
die_const_pro = 1.0d0
die_const_slv = 78.5d0
!
soften_short = .true.
vdw_type = 'LJ'
vdw_r_scale = 1.0d0 ! default: 0.9d0

! solvation
use_solv = .true.
use_SA = .true.
!
solv_w = 0.2d0   !weight changed from 0.2d0-0.16d0
SA_w = 0.02d0    !weight changed from 0.05d0
!
solv_type = SOLV_FACTS
SA_type = SOLV_FACTS
ASP_type = 0

! dDFIRE
use_ddfire = .true.
ddfire_add_scale = 1.0d0
dfire_w = 5.0d0  !weight changed from 12.0d0

! restraint
if (trim(rsr_file) /= '') use_rsr = .true.
rsr_w = 1.0d0

! hydrogen bond
use_hbond = .true.
hbond_w = 0.6d0  !weight changed from 4.0d0
use_coHbond = .true.
use_SAHbond = .true.

! ramachandran
use_rama_score = .true.
rama_w = 0.6d0   !weight changed from 1.2d0
rama_mode = 'correct'

! rotamer score
use_rotamer_score = .true.
rot_w = 1.0d0

end subroutine set_PS2loop
!-------------------------------------------------------------------------------
subroutine set_PS1tbm()
!-------------------------------------------------------------------------------
use_modeling_E = .true.

Ron = 7.0d0
Roff = 9.0d0
LRoff = 14.5d0

! molecular mechanics
use_bond = .true.
use_LJ = .true.
use_elec = .false.
!
bond_w = 1.0d0
vdw_w = 1.0d0
elec_w = 0.0d0
LJ_att_w = 1.0d0
LJ_rep_w = 1.0d0
!
rdie = .true.
truncate_by_shift = .false. ! for cases when elec_w is not 0.0
die_const_pro = 1.0d0
die_const_slv = 78.5d0
!
soften_short = .true.
vdw_type = 'LJ'
vdw_r_scale = 1.0d0 ! default: 0.9d0

! solvation
use_solv = .true.
solv_w = 0.5d0
solv_type = SOLV_EEF1

! dDFIRE
use_ddfire = .true.
ddfire_add_scale = 1.0d0
dfire_w = 20.0d0

! restraint 
if (trim(rsr_file) /= '') use_rsr = .true.
rsr_w = 1.0d0

! hydrogen bond
use_hbond = .true.
hbond_w = 5.0d0
use_coHbond = .true.
use_SAHbond = .true.

! Ramachandran
use_rama_score = .true.
rama_w = 1.0d0
rama_mode = 'correct'

! Rotamer score
use_rotamer_score = .true.
rot_w = 1.0d0

end subroutine set_PS1tbm
!-------------------------------------------------------------------------------
subroutine set_Perseus()
!-------------------------------------------------------------------------------
use_modeling_E = .true.

! molecular mechanics
use_bond = .true.
use_LJ = .false.
use_softsphere = .true.
use_elec = .false.
!
bond_w = 1.0d0
vdw_w = 1.0d0
elec_w = 0.0d0
!
vdw_soft_k = 0.0d0

! solvation
use_SA = .true.
SA_w = 2.0d0
SA_type = SOLV_HABER
ASP_type = 3

! dDFIRE 
use_ddfire = .true.
ddfire_add_scale = 0.0d0
dfire_w = 5.0d0

! restraints
if (trim(rsr_file) /= '') use_rsr = .true.
rsr_w = 1.0d0

! hydrogen bond
use_hbond = .true.
hbond_w = 15.0d0

! Ramachandran
use_rama_score = .true.
rama_w = 1.5d0
rama_mode = 'hybrid'

end subroutine set_Perseus
!-------------------------------------------------------------------------------
subroutine set_Cassiopeia()
!-------------------------------------------------------------------------------
use_modeling_E = .true.

! Molecular mechanics
use_bond = .true.
use_LJ = .true.
use_elec = .true.
!
bond_w = 1.0d0
vdw_w = 1.0d0
elec_w = 0.2d0
LJ_att_w = 1.0d0
LJ_rep_w = 1.0d0
!
rdie = .false.
die_const_pro = 1.0d0
die_const_slv = 78.5d0
!
soften_short = .true.
vdw_type = 'LJ'
vdw_r_scale = 1.0d0

! solvation
use_solv = .true.
use_SA = .true.
solv_w = 0.2d0
SA_w = 0.04d0
solv_type = SOLV_FACTS
SA_type = SOLV_FACTS
ASP_type = 0

! dDFIRE 
use_ddfire = .true.
ddfire_add_scale = 1.0d0
dfire_w = 5.0d0

! restraint
if (trim(rsr_file) /= '') use_rsr = .true.
rsr_w = 0.5d0

! hydrogen bond
use_hbond = .true.
hbond_w = 2.0d0
use_coHbond = .true.
use_SAHbond = .true.

! Ramachandran
use_rama_score = .true.
rama_w = 2.0d0
rama_mode = 'neighbor'

! rotamer score
use_rotamer_score = .true.
rot_w = 2.0d0

end subroutine set_Cassiopeia
!-------------------------------------------------------------------------------
subroutine set_Prot2016()
!-------------------------------------------------------------------------------
use_modeling_E = .true.

! Molecular mechanics
use_bond = .true.
use_LJ = .true.
use_elec = .true.
!
bond_w = 1.0d0
vdw_w = 1.0d0
elec_w = 0.2d0
LJ_att_w = 1.0d0
LJ_rep_w = 1.0d0
!
rdie = .false.
die_const_pro = 1.0d0
die_const_slv = 78.5d0
!
soften_short = .true.
vdw_type = 'LJ'
vdw_r_scale = 1.0d0

! solvation
use_solv = .true.
use_SA = .true.
solv_w = 0.2d0
SA_w = 0.04d0
solv_type = SOLV_FACTS
SA_type = SOLV_FACTS
ASP_type = 0

! kgb
use_kgb = .true.
kgb_w = 10.0d0

! restraint
if (trim(rsr_file) /= '') use_rsr = .true.
rsr_w = 0.5d0

! hydrogen bond
use_hbond = .true.
hbond_w = 2.0d0
use_coHbond = .true.
use_SAHbond = .true.

! Ramachandran
use_rama_score = .true.
rama_w = 0.6d0
rama_mode = 'neighbor'

! rotamer score
use_rotamer_score = .true.
rot_w = 2.0d0

end subroutine set_Prot2016
!-------------------------------------------------------------------------------
subroutine set_C_allH()
!-------------------------------------------------------------------------------
use_modeling_E = .true.

! molecular mechanics
use_bond = .true.
use_LJ = .true.
use_elec = .true.
!
rattle_type = 1 !bonds with hydrogen constrained
!
bond_w = 1.0d0
vdw_w = 1.0d0
elec_w = 0.15d0
LJ_att_w = 1.0d0
LJ_rep_w = 1.0d0
!
rdie = .false.
!
vdw_r_scale = 1.0d0
vdw_type = 'LJ'

! solvation
use_solv = .true.
use_SA = .true.
solv_w = 0.15d0
SA_w = 0.02d0
solv_type = SOLV_FACTS
SA_type = SOLV_FACTS
ASP_type = 0

! dDFIRE
use_ddfire = .true.
ddfire_add_scale = 1.0d0
dfire_w = 6.0d0

! restraints
if (trim(rsr_file) /= '') use_rsr = .true.
rsr_w = 0.5d0

! Hydrogen bond
use_hbond = .true.
hbond_w = 2.0d0
use_SAHbond = .true.
use_coHbond = .true.

! Ramachandran
use_rama_score = .true.
rama_w = 4.0d0
rama_mode = 'standard'

! rotamer score
use_rotamer_score = .true.
rot_w = 2.0d0

end subroutine set_C_allH
!-------------------------------------------------------------------------------
subroutine set_Cassio7TM()
!-------------------------------------------------------------------------------
use_modeling_E = .true.

! molecular mechanics
use_bond = .true.
use_LJ = .true.
use_elec = .true.
!
bond_w = 1.0d0
vdw_w = 1.0d0
elec_w = 0.4d0
LJ_att_w = 1.0d0
LJ_rep_w = 1.0d0
!
rdie = .false.
die_const_pro = 1.0d0
die_const_slv = 80.0d0
!
vdw_r_scale = 1.0d0
vdw_type = 'LJ'

! solvation
use_solv = .true.
use_SA = .true.
use_factsmem = .true.
solv_w = 0.4d0
SA_w = 0.03d0
solv_type = SOLV_FACTS
SA_type = SOLV_FACTS
ASP_type = 0

! dDFIRE
use_ddfire = .true.
ddfire_add_scale = 1.0d0
dfire_w = 5.0d0

! restraint
if (trim(rsr_file) /= '') use_rsr = .true.
rsr_w = 0.5d0

! hydrogen bond
use_hbond = .true.
hbond_w = 3.0d0
use_SAHbond = .true.
use_coHbond = .true.

! Ramachandran
use_rama_score = .true.
rama_w = 2.5d0
rama_mode = 'neighbor'

! Rotamer score
use_rotamer_score = .true.
rot_w = 2.0d0

end subroutine set_Cassio7TM
!-------------------------------------------------------------------------------
subroutine set_PPDOCK()
!-------------------------------------------------------------------------------
Ron = 8.0d0
Roff = 10.0d0
LRoff = 14.5d0

! molecular mechanics
use_bond = .true.
use_LJ = .true.
use_elec = .true.
!
bond_w = 1.0d0
vdw_w = 1.0d0
elec_w = 0.15d0
LJ_att_w = 1.0d0
LJ_rep_w = 1.0d0
truncate_by_shift = .false.
skip_fixed = .false.
pair_type_cut = 4
!
rdie = .true.
die_const_pro = 1.0d0
die_const_slv = 78.5d0
!
soften_short = .true.
vdw_type = 'LJ'
vdw_r_scale = 0.9d0 ! default: 0.9d0

! solvation
use_solv = .false.
solv_w = 0.0d0
! SASA
use_SA = .true.
use_SApp = .true.
SA_w = 4.5d0
SA_type = SOLV_HABER
ASP_type = 2
r_probe = 1.0d0

! dDFIRE
use_ddfire = .true.
ddfire_add_scale = 0.0d0
dfire_w = 8.0d0

! restraint 
if (trim(rsr_file) /= '') use_rsr = .true.
rsr_w = 1.0d0

! hydrogen bond
use_hbond = .true.
hbond_w = 6.0d0

! Rotamer score
use_rotamer_score = .true.
rot_w = 3.0d0

! Conservation score
use_conserve = .true.
conserve_w = 3.0d0

use_modeling_E = .true.
use_ppdock_E = .true.

end subroutine set_PPDOCK
!-------------------------------------------------------------------------------
subroutine set_PEPDOCK()
!-------------------------------------------------------------------------------
Ron = 8.0d0
Roff = 10.0d0
LRoff = 14.5d0

! molecular mechanics
use_bond = .false.
use_LJ = .true.
use_elec = .true.
!
bond_w = 1.0d0
vdw_w = 1.0d0
elec_w = 0.15d0
LJ_att_w = 1.0d0
LJ_rep_w = 1.0d0
truncate_by_shift = .false.
skip_fixed = .false.
pair_type_cut = 4
!
rdie = .true.
die_const_pro = 1.0d0
die_const_slv = 78.5d0
!
soften_short = .true.
vdw_type = 'LJ'
vdw_r_scale = 0.9d0 ! default: 0.9d0

! solvation
use_solv = .false.
solv_w = 0.0d0
! SASA
use_SA = .true.
use_SApp = .true.
SA_w = 4.5d0
SA_type = SOLV_HABER
ASP_type = 2
r_probe = 1.0d0

! dDFIRE
use_ddfire = .true.
ddfire_add_scale = 1.0d0
dfire_w = 8.0d0

! restraint 
if (trim(rsr_file) /= '') use_rsr = .true.
rsr_w = 1.0d0

! hydrogen bond
use_hbond = .true.
hbond_w = 6.0d0

! Rotamer score
use_rotamer_score = .true.
rot_w = 3.0d0

! Conservation score
use_conserve = .true.
conserve_w = 3.0d0

use_modeling_E = .true.
use_ppdock_E = .true.

end subroutine set_PEPDOCK
!-------------------------------------------------------------------------------
subroutine set_GalaxySite()
!-------------------------------------------------------------------------------
use_bond = .false. !There is no part for reading in 'bond_w' so this should be default
use_LJ = .false.   !LJ_type should be allocated for solv/sa even when vdw is off

use_atdk3 = .true.
atdk_w = 1.0d0
atdk_vdw_w = 1.0d0
atdk_hbond_w = 1.0d0
atdk_elec_w = 0.1146d0
atdk_solv_w = 1.0d0

atdk_int_v_w = 1.0d0
atdk_int_h_w = 1.0d0
atdk_int_e_w = 0.1146d0
atdk_int_s_w = 1.0d0
atdk_int_w = 1.0d0

infile_atdk_prm = trim(data_dir) // 'atdk_eng_par.dat'
infile_atdk3_sol_prm = trim(data_dir) // 'atdk_sol_par.dat'

use_PLP = .true.
PLP_w = 0.1d0

use_rsr = .true.
rsr_w = 0.01d0

soften_dock_E = .true.

end subroutine set_GalaxySite
!-------------------------------------------------------------------------------
subroutine set_GalaxyDock()
!-------------------------------------------------------------------------------
use_bond = .false. !There is no part for reading in 'bond_w' so this should be default
use_LJ = .false.   !LJ_type should be allocated for solv/sa even when vdw is off

use_atdk3 = .true.
atdk_w = 1.0d0
atdk_vdw_w = 1.0d0
atdk_hbond_w = 1.0d0
atdk_elec_w = 0.1146d0
atdk_solv_w = 1.0d0

atdk_int_v_w = 1.0d0
atdk_int_h_w = 1.0d0
atdk_int_e_w = 0.1146d0
atdk_int_s_w = 1.0d0
atdk_int_w = 1.0d0

infile_atdk_prm = trim(data_dir) // 'atdk_eng_par.dat'
infile_atdk3_sol_prm = trim(data_dir) // 'atdk_sol_par.dat'

use_PLP = .true.
PLP_w = 0.1d0

use_rota = .true.
ROTA_w = 0.05d0
prot_w = 0.06d0

end subroutine set_GalaxyDock
!-------------------------------------------------------------------------------
subroutine set_GalaxyDock2()
!-------------------------------------------------------------------------------
use_bond = .false. !There is no part for reading in 'bond_w' so this should be default
use_LJ = .false.   !LJ_type should be allocated for solv/sa even when vdw is off

! AutoDock4
use_atdk4 = .true.
atdk_w = 1.0d0
atdk_vdw_w = 1.0d0 * 0.1662d0      
atdk_hbond_w = 0.85d0 * 0.1209d0    
atdk_elec_w = 0.93d0 * 0.1406d0     
atdk_solv_w = 0.12d0 * 0.1322d0    

atdk_int_v_w = 0.8d0 * 0.1662d0    
atdk_int_h_w = 0.25d0 * 0.1209d0    
atdk_int_e_w = 1.35d0 * 0.1406d0    
atdk_int_s_w = 1.0d0 * 0.1322d0    
atdk_int_w = 1.0d0

infile_atdk_prm = trim(data_dir) // 'AD4_default_parameters.dat'

! DrugScore
use_drugscore = .true.
drugscore_w = 0.07d0              
infile_drugscore_prm = trim(data_dir) // 'DS_pair_potential.txt'

! X-score (HM)
use_Xscore = .true.
X_HM_w = -0.8d0
conf_w = 0.0d0

! PLP_tor
use_PLP = .true.
PLP_w = 0.01d0

end subroutine set_GalaxyDock2
!-------------------------------------------------------------------------------
subroutine set_AutoDock3()
!-------------------------------------------------------------------------------
use_bond = .false. !There is no part for reading in 'bond_w' so this should be default
use_LJ = .false.   !LJ_type should be allocated for solv/sa even when vdw is off

use_atdk3 = .true.
atdk_w = 1.0d0
atdk_vdw_w = 1.0d0
atdk_hbond_w = 1.0d0
atdk_elec_w = 0.1146d0
atdk_solv_w = 1.0d0

atdk_int_v_w = 1.0d0
atdk_int_h_w = 1.0d0
atdk_int_e_w = 0.1146d0
atdk_int_s_w = 1.0d0
atdk_int_w = 1.0d0

infile_atdk_prm = trim(data_dir) // 'atdk_eng_par.dat'
infile_atdk3_sol_prm = trim(data_dir) // 'atdk_sol_par.dat'

end subroutine set_AutoDock3
!-------------------------------------------------------------------------------
subroutine set_AutoDock4()
!-------------------------------------------------------------------------------
use_bond = .false. !There is no part for reading in 'bond_w' so this should be default
use_LJ = .false.   !LJ_type should be allocated for solv/sa even when vdw is off

use_atdk4 = .true.
atdk_w = 1.0d0
atdk_vdw_w = 0.1662d0
atdk_hbond_w = 0.1209d0
atdk_elec_w = 0.1406d0
atdk_solv_w = 0.1322d0

atdk_int_v_w = 0.1662d0
atdk_int_h_w = 0.1209d0
atdk_int_e_w = 0.1406d0
atdk_int_s_w = 0.1322d0
atdk_int_w = 1.0d0

infile_atdk_prm = trim(data_dir) // 'AD4_default_parameters.dat'

end subroutine set_AutoDock4
!-------------------------------------------------------------------------------
subroutine set_DrugScore()
!-------------------------------------------------------------------------------
use_bond = .false. !There is no part for reading in 'bond_w' so this should be default
use_LJ = .false.   !LJ_type should be allocated for solv/sa even when vdw is off

use_drugscore = .true.
drugscore_w = 1.0d0

infile_drugscore_prm = trim(data_dir) // 'DS_pair_potential.txt'

end subroutine set_DrugScore
!-------------------------------------------------------------------------------
subroutine report_weight(me, log_unit, prefix)
!-------------------------------------------------------------------------------
integer, intent(in), optional :: me
integer, intent(in), optional :: log_unit
character(len=len_log), intent(in), optional :: prefix
character(len=len_log) :: prefix_out
real(dp) :: modeling_w(n_E_component_modeling)
real(dp) :: ppdock_w  (n_E_component_ppdock)
real(dp) :: ligdock_w (n_E_component_ligdock)
integer :: i
if (use_modeling_E) then
    do i = 1,4
        modeling_w(i) = bond_w
    enddo
    modeling_w(5) = vdw_w
    modeling_w(6) = elec_w
    modeling_w(7) = solv_w
    modeling_w(8) = SA_w
    modeling_w(9) = rot_w
    modeling_w(10) = rama_w
    modeling_w(11) = hbond_w
    if (use_ddfire) then
        modeling_w(12) = dfire_w
        modeling_w(13) = dfire_w*ddfire_add_scale
    else if (use_kgb) then
        modeling_w(12) = kgb_w
        modeling_w(13) = kgb_w
    endif
    do i = 14,15
        modeling_w(i) = goap_w
    enddo
    do i = 16,20
        modeling_w(i) = rsr_w
    enddo
endif

if (use_ppdock_E) then
    ppdock_w(1) = conserve_w
    ppdock_w(2) = interev_w
    ppdock_w(3) = interev_w
    ppdock_w(4) = symm_rsr_w
endif

if (use_ligdock_E) then
    do i = 1,4
        ligdock_w(i) = atdk_w
    enddo
    do i = 5,8
        ligdock_w(i) = atdk_int_w
    enddo
    ligdock_w(9) = drugscore_w
    ligdock_w(10) = x_HM_w
    ligdock_w(11) = PLP_w
    ligdock_w(12) = ROTA_w
    ligdock_w(13) = prot_w
endif
if (present(me) .and. me /= king) return
!
if (present(prefix)) then
    prefix_out = prefix
else
    prefix_out = '-'
end if
!
if (present(log_unit)) then
    if (use_modeling_E) then
        write(log_unit, 220) trim(prefix_out), "E_WEIGHT.2.0", modeling_w(1:n_E_component_modeling)
    end if
    if (use_ppdock_E) then
        write(log_unit, 220) trim(prefix_out), "E_WEIGHT.2.1", ppdock_w(1:n_E_component_ppdock)
    end if
    if (use_ligdock_E) then
        write(log_unit, 220) trim(prefix_out), "E_WEIGHT.2.2", ligdock_w(1:n_E_component_ligdock)
    end if
else
    if (use_modeling_E) then
        write(*, 220) trim(prefix_out), "E_WEIGHT.2.0", modeling_w(1:n_E_component_modeling)
    end if
    if (use_ppdock_E) then
        write(*, 220) trim(prefix_out), "E_WEIGHT.2.1", ppdock_w(1:n_E_component_ppdock)
    end if
    if (use_ligdock_E) then
        write(*, 220) trim(prefix_out), "E_WEIGHT.2.2", ligdock_w(1:n_E_component_ligdock)
    end if
end if

220 format (A,1x,A15,20(1x,f10.3))

end subroutine report_weight
!-------------------------------------------------------------------------------
END MODULE ENERGY_WEIGHT
