import os as _os
from pathlib import Path as _Path

import importlib_resources as _ir

HOME_DIR = _Path.cwd()

MAIN_OUT_DIR = HOME_DIR / 'test'
MAIN_DATA_DIR = HOME_DIR / 'total_data'

CASF_DIR = MAIN_DATA_DIR / 'CASF-2016' / 'coreset'

POSEBUSTER_DIR = MAIN_DATA_DIR / 'posebusters'

POSEBUSTER_SET = POSEBUSTER_DIR / 'posebusters_benchmark_set_ids.txt'

POSEBUSTER_SET_DIR = POSEBUSTER_DIR / 'posebusters_benchmark_set'

PREPROCESS_SCRIPT_PATH = HOME_DIR / 'scripts' / 'gd2_preprocess_ligand.py'

CORINA_MAIN_DIR = MAIN_DATA_DIR / 'corina_ligands'

CONDA_BIN = _Path(_os.environ['CONDA_PREFIX']) / 'bin'

OBRMS_PATH = CONDA_BIN / 'obrms'
OBABEL_PATH = CONDA_BIN / 'obabel'

GD_DL_BIN_PATH: _Path = _ir.files('gd_dl.bin') / 'ligdock'
GD_DL_DATA_PATH: _Path = _ir.files('gd_dl') / 'data'
SYBYL_FILE = GD_DL_DATA_PATH / 'sybyl_lcs.types'
RESIDUE_FILE = GD_DL_DATA_PATH / 'X_score_res_lcs.prm'
assert GD_DL_BIN_PATH.exists()
assert GD_DL_DATA_PATH.is_dir()

DATA_DIR_3 = MAIN_DATA_DIR / 'data_3'

DATA_DIR = DATA_DIR_3
TRAIN_FILE = DATA_DIR / 'train_3.set'
VALID_FILE = DATA_DIR / 'valid_3.set'
CASF_2016_CORE_SET = DATA_DIR / 'core.set'
PDB_DATA_DICT = MAIN_DATA_DIR / 'total_pdb_data_dict_3.pkl'
CENTER_COORD_FILE = MAIN_DATA_DIR / 'total_center_coord_dict_3.pkl'