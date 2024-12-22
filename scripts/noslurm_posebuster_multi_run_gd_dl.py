import subprocess as sp
import sys

from pathlib import Path
from typing import List
from concurrent.futures import ProcessPoolExecutor

from gd_dl.path_setting import (POSEBUSTER_SET,
                                POSEBUSTER_SET_DIR,
                                MAIN_OUT_DIR,
                                GD_DL_BIN_PATH,
                                OBRMS_PATH,
                                PREPROCESS_SCRIPT_PATH,
                                CORINA_MAIN_DIR)
from gd_dl.lib_pdb_mol2 import MOL2
from run_gd_dl import preprocess_for_docking

DATA_SET = sys.argv[1]

STEP = 10
BOX_SIZE = 25.0

RMSD_CUTOFF = 2.0
TOPK = 1

if DATA_SET.startswith('posebuster'):
    DATASET_CODES = list(filter(None,POSEBUSTER_SET.read_text().splitlines()))
    DATASET_DIR = POSEBUSTER_SET_DIR

USE_CORINA = False

if DATA_SET.endswith('_corina'):
    USE_CORINA = True

OUT_DIR = MAIN_OUT_DIR/DATA_SET/'gd_dl'
if not OUT_DIR.exists():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
def preprocess_ligand(in_file: Path,
                      out_file: Path) -> None:
    sp.run(['chimera', '--nogui', str(PREPROCESS_SCRIPT_PATH), str(in_file), str(out_file)])

    return

def get_top_n_success_rate(rmsd_output_path, topk) -> int:
    lines = list(filter(None,rmsd_output_path.read_text().splitlines()))
    assert len(lines) == 100
    
    for i, line in enumerate(lines):
        if i >= topk:
            return 0
        
        rmsd = float(line.split()[-1])
        if rmsd < RMSD_CUTOFF:
            return 1    

def run_gd_dl(out_dir) -> None:
    sp.run([str(GD_DL_BIN_PATH), './gd_dl.in'], check=True, cwd=str(out_dir))

def single_run(pl_code) -> None:
    out_dir_target = OUT_DIR/pl_code
    if not out_dir_target.exists():
        out_dir_target.mkdir(exist_ok=True)

    for i in range(STEP):
        dest_dir = out_dir_target/f'step{i}'
        #output_path= dest_dir/'log'
        output_path= dest_dir/'GalaxyDock_fb.mol2'
        # skip if output file exists
        if output_path.exists():
            continue
        run_gd_dl(dest_dir)

def multi_run_gd_dl(mode: str,
                    n_cpu: int = None) -> None:
    if mode == 'result':
        total_success_rate = 0
        success_rate_count = 0

    if mode == 'run' and n_cpu is not None:
        with ProcessPoolExecutor(n_cpu) as executor:
            executor.map(single_run, DATASET_CODES)
    else:
        for pl_code in DATASET_CODES:
            out_dir_target = OUT_DIR/pl_code
            
            if not out_dir_target.exists():
                out_dir_target.mkdir(exist_ok=True)

            for i in range(STEP):
                dest_dir = out_dir_target/f'step{i}'
                #output_path= dest_dir/'log'
                output_path= dest_dir/'GalaxyDock_fb.mol2'
                if mode == 'run':
                    # skip if output file exists
                    if output_path.exists():
                        continue
                    run_gd_dl(dest_dir)

                elif mode == 'rmsd':
                    rmsd_output_path = dest_dir/'rmsd.info'
                    lig_fn = DATASET_DIR/pl_code/f'{pl_code}_ligand.mol2'
                    if output_path.exists():
                        with rmsd_output_path.open('w') as f:
                            sp.run([str(OBRMS_PATH), '-firstonly', str(lig_fn), str(output_path)], check=True, stdout=f)
                    else:
                        print(f'ERROR: {str(output_path)} do not exist')
                elif mode == 'result':
                    rmsd_output_path = dest_dir/'rmsd.info'
                    if not rmsd_output_path.exists():
                        continue
                    
                    total_success_rate += get_top_n_success_rate(rmsd_output_path, TOPK)
                    success_rate_count += 1
                else:
                    raise Exception('Wrong mode name')

    if mode == 'result':
        total_success_rate = total_success_rate/success_rate_count*100
        print(total_success_rate)
        print(success_rate_count)
    
    return

def create_prep_list(mode) -> List:
    assert mode == 'prep'
    prep_list = []
    for pl_code in DATASET_CODES:
        out_dir_target = OUT_DIR/pl_code
        data_dir_target = DATASET_DIR/pl_code
        
        crystal_ligand_sdf_file = data_dir_target/f'{pl_code}_ligand.sdf'
        crystal_ligand_mol2_file = crystal_ligand_sdf_file.with_suffix('.mol2')
        
        if not crystal_ligand_mol2_file.exists():
            preprocess_ligand(in_file=crystal_ligand_sdf_file,
                              out_file=crystal_ligand_mol2_file)
        
        crystal_ligand_mol2 = MOL2(crystal_ligand_mol2_file)
        crystal_ligand_mol2.read()
        assert len(crystal_ligand_mol2) == 1
        center_coord = crystal_ligand_mol2[0].get_coordinates_np_array().mean(axis=0)

        for i in range(STEP):
            dest_dir = out_dir_target/f'step{i}'

            output_path = dest_dir/'GalaxyDock_fb.mol2'
            if output_path.exists():
                continue
            
            if USE_CORINA:
                corina_ligand_mol2_file = CORINA_MAIN_DIR/'posebusters'/pl_code/f'corina_{pl_code}_ligand.mol2'
                charged_corina_ligand_mol2_file = CORINA_MAIN_DIR/'posebusters'/pl_code/f'charged_corina_{pl_code}_ligand.mol2'
                if not charged_corina_ligand_mol2_file.exists():
                    preprocess_ligand(in_file=corina_ligand_mol2_file,
                                      out_file=charged_corina_ligand_mol2_file)
                arg_i = (data_dir_target/f'{pl_code}_protein.pdb',
                charged_corina_ligand_mol2_file,
                center_coord,
                dest_dir,
                i,
                BOX_SIZE,
                False)
            else:
                arg_i = (data_dir_target/f'{pl_code}_protein.pdb',
                crystal_ligand_mol2_file,
                center_coord,
                dest_dir,
                i,
                BOX_SIZE,
                False)
                
            prep_list.append(arg_i)

    return prep_list

if __name__ == '__main__':
    mode = sys.argv[2]

    if len(sys.argv) == 4:
        n_cpu = int(sys.argv[3])
    else:
        n_cpu = None

    if mode == 'prep':
        prep_list = create_prep_list(mode)
        for args in prep_list:
            preprocess_for_docking(args)
    else:
        multi_run_gd_dl(mode,
                        n_cpu)
