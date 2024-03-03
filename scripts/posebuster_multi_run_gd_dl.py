import subprocess as sp
import sys
import time

from pathlib import Path

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

NODE_LIST = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 47, 48, 49]
SKIP_LIST = [i for i in range(1,50) if i not in NODE_LIST]
# "star": name of cpu slurm nodes
EXCLUDE_NODE_LIST = ['star%s'%(str(node_i).zfill(3)) for node_i in SKIP_LIST]
EXCLUDE_NODES = ','.join(EXCLUDE_NODE_LIST)

STEP = 10
BOX_SIZE = 25.0

RMSD_CUTOFF = 2.0
TOPK = 1

if DATA_SET.startswith('posebuster'):
    DATASET_CODES = list(filter(None,POSEBUSTER_SET.read_text().splitlines()))
    DATASET_DIR = POSEBUSTER_SET_DIR

PCPU = min(96, len(DATASET_CODES))

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

def get_top_n_success_rate(rmsd_output_path, topk):
    lines = list(filter(None,rmsd_output_path.read_text().splitlines()))
    assert len(lines) == 100
    
    for i, line in enumerate(lines):
        if i >= topk:
            return 0
        
        rmsd = float(line.split()[-1])
        if rmsd < RMSD_CUTOFF:
            return 1    

def run_gd_dl_slurm(out_dir):
    log_file = out_dir/'log'
    if log_file.exists():
        log_file.unlink()

    bash_file = Path('run_pb_slurm_.sh').resolve()
    with bash_file.open('w') as f:
        f.write('#!/bin/sh\n')
        f.write(f'#SBATCH -J gd_dl_{DATA_SET}\n')
        f.write('#SBATCH -p normal.q\n')
        f.write(f'#SBATCH --exclude={EXCLUDE_NODES}\n')
        f.write('#SBATCH -N 1\n')
        f.write('#SBATCH -n 1\n')
        f.write('#SBATCH -c 2\n')
        f.write('#SBATCH -o ./log\n')
        f.write('#SBATCH -e ./log\n')
        f.write('#SBATCH --nice=100000\n')
        f.write(f'{str(GD_DL_BIN_PATH)} ./gd_dl.in')

    sp.run(['sbatch',str(bash_file)], check=True, cwd=str(out_dir))
    
    bash_file.unlink()

def run_prep_slurm():
    for i in range(PCPU):
        bash_file = Path('run_prep_.sh').resolve()
        with bash_file.open('w') as f:
            f.write('#!/bin/sh\n')
            f.write(f'#SBATCH -J prep_{DATA_SET}\n')
            f.write('#SBATCH -p normal.q\n')
            f.write(f'#SBATCH --exclude={EXCLUDE_NODES}\n')
            f.write('#SBATCH -N 1\n')
            f.write('#SBATCH -n 1\n')
            f.write('#SBATCH -c 16\n')
            f.write(f'#SBATCH -o ./prep.log\n')
            f.write(f'#SBATCH -e ./prep.log\n')
            f.write('#SBATCH --nice=100000\n')
            f.write(f'python {__file__} {DATA_SET} prep {i}')

        sp.run(['sbatch',str(bash_file)], check=True)
        
        bash_file.unlink()

def hold_for_slurm():
    while True:
        result = sp.run(['squeue', '-n', f'gd_dl_{DATA_SET}'], stdout=sp.PIPE, text=True)
        output = result.stdout
        if output.count('PD') > 200:
            time.sleep(20)
        else:
            return

def multi_run_gd_dl(mode):
    if mode == 'result':
        total_success_rate = 0
        success_rate_count = 0

    count=0
    
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
                run_gd_dl_slurm(dest_dir)
                # In the case the number of slurm jobs should be less than a centain value
                
                if count > 1300:
                    hold_for_slurm()
                    
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
            count+=1

    if mode == 'result':
        total_success_rate = total_success_rate/success_rate_count*100
        print(total_success_rate)
        print(success_rate_count)
    
    return

def create_prep_list(mode):
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
    skip_list = SKIP_LIST
    mode = sys.argv[2]

    if mode == 'prep':
        if len(sys.argv) > 3:
            split_index = int(sys.argv[3])
        else:
            run_prep_slurm()
            exit()
        prep_list = create_prep_list(mode)
        split_prep_list = [args for i, args in enumerate(prep_list) if i%PCPU == split_index]
        for args in split_prep_list:
            preprocess_for_docking(args)
    else:
        multi_run_gd_dl(mode)
