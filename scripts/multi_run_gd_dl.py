import subprocess as sp
import sys
import time

from pathlib import Path

from gd_dl.lib_pdb_mol2 import MOL2
from gd_dl.path_setting import (CASF_2016_CORE_SET,
                                MAIN_OUT_DIR,
                                CASF_DIR,
                                GD_DL_BIN_PATH,
                                OBRMS_PATH)

from run_gd_dl import preprocess_for_docking

BOX_SIZE = 22.5

PCPU = 96

#NODE_LIST = list(range(46,56))# + list(range(37,47))
NODE_LIST = list(range(1,19)) + list(range(45,50))

SKIP_LIST = [i for i in range(1,50) if i not in NODE_LIST]
# "star": name of cpu slurm nodes
EXCLUDE_NODE_LIST = ['star%s'%(str(node_i).zfill(3)) for node_i in SKIP_LIST]
EXCLUDE_NODES = ','.join(EXCLUDE_NODE_LIST)

STEP = 10

RMSD_CUTOFF = 2.0
TOPK = 1

TEST_S = list(filter(None,CASF_2016_CORE_SET.read_text().splitlines()))

NUM_TEST_SET = len(TEST_S)

OUT_DIR = MAIN_OUT_DIR/'test'/'gd_dl'
if not OUT_DIR.exists():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

def get_top_n_success_rate(rmsd_output_path, topk):
    lines = list(filter(None,rmsd_output_path.read_text().splitlines()))
    assert len(lines) == 100
    
    count = 0
    for line in lines:
        count += 1
        if count > topk:
            return 0
        
        rmsd = float(line.split()[-1])
        if rmsd < RMSD_CUTOFF:
            return 1    

def run_gd_dl_slurm(out_dir):
    log_file = out_dir/'log'
    if log_file.exists():
        log_file.unlink()

    bash_file = Path('run_slurm_.sh').resolve()
    with bash_file.open('w') as f:
        f.write('#!/bin/sh\n')
        f.write(f'#SBATCH -J gd_dl\n')
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
            f.write(f'#SBATCH -J prep\n')
            f.write('#SBATCH -p normal.q\n')
            f.write(f'#SBATCH --exclude={EXCLUDE_NODES}\n')
            f.write('#SBATCH -N 1\n')
            f.write('#SBATCH -n 1\n')
            f.write('#SBATCH -c 4\n')
            f.write(f'#SBATCH -o ./prep.log\n')
            f.write(f'#SBATCH -e ./prep.log\n')
            f.write('#SBATCH --nice=100000\n')
            f.write(f'python {__file__} prep {i}')

        sp.run(['sbatch',str(bash_file)], check=True)
        
        bash_file.unlink()

def hold_for_slurm():
    while True:
        result = sp.run(['squeue', '-n', f'gd_dl'], stdout=sp.PIPE, text=True)
        output = result.stdout
        if output.count('PD') > 200:
            time.sleep(20)
        else:
            return

def multi_run_gd_dl(mode):
    if mode == 'result':
        total_success_rate = 0

    count=0
    success_check_count = 0
    for target in TEST_S:
        target_dir = CASF_DIR/target

        out_dir_target = OUT_DIR/target
        if not out_dir_target.exists():
            out_dir_target.mkdir(exist_ok=True)

        for i in range(STEP):
            dest_dir = out_dir_target/f'step{i}'
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
                lig_fn = target_dir/f'{target}_ligand.mol2'
                if output_path.exists():
                    with rmsd_output_path.open('w') as f:
                        sp.run([str(OBRMS_PATH), '-firstonly', str(lig_fn), str(output_path)], check=True, stdout=f)
                else:
                    print(f'ERROR: {str(output_path)} do not exist')
            elif mode == 'result':
                rmsd_output_path = dest_dir/'rmsd.info'
                assert rmsd_output_path.exists()
                
                if rmsd_output_path.exists():
                    total_success_rate += get_top_n_success_rate(rmsd_output_path, TOPK)
                    success_check_count += 1
            else:
                raise Exception('Wrong mode name')
            count+=1

    if mode == 'result':
        total_success_rate = total_success_rate/success_check_count*100
        print(total_success_rate)
        print(success_check_count)
    
    return

def create_prep_list(mode):
    assert mode == 'prep'
    prep_list = []
    for target in TEST_S:
        target_dir = CASF_DIR/target
        
        out_dir_target = OUT_DIR/target
        
        crystal_ligand_mol2_file = target_dir/f'{target}_ligand.mol2'
        
        assert crystal_ligand_mol2_file.exists()
        
        crystal_ligand_mol2 = MOL2(crystal_ligand_mol2_file)
        crystal_ligand_mol2.read()
        assert len(crystal_ligand_mol2) == 1
        center_coord = crystal_ligand_mol2[0].get_coordinates_np_array().mean(axis=0)
        
        for i in range(STEP):
            dest_dir = out_dir_target/f'step{i}'

            output_path = dest_dir/'GalaxyDock_fb.mol2'
            if output_path.exists():
                continue

            arg_i = (target_dir/f'{target}_protein.pdb',
            target_dir/f'{target}_ligand.mol2',
            center_coord,
            dest_dir,
            i,
            BOX_SIZE,
            True)
            prep_list.append(arg_i)

    return prep_list

if __name__ == '__main__':
    skip_list = SKIP_LIST
    mode = sys.argv[1]

    if mode == 'prep':
        if len(sys.argv) > 2:
            split_index = int(sys.argv[2])
        else:
            run_prep_slurm()
            exit()
        prep_list = create_prep_list(mode)
        split_prep_list = [args for i, args in enumerate(prep_list) if i%PCPU == split_index]
        for args in split_prep_list:
            preprocess_for_docking(args)
    else:
        multi_run_gd_dl(mode)
