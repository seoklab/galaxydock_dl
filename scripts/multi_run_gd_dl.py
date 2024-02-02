import subprocess as sp
import sys
import time
import pickle as pkl
from pathlib import Path

from gd_dl.path_setting import (TRAIN_FILE,
                                VALID_FILE,
                                CASF_2016_CORE_SET,
                                CENTER_COORD_FILE,
                                MAIN_OUT_DIR,
                                DATA_DIR,
                                GD_DL_BIN_PATH,
                                OBRMS_PATH)

from run_gd_dl import preprocess_for_docking

BOX_SIZE = 22.5
DATA_SET = sys.argv[1]

PCPU = 96

NODE_LIST = list(range(46,56))# + list(range(37,47))
NOT_EXIST_NODE_LIST = list(range(40,46))
SKIP_LIST = [i for i in range(1,56) if i not in NODE_LIST+NOT_EXIST_NODE_LIST]
# "star": name of cpu slurm nodes
EXCLUDE_NODE_LIST = ['star%s'%(str(node_i).zfill(3)) for node_i in SKIP_LIST]
EXCLUDE_NODES = ','.join(EXCLUDE_NODE_LIST)

STEP = 10

RMSD_CUTOFF = 2.0
TOPK = 1

TRAIN_S = list(filter(None,TRAIN_FILE.read_text().splitlines()))
VALID_S = list(filter(None,VALID_FILE.read_text().splitlines()))
TEST_S = list(filter(None,CASF_2016_CORE_SET.read_text().splitlines()))

with CENTER_COORD_FILE.open('rb') as f:
    CENTER_COORD_DICT = pkl.load(f)

NUM_TRAIN_SET = len(TRAIN_S)
NUM_VALID_SET = len(VALID_S)
NUM_TEST_SET = len(TEST_S)

if DATA_SET == 'train':
    RUN_N = NUM_TRAIN_SET
elif DATA_SET == 'valid':
    RUN_N = NUM_VALID_SET
elif DATA_SET == 'test':
    RUN_N = NUM_TEST_SET

OUT_DIR = MAIN_OUT_DIR/DATA_SET/'gd_dl'
if not OUT_DIR.exists():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

def create_gd_dl_arguments(run_i):
    if DATA_SET == 'train':
        target = TRAIN_S[run_i]
    elif DATA_SET == 'valid':
        target = VALID_S[run_i]
    elif DATA_SET == 'test':
       target = TEST_S[run_i]

    target_dir = DATA_DIR/target

    return target, target_dir

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

    bash_file = Path('run2_slurm_.sh').resolve()
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
            f.write('#SBATCH -c 4\n')
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

def multi_run_gd_dl(program_mode, mode):
    if mode == 'result':
        total_success_rate = 0

    if program_mode in ['train','valid','test']:
        count=0
        success_check_count = 0
        for run_i in range(RUN_N):            
            target, target_dir = create_gd_dl_arguments(run_i)

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
                    
                    total_success_rate += get_top_n_success_rate(rmsd_output_path, TOPK)
                    success_check_count += 1
                else:
                    raise Exception('Wrong mode name')
                count+=1

    if mode == 'result':
        total_success_rate = total_success_rate/count*100
        print(total_success_rate)
        print(count)
    
    return

def create_prep_list(program_mode, mode):
    assert mode == 'prep'
    prep_list = []
    if program_mode in ['train','valid','test']:
        for run_i in range(RUN_N):
            target, target_dir = create_gd_dl_arguments(run_i)
            out_dir_target = OUT_DIR/target

            center_coord = CENTER_COORD_DICT[target]
            for i in range(STEP):
                dest_dir = out_dir_target/f'step{i}'

                output_path = dest_dir/'GalaxyDock_fb.mol2'
                if output_path.exists():
                    continue

                arg_i = (target_dir/f'{target}_contact.pdb',
                target_dir/f'{target}_ligand.mol2',
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
        prep_list = create_prep_list(DATA_SET, mode)
        split_prep_list = [args for i, args in enumerate(prep_list) if i%PCPU == split_index]
        for args in split_prep_list:
            preprocess_for_docking(args)
    else:
        multi_run_gd_dl(DATA_SET, mode)
