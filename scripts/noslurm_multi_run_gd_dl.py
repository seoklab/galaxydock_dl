import subprocess as sp
import sys
from concurrent.futures import ProcessPoolExecutor

from gd_dl.lib_pdb_mol2 import MOL2
from gd_dl.path_setting import (CASF_2016_CORE_SET,
                                MAIN_OUT_DIR,
                                CASF_DIR,
                                GD_DL_BIN_PATH,
                                OBRMS_PATH)

from run_gd_dl import preprocess_for_docking

BOX_SIZE = 22.5

STEP = 10

RMSD_CUTOFF = 2.0
TOPK = 1

TEST_S = list(filter(None,CASF_2016_CORE_SET.read_text().splitlines()))

NUM_TEST_SET = len(TEST_S)

OUT_DIR = MAIN_OUT_DIR/'test'/'gd_dl'
if not OUT_DIR.exists():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

def get_top_n_success_rate(rmsd_output_path, topk) -> int:
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

def run_gd_dl(out_dir) -> None:
    sp.run([str(GD_DL_BIN_PATH), './gd_dl.in'], check=True, cwd=str(out_dir))

def single_run(target) -> None:
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
            run_gd_dl(dest_dir)

def multi_run_gd_dl(mode: str,
                    n_cpu: int = None) -> None:
    if mode == 'result':
        total_success_rate = 0

    if mode == 'run' and n_cpu is not None:
        with ProcessPoolExecutor(n_cpu) as executor:
            executor.map(single_run, TEST_S)
    else:
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
                    run_gd_dl(dest_dir)
                        
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

    if mode == 'result':
        total_success_rate = total_success_rate/success_check_count*100
        print(total_success_rate)
        print(success_check_count)
    
    return

def create_prep_list(mode) -> None:
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
    mode = sys.argv[1]

    if len(sys.argv) == 3:
        n_cpu = int(sys.argv[2])
    else:
        n_cpu = None

    if mode == 'prep':
        prep_list = create_prep_list(mode)
        for args in prep_list:
            preprocess_for_docking(args)
    else:
        multi_run_gd_dl(mode,
                        n_cpu)
