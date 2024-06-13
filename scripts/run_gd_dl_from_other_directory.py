import subprocess as sp
import math
from pathlib import Path
import torch
import argparse
from typing import List, Tuple, Dict
from gd_dl.lib_pdb_mol2 import PDB
from gd_dl.utils import str2bool
from gd_dl.path_setting import (
    OBABEL_PATH, GD_DL_BIN_PATH, SYBYL_FILE, RESIDUE_FILE)

if SYBYL_FILE.exists() and RESIDUE_FILE.exists():
    print('Parameter files exist')
else:
    raise Exception('Wrong global directory location or parameter files do not exist')

sybyl_lines = SYBYL_FILE.read_text().splitlines()

sybyl_types = []
for line in sybyl_lines:
    if line.startswith('ATOM') or len(line) < 3:
        continue
    sybyl_types.append(line.split()[0])

SYBYL = [i.lower() for i in sybyl_types]

sybyl_to_lig = {'c.3':(0,0,0),'c.2':(0,0,1),'c.1':(0,0,2),'c.ar':(0,0,3),'c.cat':(0,0,1),\
    'n.3':(1,0,0),'n.2':(1,0,1),'n.1':(1,0,2),'n.ar':(1,0,3),'n.am':(1,0,4),'n.pl3':(1,0,4),\
        'n.4':(1,0,0),'o.3':(2,0,0),'o.2':(2,0,1),'o.co2':(2,0,4),'s.3':(2,1,0),'s.2':(2,1,1),\
            's.o':(2,1,5),'s.o2':(2,1,5),'p.3':(1,1,0),'f':(3,0,0),'cl':(3,1,0),'br':(3,2,0),'i':(3,3,0)}

N_NODE = 29
N_LIG_FET = 1+4+4+6
N_SYBYL = 14

BOND_IDX = {'1':0,'2':1,'3':2,'ar':3,'am':4}

RES_IDX = {'ALA':0,'ARG':1,'ASN':2,'ASP':3,'CYS':4,'GLN':5,'GLU':6,'GLY':7,'HIS':8,'ILE':9, 'LEU':10,\
        'LYS':11, 'MET':12, 'PHE':13, 'PRO':14, 'SER':15, 'THR':16, 'TRP':17, 'TYR':18, 'VAL':19}

PROT_SYBYL = dict()

with RESIDUE_FILE.open('r') as f:
    for line in f:
        if line.startswith('RESI'):
            res_type = line.split()[1]
        elif line.startswith('ATOM'):
            if res_type not in RES_IDX and res_type != 'HET':
                continue
            atm_type, sybyl_type, _, vdw, _, charge, hbond, hydrophobic, solv, ring = line.split()[1:11]
            try:
                PROT_SYBYL[(res_type, atm_type)] = SYBYL.index(sybyl_type.lower())
            except:
                continue

def extract_contacts_and_cofactors(protein_pdb, out_dir):
    metal_ion_set = {'LI','NA','K','CA','MG','AL','MN','FE','NI','CD','CO','CU','ZN','HG'}
    lines = protein_pdb.read_text().splitlines()

    contact_lines = []
    cofactor_dict = {}

    for line in lines:
        if line.startswith('ATOM'):
            contact_lines.append(line)
        elif line.startswith('HETATM'):
            res_name = line[17:20].strip()
            element_name = line[76:78].strip()
            if res_name in metal_ion_set and element_name.upper() in metal_ion_set:
                contact_lines.append(line)
                if res_name not in cofactor_dict:
                    cofactor_dict[res_name] = line
                else:
                    continue

    cofactor_file_dict = {}
    if len(cofactor_dict) > 0:
        for cofac in cofactor_dict:
            cofactor_file = out_dir/f'{cofac}.pdb'
            out_mol2 = cofactor_file.with_suffix('.mol2')
            if not out_mol2.exists():
                with open(cofactor_file, 'w') as f:
                    f.write(cofactor_dict[cofac])
                sp.check_call(f'{OBABEL_PATH} {str(cofactor_file)} -O {str(out_mol2)}',shell=True)

                mol2_lines = out_mol2.read_text().splitlines()

                with out_mol2.open('w') as f:
                    read_check = False
                    for line in mol2_lines:
                        if line.startswith('@<TRIPOS>ATOM'):
                            read_check = True
                            f.write(f'{line}\n')
                        elif line.startswith('@<TRIPOS>BOND'):
                            read_check = False
                            f.write(f'{line}\n')
                        elif read_check:
                            charge = line.split()[-1]
                            new_line = '0.3000'.join(line.rsplit(charge,1))
                            if '\n' not in new_line:
                                new_line += '\n'
                            f.write(new_line)
                        else:
                            f.write(f'{line}\n')

            cofactor_file_dict[cofac] = out_mol2

    return cofactor_file_dict, contact_lines

def process_input_status(contact_lines, out_dir):
    protein_pdb_processed = out_dir/'pdb_processed.pdb'

    with protein_pdb_processed.open('w') as f:
        f.write('\n'.join(contact_lines))

    pdb = PDB(str(protein_pdb_processed))

    missing_bb = False
    for residue in pdb[0].get_residues():
        if residue.isHetatm():
            if residue.resName() == 'HOH':
                continue
            if residue.resName() == 'WAT':
                continue
        else:
            if not residue.check_bb():
                missing_bb = True
    if missing_bb:
        wrt = pdb.write()
        with protein_pdb_processed.open('w') as f:
            f.writelines(wrt)

    return protein_pdb_processed

def create_gd_dl_in_file(program_dir: Path,
                       ligand_mol2: Path,
                       contact_lines: List[str],
                       cofactor_file_dict: Dict[str, Path],
                       center_coord: Tuple[float],
                       grid_n_elem: Tuple[int],
                       grid_width: float,
                       out_dir: Path,
                       random_seed: int,
                       prep_check: bool):
    receptor_pdb = process_input_status(contact_lines, out_dir)
    output_file_name = out_dir/'gd_dl.in'
    
    str_center_coord = map(str,center_coord)
    str_center_coord = '   '.join(str_center_coord)
    
    str_grid_n_elem = map(str,grid_n_elem)
    str_grid_n_elem = '   '.join(str_grid_n_elem)
    
    input_ligand_mol2 = ligand_mol2
    if prep_check:
        charged_ligand_mol2 = out_dir/'charged_ligand.mol2'
        preprocess_script_path = program_dir/'scripts'/'gd2_preprocess_ligand.py'
        sp.run(['chimera', '--nogui', str(preprocess_script_path), str(ligand_mol2), str(charged_ligand_mol2)], check=True)
        input_ligand_mol2 = charged_ligand_mol2
    
    with output_file_name.open('w') as f:
        f.write('%-21s %s\n'%('data_directory',str(program_dir/'src'/'gd_dl'/'data')))
        f.write('%-21s %s\n'%('top_type','polarh'))
        f.write('%-21s %s\n'%('print_level','30'))
        f.write('%-21s %s\n'%('energy_print_level','30'))
        f.write('%-21s %s\n'%('fix_type','all'))
        f.write('%-21s %s\n'%('weight_type','GalaxyDock2'))
        f.write('%-21s %s\n'%('first_bank','rand'))
        f.write('%-21s %s\n'%('ligdock_prefix','GalaxyDock'))
        f.write('%-21s %s\n'%('grid_box_cntr',str_center_coord))
        f.write('%-21s %s\n'%('grid_n_elem',str_grid_n_elem))
        f.write('%-21s %s\n'%('grid_width',str(grid_width)))
        f.write('%-21s %s\n'%('infile_pdb',str(receptor_pdb)))
        f.write('%-21s %s\n'%('infile_ligand',str(input_ligand_mol2)))
        if cofactor_file_dict != {}:
            for cofac in cofactor_file_dict:
                f.write('%-21s %s\n'%('infile_mol2_topo','%s %s'%(str(cofactor_file_dict[cofac]),cofac)))
        f.write('%-21s %s\n'%('infile_pre_ML',str(out_dir/'pre_DL.pt')))
        f.write('%-21s %s\n'%('random',str(random_seed)))        

    return receptor_pdb

def prepare_pdb(receptor_pdb):
    epsilon=1e-10
    coord = []
    atom_type = []

    with receptor_pdb.open('r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atm_type = line[12:16].strip()
                res_type = line[17:20]
                try:
                    tmp_idx = PROT_SYBYL[(res_type, atm_type)]
                    tmp_node = torch.zeros(N_NODE)
                    tmp_node[N_LIG_FET+tmp_idx] = 1.0
                    atom_type.append(tmp_node)
                except:
                    continue
                else:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coord.append(torch.tensor([x, y, z]))
            elif line.startswith('HETATM'):
                atm_type = line[76:78].strip().upper() #use elemnet symbol instead of atm name
                res_type = 'HET'
                try:
                    tmp_idx = PROT_SYBYL[(res_type, atm_type)]
                    tmp_node = torch.zeros(N_NODE)
                    tmp_node[N_LIG_FET+tmp_idx] = 1.0
                    atom_type.append(tmp_node)
                except:
                    continue
                else:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coord.append(torch.tensor([x, y, z]))

    atom_type_r = torch.stack(atom_type)

    coord_r = torch.stack(coord)
    ref_vec_r = []
    ref_vec_size_r = []
    coord_r_row = coord_r.unsqueeze(1)
    coord_r_col = coord_r.unsqueeze(0)
    coord_r_dm = torch.sqrt((coord_r_row - coord_r_col).pow(2).sum(-1))
    protein_indice = torch.arange(len(coord_r))
    for i, crd in enumerate(coord_r):
        if atom_type_r[i][-1] == 1.0 or atom_type_r[i][-2] == 1.0:
            mask = torch.lt(coord_r_dm[i],2.5)
        else:
            mask = torch.lt(coord_r_dm[i],1.9)
        mask[i] = False
        neighbor_indice = protein_indice[mask]
        if True not in mask:
            ref_vec = torch.zeros(3)
            ref_vec_size = 0.0
        else:
            sum_vec = torch.zeros(3)
            for idx in neighbor_indice:
                sum_vec += coord_r[idx]
            ref_vec = -(sum_vec - len(neighbor_indice)*crd)
            ref_vec_size = torch.sqrt(ref_vec.pow(2).sum(-1))
            ref_vec = ref_vec / (ref_vec_size + epsilon)
            if ref_vec_size < 1.0:
                ref_vec_size = 0.0
                ref_vec = torch.zeros(3)
            else:
                ref_vec_size = 1.0
        ref_vec_r.append(ref_vec)
        ref_vec_size_r.append(ref_vec_size)

    ref_vec_r = torch.stack(ref_vec_r)

    receptor_zip = [coord_r,atom_type_r,ref_vec_r,ref_vec_size_r]

    if len(coord_r) != len(atom_type_r):
        raise Exception('coord_r != atom_type_r len')

    return receptor_zip

def prepare_ligand(ligand_mol2,n_edge):
    atom_type = list()

    bond = {}
    neigh_dict = {}
    read_atom = False
    read_bond = False
    hydrogen_index = []
    atm_idx_list = []
    with ligand_mol2.open('r') as f:
        for line in f:
            if '@<TRIPOS>ATOM' in line:
                read_atom = True
            elif '@<TRIPOS>BOND' in line:
                read_bond = True
                read_atom = False
            elif 'SUBSTRUCTURE' in line:
                break
            elif read_atom and len(line) > 3:
                tmp_sybyl = line.split()[5].lower()

                atm_idx = int(line.split()[0])
                if tmp_sybyl[0] == 'h':
                    hydrogen_index.append(atm_idx)
                    continue
                try:
                    group,period,sybyl_features = sybyl_to_lig[tmp_sybyl]
                except:
                    print(str(ligand_mol2))
                    exit()
                tmp_node = torch.zeros(N_NODE)
                tmp_node[0] = 1.0
                tmp_node[1+group] = 1.0
                tmp_node[5+period] = 1.0
                tmp_node[9+sybyl_features] = 1.0
                atom_type.append(tmp_node)

                atm_idx_list.append(atm_idx)
            elif read_bond:
                splitted = line.split()
                if len(splitted) < 4:
                    read_bond = False
                    continue
                
                start = int(splitted[1])
                end = int(splitted[2])
                if start in hydrogen_index or end in hydrogen_index:
                    continue
                bond[(start, end)] = BOND_IDX[splitted[3]]
                bond[(end, start)] = BOND_IDX[splitted[3]]

    n_atom = len(atom_type)
    x = torch.stack(atom_type)

    cov_edge_index_list = []
    cov_edge_attr_list = []

    for i_atom in range(n_atom):
        for j_atom in range(i_atom):
            if (atm_idx_list[j_atom], atm_idx_list[i_atom]) in bond:
                tmp_edge_index = torch.tensor([j_atom, i_atom],dtype=torch.long)
                tmp_edge_attr = torch.zeros(n_edge)
                tmp_edge_attr[0] = 1.0
                tmp_edge_attr[1+bond[(atm_idx_list[j_atom],atm_idx_list[i_atom])]] = 1.0

                cov_edge_index_list.append(tmp_edge_index)
                cov_edge_attr_list.append(tmp_edge_attr)

                if i_atom not in neigh_dict:
                    neigh_dict[i_atom] = {j_atom}
                else:
                    neigh_dict[i_atom].add(j_atom)

                if j_atom not in neigh_dict:
                    neigh_dict[j_atom] = {i_atom}
                else:
                    neigh_dict[j_atom].add(i_atom)
    #print(neigh_dict)
    if len(cov_edge_index_list) == 0 or len(cov_edge_attr_list) == 0:
        raise Exception('no covalent bond')

    ligand_zip = [x, cov_edge_index_list, cov_edge_attr_list, n_atom, atm_idx_list, bond, neigh_dict]

    return ligand_zip

def create_static_files(receptor_pdb, ligand_mol2, out_dir):
    save_fn = out_dir/'pre_DL.pt'

    ligand_zip_1 = prepare_ligand(ligand_mol2,17)
    ligand_zip_2 = prepare_ligand(ligand_mol2,19)
    receptor_zip = prepare_pdb(receptor_pdb)

    total_zip = [receptor_zip,ligand_zip_1,ligand_zip_2]

    torch.save(total_zip,str(save_fn))

    return

def preprocess_for_docking(args):
    protein_pdb, ligand_mol2, home_dir, center_coord, out_dir, i, box_size, prep = args

    assert home_dir.exists()
    program_dir = home_dir
    
    grid_width = 0.375
    
    grid_n_elem = (int(math.ceil(box_size/grid_width)) + 1 for _ in range(3))
    
    if not out_dir.exists():
        out_dir.mkdir(parents=True)

    cofactor_file_dict, contact_lines = extract_contacts_and_cofactors(protein_pdb, out_dir)
    receptor_pdb = create_gd_dl_in_file(program_dir=program_dir,
                                      ligand_mol2=ligand_mol2,
                                      contact_lines=contact_lines,
                                      cofactor_file_dict=cofactor_file_dict,
                                      center_coord=center_coord,
                                      grid_n_elem=grid_n_elem,
                                      grid_width=grid_width,
                                      out_dir=out_dir,
                                      random_seed=i,
                                      prep_check=prep)
    create_static_files(receptor_pdb, ligand_mol2, out_dir)

    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GalaxyDock2-DL')
    parser.add_argument('-p', required=True, help='protein.pdb location')
    parser.add_argument('-l', required=True, help='ligand.mol2 location')
    parser.add_argument('-d', required=True, help='main folder location')
    parser.add_argument('-x', required=True, help='Center coordinate of x')
    parser.add_argument('-y', required=True, help='Center coordinate of y')
    parser.add_argument('-z', required=True, help='Center coordinate of z')
    
    parser.add_argument('--out_dir', help='Path of output directory', default=Path.cwd())
    parser.add_argument('--random_seed', type=int, default=0, help='Random seed integer value')
    parser.add_argument('--box_size', type=float, default=22.5, help='Length of docking box')
    parser.add_argument('--prep', type=str2bool, default=True, help='Preparation for partial charge calculation')

    parse_args = parser.parse_args()

    center_coord = (float(parse_args.x), float(parse_args.y), float(parse_args.z))

    out_dir = Path(parse_args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    args = (Path(parse_args.p).resolve(True),
            Path(parse_args.l).resolve(True),
            Path(parse_args.d).resolve(True),
            center_coord,
            out_dir,
            parse_args.random_seed,
            parse_args.box_size,
            parse_args.prep)
    preprocess_for_docking(args)
    sp.run([f'{str(GD_DL_BIN_PATH)}', './gd_dl.in'], check=True, cwd=out_dir)
