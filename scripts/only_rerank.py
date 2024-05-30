from pathlib import Path
import numpy as np
import torch
import argparse
import gc

from torch_geometric.data import Data
from torch_geometric.data import Batch

from gd_dl.lib_pdb_mol2 import PDB
from gd_dl.path_setting import (
    SYBYL_FILE, RESIDUE_FILE)

from gd_dl.rerank_model import Rerank_model

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
READOUT = 'mean'
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

def extract_contacts(protein_pdb):
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

    return contact_lines

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
                    print(str(ligand_mol2))
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

def create_static_files(receptor_pdb, ligand_mol2):
    ligand_zip_1 = prepare_ligand(ligand_mol2,17)
    ligand_zip_2 = prepare_ligand(ligand_mol2,19)
    receptor_zip = prepare_pdb(receptor_pdb)

    total_zip = [receptor_zip,ligand_zip_1,ligand_zip_2]

    return total_zip

def preprocess_for_docking(args):
    protein_pdb, ligand_mol2, out_dir = args.p, args.l, args.out_dir
    
    if not out_dir.exists():
        out_dir.mkdir(parents=True)

    contact_lines = extract_contacts(protein_pdb)
    receptor_pdb = process_input_status(contact_lines, out_dir)
    all_zip = create_static_files(receptor_pdb, ligand_mol2)

    return all_zip

def get_ref_vec(i_atom, neigh_dict, coord_l):
    epsilon=1e-10
    crd_ref = coord_l[i_atom]
    ref_vec = torch.zeros(3)
    for i in neigh_dict[i_atom]:
        crd_neigh = coord_l[i]
        ref_vec += crd_neigh - crd_ref
    ref_vec = -ref_vec
    ref_vec_size = torch.sqrt(ref_vec.pow(2).sum(-1))
    ref_vec = ref_vec / (ref_vec_size + epsilon)
    if ref_vec_size < 1.0:
        ref_vec_size = 0.0
    else:
        ref_vec_size = 1.0
    
    return ref_vec, ref_vec_size

def prepare_input(all_zip, mol_fn, feature_dict):
    dist_cutoff = feature_dict['dist_cutoff']
    bucket = feature_dict['bucket']
    n_edge = feature_dict['n_edge']
    
    receptor_zip = all_zip[0]
    coord_r, atom_type_r,ref_vec_r,ref_vec_size_r = receptor_zip
    if dist_cutoff == 5.0:
        ligand_zip = all_zip[1]
    elif dist_cutoff == 6.0:
        ligand_zip = all_zip[2]
    else:
        raise Exception('wrong dist_cutoff')
    x, cov_edge_index_list, cov_edge_attr_list, n_atom, atm_idx_list, bond, neigh_dict = ligand_zip

    coord = []
    mol2_unit_list = []
    i_model = -1
    with mol_fn.open('r') as f:
        read = False
        for l in f:
            if l.startswith('@<TRIPOS>'):
                if 'MOLECULE' in l:
                    i_model += 1
                    mol2_unit_list.append(MOL2_unit())
                    tmp_mol2 = mol2_unit_list[i_model]
                    tmp_coord = []
                elif 'ATOM' in l:
                    read = True
                elif 'BOND' in l:
                    read = False
                    coord.append(tmp_coord)
            elif read:
                ll = l.split()
                tmp = ll[5]
                if tmp[0] == 'H':
                    tmp_mol2.read_line(l)
                    continue
                tmp_coord.append(ll[2:5])
            tmp_mol2.read_line(l)
    coord = np.array(coord,dtype=np.float32)
    coord_l = torch.from_numpy(coord)
    
    # generate graph
    # ligand internal graph    
    
    graph_s = list()

    protein_indice = torch.arange(len(coord_r))

    coord_r = coord_r.unsqueeze(0)

    n_model = len(coord_l)
    
    for i_model in range(n_model):
        #clash_energy = 0.0
        coord_l_i = coord_l[i_model]
        # distances
        # generate ref_vec
        ref_vec_dict = {}
        for i_atom in range(n_atom):
            ref_vec, ref_vec_size = get_ref_vec(i_atom, neigh_dict, coord_l_i)
            ref_vec_dict[i_atom] = (ref_vec, ref_vec_size)
        #
        # internal
        # unsqueese for broadcasting --> graph
        y1 = coord_l_i.unsqueeze(0)
        y2 = coord_l_i.unsqueeze(1)
        dm_ll = torch.sqrt((y1-y2).pow(2).sum(-1))
        #
        non_edge_index_list = []
        non_edge_attr_list = [] 

        for i_atom in range(n_atom):
            for j_atom in range(i_atom):
                tmp_dist = dm_ll[i_atom,j_atom]
                if tmp_dist <= dist_cutoff and (atm_idx_list[j_atom], atm_idx_list[i_atom]) not in bond:
                    tmp_edge_index = torch.tensor([j_atom, i_atom],dtype=torch.long)
                    tmp_edge_attr = torch.zeros(n_edge)
                    tmp_idx = torch.bucketize(tmp_dist,bucket)
                    tmp_edge_attr[6+tmp_idx] = 1.0

                    ref_vec_i, ref_vec_size_i = ref_vec_dict[i_atom]
                    ref_vec_j, ref_vec_size_j = ref_vec_dict[j_atom]
                    if ref_vec_size_i == 1.0 and ref_vec_size_j == 1.0:
                        cos_theta = torch.dot(ref_vec_i,ref_vec_j)
                        tmp_edge_attr[-2] = cos_theta
                        tmp_edge_attr[-1] = 1.0
                    
                    non_edge_index_list.append(tmp_edge_index)
                    non_edge_attr_list.append(tmp_edge_attr)

        # rec-lig
        dm_pl = (y2 - coord_r).pow(2).sum(-1)
        
        # receptor residues on interface
        prot_node_dict = {}
        prot_node_list = []

        pl_edge_index_list = []
        pl_edge_attr_list = []

        for i_atom in range(n_atom):
            per_dist = dm_pl[i_atom]
            mask = torch.le(per_dist,dist_cutoff**2)
            tmp_dist_s = per_dist[mask]            
            per_prot_indice = protein_indice[mask]
            for idx, prot_idx in enumerate(per_prot_indice):
                tmp_dist = tmp_dist_s[idx]
                if prot_idx not in prot_node_dict:
                    prot_node_dict[prot_idx] = n_atom + len(prot_node_dict)                    
                    prot_node_list.append(atom_type_r[prot_idx])
                tmp_edge_index = torch.tensor([i_atom, prot_node_dict[prot_idx]],dtype=torch.long)                    
                tmp_edge_attr = torch.zeros(n_edge)
                tmp_dist = tmp_dist**0.5
                tmp_idx = torch.bucketize(tmp_dist,bucket)
                tmp_edge_attr[6+tmp_idx] = 1.0

                ref_vec_i, ref_vec_size_i = ref_vec_dict[i_atom]
                ref_vec_ip, ref_vec_size_ip = ref_vec_r[prot_idx], ref_vec_size_r[prot_idx]
                if ref_vec_size_i == 1.0 and ref_vec_size_ip == 1.0:
                    cos_theta = torch.dot(ref_vec_i,ref_vec_ip)
                    tmp_edge_attr[-2] = cos_theta
                    tmp_edge_attr[-1] = 1.0
                
                pl_edge_index_list.append(tmp_edge_index)
                pl_edge_attr_list.append(tmp_edge_attr)              

        if len(prot_node_list) > 0:
            x2 = torch.stack(prot_node_list)
            tot_x = torch.cat((x,x2), dim=0)
        else:
            tot_x = x
        tot_edge_index = cov_edge_index_list + non_edge_index_list + pl_edge_index_list
        tot_edge_attr = cov_edge_attr_list + non_edge_attr_list + pl_edge_attr_list

        tot_edge_index = torch.stack(tot_edge_index)
        tot_edge_index = tot_edge_index.t().contiguous()
        tot_edge_attr = torch.stack(tot_edge_attr)
        
        row, col = tot_edge_index
        row, col = torch.cat([row, col], dim=0), torch.cat([col, row], dim=0)
        tot_edge_index = torch.stack([row, col], dim=0)
        tot_edge_attr = torch.cat([tot_edge_attr, tot_edge_attr], dim=0)
        #
        dat = Data(x=tot_x, edge_index=tot_edge_index, edge_attr=tot_edge_attr)
        graph_s.append(dat)
    
    return graph_s, n_atom, mol2_unit_list

class MOL2_unit():
    def __init__(self):
        self.line_list = []
        self.energy = None
    
    def read_line(self,line):
        self.line_list.append(line)
    def read_energy(self,energy):
        self.energy = energy

def calc_e(args):
    out_dir = args.out_dir
    mol_fn = args.out_mol2
    rerank_model = args.d/'src'/'gd_dl'/'data'/'rerank_model.pt'
    
    device = torch.device('cpu')

    model_1_feature_dict = {'dist_cutoff' : 5.0,
    'bucket' : torch.tensor([1.5,1.9,2.3,2.7,3.1,3.5,4.0,4.5,5.0]),
    'n_edge': 17,'zip':1}
    model_2_feature_dict = {'dist_cutoff' : 6.0,
    'bucket' : torch.tensor([1.5,1.9,2.3,2.7,3.1,3.5,4.0,4.5,5.0,5.5,6.0]),
    'n_edge': 19,'zip':2}

    fn_sampling_model = rerank_model.parent/'sampling_model.pt'
    fn_rerank_model = rerank_model

    fn_sampling_model_list = [str(fn_sampling_model),
                              str(fn_sampling_model)[:-3] + '_0.pt',
                              str(fn_sampling_model)[:-3] + '_1.pt']
    
    sampling_model_list = []
    
    for fn_model_i in fn_sampling_model_list:      
        tmp_model = Rerank_model(node_dim_hidden=64,
                                 edge_dim_hidden=32,
                                 readout=READOUT,
                                 ligand_only=True).to(device)
        checkpoint = torch.load(fn_model_i, map_location=device)
        tmp_model.load_state_dict(checkpoint['model_state_dict'])
        tmp_model.eval()
        
        sampling_model_list.append(tmp_model)
    
    rerank_model = Rerank_model(node_dim_hidden=64,
                                edge_dim_hidden=32,
                                edge_dim_in=model_2_feature_dict['n_edge'],
                                readout=READOUT,ligand_only=True).to(device)

    checkpoint = torch.load(str(fn_rerank_model), map_location=device)
    rerank_model.load_state_dict(checkpoint['model_state_dict'])
    rerank_model.eval()

    all_zip = preprocess_for_docking(args)
    
    energy = 0.0
    with torch.no_grad():
        graph_s_1, n_atom, mol2_unit_list = prepare_input(all_zip, mol_fn, model_1_feature_dict)
        graph_s_1 = Batch.from_data_list(graph_s_1)
        graph_s_1 = graph_s_1.to(device=device)
        for model_i in sampling_model_list:
            energy_i = model_i(graph_s_1,n_atom)
            energy += energy_i
        
        del graph_s_1
        gc.collect()

        graph_s_2, _, _ = prepare_input(all_zip, mol_fn, model_2_feature_dict)
        graph_s_2 = Batch.from_data_list(graph_s_2)
        graph_s_2 = graph_s_2.to(device=device)
        energy_rerank = rerank_model(graph_s_2,n_atom)
        energy += energy_rerank
        
        del graph_s_2
        gc.collect()
        
        energy /= len(fn_sampling_model_list) + 1

    out_fn_mol2 = out_dir/'reranked_fb.mol2'
    out_fn_energy = out_dir/'reranked_score.info'
    
    for i, energy_ in enumerate(energy):
        mol2_unit_list[i].read_energy(energy_)

    mol2_unit_list.sort(key=lambda unit : unit.energy)

    with out_fn_mol2.open('w') as f:
        for unit in mol2_unit_list:
            f.writelines(unit.line_list)
            #f.write('\n')
    
    with out_fn_energy.open('w') as f:
        for unit in mol2_unit_list:
            f.write('%10.3f\n'%(unit.energy))

    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GalaxyDock-DL')
    parser.add_argument('-p', required=True, type=Path, help='protein.pdb location')
    parser.add_argument('-l', required=True, type=Path, help='ligand.mol2 location')
    parser.add_argument('-d', required=True, type=Path, help='main folder location')
    parser.add_argument('--out_mol2', required=True, type=Path, help='output.mol2 location')
    parser.add_argument('--out_dir', type=Path, help='Path of output directory', default=Path.cwd())

    args = parser.parse_args()
    
    assert args.d.exists()

    calc_e(args)
    
