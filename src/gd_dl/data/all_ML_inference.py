import os
import time

import argparse
import numpy as np
import torch

from torch_geometric.data import Data
from torch_geometric.data import Batch

import gc

from gd_dl.rerank_model import Rerank_model

READOUT = 'mean'

device = torch.device('cpu')

torch.set_num_threads(1)
torch.set_num_interop_threads(1)

N_NODE = 29

model_1_feature_dict = {'dist_cutoff' : 5.0,
 'bucket' : torch.tensor([1.5,1.9,2.3,2.7,3.1,3.5,4.0,4.5,5.0]),
  'n_edge': 17,'zip':1}

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

def load_files(pre_ML_fn, mol_fn):
    all_zip = torch.load(pre_ML_fn)

    coord = []
    i_model = -1
    with open(mol_fn) as f:
        read = False
        for l in f:
            if l.startswith('@<TRIPOS>'):
                if 'MOLECULE' in l:
                    i_model += 1
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
                    continue
                tmp_coord.append(ll[2:5])
    coord = np.array(coord,dtype=np.float32)
    coord_l = torch.from_numpy(coord)

    return all_zip, coord_l

def calculate_energy(all_zip, coord_l, model, feature_dict):
    dist_cutoff = feature_dict['dist_cutoff']
    bucket = feature_dict['bucket']
    n_edge = feature_dict['n_edge']

    coord_r, atom_type_r,ref_vec_r,ref_vec_size_r = all_zip[0]
    x, cov_edge_index_list, cov_edge_attr_list, n_atom, atm_idx_list, bond, neigh_dict =  all_zip[feature_dict['zip']]

    energy_s = []

    # generate graph

    protein_indice = torch.arange(len(coord_r))

    coord_r = coord_r.unsqueeze(0)

    n_model = len(coord_l)

    for i_model in range(n_model):
        if i_model%50 == 0:
            graph_s = []

        coord_l_i = coord_l[i_model]

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

        if len(graph_s) == 50:
            graph_s = Batch.from_data_list(graph_s)
            graph_s = graph_s.to(device=device)

            tmp_energy = model(graph_s,n_atom)

            energy_s.append(tmp_energy)

            del graph_s
            gc.collect()

    energy_s = torch.cat(energy_s)

    return energy_s

def calc_e(args):
    pre_ML_fn = args.infile_pre_ML
    mol_fn = os.path.join(args.mol2_prefix+'.mol2')

    fn_model = args.load_model
    
    fn_model_list = [fn_model, fn_model[:-3] + '_0.pt', fn_model[:-3] + '_1.pt']
    
    model_list = []
    
    for fn_model_i in fn_model_list:      
        tmp_model = Rerank_model(node_dim_hidden=64,
                                 edge_dim_hidden=32,
                                 readout=READOUT,
                                 ligand_only=True).to(device)
        checkpoint = torch.load(fn_model_i, map_location=device)
        tmp_model.load_state_dict(checkpoint['model_state_dict'])
        tmp_model.eval()
        
        model_list.append(tmp_model)

    all_zip, coord_l = load_files(pre_ML_fn, mol_fn)

    with torch.no_grad():
        energy_s = 0.0
        for model_i in model_list:
            energy_s += calculate_energy(all_zip, coord_l, model_i, model_1_feature_dict)

    energy_s /= len(model_list)
    
    print('%s is working!'%mol_fn)

    out_fn = os.path.join(args.mol2_prefix+'.mol2.th1')
    with open(out_fn, 'w') as f:
        for energy_ in energy_s:
            f.write('%10.3f\n'%(energy_))

    return

def main(args):
    st = time.time()
    calc_e(args)
    et = time.time()
    print ("Time for inference:", round(et-st, 2), " (s)")
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile_pre_ML', type=str, required=True,
                                            help='torch_file_for_ML')
    parser.add_argument('--mol2_prefix', type=str, required=True,
                                            help='Prefix of a mol2 file')
    parser.add_argument('--load_model', type=str, required=True,
                                            help='saved_model')
    args = parser.parse_args()

    main(args)
