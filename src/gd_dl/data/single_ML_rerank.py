#!/usr/bin/env python


import os

import time

import argparse
import torch

from gd_dl.rerank_model import Rerank_model
from gd_dl.data.ML_inference import load_files, calculate_energy

READOUT = 'mean'

device = torch.device('cpu')

torch.set_num_threads(1)
#torch.set_num_interop_threads(1)

REP_DIST = 1.5
N_NODE = 29

model_1_feature_dict = {'dist_cutoff' : 5.0,
 'bucket' : torch.tensor([1.5,1.9,2.3,2.7,3.1,3.5,4.0,4.5,5.0]),
  'n_edge': 17,'zip':1}
model_2_feature_dict = {'dist_cutoff' : 6.0,
 'bucket' : torch.tensor([1.5,1.9,2.3,2.7,3.1,3.5,4.0,4.5,5.0,5.5,6.0]),
  'n_edge': 19,'zip':2}

def calc_e(args):
    pre_ML_fn = args.infile_pre_ML
    mol_fn = os.path.join(args.mol2_prefix+'.mol2')
    rerank_model = args.load_model

    fn_model_1 = rerank_model.replace('rerank_model.pt','sampling_model.pt')
    fn_model_2 = rerank_model

    model_1 = Rerank_model(node_dim_hidden=64,
    edge_dim_hidden=32,
    edge_dim_in=model_1_feature_dict['n_edge'],
    readout=READOUT,ligand_only=True).to(device)

    checkpoint = torch.load(fn_model_1, map_location=device)
    model_1.load_state_dict(checkpoint['model_state_dict'])
    model_1.eval()

    model_2 = Rerank_model(node_dim_hidden=64,
    edge_dim_hidden=32,
    edge_dim_in=model_2_feature_dict['n_edge'],
    readout=READOUT,ligand_only=True).to(device)

    checkpoint = torch.load(fn_model_2, map_location=device)
    model_2.load_state_dict(checkpoint['model_state_dict'])
    model_2.eval()

    all_zip, coord_l = load_files(pre_ML_fn, mol_fn)

    with torch.no_grad():
        energy_s_1 = calculate_energy(all_zip, coord_l, model_1, model_1_feature_dict)
        energy_s_2 = calculate_energy(all_zip, coord_l, model_2, model_2_feature_dict)

        energy = (energy_s_1 + energy_s_2)/2

    out_fn = os.path.join(args.mol2_prefix+'.mol2.th1')
    with open(out_fn, 'w') as f:
        for energy_ in energy:
            f.write('%10.3f\n'%(energy_))

    print('%s is working!'%mol_fn)

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
