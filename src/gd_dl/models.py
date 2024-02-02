import torch
import torch.nn as nn

from torch_geometric.nn.glob.glob import global_add_pool

from typing import Tuple, Union

import torch.nn.functional as F
from torch import Tensor
from torch.nn import BatchNorm1d, Linear

from torch_geometric.nn.conv import MessagePassing
from torch_geometric.typing import Adj, OptTensor, PairTensor

class Weight_CGConv(MessagePassing):
    """
    The crystal graph convolutional operator from the
    `"Crystal Graph Convolutional Neural Networks for an
    Accurate and Interpretable Prediction of Material Properties"
    <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.145301>`_
    paper
    We modifed the crystal graph convolutional operator to incorporate manual edge weights
    """
    def __init__(self, channels: Union[int, Tuple[int, int]], dim: int = 0,
                 aggr: str = 'add', batch_norm: bool = False,
                 bias: bool = True, **kwargs):
        super().__init__(aggr=aggr, **kwargs)
        self.channels = channels
        self.dim = dim
        self.batch_norm = batch_norm

        if isinstance(channels, int):
            channels = (channels, channels)

        self.lin_f = Linear(sum(channels) + dim, channels[1], bias=bias)
        self.lin_s = Linear(sum(channels) + dim, channels[1], bias=bias)
        if batch_norm:
            self.bn = BatchNorm1d(channels[1])
        else:
            self.bn = None

        self.reset_parameters()

    def reset_parameters(self):
        self.lin_f.reset_parameters()
        self.lin_s.reset_parameters()
        if self.bn is not None:
            self.bn.reset_parameters()

    def forward(self, x: Union[Tensor, PairTensor], edge_index: Adj,
                edge_attr: OptTensor = None, edge_weight: OptTensor = None) -> Tensor:

        if isinstance(x, Tensor):
            x: PairTensor = (x, x)

        # propagate_type: (x: PairTensor, edge_attr: OptTensor)
        out = self.propagate(edge_index, x=x, edge_attr=edge_attr, edge_weight=edge_weight, size=None)
        out = out if self.bn is None else self.bn(out)
        out += x[1]
        return out

    def message(self, x_i, x_j, edge_attr: OptTensor, edge_weight: OptTensor) -> Tensor:
        if edge_attr is None:
            z = torch.cat([x_i, x_j], dim=-1)
        else:
            z = torch.cat([x_i, x_j, edge_attr], dim=-1)
            if edge_weight is not None:
                return edge_weight.unsqueeze(1) * self.lin_f(z).sigmoid() * F.softplus(self.lin_s(z))

        return self.lin_f(z).sigmoid() * F.softplus(self.lin_s(z))

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.channels}, dim={self.dim})'


class Cgc_block(nn.Module):
    def __init__(self,node_dim=32,edge_dim=16):
        super(Cgc_block, self).__init__()

        self.conv = Weight_CGConv(node_dim,dim=edge_dim)
        self.norm = nn.LayerNorm(node_dim)
        self.linear = nn.Linear(node_dim, node_dim)
        self.act = nn.ELU(inplace=True)

    def forward(self,args):
        x_in, e_idx, e_attr, e_weight = args
        x = self.conv(x_in, e_idx, e_attr, e_weight)
        x = self.norm(x)
        x = self.linear(x)
        x += x_in
        out = self.act(x)

        return (out, e_idx, e_attr, e_weight)

class Rerank_model(nn.Module):
    def __init__(self,node_dim_in=29,node_dim_hidden=64,edge_dim_in=17,edge_dim_hidden=32, ligand_only=True,readout='mean'):
        super(Rerank_model, self).__init__()
        self.embed_x = nn.Sequential(nn.Linear(node_dim_in, node_dim_hidden), nn.ELU(inplace=True))
        self.embed_e = nn.Sequential(nn.Linear(edge_dim_in, edge_dim_hidden), nn.ELU(inplace=True))
        self.readout = readout
        self.ligand_only = ligand_only
        self.mlp_dropout = 0.5

        # graph convolution layer

        self.covalent_block = Cgc_block(node_dim_hidden, edge_dim_hidden)

        self.noncovalent_block = Cgc_block(node_dim_hidden, edge_dim_hidden)

        self.mlp = nn.Sequential(nn.Linear(node_dim_hidden, node_dim_hidden//2),\
                nn.ELU(inplace=True), nn.Dropout(p=self.mlp_dropout), nn.Linear(node_dim_hidden//2, node_dim_hidden//4),\
                nn.ELU(inplace=True), nn.Dropout(p=self.mlp_dropout))

        self.linear_trans = nn.Linear(node_dim_hidden//4, 1, bias = False)

    def forward(self, G, n_atom):
        x0 = G.x
        e_idx = G.edge_index
        e_attr = G.edge_attr

        # Manual edge weights
        cov_mask =  torch.where(e_attr[:,0] > 0.5, 1.0, 0.0)
        non_mask = torch.where(e_attr[:,0] > 0.5, 0.0, 1.0)

        if self.ligand_only:
            mask = torch.where(x0[:,0] > 0.5, 1.0, 0.0)
            mask = mask.view(-1,1)

        e_attr = self.embed_e(e_attr)
        x0 = self.embed_x(x0)

        gcn_out, _, _, _= self.covalent_block((x0, e_idx, e_attr,cov_mask))

        gcn_out, _, _, _ = self.noncovalent_block((gcn_out, e_idx, e_attr,non_mask))

        gcn_out = self.mlp(gcn_out)

        if self.readout in ['mean', 'sum']:
            if self.ligand_only:
                gcn_out = torch.multiply(gcn_out, mask)

            if self.readout == 'mean':
                gcn_mean = global_add_pool(gcn_out, G.batch)
                energy_ligand = self.linear_trans(gcn_mean)
                energy_ligand = energy_ligand.view(-1)
                energy_ligand = energy_ligand / n_atom
            elif self.readout == 'sum':
                gcn_sum = global_add_pool(gcn_out, G.batch)
                energy_ligand = self.linear_trans(gcn_sum)
                energy_ligand = energy_ligand.view(-1)

        return energy_ligand
