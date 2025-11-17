#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GAT Dataset Builder for ORIH Project
-----------------------------------
This module constructs PyTorch Geometric dataset objects from
variant and graph tables for training, validation, and testing.

Author: ERCSB
Date: 2025-02-11
"""

import os.path as osp
import pandas as pd
import torch
from torch_geometric.data import Data, InMemoryDataset


# ============================================================
# 1️⃣ Build gene graph edges from gene-level variant table
# ============================================================
def from_gene_ver3(data_df, graph_df):
    """
    Build gene–gene interaction edges (bidirectional)
    using the provided variant dataframe and HumanNet graph.
    """
    genes = data_df.iloc[:, 9:].columns.to_list()
    print(f"Detected {len(genes)} genes")

    mapping = {gene: i for i, gene in enumerate(genes)}

    # Filter graph edges to genes present in dataset
    graph_df_lm = graph_df[graph_df['gene_id1'].isin(genes)]
    graph_df_cut = graph_df_lm[graph_df_lm['gene_id2'].isin(genes)].reset_index(drop=True)

    # Make bidirectional edge list
    humannet_lm_a = graph_df_cut.copy()
    graph_df_cut.columns = ['gene_id2', 'gene_id1']
    humannet_lm_b = graph_df_cut[['gene_id1', 'gene_id2']]
    humannet_lm_all = pd.concat([humannet_lm_a, humannet_lm_b]).drop_duplicates().reset_index(drop=True)

    # Convert gene IDs to indices
    geneid1 = [mapping[gid] for gid in humannet_lm_all['gene_id1']]
    geneid2 = [mapping[gid] for gid in humannet_lm_all['gene_id2']]
    gene_edge = torch.tensor([geneid1, geneid2])

    print(f"Final number of edges: {len(gene_edge[0])}")
    return gene_edge


# ============================================================
# 2️⃣ Base Dataset Class
# ============================================================
class _BaseGACorihDataset(InMemoryDataset):
    """Base dataset class for GAT variant modeling."""

    def __init__(self, root, transform=None, pre_transform=None):
        super().__init__(root, transform, pre_transform)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        return [
            self.raw_variant_file,  # variant matrix
            "st12_gene_ns.csv",     # gene-gene interaction
            "Korean_ASD_WGS_v4_2186.ped",  # pedigree
            "ko_info_ac2.csv"       # clinical labels (AC_check)
        ]

    @property
    def processed_file_names(self):
        return [self.processed_filename]

    def _load_common(self):
        """Common loading and preprocessing logic shared by all splits."""
        variant_df = pd.read_csv(self.raw_paths[0], low_memory=False)
        graph_df = pd.read_csv(self.raw_paths[1], index_col=0, dtype="str")
        ped_df = pd.read_csv(self.raw_paths[2], sep=" ", header=None, low_memory=False)
        patient_df = pd.read_csv(self.raw_paths[3], index_col=0)

        ped_df.columns = ["fam_id", "s", "f_id", "m_id", "sex", "pheno"]
        ped_patient_df = pd.merge(ped_df, patient_df, how="left")

        # Fill missing and normalize AC_check (2 → 0)
        ped_patient_df["AC_check"] = ped_patient_df["AC_check"].fillna(0).replace(2, 0)
        ped_variant_df = pd.merge(ped_patient_df, variant_df)

        return ped_variant_df, graph_df

    def _process_core(self):
        """Shared core logic to create graph-based Data objects."""
        ped_variant_df, graph_df = self._load_common()
        gene_edge = from_gene_ver3(ped_variant_df, graph_df)

        data_list = []
        for idx in range(len(ped_variant_df)):
            s_data = ped_variant_df.s.iloc[idx]
            y_data = torch.tensor(ped_variant_df["AC_check"][idx]).view(1, 1).float()
            x_data = torch.FloatTensor(ped_variant_df.iloc[idx, 9:]).view(-1, 1)
            data_list.append(Data(x=x_data, y=y_data, edge_index=gene_edge, s=s_data))

        data, slices = self.collate(data_list)
        torch.save((data, slices), self.processed_paths[0])


# ============================================================
# 3️⃣ Dataset Variants: Train / Val / Test
# ============================================================
class GACorih_v3_bi_trainsp0_st12_nca_sc_deg_Dataset(_BaseGACorihDataset):
    raw_variant_file = "orih_v3_bi_trainsp0_sc_nca_ct_ent.csv"
    processed_filename = "GACorih_v3_bi_trainsp0_st12_nca_sc_deg.pt"

    def process(self):
        self._process_core()


class GACorih_v3_bi_valsp0_st12_nca_sc_deg_Dataset(_BaseGACorihDataset):
    raw_variant_file = "orih_v3_bi_valsp0_sc_nca_ct_ent.csv"
    processed_filename = "GACorih_v3_bi_valsp0_st12_nca_sc_deg.pt"

    def process(self):
        self._process_core()


class GACorih_v3_bi_testsp0_st12_nca_sc_deg_Dataset(_BaseGACorihDataset):
    raw_variant_file = "orih_v3_bi_testsp0_sc_nca_ct_ent.csv"
    processed_filename = "GACorih_v3_bi_testsp0_st12_nca_sc_deg.pt"

    def process(self):
        self._process_core()




        