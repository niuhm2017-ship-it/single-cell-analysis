# Extracted from E:/nhm/scRNA/scRNA_PY1.ipynb
# Review paths and parameters before running.

# %% Cell 0
# Core scverse libraries
import scanpy as sc
import anndata as ad
import os
import glob
import numpy as np
import pandas as pd
import scipy.io
from collections import Counter
import warnings

# Data retrieval (pooch kept for any remote fetches if needed)
import pooch

# %% Cell 1
import os
import scanpy as sc

h5_root = "/mnt/e/nhm/scRNA/h5"
sample_dirs = sorted([d for d in os.listdir(h5_root) if os.path.isdir(os.path.join(h5_root, d))])

adatas = []
for s in sample_dirs:
    candidate = os.path.join(h5_root, s, 'outs', 'filtered_feature_bc_matrix.h5')
    if not os.path.exists(candidate):
        candidate = os.path.join(h5_root, s, 'filtered_feature_bc_matrix.h5')
    if not os.path.exists(candidate):
        print(f"Skipping {s}: no file found")
        continue

    sample_adata = sc.read_10x_h5(candidate)

    for col in ('gene_symbols', 'gene_names', 'gene_name', 'symbol'):
        if col in sample_adata.var.columns:
            sample_adata.var_names = sample_adata.var[col].astype(str)
            break

    sample_adata.var_names_make_unique()
    sample_adata.obs_names = [f"{s}-{bc}" for bc in sample_adata.obs_names]
    sample_adata.obs['sample'] = s
    adatas.append(sample_adata)

if len(adatas) == 0:
    raise RuntimeError("No samples read!")

# %% Cell 3
adata = sc.concat(adatas, join='outer', fill_value=0)
adata.var_names_make_unique()

# %% Cell 5
marker_genes = ['GAPDH', 'ACTB', 'CD4', 'CD8A', 'MKI67', 'PRDM1']
present_genes = [g for g in marker_genes if g in adata.var_names]
missing_genes = [g for g in marker_genes if g not in adata.var_names]

print("Present:", present_genes)
print("Missing:", missing_genes)

# %% Cell 6
sample_name_map = {
    "Pla_HDBR10917733": "age_5week_Hrv100_1",
    "Pla_HDBR10701667": "age_6week_Hrv43_4",
    "Pla_HDBR10142769": "age_9week_Hrv46_1",
    "Pla_HDBR10917730": "age_9week_Hrv98_1",
    "Pla_HDBR10701668": "age_9week_Hrv46_3",
    "Pla_HDBR10917731": "age_8week_Hrv99_1",
    "Pla_HDBR10142770": "age_9week_Hrv46_2",
    "Pla_HDBR10701666": "age_6week_Hrv43_1",
    "Pla_HDBR10142767": "age_6week_Hrv43_2",
    "Pla_HDBR10142768": "age_6week_Hrv43_3",
    "Pla_HDBR8624431": "age_12week_H2_2",
    "Pla_HDBR8624430": "age_12week_H2_1",
    "Pla_HDBR8715514": "age_12_13week_H7a9_2",
    "Pla_HDBR8715512": "age_12_13week_H7a9_1",
}

adata.obs["name"] = adata.obs["sample"].map(sample_name_map)
print(adata.obs[["sample", "name"]].drop_duplicates().sort_values("sample"))

# %% Cell 8
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

# %% Cell 9
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

# %% Cell 10
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True)

# %% Cell 11
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# %% Cell 12
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)

# %% Cell 13
import scrublet as scr
import scipy.sparse as sp
import traceback
import sys
print(sys.executable)
import importlib.util
print(importlib.util.find_spec('scrublet'))

# %% Cell 14
counts_matrix = adata.X.copy()
scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
adata.obs['doublet_score'] = doublet_scores
adata.obs['predicted_doublet'] = predicted_doublets
adata = adata[~adata.obs['predicted_doublet'], :].copy()
print('After doublet removal, AnnData shape:', adata.shape)

# %% Cell 15
adata.layers["counts"] = adata.X.copy()

# %% Cell 16
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# %% Cell 18
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", batch_key='sample')

# %% Cell 20
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True, flavor="seurat_v3", batch_key='sample')

# %% Cell 21
genes_to_check = ["GAPDH", "ACTB", "CD4", "CD8A", "PRDM1", "MKI67"]
vn = adata.var_names.astype(str)

def find_matches(gene):
    if gene in vn:
        return [gene]
    matches = [x for x in vn if x.lower() == gene.lower()]
    if matches:
        return matches
    import re
    stripped = [re.sub(r'\.\d+$', '', x) for x in vn]
    matches = [vn[i] for i, x in enumerate(stripped) if x.lower() == gene.lower()]
    if matches:
        return matches
    return [x for x in vn if gene.lower() in x.lower()]

for g in genes_to_check:
    m = find_matches(g)
    print(g, "->", ("FOUND: " + str(m)) if m else "NOT FOUND")

# %% Cell 22
[gene for gene in adata.var_names if 'GAPDH' in gene]

# %% Cell 23
sc.pl.highly_variable_genes(adata)

# %% Cell 24
sc.tl.pca(adata)

# %% Cell 25
sc.pl.pca_variance_ratio(adata, log=True)

# %% Cell 26
sc.pl.pca(adata, color=["sample", "n_genes_by_counts", "total_counts", "pct_counts_mt"], components=["1,2", "3,4"], size=20)

# %% Cell 27
sc.pp.neighbors(adata)

# %% Cell 28
sc.tl.umap(adata)

# %% Cell 29
sc.pl.umap(adata, color="sample", size=2)

# %% Cell 30
sc.pl.umap(adata, color="name", size=2)

# %% Cell 31
sc.tl.leiden(adata, n_iterations=2)

# %% Cell 32
sc.pl.umap(adata, color=["leiden"])

# %% Cell 33
adata.obs.columns

# %% Cell 34
adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype("category")
adata.obs["doublet_score"] = adata.obs["doublet_score"].astype("float32")
sc.pl.umap(adata, color=["leiden", "predicted_doublet", "doublet_score"], wspace=0.5, size=3)

# %% Cell 35
sc.pl.umap(adata, color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"], wspace=0.5, ncols=2)

# %% Cell 36
adata.write_h5ad("/mnt/storage/yxh/nhm/scRNA/placenta/combined_processed.h5ad")

# %% Cell 39
for res in [0.02, 0.5, 2.0]:
    sc.tl.leiden(adata, key_added=f'leiden_res_{res:4.2f}', resolution=res)

# %% Cell 40
sc.pl.umap(adata, color=[f'leiden_res_0.02', "leiden_res_0.50", "leiden_res_2.00"], wspace=0.5)

# %% Cell 42
marker_genes = {
    "CD14+ Mono": ["FCN1", "CD14"],
    "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
    "cDC2": ["CST3", "COTL1", "LYZ", "DMXL2", "CLEC10A", "FCER1A"],
    "Erythroblast": ["MKI67", "HBA1", "HBB"],
    "Proerythroblast": ["CDK6", "SYNGR1", "HBM", "GYPA"],
    "NK": ["GNLY", "NKG7", "CD247", "FCER1G", "TYROBP", "KLRG1", "FCGR3A"],
    "ILC": ["ID2", "PLCG2", "GNLY", "SYNE1"],
    "Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
    "B cells": ["MS4A1", "ITGB1", "COL4A4", "PRDM1", "IRF4", "PAX5", "BCL11A", "BLK", "IGHD", "IGHM"],
    "Plasma cells": ["MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"],
    "Plasmablast": ["XBP1", "PRDM1", "PAX5"],
    "CD4+ T": ["CD4", "IL7R", "TRBC2"],
    "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH", "GZMA"],
    "T naive": ["LEF1", "CCR7", "TCF7"],
    "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],
}

# %% Cell 43
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.02", standard_scale="var")

# %% Cell 45
adata.obs["cell_type_lvl1"] = adata.obs["leiden_res_0.02"].map({"0": "B Cells", "1": "Monocytes", "2": "Erythroid"})

# %% Cell 47
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.50", standard_scale="var")

# %% Cell 49
sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")

# %% Cell 50
sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden_res_0.50", standard_scale="var", n_genes=5)

# %% Cell 51
sc.get.rank_genes_groups_df(adata, group="7").head(5)

# %% Cell 52
dc_cluster_genes = sc.get.rank_genes_groups_df(adata, group="7").head(5)["names"]
sc.pl.umap(adata, color=[*dc_cluster_genes, "leiden_res_0.50"], legend_loc="on data", frameon=False, ncols=3)

# %% Cell 55
adata.obs.columns

# %% Cell 60
adata.obs['cell_type_lvl1'].value_counts()
