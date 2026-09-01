from pathlib import Path

import pandas as pd
import scanpy as sc
import yaml


ROOT = Path(".")
PARAMS = yaml.safe_load((ROOT / "configs" / "analysis_params.yaml").read_text())
SAMPLES = pd.read_csv(ROOT / "configs" / "samples.csv")


def read_sample(row):
    adata = sc.read_10x_mtx(row.data_path, var_names="gene_symbols", cache=True)
    adata.var_names_make_unique()
    adata.obs["sample_id"] = row.sample_id
    adata.obs["group"] = row.group
    adata.obs["batch"] = row.batch
    adata.obs_names = [f"{row.sample_id}_{barcode}" for barcode in adata.obs_names]
    return adata


adatas = [read_sample(row) for row in SAMPLES.itertuples(index=False)]
adata = adatas[0].concatenate(
    *adatas[1:],
    batch_key="sample_batch",
    batch_categories=SAMPLES["sample_id"].tolist(),
    index_unique=None,
)

adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

qc = PARAMS["quality_control"]
adata = adata[
    (adata.obs["n_genes_by_counts"] >= qc["min_genes_per_cell"])
    & (adata.obs["n_genes_by_counts"] <= qc["max_genes_per_cell"])
    & (adata.obs["pct_counts_mt"] <= qc["max_mito_percent"])
].copy()

norm = PARAMS["normalization"]
sc.pp.normalize_total(adata, target_sum=norm["scale_factor"])
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=norm["n_highly_variable_genes"])
adata = adata[:, adata.var["highly_variable"]].copy()

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")

dim = PARAMS["dimension_reduction"]
sc.pp.neighbors(adata, n_pcs=dim["n_pcs"])
sc.tl.leiden(adata, resolution=dim["clustering_resolution"])
sc.tl.umap(adata, random_state=dim["random_seed"])

figures_dir = ROOT / PARAMS["outputs"]["figures_dir"]
tables_dir = ROOT / PARAMS["outputs"]["tables_dir"]
objects_dir = ROOT / PARAMS["outputs"]["objects_dir"]

for path in [figures_dir, tables_dir, objects_dir]:
    path.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = figures_dir
sc.pl.umap(adata, color=["leiden", "sample_id", "group"], save="_clusters_samples.png", show=False)

sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv(tables_dir / "scanpy_cluster_markers.csv", index=False)

adata.write_h5ad(objects_dir / "scanpy_processed.h5ad")
