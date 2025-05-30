import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# DATA

data_path = "/home/c/data/filtered_feature_bc_matrix"
adata = sc.read_10x_mtx(data_path, var_names='gene_symbols', cache=True)
adata.var_names_make_unique()

adata.obs['n_counts'] = adata.X.sum(axis=1).A1
adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1

# Add percent mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# FILTER CELLS / GENES

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs['pct_counts_mt'] < 10, :]

# NORMALIZATION & LOG TRANSFORM

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)

# DIMENSIONALITY REDUCTION

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# NEIGHBORHOOD GRAPH & CLUSTERING

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)  # Use leiden for clustering

# PLOT UMAP
sc.pl.umap(adata, color=['leiden', 'n_genes', 'pct_counts_mt'], save="_basic.png")

# FIND MARKER GENES
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, save="_markers.png")

# SAVE RESULTS
