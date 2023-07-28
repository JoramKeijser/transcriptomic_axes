# To do: fix Colquitt order. Missing subclasses

import numpy as np
import pandas as pd
import argparse
import anndata as ad
import scanpy as sc
import os
import matplotlib.pyplot as plt
import seaborn as sns
from src import constants, pca_tools
sns.set_palette("colorblind")
sns.set_context("poster")
sc.settings.figdir= "./figures/figure2/"
sc.settings.dpi_save= 300


 def subsample(adata, subclasses = ['Pvalb', 'Sst', 'Vip']):
        counts = adata.obs['Subclass'].value_counts()
        smallest_n = min(counts[subclasses])
        adata_sub = ad.AnnData()
        for i, subclass in enumerate(subclasses):
            idx = np.where(adata.obs["Subclass"] == subclass)[0]
            select = np.random.choice(idx, size = (smallest_n, ), replace=False)
            if i == 0:
                adata_sub = adata[select]
            else:
                adata_sub = ad.concat((adata_sub, adata[select]))
        return adata_sub

#    title = species[savename]
shared_genes = np.loadtxt(snakemake.input.shared_genes, dtype=str)
adata = ad.read_h5ad(snakemake.input.raw_anndata)
adata = adata[:, shared_genes] # restrict to shared genes

# Preserve the order
subclass_order = ['Pvalb', 'Sst', 'Lamp5', 'Vip', 'Sncg', 'Meis2']
adata.obs['Subclass'] = adata.obs['Subclass'].astype("category")
missing_subclasses = [subclass for subclass in subclass_order if subclass not in adata.obs['Subclass'].cat.categories]
adata.obs['Subclass'] = adata.obs['Subclass'].cat.add_categories(missing_subclasses)
adata.obs['Subclass'] = adata.obs['Subclass'].cat.reorder_categories(subclass_order)

sc.pp.normalize_total(adata, target_sum = constants.NORMALIZE_TARGET_SUM)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes = constants.NUM_HVG_GENES)
sc.pp.pca(adata, n_comps = constants.NUM_PCS)
# Choose (arbitrary) sign of PCs consistenly across datasets
if np.sum(adata.obs['Subclass'] == "Meis2") > 0:
    if adata[adata.obs['Subclass'] == "Pvalb"].obsm['X_pca'][:,0].mean() > adata[adata.obs['Subclass'] == "Meis2"].obsm['X_pca'][:,0].mean():
        adata.obsm['X_pca'][:,0] *= -1 
        adata.varm['PCs'][:,0] *= -1 
adata = pca_tools.orient_axes(adata)
print(snakemake.input.raw_anndata)
print(f"First 2 tPCs capture {adata.uns['pca']['variance_ratio'][:2].sum()*100:.1f}%")
fig, ax = plt.subplots()
sns.despine()
#if savename.startswith("tosches"):
#    size = 100 # small sample size -> larger dots
#else:
#    size = 50
fig, ax = plt.subplots()
sc.pl.pca(adata, color='Subclass', ax=ax, frameon=False, annotate_var_explained=True, show=False, #size=size, 
        legend_loc='on data', legend_fontsize=20, save=False, title = snakemake.params.species)
fig.tight_layout()
plt.savefig(snakemake.output.figure)

adata.write_h5ad(snakemake.output.anndata)
    
