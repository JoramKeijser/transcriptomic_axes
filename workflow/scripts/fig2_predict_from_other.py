# Predict state modulation from PCs of other datasets

import numpy as np
import pandas as pd
import argparse
import anndata as ad
import scanpy as sc
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

from src import constants, pca_tools,regression_tools
sns.set_palette("colorblind")
sns.set_context("poster")
sc.settings.figdir= "./figures/figure2/"
sc.settings.dpi_save= 300

def main():
    alpha = 0.6
    fontsize = 20
    marker_size = 25

    datadir = "./data/anndata/"
    bugeon = ad.read_h5ad(snakemake.input.bugeon)
    adata = ad.read_h5ad(snakemake.input.adata)
    species = snakemake.params.species
    # Restrict to shared genes
    adata = adata[:, [gene for gene in bugeon.var_names if gene in adata.var_names]]
    print(f"{adata.shape[1]} shared genes")
    
    # Project Bugeon cells (for which we have state modulation) onto tPC1 of this dataset
    bugeon = bugeon[:, adata.var_names]
    bugeon.obs['tPC1'] = np.array(((bugeon.X - bugeon.X.mean(0)) @ adata.varm['PCs'][:,0])).flatten()
    bugeon.obs['tPC1']
        
    # Now predict
    fig, ax = plt.subplots(figsize=(5, 4))
    sns.scatterplot(x='tPC1', y="state_modulation",  
            data=bugeon.obs, ax=ax, hue="Subclass", legend=None, alpha=alpha, s=marker_size)
    x = bugeon.obs['tPC1'].to_numpy()[:,None]
    y = bugeon.obs["state_modulation"].to_numpy()
    R2, x_to_plot, y_to_plot = regression_tools.leave_one_out_lr(x, y)
    print(f"R^2 subtypes: {R2:0.3f}")
    ax.plot(x_to_plot, y_to_plot, color='black', alpha=alpha)
    # Test R2
    corr_test = pearsonr(np.squeeze(x), y)
    ax.set_title(f"$R^2$ = {R2:0.2f} r = {corr_test.statistic:0.2f} ", fontsize=fontsize)
    print(corr_test)
    # Label axes
    ax.set_ylabel("State modulation", fontsize=fontsize)
    ax.set_xlabel(species + " tPC1", fontsize=fontsize)
    sns.despine()
    fig.tight_layout()
    plt.savefig(snakemake.output.regression, dpi=300)


if __name__ == "__main__":
    main()
