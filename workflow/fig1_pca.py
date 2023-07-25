"""
Preprocess Bugeon et al. data for Fig 1 analyses
Assumes we've previously run fig1_extract.py to 
put the raw data into an annotated data frame. 
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import argparse
import anndata as ad
import scanpy as sc
from scipy.stats import pearsonr
from src import regression_tools, constants, pca_tools

plt.rcParams["axes.grid"] = False
sc.settings.figdir="./figures/figure1"
sc.settings.fontsize = 10
marker_size = 100
alpha = 0.85
# Reorder colours
colorblind = sns.color_palette("colorblind", 10)
sns.set_palette(colorblind)
sns.set_context("poster")

def main(args):


    # Load data
    bugeon = ad.read_h5ad(snakemake.input.transcriptomics) #"./data/anndata/bugeon.h5ad"
    activity = pd.read_csv(snakemake.input.activity, index_col=0) #"./results/pandas/bugeon_activity.h5ad"
    print(bugeon.shape, activity.shape)
    assert np.alltrue(np.array(bugeon.obs['Subclass']) == np.array(activity['Subclass']))
    assert np.alltrue(np.array(bugeon.obs['Subtype']) == np.array(activity['Subtype']))
    bugeon.obs = activity # simply overwrite with more comprehensive metadata
    subclass_order = subclass_order = ['Pvalb', 'Sst', 'Lamp5', 'Vip', 'Sncg']
    bugeon.obs['Subclass'] = bugeon.obs['Subclass'].astype("category").cat.reorder_categories(subclass_order)

    # log-normalized and do PCA
    sc.pp.normalize_total(bugeon, target_sum=constants.NORMALIZE_TARGET_SUM)
    if snakemake.params.transform == "counts":
        transform = "counts"
        print("Don't apply log-transform but use normalized counts")
        sc.pp.pca(bugeon, n_comps=71)
        bugeon.obsm['X_pca'][:, 1] *= -1
        bugeon.varm['PCs'][:,1] *= -1
    elif snakemake.params.transform == "log":
        transform = "log"
        print("Apply log-transform")
        sc.pp.log1p(bugeon)
        sc.pp.pca(bugeon, n_comps=71)
    else:
        raise NotImplementedError(f"transform {snakemake.params.transform} not implemented")
        
    # Flip sign of first PC for consistent orientation with other datasets
    bugeon.obsm['X_pca'][:, 0] *= -1
    bugeon.varm['PCs'][:,0] *= -1
        
    print("Oriented")
    # Add to obs data frame for later use
    n_pcs = 30
    for pc in range(n_pcs):
        bugeon.obs[f'PC{pc+1}'] = bugeon.obsm['X_pca'][:, pc]

    bugeon.write_h5ad(snakemake.output.anndata) #f"./results/anndata/bugeon.h5ad"
    # Show
    """
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.despine()
    sc.pl.pca(bugeon, color='Subclass', title = "Subclass",
        legend_fontsize=15, size = marker_size, palette = colorblind,
            save=f"_subclass_{transform}.png", show=False, alpha=0.67,
            ax=ax, annotate_var_explained=True, legend_loc='on data')

    # Color by state modulation
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.despine()
    sc.pl.pca(bugeon, color='State modulation', title = "State modulation",
                vmin=-0.25, vmax=0.55, size=marker_size,
                save=f"_modulation_{transform}.png", show=False, 
                ax=ax, annotate_var_explained=True)

    """
    # Save

   
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess Bugeon data")
    parser.add_argument('--raw', help = "Apply log transform to 'counts' ",
                        action='store_const', default=False, const=True)
    args = parser.parse_args()

    main(args)
