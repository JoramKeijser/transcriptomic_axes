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
sc.settings.fontsize = 10
marker_size = 100
alpha = 0.85
# Reorder colours
colorblind = sns.color_palette("colorblind", 10)
sns.set_palette(colorblind)
sns.set_context("poster")

def main(args):

    sc.settings.figdir= "../" #snakemake.input.figdir # "figures/figure1"
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
    if snakemake.params.transform == "raw":
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

    bugeon.write_h5ad(snakemake.output.annotated) 
    # Show
    
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.despine()
    sc.pl.pca(bugeon, color='Subclass', title = "Subclass",
        legend_fontsize=15, size = marker_size, palette = colorblind,
            save=False, show=False, alpha=0.67,
            ax=ax, annotate_var_explained=True, legend_loc='on data')
    fig.tight_layout()
    plt.savefig(snakemake.output.pca_subclass)

    # Color by state modulation
    
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.despine()
    sc.pl.pca(bugeon, color='State modulation', title = "State modulation",
                vmin=-0.25, vmax=0.55, size=marker_size,
                save=False, show=False, 
                ax=ax, annotate_var_explained=True)
    fig.tight_layout()
    plt.savefig(snakemake.output.pca_modulation)
    
    # Group by subtype so we can annotate with receptors from Tasic
    df = bugeon.obs.drop("Subclass", axis="columns")
    by_subtype = bugeon.obs[['Subtype', 'State modulation', 'PC1']].groupby("Subtype").mean()
    #by_subtype = bugeon.obs.groupby("Subtype").mean()
    by_subtype['Subclass'] = [subtype.split("-")[0] for subtype in by_subtype.index]
    by_subtype['Subclass'] = by_subtype['Subclass'].replace({'Serpinf1': 'Vip'})
    by_subtype['Subclass'] = by_subtype['Subclass'].astype("category")
    subclass_order = ['Pvalb', 'Sst', 'Lamp5', 'Vip', 'Sncg']
    by_subtype['Subclass'] = by_subtype['Subclass'].cat.reorder_categories(subclass_order)
    by_subtype.to_csv(snakemake.output.by_subtype)


    print("Predicting modulation")
    inch_to_cm = 1 #2.54
    fontsize = 25
    modulation = 'State modulation'
    print(f"Predict {modulation} from tPCs")
    dropped_na = by_subtype.dropna() 
    dropped_na['Subclass'] = by_subtype['Subclass'].cat.reorder_categories(subclass_order)
    fig, ax = plt.subplots(figsize=(5/inch_to_cm, 4/inch_to_cm))
    sns.scatterplot(x='PC1', y=modulation,  
        data=by_subtype, ax=ax, hue="Subclass", legend=None, alpha=alpha, s=marker_size)
    x = dropped_na['PC1'].to_numpy()[:,None]
    y = dropped_na[modulation].to_numpy()
    R2, x_to_plot, y_to_plot = regression_tools.leave_one_out_lr(x, y)
    print(f"R^2 subtypes: {R2:0.3f}")
    ax.plot(x_to_plot, y_to_plot, color='black', alpha=alpha)
    # Test R2
    corr_test = pearsonr(np.squeeze(x), y)
    ax.set_title(f"$R^2$ = {R2:0.2f} r = {corr_test.statistic:0.2f} ", fontsize=fontsize)
    print(corr_test)
    # Label axes
    ax.set_ylabel(modulation, fontsize=fontsize)
    ax.set_xlabel("tPC1", fontsize=fontsize)
    sns.despine()
    fig.tight_layout()
    plt.savefig(snakemake.output.regression, dpi=300)


    # Individual neurons TODO: log it
    print(f"n = {np.sum(~bugeon.obs[modulation].isna())} neurons")
    idx = ~bugeon.obs[modulation].isna()
    x = bugeon.obs['PC1'].to_numpy()[idx,None]
    y = bugeon.obs[modulation].to_numpy()[idx]
    R2 = regression_tools.leave_one_out_lr(x, y)[0]
    print(f"R^2 individual neurons: {R2:0.2f}")

    #TODO: successive dimensions

   
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess Bugeon data")
    parser.add_argument('--raw', help = "Apply log transform to 'counts' ",
                        action='store_const', default=False, const=True)
    args = parser.parse_args()

    main(args)
