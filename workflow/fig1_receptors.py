"""
ACh receptor expression vs state modulation
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import argparse
import pickle

import anndata as ad
import scanpy as sc
from scipy.stats import pearsonr
from src import regression_tools 
from statsmodels.stats.multitest import fdrcorrection

sc.set_figure_params(dpi_save=300, figsize=(7,5), fontsize=7, 
    frameon=True, transparent=False, color_map="viridis")
plt.rcParams["axes.grid"] = False
marker_size = 10
alpha = 0.85
fontsize = 6.5
subclass_cmap = sns.color_palette("colorblind", 10)
sns.set_palette("colorblind", 10)
#sns.set_context("poster")
cm = 1/2.54 # centimeters in inches
mpl.rcParams.update({'font.size': 7})


def main(args):
    
    # Load data - avg by subtype because we don't have ACh expression for single cells
    subclass_order = ['Pvalb', 'Sst', 'Lamp5', 'Vip', 'Sncg']
    by_subtype = pd.read_csv(snakemake.input.by_subtype)
    by_subtype['Subclass'] = by_subtype['Subclass'].astype("category")
    by_subtype['Subclass'] = by_subtype['Subclass'].cat.reorder_categories(subclass_order)
    by_subtype.index = by_subtype['Subtype']
    # Load Tasic dataset for ACh receptor expression
    tasic = ad.read_h5ad(snakemake.input.tasic)
    # Consistent coding of subtype names: dash instead of space
    #tasic.obs['Subtype'] = [subtype.replace(" ", "-") for subtype in tasic.obs['subtype']]
    #tasic.obs['Subtype'] = tasic.obs['Subtype'].astype("category")
    # Select bugeon subtypes
    tasic = tasic[tasic.obs['Subtype'].isin(by_subtype['Subtype'])]
    sc.pp.normalize_total(tasic, target_sum=1e4) # counts -> CP10K
    
    # Add ACh receptor genes to obs data frame for convenient plotting
    receptors = ['Chrm3', 'Chrm4', 'Chrna4', 'Chrna5']
    for receptor in receptors:
        if snakemake.params.transform == "log":
            tasic.obs[receptor] = np.log1p(np.array(tasic[:,receptor].X[:,0]))
        else:
            tasic.obs[receptor] = np.array(tasic[:,receptor].X[:,0])
    tasic_by_subtype = tasic.obs.groupby("Subtype").mean()
    by_subtype = pd.concat([by_subtype, tasic_by_subtype[receptors]], axis=1, ignore_index=False)

    # Visualize the 4 receptors from Bugeon et al. 
    fontsize = 8
    modulation = 'State modulation'
    print(modulation)
    significant_receptors = []
    really_significant_receptors = []
    pvals = []
    fig, ax = plt.subplots(1, len(receptors), figsize=(18*cm * len(receptors)/5, 3.5*cm))
    for i, receptor in enumerate(receptors):
        sns.scatterplot(x=receptor, y=modulation, data=by_subtype, 
            hue='Subclass', ax=ax[i],legend=None, alpha=alpha, s=marker_size)
        if i == 0:
            ax[i].set_ylabel(modulation, fontsize=fontsize)
        else:
            ax[i].set_ylabel("")
            ax[i].set_yticklabels([])
        if snakemake.params.transform == "log":
            ax[i].set_xlabel(f"{receptor} (log CP10K)", fontsize=fontsize)
        elif snakemake.params.transform == "raw":
            ax[i].set_xlabel(f"{receptor} (CP10K)",fontsize=fontsize)
        # Show fit (train data)
        x = by_subtype[receptor].to_numpy()[:,None]
        y = by_subtype[modulation].to_numpy()
        R2, x_to_plot, y_to_plot = regression_tools.leave_one_out_lr(x, y)
        ax[i].plot(x_to_plot, y_to_plot, color='black',  alpha=alpha, lw=1)
        # Test data
        corr_test = pearsonr(by_subtype[receptor], by_subtype[modulation])
        ax[i].set_title(f"$R^2$ = {R2:0.2f} r = {corr_test.statistic:0.2f}", fontsize=fontsize)
        print(receptor, corr_test)
        if corr_test.pvalue < 0.05:
            significant_receptors.append(receptor)
        if corr_test.pvalue < 0.05 / len(receptors):
            really_significant_receptors.append(receptor)

    sns.despine()
    fig.tight_layout()
    plt.savefig(snakemake.output.figure, dpi=300)
    print(f"# Significant: {len(significant_receptors)} / {len(receptors)}, {len(really_significant_receptors)} / {len(receptors)}")
    print(really_significant_receptors)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess Bugeon data")
    parser.add_argument("--transform", default = "log", 
                        help="Logarithmize data (log) or not (counts)")
    args = parser.parse_args()

    main(args)

