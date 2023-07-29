"""
Test all receptors 
TODO
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
sc.settings.figdir="figures/figure1/"
marker_size = 10
alpha = 0.85
fontsize = 50
subclass_cmap = sns.color_palette("colorblind", 10)
sns.set_palette("colorblind", 10)
#sns.set_context("poster")
cm = 1/2.54 # centimeters in inches
mpl.rcParams.update({'font.size': 7})


def main(args):
    
    # Load data - avg by subtype because we don't have ACh expression and state mod for the same cells
    subclass_order = ['Pvalb', 'Sst', 'Lamp5', 'Vip', 'Sncg']
    by_subtype = pd.read_csv(snakemake.input.by_subtype)
    by_subtype['Subclass'] = by_subtype['Subclass'].astype("category")
    by_subtype['Subclass'] = by_subtype['Subclass'].cat.reorder_categories(subclass_order)
    by_subtype.index = by_subtype['Subtype']
    # Load Tasic dataset for ACh receptor expression
    tasic = ad.read_h5ad(anndata.tasic)
    # Consistent coding of subtype names: dash instead of space
    #tasic.obs['Subtype'] = [subtype.replace(" ", "-") for subtype in tasic.obs['subtype']]
    #tasic.obs['Subtype'] = tasic.obs['Subtype'].astype("category")
    # Select bugeon subtypes
    tasic = tasic[tasic.obs['Subtype'].isin(by_subtype['Subtype'])]
    tasic.obs['Subclass'] = tasic.obs['Subclass'].cat.reorder_categories(subclass_order)
    sc.pp.normalize_total(tasic, target_sum=1e4)

    # Sort by expression
    all_receptors = [gene for gene in tasic.var_names if gene.startswith("Chrm") or gene.startswith("Chrn")]
    expression = np.array([tasic[:, receptor].X.sum() for receptor in all_receptors])
    # Threshold?
    # now sort by correlation
    idx = np.argsort(expression)[::-1]
    all_receptors = np.array(all_receptors)[idx]
    expression = expression[idx]
   
    fname = snakemake.output.tracks.split(sc.figdir)[1] 
    sc.pl.tracksplot(tasic, all_receptors, groupby='Subclass', 
                dendrogram=False, log=False, save="_tasic.png", show=False)
    # Combine at subtype level
    obs_df = tasic.obs.copy()
    for receptor in all_receptors:
        if args.transform == "log":
            obs_df[receptor] = np.log1p(np.array(tasic[:,receptor].X[:,0]))
        else:
            obs_df[receptor] = np.array(tasic[:,receptor].X[:,0])
    tasic_by_subtype = obs_df.groupby("Subtype").mean()
    by_subtype = pd.concat([by_subtype, tasic_by_subtype[all_receptors]], axis=1, ignore_index=False)

    print("Computing significance")
    modulation = 'State modulation'
    significant_receptors = []
    R2_perm = {} 
    R2 = {}
    for i, receptor in enumerate(all_receptors):
        R2_perm[receptor] = np.zeros((args.permutations, ))
        x = by_subtype[receptor].to_numpy()[:,None]
        y = by_subtype[modulation].to_numpy()
        R2[receptor], _, _ = regression_tools.leave_one_out_lr(x, y)
        for p in range(args.permutations):
            R2_perm[receptor][p], _, _ = regression_tools.leave_one_out_lr(x, np.random.permutation(y))
        pvalue = np.mean(R2_perm[receptor] >= R2[receptor]) 
        expr = tasic[:,receptor].X.sum()
        if pvalue < 0.05:
            significant_receptors.append(receptor)
            print(f"{receptor}: p = {pvalue:0.4f}*, expr = {expr}")
        else:
            print(f"{receptor}: p = {pvalue:0.4f} (n.s.), expr = {expr}")
    print(significant_receptors)
    significant_receptors = np.array(significant_receptors)

    # Make tracksplot of significant receptors
    # Sort by 
    significant_expression = np.array([tasic[:, receptor].X.sum() for receptor in significant_receptors])
    significant_receptors = significant_receptors[np.argsort(significant_expression)[::-1]]
    not_significant = np.array([receptor for receptor in all_receptors if receptor not in significant_receptors])
    not_significant_expression = np.array([tasic[:, receptor].X.sum() for receptor in not_significant])
    not_significant = not_significant[np.argsort(not_significant_expression)[::-1]]
    all_receptors = np.concatenate((significant_receptors, not_significant))
    all_expression = np.concatenate((significant_expression, not_significant_expression))
    # Plot those with expr >= 1 CP10K
    sc.pl.tracksplot(tasic, all_receptors[all_expression >= 1.0], groupby='Subclass', 
                dendrogram=False, log=False, save="_tasic_significant_first.png", show=False)

    sns.set_context("poster")
    fontsize = 20
    # Next: plot examples of state modulation. Equal number predictive and not predictive
    n_significant = len(significant_receptors)
    for receptor_list, listname in zip([significant_receptors, not_significant[:n_significant]], ['predictive', 'unpredictive']):
        fig, ax = plt.subplots(n_significant, 2, figsize=(10, 15))
        for row, receptor in enumerate(receptor_list):
            x = by_subtype[receptor].to_numpy()[:,None]
            y = by_subtype[modulation].to_numpy()
            _, x_to_plot, y_to_plot = regression_tools.leave_one_out_lr(x, y)
            # Show fit (train data)
            ax[row, 0].plot(x_to_plot, y_to_plot, color='black',  alpha=alpha, lw=3)     
            # Expression vs modulation
            sns.scatterplot(x=receptor, y=modulation, data=by_subtype, hue='Subclass', ax=ax[row,0], 
                            legend=None, alpha=0.75, s=100)
            if row == 0:
                ax[row,0].set_ylabel(modulation,fontsize=fontsize)
                ax[row,1].set_ylabel("# shuffles", fontsize=fontsize)
                ax[row,1].set_xlabel("$R^2$", fontsize=fontsize)

                if args.transform == "log":
                    ax[row,0].set_xlabel(f"Expr. (log CP10K)", fontsize=fontsize)
                else:
                    ax[row,0].set_xlabel(f"Expr. (CP10K)", fontsize=fontsize)
            else:
                ax[row,0].set_xlabel("")
                ax[row,0].set_ylabel("")
                ax[row,1].set_xlabel("")
                ax[row,1].set_ylabel("")
               
            # Hist of R2s 
            ax[row,1].hist(R2_perm[receptor], bins = int(0.1 * args.permutations), color='gray')
            ax[row,1].vlines(R2[receptor], 0, .1*args.permutations, linestyles=":", color='tab:red', lw=3)
            ax[row,0].set_title(receptor)
            pvalue = np.mean(R2_perm[receptor] >= R2[receptor]) 
            ax[row, 1].set_title(f"$R^2$ = {R2[receptor]:0.2f}")
        sns.despine()
        fig.tight_layout()
        plt.savefig(f"./figures/figure1/receptors_supp_{args.transform}_examples_{listname}_p{args.permutations}", dpi=300)


        # Violin plot
        for i, (receptor, R2vals) in enumerate(R2_perm.items()):
            if i == 0:
                R2_df = pd.DataFrame()
                R2_df['R2'] = R2vals
                R2_df['receptor'] = receptor
            else:
                R2_df2 = pd.DataFrame()
                R2_df2['R2'] = np.array(R2vals)
                R2_df2['receptor'] = receptor
                R2_df = pd.concat((R2_df, R2_df2), axis=0)
        sns.set_context("poster")
        fig, ax = plt.subplots(figsize=(7,5))
        significant_receptors = np.sort(significant_receptors)
        R2_df_significant = R2_df[[receptor in list(significant_receptors) for receptor in R2_df['receptor']]]
        sns.stripplot(y='receptor', x='R2', color='gray', ax=ax, alpha=0.6,
                    data = R2_df_significant)

        for i, receptor in enumerate(significant_receptors):
            plt.scatter(R2[receptor], i, color=sns.color_palette()[3], marker='+', s=400)
            pval = np.mean((R2_df.loc[R2_df['receptor'] == receptor, "R2"] >= R2[receptor]))
            print(i, receptor, pval)
        ax.set_xlim([-0.3,0.3])
        ax.set_xticks([0, 0.3])
        ax.set_ylabel("")
        ax.set_yticks(np.arange(len(significant_receptors)))
        ax.set_yticklabels(significant_receptors)#, rotation=45)
        ax.set_xlabel("Cross-validated $R^2$")
        sns.despine()
        fig.tight_layout()
        plt.savefig(f"./figures/figure1/receptors_supp_{args.transform}_violin_{listname}_p{args.permutations}", dpi=300)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess Bugeon data")
    parser.add_argument("--transform", default = "log", 
                        help="Logarithmize data (log) or not (counts)")
    parser.add_argument("--permutations", default=1000, type=int, 
                        help = "Number of permutations for significance testing")
    args = parser.parse_args()

    main(args)

