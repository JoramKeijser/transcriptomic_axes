# Compute correlation of all genes with state modulation
# Could still do some GO analyses
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
# and look at conservation
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
from adjustText import adjust_text
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

sc.set_figure_params(
    dpi_save=300,
    figsize=(7, 5),
    fontsize=25,
    frameon=True,
    transparent=False,
    color_map="viridis",
)
plt.rcParams["axes.grid"] = False
marker_size = 150
alpha = 0.67
fontsize = 25
subclass_cmap = sns.color_palette("colorblind", 10)
sns.set_palette("colorblind", 10)
sns.set_context("poster")


def main():
    bugeon = ad.read_h5ad(snakemake.input.bugeon)
    tasic = ad.read_h5ad(snakemake.input.tasic)
    tasic = tasic[:, tasic.X.var(0) > 0].copy()  # can exclude constant genes

    # avg expression for each Bugeon cluster
    subtypes = bugeon.obs["Subtype"].cat.categories
    avg_expr = pd.DataFrame(columns=tasic.var_names, index=subtypes)

    for subtype in subtypes:
        avg_expr.loc[subtype] = tasic[tasic.obs["Subtype"].isin([subtype]), :].X.mean(0)

    # Correlate with state modulation
    x = (
        bugeon.obs[["Subtype", "state_modulation"]]
        .groupby("Subtype")
        .mean()["state_modulation"]
        .to_numpy()
    )
    n_genes = avg_expr.shape[1]
    corrs = np.zeros((n_genes,))
    pvals = np.zeros((n_genes,))
    tests = {}
    for i, gene in enumerate(avg_expr.columns):
        y = avg_expr.iloc[:, i]
        if np.std(y) > 0:
            tests[gene] = pearsonr(x, y)
            corrs[i] = tests[gene].statistic
            pvals[i] = tests[gene].pvalue
        else:
            corrs[i] = 0
            pvals[i] = 1.0

    # Plot
    size = 2
    idx = np.argsort(corrs)
    significant = pvals[idx] < 0.05
    fig, ax = plt.subplots()
    ax.plot(
        np.arange(pvals.shape[0])[~significant],
        corrs[idx][~significant],
        "-",
        color="gray",
        lw=2,
    )
    ax.plot(
        np.arange(pvals.shape[0])[significant],
        corrs[idx][significant],
        "o",
        markersize=size,
        alpha=0.3,
    )

    # Which genes are significant?
    pvals_adj = multipletests(pvals, method="fdr_bh")[1]
    print(avg_expr.columns[pvals_adj < 0.05].shape)
    print(avg_expr.columns[pvals < 0.05])
    np.savetxt(
        snakemake.output.significant_genes,
        np.array(avg_expr.columns[pvals_adj < 0.05]),
        fmt="%s",
    )

    order = np.argsort(corrs)
    dx = 0.7
    chr_positions = []
    chr_corrs = []
    chrs = ["Chrm3", "Chrna4", "Chrna5", "Chrna3", "Chrm4"]
    texts = []

    fontsize = 20
    for i, receptor in enumerate(chrs):
        x = np.where(avg_expr.columns[np.argsort(corrs)] == receptor)[0][0]
        y = tests[receptor].statistic
        texts.append(ax.text(float(x), y, receptor, fontsize=fontsize))
        i = np.where(avg_expr.columns == receptor)[0]
        if y > 0:
            pct = 100 * (1 - x / n_genes)
        else:
            pct = 100 * x / n_genes
        print(f"{receptor}: {x}/{n_genes} = {pct:0.1f}%")  # , p = {pvals[i]:0.4f}")
        plt.scatter(x, y, s=30, color=sns.color_palette()[3])

    adjust_text(
        texts,
        force_points=0.2,
        force_text=0.1,
        expand_points=(1, 1),
        expand_text=(1, 1),
        arrowprops=dict(arrowstyle="-", color="black", lw=1, alpha=0.8),
    )

    plt.hlines(0, 0, n_genes, color="k", linestyles=":", lw=1)
    plt.ylim([-1, 1])
    plt.xticks([1, n_genes])
    sns.despine()
    plt.ylabel("correlation w/ \n state modulation")
    plt.tight_layout()
    plt.savefig(snakemake.output.figure, dpi=300)

    return 0


if __name__ == "__main__":
    main()
