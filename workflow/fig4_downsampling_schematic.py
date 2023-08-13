# Make schematic to illustrate downsampling
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import argparse
import anndata as ad
import scanpy as sc

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
    sns.set_palette("Set2")
    print("Load data")
    names = ["tasic", "bakken"]
    datasets = {}
    seed = 1747
    rng = np.random.RandomState(seed)
    # for name in names:
    for name, file in snakemake.input.items():
        print(name)
        datasets[name] = ad.read_h5ad(file)
        datasets[name].obs["total_counts"] = datasets[name].X.sum(1)
        datasets[name].obs["n_genes"] = (datasets[name].X > 0).sum(1)
        print(
            f"# counts: {np.median(datasets[name].obs['total_counts'])}, # genes; {np.median(datasets[name].obs['n_genes'])}"
        )

    # Make the plot
    def plot_loghist(x, bins, color=None, alpha=1, histtype="bar"):
        hist, bins = np.histogram(x, bins=bins)
        logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
        if color is None:
            plt.hist(x, bins=logbins, alpha=alpha, histtype=histtype)
        else:
            plt.hist(x, bins=logbins, color=color, alpha=alpha, histtype=histtype)
        plt.xscale("log")

    _, mouse_bins = np.histogram(datasets["tasic"].obs["total_counts"], bins=100)
    plot_loghist(datasets["tasic"].obs["total_counts"], bins=mouse_bins, alpha=0.5)
    plot_loghist(
        datasets["tasic"].obs["total_counts"],
        bins=mouse_bins,
        alpha=1.0,
        histtype="step",
        color=sns.color_palette()[0],
    )
    _, bins = np.histogram(datasets["bakken"].obs["total_counts"], bins=100)
    plot_loghist(datasets["bakken"].obs["total_counts"], bins=bins, alpha=0.5)
    sns.despine()
    plt.text(
        np.median(datasets["tasic"].obs["total_counts"]) * 0.3,
        1000,
        "mouse",
        color=sns.color_palette()[0],
    )
    plt.text(
        np.median(datasets["bakken"].obs["total_counts"]) * 0.3,
        1000,
        "human",
        color=sns.color_palette()[1],
    )
    plt.arrow(
        7.5e5,
        500,
        -6.5e5,
        0,
        head_width=52,
        head_length=3.2e4,
        facecolor="black",
        width=10,
    )
    plt.text(3.5e4, 600, "downsample", fontsize=20)
    plt.yticks([0, 500, 1000])
    plt.xlabel("RNA count")
    plt.ylabel("# cells")
    sns.despine()

    # Add the downsampled mouse data
    x = np.array(datasets["tasic"].X.copy(), dtype=int)
    depth = np.median(datasets["bakken"].obs["total_counts"])
    relative_depth = depth / np.median(datasets["tasic"].obs["total_counts"])
    print(relative_depth)
    x_sub = rng.binomial(x, p=relative_depth)
    subsampled_counts = x_sub.sum(1)
    _, bins = np.histogram(subsampled_counts, bins=100)
    plot_loghist(
        subsampled_counts,
        bins=bins,
        histtype="step",
        color=sns.color_palette()[0],
        alpha=1.0,
    )
    plt.tight_layout()
    plt.savefig(snakemake.output.figure, dpi=300)


if __name__ == "__main__":
    main()
