# Directly compare principal angles & var explained
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import anndata as ad
import scanpy as sc
import pickle

sns.set_context("poster")


def main():
    colors = sns.color_palette("Set2")
    show = 5
    # Original data

    with open(snakemake.input.angles_complete, "rb") as handle:
        angles = pickle.load(handle)

    with open(snakemake.input.cv_complete, "rb") as handle:
        crossvariance = pickle.load(handle)

    with open(snakemake.input.angles_control, "rb") as handle:
        angles_condition = pickle.load(handle)

    with open(snakemake.input.cv_control, "rb") as handle:
        crossvariance_condition = pickle.load(handle)

    plt.figure()  # Do the angles
    plt.xlabel("Angles")
    for i, key in enumerate(["Human", "Zebra finch", "Turtle"]):
        plt.scatter(
            angles[key][:show],
            angles_condition[key][:show],
            label=key,
            s=100,
            alpha=1,
            color=colors[i + 1],
        )
        for i, (x, y) in enumerate(
            zip(angles[key][:show], angles_condition[key][:show])
        ):
            plt.text(x, y, i + 1, fontsize=15)
    # plt.legend(bbox_to_anchor=(1,1), fontsize = 20)
    plt.axis("equal")
    if snakemake.wildcards.control == "integrated_rpca":
        ticks = np.arange(0, 105, 15)
    else:
        ticks = np.arange(45, 105, 15)
    print("ticks", ticks)
    plt.xticks(ticks, fontsize=20)
    plt.yticks(ticks, fontsize=20)
    # identity line
    plt.plot([ticks[0], ticks[-1]], [ticks[0], ticks[-1]], color="black", linestyle=":")
    sns.despine()
    plt.ylabel("Angles \n" + snakemake.params.control)
    plt.title("Principal angles (deg.)")
    plt.tight_layout()
    plt.savefig(snakemake.output.angles, dpi=300)

    plt.figure()
    organisms = {"bakken": "Human", "colquitt": "Zebra finch", "tosches": "Turtle"}
    for i, key in enumerate(crossvariance.keys()):
        plt.scatter(
            crossvariance[key][:show] * 100,
            crossvariance_condition[key][:show] * 100,
            label=organisms[key],
            s=60,
            alpha=1,
            color=colors[i + 1],
        )
        for i, (x, y) in enumerate(
            zip(crossvariance[key][:show], crossvariance_condition[key][:show])
        ):
            plt.text(x * 100, y * 100, i + 1, fontsize=15)
    plt.axis("equal")
    if snakemake.wildcards.control == "integrated_rpca":
        ticks = np.arange(0, 100, 20)
    else:
        ticks = [0, 10, 20]
    plt.xlabel("% variance")
    plt.xticks(ticks)
    plt.yticks(ticks)
    plt.plot([ticks[0], ticks[-1]], [ticks[0], ticks[-1]], color="black", linestyle=":")
    sns.despine()
    plt.ylabel("% variance. \n" + snakemake.params.control)
    plt.title("Rel. variance explained (%)")
    plt.tight_layout()
    plt.savefig(snakemake.output.cv, dpi=300)


if __name__ == "__main__":
    main()
