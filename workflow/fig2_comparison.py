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
    show = 10
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
    plt.plot([45, 90], [45, 90], color="black", linestyle=":")
    plt.xlabel("Angles (all cells)")
    plt.xticks([50, 70, 90])
    plt.yticks([50, 70, 90])
    for i, key in enumerate(["Human", "Zebra finch", "Turtle"]):
        plt.scatter(
            angles[key][:show],
            angles_condition[key][:show],
            label=key,
            s=100,
            alpha=np.linspace(1, 0.1, show),
            color=colors[i + 1],
        )
        print(key, angles[key][:show].sum(), angles_condition[key][:show].sum())
    # plt.legend()
    sns.despine()
    plt.ylabel("Angles " + snakemake.params.control)
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
            alpha=np.linspace(1, 0, show),
            color=colors[i + 1],
        )
    else:
        plt.plot([0, 20], [0, 20], color="black", linestyle=":")
        plt.xlabel("% var (all cells)")
        plt.xticks([0, 10, 20])
        plt.yticks([0, 10, 20])
    # plt.legend()
    sns.despine()
    plt.ylabel("% var. " + snakemake.params.control)
    plt.title("Rel. variance explained (%)")
    plt.tight_layout()
    plt.savefig(snakemake.output.cv, dpi=300)


if __name__ == "__main__":
    main()
