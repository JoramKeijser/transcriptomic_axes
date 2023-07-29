import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.linalg import subspace_angles
import argparse
from src import constants
sns.set_context("poster")
sc.settings.figdir= "./figures/figure2/"
sc.settings.dpi_save= 300

species = {'colquitt': "Zebra finch", "tosches": "Turtle", 
            "tasic": "Mouse", "bugeon": "Mouse L1-3", 
            "bakken": "Human"}
def main(args):
    PC_subspace = {}
    for dataset in snakemake.input:
        name = dataset.split("/")[-1].split("h5ad")[0].split("_")[0]
        print(name)
        adata = ad.read_h5ad(dataset)
        PC_subspace[name] = adata.varm['PCs'][:, :constants.NUM_PCS]

    # Now compute principal angles wrt mouse data
    reference = 'tasic'
    remaining = np.sort([dataset for dataset in PC_subspace.keys() if dataset != reference])

    # Add random
    n_genes = PC_subspace[reference].shape[0]
    rng = np.random.RandomState(0)
    idx = rng.permutation(np.arange(n_genes))
    label = 'Chance'
    plt.plot(subspace_angles(PC_subspace[reference][idx], PC_subspace[reference])[::-1] * 180 / np.pi, 
                label = label, color='black', linestyle = ":")
            
    colors = sns.color_palette("Set2")
    for i, name in enumerate(remaining):
        plt.plot(subspace_angles(PC_subspace[name], PC_subspace[reference])[::-1] * 180 / np.pi, 
                label = species[name], color=colors[i+1])

    plt.xlabel("tPC subspace dimension")
    plt.ylabel("Principal angle (deg.)")
    #plt.ylim([45, 90])
    plt.yticks([50, 70, 90])
    plt.legend(handlelength=0.5)
    sns.despine()
    plt.tight_layout()
    plt.savefig(snakemake.output.figure, dpi=300)

    angles = {}
    for (organism, name) in zip(['Human', 'Zebra finch', 'Turtle'], ['bakken', 'colquitt', 'tosches']):
        angles[organism] = subspace_angles(PC_subspace['tasic'], PC_subspace[name])[::-1] * 180 / np.pi
    with open(snakemake.output.angles, 'wb') as handle:
        pickle.dump(angles, handle, protocol=pickle.HIGHEST_PROTOCOL)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--exclude_meis2', action='store_const', default=False, const=True)
    parser.add_argument('--match_abundance', action='store_const', default=False, const=True)
    parser.add_argument('--integrated', action='store_const', default=False, const=True)
    parser.add_argument('--downsample', action='store_const', default=False, const=True)
    args = parser.parse_args()
    main(args)