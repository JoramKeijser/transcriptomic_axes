# Mouse variance explained by other datasets
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import os
import pickle
import argparse
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from src.pca_tools import compare_variance
from src import constants
sns.set_palette("Set2")
sns.set_context("poster")
sc.settings.figdir= "./figures/figure4/"
sc.settings.dpi_save= 300


def main(args):
   
    # Do the reference
    areas = snakemake.params.areas
    reference = snakemake.params.reference 
   
    datasets = {}
    print("Loading data")
    for dataset in snakemake.input:
        name = dataset.split("/")[-1].split("h5ad")[0].split("_")[0]
        print(name)
        datasets[name] = ad.read_h5ad(dataset)

    remaining_datasets = [name for name in list(datasets.keys()) if name != reference]
    print(remaining_datasets)
        
    # PCA on reference
    # TODO: NUM_PCs and HVGs to config file
    hvgs = datasets[reference].var.highly_variable
    datasets[reference] = datasets[reference] [:, hvgs]
    pca  = PCA(n_components=constants.NUM_PCS).fit(datasets[reference].X)
    C = pca.get_covariance()
    w = pca.components_[0]
    variance = w@C@w / np.trace(C)


    # Another pass 
    print("Compute cross-covariance")
    cross_variance = {}
    for name, dataset in datasets.items():
        if name != reference:
            # Focus on HVGs of reference
            datasets[name] = datasets[name][:, hvgs]
            cross_variance[name] = compare_variance(datasets[reference], datasets[name], pca)
            print(f"{name}: {cross_variance[name][0]*100:.1f}% variance of mouse tPC1")
            #TODO: log

    variance = np.diag(pca.components_[:10]@C@pca.components_[:10].T) / np.trace(C)
    variance /= variance[0] 

    # Plot
    x = np.arange(5)*4
    n_components = 10
    # Chance level: random components
    rng = np.random.RandomState(0)
    random_components = rng.normal(size=(n_components, C.shape[0]))
    random_components /= np.linalg.norm(random_components, axis=1, keepdims=True)
    random_variance = np.diag(random_components[:10]@C@random_components[:10].T) / float(pca.components_[0]@C@pca.components_[0])
    print(f"random tPC1: {random_variance[0] * 100:0.1f}% of mouse tPC1")
    print(f"random tPC1-10: {np.sum(random_variance) * 100:0.1f}% of mouse tPC1")
    plt.figure()
    plt.hlines(random_variance[0], -1, 20, color='black', linestyles=":", lw = 3, label = 'Chance')
    print(f"tPC1-10: {variance[:10].sum()*100:0.1f}% of reference tPC1")
    plt.bar(x, variance[:5], label = areas[f'{[key for key in areas.keys()][0]}'])
    
    for i, name in enumerate(remaining_datasets): # TODO: don't hard code
        plt.bar(x+.8*(i+1), cross_variance[name][:5], label = areas[name])
        print(f"{name} tPC1-10: {cross_variance[name][:10].sum()*100:0.1f}% of reference tPC1")
    plt.xticks(x+0.4, np.arange(5)+1)
    plt.yticks([0, 1])
    sns.despine()
    plt.ylabel("Variance (norm.)")
    plt.xlabel("tPC")
    plt.ylim([-0.04, 1])
    plt.legend(handlelength=1)
    plt.tight_layout()
    # Save figure 
    plt.savefig(snakemake.output.figure, dpi=300)

    with open(snakemake.output.data, 'wb') as handle:
        pickle.dump(cross_variance, handle, protocol=pickle.HIGHEST_PROTOCOL)
   

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--human', action='store_const', default=False, const=True)
    args = parser.parse_args()
    main(args)