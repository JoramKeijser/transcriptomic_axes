"""
Subsample cells down to cell count of smallest
dataset
"""
import numpy as np

DATASETS = ["tasic", "tosches", "bakken", "colquitt"]
SEEDS = np.arange(10)  # Sample so many times at random
MIN_CELLS = 640  # number of cells from smallest dataset (could )
species = {
    "colquitt": "Zebra finch",
    "tosches": "Turtle",
    "tasic": "Mouse",
    "bugeon": "Mouse L1-3",
    "bakken": "Human",
    "yao": "mouse",
    "hodge": "Human MTG",
}


rule all:
    input:
        expand(
            "figures/subsampling/pca_{dataset}_{seed}.png",
            dataset=DATASETS,
            seed=SEEDS,
        ),
        "figures/subsampling/principal_angles_subsample.png",
        "figures/subsampling/cross_variance_subsample.png",


rule subsample_and_pca:
    input:
        raw_anndata="data/anndata/{dataset}.h5ad",
        shared_genes="results/gene_lists/shared_genes.txt",
    params:
        species=lambda wildcards: species[wildcards.dataset],
        seed=lambda wildcards: int(wildcards.seed),
        n_sample_cells=MIN_CELLS,
    output:
        figure="figures/subsampling/pca_{dataset}_{seed}.png",
        anndata="results/subsampling/{dataset}_{seed}.h5ad",
    script:
        "../scripts/subsample_and_pca.py"


rule principal_angles:
    input:
        expand(
            "results/subsampling/{dataset}_{seed}.h5ad", dataset=DATASETS, seed=SEEDS
        ),
    params:
        areas=species,
        reference="tasic",
    output:
        figure="figures/subsampling/principal_angles_subsample.png",
    script:
        "../scripts/subsample_principal_angles.py"


rule cross_variance:
    input:
        expand(
            "results/subsampling/{dataset}_{seed}.h5ad", dataset=DATASETS, seed=SEEDS
        ),
    params:
        areas=species,
        reference="tasic",
    output:
        figure="figures/subsampling/cross_variance_subsample.png",
    script:
        "../scripts/subsample_cross_variance.py"
