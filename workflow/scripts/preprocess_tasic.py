"""
Extract raw Tasic data and put it into an annotated data frame
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import anndata as ad
import scanpy as sc
from src import data_tools
from src import constants


def main():
    print("Loading data")
    counts = pd.read_csv(snakemake.input.exons, index_col=0)  # genes x cells
    intron_counts = pd.read_csv(snakemake.input.introns, index_col=0)  # genes x cells
    counts += intron_counts  # Combine exons and introns
    gene_info = pd.read_csv(snakemake.input.genes)
    cell_info = pd.read_csv(snakemake.input.cells, index_col=0)

    # Consistent naming of metadata vars
    cell_info = cell_info.rename(
        columns={"cluster": "Subtype", "subclass": "Subclass", "class": "Class"}
    )
    # Dash instead of space
    cell_info["Subtype"] = [
        subtype.replace(" ", "-") for subtype in cell_info["Subtype"]
    ]
    cell_info["Class"] = cell_info["Class"].astype("category")
    cell_info["Subtype"] = cell_info["Subtype"].astype("category")
    cell_info["Subclass"] = cell_info["Subclass"].astype("category")
    # Merge Serpinf1 in Vip as done by Bugeon
    cell_info["Subclass"] = cell_info["Subclass"].replace({"Serpinf1": "Vip"})

    print("Creating anndata")
    adata = ad.AnnData(X=counts.T, dtype=int, obs=cell_info)
    adata.var_names = gene_info["gene_symbol"]
    adata = adata[adata.obs["Class"] == "GABAergic"]

    adata = data_tools.organize_subclass_labels(adata)

    print(adata.obs["Subclass"].value_counts())

    print(f"Save data to {snakemake.output.anndata}")
    adata.write_h5ad(snakemake.output.anndata)

    print("HVGS")
    sc.pp.normalize_total(adata, target_sum=constants.NORMALIZE_TARGET_SUM)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=constants.NUM_HVG_GENES)
    hvgs = list(adata.var_names[adata.var["highly_variable"]])
    genes_path = "./results/gene_lists/"
    # if not os.path.exists(genes_path):
    #     os.mkdir(genes_path)
    fname = genes_path + f"tasic_hvgs_{constants.NUM_HVG_GENES}.txt"
    print(f"Save top {constants.NUM_HVG_GENES} HVGS to ", fname)
    np.savetxt(snakemake.output.hv_genes, hvgs, fmt="%s")


if __name__ == "__main__":
    main()
