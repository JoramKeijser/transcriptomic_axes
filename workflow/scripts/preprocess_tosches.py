"""
Put raw Tosches data into anndata
"""
import numpy as np
import pandas as pd
import anndata as ad
from src import data_tools


def main():
    print("Loading data")
    cell_info = pd.read_csv(snakemake.input.metadata, delimiter=" ")
    counts = pd.read_csv(snakemake.input.counts, delimiter=" ")  # gene x samples
    # Select samples for which we have cell info
    cell_names = cell_info.index
    counts = counts.loc[:, np.array(cell_names)]

    adata = ad.AnnData(X=counts.T, dtype=int, obs=cell_info)
    # Select ins
    adata = adata[[cluster.startswith("i") for cluster in adata.obs["clusters"]]]
    # Mouse gene names
    adata.var_names = [gene[0] + gene[1:].lower() for gene in adata.var_names]
    # Map to putative mammalian homologues (Tosches Fig 5)
    naming_map = {}
    for cluster in range(1, 7):
        naming_map[f"i0{cluster}"] = "Meis2"
    for cluster in range(7, 11):
        if cluster < 10:
            naming_map[f"i0{cluster}"] = "Sst"
        else:
            naming_map[f"i{cluster}"] = "Sst"
    for cluster in range(11, 14):
        naming_map[f"i{cluster}"] = "Pvalb"
    for cluster in range(14, 19):
        naming_map[f"i{cluster}"] = "Vip"
    # Consistent order
    adata.obs["Subclass"] = adata.obs["clusters"].map(naming_map)
    adata = data_tools.organize_subclass_labels(adata)

    print("Save to ", snakemake.output.anndata)
    adata.write_h5ad(snakemake.output.anndata)


if __name__ == "__main__":
    main()
