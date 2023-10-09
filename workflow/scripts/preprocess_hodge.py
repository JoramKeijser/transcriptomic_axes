# Put raw hodge data into AnnData
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from src import data_tools
from src import constants, pca_tools
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_palette("colorblind")
sns.set_context("poster")
sc.settings.figdir = "./figures/figure4/"
sc.settings.dpi_save = 300


def main():
    metadata = pd.read_csv(snakemake.input.cells, index_col=0)
    counts = pd.read_csv(snakemake.input.exons, index_col=0, header=0, dtype=int)
    counts += pd.read_csv(snakemake.input.introns, index_col=0, header=0, dtype=int)
    print(metadata.shape, counts.shape)
    genes = pd.read_csv(snakemake.input.genes)
    genes = [gene[0] + gene[1:].lower() for gene in genes["gene"]]
    # Put into AnnData, subset GABAergic
    adata = ad.AnnData(X=np.array(counts).T, obs=metadata)
    adata.var_names = genes
    depth = np.median(adata.X.sum(1))  # reads / cell
    num_genes = np.median((adata.X > 0).sum(1))  # genes / cell
    print(f"# counts: {depth}, # genes: {num_genes}")

    adata = adata[adata.obs["class"] == "GABAergic"]
    # Subclass assignment
    order = ["Pvalb", "Sst", "Lamp5", "Vip", "Sncg", "Meis2"]

    def camelcase(name):
        return name[0].upper() + name[1:].lower()

    # Extract subclass from cluster name
    # Exceptions: "Inh L1 Lamp5 NMMBR" from paper is coded in the data as SST
    # Pax6 correspond to Lamp5 cells; Gad1 MC4R to VIP, GAD1 GLP1R to Pvalb
    def subclass_from_cluster(cluster):
        candidate = camelcase(cluster.split(" ")[2])
        if cluster == "Inh L1 SST NMBR" or candidate == "Pax6":
            return "Lamp5"
        elif cluster == "Inh L5-6 GAD1 GLP1R":
            return "Pvalb"
        elif cluster == "Inh L1-2 GAD1 MC4R":
            return "Lamp5"
        else:
            return candidate

    adata.obs["Subclass"] = [subclass_from_cluster(cl) for cl in adata.obs["cluster"]]
    # Usual order with 2 additional human-specific types
    adata.obs["Subclass"] = adata.obs["Subclass"].astype("category")  #
    missing = [
        subclass
        for subclass in order
        if subclass not in np.unique(adata.obs["Subclass"])
    ]
    adata.obs["Subclass"] = adata.obs["Subclass"].cat.add_categories(missing)
    adata.obs["Subclass"] = adata.obs["Subclass"].cat.reorder_categories(order)

    adata.write_h5ad(snakemake.output.anndata)  # TODO: also count introns


if __name__ == "__main__":
    main()
