"""
Extract raw Bakken data and put it into an annotated data frame
"""
import pandas as pd
import anndata as ad
from src import data_tools


def main():
    cell_info = pd.read_csv(snakemake.input.metadata)
    print("Loading count data")
    counts = pd.read_csv(snakemake.input.counts)

    cell_info = cell_info.rename(
        columns={"subclass_label": "Subclass", "class_label": "Class"}
    )
    # Subset
    counts = counts[cell_info["Class"] == "GABAergic"]
    cell_info = cell_info[cell_info["Class"] == "GABAergic"]
    cell_info["Subclass"] = cell_info["Subclass"].astype("category")
    # Merge Sst Chodl into Sst, and Serpinf1 into Vip
    cell_info["Subclass"] = cell_info["Subclass"].replace(
        {"Serpinf1": "Vip", "Sst Chodl": "Sst"}
    )

    print("Put everything into annotated data frame")
    adata = ad.AnnData(X=counts.to_numpy()[:, 1:], dtype=int, obs=cell_info)
    adata.var_names = [
        gene[0] + gene[1:].lower() for gene in counts.columns[1:]
    ]  # mouse convention "Elfn1"
    # Add missing subclasses
    adata = data_tools.organize_subclass_labels(adata)
    print(adata.obs["Subclass"].value_counts())

    print("Save to ", snakemake.output.anndata)
    adata.write_h5ad(snakemake.output.anndata)


if __name__ == "__main__":
    main()
