"""
Extract raw Colquitt data and put it into an annotated data frame
"""
import pandas as pd
import anndata as ad
from src import data_tools


def main():
    print("Loading data")
    df = pd.read_csv(snakemake.input[0], index_col=0)
    # First 10 columns are cell metadata, remaining columns are counts
    n_cell_info = 10
    counts = df.iloc[:, n_cell_info:]
    cell_info = df.iloc[:, :n_cell_info]
    # Matched index
    cell_info.index = cell_info["cell"]
    counts.index = cell_info["cell"]

    print("Creating anndata")
    adata = ad.AnnData(X=counts, dtype=int, obs=cell_info)
    adata.var_names = [gene[0] + gene[1:].lower() for gene in counts.columns]
    adata = adata[(adata.obs["species"] == "zf")]  # subset species
    adata = adata[
        ["GABA" in cluster for cluster in adata.obs["cluster_int_sub2"]]
    ]  # subset GABA
    # Remove clusters that correspond to subcortical mammalian neurons
    adata = adata[
        [
            cluster not in ["GABA-7", "GABA-8", "GABA-Pre"]
            for cluster in adata.obs["cluster_int_sub2"]
        ]
    ]
    # Annotate clusters with putative mammalian homologues (Colquitt Fig. 4)
    adata.obs["Subclass"] = (
        adata.obs["cluster_int_sub2"]
        .map(
            {
                "GABA-1-1": "Meis2",
                "GABA-1-2": "Meis2",
                "GABA-2": "Sst",
                "GABA-3": "Pvalb",
                "GABA-4": "Pvalb",
                "GABA-5-1": "Vip",
                "GABA-5-2": "Vip",
                "GABA-5-3": "Vip",
                "GABA-6": "Lamp5",
            }
        )
        .astype("category")
    )
    adata = data_tools.organize_subclass_labels(adata)
    print(adata.obs["Subclass"].value_counts())

    print("Save to ", snakemake.output.anndata)
    adata.write_h5ad(snakemake.output.anndata)


if __name__ == "__main__":
    main()
