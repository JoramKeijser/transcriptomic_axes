"""
Description of each dataset
* Species
* Areas (not here?)
* Technology
* # neurons
* sequencing depth
* # genes
"""
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns


adata = ad.read_h5ad(snakemake.input.anndata)
dataset = snakemake.params.dataset
n_neurons = adata.shape[0]
depth = np.median(np.sum(adata.X, axis=1), axis=0)  # reads / cell
num_genes = np.median(np.sum(adata.X > 0, axis=1), axis=0)  # genes / cell
print(f"{dataset}: n = {n_neurons}, depth = {depth}, genes = {num_genes}")
stats = [dataset, n_neurons, depth, num_genes]
adata.obs["Sequencing depth"] = depth
adata.obs["Number of genes"] = num_genes
# CQ plot: depth, genes
fig, ax = plt.subplots(3, 1, figsize=(4, 9))
for i, feature in enumerate(["Sequencing depth", "Number of genes"]):
    sns.violinplot(
        y=feature, data=adata.obs, ax=ax[i], inner=None, linewidth=0, alpha=0.9
    )
    sns.stripplot(y=feature, data=adata.obs, ax=ax[i], s=2, color="gray", alpha=0.5)
    # One vs the other
ax[0].set_title(dataset[0].upper() + dataset[1:] + " et al.")
ax[2].scatter(
    adata.obs["Sequencing depth"],
    adata.obs["Number of genes"],
    s=1,
    color="gray",
    alpha=0.25,
)
ax[2].set_xlabel("Sequencing depth")
ax[2].set_ylabel("Number of genes")
fig.tight_layout()
sns.despine()

df = pd.DataFrame(np.array(stats)[None], columns=["name", "neurons", "depth", "genes"])
df.to_csv(snakemake.output.table, index=False)
print(np.sort(adata.var_names)[:5])
np.savetxt(snakemake.output.genes, np.sort(adata.var_names), fmt="%s")
