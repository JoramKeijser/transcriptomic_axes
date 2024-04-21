"""
Visualize receptor expression
"""
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_palette("colorblind")
sns.set_context("poster")

receptors = np.loadtxt(snakemake.input.receptors, dtype=str)
shared_genes = np.loadtxt(snakemake.input.shared_genes, dtype=str)
receptors = np.sort(list(set(receptors).intersection(shared_genes)))
print("Shared receptors", receptors)
adata = ad.read_h5ad(snakemake.input.anndata)
# Subset receptors, z-score
print("Receptors")
print(receptors)
adata = adata[:, receptors]
adata.layers["scaled"] = sc.pp.scale(adata, copy=True).X

fig, ax = plt.subplots(figsize=(8, 6))
sc.pl.dotplot(
    adata,
    receptors,
    "Subclass",
    dendrogram=False,
    swap_axes=True,
    ax=ax,
    colorbar_title="mean z-score",
    layer="scaled",
    vmin=-2,
    vmax=2,
    cmap="Blues",
    return_fig=False,
    show=False,
    var_group_rotation=45,
    title=None,
    dot_max=0.60,
)
fig.tight_layout()
fig.savefig(snakemake.output.dotplot, dpi=300)
