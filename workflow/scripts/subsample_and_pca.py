"""
Apply PCA to single datasets
"""
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from src import constants, pca_tools, data_tools

sns.set_palette("colorblind")
sns.set_context("poster")
sc.settings.dpi_save = 300


# Load data
shared_genes = np.loadtxt(snakemake.input.shared_genes, dtype=str)
adata = ad.read_h5ad(snakemake.input.raw_anndata)
# Subset cells and genes
n_total_cells = adata.shape[0]
n_sample_cells = 640
rng = np.random.RandomState(snakemake.params.seed)
sampled_cells = rng.choice(
    n_total_cells, snakemake.params.n_sample_cells, replace=False
)
adata = adata[sampled_cells, np.sort(adata.var_names.intersection(shared_genes))]

# Preprocess & PCA
sc.pp.normalize_total(adata, target_sum=constants.NORMALIZE_TARGET_SUM)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=constants.NUM_HVG_GENES)
sc.pp.pca(adata, n_comps=constants.NUM_PCS)


# Preserve the order
order = ["Pvalb", "Sst", "Lamp5", "Vip", "Sncg", "Meis2"]
adata.obs["Subclass"] = adata.obs["Subclass"].astype("category")
missing_subclasses = [
    subclass
    for subclass in order
    if subclass not in adata.obs["Subclass"].cat.categories
]
adata.obs["Subclass"] = adata.obs["Subclass"].cat.add_categories(missing_subclasses)
adata.obs["Subclass"] = adata.obs["Subclass"].cat.reorder_categories(order)
# Choose (arbitrary) sign of PCs consistenly across datasets
if np.sum(adata.obs["Subclass"] == "Meis2") > 0:
    if (
        adata[adata.obs["Subclass"] == "Pvalb"].obsm["X_pca"][:, 0].mean()
        > adata[adata.obs["Subclass"] == "Meis2"].obsm["X_pca"][:, 0].mean()
    ):
        adata.obsm["X_pca"][:, 0] *= -1
        adata.varm["PCs"][:, 0] *= -1

adata = pca_tools.orient_axes(adata)
if "colquitt" in snakemake.input.raw_anndata:
    # Flip 2nd axis
    adata.obsm["X_pca"][:, 1] *= -1
    adata.varm["PCs"][:, 1] *= -1
print(f"First 2 tPCs capture {adata.uns['pca']['variance_ratio'][:2].sum()*100:.1f}%")


fig, ax = plt.subplots()
sc.pl.pca(
    adata,
    color="Subclass",
    ax=ax,
    frameon=False,
    annotate_var_explained=True,
    show=False,  # size=size,
    legend_loc="on data",
    legend_fontsize=20,
    save=False,
    title=snakemake.params.species,
)
fig.tight_layout()
plt.savefig(snakemake.output.figure)

adata.write_h5ad(snakemake.output.anndata)
