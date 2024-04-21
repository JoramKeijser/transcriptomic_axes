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


def subsample(adata, subclasses=["Pvalb", "Sst", "Vip"]):
    counts = adata.obs["Subclass"].value_counts()
    smallest_n = min(counts[subclasses])
    adata_sub = ad.AnnData()
    for i, subclass in enumerate(subclasses):
        idx = np.where(adata.obs["Subclass"] == subclass)[0]
        select = np.random.choice(idx, size=(smallest_n,), replace=False)
        if i == 0:
            adata_sub = adata[select]
        else:
            adata_sub = ad.concat((adata_sub, adata[select]))
    return adata_sub


shared_genes = np.loadtxt(snakemake.input.shared_genes, dtype=str)
adata = ad.read_h5ad(snakemake.input.raw_anndata)
adata = adata[:, np.sort(adata.var_names.intersection(shared_genes))]

if snakemake.params.control == "meis2":
    n = np.sum(adata.obs["Subclass"] == "Meis2")
    print(f"Exclude {n} Meis2 cells")
    adata = adata[adata.obs["Subclass"] != "Meis2"]
elif snakemake.params.control == "abundance":
    adata = subsample(adata)
elif snakemake.params.control == "depth":
    fewest_reads = 3150.5  # Colquitt
    reads = np.median(adata.X.sum(1))
    relative_depth = fewest_reads / reads
    print(f"Depth: {reads}, Relative depth: {relative_depth:1.0E}")
    if relative_depth < 1:
        adata.X = np.random.binomial(np.array(adata.X, dtype=int), p=relative_depth)
elif snakemake.params.control == "bugeonabundance":
    bugeon = ad.read_h5ad(snakemake.input.bugeon)
    bugeon_number = bugeon.obs["Subclass"].value_counts()
    adata_sub = {}
    for subclass in bugeon.obs["Subclass"].cat.categories:
        adata_sub[subclass] = adata[adata.obs["Subclass"] == subclass]
        sc.pp.subsample(adata_sub[subclass], n_obs=bugeon_number[subclass])
    adata = ad.concat(list(adata_sub.values()))
    adata = data_tools.organize_subclass_labels(adata)
elif snakemake.params.control == "bugeonsst":
    # Select only upper layer Sst subtypes, as they are present in the Bugeon dataset
    bugeon = ad.read_h5ad(snakemake.input.bugeon)

    def abundance(adata, relative=True):
        # Compute relative abundance of each subclass
        abundance = adata.obs["Subclass"].value_counts()
        if relative:
            return abundance / abundance.sum()
        return abundance

    n_sst = abundance(adata, False)["Sst"]

    # Split into Sst and non-Sst cells
    sst_ix = adata.obs["Subclass"] == "Sst"
    n_sst = np.sum(sst_ix)
    adata_sst = adata[sst_ix]
    adata_other = adata[~sst_ix]
    # Select upper Sst types
    upper_subtypes = [
        subtype
        for subtype in bugeon.obs["Subtype"].cat.categories
        if subtype.startswith("Sst")
    ]
    adata_sst = adata_sst[
        [subtype in upper_subtypes for subtype in adata_sst.obs["Subtype"]]
    ]
    n_upper = adata_sst.shape[0]

    def subsample(tasic_sst, tasic_other):
        tasic_subsampled = tasic_sst[
            [subtype in upper_subtypes for subtype in tasic_sst.obs["Subtype"]]
        ]
        for subclass in tasic_other.obs["Subclass"].cat.categories:
            ix = tasic_other.obs["Subclass"] == subclass
            t2 = sc.pp.subsample(tasic_other[ix], fraction=n_upper / n_sst, copy=True)
            tasic_subsampled = ad.concat((tasic_subsampled, t2))
        return tasic_subsampled

    adata = subsample(adata_sst, adata_other)

elif snakemake.params.control == "72g":
    bugeon = ad.read_h5ad(snakemake.input.bugeon)
    shared_genes = bugeon.var_names.intersection(adata.var_names)
    adata = adata[:, shared_genes]

if "integrated" not in snakemake.params.control:
    # still need to log norm
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
