"""
Compute the transcriptomic PCs a la Bugeon et al.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.stats import pearsonr
from src import regression_tools, constants

plt.rcParams["axes.grid"] = False
sc.settings.fontsize = 10
MARKERSIZE = 100
ALPHA = 0.85
# Reorder colours
colorblind = sns.color_palette("colorblind", 10)
sns.set_palette(colorblind)
sns.set_context("poster")


# Load data
bugeon = ad.read_h5ad(snakemake.input.transcriptomics)
activity = pd.read_csv(snakemake.input.activity, index_col=0)
print(bugeon.shape, activity.shape)
assert np.alltrue(np.array(bugeon.obs["Subclass"]) == np.array(activity["Subclass"]))
assert np.alltrue(np.array(bugeon.obs["Subtype"]) == np.array(activity["Subtype"]))
bugeon.obs = activity  # simply overwrite with more comprehensive metadata
subclass_order = subclass_order = ["Pvalb", "Sst", "Lamp5", "Vip", "Sncg"]
bugeon.obs["Subclass"] = (
    bugeon.obs["Subclass"].astype("category").cat.reorder_categories(subclass_order)
)

# log-normalized and do PCA
if snakemake.params.transform == "scaled":
    print("Don't apply log-transform but use normalized counts")
    sc.pp.normalize_total(bugeon, target_sum=constants.NORMALIZE_TARGET_SUM)
    sc.pp.pca(bugeon, n_comps=71)
    bugeon.obsm["X_pca"][:, 1] *= -1
    bugeon.varm["PCs"][:, 1] *= -1
elif snakemake.params.transform == "log":
    print("Apply log-transform")
    sc.pp.normalize_total(bugeon, target_sum=constants.NORMALIZE_TARGET_SUM)
    sc.pp.log1p(bugeon)
    sc.pp.pca(bugeon, n_comps=71)
elif snakemake.params.transform == "raw":
    print("Don't scale, use data as-is")
    sc.pp.pca(bugeon, n_comps=71)
    bugeon.obsm["X_pca"][:, 0] *= -1
    bugeon.varm["PCs"][:, 0] *= -1
else:
    raise NotImplementedError(f"transform {snakemake.params.transform} not implemented")

# Flip sign of first PC for consistent orientation with other datasets
bugeon.obsm["X_pca"][:, 0] *= -1
bugeon.varm["PCs"][:, 0] *= -1

print("Oriented")
# Add to obs data frame for later use
NPCS = 30
for pc in range(NPCS):
    bugeon.obs[f"PC{pc+1}"] = bugeon.obsm["X_pca"][:, pc]

bugeon.write_h5ad(snakemake.output.annotated)
# Show

fig, ax = plt.subplots(figsize=(5, 4))
sns.despine()
sc.pl.pca(
    bugeon,
    color="Subclass",
    title="Subclass",
    legend_fontsize=15,
    size=MARKERSIZE,
    palette=colorblind,
    save=False,
    show=False,
    alpha=0.67,
    ax=ax,
    annotate_var_explained=True,
    legend_loc="on data",
)
fig.tight_layout()
plt.savefig(snakemake.output.pca_subclass)

# Color by state modulation

fig, ax = plt.subplots(figsize=(5, 4))
sns.despine()
sc.pl.pca(
    bugeon,
    color="State modulation",
    title="State modulation",
    vmin=-0.25,
    vmax=0.55,
    size=MARKERSIZE,
    save=False,
    show=False,
    ax=ax,
    annotate_var_explained=True,
)
fig.tight_layout()
plt.savefig(snakemake.output.pca_modulation)

# Group by subtype so we can annotate with receptors from Tasic
# groupby cannot handle/ignore categorical vars, so we need to do some manual work
by_subtype = bugeon.obs[["Subtype", "State modulation"]].groupby("Subtype").mean()
for pc in range(NPCS):
    by_subtype[f"PC{pc+1}"] = (
        bugeon.obs[["Subtype", f"PC{pc+1}"]].groupby("Subtype").mean()
    )

by_subtype["Subclass"] = [subtype.split("-")[0] for subtype in by_subtype.index]
by_subtype["Subclass"] = by_subtype["Subclass"].replace({"Serpinf1": "Vip"})
by_subtype["Subclass"] = by_subtype["Subclass"].astype("category")
subclass_order = ["Pvalb", "Sst", "Lamp5", "Vip", "Sncg"]
by_subtype["Subclass"] = by_subtype["Subclass"].cat.reorder_categories(subclass_order)
by_subtype.to_csv(snakemake.output.by_subtype)


print("Predicting modulation")
FONTSIZE = 22
MODULATION = "State modulation"
print(f"Predict {MODULATION} from tPCs")
dropped_na = by_subtype.dropna()
dropped_na["Subclass"] = by_subtype["Subclass"].cat.reorder_categories(subclass_order)
fig, ax = plt.subplots(figsize=(5, 4))
sns.scatterplot(
    x="PC1",
    y=MODULATION,
    data=by_subtype,
    ax=ax,
    hue="Subclass",
    legend=None,
    alpha=ALPHA,
    s=MARKERSIZE,
)
x = dropped_na["PC1"].to_numpy()[:, None]
y = dropped_na[MODULATION].to_numpy()
R2, x_to_plot, y_to_plot = regression_tools.leave_one_out_lr(x, y)
print(f"R^2 subtypes: {R2:0.3f}")
ax.plot(x_to_plot, y_to_plot, color="black", alpha=ALPHA)
# Test R2
corr_test = pearsonr(np.squeeze(x), y)
ax.set_title(f"$R^2$ = {R2:0.2f} r = {corr_test.statistic:0.2f} ", fontsize=FONTSIZE)
print(corr_test)
# Label axes
ax.set_ylabel(MODULATION, fontsize=FONTSIZE)
ax.set_xlabel("tPC1", fontsize=FONTSIZE)
sns.despine()
fig.tight_layout()
plt.savefig(snakemake.output.regression, dpi=300)


# Individual neurons TODO: log it
print(f"n = {np.sum(~bugeon.obs[MODULATION].isna())} neurons")
idx = ~bugeon.obs[MODULATION].isna()
x = bugeon.obs["PC1"].to_numpy()[idx, None]
y = bugeon.obs[MODULATION].to_numpy()[idx]
R2 = regression_tools.leave_one_out_lr(x, y)[0]
print(f"R^2 individual neurons: {R2:0.2f}")

# Now try successive dimensions
R2s = np.zeros((NPCS,))
for pc in range(1, NPCS):
    x = dropped_na[f"PC{pc+1}"].to_numpy()[:, None]
    y = dropped_na[MODULATION].to_numpy()
    R2s[pc], x_to_plot, y_to_plot = regression_tools.leave_one_out_lr(x, y)
print(f"Best PC > 1: tPC{np.argmax(R2s)+1} with R^2 = {np.max(R2s):0.2f}")
# How much variance
pc1 = bugeon.varm["PCs"][:, 0]
pc_alt = bugeon.varm["PCs"][:, np.argmax(R2s)]
C = np.cov(bugeon.X.T)
fraction_of_var = lambda w: w @ C @ w / np.trace(C) * 100
print(
    f"Its variance: {fraction_of_var(pc_alt):0.1f}% vs. {fraction_of_var(pc1):0.1f}% (tPC1)"
)
