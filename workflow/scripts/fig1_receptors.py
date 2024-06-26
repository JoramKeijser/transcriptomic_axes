"""
ACh receptor expression vs state modulation
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.stats import pearsonr
from src import regression_tools

sc.set_figure_params(
    dpi_save=300,
    figsize=(7, 5),
    fontsize=7,
    frameon=True,
    transparent=False,
    color_map="viridis",
)
plt.rcParams["axes.grid"] = False
MARKERSIZE = 10
ALPHA = 0.85
FONTSIZE = 6.5
subclass_cmap = sns.color_palette("colorblind", 10)
sns.set_palette("colorblind", 10)
# sns.set_context("poster")
CM2INCH = 1 / 2.54  # centimeters in inches
mpl.rcParams.update({"font.size": 7})


# Load data - avg by subtype because we don't have ACh expression for single cells
subclass_order = ["Pvalb", "Sst", "Lamp5", "Vip", "Sncg"]
by_subtype = pd.read_csv(snakemake.input.by_subtype)
by_subtype["Subclass"] = by_subtype["Subclass"].astype("category")
by_subtype["Subclass"] = by_subtype["Subclass"].cat.reorder_categories(subclass_order)
by_subtype.index = by_subtype["Subtype"]
# Load Tasic dataset for ACh receptor expression
tasic = ad.read_h5ad(snakemake.input.tasic)
# Select bugeon subtypes
tasic = tasic[tasic.obs["Subtype"].isin(by_subtype["Subtype"])]

# Add ACh receptor genes to obs data frame for convenient plotting
receptors = ["Chrm3", "Chrm4", "Chrna4", "Chrna5"]
for receptor in receptors:
    if snakemake.params.transform == "log":
        sc.pp.normalize_total(tasic, target_sum=1e4)  # counts -> CP10K
        tasic.obs[receptor] = np.log1p(np.array(tasic[:, receptor].X[:, 0]))
    elif snakemake.params.transform == "scaled":
        sc.pp.normalize_total(tasic, target_sum=1e4)  # counts -> CP10K
        tasic.obs[receptor] = np.array(tasic[:, receptor].X[:, 0])
    elif snakemake.params.transform == "raw":
        tasic.obs[receptor] = np.array(tasic[:, receptor].X[:, 0])
    else:
        raise NotImplementedError(
            f"transform {snakemake.params.transform} not implemented"
        )

tasic_by_subtype = tasic.obs[receptors + ["Subtype"]].groupby("Subtype").mean()
by_subtype = pd.concat([by_subtype, tasic_by_subtype], axis=1, ignore_index=False)

# Visualize the 4 receptors from Bugeon et al.
FONTSIZE = 8
MODULATION = "State modulation"
significant_receptors = []
really_significant_receptors = []
pvals = []
fig, ax = plt.subplots(
    1, len(receptors), figsize=(16 * CM2INCH * len(receptors) / 5, 3.5 * CM2INCH)
)
for i, receptor in enumerate(receptors):
    sns.scatterplot(
        x=receptor,
        y=MODULATION,
        data=by_subtype,
        hue="Subclass",
        ax=ax[i],
        legend=None,
        alpha=ALPHA,
        s=MARKERSIZE,
    )
    if i == 0:
        ax[i].set_ylabel(MODULATION, fontsize=FONTSIZE)
    else:
        ax[i].set_ylabel("")
        ax[i].set_yticklabels([])
    if snakemake.params.transform == "log":
        ax[i].set_xlabel(f"{receptor} (log CP10K)", fontsize=FONTSIZE)
    elif snakemake.params.transform == "scaled":
        ax[i].set_xlabel(f"{receptor} (CP10K)", fontsize=FONTSIZE)
    elif snakemake.params.transform == "raw":
        ax[i].set_xlabel(f"{receptor} (counts)", fontsize=FONTSIZE)
    # Show fit (train data)
    x = by_subtype[receptor].to_numpy()[:, None]
    y = by_subtype[MODULATION].to_numpy()
    R2, x_to_plot, y_to_plot = regression_tools.leave_one_out_lr(x, y)
    ax[i].plot(x_to_plot, y_to_plot, color="black", alpha=ALPHA, lw=1)
    # Test data
    corr_test = pearsonr(by_subtype[receptor], by_subtype[MODULATION])
    ax[i].set_title(
        f"$R^2$ = {R2:0.2f} r = {corr_test.statistic:0.2f}", fontsize=FONTSIZE
    )
    print(receptor, corr_test)
    if corr_test.pvalue < 0.05:
        significant_receptors.append(receptor)
    if corr_test.pvalue < 0.05 / len(receptors):
        really_significant_receptors.append(receptor)

sns.despine()
fig.tight_layout()
plt.savefig(snakemake.output.figure, dpi=300)
print(
    f"# Significant: {len(significant_receptors)} / {len(receptors)}, \
    {len(really_significant_receptors)} / {len(receptors)}"
)
print(really_significant_receptors)
