"""
Test all receptors
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import anndata as ad
import scanpy as sc
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
sc.settings.figdir = snakemake.params.figdir
ALPHA = 0.85
FONTSIZE = 50
subclass_cmap = sns.color_palette("colorblind", 10)
sns.set_palette("colorblind", 10)
# sns.set_context("poster")
mpl.rcParams.update({"font.size": 7})


# Load data - avg by subtype because we don't have ACh expression and state mod for the same cells
subclass_order = ["Pvalb", "Sst", "Lamp5", "Vip", "Sncg"]
by_subtype = pd.read_csv(snakemake.input.by_subtype)
by_subtype["Subclass"] = by_subtype["Subclass"].astype("category")
by_subtype["Subclass"] = by_subtype["Subclass"].cat.reorder_categories(subclass_order)
by_subtype.index = by_subtype["Subtype"]
# Load Tasic dataset for ACh receptor expression
tasic = ad.read_h5ad(snakemake.input.tasic)
# Consistent coding of subtype names: dash instead of space
# tasic.obs['Subtype'] = [subtype.replace(" ", "-") for subtype in tasic.obs['subtype']]
# tasic.obs['Subtype'] = tasic.obs['Subtype'].astype("category")
# Select bugeon subtypes
tasic = tasic[tasic.obs["Subtype"].isin(by_subtype["Subtype"])]
tasic.obs["Subclass"] = tasic.obs["Subclass"].cat.reorder_categories(subclass_order)
obs_df = tasic.obs["Subtype"].to_frame()
sc.pp.normalize_total(tasic, target_sum=1e4)

# Sort by expression
all_receptors = [
    gene
    for gene in tasic.var_names
    if gene.startswith("Chrm") or gene.startswith("Chrn")
]
expression = np.array([tasic[:, receptor].X.sum() for receptor in all_receptors])
# Threshold?
# now sort by correlation
idx = np.argsort(expression)[::-1]
all_receptors = np.array(all_receptors)[idx]
expression = expression[idx]

savename = snakemake.output.tracksplot.split(snakemake.params.figdir + "tracksplot")[1]
sc.pl.tracksplot(
    tasic,
    all_receptors,
    groupby="Subclass",
    dendrogram=False,
    logTODO=False,
    save=savename,
    show=False,
)
# Combine at subtype level
for receptor in all_receptors:
    obs_df[receptor] = np.log1p(np.array(tasic[:, receptor].X[:, 0]))
tasic_by_subtype = obs_df.groupby("Subtype").mean()
by_subtype = pd.concat(
    [by_subtype, tasic_by_subtype[all_receptors]], axis=1, ignore_index=False
)

print("Computing significance")
MODULATION = "State modulation"
significant_receptors = []
R2_perm = {}
R2 = {}
for i, receptor in enumerate(all_receptors):
    R2_perm[receptor] = np.zeros((snakemake.params.n_permutations,))
    x = by_subtype[receptor].to_numpy()[:, None]
    y = by_subtype[MODULATION].to_numpy()
    R2[receptor], _, _ = regression_tools.leave_one_out_lr(x, y)
    for p in range(snakemake.params.n_permutations):
        R2_perm[receptor][p], _, _ = regression_tools.leave_one_out_lr(
            x, np.random.permutation(y)
        )
    pvalue = np.mean(R2_perm[receptor] >= R2[receptor])
    expr = tasic[:, receptor].X.sum()
    if pvalue < 0.05:
        significant_receptors.append(receptor)
        print(f"{receptor}: p = {pvalue:0.4f}*, expr = {expr}")
    else:
        print(f"{receptor}: p = {pvalue:0.4f} (n.s.), expr = {expr}")
print(significant_receptors)
significant_receptors = np.array(significant_receptors)

# Make tracksplot of significant receptors
savename = snakemake.output.tracksplot_significant.split(
    snakemake.params.figdir + "tracksplot"
)[1]
significant_expression = np.array(
    [tasic[:, receptor].X.sum() for receptor in significant_receptors]
)
significant_receptors = significant_receptors[np.argsort(significant_expression)[::-1]]
not_significant = np.array(
    [receptor for receptor in all_receptors if receptor not in significant_receptors]
)
not_significant_expression = np.array(
    [tasic[:, receptor].X.sum() for receptor in not_significant]
)
not_significant = not_significant[np.argsort(not_significant_expression)[::-1]]
all_receptors = np.concatenate((significant_receptors, not_significant))
all_expression = np.concatenate((significant_expression, not_significant_expression))
# Plot those with expr >= 1 CP10K
sc.pl.tracksplot(
    tasic,
    all_receptors[all_expression >= 1.0],
    groupby="Subclass",
    dendrogram=False,
    log=False,
    save=savename,
    show=False,
)

sns.set_context("poster")
FONTSIZE = 20
# Next: plot examples of state modulation. Equal number predictive and not predictive
n_significant = len(significant_receptors)
for receptor_list, figname in zip(
    [significant_receptors, not_significant[:n_significant]],
    [snakemake.output.hist_sig, snakemake.output.hist_notsig],
):
    fig, ax = plt.subplots(n_significant, 2, figsize=(10, 15))
    for row, receptor in enumerate(receptor_list):
        x = by_subtype[receptor].to_numpy()[:, None]
        y = by_subtype[MODULATION].to_numpy()
        _, x_to_plot, y_to_plot = regression_tools.leave_one_out_lr(x, y)
        # Show fit (train data)
        ax[row, 0].plot(x_to_plot, y_to_plot, color="black", alpha=ALPHA, lw=3)
        # Expression vs modulation
        sns.scatterplot(
            x=receptor,
            y=MODULATION,
            data=by_subtype,
            hue="Subclass",
            ax=ax[row, 0],
            legend=None,
            alpha=ALPHA,
            s=100,
        )
        if row == 0:
            ax[row, 0].set_ylabel(MODULATION, fontsize=FONTSIZE)
            ax[row, 1].set_ylabel("# shuffles", fontsize=FONTSIZE)
            ax[row, 1].set_xlabel("$R^2$", fontsize=FONTSIZE)

            ax[row, 0].set_xlabel("Expr. (log CP10K)", fontsize=FONTSIZE)
        else:
            ax[row, 0].set_xlabel("")
            ax[row, 0].set_ylabel("")
            ax[row, 1].set_xlabel("")
            ax[row, 1].set_ylabel("")

        # Hist of R2s
        ax[row, 1].hist(
            R2_perm[receptor],
            bins=int(0.1 * snakemake.params.n_permutations),
            color="gray",
        )
        ax[row, 1].vlines(
            R2[receptor],
            0,
            0.1 * snakemake.params.n_permutations,
            linestyles=":",
            color="tab:red",
            lw=3,
        )
        ax[row, 0].set_title(receptor)
        pvalue = np.mean(R2_perm[receptor] >= R2[receptor])
        ax[row, 1].set_title(f"$R^2$ = {R2[receptor]:0.2f}")
    sns.despine()
    fig.tight_layout()
    plt.savefig(figname, dpi=300)
