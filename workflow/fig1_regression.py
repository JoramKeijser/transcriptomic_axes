"""
Regression by subclass
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.stats import pearsonr
from src import regression_tools


sc.set_figure_params(
    dpi_save=300,
    figsize=(7, 5),
    fontsize=25,
    frameon=True,
    transparent=False,
    color_map="viridis",
)
plt.rcParams["axes.grid"] = False
MARKERSIZE = 150
subclass_cmap = sns.color_palette("colorblind", 10)
sns.set_palette("colorblind", 10)
sns.set_context("poster")


bugeon = ad.read_h5ad(snakemake.input[0])  # "./results/anndata/bugeon.h5ad"
bugeon.obs["tPC1"] = bugeon.obsm["X_pca"][:, 0]  # Add tPC1 score for easy plotting
print(bugeon)
# Plot all
fig, ax = plt.subplots()
sns.scatterplot(
    bugeon.obs,
    x="tPC1",
    y="State modulation",
    hue="Subclass",
    ax=ax,
    linewidth=0.33,
    alpha=0.67,
)
ax.legend().remove()
sns.despine()
x = np.array(bugeon.obs["tPC1"])[:, None]
y = np.array(bugeon.obs["State modulation"])
ix = np.isnan(y)
print(f"{np.sum(ix)} NaNs")
x, y = x[~ix], y[~ix]
R2, x_to_plot, y_to_plot = regression_tools.leave_one_out_lr(x, y)
ax.plot(x_to_plot, y_to_plot, color="black")
ax.set_title(f"All: $R^2$ = {R2:0.2f}")
ax.set_ylabel("State modulation")
test = pearsonr(x.flatten(), y)
print(f"All: p = {test.pvalue:0.4f}, corr = {test.statistic:0.2f}")
fig.tight_layout()
plt.savefig(snakemake.output[0], dpi=300)

# Plot subclasses
subclasses = bugeon.obs["Subclass"].cat.categories
print(subclasses)
stats = {}
for i, subclass in enumerate(subclasses):
    fig, ax = plt.subplots()
    sns.despine()
    sns.scatterplot(
        bugeon.obs,
        x="tPC1",
        y="State modulation",
        hue="Subclass",
        ax=ax,
        alpha=0.1,
        size=1,
        legend=False,
        linewidth=0.33,
    )
    obs = bugeon[bugeon.obs["Subclass"] == subclass].obs
    sns.scatterplot(
        obs,
        x="tPC1",
        y="State modulation",
        color=sns.color_palette()[i],
        ax=ax,
        alpha=0.67,
        size=1,
        legend=False,
        linewidth=0.33,
    )
    x = np.array(obs["tPC1"])[:, None]
    y = obs["State modulation"]
    ix = np.isnan(y)
    print(f"{np.sum(ix)} NaNs")
    x, y = x[~ix], y[~ix]
    R2, x_to_plot, y_to_plot = regression_tools.leave_one_out_lr(x, y)
    stats[subclass] = (R2, obs.shape[0])  # record R^2 and sample size
    ax.plot(x_to_plot, y_to_plot, color="black")
    ax.set_title(f"{subclass}: $R^2$ = {R2:0.2f}")
    test = pearsonr(x.flatten(), y)
    print(f"{subclass}: p = {test.pvalue:0.4f}, corr = {test.statistic:0.2f}")
    ax.set_ylabel("State modulation")
    fig.tight_layout()
    plt.savefig(snakemake.output[i + 1], dpi=300)

# Now test whether drop in performance is due sample size
# Observation
R2s = [value[0] for value in stats.values()]
ns = [value[1] for value in stats.values()]
plt.figure()
plt.scatter(ns, R2s)
for i, subclass in enumerate(subclasses):
    plt.scatter(ns[i], R2s[i], alpha=0.75, s=200)
    plt.text(ns[i] - 20, R2s[i], subclass)  # annote, center x-coord
sns.despine()
plt.xlabel("# cells")
plt.ylabel("$R^2$")
plt.tight_layout()
plt.savefig(snakemake.output[-2], dpi=300)

subclasses = bugeon.obs["Subclass"].cat.categories
n_subsets = snakemake.params.n_subsets
R2_sub = {}
for i, subclass in enumerate(subclasses):
    obs = bugeon[bugeon.obs["Subclass"] == subclass].obs
    n_subclass = obs.shape[0]
    # Randomly sample from bugeon
    R2_sub[subclass] = np.zeros((n_subsets, 1))
    for s in range(n_subsets):
        bugeon_subsample = sc.pp.subsample(
            bugeon, n_obs=n_subclass, copy=True, random_state=s
        )
        x = np.array(bugeon_subsample.obs["tPC1"])[:, None]
        y = bugeon_subsample.obs["State modulation"]
        ix = np.isnan(y)
        x, y = x[~ix], y[~ix]
        (
            R2_sub[subclass][s],
            x_to_plot,
            y_to_plot,
        ) = regression_tools.leave_one_out_lr(x, y)
    print(
        f"{subclass}: n = {n_subclass}, p = {np.mean(R2_sub[subclass] <= R2s[i]):0.4f}"
    )
    # Add to df
    data = np.concatenate(
        (R2_sub[subclass], np.array([subclass] * n_subsets)[:, None]), axis=1
    )
    if i == 0:
        df = pd.DataFrame(data, columns=["R2", "Subclass"])
    else:
        df2 = pd.DataFrame(data, columns=["R2", "Subclass"])
        df = pd.concat((df, df2), axis=0)
df["R2"] = df["R2"].astype(float)

# Plot this one
fig, ax = plt.subplots()
sns.violinplot(df, x="R2", y="Subclass", inner="box", linewidth=1, ax=ax, alpha=0.5)
ax.set_xlim([-0.2, 0.6])
for i, subclass in enumerate(subclasses):
    ax.vlines(stats[subclass][0], i - 0.33, i + 0.33, color="gray")
    obs = bugeon[bugeon.obs["Subclass"] == subclass].obs
    n_subclass = obs.shape[0]
    ax.text(0.65, i, n_subclass)
ax.text(0.6, -1, "# cells")
sns.despine()
ax.set_xlabel("$R^2$")
ax.set_ylabel("")
plt.tight_layout()
plt.savefig(snakemake.output[-1], dpi=300)
