"""
Mouse variance explained by other datasets
"""
import numpy as np
import anndata as ad
import scanpy as sc
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from src.pca_tools import compare_variance
from src import constants

sns.set_palette("Set2")
sns.set_context("poster")
sc.settings.dpi_save = 300


colors = {
    "tasic": sns.color_palette("Set2")[0],
    "bakken": sns.color_palette("Set2")[1],
    "colquitt": sns.color_palette("Set2")[2],
    "tosches": sns.color_palette("Set2")[3],
}
species = {
    "colquitt": "Zebra finch",
    "tosches": "Turtle",
    "tasic": "Mouse",
    "bakken": "Human",
}
# Load data
reference = snakemake.params.reference
datasets = {}
seeds = []
print("Loading data")
for dataset in snakemake.input:
    name = dataset.split("/")[-1].split("h5ad")[0].split("_")[0]
    seed = int(dataset.split("/")[-1].split(".h5ad")[0].split("_")[-1])
    if name == snakemake.params.reference:
        seeds.append(seed)
    if name not in datasets.keys():
        datasets[name] = {}
    datasets[name][seed] = ad.read_h5ad(dataset)

print(seeds)
cross_variance = {}  # dist of lists, one for each dataset
variance = []  # one for each seed
print("Compute cross-covariance")
for seed in seeds:
    # PCA on reference
    hvgs = datasets[reference][seed].var.highly_variable
    datasets[reference][seed] = datasets[reference][seed][:, hvgs]
    pca = PCA(n_components=constants.NUM_PCS).fit(datasets[reference][seed].X)
    C = pca.get_covariance()
    w = pca.components_[0]
    # Normalized variance for reference
    unnormed_var = np.diag(
        pca.components_[:10] @ C @ pca.components_[:10].T
    ) / np.trace(C)
    variance.append(unnormed_var / unnormed_var[0])

    # Another pass, now across other datasets
    seed_datasets = [key for key in datasets.keys()]
    for name in seed_datasets:
        dataset = datasets[name][seed]
        if name == reference:
            continue
        # Focus on HVGs of reference
        datasets[name][seed] = datasets[name][seed][:, hvgs]
        if name not in cross_variance.keys():
            cross_variance[name] = []
        cross_variance[name].append(
            compare_variance(datasets[reference][seed], datasets[name][seed], pca)
        )

# Nor plt them all
x = np.arange(5) * 4
plt.scatter(
    x, np.array(variance)[:, :5].mean(0), color=colors[reference], marker="d", s=100
)
for i_s, seed in enumerate(seeds):
    plt.scatter(
        x, np.array(variance)[:, :5].mean(0), color=colors[reference], s=20, alpha=0.2
    )


for i, name in enumerate(["bakken", "colquitt", "tosches"]):
    plt.scatter(
        x + 0.8 * (i + 1),
        np.array(cross_variance[name])[:, :5].mean(0),
        color=colors[name],
        marker="d",
        s=100,
    )

    for i_s, seed in enumerate(seeds):
        jitter = np.random.uniform(-0.3, 0.3)
        plt.scatter(
            x + 0.8 * (i + 1) + jitter,
            np.array(cross_variance[name])[i_s, :5],
            color=colors[name],
            s=20,
            alpha=0.2,
        )

plt.xticks(x + 1.0, np.arange(5) + 1)
plt.yticks([0, 0.5, 1])
sns.despine()
plt.ylabel("Variance (norm.)")
plt.xlabel("tPC")
plt.ylim([-0.04, 1.1])
plt.tight_layout()
# Save figure
plt.savefig(snakemake.output.figure, dpi=300)
