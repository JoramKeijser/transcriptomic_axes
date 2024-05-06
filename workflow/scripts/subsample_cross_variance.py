"""
Mouse variance explained by other datasets
"""
import pickle
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
sc.settings.figdir = "./figures/figure2/"
sc.settings.dpi_save = 300


def main():
    plt.figure()
    colors = {
            "tasic": sns.color_palette('Set2')[0],
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
            datasets[name] = []
        datasets[name].append(ad.read_h5ad(dataset))

    # PCA on reference
    for seed in seeds:
        cross_variance = {}
        hvgs = datasets[reference][seed].var.highly_variable
        datasets[reference][seed] = datasets[reference][seed][:, hvgs]
        pca = PCA(n_components=constants.NUM_PCS).fit(datasets[reference][seed].X)
        C = pca.get_covariance()
        w = pca.components_[0]
        variance = w @ C @ w / np.trace(C)

        # Another pass
        print("Compute cross-covariance")
        seed_datasets = [key for key in datasets.keys()]
        for name in seed_datasets:
            dataset = datasets[name][seed]
            if name == reference:
                continue
            # Focus on HVGs of reference
            datasets[name][seed] = datasets[name][seed][:, hvgs]
            cross_variance[name] = compare_variance(
                datasets[reference][seed], datasets[name][seed], pca
            )

        variance = np.diag(pca.components_[:10] @ C @ pca.components_[:10].T) / np.trace(C)
        variance /= variance[0]

        # Plot
        x = np.arange(5) * 4
        # Chance level: random components
        rng = np.random.RandomState(0)
        random_components = rng.normal(size=(10, C.shape[0]))
        random_components /= np.linalg.norm(random_components, axis=1, keepdims=True)
        random_variance = np.diag(
            random_components[:10] @ C @ random_components[:10].T
        ) / float(pca.components_[0] @ C @ pca.components_[0])
        print(f"random tPC1: {random_variance[0] * 100:0.1f}% of mouse tPC1")
        print(f"random tPC1-10: {np.sum(random_variance) * 100:0.1f}% of mouse tPC1")
        plt.hlines(
            random_variance[0], -1, 20, color="black", linestyles=":", lw=3,
        )
        print(f"tasic tPC1-10: {variance[:10].sum()*100:0.1f}% of mouse tPC1")
        jitter = np.random.uniform(-0.3, 0.3)
        plt.scatter(x + jitter, variance[:5], s=30, alpha=0.3, color=colors[reference])

        for i, name in enumerate(["bakken", "colquitt", "tosches"]):
            jitter = np.random.uniform(-0.3, 0.3)
            plt.scatter(x + 0.8 * (i + 1) + jitter, cross_variance[name][:5], s=30, alpha=0.3,
                        color=colors[name])
            print(
                f"{name} tPC1-10: {cross_variance[name][:10].sum()*100:0.1f}% of mouse tPC1"
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


if __name__ == "__main__":
    main()
