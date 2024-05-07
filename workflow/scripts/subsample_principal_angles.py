"""
Compute and visualize principal angles
"""
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.linalg import subspace_angles
from src import constants

# Settings for plotting
sns.set_context("poster")
sc.settings.dpi_save = 300
MARKERSIZE = 20
ALPHA = 0.2
LINEWIDTH = 3

print("Loading PCs")
PC_subspace = {}
seeds = []
for dataset in snakemake.input:
    name = dataset.split("/")[-1].split("h5ad")[0].split("_")[0]
    seed = int(dataset.split("/")[-1].split(".h5ad")[0].split("_")[-1])
    if name == snakemake.params.reference:
        seeds.append(seed)
    adata = ad.read_h5ad(dataset)
    if name not in PC_subspace.keys():
        PC_subspace[name] = []
    PC_subspace[name].append(adata.varm["PCs"][:, : constants.NUM_PCS])

for name in PC_subspace.keys():
    PC_subspace[name] = np.array(PC_subspace[name])

print("Plotting principal angles")
reference = snakemake.params.reference
remaining = np.sort([dataset for dataset in PC_subspace.keys() if dataset != reference])


# Add random
n_genes = PC_subspace[reference][0].shape[0]
n_pcs = PC_subspace[reference][seeds[0]].shape[1]
rng = np.random.RandomState(0)
angles = np.zeros((len(seeds), n_pcs))
for i_s, seed in enumerate(seeds):
    idx = rng.permutation(np.arange(n_genes))
    angles[i_s] = (
        subspace_angles(
            PC_subspace[reference][seed][idx], PC_subspace[reference][seed]
        )[::-1]
        * 180
        / np.pi
    )
    plt.scatter(np.arange(n_pcs), angles[i_s], color="gray", s=MARKERSIZE, alpha=ALPHA)
plt.plot(np.arange(n_pcs), angles.mean(0), label="Chance", color="gray", lw=LINEWIDTH)

colors = {
    "bakken": sns.color_palette("Set2")[1],
    "colquitt": sns.color_palette("Set2")[2],
    "tosches": sns.color_palette("Set2")[3],
}
# Compare reference
print("Seeds: ", seeds)
for i, name in enumerate(remaining):
    angles = np.zeros((len(seeds), n_pcs))
    for i_s, seed in enumerate(seeds):
        angles[i_s] = (
            subspace_angles(PC_subspace[name][seed], PC_subspace[reference][seed])[::-1]
            * 180
            / np.pi
        )
        plt.scatter(
            np.arange(n_pcs),
            angles[i_s],
            color=colors[name],
            alpha=ALPHA,
            s=MARKERSIZE,
        )
    # Plot mean
    plt.plot(
        np.arange(n_pcs),
        angles.mean(0),
        label=snakemake.params.areas[name],
        color=colors[name],
        lw=LINEWIDTH,
    )

plt.xlabel("tPC subspace dimension")
plt.ylabel("Principal angle (deg.)")
plt.yticks([50, 70, 90])
plt.legend(handlelength=0.5, framealpha=0.2)
sns.despine()
plt.tight_layout()
plt.savefig(snakemake.output.figure, dpi=300)
