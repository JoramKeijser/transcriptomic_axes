"""
Compute and visualize principal angles
"""
import pickle
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.linalg import subspace_angles
from src import constants

sns.set_context("poster")
sc.settings.dpi_save = 300

print("Loading PCs")
PC_subspace = {}
seeds = []
for dataset in snakemake.input:
    print(dataset)
    name = dataset.split("/")[-1].split("h5ad")[0].split("_")[0]
    seed = int(dataset.split("/")[-1].split(".h5ad")[0].split("_")[-1])
    if name == snakemake.params.reference:
        seeds.append(seed)
    print(dataset, name, seed)
    adata = ad.read_h5ad(dataset)
    if name not in PC_subspace.keys():
        PC_subspace[name] = []
    PC_subspace[name].append(adata.varm["PCs"][:, : constants.NUM_PCS])

for name in PC_subspace.keys():
    PC_subspace[name] = np.array(PC_subspace[name])
    print(PC_subspace[name].shape)

print("Plotting principal angles")
reference = snakemake.params.reference
remaining = np.sort([dataset for dataset in PC_subspace.keys() if dataset != reference])


# Add random
n_genes = PC_subspace[reference][0].shape[0]
rng = np.random.RandomState(0)
for seed in seeds:
    idx = rng.permutation(np.arange(n_genes))
    angles = (
        subspace_angles(
            PC_subspace[reference][seed][idx], PC_subspace[reference][seed]
        )[::-1]
        * 180
        / np.pi
    )
    if seed == min(seeds):  # add label
        plt.plot(angles, label="Chance", color="black", linestyle=":", lw=1, alpha=0.5)
    else:
        plt.plot(angles, color="black", linestyle=":", lw=1, alpha=0.5)

colors = {
    "bakken": sns.color_palette("Set2")[1],
    "colquitt": sns.color_palette("Set2")[2],
    "tosches": sns.color_palette("Set2")[3],
}
# Compare reference
print("Seeds: ", seeds)
for seed in seeds:
    for i, name in enumerate(remaining):
        angles = (
            subspace_angles(PC_subspace[name][seed], PC_subspace[reference][seed])[::-1]
            * 180
            / np.pi
        )
        print(f"Seed {seed}: {angles[0]:0.2f}")
        if seed == min(seeds):  # Label once
            plt.plot(
                angles,
                label=snakemake.params.areas[name],
                color=colors[name],
                lw=1,
                alpha=0.5,
            )
        else:
            print("plot seed")
            plt.plot(angles, color=colors[name], lw=1, alpha=0.5)

plt.xlabel("tPC subspace dimension")
plt.ylabel("Principal angle (deg.)")
plt.yticks([50, 70, 90])
plt.legend(handlelength=0.5)
sns.despine()
plt.tight_layout()
plt.savefig(snakemake.output.figure, dpi=300)
