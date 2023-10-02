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

areas = snakemake.params.areas
print("Loading PCs")
PC_subspace = {}
for dataset in snakemake.input.same_species:
    name = dataset.split("/")[-1].split("h5ad")[0].split("_")[0]
    print(name)
    adata = ad.read_h5ad(dataset)
    PC_subspace[name] = adata.varm["PCs"][:, : constants.NUM_PCS]

print("Plotting principal angles")
reference = snakemake.params.reference
remaining = np.sort([dataset for dataset in areas.keys() if dataset != reference])

# Add random
n_genes = PC_subspace[reference].shape[0]
rng = np.random.RandomState(0)
idx = rng.permutation(np.arange(n_genes))
plt.plot(
    subspace_angles(PC_subspace[reference][idx], PC_subspace[reference])[::-1]
    * 180
    / np.pi,
    label="Chance",
    color="black",
    linestyle=":",
)

colors = {
    "yao": sns.color_palette("Set2")[1],
    "bugeon": sns.color_palette("Set2")[2],
    "hodge": sns.color_palette("Set2")[1],
}
for i, name in enumerate(remaining):
    plt.plot(
        subspace_angles(PC_subspace[name], PC_subspace[reference])[::-1] * 180 / np.pi,
        label=areas[name],
        color=colors[name],
    )

# Add between-species baseline
if snakemake.params.control != "72g": 
    with open(snakemake.input.baseline_angles, 'rb') as handle:
        d = pickle.load(handle)
    for i, angles in enumerate(d.values()):
        if i == 0:
            plt.plot(angles, color="gray", label = "Other species")
        else:
            plt.plot(angles, color="gray")

plt.xlabel("tPC subspace dimension")
plt.ylabel("Principal angle (deg.)")
plt.yticks([0, 45, 90])
plt.legend(handlelength=0.5)
sns.despine()
plt.tight_layout()
plt.savefig(snakemake.output.figure, dpi=300)

print("Saving data")
angles = {}
for name in areas.keys():
    angles[name] = (
        subspace_angles(PC_subspace[reference], PC_subspace[name])[::-1] * 180 / np.pi
    )

savename = snakemake.output.angles
print("Save as", savename)
with open(savename, "wb") as handle:
    pickle.dump(angles, handle, protocol=pickle.HIGHEST_PROTOCOL)
