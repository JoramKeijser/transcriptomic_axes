import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.linalg import subspace_angles
import argparse
from src import constants
sns.set_context("poster")
sc.settings.dpi_save= 300

areas = snakemake.params.areas

print("Loading PCs")
PC_subspace = {}
for dataset in snakemake.input:
    name = dataset.split("/")[-1].split("h5ad")[0].split("_")[0]
    print(name)
    adata = ad.read_h5ad(dataset)
    PC_subspace[name] = adata.varm['PCs'][:, :constants.NUM_PCS]

print("Plotting principal angles")
reference = 'tasic'
remaining = np.sort([dataset for dataset in areas.keys() if dataset != reference])

# Add random
n_genes = PC_subspace[reference].shape[0]
rng = np.random.RandomState(0)
idx = rng.permutation(np.arange(n_genes))
label = 'Chance'
plt.plot(subspace_angles(PC_subspace[reference][idx], PC_subspace[reference])[::-1] * 180 / np.pi, 
            label = label, color='black', linestyle = ":")
        
colors = {'yao': sns.color_palette("Set2")[1], 'bugeon': sns.color_palette("Set2")[2]}
for i, name in enumerate(remaining):
    plt.plot(subspace_angles(PC_subspace[name], PC_subspace[reference])[::-1] * 180 / np.pi, 
            label = areas[name], color=colors[name])

plt.xlabel("tPC subspace dimension")
plt.ylabel("Principal angle (deg.)")
plt.yticks([45, 90])
plt.legend(handlelength=0.5)
sns.despine()
plt.tight_layout()
plt.savefig(snakemake.output.figure, dpi=300)

print("Saving data")
angles = {}
for name in areas.keys():
    angles[name] = subspace_angles(PC_subspace['tasic'], PC_subspace[name])[::-1] * 180 / np.pi

savename = snakemake.output.angles
print("Save as", savename)
with open(savename, 'wb') as handle:
    pickle.dump(angles, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
