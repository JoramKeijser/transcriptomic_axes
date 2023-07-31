# Find significant receptors
import numpy as np
import anndata as ad
import pandas as pd
import os
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
sns.set_palette("colorblind")
sns.set_context("poster")


rng = np.random.RandomState(123)
# Load data
bugeon = ad.read_h5ad(snakemake.input.bugeon)
tasic = ad.read_h5ad(snakemake.input.tasic)
# Find ACh receptors, create df with their subtype-averaged expression
receptors = np.sort([var for var in tasic.var_names if var.startswith("Chrm")  or var.startswith("Chrna")])
tasic_obs = tasic.obs[['Subtype']]
tasic_obs[receptors] = np.array(tasic[:, receptors].X)
tasic_subtypes =  tasic_obs.groupby("Subtype").mean()
# Combine with subtype-averaged state modulation
bugeon_subtypes =  bugeon.obs[['state_modulation', 'Subtype']].groupby("Subtype").mean()
bugeon_subtypes.loc[:, receptors] = tasic_subtypes.loc[bugeon_subtypes.index, receptors]

# Correlate
def get_corr(shuffle = False):
    bugeon_subtypes = bugeon.obs[['state_modulation', 'Subtype']].groupby("Subtype").mean()
    n_subtypes = bugeon_subtypes.shape[0]
    bugeon_subtypes.loc[:, receptors] = tasic_subtypes.loc[bugeon_subtypes.index, receptors]
    if shuffle:
        idx = np.random.choice(n_subtypes, n_subtypes)
        bugeon_subtypes.loc[:, receptors] = bugeon_subtypes.iloc[idx].loc[:, receptors].to_numpy()
    C = bugeon_subtypes.corr()['state_modulation']
    C.drop("state_modulation", inplace=True)
    return C


C = get_corr(False)
random_corrs = pd.DataFrame()
# Take care of occasional NaNs
rep = 0
while rep < snakemake.params.permutations:
    new_df = get_corr(True)
    if np.sum(new_df.isna()) == 0:
        random_corrs = pd.concat((random_corrs, new_df), axis=1)
        rep += 1

# Plot
print(random_corrs.shape)
plt.figure(figsize=(7,5))
idx = np.argsort(C)
significant_receptors = []
short_names = [receptor.split("Chr")[1] for receptor in C[idx].index]
plt.xticks(np.arange(len(idx)), short_names, rotation=90)
plt.ylim([-.7, .7])
plt.hlines(0, 0, len(idx)-1, linestyles=':', color='black')
print("receptor, p-value:")
for i, receptor in enumerate(receptors):
    if receptor in ['Chrm2', 'Chrm4']:
        color = 'tab:blue'
    else:
        color = 'tab:red'
    x = np.where(C[idx].index == receptor)[0][0]
    plt.scatter(x, C.loc[receptor], color=color, marker='o', s=40)
    p = np.mean(random_corrs.loc[receptor] > np.abs(C.loc[receptor]))
    if p < 0.05:
        plt.text(x, .65, "*", color='black')
        significant_receptors.append(receptor)
        print(receptor, p)
sns.despine()
m = random_corrs.iloc[idx].mean(1)
sd = random_corrs.iloc[idx].std(1) 
plt.fill_between(np.arange(len(idx)), m - sd, m + sd, color='gray', alpha=0.33)
plt.yticks([-0.6, -0.3, 0, 0.3, 0.6]);
plt.ylabel("corr. w/ \n state modulation")
plt.tight_layout()
plt.savefig(snakemake.output.figure, dpi=300)

print("Significant receptors:")
print(significant_receptors)
# To do: intersect with shared genes
# if not (os.path.exists(fname) and os.path.exists(fname)):
    
#     np.savetxt(results_dir + "corrs.txt", C)
#     np.savetxt(results_dir + f"corrs_shuffle_{args.permutations}.txt", random_corrs)

np.savetxt(snakemake.output.receptors, significant_receptors, fmt="%s")
