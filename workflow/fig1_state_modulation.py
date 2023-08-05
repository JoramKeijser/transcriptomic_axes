"""
Compute state modulation
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import anndata as ad
from scipy.stats import mannwhitneyu
from src import data_tools

sns.set_palette("colorblind")
sns.set_context("poster")


basedir = snakemake.input.data_dir + "/"
stimulus = snakemake.params.stimulus
all_adatas = ad.AnnData()
adata_dict = {}
exp_folders = [d for d in os.listdir(basedir) if os.path.isdir(basedir + d)]
# Grab all experiments with the right stiulus
for exp_folder in exp_folders:
    sessions = [d for d in os.listdir(basedir + exp_folder) if os.path.isdir(basedir + exp_folder)]
    for session in sessions:
        session_path = basedir + exp_folder + "/" + session + "/"
        adata = data_tools.load_activity(basedir, session_path, session_path + stimulus)
        adata.var['experiment'] = exp_folder + "/" + session
        adata_dict[exp_folder + "/" + session] = adata.copy()
        if all_adatas.shape[0] == 0: # No activity to add
            all_adatas = adata.copy()
        else:
            all_adatas = ad.concat([adata, all_adatas], axis=1)

var_df = pd.concat([adata.var for adata in adata_dict.values()])
var_df['uid'] = var_df.index
# Same neuron recorded twice - only keep first recording:
var_df = var_df.drop_duplicates('uid', keep = 'first')
var_df.drop("uid", inplace=True, axis=1)
var_df.index = var_df.index.astype("str")
# Drop cells that belong to t types with too few counts
few_counts = var_df.value_counts("Subtype")[np.where(var_df.value_counts("Subtype") < 3 )[0]].index
var_df = var_df[~var_df['Subtype'].isin(few_counts)]

nan_cell = var_df["State modulation"].isna()
print(f'{np.sum(nan_cell)} cells with missing State modulation, {np.sum(~nan_cell)} left')

for variable in ["Running", "Stationary Desynchronized", "State modulation"]:
    var_df[variable] = var_df[variable].astype(float)
var_df.to_csv(snakemake.output.data)

# Visualize
FONTSIZE = 20
MODULATION = 'State modulation'
fig, ax = plt.subplots(figsize=(5,4))
sns.violinplot(x = 'Subclass', y = MODULATION, data=var_df, hue='Subclass', legend=None, inner=None,
linewidth=.0, dodge=0, saturation=.9, ax=ax)
ax.legend().remove()
sns.swarmplot(x = 'Subclass', y = MODULATION, data=var_df,  legend=False,
    size=1., color = 'black', alpha=0.67)
sns.despine()
ax.hlines(0, -.5, 4.5, linestyles=":", color='black', lw=1)
ax.set_ylim([-1.2, 1.2])
ax.set_yticks([-1, 0,  1])
ax.set_xlabel("")
ax.set_ylabel(MODULATION, fontsize=FONTSIZE)
ax.set_xticklabels(var_df['Subclass'].cat.categories, rotation=45, fontsize=FONTSIZE)
plt.tight_layout()
print(f"{MODULATION}: n = {var_df[~var_df[MODULATION].isna()].shape[0]} interneurons:")
print(var_df[~var_df[MODULATION].isna()].value_counts("Subclass"))
print(var_df.columns)
print(var_df[[MODULATION, "Subclass"]].groupby("Subclass").mean())
# Add significance
# Add pairwise significance
var_df.dropna(inplace=True, subset="State modulation")
y, h = 1.1, .05
for i, subclass in enumerate(var_df['Subclass'].cat.categories):
    mod1 = np.array(var_df[var_df['Subclass'] == subclass][MODULATION]).flatten()
    for j, subclass2 in enumerate(var_df['Subclass'].cat.categories):
        mod2 = np.array(var_df[var_df['Subclass'] == subclass2][MODULATION]).flatten()
        if j > i:
            pval = mannwhitneyu(mod1, mod2).pvalue
            if j == i + 1:
                # Add it
                if pval < 0.05:
                    plt.plot([i, i, j, j], [y, y+h, y+h, y], lw=.8, color='black')
                    plt.text((i+j)*.5, y+h, "*",
                             ha='center', va='bottom', color='black', fontsize=FONTSIZE)
plt.tight_layout()
plt.savefig(snakemake.output.figure, dpi=300)
