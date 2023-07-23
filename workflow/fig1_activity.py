"""
Visualize neural activity
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.gridspec import GridSpec
import seaborn as sns
import pandas as pd
import anndata as ad
from src import data_tools
from sklearn.decomposition import PCA
from scipy.stats import ttest_ind, mannwhitneyu
from scipy.ndimage import convolve1d, gaussian_filter1d
sns.set_palette("colorblind")
sns.set_context("poster")


#save_dir = "./results/pandas/"
basedir = snakemake.input[0] + "/" 
stimulus = "Blank"  # stimulus - only take gray screen
all_adatas = ad.AnnData()
adata_dict = {}
exp_folders = [d for d in os.listdir(basedir) if os.path.isdir(basedir + d)]
#exp_folder = exp_folders[0]
for exp_folder in exp_folders:
    sessions = [d for d in os.listdir(basedir + exp_folder) if os.path.isdir(basedir + exp_folder)]
    for session in sessions:
        session_path = basedir + exp_folder + "/" + session + "/"
        adata = data_tools.load_activity(basedir, session_path, session_path + stimulus)
        adata.var['experiment'] = exp_folder + "/" + session   
        adata_dict[exp_folder + "/" + session] = adata.copy()
        if all_adatas.shape[0] == 0:
            all_adatas = adata.copy()
        else:
            all_adatas = ad.concat([adata, all_adatas], axis=1)

var_df = pd.concat([adata.var for adata in adata_dict.values()])
print("var_df:", var_df.shape)
var_df['uid'] = var_df.index 
var_df = var_df.drop_duplicates('uid', keep = 'first') # Same neuron recorded twice - only keep first recording
var_df.drop("uid", inplace=True, axis=1) 
var_df.index = var_df.index.astype("str")
# Drop cells that belong to t types with too few counts
few_counts = var_df.value_counts("Subtype")[np.where(var_df.value_counts("Subtype") < 3 )[0]].index
var_df = var_df[~var_df['Subtype'].isin(few_counts)]
print("var_df:", var_df.shape)
print(var_df.head())

print(var_df["State modulation"])
nan_cell = var_df["State modulation"].isna()
print(f'{np.sum(nan_cell)} cells with missing State modulation, {np.sum(~nan_cell)} left')

for variable in ["Running", "Stationary Desynchronized", "State modulation"]:
    var_df[variable] = var_df[variable].astype(float) 
var_df.to_csv(snakemake.output[0]) #f"{save_dir}bugeon_activity.h5ad"

    # Visualize 
fontsize = 20
modulation = 'State modulation'
fig, ax = plt.subplots(figsize=(5,4))
sns.violinplot(x = 'Subclass', y = modulation, data=var_df, hue='Subclass', legend=None, inner=None, 
linewidth=.0, dodge=0, saturation=.9, ax=ax)
ax.legend().remove()
sns.swarmplot(x = 'Subclass', y = modulation, data=var_df,  legend=False,
    size=1., color = 'black', alpha=0.67)
sns.despine()
ax.hlines(0, -.5, 4.5, linestyles=":", color='black', lw=1)
ax.set_ylim([-1.2, 1.2]) 
ax.set_yticks([-1, 0,  1])
ax.set_xlabel("")
ax.set_ylabel(modulation, fontsize=fontsize)
ax.set_xticklabels(var_df['Subclass'].cat.categories, rotation=45, fontsize=fontsize)
plt.tight_layout()
print(f"{modulation}: n = {var_df[~var_df[modulation].isna()].shape[0]} interneurons:")
print(var_df[~var_df[modulation].isna()].value_counts("Subclass"))                    
plt.tight_layout()
plt.savefig(snakemake.output[1], dpi=300)


