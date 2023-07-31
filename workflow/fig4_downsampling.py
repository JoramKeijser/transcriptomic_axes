## Test effect of varying sequencing depth
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import argparse
import pickle
import anndata as ad
import scanpy as sc

sc.set_figure_params(dpi_save=300, figsize=(7,5), fontsize=25, 
    frameon=True, transparent=False, color_map="viridis")
plt.rcParams["axes.grid"] = False
marker_size = 150
alpha = 0.67
fontsize = 25
subclass_cmap = sns.color_palette("colorblind", 10)
sns.set_palette("colorblind", 10)
sns.set_context("poster")

rng = np.random.RandomState(123)


datasets = {}
datasets['reference'] = ad.read_h5ad(snakemake.input.dataset)
datasets['shallow'] = ad.read_h5ad(snakemake.input.shallow)
# Which receptors to compare?
#receptors = np.loadtxt(snakemake.input.receptors, dtype=str)
#shared_genes = np.loadtxt(snakemake.input.shared_genes, dtype=str)
receptors = np.loadtxt(snakemake.input.receptors, dtype=str)
shared_genes = np.loadtxt(snakemake.input.shared_genes, dtype=str)
genes = datasets['reference'].var_names # TODO 
genes = set(datasets['reference'].var_names).intersection(datasets['shallow'].var_names)
shared_genes = np.sort(list(set(receptors).intersection(genes)))
# Subset to shared genes - compute depth
for dataset in ['reference', 'shallow']:
    datasets[dataset] = datasets[dataset][:, shared_genes].copy()
    datasets[dataset].obs['total_counts'] = datasets[dataset].X.sum(1)
genes = np.sort(list(set(receptors).intersection(shared_genes)))

depth = np.median(datasets['shallow'].obs['total_counts']) 
relative_depth = float(depth / np.median(datasets['reference'].obs['total_counts']) )
print(f"relative depth: {relative_depth:0.3f}")
assert relative_depth <= 1.0
datasets['subsampled'] = datasets['reference'].copy() 
# Make sure shallower is not larger
if datasets['shallow'].shape[0] > datasets['reference'].shape[0]:
    print("Subsampling shallower dataset")
    _, idx = sc.pp.subsample(datasets['shallow'].X, n_obs = datasets['reference'].shape[0], copy=True)
    datasets['shallow'] = datasets['shallow'][idx]

all_samples = {gene: [] for gene in genes}
control = {gene: datasets['shallow'][:, gene].X.sum() for gene in genes}
subsampled_counts = {gene: np.zeros((snakemake.params.num_samples, )) for gene in genes}
for sample in range(snakemake.params.num_samples):
    _, idx = sc.pp.subsample(datasets['reference'].X, n_obs = datasets['shallow'].shape[0], copy=True)
    datasets['subsampled'] = datasets['reference'][idx].copy()
    datasets['subsampled'].X = rng.binomial(np.array(datasets['subsampled'].X, dtype=int), p = relative_depth)
    for gene in genes:
        subsampled_counts[gene][sample] = datasets['subsampled'][:, gene].X.sum()

# Put result into data frame
df = pd.DataFrame(columns = ['statistic', 'gene'])
pseudo_count = 1.0
for gene in subsampled_counts.keys():
    statistic = np.log2(float(control[gene] + pseudo_count) / (subsampled_counts[gene]  + pseudo_count))
    gene_label = np.array([gene] * len(statistic))[:,None]
    name_label = np.array([snakemake.params.dataset] * len(statistic))[:,None]
    data = np.concatenate((statistic[:,None], gene_label, name_label), -1)
    df = pd.concat((df, pd.DataFrame(data, columns = ['statistic', 'gene', 'name'])))
df['statistic'] = df['statistic'].astype("float")

fig, ax = plt.subplots(figsize=(6,5))
sns.violinplot(df, x='statistic', y='gene', inner='box', 
    color = sns.color_palette("colorblind")[2],alpha=0.2,
    saturation=1, linewidth=1, ax=ax, scale='width')
plt.xlabel("$\Delta$ expression (log2fc)")
plt.ylabel("")
plt.vlines(0, -.5, 4.5, color='gray', linestyle=":", lw=1)
plt.fill_betweenx([-.5, 4.5], -2, 2, color='gray', alpha=.1)
sns.despine()
plt.title(snakemake.params.species)
plt.tight_layout()
plt.savefig(snakemake.output.vln, dpi=300)
#plt.xlim([-7, 7])
#plt.tight_layout()
#plt.savefig(f"./figures/figure5/subsample_vln_{args.shallow_data}_{args.num_samples}_sharedlim.png", dpi=300)

#df.to_csv(savedir + f"{args.shallow_data}_{args.num_samples}.csv")

