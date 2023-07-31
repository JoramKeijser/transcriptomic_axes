import pandas as pd
import anndata as ad
import numpy as np
import h5py


START = snakemake.params.start_row
END = START + snakemake.params.num_rows
# Read & rename all cell info aka meta data
cell_info = pd.read_csv(snakemake.input.metadata, index_col=0, header=0) #(1169213, 56)
n_cells = cell_info.shape[0]
cell_info = cell_info.rename(columns={'subclass_label':'Subclass', 'class_label': 'Class'})
n_cells = cell_info.shape[0]
# Subset gaba
cell_info = cell_info[cell_info['Class'] == "GABAergic"]
gaba_ids = np.array(cell_info.index).flatten()
# How many GABA cells?
n_gaba = np.sum(cell_info.shape[0])
cell_info['Subclass'] = cell_info['Subclass'].replace({'Serpinf1': 'Vip', 'Sst Chodl': 'Sst'})

# Load just the gene names counts
f = h5py.File(snakemake.input.counts, "r")
gene_names = np.array(f['data']['gene'], dtype=str)
sample_names = np.array(f['data']['samples'], dtype=str)
# Load actual counts
f = h5py.File(snakemake.input.counts, "r")
counts = f['data']['counts'][:, START:END]
counts = counts.T
counts = pd.DataFrame(counts, index = sample_names[START:END], columns = gene_names)
# Subset GABA
shared_ids = list(set(gaba_ids).intersection(counts.index))
counts = counts.loc[shared_ids]
cell_info = cell_info.loc[shared_ids]

# Create adata with current sub
adata = ad.AnnData(X = counts, dtype=int, obs = cell_info)
# To do: order
subclass_order = ['Pvalb', 'Sst', 'Lamp5', 'Vip', 'Sncg', 'Meis2']
adata.obs['Subclass'] = adata.obs['Subclass'].astype("category")

missing_subclasses = [subclass for subclass in subclass_order if subclass not in adata.obs['Subclass'].cat.categories]
adata.obs['Subclass'] = adata.obs['Subclass'].cat.add_categories(missing_subclasses)
adata.obs['Subclass'] = adata.obs['Subclass'].cat.reorder_categories(subclass_order)

adata.write_h5ad(snakemake.output.anndata)

