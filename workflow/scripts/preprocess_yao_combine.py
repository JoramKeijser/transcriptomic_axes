"""
Contruct all GABAergic cells from the yao datasets
using the indidual subsets (see preprocess_yao_partition.py)
"""
import anndata as ad
from src import data_tools

files = snakemake.input.files
adata = ad.read_h5ad(files.pop(), backed="r")
for i, file in enumerate(files):
    print(f"{i+1}/{len(files)}")
    next_adata = ad.read_h5ad(file, backed="r")
    # We know next_data has same genes as all others
    adata = ad.concat((adata, next_adata))
    del next_adata

# Add missing subclasses
adata = data_tools.organize_subclass_labels(adata)
print("Final:", adata.shape)
adata.write_h5ad(snakemake.output.anndata)
