# Combine partitions
import numpy as np
import anndata as ad
import pandas as pd
import os
import re
import argparse
from src import data_tools

# TODO: all genes?
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

# Check that we didn't miss any interneurons TODO
