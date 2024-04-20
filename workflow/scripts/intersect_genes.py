"""
Find one-to-one orthologs that are shared between datasets
"""
import numpy as np

for i, genelist in enumerate(snakemake.input):
    if i == 0:
        shared_genes = np.loadtxt(genelist, dtype=str)
    else:
        new_genes = np.loadtxt(genelist, dtype=str)
        shared_genes = set(shared_genes).intersection(new_genes)

print(f"Found {len(shared_genes)} shared genes")
np.savetxt(snakemake.output.shared_genes, list(shared_genes), fmt="%s")
