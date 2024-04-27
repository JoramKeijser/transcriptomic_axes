"""
Save gene names from an AnnData file
"""
import numpy as np
import anndata as ad

adata = ad.read_h5ad(snakemake.input.adata)
genes = np.array(np.sort(adata.var_names), dtype=str)
np.savetxt(snakemake.output.gene_list, genes, fmt="%s")
