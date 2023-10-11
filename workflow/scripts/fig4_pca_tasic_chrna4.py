# Color PCA plot by Chrna4 expression
import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc
import seaborn as sns

sns.set_context("poster")

sc.set_figure_params(dpi_save=300)
fig, ax = plt.subplots(figsize=(4, 3))
adata = ad.read_h5ad(snakemake.input[0])
sc.pl.pca(adata, color="Chrna4", frameon=False, ax=ax)
fig.tight_layout()
plt.savefig(snakemake.output[0])
