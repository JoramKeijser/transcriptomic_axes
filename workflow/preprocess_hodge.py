# Put raw hodge data into AnnData
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from src import data_tools 
from src import constants, pca_tools
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_palette("colorblind")
sns.set_context("poster")
sc.settings.figdir= "./figures/figure4/"
sc.settings.dpi_save = 300

def main():
   
    metadata = pd.read_csv(snakemake.input.cells, index_col=0)
    introns = pd.read_csv(snakemake.input.introns, 
                      index_col=0, header=0)
    exons = pd.read_csv(snakemake.input.exons, 
                      index_col=0, header=0)
    counts = introns + exons
    del introns, exons
    genes = pd.read_csv(snakemake.input.genes)
    genes = [gene[0] + gene[1:].lower() for gene in genes['gene']]

    # Put into AnnData, subset GABAergic
    adata = ad.AnnData(X = np.array(counts).T, obs=metadata)
    adata.var_names = genes
    depth = np.median(adata.X.sum(1)) # reads / cell
    num_genes = np.median((adata.X>0).sum(1)) # genes / cell
    print(f"# counts: {depth}, # genes: {num_genes}")

    adata = adata[adata.obs['class'] == "GABAergic"]
    #adata.obs['Subclass'] = [cluster.split(" ")[2] for cluster in adata.obs['cluster']]
    #adata.obs['Subclass'] = [subclass[0] + subclass[1:].lower() for subclass in adata.obs['Subclass']]  
    adata = data_tools.organize_subclass_labels(adata)
    adata.write_h5ad(snakemake.output.anndata)

    """
    print("Saving data as", savedir + "hodge.h5ad")
    adata.write_h5ad(savedir + "hodge.h5ad")

    print("PCA")
    shared_genes = np.loadtxt(snakemake.input.shared_genes, dtype = str)
    adata = adata[:, shared_genes]
    sc.pp.normalize_total(adata, target_sum=constants.NORMALIZE_TARGET_SUM)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes = constants.NUM_HVG_GENES)
    sc.pp.pca(adata, constants.NUM_PCS)
    # Consistent (but arbitrary) orientation along PC1,2 for visual comparison
    adata = pca_tools.orient_axes(adata)
    # Assign clusters to the right subclass based on PCA
    print("Clustering")
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=.1)
    sc.pl.pca(adata, color='leiden', legend_loc='on data', save="_hodge_leiden.png")
    equiv = {"0":"Pvalb", "8":"Pvalb", 
         "2":"Sst", "7":"Sst", 
         "3": "Lamp5", "4": "Lamp5", "6":"Lamp5", 
         "1": "Vip", "5": "Vip"}
    sc.pl.pca(adata, color='leiden', legend_loc='on data', save="_hodge_subclass.png")

    # Plot it
    fig, ax = plt.subplots()
    sns.despine()
    sc.pl.pca(adata, color='Subclass', ax=ax, frameon=False, annotate_var_explained=True, show=False,
                title = 'Human MTG', save= "_hodge.png",  legend_loc='on data')
    
    sc.pl.pca(adata, color='cluster', ax=ax, frameon=False, annotate_var_explained=True, show=False,
                title = 'Human MTG', save= "_hodge.png",  legend_loc='on data')
    
    print("Saving PCA'd data as", pca_savedir + "hodge.h5ad")
    """

if __name__=="__main__":
    main()