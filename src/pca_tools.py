import numpy as np
from sklearn.decomposition import PCA
import scanpy as sc
from src import constants


def do_pca(adata):
    """
    Perform PCA after preprocessing
    """
    sc.pp.normalize_total(adata, target_sum = constants.NORMALIZE_TARGET_SUM)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes = constants.NUM_HVG_GENES)
    sc.pp.pca(adata, n_comps = constants.NUM_PCS)
    return adata

def compare_variance(adata_ref, adata, pca_ref, n_pcs = 10):
        """
        How much variance is explained by PCs of a second datasets?
        Arguments:
            adata_ref (anndata): reference dataset
            adata (anndata): the other dataset
            pca_ref (PCA): fitted PCA object
            n_pcs (int): number of PCs to compare
        Returns:
            normalized cross_variance "(n_pcs)
        """
        # get PCs from reference
        C = pca_ref.get_covariance()
        variance = pca_ref.components_[0]@C@pca_ref.components_[0]
        # PCA on other dataset
        pca = PCA().fit(adata.X)
        PCs = pca.components_
        cross_variance = np.diag(pca.components_[:n_pcs]@C@pca.components_[:n_pcs].T)
        # Normalize by first ref. PC
        return cross_variance / variance


def orient_axes(adata):
    """
    Orient each of the first 2 PCs for consistent visualization
    Arguments:
        adata (anndata): AnnData set
    Returns:
        AnnData with Pv cells < Vip on PC1, Pv < Lamp5 on PC2
    """
    assert "X_pca" in adata.obsm
    if adata[adata.obs['Subclass'] == "Pvalb"].obsm['X_pca'][:,0].mean() > adata[adata.obs['Subclass'] == "Vip"].obsm['X_pca'][:,0].mean():
        # Orient to have Pvalb left, Vip right
        print("Flip PC1")
        adata.obsm['X_pca'][:,0] *= -1
        adata.varm['PCs'][:,0] *= -1
    if adata[adata.obs['Subclass'] == "Pvalb"].obsm['X_pca'][:,1].mean() > adata[adata.obs['Subclass'] == "Lamp5"].obsm['X_pca'][:,1].mean():
        # Pvalb top, Sst bottom
        print("Flip PC2")
        adata.obsm['X_pca'][:,1] *= -1
        adata.varm['PCs'][:,1] *= -1

    if np.sum(adata.obs['Subclass'] == "Meis2") > 0:
        if adata[adata.obs['Subclass'] == "Pvalb"].obsm['X_pca'][:,0].mean() > adata[adata.obs['Subclass'] == "Meis2"].obsm['X_pca'][:,0].mean():
            adata.obsm['X_pca'][:,0] *= -1
            adata.varm['PCs'][:,0] *= -1
    return adata
