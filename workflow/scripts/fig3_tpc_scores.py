# Compute subclass-specific tPC1 score
# for different versions of the L1-3 and L1-6 data
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_palette("colorblind")


def subclass_pc_scores(adata, name, subclass, pc_num=1):
    """
    Compute tPCn scorees of a particular subclass

    Arguments:
        adata (anndata): annotated data frame
        name (str): description of the data
        subclass (str): Pvalb, Sst etc.
        pc_num (int): PC number

    Returns:
        DataFrame:

    """
    # PC1 projection of a particular subclass
    assert subclass in adata.obs["Subclass"].cat.categories
    df = pd.DataFrame()
    ix = adata.obs["Subclass"] == subclass
    pc_score = np.array(adata[ix].obsm["X_pca"][:, pc_num - 1])
    df2 = pd.DataFrame()
    df2["tPC"] = pc_score
    df2["Subclass"] = [subclass] * len(pc_score)
    df2["Dataset"] = [name] * len(pc_score)
    return df2


descriptions = {
    "bugeon_log": "L1-3",
    "complete": "L1-6",
    "bugeonsst": "Match subtypes",
    "72g": "Match genes",
    "bugeonabundance": "Match abundance",
}

for file in snakemake.input:
    adata = ad.read_h5ad(file)
    # Extract description from file name
    name = name.split("_")[1].split(".h5ad")[0]
    df_dataset = pd.DataFrame()
    for subclass in ["Pvalb", "Sst"]:
        df_subclass = subclass_pc_scores(adata, name, subclass, snakemake.params.pc_num)
        # Normalize by Pvalb stats
        if subclass == "Pvalb":
            m = df_subclass["tPC"].median()
            sd = df_subclass["tPC"].std()
        df_dataset = pd.concat((df_dataset, df_subclass))
    df_dataset["tPC"] -= m
    df_dataset["tPC"] /= sd
    df = pd.concat((df, df_dataset))


fig, ax = plt.subplots(figsize=(8, 5))
sns.boxplot(
    df,
    x="tPC",
    y="Dataset",
    hue="Subclass",
    fliersize=1,
    ax=ax,
    linewidth=1,
    saturation=0.9,
)
sns.despine()
ax.legend().remove()
ax.set_ylabel("")
ax.set_xlabel("tPC1 score (norm.)")
fig.tight_layout()
plt.savefig(snakemake.output[0], dpi=300)
