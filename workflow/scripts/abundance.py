"""
Count cell types
"""
import anndata as ad
import pandas as pd

# Create DF with cell type counts for each dataset
subclass_order = ["Pvalb", "Sst", "Lamp5", "Vip", "Sncg", "Meis2"]
names = [name[0].upper() + name[1:] for name in snakemake.params.names]
df = pd.DataFrame(columns=names)
for name, file in zip(names, snakemake.input):
    print(file)
    adata = ad.read_h5ad(file)
    counts = adata.obs['Subclass'].value_counts()
    # Add 0 counts 
    counts = counts.reindex(subclass_order, fill_value=0)
    percentages = counts / counts.sum() * 100
    df[name] = percentages

#df = df.reindex(subclass_order)
# Save as csv and tex
df.to_csv(snakemake.output.csv, index=True)
df.to_latex(snakemake.output.tex, index=True, float_format="%.1f")
print(df.round(2))



