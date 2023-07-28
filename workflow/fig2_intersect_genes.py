"""
Description of each dataset
* Species
* Areas (not here?)
* Technology
* # neurons
* sequencing depth
* # genes
"""
import numpy as np
import pandas as pd
import anndata as ad
import os
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_palette("colorblind")

def main():

    for i, genelist in enumerate(snakemake.input):
        if i == 0:
            shared_genes = np.loadtxt(genelist, dtype=str)
        else:
            new_genes = np.loadtxt(genelist, dtype=str)
            shared_genes = set(shared_genes).intersection(new_genes)

        
    print(f"Found {len(shared_genes)} shared genes")
    np.savetxt(snakemake.output.shared_genes, 
               list(shared_genes), fmt="%s")
    
if __name__ == "__main__":
    main()