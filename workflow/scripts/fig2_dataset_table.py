"""
Description of each dataset
* Species
* Areas (not here?)
* Technology
* # neurons
* sequencing depth
* # genes
"""
import pandas as pd
import seaborn as sns


def main():
    for i, dataset in enumerate(snakemake.input):
        if i == 0:
            df = pd.read_csv(dataset)
        else:
            df2 = pd.read_csv(dataset)
            df = pd.concat((df, df2), axis=0)
    df.to_csv(snakemake.output.table, index=False)
    df.to_latex(snakemake.output.latex, index=False, float_format="{:.0f}".format)


if __name__ == "__main__":
    main()
