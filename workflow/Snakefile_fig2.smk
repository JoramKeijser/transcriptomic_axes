# Figu 2: PCA
DATASETS = ["bakken", "tosches", "tasic", "colquitt"]
REFERENCE = 'tasic'
OTHERS = ["bakken", "tosches", "colquitt"]
# TODO: save e.g. gene lists to scratchpath
species = {'colquitt': "Zebra finch", "tosches": "Turtle",
            "tasic": "Mouse", "bugeon": "Mouse L1-3",
            "bakken": "Human", "yao": "mouse", "hodge": "Human MTG"}
CONTROLS = ['complete', 'meis2', 'abundance', 'depth',
    "integrated_rpca", "integrated_cca"]
DESCRIPTION = {"complete": "", "meis2": "No Meis2",
    "integrated_rpca": "Integrated rPCA",
    "integrated_cca": "Integrated CCA",
    "abundance": "Matched Abundance", "depth": "Matched Depth"}
# TODO: flag to avoid naming "complete"

rule all:
    input:
        expand("figures/figure2/pca_{dataset}_{control}.png",
        dataset=DATASETS, control=CONTROLS),
        expand("figures/figure2/principal_angles_{control}.png",
        control=CONTROLS),
        expand("figures/figure2/cross_variance_{control}.png",
        control=CONTROLS),
        expand("figures/figure2/compare_angles_{control}.png",
        control=['meis2', 'abundance', 'depth'] )


rule integrate:
    # TODO: Name input files - or not
    # TODO: still need to do PCA on each separately.
    input:
        reference = f"data/anndata/{REFERENCE}.h5ad",
        others = expand("data/anndata/{dataset}.h5ad",
            dataset=OTHERS)
    params:
        method = lambda wildcards : wildcards.method
    output:
        expand(
            "results/anndata/{dataset}_integrated_{{method}}.h5ad",
            dataset = DATASETS
        )
    script:
        "fig2_integrate.R"


rule fig2_dag:
    input:
        "figures/figure2/principal_angles_complete.png"
    output:
        dag = "figures/dags/dag_fig2.png"
    shell:
        """
        snakemake {input} --dag | dot -Tpng -Gdpi=300 > {output}
        """

rule comparison:
    input:
        cv_complete = 'results/pc_comparison/cross_variance_complete.pickle',
        angles_complete = 'results/pc_comparison/principal_angles_complete.pickle',
        angles_control = 'results/pc_comparison/principal_angles_{control}.pickle',
        cv_control = 'results/pc_comparison/cross_variance_{control}.pickle'
    params:
        control = lambda wildcards : DESCRIPTION[wildcards.control]
    output:
        angles = "figures/figure2/compare_angles_{control}.png",
        cv = "figures/figure2/compare_variance_{control}.png"
    script:
        "fig2_comparison.py"


rule principal_angles:
    #TODO: use same script as fig3
    input:
        expand("results/anndata/{dataset}_{{control}}.h5ad",
            dataset=DATASETS)
    params:
        areas = species,
        reference = "tasic"
    output:
        figure= "figures/figure2/principal_angles_{control}.png",
        angles = 'results/pc_comparison/principal_angles_{control}.pickle'
    script:
        "fig2_principal_angles.py"

rule cross_variance:
    input:
        expand("results/anndata/{dataset}_{{control}}.h5ad",
            dataset=DATASETS)
    params:
        areas = species,
        reference = "tasic"
    output:
        figure= "figures/figure2/cross_variance_{control}.png",
        data = 'results/pc_comparison/cross_variance_{control}.pickle'
    script:
        "fig3_cross_variance.py"

rule pca:
    input:
        raw_anndata = "data/anndata/{dataset}.h5ad",
        shared_genes="results/gene_lists/shared_genes.txt"
    params:
        species = lambda wildcards : species[wildcards.dataset],
        control = lambda wildcards : wildcards.control
    output:
        figure = "figures/figure2/pca_{dataset}_{control}.png",
        anndata = "results/anndata/{dataset}_{control}.h5ad",
    script:
        "fig2_pca.py"

rule intersect_genes:
    input:
        expand("results/gene_lists/genes_{dataset}.csv",
        dataset=DATASETS),
    output:
        shared_genes="results/gene_lists/shared_genes.txt",
    script:
        "fig2_intersect_genes.py"


rule datasets_table:
    input:
        expand("results/pandas/overview_{dataset}.csv", dataset=DATASETS),
    output:
        table="results/pandas/overview.csv",
        latex="results/pandas/overview.tex",
    script:
        "fig2_dataset_table.py"

rule datasets:
    input:
        anndata="data/anndata/{dataset}.h5ad",
    output:
        figure="figures/figure2/QC_{dataset}.png",
        table="results/pandas/overview_{dataset}.csv",
        genes="results/gene_lists/genes_{dataset}.csv",
    params:
        dataset=lambda wildcards: wildcards.dataset,
    script:
        "fig2_dataset_stats.py"
