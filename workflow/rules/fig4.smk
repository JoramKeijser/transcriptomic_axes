DATASETS = ["yao", "tasic", "tosches", "bakken", "colquitt", "hodge"]
PERMUTATIONS = 10000
# Estimating null distribution of corr. between modulation and expression
SAMPLES = 100
# Estimating RNA count distribution by downsampling
species = {
    "colquitt": "Zebra finch",
    "tosches": "Turtle",
    "tasic": "Mouse VISp L1-6",
    "bugeon": "Mouse VISp L1-3",
    "bakken": "Human M1",
    "hodge": "Human MTG",
    "yao": "Mouse Ctx & Hpc",
}
MEM = 20000


def script_path(x):
    return f"../scripts/{x}"


rule all:
    input:
        expand("figures/figure4/dotplot_{dataset}_1000.png", dataset=DATASETS),
        "figures/figure4/schematic.png",
        expand(
            "figures/figure4/sub_{dataset}_tasic_n100_p1000.png",
            dataset=[d for d in DATASETS if d not in ["hodge", "tasic"]],
        ),
        # Also use hodge as reference for other human dataset
        "figures/figure4/sub_bakken_hodge_n100_p1000.png",
        "figures/figure4/all_genes.svg",
        #"figures/figure4/pca_tasic_chrna4_bugeon_abundance.png",


rule downsampling_schematic:
    input:
        tasic="data/anndata/tasic.h5ad",
        bakken="data/anndata/bakken.h5ad",
    output:
        figure="figures/figure4/schematic.png",
    resources:
        mem_mb=32000,
    script:
        script_path("fig4_downsampling_schematic.py")


rule downsampling:
    input:
        dataset="data/anndata/{reference}.h5ad",  # deep reference dataset
        shallow="data/anndata/{dataset}.h5ad",  # shallower dataset
    resources:
        mem_mb=64000,
    output:
        vln="figures/figure4/sub_{dataset}_{reference}_n{num_samples}_p{permutations}.png",
    params:
        num_samples=lambda wildcards: int(wildcards.num_samples),
        dataset=lambda wildcards: wildcards.dataset,
        species=lambda wildcards: species[wildcards.dataset],
    script:
        script_path("fig4_downsampling.py")


rule downsampling_human:
    # Downsample deeper Hodge dataset to compare with Bakken
    # Could combine this with the previous rule
    input:
        dataset="data/anndata/hodge.h5ad",
        shallow="data/anndata/bakken.h5ad",
    resources:
        mem_mb=64000,
    output:
        vln="figures/figure4/subsample_vln_hodge_n{num_samples}_p{permutations}_human.png",
    params:
        num_samples=lambda wildcards: int(wildcards.num_samples),
        dataset=lambda wildcards: wildcards.dataset,
        species=lambda wildcards: species[wildcards.dataset],
    script:
        script_path("fig4_downsampling.py")


rule dotplot:
    input:
        anndata="data/anndata/{dataset}.h5ad",
        shared_genes="results/gene_lists/shared_genes.txt",  #TODO: receptors only
        receptors="results/gene_lists/significant_receptors_{permutations}.txt",
    resources:
        mem_mb=32000,
    output:
        dotplot="figures/figure4/dotplot_{dataset}_{permutations}.png",
    script:
        script_path("fig4_dotplots.py")


rule find_receptors:
    input:
        bugeon="data/anndata/bugeon.h5ad",
        tasic="data/anndata/tasic.h5ad",
        bugeon_grouped="results/pandas/bugeon_by_subtype_log.csv",
    params:
        permutations=lambda wildcards: int(wildcards.permutations),
    resources:
        mem_mb=8000,
    output:
        figure="figures/figure4/corr_{permutations}.png",
        receptors="results/gene_lists/significant_receptors_{permutations}.txt",
    script:
        script_path("fig4_find_receptors.py")


rule allgenes:
    input:
        bugeon="data/anndata/bugeon.h5ad",
        tasic="data/anndata/tasic.h5ad",
    resources:
        mem_mb=8000,
    output:
        png_fig="figures/figure4/all_genes.png",
        svg_fig="figures/figure4/all_genes.svg",
        significant_genes="results/gene_lists/significant_genes.txt",
    script:
        script_path("fig4_allgenes.py")


""" rule receptor_viz:
    input:
        "results/anndata/tasic_bugeonabundance.h5ad",
    output:
        "figures/figure4/pca_tasic_chrna4_bugeon_abundance.png",
    script:
        script_path("fig4_pca_tasic_chrna4.py") """
