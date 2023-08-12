DATASETS = ["yao", "tasic", "tosches", "bakken", "colquitt"]
# TODO: add 'hodge' - also
PERMUTATIONS = 10
SAMPLES = 10
species = {
    "colquitt": "Zebra finch",
    "tosches": "Turtle",
    "tasic": "Mouse VISp L1-6",
    "bugeon": "Mouse VISp L1-3",
    "bakken": "Human",
    "yao": "Mouse Ctx & Hpc",
}
MEM = 20000


rule all:
    input:
        expand("figures/figure4/dotplot_{dataset}_1000.png", dataset=DATASETS),
        "figures/figure4/schematic.png",
        expand(
            "figures/figure4/subsample_vln_{dataset}_n1000_p1000.png", dataset=DATASETS
        ),


rule downsampling_schematic:
    input:
        tasic="data/anndata/tasic.h5ad",
        bakken="data/anndata/bakken.h5ad",
    output:
        figure="figures/figure4/schematic.png",
    resources:
        mem_mb=32000,
    script:
        "fig4_downsampling_schematic.py"


rule downsampling:
    input:
        dataset="data/anndata/tasic.h5ad",  # deep reference dataset TODO: use yao for hodge
        shallow="data/anndata/{dataset}.h5ad",  # shallower dataset
        shared_genes="results/gene_lists/shared_genes.txt",  #TODO: shared receptors
        receptors="results/gene_lists/significant_receptors_{permutations}.txt",
    resources:
        mem_mb=32000,
    output:
        vln="figures/figure4/subsample_vln_{dataset}_n{num_samples}_p{permutations}.png",
    params:
        num_samples=lambda wildcards: int(wildcards.num_samples),
        dataset=lambda wildcards: wildcards.dataset,
        species=lambda wildcards: species[wildcards.dataset],
    script:
        "fig4_downsampling.py"


rule dotplot:
    input:
        anndata="data/anndata/{dataset}.h5ad",
        shared_genes="results/gene_lists/shared_genes.txt",  #TODO: receptors only
        receptors="results/gene_lists/significant_receptors_{permutations}.txt",
    resources:
        mem_mb=16000,
    output:
        dotplot="figures/figure4/dotplot_{dataset}_{permutations}.png",
    script:
        "fig4_dotplots.py"


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
        "fig4_find_receptors.py"


# TODO: parallelize
rule allgenes:
    input:
        bugeon="data/anndata/bugeon.h5ad",
        tasic="data/anndata/tasic.h5ad",
    resources:
        mem_mb=8000,
    output:
        figure="figures/figure4/all_genes.png",
        significant_genes="results/gene_lists/significant_genes.txt",
    script:
        "fig4_allgenes.py"
