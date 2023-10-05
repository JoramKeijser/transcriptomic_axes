DATASETS = ["yao", "tasic", "tosches", "bakken", "colquitt"]
PERMUTATIONS = 10000 
# Estimating null distribution of corr. between modulation and expression
SAMPLES = 100
# Estimating RNA count distribution by downsampling 
species = {
    "colquitt": "Zebra finch",
    "tosches": "Turtle",
    "tasic": "Mouse VISp L1-6",
    "bugeon": "Mouse VISp L1-3",
    "bakken": "Human",
    "yao": "Mouse Ctx & Hpc",
}
MEM = 20000


def script_path(x):
    return f"../scripts/{x}"


rule all:
    input:
        "figures/figure4/corr_1000.png",
        "figures/figure4/all_genes.png",
        expand("figures/figure4/dotplot_{dataset}_1000.png", dataset=DATASETS),
        "figures/figure4/schematic.png",
        expand(
            "figures/figure4/subsample_vln_{dataset}_n100_p1000.png", dataset=DATASETS
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
        script_path("fig4_downsampling_schematic.py")


rule downsampling:
    input:
        dataset="data/anndata/tasic.h5ad",  # deep reference dataset TODO: use yao for hodge
        shallow="data/anndata/{dataset}.h5ad",  # shallower dataset
    resources:
        mem_mb=64000,
    output:
        vln="figures/figure4/subsample_vln_{dataset}_n{num_samples}_p{permutations}.png",
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
        figure="figures/figure4/all_genes.png",
        significant_genes="results/gene_lists/significant_genes.txt",
    script:
        script_path("fig4_allgenes.py")
