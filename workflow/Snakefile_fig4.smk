DATASETS = ['yao', 'tasic', 'tosches', "bakken", "colquitt"]
#TODO: add 'hodge' - also
PERMUTATIONS = 1000
species = {'colquitt': "Zebra finch", "tosches": "Turtle",
            "tasic": "Mouse VISp L1-6", "bugeon": "Mouse VISp L1-3",
            "bakken": "Human", "yao": "Mouse Ctx & Hpc"}

rule all:
    input:
        expand("figures/figure4/dotplot_{dataset}.png",
        dataset=DATASETS
        ),
        expand("figures/figure4/subsample_vln_{dataset}_100.png",
        dataset=DATASETS # TODO: no tasic
        ),
        "figures/figure4/schematic.png"

rule downsampling_schematic:
    input:
        tasic = "data/anndata/tasic.h5ad",
        bakken = "data/anndata/bakken.h5ad"
    output:
        figure = "figures/figure4/schematic.png"
    script:
        "fig4_downsampling_schematic.py"

rule downsampling:
    input:
        dataset = "data/anndata/tasic.h5ad", # deep reference dataset TODO: use yao for hodge
        shallow = "data/anndata/{dataset}.h5ad", # shallower dataset
        shared_genes = "results/gene_lists/shared_genes.txt", #TODO: shraed receptors
        receptors = "results/gene_lists/significant_receptors_1000.txt"
    output:
        vln = "figures/figure4/subsample_vln_{dataset}_{num_samples}.png"
    params:
        num_samples = lambda wildcards : int(wildcards.num_samples),
        dataset = lambda wildcards : wildcards.dataset,
        species = lambda wildcards : species[wildcards.dataset]
    script:
        "fig4_downsampling.py"

rule dotplot:
    input:
        anndata = "data/anndata/{dataset}.h5ad",
        shared_genes = "results/gene_lists/shared_genes.txt", #TODO: shraed receptors
        receptors = "results/gene_lists/significant_receptors_1000.txt"
    output:
        dotplot = "figures/figure4/dotplot_{dataset}.png"
    script:
        "fig4_dotplots.py"


rule find_receptors:
    input:
        bugeon = "data/anndata/bugeon.h5ad",
        tasic = "data/anndata/tasic.h5ad",
        bugeon_grouped = "results/pandas/bugeon_by_subtype_log.csv"
    params:
        permutations = lambda wildcards : int(wildcards.permutations)
    output:
        figure = "figures/figure4/corr_{permutations}.png",
        receptors = "results/gene_lists/significant_receptors_{permutations}.txt"
    script:
        "fig4_find_receptors.py"

# TODO: parallelize
rule allgenes:
    input:
        bugeon = "data/anndata/bugeon.h5ad",
        tasic = "data/anndata/tasic.h5ad",
    output:
        figure = "figures/figure4/all_genes.png",
        significant_genes = "results/gene_lists/significant_genes.txt"
    script:
        "fig4_allgenes.py"
