"""
Mouse datasets
"""
DATASETS = ["tasic", "yao", "bugeon"]
species = {"tasic": "VISp L1-6", "bugeon": "VISp L1-3", 
             "yao": "Ctx & Hpc"}
CONTROLS = ['complete', '72g']
#TODO: bugeon_log -> bugeon

def genelist_from_label(wildcards):
    lists = {"complete": "results/gene_lists/shared_mouse_genes.txt",
            "72g": "data/bugeon/genes.names.txt",
            "bugeonabundance": "results/gene_lists/shared_mouse_genes.txt"}
    return lists[wildcards.control]


def datasets_from_condition(wildcards):
    lists = {"complete": ["tasic", "yao"],
            "72g": ["tasic", "yao", "bugeon"] }
    return expand(
        expand("results/anndata/{dataset}_{control}.h5ad",
            dataset=lists[wildcards.control], control=wildcards.control)
    )

def areas_from_condition(wildcards):
    areas = {"complete": {"tasic": "VISp L1-6", 
             "yao": "Ctx & Hpc"},
             "72g": {"tasic": "VISp L1-6", "bugeon": "VISp L1-3", 
             "yao": "Ctx & Hpc"}}
    return areas[wildcards.control]


rule cross_variance:
    input:
        datasets_from_condition,
    params:
        areas = areas_from_condition,
        reference = "tasic"
    output:
        figure= "figures/figure3/cross_variance_{control}.png",
        data = 'results/pc_comparison/cross_variance_mouse_{control}.pickle'
    script:
        "fig3_cross_variance.py"

# TODO: same script for fig2 and fig3 principal angles
rule principal_angles:
    input:
       datasets_from_condition,
    params:
        areas = areas_from_condition
    output:
        figure = "figures/figure3/principal_angles_{control}.png",
        angles = 'results/pc_comparison/principal_angles_mouse_{control}.pickle'
    script:
        "fig3_principal_angles.py"

rule pca:
    input:
        raw_anndata = "data/anndata/{dataset}.h5ad",
        shared_genes=genelist_from_label,
        bugeon = "data/anndata/bugeon.h5ad"
    params:
        species = lambda wildcards : species[wildcards.dataset],
        control = lambda wildcards : wildcards.control
    output:
        figure = "figures/figure3/pca_{dataset}_{control}.png",
        anndata = "results/anndata/{dataset}_{control}.h5ad",
    script:
        "fig2_pca.py"


rule intersect_mouse_genes:
    input:
        expand("results/gene_lists/genes_{dataset}.csv",
        dataset=["tasic", "yao"])
    output:
        shared_genes="results/gene_lists/shared_mouse_genes.txt"
    script:
        "fig2_intersect_genes.py" 
