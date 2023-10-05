"""
Mouse datasets
"""
areas = {
    "tasic": "VISp L1-6",
    "bugeon": "VISp L1-3",
    "yao": "Ctx & Hpc",
    "bakken": "Human M1",
    "hodge": "Human MTG",
    "colquitt": "Zebra Finch",
    "tosches": "Turtle",
}
CONTROLS = ["complete", "72g", "human"]
# TODO: bugeon_log -> bugeon
MEM = 20000


def script_path(x):
    return f"../scripts/{x}"


def genelist_from_label(wildcards):
    lists = {
        "complete": "results/gene_lists/shared_mouse_genes.txt",
        "72g": "data/bugeon/genes.names.txt",
        "log": "data/bugeon/genes.names.txt",
        "bugeonabundance": "results/gene_lists/shared_mouse_genes.txt",
        "bugeonsst": "results/gene_lists/shared_mouse_genes.txt",
        "human": "results/gene_lists/shared_mouse_genes.txt",
    }
    return lists[wildcards.control]


def datasets_from_condition(wildcards):
    lists = {
        "complete": ["tasic", "yao"],
        "72g": ["tasic", "yao", "bugeon"],
        "human": ["bakken", "hodge"],
        "bugeonabundance": ["tasic"],
        "bugeonsst": ["tasic"],
    }
    return expand(
        expand(
            "results/anndata/{dataset}_{control}.h5ad",
            dataset=lists[wildcards.control],
            control=wildcards.control,
        )
    )


def areas_from_condition(wildcards):
    areas = {
        "complete": {"tasic": "VISp L1-6", "yao": "Ctx & Hpc"},
        "72g": {"tasic": "VISp L1-6", "bugeon": "VISp L1-3", "yao": "Ctx & Hpc"},
        "human": {"bakken": "Human M1", "hodge": "Human MTG"},
        "bugeonsst": {"tasic": "VISp L1-6", "bugeon": "VISp L1-3"},
        "bugeonabundance": {"tasic": "VISp L1-6", "bugeon": "VISp L1-3"},
    }
    return areas[wildcards.control]


def reference_from_condition(wildcards):
    if wildcards.control == "human":
        return "bakken"
    else:
        return "tasic"


rule all:
    input:
        expand("figures/figure3/cross_variance_{control}.png", control=CONTROLS),
        expand("figures/figure3/principal_angles_{control}.png", control=CONTROLS),


rule cross_variance:
    input:
        datasets_from_condition,
    params:
        areas=areas_from_condition,
        reference=reference_from_condition,
    resources:
        mem_mb=MEM * 4,
    output:
        figure="figures/figure3/cross_variance_{control}.png",
        data="results/pc_comparison/cross_variance_mouse_{control}.pickle",
    script:
        script_path("fig3_cross_variance.py")


rule principal_angles:
    input:
        same_species=datasets_from_condition,
        baseline_angles="results/pc_comparison/principal_angles_complete.pickle",
    params:
        control=lambda wildcards: wildcards.control,
        areas=areas_from_condition,
        reference=reference_from_condition,
    resources:
        mem_mb=MEM,
    output:
        figure="figures/figure3/principal_angles_{control}.png",
        angles="results/pc_comparison/principal_angles_mouse_{control}.pickle",
    script:
        script_path("fig3_principal_angles.py")


rule pca:
    input:
        raw_anndata="data/anndata/{dataset}.h5ad",
        shared_genes=genelist_from_label,
        bugeon="data/anndata/bugeon.h5ad",
    params:
        species=lambda wildcards: areas[wildcards.dataset],
        control=lambda wildcards: wildcards.control,
    output:
        figure="figures/figure3/pca_{dataset}_{control}.png",
        anndata="results/anndata/{dataset}_{control}.h5ad",
    resources:
        mem_mb=MEM * 4,
    script:
        script_path("fig2_pca.py")


rule intersect_mouse_genes:
    input:
        expand("results/gene_lists/genes_{dataset}.csv", dataset=["tasic", "yao"]),
    output:
        shared_genes="results/gene_lists/shared_mouse_genes.txt",
    script:
        script_path("fig2_intersect_genes.py")
