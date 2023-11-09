# Figu 2: PCA
DATASETS = ["bakken", "tosches", "tasic", "colquitt"]
REFERENCE = "tasic"
OTHERS = ["bakken", "tosches", "colquitt"]
# TODO: save e.g. gene lists to scratchpath
species = {
    "colquitt": "Zebra finch",
    "tosches": "Turtle",
    "tasic": "Mouse",
    "bugeon": "Mouse L1-3",
    "bakken": "Human",
    "yao": "mouse",
    "hodge": "Human MTG",
}
CONTROLS = ["complete", "meis2", "abundance", "depth", "integrated_rpca"]
# "integrated_rpca", "integrated_cca"]
DESCRIPTION = {
    "complete": "",
    "meis2": "No Meis2",
    "integrated_rpca": "Integrated rPCA",
    "integrated_cca": "Integrated CCA",
    "abundance": "Matched Abundance",
    "depth": "Matched Depth",
}
MEM = 20000


def script_path(x):
    return f"../scripts/{x}"


rule all:
    input:
        expand(
            "figures/figure2/pca_{dataset}_{control}.png",
            dataset=DATASETS,
            control=[c for c in CONTROLS if "integrated" not in c],
        ),
        expand(
            "figures/figure2/pca_{dataset}_integrated_rpca.png",
            dataset=DATASETS,
        ),
        expand("figures/figure2/principal_angles_{control}.png", control=CONTROLS),
        expand("figures/figure2/cross_variance_{control}.png", control=CONTROLS),
        expand(
            "figures/figure2/compare_angles_{control}.png",
            control=CONTROLS,
        ),
        expand("figures/figure2/predict_from_{dataset}.png", dataset=DATASETS),

rule schematics:
    output:
        pc_sketch="figures/figure2/pcs_schematic.png",
        pc_gradient="figures/figure2/pcs_gradient.png",
    script:
        script_path("fig2_schematics.py")


rule integrate:
    input:
        reference=f"data/anndata/{REFERENCE}.h5ad",
        others=expand("data/anndata/{dataset}.h5ad", dataset=OTHERS),
    params:
        method=lambda wildcards: wildcards.method,
    output:
        expand(
            "results/anndata_integrated/{dataset}_integrated_{{method}}.h5ad",
            dataset=DATASETS,
        ),
    resources:
        mem_mb=64000,
    script:
        script_path("fig2_integrate.R")


rule fig2_dag:
    input:
        "figures/figure2/principal_angles_complete.png",
    output:
        dag="figures/dags/dag_fig2.png",
    shell:
        """
        snakemake {input} --dag | dot -Tpng -Gdpi=300 > {output}
        """


rule comparison:
    input:
        cv_complete="results/pc_comparison/cross_variance_complete.pickle",
        angles_complete="results/pc_comparison/principal_angles_complete.pickle",
        angles_control="results/pc_comparison/principal_angles_{control}.pickle",
        cv_control="results/pc_comparison/cross_variance_{control}.pickle",
    params:
        control=lambda wildcards: DESCRIPTION[wildcards.control],
    resources:
        mem_mb=MEM,
    output:
        angles="figures/figure2/compare_angles_{control}.png",
        cv="figures/figure2/compare_variance_{control}.png",
    script:
        script_path("fig2_comparison.py")


rule principal_angles:
    # TODO: use same script as fig3
    input:
        expand("results/anndata/{dataset}_{{control}}.h5ad", dataset=DATASETS),
    params:
        areas=species,
        reference="tasic",
    resources:
        mem_mb=MEM,
    output:
        figure="figures/figure2/principal_angles_{control}.png",
        angles="results/pc_comparison/principal_angles_{control}.pickle",
    script:
        script_path("fig2_principal_angles.py")


rule cross_variance:
    input:
        expand("results/anndata/{dataset}_{{control}}.h5ad", dataset=DATASETS),
    params:
        areas=species,
        reference="tasic",
    resources:
        mem_mb=MEM,
    output:
        figure="figures/figure2/cross_variance_{control}.png",
        data="results/pc_comparison/cross_variance_{control}.pickle",
    script:
        script_path("fig3_cross_variance.py")


rule pca_on_integrated:
    input:
        raw_anndata="results/anndata_integrated/{dataset}_integrated_{method}.h5ad",
        shared_genes="results/gene_lists/shared_genes.txt",
    params:
        species=lambda wildcards: species[wildcards.dataset],
        control=lambda wildcards: f"integrated_{wildcards.method}",
    resources:
        mem_mb=MEM,
    output:
        figure="figures/figure2/pca_{dataset}_integrated_{method}.png",
        anndata="results/anndata/{dataset}_integrated_{method}.h5ad",
    script:
        script_path("fig2_pca.py")


rule pca:
    # Input function based on contol (integrated or not)
    # if integrated - no need for shared genes
    input:
        raw_anndata="data/anndata/{dataset}.h5ad",
        shared_genes="results/gene_lists/shared_genes.txt",
        bugeon="data/anndata/bugeon.h5ad",
    params:
        species=lambda wildcards: species[wildcards.dataset],
        control=lambda wildcards: wildcards.control,
    resources:
        mem_mb=MEM,
    output:
        figure="figures/figure2/pca_{dataset}_{control}.png",
        anndata="results/anndata/{dataset}_{control}.h5ad",
    script:
        script_path("fig2_pca.py")


rule predict_across_species:
    # Predict state modulation in mouse (Bugeon) data
    # from tPC1 from other species
    input:
        bugeon="results/anndata/bugeon_log.h5ad",
        adata="results/anndata/{dataset}_72g.h5ad",
    params:
        species=lambda wildcards: species[wildcards.dataset],
    output:
        regression="figures/figure2/predict_from_{dataset}.png",
    script:
        script_path("fig2_predict_from_other.py")


rule datasets_table:
    input:
        expand("results/pandas/overview_{dataset}.csv", dataset=DATASETS),
    output:
        table="results/pandas/overview.csv",
        latex="results/pandas/overview.tex",
    script:
        script_path("fig2_dataset_table.py")
