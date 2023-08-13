SUBCLASSES = ["Pvalb", "Sst", "Lamp5", "Vip", "Sncg"]
TRANSFORMS = ["log", "raw"]
NSUBSETS = 1000
NPERMUTATIONS = 1000


def script_path(x):
    return f"workflow/scripts/{x}"


rule all:
    input:
        "results/pandas/activity.h5ad",
        "figures/figure1/state_modulation.png",
        "figures/figure1/example_trial.png",
        "data/anndata/bugeon.h5ad",
        "results/pandas/activity.h5ad",
        expand("results/anndata/bugeon_{transform}.h5ad", transform=TRANSFORMS),
        expand("figures/figure1/receptors_{transform}.png", transform=TRANSFORMS),
        "figures/figure1/regression_All.png",
        "figures/figure1/tracksplot_tasic.png",


rule state_modulation:
    input:
        data_dir="data/bugeon/",
    params:
        stimulus="Blank",
    output:
        data="results/pandas/activity.h5ad",
        figure="figures/figure1/state_modulation.png",
    script:
        script_path("fig1_state_modulation.py")


rule example_trial:
    input:
        home_dir="data/bugeon/",
        experiment_dir="data/bugeon/SB026/2019-10-16/",
    params:
        experiment="SB026/2019-10-16/",
        stimulus="Blank",
    output:
        figure="figures/figure1/example_trial.png",
    script:
        script_path("fig1_example_trial.py")


rule pca:
    input:
        transcriptomics="data/anndata/bugeon.h5ad",
        activity="results/pandas/activity.h5ad",
    params:
        transform=lambda wildcards: wildcards.transform,
    output:
        annotated="results/anndata/bugeon_{transform}.h5ad",
        by_subtype="results/pandas/bugeon_by_subtype_{transform}.csv",
        pca_subclass="figures/figure1/pca_subclass_{transform}.png",
        pca_modulation="figures/figure1/pca_modulation_{transform}.png",
        regression="figures/figure1/predict_mod_{transform}.png",
    script:
        script_path("fig1_pca.py")


rule receptors:
    input:
        by_subtype="results/pandas/bugeon_by_subtype_{transform}.csv",
        tasic="data/anndata/tasic.h5ad",
    output:
        figure="figures/figure1/receptors_{transform}.png",
    params:
        transform=lambda wildcards: wildcards.transform,
    resources:
        mem_mb=8000,
    script:
        script_path("fig1_receptors.py")


# TODO: add receptors_supp
rule all_receptors:
    input:
        by_subtype="results/pandas/bugeon_by_subtype_log.csv",
        tasic="data/anndata/tasic.h5ad",
    params:
        n_permutations=NPERMUTATIONS,
        figdir="figures/figure1/",
    resources:
        mem_mb=8000,
    output:
        tracksplot="figures/figure1/tracksplot_tasic.png",
        tracksplot_significant="figures/figure1/tracksplot_tasic_significant_first.png",
        hist_sig="figures/figure1/receptors_examples_significant.png",
        hist_notsig="figures/figure1/receptors_examples_not_significant.png",
    script:
        script_path("fig1_receptors_supp.py")


rule regression:
    input:
        "results/anndata/bugeon_log.h5ad",
    output:
        "figures/figure1/regression_All.png",
        expand(
            "figures/figure1/regression_{subclass}.png",
            subclass=SUBCLASSES,
            transform=["log"],
        ),
        "figures/figure1/ncell_corr.png",
        "figures/figure1/corr_subsample.png",
    params:
        n_subsets=NSUBSETS,
    script:
        script_path("fig1_regression.py")
