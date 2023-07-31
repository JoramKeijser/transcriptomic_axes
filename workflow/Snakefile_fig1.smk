SUBCLASSES = ["Pvalb", "Sst", "Lamp5", "Vip", "Sncg"]
DATASETS = ["bugeon", "bakken", "tosches", "tasic", "colquitt"]  # , "yao", , ,
TRANSFORMS = ["log", "raw"]
#TODO: bugeon_log -> bugeon in pca rule

rule figure1:
    input:
        "results/pandas/activity.h5ad",
        "figures/figure1/state_modulation.png",
        "figures/figure1/example_trial.png",
        "data/anndata/bugeon.h5ad",
        expand("results/anndata/bugeon_{transform}.h5ad", transform=TRANSFORMS),
        expand(
            "figures/figure1/regression_{subclass}.png",
            subclass=SUBCLASSES,
            transform=["log"],
        ),
        expand("/figures/figure1/receptors_{transform}.png", transform=TRANSFORMS),


rule state_modulation:
    input:
        data_dir="data/bugeon/",
    params:
        stimulus="Blank",
    output:
        data="results/pandas/activity.h5ad",
        figure="figures/figure1/state_modulation.png",
    script:
        "fig1_state_modulation.py"


rule example_trial:
    input:
        "data/bugeon/SB026/2019-10-16/",
    params:
        experiment="SB026/2019-10-16/",
        stimulus="Blank",
    output:
        figure="figures/figure1/example_trial.png",
    script:
        "fig1_example_trial.py"


rule pca:
    input:
        figdir="figures/figure1/",
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
        "fig1_pca.py"


# TODO: fix problem that we now only have bugeon_log


rule regression:
    # TODO: fix num_subsets hard coding
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
        n_subsets=100,  #TODO mbda wildcards : wildcards.n_subsets
    script:
        "fig1_regression.py"


rule receptors:
    # TODO: tasic
    input:
        by_subtype="results/pandas/bugeon_by_subtype_{transform}.csv",
        tasic="results/anndata/tasic.h5ad",
    output:
        figure="/figures/figure1/receptors_{transform}.png",
    params:
        transform=lambda wildcards: wildcards.transform,
    script:
        "fig1_receptors.py"


# TODO: add receptors_supp
