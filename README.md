# Transcriptomic axes of GABAergic interneuron diversity

Neural activity is influenced not only by external stimuli but also by an animal's internal state. Much is known about how this happens in the mouse brain, but it is unclear if the same principles
apply to other species, including humans. To bridge this gap, we compare gene expression patterns in neurons from different species. The following
[preprint](https://doi.org/10.1101/2023.12.04.569849) describes our findings; this repository contains the corresponding code.

>Keijser, J., Hertäg, L., & Sprekeler, H. (2023). Transcriptomic correlates of state modulation in GABAergic interneurons:
>A cross-species analysis. bioRxiv, 2023-12.

<p align="center">
  <img width="800" src="./figures/tpc_fig0.png">
</p>

The actual analysis code is combined into a pipeline using [Snakemake](https://snakemake.readthedocs.io/en/stable/), which provides a reproducible and scalable way of running large-scale, interdependent analyses across multiple datasets. Briefly, a Snakemake pipeline consists of multiple input/output rules, forming a directed acyclic graph (DAG) through which data flows --- from raw input to the final result. The DAG of the first supplementary figure is shown below.

<p align="center">
  <img width="600" src="./figures/dags/fig1.png">
</p>

After installation, this figure can be recreated by running:
```
snakemake -s workflow/rules/fig1.smk -f --dag | dot -Tpng -Gdpi=300 > figures/dags/fig1.png
```
## Organization

The repository is organized as follows:
```
├── data: raw data, will be created by running the pipeline
├── env: Snakemake configuration file(s) for running on a computing cluster
├── figures: will be created when running the pipeline
├── log: will be created when running the pipeline on a cluster
├── profile: Snakemake configuration file(s) for running locally
├── renv: files to reproduce the R environment
├── results: will contain output data after running the pipeline
├── src: source code
├── src.egg-info: project package
├── workflow
  ├── rules: Snakemake rules defining data flow
  └── scripts: Python/R/Bash scripts 
```

## Installation

To start, clone the repository:
```
git clone github.com/JoramKeijser/transcriptomic_axes/
```
Then recreate the [conda/mambda](https://github.com/mamba-org/mamba) environment with the required Python packages:
```
cd transcriptomic_axes/
mamba env create --name tpc --file environment.yml
mamba activate tpc
```
Install the project package:
```
pip install -e .
```

Finally, recreate the R environment using [renv](https://rstudio.github.io/renv/index.html):
```
R -e 'renv::init()'
```


## Running the code
The pipeline is distributed across multiple rules, one for each (main) figure. There are also two preliminary rules for downloading and preprocessing the datasets. Individual rules are not connected; the `download` and `preprocess` workflows should be executed before any others.

To generate all results (figures and data output) corresponding to, say, figure 1, run:

```
snakemake -s workflow/rules/fig1.smk --cores
```
The `cores` argument tells Snakemake to use all available cores. Executed this way, the code will run on your local machine, but it's advisable to use a computing cluster to maximally benefit from Snakemake's parallelization power (and to have sufficient RAM for handling the larger datasets). The `run.sh` script submits jobs to a [SLURM](https://slurm.schedmd.com/overview.html) cluster:
```
sbatch run.sh -s workflow/rules/fig1.smk
```
In this case, the output from rule `all` will not be printed to the command line but to a file in the `log/` directory; the output of individual jobs will be printed to files in `log/jobs/`. 

## Acknowledgements

We are grateful to everyone who previously collected & published the data analysed here:
* Bakken, Trygve E., et al. "Comparative cellular analysis of motor cortex in human, marmoset and mouse." Nature 598.7879 (2021): 111-119.
* Bugeon, Stephane, et al. "A transcriptomic axis predicts state modulation of cortical interneurons." Nature 607.7918 (2022): 330-338.
* Colquitt, Bradley M., et al. "Cellular transcriptomics reveals evolutionary identities of songbird vocal circuits." Science 371.6530 (2021): eabd9704.
* Hodge, Rebecca D., et al. "Conserved cell types with divergent features in human versus mouse cortex." Nature 573.7772 (2019): 61-68.
* Tasic, Bosiljka, et al. "Shared and distinct transcriptomic cell types across neocortical areas." Nature 563.7729 (2018): 72-78.
* Tosches, Maria Antonietta, et al. "Evolution of pallium, hippocampus, and cortical cell types revealed by single-cell transcriptomics in reptiles." Science 360.6391 (2018): 881-888.
* Yao, Zizhen, et al. "A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation." Cell 184.12 (2021): 3222-3241.
