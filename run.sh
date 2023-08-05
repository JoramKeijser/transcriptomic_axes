#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=3850
#SBATCH --time=01:00:00
#SBATCH --output=log/main_%j

source ~/mambaforge/bin/activate
#module load mamba/4.14.0
eval "$(conda shell.bash hook)"
conda activate snakemake_cluster

# set up proxy (no internet on nodes)
export https_proxy=http://frontend01:3128/
export http_proxy=http://frontend01:3128/ 


snakemake --profile ./env/slurm "$@"

