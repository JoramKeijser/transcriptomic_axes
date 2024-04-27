#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=3850
#SBATCH --time=4:00:00 #to allow queueing - 1h suffices for compute
#SBATCH --output=log/main_%j
#SBATCH --partition=TestAndBuild # alternative: standard

snakemake --profile env/slurm "$@" --cores 
