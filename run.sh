#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=3850
#SBATCH --time=01:00:00
#SBATCH --output=log/main_%j
#SBATCH --partition=standard

snakemake --cores --profile env/slurm "$@" --cores 
