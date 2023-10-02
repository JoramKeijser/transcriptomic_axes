# tpc_pipeline

Run a particular rule (using all cores): `snakemake -s workflow/rules/fig2.smk --cores`. 

Create a particular output: `snakemake -s workflow/rules/fig1.smk --cores figures/figure2/principal_angles_complete.png`.

Submit to SLURM cluster `sbatch run.sh -s workflow/rules/fig1.smk --cores figures/figure2/principal_angles_complete.png`.

