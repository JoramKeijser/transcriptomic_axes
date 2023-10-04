# tpc_pipeline

Run a particular rule (using all cores): `snakemake -s workflow/rules/fig2.smk --cores`. 

Create a particular output: `snakemake -s workflow/rules/fig2.smk --cores figures/figure2/principal_angles_complete.png`.

Submit to SLURM cluster `sbatch run.sh -s workflow/rules/fig2.smk --cores`. In this case, the output from rule `all` will not be printed to the command line, but to a file in the `log/` directory; the output of individual jobs will be printed to files in `log/jobs/`.  
