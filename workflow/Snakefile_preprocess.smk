# Raw data to anndata
import numpy as np

DATASETS = ["bakken", "bugeon", "colquitt", "hodge", "tasic", "tosches", "yao"]
NUM_JOBS = 100
NUM_CELLS = 1169213 #int(1169213 / 2)  # 1169213
ROWS_PER_JOB = int(NUM_CELLS / NUM_JOBS)
START_ROWS = np.arange(0, NUM_CELLS, ROWS_PER_JOB)
SMALLMEM = 4000
MEM = 60000
LARGEMEM = 240000


#rule num_rows:
#    input:
#        metadata="data/yao/metadata.csv",
#    output:
#        "results/pandas/num_jobs.txt",
#    shell:
#        """
#        num_cells=$(wc -l < {input})
#        num_jobs=100
#        rows_per_job=$((num_cells/{{NUM_JOBS}}+1)) # >1/100th of the 1.16 m cells
#        echo $rows_per_job > {output}
#        """


rule all:
    input:
        expand("data/anndata/{dataset}.h5ad",
 	dataset=DATASETS),


#rule intersect_genes:
#    input:
#        expand("results/pandas/genes_{dataset}.csv",
#        dataset=DATASETS),
#    output:
#        shared_genes="results/gene_lists/shared_genes.txt",
#    script:
#        "fig2_intersect_genes.py"


rule pp_yao:  # Merge partitions. Could use rows per job here
   input:
       files=expand(
           "data/scratch/yao_{start}_{num}.h5ad", start=START_ROWS, num=ROWS_PER_JOB
       ),
   resources:
       mem_mb=LARGEMEM,
   output:
       anndata="data/anndata/yao.h5ad",
   script:
       "preprocess_yao_combine.py"


rule partition:
    input:
        metadata="data/yao/metadata.csv",
        counts="data/yao/expression_matrix.hdf5",
        shared_genes="results/gene_lists/shared_genes.txt",
        bugeon_genes="data/bugeon/genes.names.txt",
    output:
        anndata="data/scratch/yao_{start_row}_{num_rows}.h5ad",
    resources:
        mem_mb=8000,  
    params:
        start_row=lambda wildcards: int(wildcards.start_row),
        num_rows=lambda wildcards: int(wildcards.num_rows),
    script:
        "preprocess_yao_partition.py"


rule pp_bugeon:
    input:
        data_dir="data/bugeon/",
    output:
        anndata="data/anndata/bugeon.h5ad",
    resources:
        mem_mb=SMALLMEM,
    script:
        "preprocess_bugeon.py"


rule pp_hodge:
    input:
        genes="data/hodge/human_MTG_2018-06-14_genes-rows.csv",
        exons="data/hodge/human_MTG_2018-06-14_exon-matrix.csv",
        introns="data/hodge/human_MTG_2018-06-14_intron-matrix.csv",
        cells="data/hodge/human_MTG_2018-06-14_samples-columns.csv",
    resources:
        mem_mb=MEM,
    output:
        anndata="data/anndata/hodge.h5ad",
    script:
        "preprocess_hodge.py"


rule pp_tasic:
    input:
        exons="data/tasic/mouse_VISp_2018-06-14_exon-matrix.csv",
        introns="data/tasic/mouse_VISp_2018-06-14_intron-matrix.csv",
        genes="data/tasic/mouse_VISp_2018-06-14_genes-rows.csv",
        cells="data/tasic/mouse_VISp_2018-06-14_samples-columns.csv",
    resources:
        mem_mb=MEM,
    output:
        anndata="data/anndata/tasic.h5ad",
        hv_genes="results/gene_lists/tasic_hvgs.txt",
    script:
        "preprocess_tasic.py"


rule extract_tosches:
    input:
        robj="data/tosches/turtle.neurons.Robj",
    output:
        counts="data/tosches/counts.csv",
        metadata="data/tosches/metadata.csv",
    resources:
        mem_mb=SMALLMEM,
    script:
        "extract_tosches.R"


rule pp_tosches:
    input:
        counts="data/tosches/counts.csv",
        metadata="data/tosches/metadata.csv",
    output:
        anndata="data/anndata/tosches.h5ad",
    resources:
        mem_mb=SMALLMEM,
    script:
        "preprocess_tosches.py"


rule pp_colquitt:
    input:
        data="data/colquitt/HVC_RA_RNA_counts.csv",
    output:
        anndata="data/anndata/colquitt.h5ad",
    resources:
        mem_mb=MEM,
    script:
        "preprocess_colquitt.py"


rule pp_bakken:
    input:
        counts="data/bakken/matrix.csv",
        metadata="data/bakken/metadata.csv",
    output:
        anndata="data/anndata/bakken.h5ad",
    resources:
        mem_mb=LARGEMEM,
    script:
        "preprocess_bakken.py"
