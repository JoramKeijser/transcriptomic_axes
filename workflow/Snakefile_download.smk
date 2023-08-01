SUBCLASSES = ["Pvalb", "Sst", "Lamp5", "Vip", "Sncg"]
DATASETS = ["bugeon", "bakken", "tosches", "tasic", "colquitt",
            "yao", "hodge"]  # , "yao", , ,


rule all:
    input:
        expand("data/{dataset}", dataset=DATASETS),


rule download_bugeon:
    output:
        directory("data/bugeon"),
    shell:
        """
        wget https://figshare.com/ndownloader/files/38046633
        unzip -q 38046633 -d {output}
        rm 38046633
        """

rule download_hodge:
    output:
        directory("data/hodge")
    shell:
        """
        wget -P {output} https://celltypes.brain-map.org/api/v2/well_known_file_download/694416044
        7z x {output}/694416044
        mv human*.{{txt,csv}} {output}
        rm -r {output}/694413985
        """

# TODO: one rm line. Fix.
rule download_tasic:
    output:
        directory("data/tasic"),
    shell:
        """
        mkdir {output}
        wget -P {output} https://celltypes.brain-map.org/api/v2/well_known_file_download/694413985
        7z x {output}/694413985
        mv mouse*.{{txt,csv}} {output}
        rm -r {output}/694413985
        """


rule download_bakken:
    output:
        directory("data/bakken"),
    shell:
        """
        wget -P {output} https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv
        wget -P {output} https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv
        """


rule download_colquitt:
    output:
        directory("data/colquitt"),
    shell:
        """
        mkdir {output}
        wget -O {output}/HVC_RA_RNA_counts.csv 'https://cloud.biohpc.swmed.edu/index.php/s/nLicEtkmjGGmRF8/download?path=%2FHVC_RA&files=HVC_RA_RNA_counts.csv&downloadStartSecret=6cncz7uwc0g'
        """


rule download_tosches:
    output:
        directory("data/tosches"),
    shell:
        """
        mkdir {output}
        wget -P {output} https://public.brain.mpg.de/Laurent/ReptilePallium2018/turtle.neurons.Robj
        """


rule download_yao:
    output:
        directory("data/yao") # TODO: actual files
    shell:
        """
        echo "Download dataset 6/6: Yao et al. (mouse Ctx & Hpc)"
        mkdir {output}
        wget -P {output} https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_10x/metadata.csv
        wget -P {output} https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_10x/expression_matrix.hdf5
        """
