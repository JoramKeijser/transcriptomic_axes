onstart:
    print("### Complete workflow ###")
    shell("mkdir -p log/jobs")


module download:
    snakefile:
        "download.smk"


use rule * from download as download_*


module preprocess:
    snakefile:
        "preprocess.smk"


use rule * from preprocess as preprocess_*


module figure1:
    snakefile:
        "fig1.smk"


use rule * from figure1 as figure1_*


module figure2:
    snakefile:
        "fig2.smk"


use rule * from figure2 as figure2_*


module figure3:
    snakefile:
        "fig3.smk"


use rule * from figure3 as figure3_*


module figure4:
    snakefile:
        "fig4.smk"


use rule * from figure4 as figure4_*


module figure5:
    snakefile:
        "fig5.smk"


use rule * from figure5 as figure5_*


rule all:
    input:
        rules.download_all.input,
        rules.preprocess_all.input,
        rules.figure1_all.input,
        rules.figure2_all.input,
        rules.figure3_all.input,
        rules.figure4_all.input,
        rules.figure5_all.input,


ruleorder: figure1_pca > figure2_pca
