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


rule all:
    input:
        rules.download_all.input,
        rules.preprocess_all.input,
