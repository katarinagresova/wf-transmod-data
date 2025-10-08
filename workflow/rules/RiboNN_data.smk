rule prepare_RiboNN_data:
    output:
        os.path.join(DATA_DIR, "{dataset}.csv")
    params:
        url=lambda wildcards: RIBONN_DATASETS[wildcards.dataset]['url'],
        sheet_name=lambda wildcards: RIBONN_DATASETS[wildcards.dataset]['sheet_name'],
        add_seqs=lambda wildcards: RIBONN_DATASETS[wildcards.dataset]['add_seqs']
    conda:
        "../envs/RiboNN.yaml"
    log:
        os.path.join(LOG_DIR, "prepare_RiboNN_data_{dataset}.log")
    script:
        "../scripts/prepare_RiboNN_data.py"