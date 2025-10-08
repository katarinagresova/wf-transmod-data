rule prepare_RiboNN_data:
    output:
        os.path.join(DATA_DIR, "{dataset}.csv")
    params:
        url=lambda wildcards: RIBONN_DATASETS[wildcards.dataset]['url'],
        sheet_name=lambda wildcards: RIBONN_DATASETS[wildcards.dataset]['sheet_name']
    conda:
        "../envs/RiboNN.yaml"
    log:
        os.path.join(LOG_DIR, "prepare_RiboNN_data_{dataset}.log")
    shell:
        """
        python {workflow.basedir}/scripts/prepare_RiboNN_data.py \
            --url {params.url} \
            --sheet_name {params.sheet_name} \
            --output {output} \
            --log_file {log} \
            --log_level DEBUG
        """