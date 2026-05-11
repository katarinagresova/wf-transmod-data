def get_ribonn_param(wildcards, key):
    if wildcards.dataset not in RIBONN_DATASETS:
        available_datasets = ", ".join(sorted(RIBONN_DATASETS))
        raise ValueError(
            f"Unknown RiboNN dataset '{wildcards.dataset}'. Available datasets: {available_datasets}"
        )
    return RIBONN_DATASETS[wildcards.dataset][key]


rule prepare_RiboNN_data:
    output:
        os.path.join(DATA_DIR, "RiboNN_{dataset}.csv")
    params:
        url=lambda wildcards: get_ribonn_param(wildcards, 'url'),
        sheet_name=lambda wildcards: get_ribonn_param(wildcards, 'sheet_name'),
        add_seqs=lambda wildcards: get_ribonn_param(wildcards, 'add_seqs')
    conda:
        "../envs/RiboNN.yaml"
    log:
        os.path.join(LOG_DIR, "prepare_RiboNN_data_{dataset}.log")
    script:
        "../scripts/prepare_RiboNN_data.py"