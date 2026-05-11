rule prepare_HCT116_isoforms:
    input:
        gtf=HCT116_ISOFORMS['gtf'],
        fasta=HCT116_ISOFORMS['fasta']
    output:
        os.path.join(DATA_DIR, "HCT116_isoforms_metadata.csv")
    params:
        emit_utrs_without_cds=HCT116_ISOFORMS['emit_utrs_without_cds'],
        filter_protein_coding=HCT116_ISOFORMS['filter_protein_coding'],
        log_level=LOG_LEVEL
    conda:
        "../envs/HCT116.yaml"
    log:
        os.path.join(LOG_DIR, "prepare_HCT116_isoforms_metadata.log")
    script:
        "../scripts/prepare_HCT116_isoforms_metadata.py"


rule prepare_HCT116_isoforms_TE:
    input:
        logTE=HCT116_ISOFORMS_TE['logTE']
    output:
        os.path.join(DATA_DIR, "HCT116_isoforms_TE_{condition}_{timepoint}.csv")
    params:
        log_level=LOG_LEVEL,
        ribo_tpm_threshold=HCT116_ISOFORMS_TE.get('ribo_tpm_threshold', 5)
    conda:
        "../envs/HCT116.yaml"
    log:
        os.path.join(LOG_DIR, "prepare_HCT116_isoforms_TE_{condition}_{timepoint}.log")
    script:
        "../scripts/prepare_HCT116_isoforms_TE.py"