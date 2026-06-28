rule download_uORF_source:
    output: config["UORF_SOURCE"]["excel_path"]
    params: url=config["UORF_SOURCE"]["url"]
    log: os.path.join(LOG_DIR, "download_uORF_source.log")
    shell:
        """
        python ../scripts/download_uORF_source.py {log} {output} > {log} 2>&1
        """


rule map_uORF_to_transcripts:
    input:
        uorf_excel=config["UORF_SOURCE"]["excel_path"],
        gtf=HCT116_ISOFORMS["gtf"],
        hct116_metadata=os.path.join(DATA_DIR, "HCT116_isoforms_metadata.csv")
    output:
        os.path.join(DATA_DIR, "HCT116_isoforms_uORF_transcript_mapping.csv")
    params:
        log_level=LOG_LEVEL
    conda:
        "../envs/HCT116.yaml"
    log:
        os.path.join(LOG_DIR, "map_uORF_to_transcripts.log")
    script:
        "../scripts/map_uORF_to_transcripts.py"
