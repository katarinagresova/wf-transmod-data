# Workflow Overview

This directory contains the Snakemake workflow used to build the TransMod data tables from the bundled HCT116 transcriptome, the EIF perturbation logTE matrix, the TargetScan summary counts archive, and the RiboNN supplemental spreadsheets.

The workflow produces three main data products:

1. HCT116 isoform metadata derived from the provided GTF/FASTA pair.
2. HCT116 translation efficiency tables for each condition and timepoint in the EIF dataset.
3. miRNA binding score matrices derived from TargetScan Summary_Counts.default_predictions.

## Description

The workflow is driven by [workflow/Snakefile](workflow/Snakefile) and organized into three rule groups:

- `prepare_HCT116_isoforms`: builds transcript metadata from the HCT116 transcriptome annotation.
- `prepare_HCT116_isoforms_TE`: extracts TE measurements for each condition/timepoint pair from the EIF logTE file.
- `prepare_RiboNN_data`: downloads and reshapes the RiboNN Human and Mouse supplemental Excel sheets.
- `download_targetscan_summary_counts` and `prepare_miRNA_binding_scores`: download TargetScan summary counts and pivot the selected miRNA score types into transcript-level matrices.

## How To Run

From the repository root, run either the wrapper script:

```bash
bash snakemake.sh
```

or invoke Snakemake directly:

```bash
snakemake --use-conda -s workflow/Snakefile --cores 1
```

The bundled wrapper activates the `snake` conda environment before calling Snakemake. The workflow itself uses the per-rule environments under [workflow/envs](workflow/envs).

## Options

Runtime behavior is controlled by [config/config.yaml](config/config.yaml):

- `DATA_DIR`: output directory for generated CSV files, currently `data`.
- `LOG_DIR`: directory for Snakemake script logs, currently `logs`.
- `LOG_LEVEL`: logging verbosity passed to the scripts.
- `RIBONN_DATASETS`: RiboNN dataset definitions, including download URL, sheet name, and whether to add sequence columns.
- `HCT116_ISOFORMS`: input FASTA/GTF paths and transcript filtering options.
- `HCT116_ISOFORMS_TE`: EIF logTE input plus the condition and timepoint values expanded by the workflow.
- `MIRNA_TARGETS`: TargetScan download URL, species ID, and the output score column mapping used for the miRNA matrices.

Useful workflow-specific switches are embedded in the config rather than exposed as CLI flags:

- `HCT116_ISOFORMS.emit_utrs_without_cds`: keep transcripts without CDS annotation.
- `HCT116_ISOFORMS.filter_protein_coding`: restrict isoform metadata to protein-coding transcripts.
- `HCT116_ISOFORMS_TE.ribo_tpm_threshold`: minimum `ribo_tpm_minus` threshold used to filter TE rows per condition/timepoint.

## Outputs

The workflow writes final tables to [data](data) and log files to [logs](logs).

Generated data files include:

- `data/RiboNN_Human_TE.csv`
- `data/RiboNN_Mouse_TE.csv`
- `data/HCT116_isoforms_metadata.csv`
- `data/HCT116_isoforms_TE_{condition}_{timepoint}.csv`
- `data/HCT116_isoforms_miRNA_binding_scores_conserved_sites.csv`
- `data/HCT116_isoforms_miRNA_binding_scores_context_score.csv`
- `data/HCT116_isoforms_miRNA_binding_scores_weighted_context_score.csv`

Supporting files produced by intermediate rules include:

- `resources/Summary_Counts.default_predictions.txt`

Logs are written alongside the pipeline outputs in `logs/`, one file per rule invocation.

## Notes

- The TargetScan downloader expects internet access on first run.
- The RiboNN step also downloads the Excel source files at runtime and removes the temporary download after reading it.
- The workflow assumes the input resources shipped with the repository are present at the paths listed in `config/config.yaml`.