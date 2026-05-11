import argparse
import logging
import sys

import pandas as pd

_log = logging.getLogger("prepare_targetscan_summary_counts_default_predictions_scores")

TARGETSCAN_GENE_SYMBOL_COL = "Gene Symbol"
TARGETSCAN_TRANSCRIPT_ID_COL = "Transcript ID"
TARGETSCAN_MIRNA_FAMILY_COL = "miRNA family"
TARGETSCAN_SPECIES_ID_COL = "Species ID"
INPUT_TRANSCRIPT_ID_COL = "tx_id"
INPUT_GENE_SYMBOL_COL = "gene_name"
INPUT_BIOTYPE_COL = "transcript_biotype"


def validate_targetscan_schema(columns, score_column):
    required_columns = [
        TARGETSCAN_GENE_SYMBOL_COL,
        TARGETSCAN_TRANSCRIPT_ID_COL,
        TARGETSCAN_MIRNA_FAMILY_COL,
        TARGETSCAN_SPECIES_ID_COL,
        score_column,
    ]
    for column in required_columns:
        if column not in columns:
            raise ValueError(f"Column '{column}' not found in TargetScan TSV")


def set_logging(log_file, log_level):
    _log.setLevel(log_level)
    fh = logging.FileHandler(log_file)
    ch = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    _log.addHandler(ch)
    _log.addHandler(fh)


def main():
    parser = argparse.ArgumentParser(
        description="Extract one score column from the TargetScan Summary_Counts file and aggregate by family.")
    parser.add_argument("--input_csv", "-i", required=True,
                        help="Input CSV file with transcript_id and gene_symbol columns.")
    parser.add_argument("--targetscan_tsv", "-t", required=True,
                        help="TargetScan Summary_Counts TSV file.")
    parser.add_argument("--output", "-o", required=True,
                        help="Output CSV file.")
    parser.add_argument("--score_column", "-s", default="Total context++ score",
                        help="Score column name from the TargetScan Summary_Counts TSV.")
    parser.add_argument("--species_id", default=9606, type=int,
                        help="TargetScan Species ID to keep before joining (9606 for Human).")
    parser.add_argument("--log_file", "-l", default="prepare_targetscan_summary_counts_scores.log",
                        help="Path to the log file.")
    parser.add_argument("--log_level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging level.")
    args = parser.parse_args()

    set_logging(args.log_file, args.log_level)

    _log.info(f"Reading input CSV from {args.input_csv}")
    try:
        input_df = pd.read_csv(args.input_csv)
        _log.info(f"Loaded {len(input_df)} rows from input CSV")
    except Exception as exc:
        _log.error(f"Error reading input CSV: {exc}")
        sys.exit(1)

    if INPUT_TRANSCRIPT_ID_COL not in input_df.columns:
        _log.error(f"Column '{INPUT_TRANSCRIPT_ID_COL}' not found in input CSV")
        sys.exit(1)
    if INPUT_GENE_SYMBOL_COL not in input_df.columns:
        _log.error(f"Column '{INPUT_GENE_SYMBOL_COL}' not found in input CSV")
        sys.exit(1)
    if INPUT_BIOTYPE_COL not in input_df.columns:
        _log.error(f"Column '{INPUT_BIOTYPE_COL}' not found in input CSV")
        sys.exit(1)

    # filter for protein_coding biotype
    input_df = input_df[input_df[INPUT_BIOTYPE_COL] == 'protein_coding']
    _log.info(f"Retained {len(input_df)} rows after filtering for protein_coding biotype")

    all_tx_ids = input_df[INPUT_TRANSCRIPT_ID_COL].drop_duplicates()
    _log.info(f"Input unique {INPUT_TRANSCRIPT_ID_COL} count: {len(all_tx_ids)}")

    _log.info(f"Reading TargetScan TSV from {args.targetscan_tsv}")
    try:
        targetscan_df = pd.read_csv(args.targetscan_tsv, sep='\t')
        _log.info(f"Loaded {len(targetscan_df)} rows from TargetScan TSV")
    except Exception as exc:
        _log.error(f"Error reading TargetScan TSV: {exc}")
        sys.exit(1)

    try:
        validate_targetscan_schema(targetscan_df.columns, args.score_column)
    except ValueError as exc:
        _log.error(str(exc))
        sys.exit(1)

    _log.info(f"Filtering TargetScan rows for Species ID {args.species_id}")
    targetscan_df = targetscan_df[targetscan_df[TARGETSCAN_SPECIES_ID_COL] == args.species_id]
    _log.info(f"Retained {len(targetscan_df)} rows after species filter")

    _log.info("Joining input CSV on gene_name to TargetScan on Gene Symbol")
    merged_df = input_df.merge(
        targetscan_df[[TARGETSCAN_GENE_SYMBOL_COL, TARGETSCAN_TRANSCRIPT_ID_COL, TARGETSCAN_MIRNA_FAMILY_COL, args.score_column]],
        left_on=INPUT_GENE_SYMBOL_COL,
        right_on=TARGETSCAN_GENE_SYMBOL_COL,
        how="left",
    )
    _log.info(f"After join: {len(merged_df)} rows")

    merged_df[args.score_column] = pd.to_numeric(merged_df[args.score_column], errors='coerce')

    pivot_df = merged_df.pivot_table(
        index=INPUT_TRANSCRIPT_ID_COL,
        columns=TARGETSCAN_MIRNA_FAMILY_COL,
        values=args.score_column,
        fill_value=0,
    )
    pivot_df = pivot_df.reindex(all_tx_ids, fill_value=0).reset_index()

    output_unique_tx = pivot_df[INPUT_TRANSCRIPT_ID_COL].nunique(dropna=False)
    _log.info(f"Output unique {INPUT_TRANSCRIPT_ID_COL} count: {output_unique_tx}")

    _log.info(f"Output shape: {pivot_df.shape}, (should match input CSV)")
    _log.debug(f"Output columns: {', '.join(pivot_df.columns.astype(str))}")

    try:
        pivot_df.to_csv(args.output, index=False)
        _log.info(f"Saved output to {args.output}")
    except Exception as exc:
        _log.error(f"Error writing output CSV: {exc}")
        sys.exit(1)


if __name__ == "__main__":
    snakemake = globals().get("snakemake")
    if snakemake is not None:
        args = [
            "--input_csv", snakemake.input.input_csv,
            "--targetscan_tsv", snakemake.input.targetscan_tsv,
            "--output", snakemake.output[0],
            "--log_file", snakemake.log[0],
        ]
        if hasattr(snakemake.params, "score_column"):
            args += ["--score_column", snakemake.params.score_column]
        if hasattr(snakemake.params, "species_id"):
            args += ["--species_id", str(snakemake.params.species_id)]
        if hasattr(snakemake.params, "log_level"):
            args += ["--log_level", snakemake.params.log_level]
        sys.argv[1:] = args
        main()
    else:
        main()