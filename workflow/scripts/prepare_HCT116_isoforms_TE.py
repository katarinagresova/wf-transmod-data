import pandas as pd
import numpy as np
import re
import argparse
import logging
import sys

_log = logging.getLogger("prepare_HCT116_isoforms_TE")

def set_logging(log_file, log_level):
    _log.setLevel(log_level)
    # create file handler that logs debug and higher level messages
    fh = logging.FileHandler(log_file)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    # create formatter and add it to the handlers
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    # add the handlers to logger
    _log.addHandler(ch)
    _log.addHandler(fh)

def apply_pair_ribo_tpm_filter(df, condition, timepoint, threshold, logger):
    pair = f"{condition}_{timepoint}"
    ribo_tpm_cols = [
        col for col in df.columns
        if col.startswith('ribo_tpm_') and col.endswith(f'_{pair}')
    ]
    if not ribo_tpm_cols:
        logger.warning(f"No ribo_tpm columns found for pair '{pair}'; skipping filter")
        return df

    pair_cols = [col for col in df.columns if col.endswith(f'_{pair}')]
    if not pair_cols:
        logger.warning(f"No measurement columns found for pair '{pair}'; skipping filter")
        return df

    max_ribo_tpm = df[ribo_tpm_cols].max(axis=1)
    keep_mask = max_ribo_tpm.notna() & (max_ribo_tpm >= threshold)
    filtered_pairs = int((~keep_mask).sum())
    df = df.loc[keep_mask].copy()

    logger.info(
        f"Applied ribo_tpm filter for pair '{pair}' (threshold={threshold}): "
        f"{filtered_pairs} transcript rows removed"
    )
    return df


def select_measurement_columns(df):
    prefixes = ('log2TE', 'log2FC', 'padj', 'ribo_tpm', 'rna_tpm')
    return [col for col in df.columns if col.startswith(prefixes)]


def build_measurement_frame(logte_df, condition, timepoint):
    pair = f"{condition}_{timepoint}"
    pair_df = logte_df[logte_df['factor_time'] == pair].copy()
    if pair_df.empty:
        _log.warning(f"No rows found in logTE for pair '{pair}'")

    measurement_cols = select_measurement_columns(pair_df)
    if not measurement_cols:
        return pd.DataFrame(columns=['tx_id', 'gene_name'])

    keep_cols = ['Name', 'gene_name'] + measurement_cols
    pair_df = pair_df[keep_cols].copy()
    
    # Verify that Name is unique
    duplicate_names = pair_df[pair_df.duplicated(subset=['Name'], keep=False)]
    if not duplicate_names.empty:
        dup_list = duplicate_names['Name'].unique().tolist()
        raise ValueError(
            f"Found {len(duplicate_names)} rows with duplicate 'Name' values in pair '{pair}'. "
            f"'Name' should be unique. Duplicates: {dup_list}"
        )

    rename_map = {col: f"{col}_{pair}" for col in measurement_cols}
    pair_df.rename(columns={'Name': 'tx_id', **rename_map}, inplace=True)
    pair_df = pair_df.sort_values('tx_id').reset_index(drop=True)
    return pair_df


def main():
    parser = argparse.ArgumentParser(
        description="Extract measurement values for one condition and one timepoint from logTE")
    parser.add_argument("--logTE", required=True, help="Input logTE file")
    parser.add_argument("--condition", required=True, help="Condition to filter on (e.g. eif3d)")
    parser.add_argument("--timepoint", required=True, help="Timepoint to filter on (e.g. 8h)")
    parser.add_argument("--out_tsv", required=True, help="Output CSV file with measurements")
    parser.add_argument("--log_file", default="prepare_HCT116_isoforms_TE.log", help="Log file path")
    parser.add_argument(
        "--log_level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level."
    )
    parser.add_argument("--ribo-tpm-threshold", dest="ribo_tpm_threshold", type=float, default=5,
                        help="Minimum ribo_tpm threshold for per-condition filtering (default: 5)")
    args = parser.parse_args()

    set_logging(args.log_file, getattr(logging, args.log_level.upper(), logging.INFO))

    sep = '\t' if args.logTE.endswith('.tsv') else ','
    log2TE_df = pd.read_csv(args.logTE, index_col=0, sep=sep)
    _log.info(f"Loaded log2TE data with {len(log2TE_df)} rows from {args.logTE}")

    if 'Name' not in log2TE_df.columns or 'factor_time' not in log2TE_df.columns:
        raise ValueError("logTE file must contain 'Name' and 'factor_time' columns")

    out_df = build_measurement_frame(log2TE_df, args.condition, args.timepoint)
    if args.ribo_tpm_threshold > 0 and not out_df.empty:
        out_df = apply_pair_ribo_tpm_filter(
            out_df,
            args.condition,
            args.timepoint,
            args.ribo_tpm_threshold,
            _log,
        )

    _log.info(f"Writing output CSV to {args.out_tsv} with {len(out_df)} rows")
    out_df.to_csv(args.out_tsv, index=False)


if __name__ == "__main__":
    if 'snakemake' in globals():
        args = [
            "--logTE", snakemake.input.logTE,
            "--condition", snakemake.wildcards.condition,
            "--timepoint", snakemake.wildcards.timepoint,
            "--out_tsv", snakemake.output[0],
            "--log_file", snakemake.log[0]
        ]
        if hasattr(snakemake.params, 'log_level'):
            args += ["--log_level", snakemake.params.log_level]
        if hasattr(snakemake.params, 'ribo_tpm_threshold'):
            args += ["--ribo-tpm-threshold", str(snakemake.params.ribo_tpm_threshold)]
        sys.argv[1:] = args
        main()
    else:
        main()