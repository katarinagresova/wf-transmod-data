#!/usr/bin/env python3
"""Map uORFs to HCT116 transcripts using 5' UTR exon coordinates."""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import pyranges as pr

_log = logging.getLogger("map_uORF_to_transcripts")


def set_logging(log_file, log_level):
    _log.setLevel(log_level)
    fh = logging.FileHandler(log_file)
    ch = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    _log.addHandler(ch)
    _log.addHandler(fh)


def gtf_to_pr_5utr(gtf_path):
    """Parse GTF and identify 5' UTR exons per transcript."""
    df = pd.read_csv(
        gtf_path, sep='\t', comment='#', header=None,
        names=['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score',
               'Strand', 'Frame', 'Attributes']
    )

    def get_attr(attr_str, key):
        val = attr_str.split(f'{key} "')[1].split('"')[0] if f'{key} "' in attr_str else None
        return val

    df['transcript_id'] = df['Attributes'].apply(lambda x: get_attr(x, 'transcript_id'))

    has_five_prime_utr = (df['Feature'] == 'five_prime_utr').any()

    if has_five_prime_utr:
        _log.info("Using annotated five_prime_utr features from GTF")
        utr5_df = df[df['Feature'] == 'five_prime_utr'][
            ['Chromosome', 'Start', 'End', 'Strand', 'transcript_id']
        ].copy()
        utr5_df['Start'] -= 1
    else:
        _log.info("No five_prime_utr annotations found; computing from exon + CDS")
        exons = df[df['Feature'] == 'exon'].copy()
        cdss = df[df['Feature'] == 'CDS'].copy()

        if cdss.empty:
            _log.warning("No CDS features found in GTF; returning empty PyRanges")
            return pr.PyRanges(pd.DataFrame(columns=[
                'Chromosome', 'Start', 'End', 'Strand', 'transcript_id'
            ]))

        cds_bounds = cdss.groupby('transcript_id').agg({
            'Start': 'min',
            'End': 'max',
            'Strand': 'first'
        }).reset_index()

        def get_cds_start(row):
            return row['Start'] if row['Strand'] == '+' else row['End']

        cds_bounds['cds_start'] = cds_bounds.apply(get_cds_start, axis=1)

        exons_with_cds = exons.merge(cds_bounds[['transcript_id', 'cds_start']], on='transcript_id')

        utr5_parts = []
        for _, row in exons_with_cds.iterrows():
            exon_start, exon_end = row['Start'], row['End']
            cds_start = row['cds_start']
            strand = row['Strand']

            if strand == '+':
                if exon_end <= cds_start:
                    utr5_parts.append({
                        'Chromosome': row['Chromosome'],
                        'Start': exon_start,
                        'End': exon_end,
                        'Strand': strand,
                        'transcript_id': row['transcript_id']
                    })
                elif exon_start < cds_start:
                    utr5_parts.append({
                        'Chromosome': row['Chromosome'],
                        'Start': exon_start,
                        'End': cds_start,
                        'Strand': strand,
                        'transcript_id': row['transcript_id']
                    })
            else:
                if exon_start >= cds_start:
                    utr5_parts.append({
                        'Chromosome': row['Chromosome'],
                        'Start': exon_start,
                        'End': exon_end,
                        'Strand': strand,
                        'transcript_id': row['transcript_id']
                    })
                elif exon_end > cds_start:
                    utr5_parts.append({
                        'Chromosome': row['Chromosome'],
                        'Start': cds_start,
                        'End': exon_end,
                        'Strand': strand,
                        'transcript_id': row['transcript_id']
                    })

        if not utr5_parts:
            _log.warning("No 5' UTR exons found; returning empty PyRanges")
            return pr.PyRanges(pd.DataFrame(columns=[
                'Chromosome', 'Start', 'End', 'Strand', 'transcript_id'
            ]))

        utr5_exons = pd.DataFrame(utr5_parts)
        utr5_exons['Start'] -= 1
        utr5_df = utr5_exons

    utr5_df = utr5_df[utr5_df['transcript_id'].notna()].reset_index(drop=True)
    _log.info(f"Found {len(utr5_df)} 5' UTR exon regions across {utr5_df['transcript_id'].nunique()} transcripts")

    return pr.PyRanges(utr5_df)


def gtf_to_pr_exons(gtf_path):
    """Parse GTF and extract all exon coordinates per transcript."""
    df = pd.read_csv(
        gtf_path, sep='\t', comment='#', header=None,
        names=['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score',
               'Strand', 'Frame', 'Attributes']
    )

    def get_attr(attr_str, key):
        val = attr_str.split(f'{key} "')[1].split('"')[0] if f'{key} "' in attr_str else None
        return val

    df['transcript_id'] = df['Attributes'].apply(lambda x: get_attr(x, 'transcript_id'))

    exons = df[df['Feature'] == 'exon'][
        ['Chromosome', 'Start', 'End', 'Strand', 'transcript_id']
    ].copy()
    exons['Start'] -= 1

    return pr.PyRanges(exons)


def main():
    parser = argparse.ArgumentParser(
        description="Map uORFs to HCT116 transcripts using 5' UTR exon coordinates.")
    parser.add_argument("--uorf_excel", "-i", required=True,
                        help="Path to input Excel file with uORF data.")
    parser.add_argument("--gtf", "-g", required=True,
                        help="Path to GTF file.")
    parser.add_argument("--output", "-o", required=True,
                        help="Path to output CSV file.")
    parser.add_argument("--log_file", "-l", default="map_uORF_to_transcripts.log",
                        help="Path to the log file.")
    parser.add_argument("--log_level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging level.")
    args = parser.parse_args()

    set_logging(args.log_file, args.log_level)

    excel_path = Path(args.uorf_excel)
    gtf_path = Path(args.gtf)
    output_csv = Path(args.output)

    # Load uORF data
    _log.info(f"Loading uORF data from {excel_path}")
    uORFs = pd.read_excel(excel_path, sheet_name="S2. PHASE I Ribo-seq ORFs")
    _log.info(f"Loaded {len(uORFs)} uORFs")

    # Filter out lncRNA and processed_transcript
    before = len(uORFs)
    uORFs = uORFs[(uORFs['orf_biotype'] != 'lncRNA') &
                  (uORFs['orf_biotype'] != 'processed_transcript')]
    uORFs.reset_index(drop=True, inplace=True)
    _log.info(f"Filtered to {len(uORFs)} uORFs (removed {before - len(uORFs)} lncRNAs/processed_transcripts)")

    # Load GTF-derived 5' UTR exons
    _log.info(f"Parsing GTF: {gtf_path}")
    hct116_exons = gtf_to_pr_exons(gtf_path)
    _log.info(f"Loaded {len(hct116_exons)} exons for {hct116_exons.transcript_id.nunique()} transcripts")

    hct116_exons_5utr = gtf_to_pr_5utr(gtf_path)
    _log.info(f"Loaded {len(hct116_exons_5utr)} 5' UTR exons for {hct116_exons_5utr.transcript_id.nunique()} transcripts")

    # Build uORF PyRanges
    _log.info("Building uORF PyRanges")
    rows = []
    for _, row in uORFs.iterrows():
        starts = [int(s) for s in str(row['starts']).split(';')]
        ends = [int(e) for e in str(row['ends']).split(';')]
        rows.append({
            'orf_name': row['orf_name'],
            'Chromosome': f"chr{row['chrm']}",
            'Start': starts[0],
            'End': ends[0],
            'Strand': row['strand'],
            'uORF_gene': row['gene_name'],
            'orf_biotype': row['orf_biotype'],
        })

    uORFs_pr = pr.PyRanges(pd.DataFrame(rows))
    _log.info(f"Created {len(uORFs_pr)} uORF entries")

    # Find overlaps
    _log.info("Finding overlaps between uORFs and 5' UTR exons")
    overlaps = uORFs_pr.join(hct116_exons, how="containment", suffix="_tx")
    overlaps_5utr = overlaps.join(hct116_exons_5utr, suffix="_tx_5utr")

    # Convert to DataFrame
    mapping_df = overlaps_5utr.df[[
        'orf_name', 'Chromosome', 'Start', 'End', 'Strand',
        'orf_biotype', 'transcript_id'
    ]]
    mapping_df.columns = [
        'orf_name', 'uORF_chrom', 'uORF_start', 'uORF_end', 'uORF_strand',
        'uORF_biotype', 'transcript_id'
    ]
    mapping_df = mapping_df[[
        'transcript_id', 'orf_name', 'uORF_chrom', 'uORF_start', 'uORF_end',
        'uORF_strand', 'uORF_biotype'
    ]]

    _log.info(f"\n=== Mapping Results ===")
    _log.info(f"Total uORF-transcript overlaps: {len(mapping_df)}")
    _log.info(f"Unique uORFs mapped: {mapping_df['orf_name'].nunique()}")
    _log.info(f"Unique transcripts with uORFs: {mapping_df['transcript_id'].nunique()}")

    # Save output
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    mapping_df.to_csv(output_csv, index=False)
    _log.info(f"Saved mapping to {output_csv}")


if __name__ == "__main__":
    if 'snakemake' in globals():
        from snakemake.script import snakemake
        args = [
            "--uorf_excel", snakemake.input.uorf_excel,
            "--gtf", snakemake.input.gtf,
            "--output", snakemake.output[0],
            "--log_file", snakemake.log[0]
        ]
        if hasattr(snakemake.params, 'log_level'):
            args += ["--log_level", snakemake.params.log_level]
        sys.argv[1:] = args
        main()
    else:
        main()
