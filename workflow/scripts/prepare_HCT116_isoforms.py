import pyranges as pr
import pandas as pd
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
import argparse
import logging
import sys

_log = logging.getLogger("prepare_HCT116_isoforms")

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


def load_transcriptome_fasta(fasta_path):
    """Load FASTA and return dict: id -> sequence (str)."""
    fasta_dict = SeqIO.to_dict(SeqIO.parse(open(fasta_path, 'r'), 'fasta'))
    return {str(k): str(v.seq) for k, v in fasta_dict.items()} 


def collect_annotations(df, tx_col='transcript_id', features_of_interest=None):
    """Group GTF annotations by transcript and return exon dictionaries and CDS bounds.

    Args:
        df (pd.DataFrame): GTF annotations dataframe.
        tx_col (str): Name of the transcript id column.
        features_of_interest (list): List of GTF features to include (e.g. ['exon', 'cds']).

    Returns: cds_by_tx, missing_tx_rows_count
    """
    if features_of_interest is None:
        features_of_interest = ['exon', 'cds', 'transcript']

    if tx_col not in df.columns:
        _log.warning(f"Transcript id column '{tx_col}' not found in GTF dataframe. Expected it to be present.")

    df_features = df[df['Feature'].str.lower().isin(features_of_interest)].copy()
    missing_tx_mask = df_features[tx_col].isna()
    missing_tx_count = int(missing_tx_mask.sum())
    if missing_tx_count:
        _log.warning(f"{missing_tx_count} annotation rows are missing '{tx_col}' and will be ignored.")
    df_features = df_features[~missing_tx_mask].copy()

    cds_by_tx = {}

    for tx_id, grp in df_features.groupby(tx_col):
        strand_vals = grp['Strand'].dropna().unique()
        # if multiple strands, warn and pick the first
        if len(strand_vals) > 1:
            _log.Warning(f"Transcript {tx_id} has multiple strands {strand_vals}; using the first.")
        strand = strand_vals[0] if len(strand_vals) > 0 else '+'

        # collect and sort exons by genomic start
        ex_grp = grp[grp['Feature'].str.lower() == 'exon']
        exons = []
        if not ex_grp.empty:
            for s, e in ex_grp[['Start', 'End']].itertuples(index=False, name=None):
                exons.append({'start': int(s), 'end': int(e)})
        exons = sorted(exons, key=lambda x: x['start'])

        # collect and sort CDS by genomic start
        cds_grp = grp[grp['Feature'].str.lower() == 'cds']
        cdss = []
        if not cds_grp.empty:
            for s, e in cds_grp[['Start', 'End']].itertuples(index=False, name=None):
                cdss.append({'start': int(s), 'end': int(e)})
        cdss = sorted(cdss, key=lambda x: x['start'])

        # compute transcript offsets (tx_start, tx_end) for each exon in transcript order
        tpos = 0
        for exon in exons:
            gs, ge = exon['start'], exon['end']
            exon_len = ge - gs + 1
            exon['tx_start'] = tpos
            exon['tx_end'] = tpos + exon_len  # exclusive
            tpos += exon_len

        # for negative strand, reverse exon ordering and adapt tx_start/tx_end
        if strand == '-':
            exons = exons[::-1]
            tlen = tpos
            for exon in exons:
                old_start, old_end = exon['tx_start'], exon['tx_end']
                exon['tx_start'] = tlen - old_end
                exon['tx_end'] = tlen - old_start

        # iterate over CDSs, find exon their belong to (one exon should contain both cds start and end), and compute tx_start/tx_end for CDS
        cds_list = []
        for cds in cdss:
            cds_gstart, cds_gend = cds['start'], cds['end']
            # find exon that contains both cds_gstart and cds_gend
            containing_exon = None
            for exon in exons:
                if cds_gstart >= exon['start'] and cds_gend <= exon['end']:
                    containing_exon = exon
                    break
            if containing_exon is None:
                _log.warning(f"CDS {cds_gstart}-{cds_gend} for transcript {tx_id} does not fully fit within a single exon; skipping this CDS.")
                continue
            # compute transcript coordinates for both genomic endpoints in a strand-aware way.
            # After the exon normalization above, containing_exon['tx_start'] and ['tx_end']
            # already describe the exon span in transcript coordinates (0-based, end-exclusive),
            # while containing_exon['start']/['end'] are genomic coordinates (low->high).
            # For positive strand exons the offset within exon is (genomic_pos - exon.start).
            # For negative strand exons transcript coordinate increases as genomic coordinate decreases,
            # so the offset within exon is (exon.end - genomic_pos).
            if strand == '+':
                tx_a = containing_exon['tx_start'] + (cds_gstart - containing_exon['start'])
                tx_b = containing_exon['tx_start'] + (cds_gend - containing_exon['start'])
            else:
                # negative strand: map genomic positions to transcript coords using exon end as origin
                tx_a = containing_exon['tx_start'] + (containing_exon['end'] - cds_gstart)
                tx_b = containing_exon['tx_start'] + (containing_exon['end'] - cds_gend)

            # ensure ordering: transcript coordinates should have start <= end; use min/max
            cds_tx_start = int(min(tx_a, tx_b))
            cds_tx_end = int(max(tx_a, tx_b)) + 1  # make end exclusive
            cds_list.append({'start': cds_gstart, 'end': cds_gend, 'tx_start': cds_tx_start, 'tx_end': cds_tx_end})
        cds_by_tx[tx_id] = cds_list

        # compute overall CDS bounds for transcript (min start, max end)
        if cds_list:
            cds_start = min(cds['tx_start'] for cds in cds_list)
            cds_end = max(cds['tx_end'] for cds in cds_list)
            cds_by_tx[tx_id] = (cds_start, cds_end)
        else:
            cds_by_tx[tx_id] = (None, None)

    return cds_by_tx, missing_tx_count


def build_transcript_row(tx_id, sequences, cds_by_tx):
    """Build output dict for one transcript using transcript-level CDS bounds.

    We map cds_genomic_start/end to transcript indices using exon tx_start offsets. Then we extract
    the CDS sequence as tx_seq[cds_tx_start:cds_tx_end_exclusive]. If no CDS (None bounds), we set
    utr5 to empty and utr3 to the full transcript sequence (per your request).
    """
    tx_seq = sequences.get(str(tx_id))
    if tx_seq is None:
        return None

    cds_start, cds_end = cds_by_tx.get(tx_id, (None, None))

    if cds_start is None or cds_end is None:
        # no CDS annotated: whole transcript becomes utr3
        return {
            'transcript_id': tx_id,
            'tx_sequence': tx_seq,
            'utr3_sequence': tx_seq,
            'cds_sequence': '',
            'utr5_sequence': '',
            'tx_length': len(tx_seq),
            'utr3_length': len(tx_seq),
            'cds_length': 0,
            'utr5_length': 0
        }

    # convert inclusive end to exclusive end for python slicing
    cds_end = cds_end + 1

    cds_seq = tx_seq[cds_start:cds_end]
    utr5_seq = tx_seq[:cds_start]
    utr3_seq = tx_seq[cds_end:]

    return {
        'transcript_id': tx_id,
        'tx_sequence': tx_seq,
        'utr3_sequence': utr3_seq,
        'cds_sequence': cds_seq,
        'utr5_sequence': utr5_seq,
        'tx_length': len(tx_seq),
        'utr3_length': len(utr3_seq),
        'cds_length': len(cds_seq),
        'utr5_length': len(utr5_seq)
    }


def main():
    parser = argparse.ArgumentParser(
        description="Enrich GTF with transcript coordinates and add utr5/utr3 features")
    parser.add_argument("--gtf", required=True, help="Input GTF file")
    parser.add_argument("--fasta", required=True, help="Input transcriptome FASTA file")
    parser.add_argument("--out_tsv", required=True, help="Output TSV file with transcript sequences ")
    parser.add_argument("--log_file", required=False, help="Log file path", default="prepare_HCT116_isoforms.log")
    parser.add_argument("--log_level", required=False, default="INFO")
    parser.add_argument("--emit-utrs-without-cds", dest="emit_utrs_without_cds", action="store_true",
                        help="Emit UTR feature(s) for transcripts that have no CDS (creates a single 'utr' covering the transcript)")
    args = parser.parse_args()

    set_logging(args.log_file, getattr(logging, args.log_level.upper(), logging.INFO))

    _log.info(f"Reading GTF from {args.gtf}")
    gtf = pr.read_gtf(args.gtf).as_df()
    _log.info(f"Loaded {len(gtf)} GTF annotations")

    _log.info("Collecting transcript sequences from FASTA")
    sequences = load_transcriptome_fasta(args.fasta)
    _log.info(f"Loaded {len(sequences)} transcript sequences from FASTA")

    _log.info("Collecting CDS bounds from GTF annotations")
    cds_bounds_by_tx, missing_tx_rows = collect_annotations(gtf, tx_col='transcript_id')
    _log.info(f"Collected CDS bounds for {len(cds_bounds_by_tx)} transcripts")

    if missing_tx_rows > 0:
        _log.warning(f"{missing_tx_rows} annotation rows are missing 'transcript_id' and were ignored.")

    tx_in_gtf = set(cds_bounds_by_tx.keys())
    tx_in_fasta = set(sequences.keys())
    transcripts_to_process = sorted(tx_in_gtf & tx_in_fasta)
    missing_in_fasta = sorted(tx_in_gtf - tx_in_fasta)
    if missing_in_fasta:
        _log.warning(f"Transcripts with exon annotations missing in FASTA: {len(missing_in_fasta)}")

    transcript_meta = gtf[gtf['Feature'] == 'transcript'].set_index('transcript_id')[[
        'GENCODE_gene_id', 'GENCODE_transcript_id', 'gene_name'
    ]]
    rows = []
    for tx_id in transcripts_to_process:
        row = build_transcript_row(tx_id, sequences, cds_bounds_by_tx)
        if row is None:
            # sequence missing in FASTA or cannot build
            continue
        # Attach metadata if available (use .loc only when present to avoid KeyError)
        if tx_id in transcript_meta.index:
            meta = transcript_meta.loc[tx_id]
            row['GENCODE_gene_id'] = meta.get('GENCODE_gene_id', None)
            row['GENCODE_transcript_id'] = meta.get('GENCODE_transcript_id', None)
            row['gene_name'] = meta.get('gene_name', None)
        else:
            row['GENCODE_gene_id'] = None
            row['GENCODE_transcript_id'] = None
            row['gene_name'] = None
        rows.append(row)

    out_df = pd.DataFrame(
        rows, 
        columns=[
            'transcript_id', 
            'tx_sequence', 
            'utr3_sequence', 
            'cds_sequence', 
            'utr5_sequence', 
            'tx_length', 
            'utr3_length', 
            'cds_length', 
            'utr5_length', 
            'GENCODE_gene_id', 
            'GENCODE_transcript_id', 
            'gene_name'
        ]
    )

    _log.info(f"Writing output TSV to {args.out_tsv} with {len(out_df)} transcripts")
    out_df.to_csv(args.out_tsv, index=False)


if __name__ == "__main__":
    # handle both script and snakemake execution
    if 'snakemake' in globals():
        # snakemake execution
        from snakemake.script import snakemake
        args = [
            "--gtf", snakemake.input.gtf,
            "--fasta", snakemake.input.fasta,
            "--out_tsv", snakemake.output[0],
            "--log_file", snakemake.log[0]
        ]
        if hasattr(snakemake.params, 'log_level'):
            args += ["--log_level", snakemake.params.log_level]
        if hasattr(snakemake.params, 'emit_utrs_without_cds') and snakemake.params.emit_utrs_without_cds:
            args.append("--emit-utrs-without-cds")
        sys.argv[1:] = args
        main()
    else:
        # command line execution
        main()