#!/usr/bin/env python3
"""Download uORF source data from Springer."""

import argparse
import logging
import sys
from pathlib import Path
import urllib.request

_log = logging.getLogger("download_uORF_source")


def set_logging(log_file, log_level):
    _log.setLevel(log_level)
    fh = logging.FileHandler(log_file)
    ch = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    _log.addHandler(ch)
    _log.addHandler(fh)


def main():
    parser = argparse.ArgumentParser(
        description="Download uORF source data from Springer.")
    parser.add_argument("--url", "-i", required=True,
                        help="URL of the input file to download.")
    parser.add_argument("--output", "-o", required=True,
                        help="Path to the output file.")
    parser.add_argument("--log_file", "-l", default="download_uORF_source.log",
                        help="Path to the log file.")
    parser.add_argument("--log_level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging level.")
    args = parser.parse_args()

    set_logging(args.log_file, args.log_level)

    output = Path(args.output)

    output.parent.mkdir(parents=True, exist_ok=True)

    _log.info(f"Downloading uORF source data from {args.url}")
    urllib.request.urlretrieve(args.url, output)
    _log.info(f"Saved to {output}")


if __name__ == "__main__":
    if 'snakemake' in globals():
        from snakemake.script import snakemake
        args = [
            "--url", snakemake.params.url,
            "--output", snakemake.output[0],
            "--log_file", snakemake.log[0]
        ]
        if hasattr(snakemake.params, 'log_level'):
            args += ["--log_level", snakemake.params.log_level]
        sys.argv[1:] = args
        main()
    else:
        main()
