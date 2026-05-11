import argparse
import io
import logging
import os
import sys
import urllib.request
import zipfile

_log = logging.getLogger("download_targetscan_summary_counts_default_predictions")

TARGETSCAN_OUTPUT_NAME = "Summary_Counts.default_predictions.txt"


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
        description="Download and extract the default TargetScan Summary_Counts file.")
    parser.add_argument("--url", required=True, help="TargetScan archive URL.")
    parser.add_argument("--output", required=True, help="Path to the extracted TSV output.")
    parser.add_argument("--log_file", default="download_targetscan_summary_counts.log", help="Path to the log file.")
    parser.add_argument(
        "--log_level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level.",
    )
    args = parser.parse_args()

    set_logging(args.log_file, args.log_level)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    _log.info(f"Downloading TargetScan archive from {args.url}")
    try:
        with urllib.request.urlopen(args.url) as response:
            archive_bytes = response.read()
    except Exception as exc:
        _log.error(f"Error downloading TargetScan archive: {exc}")
        sys.exit(1)

    _log.info("Extracting TargetScan TSV from archive")
    try:
        with zipfile.ZipFile(io.BytesIO(archive_bytes)) as archive:
            member_name = None
            for candidate in archive.namelist():
                if candidate.endswith(TARGETSCAN_OUTPUT_NAME):
                    member_name = candidate
                    break
            if member_name is None:
                raise ValueError(f"Could not find {TARGETSCAN_OUTPUT_NAME} in the TargetScan archive")
            with archive.open(member_name) as source, open(args.output, "wb") as target:
                target.write(source.read())
    except Exception as exc:
        _log.error(f"Error extracting TargetScan TSV: {exc}")
        sys.exit(1)

    _log.info(f"Saved TargetScan TSV to {args.output}")


if __name__ == "__main__":
    snakemake = globals().get("snakemake")
    if snakemake is not None:
        args = [
            "--url", snakemake.params.url,
            "--output", snakemake.output[0],
            "--log_file", snakemake.log[0],
        ]
        if hasattr(snakemake.params, "log_level"):
            args += ["--log_level", snakemake.params.log_level]
        sys.argv[1:] = args
        main()
    else:
        main()