import pandas as pd
import requests
import argparse
import logging
import os
import sys
import snakemake

_log = logging.getLogger("prepare_RiboNN_data")

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

def main():
    parser = argparse.ArgumentParser(
        description="Prepare RiboNN data from Excel files.")
    parser.add_argument("--url", "-i", required=True,
                        help="URL of the input Excel file.")
    parser.add_argument("--sheet_name", "-s", required=True,
                        help="Sheet name in the Excel file.")
    parser.add_argument("--output", "-o", required=True,
                        help="Path to the output CSV file.")
    parser.add_argument("--log_file", "-l", default="prepare_RiboNN_data.log",
                        help="Path to the log file.")
    parser.add_argument("--log_level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging level.")
    args = parser.parse_args()

    set_logging(args.log_file, args.log_level)

    # download the excel file from the url
    file_path = args.url.split("/")[-1]
    _log.info(f"Downloading data from {args.url} to {file_path}")
    try:
        response = requests.get(args.url)
        response.raise_for_status()  # Raise an error for bad status codes
        with open(file_path, 'wb') as f:
            f.write(response.content)
        _log.info(f"Successfully downloaded file to {file_path}")
    except Exception as e:
        _log.error(f"Error downloading file: {e}")
        return


    _log.info(f"Reading data from {file_path}, sheet: {args.sheet_name}")
    try:
        df = pd.read_excel(file_path, sheet_name=args.sheet_name, engine='openpyxl')
        _log.info(f"Successfully read {len(df)} rows.")
    except Exception as e:
        _log.error(f"Error reading Excel file: {e}")
        return
    
    # delete downloaded excel file
    try:
        os.remove(file_path)
        _log.info(f"Deleted temporary file {file_path}.")
    except Exception as e:
        _log.warning(f"Could not delete temporary file {file_path}: {e}")

    # remove TE prefix from column names
    df.columns = [col.replace("TE_", "") for col in df.columns]
    _log.info("Removed 'TE_' prefix from column names.")

    # rename gene symbol column: SYMBOL -> gene_symbol
    if "SYMBOL" in df.columns:
        df = df.rename(columns={"SYMBOL": "gene_symbol"})
        _log.info("Renamed 'SYMBOL' column to 'gene_symbol'.")
    else:
        _log.warning("'SYMBOL' column not found in the data.")

    # chceck if there are any missing values in utr5_size, cds_size or utr3_size columns
    for col in ["utr5_size", "cds_size", "utr3_size"]:
        if col in df.columns:
            missing_count = df[col].isna().sum()
            if missing_count > 0:
                _log.warning(f"Column '{col}' has {missing_count} missing values.")
        else:
            _log.warning(f"Column '{col}' not found in the data.")

    # save dataframe to csv
    try:
        df.to_csv(args.output, index=False)
        _log.info(f"Saved processed data to {args.output}.")
    except Exception as e:
        _log.error(f"Error saving CSV file: {e}")

if __name__ == "__main__":
    # handle both script and snakemake execution
    if 'snakemake' in sys.modules:
        # snakemake execution
        args = [
            "--url", snakemake.params.url,
            "--sheet_name", snakemake.params.sheet_name,
            "--output", snakemake.output[0],
            "--log_file", snakemake.log[0]
        ]
        if hasattr(snakemake.params, 'log_level'):
            args += ["--log_level", snakemake.params.log_level]
        sys.argv[1:] = args
        main()
    else:
        # script execution
        main()