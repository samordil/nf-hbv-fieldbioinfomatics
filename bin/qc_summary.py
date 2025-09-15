#!/usr/bin/env python3
"""
qc_summary.py

Summarize one or more ARTIC-style QC report files (*.qc.report.tsv)
into a single combined TSV file.

Each input file has the format (4 lines, tab-separated values):

    total-reads:     <int>
    mapped-reads:    <int>
    mapped-percent:  <float>
    coverage-percent:<float>

The script extracts the sample ID from the filename prefix
(the part before ".qc.report.tsv") and writes one row per file
with the following columns:

    sample-id    total-reads    mapped-reads    mapped-percentage    coverage-percent
"""

import argparse
import os
import csv

def parse_report_file(filepath):
    """
    Parse one qc.report.tsv file.

    Parameters
    ----------
    filepath : str
        Path to the input file.

    Returns
    -------
    dict
        A dictionary with keys:
        - sample-id
        - total-reads
        - mapped-reads
        - mapped-percentage
        - coverage-percent
    """
    # Extract sample ID from filename (prefix before ".qc.report.tsv")
    sample_id = os.path.basename(filepath).replace(".qc.report.tsv", "")
    data = {"sample-id": sample_id}

    # Read the QC report line by line
    with open(filepath, "r") as f:
        for line in f:
            if line.strip():  # skip empty lines
                key, value = line.strip().split("\t")
                if key == "total-reads:":
                    data["total-reads"] = value
                elif key == "mapped-reads:":
                    data["mapped-reads"] = value
                elif key == "mapped-percent:":
                    data["mapped-percentage"] = value
                elif key == "coverage-percent:":
                    data["coverage-percent"] = value

    return data

def main():
    """
    Parse command line arguments, process input files, and
    write the combined summary TSV file.
    """
    parser = argparse.ArgumentParser(
        description="Summarize *.qc.report.tsv files into a single TSV."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output summary TSV file"
    )
    parser.add_argument(
        "-t", "--tsv_file",
        nargs="+",
        required=True,
        help="Input .qc.report.tsv files (one or more)"
    )

    args = parser.parse_args()

    # Open the output file for writing
    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=["sample-id", "total-reads", "mapped-reads", "mapped-percentage", "coverage-percent"],
            delimiter="\t"
        )
        writer.writeheader()

        # Process each input file in order and write a row per file
        for infile in args.tsv_file:
            row = parse_report_file(infile)
            writer.writerow(row)

if __name__ == "__main__":
    main()
