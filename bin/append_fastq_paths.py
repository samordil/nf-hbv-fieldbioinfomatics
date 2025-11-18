#!/usr/bin/env python
"""
append_fastq_paths.py

Match sample IDs from a CSV metadata file to FASTQ files in a given directory.

The script looks for filenames in the FASTQ directory that contain the
sample ID (from a specified column in the CSV). It then appends a new
column (`fastq_path`) to the CSV file, containing either the absolute path
to the matched FASTQ file or "NOT_FOUND" if no match was found.

Output is written as a new CSV file with the same columns as the input
plus the additional `fastq_path` column.
"""

import csv
import argparse
from pathlib import Path


def match_fastq_paths(csv_file, fastq_dir, output_file, match_column):
    """
    Match sample IDs from a CSV file to FASTQ files in a directory.

    Parameters
    ----------
    csv_file : str
        Path to the input CSV metadata file.
    fastq_dir : str
        Directory containing FASTQ files (.fastq or .fastq.gz).
    output_file : str
        Path to the output CSV file with appended FASTQ paths.
    match_column : str
        Column name in the CSV file used to match against FASTQ filenames.

    Returns
    -------
    None
        Writes a new CSV file with an additional `fastq_path` column.
    """
    # Convert to Path objects
    csv_path = Path(csv_file)
    fastq_path = Path(fastq_dir)
    output_path = Path(output_file)

    # Verify FASTQ directory exists
    if not fastq_path.exists():
        raise ValueError(f"FASTQ directory does not exist: {fastq_path}")
    if not fastq_path.is_dir():
        raise ValueError(f"FASTQ path is not a directory: {fastq_path}")

    # Open and read the input CSV into memory
    with csv_path.open(newline='') as infile:
        reader = csv.DictReader(infile)
        rows = list(reader)

    # Collect all FASTQ files in the directory with debugging
    fastq_files = []
    for f in fastq_path.iterdir():
        if f.is_file() and (f.suffix in ('.fastq', '.gz') or f.name.endswith(('.fastq.gz', '.fq.gz'))):
            fastq_files.append(f)
    
    print(f"Found {len(fastq_files)} FASTQ files in {fastq_path}")
    print("FASTQ files found:")
    for f in fastq_files:
        print(f"  - {f.name}")

    # Match rows to FASTQ files with detailed debugging
    for row in rows:
        key = row[match_column]
        matched_file = None
        
        print(f"\nLooking for sample ID: '{key}'")
        
        # Find the first FASTQ file containing the key in its name
        for f in fastq_files:
            if key in f.name:
                matched_file = f.resolve()  # Get absolute path
                print(f"  MATCH: '{key}' found in '{f.name}'")
                break
            else:
                print(f"  No match: '{key}' not in '{f.name}'")
        
        # Add the result (absolute path or "NOT_FOUND") as a new column
        row["fastq_path"] = str(matched_file) if matched_file else "NOT_FOUND"
        print(f"  Result: {row['fastq_path']}")

    # Write the updated data to a new CSV file
    fieldnames = list(rows[0].keys())
    with output_path.open("w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    
    print(f"\nOutput written to: {output_path}")


def main():
    """
    Parse command-line arguments and run the FASTQ matching workflow.
    """
    parser = argparse.ArgumentParser(
        description="Match CSV strain/sample IDs to FASTQ files and append file paths."
    )
    parser.add_argument("-c", "--csv_file", required=True, help="Input CSV metadata file")
    parser.add_argument("-d", "--fastq_dir", required=True, help="Directory containing FASTQ files")
    parser.add_argument("-o", "--output_file", required=True, help="Output CSV file")
    parser.add_argument("-m", "--match_column", required=True, help="Column name in CSV used for matching against FASTQ filenames")

    args = parser.parse_args()

    match_fastq_paths(args.csv_file, args.fastq_dir, args.output_file, args.match_column)


if __name__ == "__main__":
    main()

