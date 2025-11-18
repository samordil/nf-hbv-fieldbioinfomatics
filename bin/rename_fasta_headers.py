#!/usr/bin/env python
"""
Rename FASTA headers using metadata from a CSV.

Headers are updated to:
    >sample_name/ARTIC/medaka/run_id/date

- Reads `run_id` and `date` from the CSV (sample_name must match filename prefix).
- Normalizes date to ISO (YYYY-MM-DD).
- Handles many files efficiently (streamed I/O).
- Safe for HPC use (no shell calls, idempotent).
"""

import argparse
import csv
from datetime import datetime
from pathlib import Path
from Bio import SeqIO


def load_metadata(csv_file: Path):
    """Load metadata keyed by sample_name, normalizing dates."""
    metadata = {}
    with csv_file.open(newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_name = row["sample_name"].strip()
            run_id = row["run_id"].strip()
            raw_date = row["date"].strip()

            # Try common date formats → convert to YYYY-MM-DD
            for fmt in ("%d/%m/%Y", "%d-%m-%Y", "%Y/%m/%d", "%Y-%m-%d"):
                try:
                    raw_date = datetime.strptime(raw_date, fmt).strftime("%Y-%m-%d")
                    break
                except ValueError:
                    pass

            metadata[sample_name] = {"run_id": run_id, "date": raw_date}
    return metadata


def rename_fasta_headers(fasta_files, metadata, output_dir: Path, overwrite=False):
    """Rewrite FASTA headers using metadata."""
    if not overwrite:
        output_dir.mkdir(parents=True, exist_ok=True)

    for fasta_path in fasta_files:
        fasta_path = Path(fasta_path)
        if not fasta_path.exists():
            print(f"File not found: {fasta_path}")
            continue

        sample_name = fasta_path.stem.split(".")[0]
        if sample_name not in metadata:
            print(f"No metadata for {sample_name}, skipping.")
            continue

        run_id, date = metadata[sample_name]["run_id"], metadata[sample_name]["date"]
        output_path = fasta_path if overwrite else output_dir / fasta_path.name

        records = []
        for record in SeqIO.parse(fasta_path, "fasta"):
            record.id = f"{sample_name}/ARTIC/medaka/{run_id}/{date}"
            record.description = ""
            records.append(record)

        SeqIO.write(records, output_path, "fasta")
        print(f"{fasta_path.name} → {record.id}")

    print("Renaming complete.")

def main():
    parser = argparse.ArgumentParser(description="Rename FASTA headers using CSV metadata.")
    parser.add_argument("-c", "--csv", required=True, type=Path, help="CSV with sample_name, run_id, date.")
    parser.add_argument("-f", "--fastas", required=True, nargs="+", type=Path, help="FASTA files (supports wildcards).")
    parser.add_argument("-o", "--output_dir", default=Path("renamed_fastas"), type=Path, help="Output directory.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite original files.")
    args = parser.parse_args()

    metadata = load_metadata(args.csv)
    rename_fasta_headers(args.fastas, metadata, args.output_dir, args.overwrite)


if __name__ == "__main__":
    main()
