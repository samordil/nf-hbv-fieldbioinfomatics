import json
from Bio import SeqIO
import argparse

def create_qc(mappingstats: str, consensus: str, output: str):
    mapping_data = json.load(open(mappingstats))
    qc_pass = mapping_data["QC-passed reads"]

    total_reads = qc_pass["total"]
    mapped_reads = qc_pass["mapped"]
    percent_mapped = qc_pass["mapped %"]

    # percentage coverage
    consensus_genome = SeqIO.read(consensus, "fasta")
    coverage_percent = round(((len(consensus_genome) - consensus_genome.seq.count("N")) / len(consensus_genome)) * 100, 2)

    
    with open(output, "w") as f:
        f.write(f"total-reads:\t{total_reads}\n")
        f.write(f"mapped-reads:\t{mapped_reads}\n")
        f.write(f"mapped-percent:\t{percent_mapped}\n")
        f.write(f"coverage-percent:\t{coverage_percent}\n")



def main():
    parser = argparse.ArgumentParser(description="Collate QC data about the run")
    parser.add_argument('--mappingstats', type=str, required=True, help="Mapping stats file")
    parser.add_argument('--consensus', type=str, required=True, help="Mapping stats file")
    parser.add_argument('--output', type=str, required=True, help="Output file location")

    args = parser.parse_args()
    create_qc(args.mappingstats, args.consensus, args.output)

if __name__ == "__main__":
    main()

