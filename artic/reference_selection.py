import pathlib
from Bio import SeqIO, SeqRecord, Seq
import argparse
import pandas as pd
import sys

from .circular import BedLine


def dealign_ref(ref: pathlib.Path, outpath: pathlib.Path):
    """
    ref: pathlib.Path
        Path to the reference fasta file
    outpath: pathlib.Path
        Path to the output fasta file
    """
    with open(outpath, "w") as out:
        for record in SeqIO.parse(ref, "fasta"):
            record.seq = record.seq.replace("-", "")
            SeqIO.write(record, out, "fasta")


def parse_samtools_coverage(cov_file: pathlib.Path) -> str:
    """
    cov_file: pathlib.Path
        Path to the samtools coverage file
    """
    df = pd.read_csv(cov_file, sep="\t")
    return df.iloc[df["numreads"].idxmax()]["#rname"]  # type: ignore


# These functions are used to map a primerscheme to a new reference genome
def read_bedlines(bedfile_path: pathlib.Path) -> list[BedLine]:
    bedlines = []
    with open(bedfile_path, "r") as f:
        for line in f:
            if not line.startswith("#"):
                data = line.strip().split("\t")
                bedlines.append(
                    BedLine(
                        data[0],
                        int(data[1]),
                        int(data[2]),
                        data[3],
                        int(data[4]),
                        data[5],
                        data[6],
                    )
                )
    return bedlines


def remap(
    bedfile_path: pathlib.Path,
    msa_path: pathlib.Path,
    id_to_remap_to: str,
    output_primer: pathlib.Path,
    output_reference: pathlib.Path,
):
    # Read in the bedlines
    bedlines = read_bedlines(bedfile_path)

    # Read in the MSA
    msa_dict = SeqIO.index(str(msa_path), "fasta")

    # Check primer's reference genome is in the msa
    chroms = set([primer.chrom for primer in bedlines])
    if len(chroms) > 1:
        # TODO spesify the primer chrom name to remap
        print("More than one reference genome in the primer.bed")
        sys.exit(1)
    remap_chrom = next(iter(chroms))
    if remap_chrom not in msa_dict:
        print(
            f"{remap_chrom} ID from primer.bed is not found in the MSA. {msa_dict}",
        )
        sys.exit(1)
    else:
        print(
            f"{remap_chrom} ID from primer.bed is found in the MSA",
        )

    # Check the new reference genome is in the msa
    if id_to_remap_to not in msa_dict:
        print(
            f"{id_to_remap_to} remapping ID is not found in the MSA",
        )
        sys.exit(1)
    else:
        print(
            f"{id_to_remap_to} remapping ID is found in the MSA",
        )

    # The the primer's reference genome to the MSA index
    primer_to_msa: dict[int, int] = {}
    ref_index = 0
    for msa_index, ref_base in enumerate(msa_dict[remap_chrom]):  # type: ignore
        if ref_base not in {"", "-"}:
            primer_to_msa[ref_index] = msa_index
            ref_index += 1

    # Create a dict that can map MSA indexes to the new reference genome
    msa_to_new_ref: dict[int, int] = {}
    new_index = 0
    for msa_index, ref_base in enumerate(msa_dict[id_to_remap_to]):  # type: ignore
        if ref_base not in {"", "-"}:
            msa_to_new_ref[msa_index] = new_index
            new_index += 1

    # Grab genome length
    msa_length = len(msa_dict.get(id_to_remap_to))  # type: ignore

    for primer in bedlines:
        # Handle the forward primer
        if primer.strand == "+":
            fp_msa = primer_to_msa[primer.end]
            fp_newref = msa_to_new_ref.get(fp_msa)

            if fp_newref is None:
                print(
                    f"Gap preventing direct mapping of {primer.primerid}_LEFT {remap_chrom}:{primer.end} -> {id_to_remap_to}"
                )
                # Walk to next valid index
                while fp_newref is None and fp_msa < msa_length:
                    fp_msa += 1
                    fp_newref = msa_to_new_ref.get(fp_msa)

                # Check fixed
                if fp_newref is None:
                    print(
                        f"Could not find or repair a valid index for {primer.primerid}_LEFT: {primer.end}"
                    )
                    sys.exit(1)
                else:
                    print(
                        f"Fixed with non-direct mapping {remap_chrom}:{primer.end} -> {id_to_remap_to}:{fp_newref}"
                    )

            primer.end = fp_newref
            primer.start = primer.end - len(primer.seq)
            primer.chrom = id_to_remap_to

        # Handle the reverse primer
        else:
            rp_msa = primer_to_msa[primer.start]
            rp_newref = msa_to_new_ref.get(rp_msa)

            if rp_newref is None:
                print(
                    f"Gap preventing direct mapping of {primer.primerid}_RIGHT {remap_chrom}:{primer.start} -> {id_to_remap_to}"
                )
                # Walk left to next valid index
                while rp_newref is None and rp_msa > 0:
                    rp_msa -= 1
                    rp_newref = msa_to_new_ref.get(rp_msa)

                if rp_newref is None:
                    print(
                        f"Could not find or repair a valid index for {primer.primerid}_RIGHT: {primer.start}"
                    )
                    sys.exit(1)
                else:
                    print(
                        f"Fixed with non-direct mapping {remap_chrom}:{primer.start} -> {id_to_remap_to}:{rp_newref}"
                    )

            primer.start = rp_newref
            primer.end = primer.start + len(primer.seq)
            primer.chrom = id_to_remap_to

    # Write the new primer.bed file
    with open(output_primer, "w") as f:
        f.write("\n".join([p.to_bed() for p in bedlines]) + "\n")

    # Write the new reference genome out
    with open(output_reference, "w") as f:
        records = [
            SeqRecord.SeqRecord(
                Seq.Seq(
                    str(msa_dict[id_to_remap_to].seq).strip().replace("-", "").upper(),  # type: ignore
                ),
                id=id_to_remap_to,
                description="",
            )
        ]
        SeqIO.write(
            records,
            f,
            "fasta",
        )


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="sub-command help", required=True)

    # Create the parser for the "dealign" command
    parser_dealign = subparsers.add_parser("dealign", help="Remove gaps from an MSA")

    parser_dealign.add_argument(
        "--ref",
        type=pathlib.Path,
        required=True,
        help="Path to the reference fasta file",
    )
    parser_dealign.add_argument(
        "--outpath",
        type=pathlib.Path,
        required=True,
        help="Path to the output fasta file",
    )
    parser_dealign.set_defaults(func=dealign_ref)

    # Create the parser for the "remap" command
    parser_remap = subparsers.add_parser(
        "remap", help="Remap a primerscheme to a new reference genome"
    )
    parser_remap.add_argument(
        "--bedfile_path",
        type=pathlib.Path,
        required=True,
        help="Path to the primer.bed file",
    )
    parser_remap.add_argument(
        "--msa_path",
        type=pathlib.Path,
        required=True,
        help="Path to the MSA file",
    )
    parser_remap.add_argument(
        "--samtools_coverage",
        type=str,
        help="The output file from samtools coverage. Used to calcuate the ID of the reference genome to remap to",
    )
    parser_remap.add_argument(
        "--id_to_remap_to",
        type=str,
        required=False,
        help="The ID of the reference genome to remap to",
    )
    parser_remap.add_argument(
        "--output_primer",
        type=pathlib.Path,
        required=True,
        help="Where to save the remapped primer.bed file",
    )
    parser_remap.add_argument(
        "--output_reference",
        type=pathlib.Path,
        required=True,
        help="Where to save the remapped reference.fasta file",
    )
    parser_remap.set_defaults(func=remap)

    args = parser.parse_args()

    if args.func == dealign_ref:
        dealign_ref(args.ref, args.outpath)
    elif args.func == remap:
        # if the id_to_remap_to is provided use it
        if args.id_to_remap_to:
            id_to_remap_to = args.id_to_remap_to
        elif args.samtools_coverage:
            id_to_remap_to = parse_samtools_coverage(
                pathlib.Path(args.samtools_coverage)
            )
        else:
            print("Please provide either --id_to_remap_to or --samtools_coverage")
            sys.exit(1)
        remap(
            args.bedfile_path,
            args.msa_path,
            id_to_remap_to,
            args.output_primer,
            args.output_reference,
        )


if __name__ == "__main__":
    main()
