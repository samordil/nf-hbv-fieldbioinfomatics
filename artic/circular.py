# Writen for the HBV analysis of circular genomes
# Chris Kent

from Bio import SeqIO
import vcf
import argparse
import pathlib


class BedLine:
    chrom: str
    start: int
    end: int
    primerid: str
    pool: int
    strand: str
    seq: str

    def __init__(self, chrom, start, end, primerid, pool, strand, seq):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.primerid = primerid
        self.pool = int(pool)
        self.strand = strand
        self.seq = seq

    def to_bed(self):
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.primerid}\t{self.pool}\t{self.strand}\t{self.seq}"


def generate_amplicons(bedfile) -> dict[str, list[list[BedLine]]]:
    # Read in the bedlines as lists
    bedlines: list[BedLine] = []
    with open(bedfile, "r") as infile:
        for line in infile:
            splitline = line.strip().split("\t")
            if splitline:
                bedlines.append(
                    BedLine(
                        chrom=splitline[0],
                        start=splitline[1],
                        end=splitline[2],
                        primerid=splitline[3],
                        pool=splitline[4],
                        strand=splitline[5],
                        seq=splitline[6],
                    )
                )

    # group bedlines by Amplicon_number
    amplicons = {}
    for bedline in bedlines:
        ampliconID = "_".join(bedline.primerid.split("_")[:2])
        direction = bedline.primerid.split("_")[-1]

        if ampliconID not in amplicons:
            amplicons[ampliconID] = [[], []]

        if direction == "LEFT":
            amplicons[ampliconID][0].append(bedline)
        elif direction == "RIGHT":
            amplicons[ampliconID][1].append(bedline)

    return amplicons


def create_or_find_circular_scheme_cli(args):
    create_or_find_circular_scheme(
        args.bedfile,
        args.ref,
        args.output_primer,
        args.output_reference,
    )


def create_or_find_circular_scheme(
    bedfile_path: pathlib.Path,
    ref_path: pathlib.Path,
    output_primer: pathlib.Path,
    output_reference: pathlib.Path,
) -> tuple[pathlib.Path, pathlib.Path, int]:
    """
    Returns:
        cbed: pathlib.Path, path to the circular bed file
        cref: pathlib.Path, path to the circular reference file
        reflen: int, orginal length of the reference genome
    """

    # Check if orignal bed and ref and found locally

    if not bedfile_path.is_file():
        raise FileNotFoundError(f"Scheme not found at {bedfile_path}")
    if not ref_path.is_file():
        raise FileNotFoundError(f"Reference not found at {ref_path}")

    # Read in linear bedfile
    amplicons = generate_amplicons(bedfile_path)
    # Find circular amplicons
    circular_amplicons = []
    for ampliconID, (lprimers, rprimers) in amplicons.items():
        # If left primer > right primer
        if int(lprimers[0].end) > int(rprimers[0].start):
            circular_amplicons.append((lprimers, rprimers))
    if not circular_amplicons:
        raise ValueError("No circular amplicons found")

    # Write circular ref
    cref = output_reference
    refrecord = SeqIO.read(ref_path, "fasta")
    reflen = len(refrecord)
    furthest_right = max(
        [rp.end for rp in (rps for lps, rps in circular_amplicons) for rp in rp]
    )

    # Append the circular region to the end of the reference genome
    refrecord.seq = refrecord.seq + refrecord.seq[:furthest_right]
    refrecord.id = refrecord.id + "_circular"
    refrecord.description = ""

    with open(cref, "w") as f:
        SeqIO.write(refrecord, f, "fasta")

    # Write circular bed
    cbed = output_primer

    bedfile_str = []
    for ampliconID, (lprimers, rprimers) in amplicons.items():
        # If left primer > right primer
        if int(lprimers[0].end) > int(rprimers[0].start):
            for p in rprimers:
                p.start += reflen
                p.end += reflen

        for lp in lprimers:
            lp.chrom = refrecord.id
            bedfile_str.append(lp.to_bed())
        for rp in rprimers:
            rp.chrom = refrecord.id
            bedfile_str.append(rp.to_bed())

    with open(cbed, "w") as f:
        f.write("\n".join(bedfile_str) + "\n")

    return cbed, cref, reflen


def decirc_vcf(args):
    # Load the reference genome
    ref = SeqIO.read(args.reference, "fasta")
    ref_len = len(ref)

    # Load the circular consensus genome vcf
    ccvcf = vcf.Reader(open(args.vcf, "r"), filename=args.vcf)
    # Open a new vcf file to write to
    vcf_writer = vcf.Writer(
        open(args.output, "w"),
        ccvcf,
        lineterminator="\n",
    )

    records = []
    for record in ccvcf:
        if record.ALT[0] != ".":
            record.POS = (int(record.POS - 1) % int(ref_len)) + 1
            record.CHROM = ref.id.strip()

        records.append(record)

    records.sort(key=lambda x: x.POS)

    for record in records:
        vcf_writer.write_record(record)


def parse_fail_vcf(args):
    # Load the reference genome
    lref = SeqIO.read(args.lref, "fasta")
    lref_len = len(lref)

    # Load the linear pass vcf
    pass_vcf = vcf.Reader(open(args.pass_vcf, "r"), filename=args.pass_vcf)
    pass_records = {}
    for record in pass_vcf:
        for i, _base in enumerate(record.REF):
            pass_records[(record.CHROM, record.POS + i)] = record

    # Read in the fail vcf
    fail_vcf = vcf.Reader(open(args.fail_vcf, "r"), filename=args.fail_vcf)
    fail_vcf_writer = vcf.Writer(
        open(args.output, "w"),
        fail_vcf,
        lineterminator="\n",
    )
    for failrecord in fail_vcf:
        write_fail = True
        # Remap the fail record to the linear reference genome
        if failrecord.ALT[0] != ".":
            failrecord.POS = (int(failrecord.POS - 1) % int(lref_len)) + 1
            failrecord.CHROM = lref.id.strip()

        # Check if the position is in the pass vcf

        for i, _ in enumerate(str(failrecord.ALT[0])):
            if (failrecord.CHROM, failrecord.POS + i) in pass_records:
                write_fail = False
                break

        if write_fail:
            fail_vcf_writer.write_record(failrecord)


# This takes the consensus genome generated from the circular reference genome and maps the circular coordinates back to the linear reference genome
def main():
    global_parser = argparse.ArgumentParser(
        description="Tools for handling circular genomes",
    )
    subparsers = global_parser.add_subparsers(
        title="subcommands", help="scheme types", required=True
    )
    parse_vcf_parser = subparsers.add_parser(
        "parse-vcf", help="map a circular genome back to linear coords"
    )
    parse_vcf_parser.add_argument(
        "--vcf",
        type=str,
        help="The vcf generated from the circular reference genome",
        required=True,
    )
    parse_vcf_parser.add_argument(
        "--reference",
        "-r",
        type=str,
        help="The reference genome in fasta format",
        required=True,
    )
    parse_vcf_parser.add_argument(
        "-o", "--output", type=str, help="The output vcf file", required=True
    )
    parse_vcf_parser.set_defaults(func=decirc_vcf)

    # Parse and De-dupe fail/vcf
    de_dupe_vcf_parser = subparsers.add_parser(
        "dedupe-vcf", help="Removed fail positions in a valid vcf"
    )
    de_dupe_vcf_parser.add_argument("--pass-vcf", type=str, help="The linear pass.vcf")
    de_dupe_vcf_parser.add_argument("--fail-vcf", type=str, help="The fail.vcf")
    de_dupe_vcf_parser.add_argument(
        "--lref", type=str, help="The linear referance genome"
    )
    de_dupe_vcf_parser.add_argument(
        "-o", "--output", type=str, help="The output vcf file", required=True
    )
    de_dupe_vcf_parser.set_defaults(func=parse_fail_vcf)

    # Create or find circular scheme
    create_circular_parser = subparsers.add_parser(
        "create-circular-scheme", help="Create a circular scheme"
    )
    create_circular_parser.add_argument(
        "--bedfile",
        type=pathlib.Path,
        help="The bed file of the circular scheme",
        required=True,
    )
    create_circular_parser.add_argument(
        "--ref",
        type=pathlib.Path,
        help="The reference genome in fasta format",
        required=True,
    )
    create_circular_parser.add_argument(
        "--output-primer",
        type=pathlib.Path,
        help="The output bed file",
        required=True,
    )
    create_circular_parser.add_argument(
        "--output-reference",
        type=pathlib.Path,
        help="The output reference genome",
        required=True,
    )
    create_circular_parser.set_defaults(func=create_or_find_circular_scheme)

    # Run
    args = global_parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
