# Written by Nick Loman (@pathogenomenick)

from clint.textui import colored
import os
import sys
import time
import requests
import hashlib
import pathlib
from .vcftagprimersites import read_bed_file
from .circular import create_or_find_circular_scheme


def check_scheme_hashes(filepath, manifest_hash):
    with open(filepath, "rb") as fh:
        data = fh.read()
        hash_sha256 = hashlib.sha256(data).hexdigest()
    if hash_sha256 != manifest_hash:
        print(
            colored.yellow(f"sha256 hash for {str(filepath)} does not match manifest"),
            file=sys.stderr,
        )
        raise SystemExit(1)


def run_commands(cmds: list[str], dry_run: bool, logfile: pathlib.Path):
    with open(logfile, "a") as logfh:
        for cmd in cmds:
            print(colored.green("Running: ") + cmd, file=sys.stderr)
            if not dry_run:
                timerStart = time.perf_counter()
                retval = os.system(cmd)
                if retval != 0:
                    print(colored.red("Command failed:") + cmd, file=sys.stderr)
                    raise SystemExit(20)
                timerStop = time.perf_counter()

                ## print the executed command and the runtime to the log file
                print("{}\t{}".format(cmd, timerStop - timerStart), file=logfh)

            ## if it's a dry run, print just the command
            else:
                print(cmd, file=logfh)


def get_scheme(
    scheme_directory: pathlib.Path, scheme_name: str, scheme_version: str
) -> tuple[pathlib.Path, pathlib.Path, str]:
    """Get and check the ARTIC primer scheme.
    When determining a version, the original behaviour (parsing the scheme_name and
    separating on /V ) is used over a specified scheme_version. If neither are
    provided, the version defaults to 1.
    If 0 is provided as version, the latest scheme will be downloaded.

    Parameters
    ----------
    scheme_name : str
        The primer scheme name
    scheme_directory : pathlib.Path
        The directory containing the primer scheme and reference sequence
    scheme_version : str
        The primer scheme version (optional)
    Returns
    -------
    pathlib.Path
        The location of the checked primer scheme
    pathlib.Path
        The location of the checked reference sequence
    str
        The version being used
    """
    # create the filenames and check they exist
    bed: pathlib.Path = (
        scheme_directory / scheme_name / scheme_version / f"{scheme_name}.scheme.bed"
    )
    ref: pathlib.Path = (
        scheme_directory
        / scheme_name
        / scheme_version
        / f"{scheme_name}.reference.fasta"
    )

    if bed.is_file() and ref.is_file():
        return bed, ref, scheme_version

    # Dont download

    print(list(ref.parent.glob("*")))
    raise SystemExit(2)

    # if they don't exist, try downloading them to the current directory
    print(
        colored.yellow(
            "could not find primer scheme and reference sequence, downloading"
        ),
        file=sys.stderr,
    )

    try:
        manifest = requests.get(
            "https://raw.githubusercontent.com/artic-network/primer-schemes/master/schemes_manifest.json"
        ).json()
    except requests.exceptions.RequestException as error:
        print("Manifest Exception:", error)
        raise SystemExit(2)

    for scheme, scheme_contents in dict(manifest["schemes"]).items():
        if (
            scheme == scheme_name.lower()
            or scheme_name.lower() in scheme_contents["aliases"]
        ):
            print(
                colored.yellow(
                    f"\tfound requested scheme:\t{scheme} (using alias {scheme_name})"
                ),
                file=sys.stderr,
            )
            if scheme_version == 0:
                print(
                    colored.yellow(
                        f"Latest version for scheme {scheme} is -> {scheme_contents['latest_version']}"
                    ),
                    file=sys.stderr,
                )
                scheme_version = scheme_contents["latest_version"]
            elif scheme_version not in dict(scheme_contents["primer_urls"]).keys():
                print(
                    colored.yellow(
                        f"Requested scheme version {scheme_version} not found; using latest version ({scheme_contents['latest_version']}) instead"
                    ),
                    file=sys.stderr,
                )
                scheme_version = scheme_contents["latest_version"]
                bed: pathlib.Path = (
                    scheme_directory
                    / scheme_name
                    / scheme_version
                    / f"{scheme_name}.scheme.bed"
                )
                ref: pathlib.Path = (
                    scheme_directory
                    / scheme_name
                    / scheme_version
                    / f"{scheme_name}.reference.fasta"
                )

            bed.parent.mkdir(parents=True, exist_ok=True)
            with requests.get(scheme_contents["primer_urls"][scheme_version]) as fh:
                open(bed, "wt").write(fh.text)

            os.makedirs(os.path.dirname(ref), exist_ok=True)
            with requests.get(scheme_contents["reference_urls"][scheme_version]) as fh:
                open(ref, "wt").write(fh.text)

            check_scheme_hashes(
                bed, scheme_contents["primer_sha256_checksums"][scheme_version]
            )
            check_scheme_hashes(
                ref, scheme_contents["reference_sha256_checksums"][scheme_version]
            )

            return bed, ref, scheme_version

    print(
        colored.yellow(
            f"\tRequested scheme:\t{scheme_name} could not be found, exiting"
        ),
        file=sys.stderr,
    )
    raise SystemExit(1)


def run(parser, args):
    # check for medaka-model
    if args.medaka and (args.medaka_model is None):
        print(
            colored.red("Must specify --medaka-model if using the --medaka workflow.")
        )
        raise SystemExit(1)

    ## create a holder to keep the pipeline commands in
    cmds = []

    # 1) check the parameters and set up the filenames
    ## find the primer scheme, reference sequence and confirm scheme version
    bed, ref, _ = get_scheme(
        args.scheme_directory, args.scheme_name, args.scheme_version
    )

    # 2) check the input and output files
    bed = bed.absolute()
    ref = ref.absolute()

    ## set up the read file
    read_files_str = " ".join([str(x.absolute()) for x in args.read_files])
    output_dir: pathlib.Path = args.output.absolute()
    log = output_dir / f"{args.sample}.minion.log.txt"

    try:
        output_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print(
            colored.red("output already exists: {}".format(output_dir)),
            file=sys.stderr,
        )
        raise SystemExit(1)

    # Check for the best match is select-ref-file is provided
    if args.select_ref_file is not None:
        select_ref_cmds = []
        if not args.select_ref_file.is_file():
            print(
                colored.red(
                    f"The provided select-ref-file does not exist: {args.select_ref_file}"
                ),
                file=sys.stderr,
            )
            raise SystemExit(1)

        # Create the degapped reference file
        degapped_fasta = output_dir / f"{args.sample}.degapped.fasta"
        select_ref_cmds.append(
            f"artic_reference_selection dealign --ref {args.select_ref_file} --outpath {degapped_fasta}"
        )
        # Map all reads to the degapped reference
        all_ref_sorted_bam = output_dir / f"{args.sample}.all_ref.sorted.bam"
        select_ref_cmds.append(
            f"minimap2 -a -x map-ont -t {args.threads} {degapped_fasta} {read_files_str} | samtools sort -o {all_ref_sorted_bam}"
        )
        # Run samtools coverage on the degapped reference
        reference_selection_coverage = (
            output_dir / f"{args.sample}.reference_selection.coverage"
        )
        select_ref_cmds.append(
            f"samtools coverage {all_ref_sorted_bam} > {reference_selection_coverage}"
        )
        # Remap the primer scheme to the best reference
        remapped_bed = output_dir / f"{args.sample}.remapped.scheme.bed"
        remapped_ref = output_dir / f"{args.sample}.remapped.reference.fasta"
        select_ref_cmds.append(
            f"artic_reference_selection remap --bedfile_path {bed} --msa_path {args.select_ref_file} --output_primer {remapped_bed} --output_reference {remapped_ref} --samtools_coverage {reference_selection_coverage}"
        )
        bed = remapped_bed
        ref = remapped_ref
        # Remove the degapped reference and the all_ref_sorted_bam
        select_ref_cmds.append(f"rm {all_ref_sorted_bam} {degapped_fasta}")

        # Run the select-ref-file commands
        run_commands(select_ref_cmds, args.dry_run, log)

    if args.circular:
        lbed = bed
        lref = ref

        bed = output_dir / f"{args.sample}.circular.scheme.bed"
        ref = output_dir / f"{args.sample}.circular.reference.fasta"

        # if the circular flag is set, create the circular scheme and overwrite the bed and ref
        _, _, reflen = create_or_find_circular_scheme(
            bedfile_path=lbed, ref_path=lref, output_primer=bed, output_reference=ref
        )

    ## if in strict mode, validate the primer scheme
    if args.strict:
        checkScheme = "artic-tools validate_scheme %s" % (bed)
        print(colored.green("Running: "), checkScheme, file=sys.stderr)
        if os.system(checkScheme) != 0:
            print(colored.red("primer scheme failed strict checking"), file=sys.stderr)
            raise SystemExit(1)

    ## collect the primer pools
    pools = set([row["PoolName"] for row in read_bed_file(bed)])

    # 3) index the ref & align with minimap
    all_sorted_bam = output_dir / f"{args.sample}.all.sorted.bam"  # All mapped reads
    mappingstats = output_dir / f"{args.sample}.mappingstats.json"  # Mapping stats
    sorted_bam = output_dir / f"{args.sample}.sorted.bam"  # Mapped reads

    # aligntrim
    alignreport_txt = output_dir / f"{args.sample}.alignreport.txt"
    alignreport_er = output_dir / f"{args.sample}.alignreport.er"
    trimmed_rg_sorted_bam = output_dir / f"{args.sample}.trimmed.rg.sorted.bam"
    primertrimmed_rg_sorted_bam = (
        output_dir / f"{args.sample}.primertrimmed.rg.sorted.bam"
    )

    cmds.append(
        f"minimap2 -a -x map-ont -t {args.threads} {ref} {read_files_str} | samtools sort -o {all_sorted_bam}"
    )
    # get mapping stats
    cmds.append(f"samtools flagstat {all_sorted_bam} -O json > {mappingstats}")
    # Filter out unmapped reads then delete all.bam
    cmds.append(
        f"samtools view -bS -F 4 {all_sorted_bam} > {sorted_bam} && rm {all_sorted_bam}"
    )
    cmds.append(f"samtools index {sorted_bam}")

    # 4) trim the alignments to the primer start sites and normalise the coverage to save time
    if args.normalise:
        normalise_string = f"--normalise {args.normalise}"
    else:
        normalise_string = ""
    cmds.append(
        f"align_trim {normalise_string} {bed} --start --remove-incorrect-pairs --report {alignreport_txt} < {sorted_bam} 2> {alignreport_er} | samtools sort -T {args.sample} - -o {trimmed_rg_sorted_bam}"
    )
    cmds.append(
        f"align_trim {normalise_string} {bed} --remove-incorrect-pairs --report {alignreport_txt} < {sorted_bam} 2> {alignreport_er} | samtools sort -T %s - -o {primertrimmed_rg_sorted_bam}"
    )
    cmds.append(f"samtools index {trimmed_rg_sorted_bam}")
    cmds.append(f"samtools index {primertrimmed_rg_sorted_bam}")

    # 6) medaka. do variant calling on each read group, either using medaka
    merge_vcf_files = []
    for p in pools:
        _rg_hdf = output_dir / f"{args.sample}.{p}.hdf"
        _rg_vcf = output_dir / f"{args.sample}.{p}.vcf"

        if _rg_hdf.exists():
            _rg_hdf.unlink()
        cmds.append(
            f"medaka consensus --model {args.medaka_model} --threads {args.threads} --chunk_len 800 --chunk_ovlp 400 --RG {p} {primertrimmed_rg_sorted_bam} {_rg_hdf}"
        )
        medaka_cmd = "snp" if args.no_indels else "variant"
        cmds.append(f"medaka {medaka_cmd} {ref} {_rg_hdf} {_rg_vcf}")

        ## if not using longshot, annotate VCF with read depth info etc. so we can filter it
        if args.no_longshot:
            cmds.append(
                f"medaka tools annotate --pad 25 --RG {p} {_rg_vcf} {ref} {primertrimmed_rg_sorted_bam} {output_dir}/tmp.medaka-annotate.vcf"
            )
            cmds.append(f"mv {output_dir}/tmp.medaka-annotate.vcf {_rg_vcf}")

        merge_vcf_files.append(f"{p}:{_rg_vcf}")

    # 7) merge the called variants for each read group
    primersitereport_txt = output_dir / f"{args.sample}.primersitereport.txt"
    merged_vcf = output_dir / f"{args.sample}.merged.vcf"
    primer_vcf = output_dir / f"{args.sample}.primer.vcf"

    merged_filtered_vcf = output_dir / f"{args.sample}.merged.filtered.vcf"
    vcfreport_txt = output_dir / f"{args.sample}.vcfreport.txt"

    merge_vcf_cmd = f"artic_vcf_merge --output_merged {merged_vcf} --output_primer {primer_vcf} {args.sample} {bed} 2> {primersitereport_txt} {' '.join([str(x) for x in merge_vcf_files])}"
    cmds.append(merge_vcf_cmd)

    # 8) check and filter the VCFs
    ## if using strict, run the vcf checker to remove vars present only once in overlap regions (this replaces the original merged vcf from the previous step)
    if args.strict:
        cmds.append(f"bgzip -f {merged_vcf}")
        cmds.append(f"tabix -p vcf {merged_vcf}.gz")
        cmds.append(
            f"artic-tools check_vcf --dropPrimerVars --dropOverlapFails --vcfOut {merged_filtered_vcf} {merged_vcf}.gz {bed} 2> {vcfreport_txt}"
        )
        cmds.append(f"mv {merged_filtered_vcf} {merged_vcf}")

    ## if doing the medaka workflow and longshot required, do it on the merged VCF
    if args.medaka and not args.no_longshot:
        cmds.append(f"bgzip -f {merged_vcf}")
        cmds.append(f"tabix -f -p vcf {merged_vcf}.gz")
        cmds.append(
            f"longshot -P 0 -F -A --no_haps --bam {primertrimmed_rg_sorted_bam} --ref {ref} --out {merged_vcf} --potential_variants {merged_vcf}.gz"
        )

    ## set up some name holder vars for ease

    method = "medaka"

    pass_vcf = output_dir / f"{args.sample}.pass.vcf"
    fail_vcf = output_dir / f"{args.sample}.fail.vcf"

    ## filter the variants to produce PASS and FAIL lists, then index them
    if args.no_frameshifts and not args.no_indels:
        cmds.append(
            f"artic_vcf_filter --{method} --no-frameshifts {merged_vcf} {pass_vcf} {fail_vcf}"
        )
    else:
        cmds.append(f"artic_vcf_filter  --{method} {merged_vcf} {pass_vcf} {fail_vcf}")

    # Modify the vcf.pass
    if args.circular:
        pass_mod_vcf = pass_vcf.with_suffix(".mod.vcf")
        fail_mod_vcf = fail_vcf.with_suffix(".mod.vcf")
        cmds.append(
            f"artic_circular parse-vcf --vcf {pass_vcf} --reference {lref} --output {pass_mod_vcf}"
        )  # Writes vcf
        pass_vcf = pass_mod_vcf
        cmds.append(
            f"artic_circular dedupe-vcf --pass-vcf {pass_vcf} --fail-vcf {fail_vcf} --lref {lref} --output {fail_mod_vcf}"
        )  # Writes vcf
        fail_vcf = fail_mod_vcf

    cmds.append(f"bgzip -f {pass_vcf}")
    cmds.append(f"tabix -p vcf {pass_vcf}.gz")

    # 9) get the depth of coverage for each readgroup, create a coverage mask and plots, and add failed variants to the coverage mask (artic_mask must be run before bcftools consensus)
    coverage_mask_txt = output_dir / f"{args.sample}.coverage_mask.txt"
    preconsensus_fasta = output_dir / f"{args.sample}.preconsensus.fasta"

    if args.circular:
        # modulo the depth mask to turn circular indexing back into linear
        cmds.append(
            f"artic_make_depth_mask --store-rg-depths {ref} --de-circ-reflength {reflen} {primertrimmed_rg_sorted_bam} {coverage_mask_txt}"
        )
        cmds.append(
            f"artic_mask {lref} {coverage_mask_txt} {fail_vcf} {preconsensus_fasta}"
        )

    else:
        cmds.append(
            f"artic_make_depth_mask --store-rg-depths {ref}  {primertrimmed_rg_sorted_bam} {coverage_mask_txt}"
        )
        cmds.append(
            f"artic_mask {ref} {coverage_mask_txt} {fail_vcf} {preconsensus_fasta}"
        )

    # Create QC depth plots
    cmds.append(
        f"artic_make_depth_plot --depth {coverage_mask_txt}.depths --min-depth 20 --output {coverage_mask_txt}.depths"
    )

    # 10) generate the consensus sequence
    consensus_fasta = output_dir / f"{args.sample}.consensus.fasta"
    cmds.append(
        f"bcftools consensus -f {preconsensus_fasta} {pass_vcf}.gz -m {coverage_mask_txt} -o {consensus_fasta}"
    )

    # 11) apply the header to the consensus sequence and run alignment against the reference sequence
    muscle_in_fasta = output_dir / f"{args.sample}.muscle.in.fasta"
    muscle_out_fasta = output_dir / f"{args.sample}.muscle.out.fasta"

    fasta_header = f"{args.sample}/ARTIC/{method}"
    cmds.append(f"artic_fasta_header {consensus_fasta} '{fasta_header}'")
    cmds.append(f"cat {consensus_fasta} {lref} > {muscle_in_fasta}")
    cmds.append(f"muscle -in {muscle_in_fasta} -out {muscle_out_fasta}")

    # 12) get some QC stats
    if args.strict:
        cmds.append(
            f"artic_get_stats --scheme {bed} --align-report {alignreport_txt} --vcf-report {vcfreport_txt} {args.sample}"
        )

    # Create a qc tsv
    qc_report = output_dir / f"{args.sample}.qc.report.tsv"
    cmds.append(
        f"artic_qc_report --mappingstats {mappingstats} --consensus {consensus_fasta} --output {qc_report}"
    )

    # 13) setup the log file and run the pipeline commands
    run_commands(cmds, args.dry_run, log)
