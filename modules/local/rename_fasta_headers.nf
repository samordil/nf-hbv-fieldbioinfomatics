process RENAME_FASTA_HEADERS {
    tag "renaning fasta headers"
    label 'process_single'
    label 'error_ignore'

    input:
        path fasta_file
        path metadata_csv

    output:
        path "fasta_dir/*.fasta"    , emit: fasta

    script:     // This script is bundled with the pipeline, in bin
    """
    rename_fasta_headers.py \\
        --fastas ${fasta_file.join(' ')} \\
        --csv $metadata_csv \\
        --output_dir fasta_dir
    """
}