process GENERATE_SAMPLESHEET {
    publishDir "${params.outDir}/assembly", mode:'copy'
    tag "samplesheet generation"
    label 'process_single'
    label 'error_ignore'

    input:
        path raw_data_dir
        path metadata_csv
        val column_to_match

    output:
        path "samplesheet_with_data_path.csv"    , emit: samplesheet


    script:     // This script is bundled with the pipeline, in bin
    def column  = column_to_match ? "$column_to_match" : ""

    """
    append_fastq_paths.py \\
        --csv_file $metadata_csv \\
        --fastq_dir $raw_data_dir \\
        --match_column $column \\
        --output_file samplesheet_with_data_path.csv
    """
}