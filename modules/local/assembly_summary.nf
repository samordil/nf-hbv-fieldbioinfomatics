process ASSEMBLY_STATS {
    publishDir "${params.outDir}/assembly", mode:'copy'
    tag "Assembly stats"
    label 'process_single'
    label 'error_ignore'

    input:
        path tsv_file

    output:
        path "assembly.stats.tsv"    , emit: tsv


    script:     // This script is bundled with the pipeline, in bin
    """
    qc_summary.py --tsv_file $tsv_file --output assembly.stats.tsv

    """
}