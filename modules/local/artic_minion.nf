process ARTIC_MINION {
    tag "assembling $sample_id"
    label 'process_medium'
    errorStrategy 'ignore'

    input:
    path scheme_dir
    val scheme_name_str
    val scheme_version_str
    val normalise_int
    val medaka_model_str
    path msa_ref_fasta
    tuple val(sample_id), path(rawReadsFastq)

    output:
    tuple val(sample_id), path("${sample_id}_output/${sample_id}.sorted.bam")                      , emit: bam
    tuple val(sample_id), path("${sample_id}_output/${sample_id}.sorted.bam.bai")                  , emit: bai
    tuple val(sample_id), path("${sample_id}_output/${sample_id}.trimmed.rg.sorted.bam")           , emit: bam_trimmed
    tuple val(sample_id), path("${sample_id}_output/${sample_id}.trimmed.rg.sorted.bam.bai")       , emit: bai_trimmed
    tuple val(sample_id), path("${sample_id}_output/${sample_id}.primertrimmed.rg.sorted.bam")     , emit: bam_primertrimmed
    tuple val(sample_id), path("${sample_id}_output/${sample_id}.primertrimmed.rg.sorted.bam.bai") , emit: bai_primertrimmed
    tuple val(sample_id), path("${sample_id}_output/${sample_id}.consensus.fasta")                 , emit: fasta
    tuple val(sample_id), path("${sample_id}_output/${sample_id}.pass.vcf")                        , emit: vcf
    tuple val(sample_id), path("${sample_id}_output/${sample_id}.pass.mod.vcf.gz.tbi")             , emit: tbi
    tuple val(sample_id), path("${sample_id}_output/${sample_id}.qc.report.tsv")                   , emit: tsv
    tuple val(sample_id), path("${sample_id}_output/${sample_id}.coverage_mask.txt.depths.png")    , emit: png
    tuple val(sample_id), path("${sample_id}_output/*.json", optional: true)  

    script:
    def multi_ref_file = msa_ref_fasta ? "--select-ref-file $msa_ref_fasta" : ''

    """
    artic minion \\
        --circular \\
        $multi_ref_file \\
        --medaka \\
        --normalise $normalise_int \\
        --threads $task.cpus \\
        --scheme-directory $scheme_dir \\
        --read-file $rawReadsFastq \\
        --medaka-model $medaka_model_str \\
        --scheme-name $scheme_name_str \\
        --scheme-version $scheme_version_str \\
        --output ${sample_id}_output \\
        $sample_id
    """
}
