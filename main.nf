#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { GENERATE_SAMPLESHEET   } from "./modules/local/generate_samplesheet"
include { ARTIC_MINION           } from "./modules/local/artic_minion.nf"
include { ASSEMBLY_STATS         } from "./modules/local/assembly_summary"
include { RENAME_FASTA_HEADERS   } from "./modules/local/rename_fasta_headers"

// Define input channel for a csv samplesheet file
channel.fromPath(params.samplesheet_csv)
        .set {ch_samplesheet_csv}

// Define an optional non-consumable input channel for a multirefence file
ch_multi_ref_file = params.multi_ref_file ?
    Channel.value(file(params.multi_ref_file)) :       // one File value
    Channel.value(null)   

// Define primerscheme directory and make sure is not consumed
Channel.value(file(params.scheme_directory))
        .set { scheme_directory_ch }

// Define fastq_pass directory channel
    Channel                                                     // Get raw fastq directory
        .fromPath(params.fastq_dir, type: 'dir', maxDepth: 1)
        .set { ch_fastq_data_dir }

// Pipeline entry point
workflow {

        // MODULE: Run bin/samplesheet_generator.py to generate samplesheet
        GENERATE_SAMPLESHEET (
        ch_fastq_data_dir,               // raw reads directory channel
        ch_samplesheet_csv,            // tsv metadata channel (can be an empty channel)
        params.matching_column
        )

        // Get the sample id with the path to the fastq file
        GENERATE_SAMPLESHEET.out.samplesheet
        .splitCsv(header:true)
        .map { row-> tuple(row.sample_name, file(row.fastq_path)) }
        .set { ch_fastq_samplesheet }

    // Run ARTIC_MINION
    ARTIC_MINION (
        scheme_directory_ch,
        params.scheme_name,
        params.scheme_version,
        params.normalize,
        params.medaka_model,
	ch_multi_ref_file,
        ch_fastq_samplesheet      // [sample_id, fastq_path]
        )

    // Add python script to rename fasta headers
    RENAME_FASTA_HEADERS (
        ARTIC_MINION.out.tsv.map{ it[1]}.collect(),
        ch_samplesheet_csv
    )

    // Generate assembly stats
    ASSEMBLY_STATS (
        ARTIC_MINION.out.tsv.map{ it[1]}.collect()
        )
}
