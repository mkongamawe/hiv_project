#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import required processes
include { STAGE_FASTQ; CELLRANGER_COUNT } from './modules/cellranger.nf'

workflow {

    // Create a QC header
    log.info """
    ======================================
    S C R N A - S E Q  P I P E L I N E
    ======================================
    Cohort        : ${params.cohort}
    Input Reads   : ${params.reads}
    Transcriptome : ${params.transcriptome}
    Output Dir    : ${params.outdir}/${params.cohort}
    ======================================
    """.stripIndent()

    // 1. Create input channel grouping files by SRR ID
    reads_ch = channel
        .fromPath(params.reads, checkIfExists: true)
        .map { file -> 
            // Strip the extensions _1, _2 to group files by SRR accession}
            def sample_id = file.name.replaceAll(/(_[123])?\.f(ast)?q\.gz$/, '')
            tuple(sample_id, file)
        }
        .groupTuple(sort: true)

    // 2. Load the 10x reference directory as a value channel
    transcriptome_ch = channel.value(file(params.transcriptome, type: 'dir', checkIfExists: true))

    // 3. Rename files to match filenames
    STAGE_FASTQ(reads_ch)

    // 4. Run Cell Ranger count
    CELLRANGER_COUNT(STAGE_FASTQ.out.staged_reads, transcriptome_ch)
}