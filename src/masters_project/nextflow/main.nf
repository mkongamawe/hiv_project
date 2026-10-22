#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import required processes
include { runFastQc; runMultiQc } from './modules/data_qc.nf'

workflow {

    // Create a QC header
    log.info """
    ======================================
    Q C  P I P E L I N E
    ======================================
    Input Reads: ${params.reads}
    Output Dir : ${params.outdir}
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

    // 2. Run fastqc
    runFastQc(reads_ch)

    // 3. Collect all fastqc files and pass to multiqc
    runMultiQc(runFastQc.out.qc_files.collect())
}