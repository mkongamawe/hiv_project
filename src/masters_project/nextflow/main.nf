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
    Metadata CSV  : ${params.metadata}
    Input Reads   : ${params.reads}
    Transcriptome : ${params.transcriptome}
    Output Dir    : ${params.outdir}/${params.cohort}
    ======================================
    """.stripIndent()

    // 1. Read the metadata CSV and extract: SRR, Biosample
    meta_ch = channel
        .fromPath(params.metadata, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.Run, row.BioSample) }

    // 2. Create input channel grouping files by SRR ID
    reads_ch = channel
        .fromPath(params.reads, checkIfExists: true, type: 'file')
        .map { file -> 
            // Strip the extensions _1, _2 to group files by SRR accession}
            def srr = file.name.replaceAll(/(_[123])?\.f(ast)?q\.gz$/, '')
            tuple(srr, file)
        }
        .groupTuple(sort: true)

    // 3. Join metadata with files on the SRR key, and then regroup by BioSample
    biosample_ch = meta_ch.join(reads_ch)
        
        // Format after join: [SRR, Biosample, [read1, read2 , ...]]
        .map { srr, biosample, reads -> tuple(biosample, reads) }

        // Group by Biosample -> [BioSample, [[read1, read 2]]]
        .groupTuple()

        // Flatten the nested file list
        .map { biosample, nested_reads -> tuple(biosample, nested_reads.flatten()) }

    // 4. Load the 10x reference directory as a value channel
    transcriptome_ch = channel.value(file(params.transcriptome, type: 'dir', checkIfExists: true))

    // 5. Stage the samples
    STAGE_FASTQ(biosample_ch)

    // 4. Run Cell Ranger count
    CELLRANGER_COUNT(STAGE_FASTQ.out.staged_reads, transcriptome_ch)
}