#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { INPUT_PREPROCESSING  } from './subworkflows/input_preprocessing.nf'
include { CELLRANGER_ALIGNMENT } from './subworkflows/cellranger_alignment.nf'
include { FASTQC; MULTIQC      } from './modules/qc.nf'

workflow {

    // Create a QC header
    log.info """
    ======================================
    S C R N A - S E Q  P I P E L I N E
    ======================================
    Mode           : ${params.mode}
    Cohort         : ${params.cohort}
    BCL Input Dir  : ${params.bcl_dir ?: 'Skipped (Using FASTQs)'}
    Input Reads    : ${params.reads}
    SRA Metadata   : ${params.metadata ?: 'None'}
    Multi Mapping  : ${params.multi_samples ?: 'N/A'}
    Transcriptome  : ${params.transcriptome}
    Probe Set      : ${params.mode == 'multi' ? params.probe_set : 'N/A'}
    Skip FastQC    : ${params.skip_fastqc}
    Output Dir     : ${params.outdir}/${params.cohort}
    ==================================================
    """.stripIndent()

    // 1. Preprocess Inputs (BCL conversion if requested + Regex FASTQ staging)
    INPUT_PREPROCESSING(
        params.bcl_dir,
        params.bcl_samplesheet,
        params.reads,
        params.metadata
    )

    // 2. Optional FastQC on clean/staged reads
    fastqc_zip_ch = channel.empty()
    if (!params.skip_fastqc) {
        FASTQC(INPUT_PREPROCESSING.out.reads)
        fastqc_zip_ch = FASTQC.out.zip.flatten()
    }

    // 3. Prepare Reference Value Channels
    transcriptome_ch = channel.value(file(params.transcriptome, type: 'dir', checkIfExists: true))
    probe_set_ch     = params.mode == 'multi' 
        ? channel.value(file(params.probe_set, type: 'file', checkIfExists: true)) 
        : channel.empty()

    feature_ref_ch = params.feature_ref
        ? channel.value(file(params.feature_ref, type: 'file', checkIfExists: true))
        : channel.value([])

    // 4. Run Cell Ranger (Count or Multi)
    CELLRANGER_ALIGNMENT(
        INPUT_PREPROCESSING.out.reads,
        params.mode,
        transcriptome_ch,
        probe_set_ch,
        feature_ref_ch, 
        params.multi_samples
    )

    // 5. Aggregate all QC outputs (BCL reports + FastQC + Cell Ranger summaries) into MultiQC
    if (!params.skip_multiqc) {
        all_qc_ch = INPUT_PREPROCESSING.out.bcl_reports
            .mix(fastqc_zip_ch)
            .mix(CELLRANGER_ALIGNMENT.out.qc_metrics)
            .collect()

        MULTIQC(all_qc_ch)
    }
}