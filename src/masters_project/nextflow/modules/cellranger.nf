/*
 * Cellranger requires strict naming convention.
 * The first process converts the names to the ones
 * required by cellranger for procesisng
 */

process STAGE_FASTQ {
    tag "$sample_id"
    label 'process_single'
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_S1_L001_R*.fastq.gz"), emit: staged_reads

    script:
    """
    ln -s \$(realpath ${reads[0]}) ${sample_id}_S1_L001_R1_001.fastq.gz
    ln -s \$(realpath ${reads[1]}) ${sample_id}_S1_L001_R1_002.fastq.gz
    """

    stub:
    """
    touch ${sample_id}_S1_L001_R1_001.fastq.gz
    touch ${sample_id}_S1_L001_R2_001.fastq.gz
    """
}

process CELLRANGER_COUNT {
    tag "sample_id"
    label "process_high"

    // Load the module
    module 'cellranger/10.0.0'

    publishDir "${params.outdir}/{params.cohort}/cellranger", mode: 'copy'

    input:
    tuple val(sample_id), path(fastqs)
    path transcriptome

    output:
    path "${sample_id}/outs/*"                              , emit: outs
    path "${sample_id}/outs/web_summary.html"               , emit: html
    path "${sample_id}/outs/filtered_feature_bc_matrix"     , emit: matrix

    script:
    """
    cellranger count \\
        --id=${sample_id} \\
        --transcriptome=${transcriptome} \\
        --fastqs=. \\
        --sample=${sample_id} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()}
    """

    stub:
    """
    mkdir -p ${sample_id}/outs/filtered_feature_bc_matrix
    touch ${sample_id}/outs/web_summary.html
    touch ${sample_id}/outs/dummy_data.txt

    echo "Dry run complete for ${sample_id} using ${task.cpus} CPUs and ${task.memory.toGiga()} GB RAM"
    """
}