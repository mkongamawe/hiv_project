/*
 * The purpose of this module is to convert BCL files to
 * FASTQ files using Illumina bcl_convert
*/

process BCL_CONVERT {
    tag "$bcl_dir.simple.Name"
    label 'process_high'

    publishDir "${params.outdir}/${params.cohort}/fastq_demux", mode: 'copy', pattern: "*.fastq.gz"
    publishDir "${params.outdir}/${params.cohort}/bcl_reports", mode: 'copy', pattern: "Reports"

    input:
    path bcl_dir
    path samplesheet

    output:
    // Emit only valid FASTQs
    path "fastq_out/[!Undetermined]*.fastq.gz", emit: fastqs
    // Emit Illumina demultiplexing stats
    path "fastq_out/Reports"                  , emit: reports
    path "fastq_out/Logs"                     , emit: logs

    script:
    """
    bcl-convert \\
        --bcl-input-directory ${bcl_dir} \\
        --output-directory fastq_out \\
        --sample-sheet ${samplesheet} \\
        --no-lane-splitting false \\
        --bcl-num-conversion-threads ${task.cpus} \\
        --bcl-num-compression-threads ${task.cpus} \\
        --force
    """

    stub:
    """
    mkdir -p fastq_out/Reports fastq_out/Logs
    touch fastq_out/Sample1_S1_L001_R1_001.fastq.gz
    touch fastq_out/Sample1_S1_L001_R2_001.fastq.gz
    touch fastq_out/Undetermined_S0_L001_R1_001.fastq.gz
    touch fastq_out/Reports/Demultiplex_Stats.csv
    """
}