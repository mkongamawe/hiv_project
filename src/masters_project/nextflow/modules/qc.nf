/*
 * The purpose of this module is to generate BCL states,  
 * FASTQC reports and cell ranger metrics into a
 * single interactive HTML dashboard
*/

process FASTQC {
    tag "$sample_id"
    label 'process_low'

    publishDir "${params.outdir}/${params.cohort}/qc/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.{html,zip}", emit: zip

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --quiet \\
        ${reads}
    """

    stub:
    """
    touch ${sample_id}_fastqc.html
    touch ${sample_id}_fastqc.zip
    """
}

process MULTIQC {
    tag "$params.cohort"
    label 'process_single'

    publishDir "${params.outdir}/${params.cohort}/qc/multiqc", mode: 'copy'

    input:
    path qc_files, stageAs: "?/*"

    output:
    // With --title, MultiQC prefixes the file names with the title
    // (e.g. <title>_multiqc_report.html), so match any prefix
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data

    script:
    """
    multiqc . \\
        --title "${params.cohort} scRNA-seq QC Report" \\
        --force
    """

    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_report.html
    """
}