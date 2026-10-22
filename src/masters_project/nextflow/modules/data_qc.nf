// A script that runs both fastqc and multiqc on samples

process runFastQc {
   tag "$sample_id"
   label 'process_low'

   // Using pre-built containers
   container 'https://depot.galaxyproject.org/singularity/fastqc:0.13.0--hdfd78af_0'
   conda     'bioconda::fastqc=0.13.1'

   publishDir "${params.outdir}/fastqc", mode: 'copy', pattern: "*.{html,zip}"

   input:
   tuple val(sample_id), path(reads)

   output:
   path "*.{html,zip}", emit: qc_files

   script:
   """
   fastqc \\
        --quiet \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${reads}
   """
}

process runMultiQc {
    tag "Aggregating QC across all samples"
    label 'process_single'

    container 'https://depot.galaxyproject.org/singularity/multiqc:1.35--pyhdfd78af_2'
    conda     'bioconda::multiqc=1.35'

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data"       , emit: data

    script:
    """
    multiqc . \\
        --force \\
        --filename multiqc_report.html
    """
}