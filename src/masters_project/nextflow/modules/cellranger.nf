/*
 * Cellranger requires strict naming convention.
 * The first process converts the names to the ones
 * required by cellranger for procesisng
 */

process STAGE_FASTQ {
    tag "$biosample"
    label 'process_single'
    
    input:
    tuple val(biosample), path(reads)

    output:
    tuple val(biosample), path("${biosample}_S1_L*_R*.fastq.gz"), emit: staged_reads

    script:
    """
    # Sort the files so reads 1 and 2 stay paired
    sorted_reads=(\$(echo "${reads}" | tr ' ' '\\n' | sort))

    count=1
    # Loop through the paired files at a time
    for (( i = 0; i<\${#sorted_reads[@]}; i+=2 )); do
        r1=\${sorted_reads[i]}
        r2=\${sorted_reads[i+1]}

        # Format a fake lane number
        lane=\$(printf "L%03d" \$count)

        ln -s \$(realpath \$r1) ${biosample}_S1_\${lane}_R1_001.fastq.gz
        ln -s \$(realpath \$r2) ${biosample}_S1_\${lane}_R2_001.fastq.gz

        ((count++))
    done
    """

    stub:
    """
    touch ${biosample}_S1_L001_R1_001.fastq.gz
    touch ${biosample}_S1_L001_R2_001.fastq.gz
    """
}

process CELLRANGER_COUNT {
    tag "$biosample"
    label "process_high"

    // Load the module
    module 'cellranger/10.0.0'

    publishDir "${params.outdir}/${params.cohort}/cellranger", mode: 'copy'

    input:
    tuple val(biosample), path(fastqs)
    path transcriptome

    output:
    path "${biosample}/outs/*"                              , emit: outs
    path "${biosample}/outs/web_summary.html"               , emit: html
    path "${biosample}/outs/filtered_feature_bc_matrix"     , emit: matrix

    script:
    """
    cellranger count \\
        --id=${biosample} \\
        --transcriptome=${transcriptome} \\
        --fastqs=. \\
        --create-bam=false \\
        --sample=${biosample} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()}
    """

    stub:
    """
    mkdir -p ${biosample}/outs/filtered_feature_bc_matrix
    touch ${biosample}/outs/web_summary.html
    touch ${biosample}/outs/dummy_data.txt

    echo "Dry run complete for ${biosample} using ${task.cpus} CPUs and ${task.memory.toGiga()} GB RAM"
    """
}