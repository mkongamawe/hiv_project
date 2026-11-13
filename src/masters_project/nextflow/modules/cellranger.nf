/*
 * This module runs cellranger on the samples
 * It has 3 processes:
 * STAGE_FASTQ - Renames non-compliant FASTQ files into
 *             - the format that cellranger requires
 *             - e.g.(_S1_L001_R1_001.fastq.gz)
 * CELLRANGER_COUNT - process for stanfard 3' single-cell workflow
 * CELLRANGER_MULTI - This is for Flex multiplexed workflows
 *                  - It dynamically constructs the multi_config.csv inside the task directory
 */

process STAGE_FASTQ {
    tag "$sample_id"
    label 'process_single'
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_S1_L*_R*.fastq.gz"), emit: staged_reads

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

        ln -s \$(realpath \$r1) ${sample_id}_S1_\${lane}_R1_001.fastq.gz
        ln -s \$(realpath \$r2) ${sample_id}_S1_\${lane}_R2_001.fastq.gz

        ((count++))
    done
    """

    stub:
    """
    touch ${sample_id}_S1_L001_R1_001.fastq.gz
    touch ${sample_id}_S1_L001_R2_001.fastq.gz
    """
}

process CELLRANGER_COUNT {
    tag "$sample_id"
    label "process_high"

    publishDir "${params.outdir}/${params.cohort}/cellranger", mode: 'copy'

    input:
    tuple val(sample_id), path(fastqs)
    path transcriptome

    output:
    path "${sample_id}/outs/*"                        , emit: outs
    path "${sample_id}/outs/web_summary.html"         , emit: html
    path "${sample_id}/outs/metrics_summary.csv"      , emit: metrics
    path "${sample_id}/outs/filtered_feature_bc_matrix", emit: matrix

    script:
    """
    cellranger count \\
        --id=${sample_id} \\
        --transcriptome=${transcriptome} \\
        --fastqs=. \\
        --create-bam=${params.create_bam} \\
        --sample=${sample_id} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()}
    """

    stub:
    """
    mkdir -p ${sample_id}/outs/filtered_feature_bc_matrix
    touch ${sample_id}/outs/web_summary.html
    touch ${sample_id}/outs/dummy_data.txt
    """
}

process CELLRANGER_MULTI {
    tag "$pool_id"
    label 'process_high'

    publishDir "${params.outdir}/{params.cohort}/cellranger_multi", mode: 'copy'

    input:
    tuple val(pool_id), path(fastqs), val(sample_rows)
    path transcriptome
    path probe_set

    output:
    path "${pool_id}/outs/*"                                                        , emit: outs
    path "${pool_id}/outs/per_sample_outs/*/web_summary.html"                       , emit: html
    path "${pool_id}/outs/per_sample_outs/*/metrics_summary.csv"                    , emit: metrics
    path "${pool_id}/outs/per_sample_outs/*/count/sample_filtered_feature_bc_matrix", emit: matrices

    script:
    // Format the sample block from the grouped metadata
    def samples_block = sample_rows.join('\n')
    """
    cat <<EOF > multi_config.csv
[gene-expression]
reference,\$(realpath ${transcriptome})
probe-set,\$(realpath ${probe_set})
create-bam,${params.create_bam}

[libraries]
fastq_id,fastqs,feature_types
${pool_id},\$(pwd),Gene Expression

[samples]
sample_id,probe_barcode_ids
${samples_block}
EOF

    cellranger multi \\
        --id=${pool_id} \\
        --csv=multi_config.csv \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()}
    """

    stub:
    """
    mkdir -p ${pool_id}/outs/per_sample_outs/Sample1/count/sample_filtered_feature_bc_matrix
    touch ${pool_id}/outs/per_sample_outs/Sample1/web_summary.html
    touch ${pool_id}/outs/per_sample_outs/Sample1/metrics_summary.csv
    """
}