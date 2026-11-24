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
    # Run name from a file name, e.g. SRR123_1.fastq.gz -> SRR123
    # (must match the regex in input_preprocessing.nf)
    run_of() { echo "\$1" | sed -E 's/(_[RI]?[123](_001)?)?\\.f(ast)?q\\.gz\$//'; }

    # Length of the first read in a gzipped FASTQ
    read_len() { gzip -cd "\$1" 2>/dev/null | awk 'NR==2 { print length(\$0); exit }'; }

    mapfile -t inputs < <(printf '%s\\n' ${reads} | sort)
    runs=\$(for f in "\${inputs[@]}"; do run_of "\$f"; done | sort -u)

    lane=1
    for run in \$runs; do
        # This run's files, shortest read first (sort -s keeps _1/_2 order on ties)
        mapfile -t ordered < <(
            for f in "\${inputs[@]}"; do
                [[ "\$(run_of "\$f")" == "\$run" ]] && echo "\$(read_len "\$f") \$f"
            done | sort -s -n -k1,1 | cut -d' ' -f2-
        )

        case \${#ordered[@]} in
            2) r1=\${ordered[0]}; r2=\${ordered[1]} ;;
            3) r1=\${ordered[1]}; r2=\${ordered[2]} ;;   # shortest = index read, skipped
            *) echo "ERROR: run \$run has \${#ordered[@]} FASTQ file(s): \${ordered[*]}" >&2; exit 1 ;;
        esac

        L=\$(printf "L%03d" \$lane)
        echo "\$run -> \$L : R1=\$r1 (\$(read_len "\$r1") bp)  R2=\$r2 (\$(read_len "\$r2") bp)"
        ln -s "\$(realpath "\$r1")" ${sample_id}_S1_\${L}_R1_001.fastq.gz
        ln -s "\$(realpath "\$r2")" ${sample_id}_S1_\${L}_R2_001.fastq.gz
        lane=\$((lane + 1))
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
        --localmem=${task.memory.toGiga() - 4}
    """

    stub:
    """
   mkdir -p ${sample_id}/outs/filtered_feature_bc_matrix
   touch ${sample_id}/outs/web_summary.html
   touch ${sample_id}/outs/metrics_summary.csv
    """
}

process CELLRANGER_MULTI {
    tag "$pool_id"
    label 'process_high'

    publishDir "${params.outdir}/${params.cohort}/cellranger_multi", mode: 'copy'

    input:
    tuple val(pool_id), path(fastqs), val(library_rows), val(sample_rows)
    path transcriptome
    path probe_set
    path feature_ref

    output:
    path "${pool_id}/outs/*"                                                        , emit: outs
    path "${pool_id}/outs/per_sample_outs/*/web_summary.html"                       , emit: html
    path "${pool_id}/outs/per_sample_outs/*/metrics_summary.csv"                    , emit: metrics
    path "${pool_id}/outs/per_sample_outs/*/count/sample_filtered_feature_bc_matrix", emit: matrices

    script:
    def samples_block   = sample_rows.join('\n')
    def libraries_block = library_rows.join('\n').replace('FASTQ_DIR', '$(pwd)')
    def feature_block   = feature_ref ? "[feature]\nreference,\$(realpath ${feature_ref})\n" : ''
    """
    cat <<EOF > multi_config.csv
[gene-expression]
reference,\$(realpath ${transcriptome})
probe-set,\$(realpath ${probe_set})
create-bam,${params.create_bam}

${feature_block}
[libraries]
fastq_id,fastqs,feature_types
${libraries_block}

[samples]
sample_id,probe_barcode_ids
${samples_block}
EOF

    cellranger multi \\
        --id=${pool_id} \\
        --csv=multi_config.csv \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga() - 4}
    """

    stub:
    """
    mkdir -p ${pool_id}/outs/per_sample_outs/Sample1/count/sample_filtered_feature_bc_matrix
    touch ${pool_id}/outs/per_sample_outs/Sample1/web_summary.html
    touch ${pool_id}/outs/per_sample_outs/Sample1/metrics_summary.csv
    """
}