/*
 * The purpose of this subworkflow is to route the FASTQ files
 * to either CELLRANGER_MULTI or CELLRANGER_COUNT
*/

include { CELLRANGER_COUNT; CELLRANGER_MULTI } from '../modules/cellranger.nf'

workflow CELLRANGER_ALIGNMENT {
    take:
    reads_ch         // tuple(id, [fastqs])
    mode             // 'count' or 'multi'
    transcriptome_ch // value channel (dir)
    probe_set_ch     // value channel (file) or empty
    multi_samples    // path to CSV (pool_id,sample_id,probe_barcode_ids)

    main:
    qc_metrics_ch = channel.empty()
    matrices_ch   = channel.empty()

    if (mode == 'count') {
        CELLRANGER_COUNT(reads_ch, transcriptome_ch)

        qc_metrics_ch = CELLRANGER_COUNT.out.metrics.mix(CELLRANGER_COUNT.out.html)
        metrices_ch   = CELLRANGER_COUNT.out.matrix
    
    } else if (mode == 'multi') {

        // Parse the multi sample mapping sheet
        pool_map_ch = channel
            .fromPath(multi_samples, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> 
                def sample_line = "${row.sample_id}, ${row.probe_barcode_ids}"
                tuple(row.pool_id, sample_line)
            }
            .groupTuple() // -> tuple(pool_id, ["Sample1,BC001", "Sample2,BC002", ...])

        // Join the FASTQ files for each pool_id
        multi_input_ch = reads_ch.join(pool_map_ch)

        CELLRANGER_MULTI(multi_input_ch, transcriptome_ch, probe_set_ch)

        qc_metrics_ch = CELLRANGER_MULTI.out.metrics.mix(CELLRANGER_MULTI.out.html).flatten()
        matrices_ch   = CELLRANGER_MULTI.out.matrices

    } else {
        error "Invalid mode '${mode}'. Please use --mode 'count' or --mode 'multi'."
    }

    emit:
    qc_metrics = qc_metrics_ch
    matrices   = matrices_ch

}