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
    feature_ref_ch   // value channel (file) or [] when no Antibody Capture library
    multi_samples    // path to CSV (pool_id,sample_id,probe_barcode_ids)

    main:
    qc_metrics_ch = channel.empty()
    matrices_ch   = channel.empty()

    if (mode == 'count') {
        CELLRANGER_COUNT(reads_ch, transcriptome_ch)

        qc_metrics_ch = CELLRANGER_COUNT.out.metrics.mix(CELLRANGER_COUNT.out.html)
        matrices_ch   = CELLRANGER_COUNT.out.matrix
    
    } else if (mode == 'multi') {

        // Parse the multi sample mapping sheet
        pool_map_ch = channel
            .fromPath(multi_samples, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> 
                def sample_line = "${row.sample_id.trim()},${row.probe_barcode_ids.trim()}"
                tuple(row.pool_id, sample_line)
            }
            .groupTuple() // -> tuple(pool_id, ["Sample1,BC001", "Sample2,BC002", ...])

        // Join the FASTQ files for each pool_id
        // Merge <pool>_gex and <pool>_ab FASTQ sets into one pool with a library row each
        def lib_types = [gex: 'Gene Expression', ab: 'Antibody Capture']
        libs_ch = reads_ch
            .map { fastq_id, fqs ->
                def m = (fastq_id =~ /^(.+)_(gex|ab)$/)
                if (!m.matches()) error "Multi mode: FASTQ id '${fastq_id}' must end in _gex or _ab"
                tuple(m.group(1), "${fastq_id},FASTQ_DIR,${lib_types[m.group(2)]}".toString(), fqs)
            }
            .groupTuple()
            .map { pool, lib_rows, fqs -> tuple(pool, fqs.flatten(), lib_rows) }

        // -> tuple(pool_id, [fastqs], [library rows], [sample rows])
        multi_input_ch = libs_ch.join(pool_map_ch, failOnMismatch: true)

        CELLRANGER_MULTI(multi_input_ch, transcriptome_ch, probe_set_ch, feature_ref_ch)

        qc_metrics_ch = CELLRANGER_MULTI.out.metrics.mix(CELLRANGER_MULTI.out.html).flatten()
        matrices_ch   = CELLRANGER_MULTI.out.matrices

    } else {
        error "Invalid mode '${mode}'. Please use --mode 'count' or --mode 'multi'."
    }

    emit:
    qc_metrics = qc_metrics_ch
    matrices   = matrices_ch

}