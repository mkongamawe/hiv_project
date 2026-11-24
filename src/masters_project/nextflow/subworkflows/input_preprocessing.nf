/*
 * The purpose of this subworkflow is to handle upstream complexities
 * 1. It runs BCL_CONVERT if --bcl-dir is provided otherwise utilises --reads
 * 2. It contains a regex routher that inspects every FASTQ filename so that:
 *    a) Compliant files `Sample_S1_L001_R1_001.fastq.gz` bypass the STAGE_FASTQ process
 *    b) Non-compliant files `SRR12623961_1.fastq.gz` are first pooled if need be
         and then sent to through STAGE_FASTQ
*/

include { BCL_CONVERT } from '../modules/bcl_convert.nf'
include { STAGE_FASTQ } from '../modules/cellranger.nf'

workflow INPUT_PREPROCESSING {
    take:
    bcl_dir        // Can be a path or null
    samplesheet    // Can be a path or null
    reads_pattern  // Can be a path or null
    metadata_csv   // Can be a path or null

    main:
    bcl_reports_ch = Channel.empty()

    // 1. Get the raw FASTQ files - could be from BCL_CONVERT of from files
    if (bcl_dir) {
        if (!samplesheet) error "--bcl_samplesheet is required when --bcl_dir is set"
        bcl_ch   = channel.fromPath(bcl_dir, checkIfExists: true, type: 'dir')
        sheet_ch = channel.fromPath(samplesheet, checkIfExists: true, type: 'file')

        BCL_CONVERT(bcl_ch, sheet_ch)
        raw_fastqs_ch  = BCL_CONVERT.out.fastqs
            .flatten()
            .filter { fq -> !fq.name.startsWith('Undetermined') }
        bcl_reports_ch = BCL_CONVERT.out.reports 
    } else {
        raw_fastqs_ch = channel.fromPath(reads_pattern, checkIfExists: true, type: 'file')
    }

    // 2. Use regex to check if the naming convention is appropriate for cellranger
    // Pattern: <SampleName>_S<SampleNum>_L<LaneNum>_<R1/R2/I1/I2>_001.fastq.gz
    def illumina_regex = ~/^(.+)_S\d+_L\d{3}_[RI][123]_001\.f(ast)?q\.gz$/

    branched_fastqs = raw_fastqs_ch.branch { file -> 
        compliant: file.name ==~ illumina_regex
        non_compliant: true
    }

    // 3A. Handle compliant FASTQs
    compliant_grouped_ch = branched_fastqs.compliant
        .map { file -> 
            def matcher = (file.name =~ illumina_regex)
            def sample_id = matcher[0][1]
            tuple(sample_id, file)
        }
        .groupTuple(sort: true)

    // 3B. Handle non-compliant FASTQs
    non_compliant_by_prefix = branched_fastqs.non_compliant
        .map { file ->
            def prefix = file.name.replaceAll(/(_[RI]?[123](_001)?)?\.f(ast)?q\.gz$/, '')
            tuple(prefix, file)
        }
        .groupTuple(sort: true)

    // If SRA metadata is provided
    if (metadata_csv) {
        meta_ch = channel
            .fromPath(metadata_csv, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> tuple(row.Run, row.BioSample) } // Change the column here to group by e.g. Experiment

        // remainder: true keeps unmatched entries (the missing side is null) so they can be reported
        joined_ch = meta_ch
            .join(non_compliant_by_prefix, remainder: true)
            .branch { row ->
                no_fastq:    row[2] == null
                no_metadata: row[1] == null
                matched:     true
            }

        joined_ch.no_fastq.view    { row -> "WARNING: ${row[0]} (${row[1]}) is in the metadata but has no FASTQ files - skipped" }
        joined_ch.no_metadata.view { row -> "WARNING: ${row[0]} has FASTQ files but is not in the metadata - skipped" }

        to_stage_ch = joined_ch.matched
            .map { srr, biosample, reads -> tuple(biosample, reads) }
            .groupTuple()
            .map { biosample, nested_reads -> tuple(biosample, nested_reads.flatten()) }
    } else {
        to_stage_ch = non_compliant_by_prefix
    }

    // Rename non-compliant files via symlinks
    STAGE_FASTQ(to_stage_ch)

    // 4. Merge compliant and newly-staged reads into one stream
    clean_reads_ch = compliant_grouped_ch.mix(STAGE_FASTQ.out.staged_reads)

    emit:
    reads       = clean_reads_ch // Format: tuple(sample_or_pool_id, [fastq_files])
    bcl_reports = bcl_reports_ch
}