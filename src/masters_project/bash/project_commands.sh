# To download the ucsd_data from NCBI run:

THREADS=8 bash src/masters_project/bash/sra_download.sh \
    data/ucsd_cohort/metadata/SRR_Acc_List.txt \
    /scratch/cmwagwabi/ucsd_cohort/fastq

# To run the cellranger stub for multi
nextflow run main.nf -profile laptop -stub --mode multi \
  --reads '<fastq folder>/*/*.fastq.gz' --multi_samples <multi_samples.csv> \
  --feature_ref <feature_reference.csv> --probe_set <probe set csv> --transcriptome <ref dir>

nextflow run src/masters_project/nextflow/main.nf -profile laptop -stub --mode multi \
  --reads "$BASE/fastq/16k_Human_PBMCs_TotalseqC_4plex_fastqs/*/*.fastq.gz" \
  --multi_samples $BASE/metadata/multi_samples.csv \
  --feature_ref $BASE/metadata/feature_reference.csv \
  --probe_set <probe set csv> --transcriptome <ref dir> \
  --outdir $BASE/stub_results -work-dir $BASE/stub_work

# To run the actual cellranger multi for the data
nextflow run src/masters_project/nextflow/main.nf -profile hpc \
  --mode multi \
  --cohort flexv2_pbmc \
  --reads "$BASE/fastq/16k_Human_PBMCs_TotalseqC_4plex_fastqs/*/*.fastq.gz" \
  --multi_samples $BASE/metadata/multi_samples.csv \
  --feature_ref $BASE/metadata/feature_reference.csv \
  --probe_set /path/to/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv \
  --transcriptome /path/to/refdata-gex-GRCh38-2024-A \
  --outdir $BASE/results \
  -work-dir $BASE/work