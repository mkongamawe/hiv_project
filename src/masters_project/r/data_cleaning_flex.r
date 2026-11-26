#!/usr/bin/env Rscript
# =============================================================================
# data_cleaning_flex.r — QC and filtering for 10x Flex (cellranger multi) data
#
# Input : <multi_dir>/<pool>/outs/per_sample_outs/<sample>/sample_filtered_feature_bc_matrix
#         (the `cellranger_multi/` folder published by the Nextflow multi mode)
# Output: <out_dir>/objects/<cohort>_qc.rds      cleaned, merged Seurat v5 object
#                                                (RNA + ADT if antibody capture)
#         <out_dir>/tables/qc/*.csv(.gz)          per-cell log, retention, cut-offs
#         <out_dir>/figures/qc/*.png              QC plots
#
# How Flex differs from 3' here:
#   * Probe-based on fixed cells: ambient RNA is low and there is no per-sample
#     raw matrix, so decontX (filtered matrix only) is used instead of SoupX.
#   * The probe set may target few/no mitochondrial genes, and fixation changes
#     mito %. If fewer than --min_mt_genes MT- genes are measured, the mito
#     filter is switched off and only reported.
#   * Antibody Capture (TotalSeq-C) is kept as an ADT assay; ADT UMIs are
#     reported and can be filtered with min_adt_counts in --thresholds.
#   * Several samples share one GEM well (pool). pool_id is stored, and used
#     as the default `batch`.
#
# Sample metadata: --metadata CSV keyed by sample_id. multi_samples.csv from
# the pipeline works as-is; add columns such as donor_id, condition and
# sort_population (B / T / Mono / Other) to it for sorted libraries.
#
# Example:
#   uvr run src/masters_project/r/data_cleaning_flex.r -- \
#     --multi_dir $BASE/results/flexv2_pbmc/cellranger_multi \
#     --metadata  $BASE/metadata/multi_samples.csv \
#     --cohort flexv2_pbmc \
#     --out_dir $BASE/results/flexv2_pbmc/seurat
# =============================================================================

suppressPackageStartupMessages(library(optparse))

script_dir <- local({
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f))) else getwd()
})
source(file.path(script_dir, "qc_utils.r"))

opts <- parse_args(OptionParser(option_list = list(
  make_option("--multi_dir", help = "Folder holding <pool>/outs/per_sample_outs [required]"),
  make_option("--metadata", default = NULL, help = "Sample metadata CSV keyed by sample_id"),
  make_option("--condition_col", default = NULL, help = "Metadata column copied to `condition`"),
  make_option("--donor_col", default = NULL, help = "Metadata column copied to `donor_id`"),
  make_option("--sort_col", default = NULL, help = "Metadata column copied to `sort_population`"),
  make_option("--cohort", default = "flex", help = "Cohort name, used in file names [%default]"),
  make_option("--out_dir", default = "results/seurat", help = "Output folder [%default]"),
  make_option("--thresholds", default = NULL, help = "Optional per-sort_population thresholds CSV"),
  make_option("--ambient", default = "decontx", help = "decontx | none [%default]"),
  make_option("--min_mt_genes", type = "integer", default = 5,
              help = "Minimum MT- genes measured to use % mito as a filter [%default]"),
  make_option("--dbr", type = "double", default = NULL, help = "Fixed doublet rate (default: scDblFinder's)"),
  make_option("--no_uncorrected", action = "store_true", default = FALSE,
              help = "Do not keep the pre-decontX counts as assay RNA_uncorrected"),
  make_option("--seed", type = "integer", default = 42)
)))
if (is.null(opts$multi_dir)) stop("--multi_dir is required")
if (!opts$ambient %in% c("decontx", "none")) {
  stop("--ambient must be decontx or none for Flex (no per-sample raw matrix for SoupX)")
}

set.seed(opts$seed)
options(future.globals.maxSize = 8 * 1024^3)
dirs <- make_out_dirs(opts$out_dir)

# Defaults for Flex PBMC. These are starting points: check the violins and
# the cut-offs table after the first run and tune via --thresholds.
default_thr <- list(
  min_counts = 500, min_features = 200,
  max_counts = NA, max_features = NA,
  max_percent_mt = 15, mt_cut_floor = 5,
  max_percent_hb = NA, min_complexity = 0.8,
  nmads_counts = 5, nmads_mt = 3, min_adt_counts = 0
)
thr_table <- load_thresholds(opts$thresholds, default_thr)

# -----------------------------------------------------------------------------
# 1. Sample table
# -----------------------------------------------------------------------------
mat_dirs <- Sys.glob(file.path(opts$multi_dir, "*", "outs", "per_sample_outs", "*",
                               "sample_filtered_feature_bc_matrix"))
if (!length(mat_dirs)) {
  stop("No <pool>/outs/per_sample_outs/<sample>/sample_filtered_feature_bc_matrix under ",
       opts$multi_dir)
}
discovered <- tibble(
  mat_dir = mat_dirs,
  sample_id = basename(dirname(mat_dirs)),
  pool_id = basename(dirname(dirname(dirname(dirname(mat_dirs)))))
)
if (anyDuplicated(discovered$sample_id)) {
  stop("The same sample_id appears in more than one pool: ",
       toString(unique(discovered$sample_id[duplicated(discovered$sample_id)])))
}

meta <- if (!is.null(opts$metadata)) {
  m <- read_csv(opts$metadata, show_col_types = FALSE, col_types = cols(.default = "c"))
  if (!"sample_id" %in% names(m)) stop("--metadata needs a sample_id column")
  m
} else {
  select(discovered, sample_id)
}

copy_col <- function(meta, from, to) {
  if (is.null(from)) return(meta)
  if (!from %in% names(meta)) stop("Column '", from, "' not found in sample metadata")
  meta[[to]] <- meta[[from]]
  meta
}
meta <- meta |>
  copy_col(opts$condition_col, "condition") |>
  copy_col(opts$donor_col, "donor_id") |>
  copy_col(opts$sort_col, "sort_population") |>
  select(-any_of("pool_id")) |>  # pool comes from the folder layout
  left_join(select(discovered, sample_id, pool_id, mat_dir), by = "sample_id")

missing <- meta$sample_id[is.na(meta$mat_dir)]
if (length(missing)) message("WARNING: no Cell Ranger output for: ", toString(missing), " — skipped")
unlisted <- setdiff(discovered$sample_id, meta$sample_id)
if (length(unlisted)) message("WARNING: Cell Ranger output not in metadata: ", toString(unlisted), " — skipped")
meta <- filter(meta, !is.na(mat_dir))

meta <- standardise_sample_meta(meta, chemistry = "10x_flex", cohort = opts$cohort)
write_csv(meta, file.path(dirs$tab, paste0(opts$cohort, "_samples.csv")))
log_msg("Processing ", nrow(meta), " samples from ", n_distinct(meta$pool_id), " pool(s)")
print(count(meta, pool_id, sort_population, condition))

# -----------------------------------------------------------------------------
# 2. Per-sample cleaning
# -----------------------------------------------------------------------------
results <- map(seq_len(nrow(meta)), function(i) {
  s <- meta[i, ]
  m <- read_10x_dir(s$mat_dir, prefix = s$sample_id)
  n_mt <- n_mt_genes(m$gex)
  use_mt <- n_mt >= opts$min_mt_genes
  if (!use_mt) {
    log_msg(s$sample_id, ": only ", n_mt, " MT- genes measured; % mito reported but not filtered")
  }
  if (!is.null(m$adt)) log_msg(s$sample_id, ": ", nrow(m$adt), " antibody features")
  clean_one_sample(
    filt = m$gex, raw = NULL, adt = m$adt,
    sample_meta = as.list(select(s, -mat_dir)),
    thr = thresholds_for(thr_table, s$sort_population),
    ambient = opts$ambient, seed = opts$seed, dbr = opts$dbr,
    keep_uncorrected = !opts$no_uncorrected, use_mt = use_mt
  )
})

# -----------------------------------------------------------------------------
# 3. Tables, plots, merged object
# -----------------------------------------------------------------------------
cell_log <- map(results, "cell_log") |> bind_rows()
cuts <- map(results, "cuts") |> bind_rows()
retention <- retention_table(cell_log)
if ("nCount_ADT" %in% names(cell_log)) {
  retention <- cell_log |>
    filter(qc_pass) |>
    group_by(sample_id) |>
    summarise(median_adt_counts = median(nCount_ADT), .groups = "drop") |>
    right_join(retention, by = "sample_id") |>
    relocate(median_adt_counts, .after = last_col())
}

write_csv(cell_log, file.path(dirs$tab, paste0(opts$cohort, "_cell_qc_log.csv.gz")))
write_csv(cuts, file.path(dirs$tab, paste0(opts$cohort, "_qc_cutoffs.csv")))
write_csv(retention, file.path(dirs$tab, paste0(opts$cohort, "_qc_retention.csv")))
write_csv(thr_table, file.path(dirs$tab, paste0(opts$cohort, "_qc_thresholds_input.csv")))
print(select(retention, sample_id, n_cellranger, n_doublets, n_retained, pct_retained, ambient_method))

qc_plots(cell_log, cuts, dirs$fig, opts$cohort)

if ("nCount_ADT" %in% names(cell_log)) {
  p_adt <- ggplot(cell_log, aes(nCount_RNA, nCount_ADT, colour = qc_pass)) +
    geom_point(size = 0.2, alpha = 0.4) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = c(`TRUE` = "#4C72B0", `FALSE` = "grey70")) +
    facet_wrap(~sample_id) +
    labs(title = "RNA vs ADT UMIs", colour = "Retained") +
    theme_bw(base_size = 9)
  save_plot(p_adt, dirs$fig, paste0(opts$cohort, "_rna_vs_adt.png"), w = 11, h = 8)
}

merged <- merge_samples(map(results, "seu"), project = opts$cohort)
obj_path <- file.path(dirs$obj, paste0(opts$cohort, "_qc.rds"))
saveRDS(merged, obj_path)
log_msg("Saved ", ncol(merged), " cells x ", nrow(merged), " genes (assays: ",
        toString(Assays(merged)), ") to ", obj_path)

writeLines(capture.output(sessionInfo()), file.path(dirs$tab, paste0(opts$cohort, "_sessionInfo.txt")))
