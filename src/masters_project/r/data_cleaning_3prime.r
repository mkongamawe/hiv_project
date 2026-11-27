#!/usr/bin/env Rscript
# =============================================================================
# data_cleaning_3prime.r — QC and filtering for 10x 3' (cellranger count) data
#
# Input : <cellranger_dir>/<sample>/outs/{filtered,raw}_feature_bc_matrix
#         (the `cellranger/` folder published by the Nextflow count mode)
# Output: <out_dir>/objects/<cohort>_qc.rds      cleaned, merged Seurat v5 object
#         <out_dir>/tables/qc/*.csv(.gz)          per-cell log, retention, cut-offs
#         <out_dir>/figures/qc/*.png              QC plots
#
# Steps per sample (see qc_utils.r): SoupX ambient correction (needs raw
# matrix; falls back to decontX) -> QC metrics -> hard floors -> scDblFinder ->
# MAD outlier flags capped by hard limits -> remove flagged cells + doublets.
#
# Sample metadata, in order of precedence:
#   --samplesheet  CSV with sample_id (+ optional outs_dir, donor_id,
#                  condition, sort_population, batch, any other columns)
#   --sra_table    SRA Run Selector table; runs are collapsed to BioSample
#                  (the Cell Ranger sample name used by the pipeline) and any
#                  column constant within a BioSample is kept.
#
# Example (UCSD cohort, stratified by viral load):
#   uvr run src/masters_project/r/data_cleaning_3prime.r -- \
#     --cellranger_dir /scratch/cmwagwabi/ucsd_cohort/results/ucsd_cohort/cellranger \
#     --sra_table data/ucsd_cohort/metadata/SraRunTable.csv \
#     --condition_col ViralLoad \
#     --cohort ucsd_cohort \
#     --out_dir /scratch/cmwagwabi/ucsd_cohort/results/ucsd_cohort/seurat
# =============================================================================

suppressPackageStartupMessages(library(optparse))

script_dir <- local({
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f))) else getwd()
})
source(file.path(script_dir, "qc_utils.r"))

opts <- parse_args(OptionParser(
  option_list = list(
    make_option(
      "--cellranger_dir",
      help = "Folder holding one <sample>/outs per sample [required]"
    ),
    make_option(
      "--samplesheet",
      default = NULL,
      help = "Sample metadata CSV (sample_id, ...)"
    ),
    make_option(
      "--sra_table",
      default = NULL,
      help = "SraRunTable.csv, used if no --samplesheet"
    ),
    make_option(
      "--sample_col",
      default = "BioSample",
      help = "SRA column naming the Cell Ranger sample [%default]"
    ),
    make_option(
      "--condition_col",
      default = NULL,
      help = "Metadata column copied to `condition` (e.g. ViralLoad)"
    ),
    make_option(
      "--donor_col",
      default = NULL,
      help = "Metadata column copied to `donor_id`"
    ),
    make_option(
      "--sort_col",
      default = NULL,
      help = "Metadata column copied to `sort_population` (B/T/Mono/Other)"
    ),
    make_option(
      "--cohort",
      default = "cohort",
      help = "Cohort name, used in file names [%default]"
    ),
    make_option(
      "--out_dir",
      default = "results/seurat",
      help = "Output folder [%default]"
    ),
    make_option(
      "--thresholds",
      default = NULL,
      help = "Optional per-sort_population thresholds CSV"
    ),
    make_option(
      "--ambient",
      default = "soupx",
      help = "soupx | decontx | none [%default]"
    ),
    make_option(
      "--dbr",
      type = "double",
      default = NULL,
      help = "Fixed doublet rate (default: scDblFinder's ~1%/1000 cells)"
    ),
    make_option(
      "--no_uncorrected",
      action = "store_true",
      default = FALSE,
      help = "Do not keep the pre-ambient-correction counts as assay RNA_uncorrected"
    ),
    make_option("--seed", type = "integer", default = 42)
  )
))
if (is.null(opts$cellranger_dir)) {
  stop("--cellranger_dir is required")
}

set.seed(opts$seed)
options(future.globals.maxSize = 8 * 1024^3)
dirs <- make_out_dirs(opts$out_dir)

# Defaults for 10x 3' PBMC. Override per sort_population with --thresholds.
default_thr <- list(
  min_counts = 500,
  min_features = 200,
  max_counts = NA,
  max_features = NA,
  max_percent_mt = 15,
  mt_cut_floor = 5,
  max_percent_hb = NA,
  min_complexity = 0.8,
  nmads_counts = 5,
  nmads_mt = 3,
  min_adt_counts = 0
)
thr_table <- load_thresholds(opts$thresholds, default_thr)

# -----------------------------------------------------------------------------
# 1. Sample table
# -----------------------------------------------------------------------------
discovered <- tibble(
  outs_dir = list.dirs(opts$cellranger_dir, recursive = FALSE)
) |>
  mutate(
    outs_dir = file.path(outs_dir, "outs"),
    sample_id = basename(dirname(outs_dir))
  ) |>
  filter(dir.exists(file.path(outs_dir, "filtered_feature_bc_matrix")))
if (!nrow(discovered)) {
  stop(
    "No <sample>/outs/filtered_feature_bc_matrix under ",
    opts$cellranger_dir
  )
}

collapse_sra <- function(path, sample_col) {
  sra <- read_csv(
    path,
    show_col_types = FALSE,
    col_types = cols(.default = "c")
  )
  if (!sample_col %in% names(sra)) {
    stop("Column '", sample_col, "' not in ", path)
  }
  # Keep only columns that are constant within a sample (drops Run, Bytes, ...).
  constant <- sra |>
    group_by(.data[[sample_col]]) |>
    summarise(across(everything(), ~ n_distinct(.x) == 1)) |>
    select(-all_of(sample_col)) |>
    summarise(across(everything(), all)) |>
    unlist()
  keep <- c(sample_col, names(constant)[constant])
  sra |>
    select(all_of(keep)) |>
    distinct() |>
    rename(sample_id = all_of(sample_col)) |>
    mutate(n_runs = as.character(table(sra[[sample_col]])[sample_id]))
}

meta <- if (!is.null(opts$samplesheet)) {
  read_csv(
    opts$samplesheet,
    show_col_types = FALSE,
    col_types = cols(.default = "c")
  )
} else if (!is.null(opts$sra_table)) {
  collapse_sra(opts$sra_table, opts$sample_col)
} else {
  select(discovered, sample_id)
}

copy_col <- function(meta, from, to) {
  if (is.null(from)) {
    return(meta)
  }
  if (!from %in% names(meta)) {
    stop("Column '", from, "' not found in sample metadata")
  }
  meta[[to]] <- meta[[from]]
  meta
}
meta <- meta |>
  copy_col(opts$condition_col, "condition") |>
  copy_col(opts$donor_col, "donor_id") |>
  copy_col(opts$sort_col, "sort_population")

if (!"outs_dir" %in% names(meta)) {
  meta <- left_join(
    meta,
    select(discovered, sample_id, outs_dir),
    by = "sample_id"
  )
}
missing <- meta$sample_id[is.na(meta$outs_dir)]
if (length(missing)) {
  message(
    "WARNING: no Cell Ranger output for: ",
    toString(missing),
    " — skipped"
  )
}
unlisted <- setdiff(discovered$sample_id, meta$sample_id)
if (length(unlisted)) {
  message(
    "WARNING: Cell Ranger output not in metadata: ",
    toString(unlisted),
    " — skipped"
  )
}
meta <- filter(meta, !is.na(outs_dir))

meta <- standardise_sample_meta(
  meta,
  chemistry = "10x_3prime",
  cohort = opts$cohort
)
write_csv(meta, file.path(dirs$tab, paste0(opts$cohort, "_samples.csv")))
log_msg("Processing ", nrow(meta), " samples")
if (all(is.na(meta$condition))) {
  message(
    "Note: no `condition` set — pass --condition_col (e.g. ViralLoad) for DE later"
  )
} else {
  print(count(meta, condition, sort_population))
}

# -----------------------------------------------------------------------------
# 2. Per-sample cleaning
# -----------------------------------------------------------------------------
results <- map(seq_len(nrow(meta)), function(i) {
  s <- meta[i, ]
  filt <- read_10x_dir(
    file.path(s$outs_dir, "filtered_feature_bc_matrix"),
    prefix = s$sample_id
  )
  raw_dir <- file.path(s$outs_dir, "raw_feature_bc_matrix")
  raw <- if (opts$ambient == "soupx" && dir.exists(raw_dir)) {
    read_10x_dir(raw_dir, prefix = s$sample_id)$gex
  }
  out <- clean_one_sample(
    filt = filt$gex,
    raw = raw,
    adt = filt$adt,
    sample_meta = as.list(select(s, -outs_dir)),
    thr = thresholds_for(thr_table, s$sort_population),
    ambient = opts$ambient,
    seed = opts$seed,
    dbr = opts$dbr,
    keep_uncorrected = !opts$no_uncorrected
  )
  rm(raw)
  gc(verbose = FALSE)
  out
})

# -----------------------------------------------------------------------------
# 3. Tables, plots, merged object
# -----------------------------------------------------------------------------
cell_log <- map(results, "cell_log") |> bind_rows()
cuts <- map(results, "cuts") |> bind_rows()
retention <- retention_table(cell_log)

write_csv(
  cell_log,
  file.path(dirs$tab, paste0(opts$cohort, "_cell_qc_log.csv.gz"))
)
write_csv(cuts, file.path(dirs$tab, paste0(opts$cohort, "_qc_cutoffs.csv")))
write_csv(
  retention,
  file.path(dirs$tab, paste0(opts$cohort, "_qc_retention.csv"))
)
write_csv(
  thr_table,
  file.path(dirs$tab, paste0(opts$cohort, "_qc_thresholds_input.csv"))
)
print(select(
  retention,
  sample_id,
  n_cellranger,
  n_doublets,
  n_retained,
  pct_retained,
  ambient_method
))

qc_plots(cell_log, cuts, dirs$fig, opts$cohort)

merged <- merge_samples(map(results, "seu"), project = opts$cohort)
obj_path <- file.path(dirs$obj, paste0(opts$cohort, "_qc.rds"))
saveRDS(merged, obj_path)
log_msg(
  "Saved ",
  ncol(merged),
  " cells x ",
  nrow(merged),
  " genes to ",
  obj_path
)

writeLines(
  capture.output(sessionInfo()),
  file.path(dirs$tab, paste0(opts$cohort, "_sessionInfo.txt"))
)
