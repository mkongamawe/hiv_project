# =============================================================================
# qc_utils.r — shared helpers for the scRNA-seq data-cleaning scripts
#
# Sourced by data_cleaning_3prime.r and data_cleaning_flex.r. Nothing here
# reads or writes files on its own except the plotting/saving helpers.
#
# Per-sample QC order used by both scripts:
#   1. read Cell Ranger matrix (+ raw matrix for SoupX)
#   2. ambient RNA correction (SoupX for 3', decontX for Flex)
#   3. QC metrics (UMIs, genes, % mito, % ribo, % haemoglobin, complexity)
#   4. hard floors -> drop empty / broken barcodes before doublet calling
#   5. scDblFinder on the uncorrected counts
#   6. MAD-based outlier flags, with hard caps, per sample
#   7. drop flagged cells and doublets, keep a per-cell log of why
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(tidyverse)
  library(patchwork)
})

# Standard metadata columns every cleaned object carries, so 3' and Flex
# objects can be merged/integrated later without renaming anything.
standard_meta_cols <- c(
  "sample_id", "donor_id", "condition", "sort_population",
  "pool_id", "batch", "chemistry", "cohort"
)

# Recognised sort labels. Anything else is kept as-is but triggers a message,
# so typos ("Bcell" vs "B") are noticed early.
known_sort_populations <- c("B", "T", "Mono", "Other", "unsorted")

log_msg <- function(...) {
  message(format(Sys.time(), "[%H:%M:%S] "), ...)
}

# -----------------------------------------------------------------------------
# Reading
# -----------------------------------------------------------------------------

# Read a 10x matrix directory and split it by feature type.
# Returns list(gex = dgCMatrix, adt = dgCMatrix or NULL).
read_10x_dir <- function(path, prefix = NULL) {
  if (!dir.exists(path)) stop("Matrix directory not found: ", path)
  m <- Read10X(path, gene.column = 2, strip.suffix = FALSE)
  if (is.list(m)) {
    gex <- m[["Gene Expression"]]
    adt <- m[["Antibody Capture"]]
    if (is.null(gex)) stop("No 'Gene Expression' features in ", path)
  } else {
    gex <- m
    adt <- NULL
  }
  if (!is.null(prefix)) {
    colnames(gex) <- paste0(prefix, "_", colnames(gex))
    if (!is.null(adt)) colnames(adt) <- paste0(prefix, "_", colnames(adt))
  }
  list(gex = as(gex, "CsparseMatrix"), adt = adt)
}

# -----------------------------------------------------------------------------
# Thresholds
# -----------------------------------------------------------------------------
# Thresholds can differ per sorted population (monocytes carry more RNA and
# mitochondria than resting T/B cells). A thresholds CSV has one row per
# sort_population plus a "default" row; blank cells inherit from defaults.
#
#   sort_population,min_counts,min_features,max_percent_mt,...
#   default,500,200,15,...
#   Mono,,,20,...
#
# Columns:
#   min_counts, min_features     hard floors (also applied before doublets)
#   max_counts, max_features     hard caps (NA = none; doublets handled by
#                                scDblFinder)
#   max_percent_mt               hard cap on % mitochondrial
#   mt_cut_floor                 MAD-derived mito cut is never set below this
#   max_percent_hb               hard cap on % haemoglobin (NA = no filter)
#   min_complexity               log10(genes)/log10(UMIs) floor
#   nmads_counts, nmads_mt       MADs for the adaptive outlier flags
#   min_adt_counts               Flex/CITE only: ADT UMI floor (0 = none)
# -----------------------------------------------------------------------------

threshold_cols <- c(
  "min_counts", "min_features", "max_counts", "max_features",
  "max_percent_mt", "mt_cut_floor", "max_percent_hb", "min_complexity",
  "nmads_counts", "nmads_mt", "min_adt_counts"
)

load_thresholds <- function(path, defaults) {
  base <- as_tibble(c(list(sort_population = "default"), defaults))
  if (is.null(path) || !nzchar(path)) return(base)
  user <- read_csv(path, show_col_types = FALSE, col_types = cols(.default = "c"))
  if (!"sort_population" %in% names(user)) {
    stop("Thresholds CSV needs a 'sort_population' column")
  }
  unknown <- setdiff(names(user), c("sort_population", threshold_cols))
  if (length(unknown)) stop("Unknown threshold columns: ", toString(unknown))
  user <- user |> mutate(across(-sort_population, as.numeric))

  # Fill a user "default" row first, then populations inherit from it.
  user_default <- filter(user, sort_population == "default")
  if (nrow(user_default)) {
    for (col in intersect(names(user_default), threshold_cols)) {
      if (!is.na(user_default[[col]])) base[[col]] <- user_default[[col]]
    }
  }
  pops <- filter(user, sort_population != "default")
  if (!nrow(pops)) return(base)
  rows <- map(seq_len(nrow(pops)), function(i) {
    row <- base
    row$sort_population <- pops$sort_population[i]
    for (col in intersect(names(pops), threshold_cols)) {
      v <- pops[[col]][i]
      if (!is.na(v)) row[[col]] <- v
    }
    row
  })
  bind_rows(base, rows)
}

thresholds_for <- function(thr_table, population) {
  hit <- filter(thr_table, sort_population == population)
  if (!nrow(hit)) hit <- filter(thr_table, sort_population == "default")
  as.list(hit[1, ])
}

# -----------------------------------------------------------------------------
# Metadata helpers
# -----------------------------------------------------------------------------

# Make sure the standard columns exist and have sensible fallbacks.
standardise_sample_meta <- function(meta, chemistry, cohort) {
  if (!"sample_id" %in% names(meta)) stop("Sample metadata needs 'sample_id'")
  if (anyDuplicated(meta$sample_id)) {
    stop("Duplicated sample_id values: ",
         toString(unique(meta$sample_id[duplicated(meta$sample_id)])))
  }
  if (any(grepl("_", meta$sample_id))) {
    # Barcodes are prefixed "<sample_id>_<barcode>"; underscores are allowed,
    # but warn because they make parsing cell names back ambiguous.
    message("Note: sample_id contains '_' — use the metadata column, not ",
            "the cell name, to recover sample identity.")
  }
  fill <- function(col, value) {
    if (!col %in% names(meta)) meta[[col]] <<- value
    meta[[col]][is.na(meta[[col]]) | meta[[col]] == ""] <<- value[1]
  }
  if (!"donor_id" %in% names(meta)) meta$donor_id <- meta$sample_id
  if (!"condition" %in% names(meta)) meta$condition <- NA_character_
  fill("sort_population", "unsorted")
  if (!"pool_id" %in% names(meta)) meta$pool_id <- meta$sample_id
  if (!"batch" %in% names(meta)) meta$batch <- meta$pool_id
  meta$chemistry <- chemistry
  meta$cohort <- cohort

  odd <- setdiff(unique(meta$sort_population), known_sort_populations)
  if (length(odd)) {
    message("Note: non-standard sort_population values: ", toString(odd),
            " (expected one of ", toString(known_sort_populations), ")")
  }
  meta |>
    mutate(across(all_of(standard_meta_cols), as.character)) |>
    relocate(all_of(standard_meta_cols))
}

# -----------------------------------------------------------------------------
# Ambient RNA
# -----------------------------------------------------------------------------

# Fast throw-away clustering, only used to tell SoupX/decontX which cells are
# alike. Not the clustering used downstream.
quick_clusters <- function(counts, seed = 42, resolution = 0.8) {
  so <- CreateSeuratObject(counts, min.cells = 3)
  so <- NormalizeData(so, verbose = FALSE)
  so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
  so <- ScaleData(so, verbose = FALSE)
  npcs <- min(30, ncol(so) - 1)
  so <- RunPCA(so, npcs = npcs, verbose = FALSE, seed.use = seed)
  so <- FindNeighbors(so, dims = seq_len(min(20, npcs)), verbose = FALSE)
  so <- FindClusters(so, resolution = resolution, verbose = FALSE,
                     random.seed = seed)
  setNames(as.character(so$seurat_clusters), colnames(so))
}

# SoupX: needs the raw (all-barcode) matrix to estimate the soup profile.
# Falls back to decontX if automatic contamination estimation fails, which
# happens on low-heterogeneity samples such as a pure sorted population.
run_soupx <- function(raw, filt, clusters, seed = 42) {
  common <- intersect(rownames(raw), rownames(filt))
  raw <- raw[common, , drop = FALSE]
  filt <- filt[common, , drop = FALSE]
  sc <- SoupX::SoupChannel(tod = raw, toc = filt, calcSoupProfile = TRUE)
  sc <- SoupX::setClusters(sc, clusters[colnames(filt)])
  est <- tryCatch(
    SoupX::autoEstCont(sc, doPlot = FALSE, verbose = FALSE),
    error = function(e) {
      message("  SoupX autoEstCont failed (", conditionMessage(e),
              "); falling back to decontX")
      NULL
    }
  )
  if (is.null(est)) return(run_decontx(filt, clusters, seed))
  rho <- unname(est$metaData$rho[1])
  corrected <- SoupX::adjustCounts(est, roundToInt = TRUE, verbose = 0)
  list(
    counts = as(corrected, "CsparseMatrix"),
    contamination = setNames(rep(rho, ncol(filt)), colnames(filt)),
    method = "SoupX"
  )
}

# decontX: works from the filtered matrix alone; contamination is per cell.
run_decontx <- function(filt, clusters, seed = 42) {
  res <- decontX::decontX(
    x = filt, z = clusters[colnames(filt)], seed = seed, verbose = FALSE
  )
  corrected <- res$decontXcounts
  corrected@x <- round(corrected@x)
  corrected <- drop0(corrected)
  list(
    counts = as(corrected, "CsparseMatrix"),
    contamination = setNames(res$contamination, colnames(filt)),
    method = "decontX"
  )
}

correct_ambient <- function(filt, raw = NULL, method = c("soupx", "decontx", "none"),
                            seed = 42) {
  method <- match.arg(method)
  if (method == "none") {
    return(list(counts = filt,
                contamination = setNames(rep(NA_real_, ncol(filt)), colnames(filt)),
                method = "none"))
  }
  clusters <- quick_clusters(filt, seed = seed)
  if (method == "soupx") {
    if (is.null(raw)) {
      message("  No raw matrix available; using decontX instead of SoupX")
      return(run_decontx(filt, clusters, seed))
    }
    return(run_soupx(raw, filt, clusters, seed))
  }
  run_decontx(filt, clusters, seed)
}

# -----------------------------------------------------------------------------
# QC metrics
# -----------------------------------------------------------------------------

hb_genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM",
              "HBQ1", "HBZ")

add_qc_metrics <- function(seu) {
  genes <- rownames(seu)
  mt <- grep("^MT-", genes, value = TRUE)
  seu$percent_mt <- if (length(mt)) PercentageFeatureSet(seu, features = mt) else NA_real_
  seu$percent_ribo <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  hb <- intersect(hb_genes, genes)
  seu$percent_hb <- if (length(hb)) PercentageFeatureSet(seu, features = hb) else 0
  seu$log10_genes_per_umi <- log10(seu$nFeature_RNA) / log10(pmax(seu$nCount_RNA, 2))
  seu
}

# Number of mitochondrial genes the matrix actually measures. Flex probe sets
# may target few or none, in which case % mito is not a usable QC metric.
n_mt_genes <- function(counts) sum(grepl("^MT-", rownames(counts)))

# -----------------------------------------------------------------------------
# Doublets
# -----------------------------------------------------------------------------

# scDblFinder on uncorrected counts of one sample. `dbr` NULL lets
# scDblFinder pick the 10x-style rate (~1% per 1000 cells).
run_scdblfinder <- function(counts, seed = 42, dbr = NULL) {
  set.seed(seed)
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  args <- list(sce = sce, verbose = FALSE)
  if (!is.null(dbr)) args$dbr <- dbr
  sce <- do.call(scDblFinder::scDblFinder, args)
  tibble(
    cell = colnames(sce),
    doublet_score = sce$scDblFinder.score,
    doublet_class = as.character(sce$scDblFinder.class)
  )
}

# -----------------------------------------------------------------------------
# Filtering
# -----------------------------------------------------------------------------

# Lower/upper bounds at median +/- n MADs (on log scale for count metrics).
mad_bounds <- function(x, nmads, log = TRUE) {
  x <- x[is.finite(x)]
  if (log) x <- log1p(x)
  med <- median(x)
  dev <- mad(x)
  b <- c(lower = med - nmads * dev, upper = med + nmads * dev)
  if (log) b <- expm1(b)
  b
}

# Flag cells of ONE sample. Returns the metadata with qc_* flag columns and a
# `qc_reason` string, plus the numeric cut-offs actually used.
flag_low_quality <- function(meta, thr, use_mt = TRUE) {
  cnt_b <- mad_bounds(meta$nCount_RNA, thr$nmads_counts)
  feat_b <- mad_bounds(meta$nFeature_RNA, thr$nmads_counts)

  lo_counts <- max(thr$min_counts, cnt_b[["lower"]])
  lo_feat <- max(thr$min_features, feat_b[["lower"]])
  hi_counts <- min(c(thr$max_counts, cnt_b[["upper"]]), na.rm = TRUE)
  hi_feat <- min(c(thr$max_features, feat_b[["upper"]]), na.rm = TRUE)

  mt_cut <- NA_real_
  if (use_mt && any(is.finite(meta$percent_mt))) {
    mt_b <- mad_bounds(meta$percent_mt, thr$nmads_mt, log = FALSE)
    mt_cut <- min(thr$max_percent_mt, max(thr$mt_cut_floor, mt_b[["upper"]]))
  }

  meta <- meta |>
    mutate(
      qc_low_counts = nCount_RNA < lo_counts,
      qc_low_features = nFeature_RNA < lo_feat,
      qc_high_counts = nCount_RNA > hi_counts,
      qc_high_features = nFeature_RNA > hi_feat,
      qc_high_mt = if (is.na(mt_cut)) FALSE else percent_mt > mt_cut,
      qc_high_hb = if (is.na(thr$max_percent_hb)) FALSE else percent_hb > thr$max_percent_hb,
      qc_low_complexity = log10_genes_per_umi < thr$min_complexity
    )
  if ("nCount_ADT" %in% names(meta) && thr$min_adt_counts > 0) {
    meta$qc_low_adt <- meta$nCount_ADT < thr$min_adt_counts
  }

  cuts <- tibble(
    min_counts_used = lo_counts, min_features_used = lo_feat,
    max_counts_used = hi_counts, max_features_used = hi_feat,
    max_percent_mt_used = mt_cut
  )
  list(meta = meta, cuts = cuts)
}

# Collapse the flag columns into one readable reason per cell.
summarise_reasons <- function(meta) {
  flag_cols <- grep("^qc_", names(meta), value = TRUE)
  flag_cols <- setdiff(flag_cols, c("qc_reason", "qc_pass", "qc_floor_fail"))
  flags <- as.matrix(meta[, flag_cols, drop = FALSE])
  flags[is.na(flags)] <- FALSE
  reason <- apply(flags, 1, function(r) {
    hit <- sub("^qc_", "", flag_cols[r])
    if (length(hit)) paste(hit, collapse = ";") else ""
  })
  if ("doublet_class" %in% names(meta)) {
    is_dbl <- !is.na(meta$doublet_class) & meta$doublet_class == "doublet"
    reason <- ifelse(is_dbl, ifelse(reason == "", "doublet", paste0(reason, ";doublet")), reason)
  }
  meta$qc_reason <- ifelse(reason == "", "pass", reason)
  meta$qc_pass <- meta$qc_reason == "pass"
  meta
}

# -----------------------------------------------------------------------------
# Outputs
# -----------------------------------------------------------------------------

make_out_dirs <- function(out_dir) {
  dirs <- list(
    fig = file.path(out_dir, "figures", "qc"),
    tab = file.path(out_dir, "tables", "qc"),
    obj = file.path(out_dir, "objects")
  )
  walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  dirs
}

save_plot <- function(p, dir, name, w = 10, h = 7) {
  ggsave(file.path(dir, name), p, width = w, height = h, dpi = 300, bg = "white")
}

# Per-sample cell retention table: Cell Ranger cells -> passed floors ->
# passed QC -> passed doublets.
retention_table <- function(cell_log) {
  cell_log |>
    group_by(sample_id, sort_population, condition) |>
    summarise(
      n_cellranger = n(),
      n_pass_floors = sum(!qc_floor_fail),
      n_doublets = sum(doublet_class == "doublet", na.rm = TRUE),
      # Removed by the MAD/hard-cap flags (doublet-only and below-floor cells
      # are counted in their own columns).
      n_low_quality = sum(!qc_floor_fail & !qc_pass & qc_reason != "doublet"),
      n_retained = sum(qc_pass),
      pct_retained = round(100 * n_retained / n_cellranger, 1),
      median_counts = median(nCount_RNA[qc_pass]),
      median_features = median(nFeature_RNA[qc_pass]),
      median_percent_mt = median(percent_mt[qc_pass]),
      contamination = median(ambient_contamination, na.rm = TRUE),
      ambient_method = first(ambient_method),
      .groups = "drop"
    )
}

qc_plots <- function(cell_log, cuts, fig_dir, prefix) {
  cell_log <- cell_log |>
    mutate(status = factor(ifelse(qc_pass, "retained", "removed"),
                           c("retained", "removed")))
  has_mt <- any(is.finite(cell_log$percent_mt))

  # Violins of the main metrics, removed cells in grey on top.
  metrics <- c("nCount_RNA", "nFeature_RNA", if (has_mt) "percent_mt",
               "percent_ribo", "percent_hb", "log10_genes_per_umi")
  long <- cell_log |>
    select(sample_id, status, all_of(metrics)) |>
    pivot_longer(all_of(metrics), names_to = "metric")
  p_vln <- ggplot(long, aes(sample_id, value, fill = status)) +
    geom_violin(scale = "width", linewidth = 0.2, position = position_dodge(0.8)) +
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c(retained = "#4C72B0", removed = "grey70")) +
    scale_y_continuous(trans = "log1p") +
    labs(x = NULL, y = NULL, fill = NULL,
         title = "QC metrics per sample (log1p scale)") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_vln, fig_dir, paste0(prefix, "_qc_violins.png"), w = 12, h = 10)

  # UMIs vs genes, coloured by mito, one panel per sample, with cut lines.
  colour_by <- if (has_mt) "percent_mt" else "log10_genes_per_umi"
  p_sc <- ggplot(cell_log, aes(nCount_RNA, nFeature_RNA, colour = .data[[colour_by]])) +
    geom_point(size = 0.2, alpha = 0.5) +
    geom_vline(data = cuts, aes(xintercept = min_counts_used), linetype = 2) +
    geom_hline(data = cuts, aes(yintercept = min_features_used), linetype = 2) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_viridis_c(option = "plasma") +
    facet_wrap(~sample_id) +
    labs(title = "UMIs vs genes (dashed = lower cut-offs used)") +
    theme_bw(base_size = 9)
  save_plot(p_sc, fig_dir, paste0(prefix, "_counts_vs_genes.png"), w = 12, h = 10)

  # Doublet scores.
  p_dbl <- cell_log |>
    filter(!is.na(doublet_score)) |>
    ggplot(aes(doublet_score, fill = doublet_class)) +
    geom_histogram(bins = 50) +
    facet_wrap(~sample_id, scales = "free_y") +
    scale_fill_manual(values = c(singlet = "#4C72B0", doublet = "#DD8452")) +
    labs(title = "scDblFinder scores") +
    theme_bw(base_size = 9)
  save_plot(p_dbl, fig_dir, paste0(prefix, "_doublet_scores.png"), w = 12, h = 8)

  # Why cells were removed.
  p_why <- cell_log |>
    mutate(reason = ifelse(qc_floor_fail, "below_floors", qc_reason)) |>
    count(sample_id, reason) |>
    ggplot(aes(sample_id, n, fill = reason)) +
    geom_col() +
    labs(x = NULL, y = "Cells", title = "Cell fate per sample") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_why, fig_dir, paste0(prefix, "_cell_fate.png"), w = 11, h = 6)

  # Ambient contamination.
  if (any(is.finite(cell_log$ambient_contamination))) {
    p_amb <- ggplot(cell_log, aes(sample_id, ambient_contamination)) +
      geom_boxplot(outlier.size = 0.2) +
      labs(x = NULL, y = "Estimated contamination fraction",
           title = "Ambient RNA contamination") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(p_amb, fig_dir, paste0(prefix, "_ambient_contamination.png"), w = 9, h = 5)
  }
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# One sample, end to end
# -----------------------------------------------------------------------------
# filt / raw / adt are matrices with cells already prefixed by sample_id.
# Returns list(seu = cleaned Seurat object or NULL, cell_log, cuts).
clean_one_sample <- function(filt, raw = NULL, adt = NULL, sample_meta, thr,
                             ambient = "soupx", seed = 42, dbr = NULL,
                             keep_uncorrected = TRUE, use_mt = TRUE) {
  sid <- sample_meta$sample_id
  log_msg(sid, ": ", ncol(filt), " Cell Ranger cells, ambient = ", ambient)

  amb <- correct_ambient(filt, raw = raw, method = ambient, seed = seed)
  log_msg(sid, ": ambient correction via ", amb$method)

  seu <- CreateSeuratObject(amb$counts, project = sid, min.cells = 0, min.features = 0)
  if (keep_uncorrected && amb$method != "none") {
    seu[["RNA_uncorrected"]] <- CreateAssay5Object(counts = filt[rownames(amb$counts), colnames(seu)])
  }
  if (!is.null(adt)) {
    seu[["ADT"]] <- CreateAssay5Object(counts = adt[, colnames(seu), drop = FALSE])
    iso <- grep("isotype|IgG", rownames(adt), ignore.case = TRUE, value = TRUE)
    if (length(iso)) {
      seu$percent_isotype <- 100 * colSums(adt[iso, colnames(seu), drop = FALSE]) /
        pmax(seu$nCount_ADT, 1)
    }
  }
  seu <- add_qc_metrics(seu)
  seu$ambient_contamination <- unname(amb$contamination[colnames(seu)])
  seu$ambient_method <- amb$method
  for (col in names(sample_meta)) seu@meta.data[[col]] <- as.character(sample_meta[[col]])

  meta <- seu@meta.data |> rownames_to_column("cell") |> as_tibble()

  # Hard floors first: these barcodes are not cells worth simulating doublets from.
  meta$qc_floor_fail <- meta$nCount_RNA < thr$min_counts | meta$nFeature_RNA < thr$min_features
  keep_for_dbl <- meta$cell[!meta$qc_floor_fail]
  if (length(keep_for_dbl) < 50) {
    warning(sid, ": only ", length(keep_for_dbl), " cells pass the hard floors; skipping")
    meta$doublet_score <- NA_real_
    meta$doublet_class <- NA_character_
    meta$qc_reason <- "below_floors"
    meta$qc_pass <- FALSE
    return(list(seu = NULL, cell_log = meta, cuts = tibble(sample_id = sid)))
  }

  dbl <- run_scdblfinder(filt[, keep_for_dbl], seed = seed, dbr = dbr)
  meta <- left_join(meta, dbl, by = "cell")
  log_msg(sid, ": ", sum(dbl$doublet_class == "doublet"), " doublets called")

  # MAD cut-offs are computed on cells above the floors only, so empty-ish
  # barcodes do not drag the median down.
  cell_order <- meta$cell
  flagged <- flag_low_quality(filter(meta, !qc_floor_fail), thr, use_mt = use_mt)
  meta <- bind_rows(flagged$meta, filter(meta, qc_floor_fail))
  meta <- meta[match(cell_order, meta$cell), ]
  meta <- summarise_reasons(meta)
  meta$qc_reason[meta$qc_floor_fail] <- "below_floors"
  meta$qc_pass[meta$qc_floor_fail] <- FALSE

  keep <- meta$cell[meta$qc_pass]
  log_msg(sid, ": retained ", length(keep), " / ", nrow(meta), " cells")

  seu <- subset(seu, cells = keep)
  add <- meta |>
    filter(qc_pass) |>
    select(cell, doublet_score, doublet_class) |>
    column_to_rownames("cell")
  seu <- AddMetaData(seu, add)

  list(seu = seu, cell_log = meta,
       cuts = mutate(flagged$cuts, sample_id = sid, .before = 1))
}

# Merge per-sample objects into one, layers joined, ready for integration.
merge_samples <- function(objs, project) {
  objs <- compact(objs)
  if (!length(objs)) stop("No samples survived QC")
  merged <- if (length(objs) == 1) objs[[1]] else merge(objs[[1]], objs[-1], project = project)
  for (a in Assays(merged)) merged[[a]] <- JoinLayers(merged[[a]])
  DefaultAssay(merged) <- "RNA"
  merged$orig.ident <- merged$sample_id
  merged
}
