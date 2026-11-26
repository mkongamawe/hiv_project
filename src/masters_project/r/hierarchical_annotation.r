#!/usr/bin/env Rscript

# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(harmony)
})
set.seed(42)

out_dir <- "results/training_ucsd"
fig_dir <- file.path(out_dir, "figures", "annotation")
tab_dir <- file.path(out_dir, "tables", "annotation")
obj_dir <- file.path(out_dir, "objects")
walk(c(fig_dir, tab_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w = 9, h = 7) {
  ggsave(file.path(fig_dir, name), p, width = w, height = h, dpi = 300,
         bg = "white")
}

pbmc <- readRDS(file.path(obj_dir, "pbmc_integrated.rds"))
cluster_col <- "RNA_snn_res.0.5"
Idents(pbmc) <- cluster_col