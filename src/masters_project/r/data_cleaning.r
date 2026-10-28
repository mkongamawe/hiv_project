#!/usr/bin/env Rscript

# Loading libraries
packages <- c("Seurat", "tidyverse", "harmony", "presto", "BiocManager")

bioc_packages <- c("scDblFinder")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg %in% bioc_packages) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", dependencies = TRUE)
      }
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}

set.seed(42)
options(future.globals.maxSize = 8 * 1024^3)

# Define directories
base_dir <- "results/training_ucsd/cellranger"
out_dir <- "results/training_ucsd"
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
obj_dir <- file.path(out_dir, "objects")
walk(
  c(fig_dir, tab_dir, obj_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------------------------
# 1. Read cell Ranger matrices, compute QC metrics, detect doublets per sample
# ------------------------------------------------------------------------------
