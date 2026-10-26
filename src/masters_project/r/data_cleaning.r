#!usr/bin/env Rscript

# Loading libraries
packages <- c("Seurat", "tidyverse", "harmony")

for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}
