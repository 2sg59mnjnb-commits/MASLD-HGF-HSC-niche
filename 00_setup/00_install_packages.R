
# Install all R packages required for the analysis pipeline
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c("limma", "edgeR", "WGCNA", "impute", "preprocessCore",
               "org.Hs.eg.db", "AnnotationDbi", "hgu133plus2.db",
               "clusterProfiler", "decoupleR", "dorothea",
               "GenomicRanges", "rtracklayer", "ComplexHeatmap", "UCell")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

cran_pkgs <- c("tidyverse", "data.table", "Matrix",
               "Seurat", "SeuratObject", "patchwork",
               "glmnet", "randomForest", "e1071", "kernlab", "xgboost", "caret",
               "pROC", "dcurves", "logistf", "brglm2", "clinfun", "diptest",
               "ggplot2", "ggrepel", "UpSetR", "VennDiagram", "reshape2",
               "TwoSampleMR", "GEOquery", "httr", "jsonlite", "future")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
github_pkgs <- c("jinworks/CellChat", "saeyslab/nichenetr",
                 "georgheinze/logistf")
for (pkg in github_pkgs) {
  pkg_name <- sub(".*/", "", pkg)
  if (!requireNamespace(pkg_name, quietly = TRUE))
    remotes::install_github(pkg)
}

cat("\n✅ All packages installed.\n")

