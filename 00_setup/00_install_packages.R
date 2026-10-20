# ============================================================================
# 00_install_packages.R -- install the R packages used in this repository
#
# Run once from the repository root:
#     Rscript 00_setup/00_install_packages.R
#
# Bioconductor packages are installed through BiocManager. Two packages that
# are not on CRAN/Bioconductor (CellChat) are installed from GitHub and are
# only needed for 03_spatial/22_CellChat_spatial.R.
# ============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))     install.packages("remotes")

cran <- c(
  # data handling / plotting
  "tidyverse", "data.table", "reshape2", "patchwork", "ggrepel", "scales",
  # modelling
  "glmnet", "randomForest", "e1071", "xgboost", "caret", "pROC", "dcurves",
  "logistf", "pwr", "simr", "lme4", "lmerTest", "sandwich", "lmtest",
  "nnls", "igraph", "ggraph",
  # single-cell / spatial
  "Seurat", "SeuratObject", "Matrix", "RANN", "harmony",
  # enrichment / annotation
  "msigdbr", "GSVA", "metafor", "clinfun", "diptest",
  # web
  "httr", "jsonlite"
)

bioc <- c(
  "limma", "edgeR", "WGCNA", "impute", "preprocessCore",
  "org.Hs.eg.db", "AnnotationDbi", "clusterProfiler",
  "decoupleR", "dorothea", "JASPAR2020", "TFBSTools", "Biostrings",
  "GenomicRanges", "UCell", "celda",
  "GEOquery", "hgu133plus2.db", "illuminaHumanv4.db"
)

for (p in cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing CRAN package: ", p)
    try(install.packages(p), silent = TRUE)
  }
}

for (p in bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing Bioconductor package: ", p)
    try(BiocManager::install(p, ask = FALSE, update = FALSE), silent = TRUE)
  }
}

# CellChat is GitHub-only and is required solely for the spatially constrained
# communication analysis. Everything else runs without it.
if (!requireNamespace("CellChat", quietly = TRUE)) {
  message("Installing CellChat from GitHub (optional)")
  try(remotes::install_github("jinworks/CellChat", upgrade = "never"),
      silent = TRUE)
}

# ---- Report ----------------------------------------------------------------
want <- c(cran, bioc, "CellChat")
status <- data.frame(
  package = want,
  installed = vapply(want, function(p) requireNamespace(p, quietly = TRUE), logical(1))
)
cat(sprintf("\n%d / %d packages available\n",
            sum(status$installed), nrow(status)))
missing <- status$package[!status$installed]
if (length(missing)) {
  cat("Still missing (install manually):\n  ", paste(missing, collapse = ", "), "\n")
}
