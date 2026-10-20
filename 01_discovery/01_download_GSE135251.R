# ============================================================================
# 01_download_GSE135251.R
#
# Downloads the training cohort: 216 per-sample count files from GSE135251.
# Also saves the sample phenotype table, which carries the METAVIR fibrosis
# stage and the disease label used to define groups.
#
# Output
#   data/GSE135251_counts/*.counts.txt.gz   216 files
#   data/GSE135251_pheno.csv                sample metadata
#
# Run time: 10-30 min depending on connection. Re-running skips existing files.
# ============================================================================

source("R/config.R")

suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
})

options(timeout = 3600)
count_dir <- file.path(DATA_DIR, "GSE135251_counts")
dir.create(count_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Series-level metadata -------------------------------------------------
cat("Fetching GSE135251 series matrix ...\n")
gse <- getGEO("GSE135251", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gse) == 0) stop("Could not retrieve GSE135251 from GEO.")
eset  <- gse[[1]]
pheno <- pData(eset)

cat(sprintf("Samples in series matrix: %d\n", nrow(pheno)))

pheno_out <- data.frame(
  gsm            = pheno$geo_accession,
  disease        = pheno[["disease:ch1"]],
  fibrosis_stage = as.numeric(as.character(pheno[["fibrosis stage:ch1"]])),
  nas_score      = pheno[["nas score:ch1"]],
  supplementary  = pheno[["supplementary_file_1"]],
  stringsAsFactors = FALSE
)
write.csv(pheno_out, file.path(DATA_DIR, "GSE135251_pheno.csv"), row.names = FALSE)

cat("Disease distribution:\n"); print(table(pheno_out$disease, useNA = "ifany"))
cat("Fibrosis stage distribution:\n"); print(table(pheno_out$fibrosis_stage, useNA = "ifany"))

# ---- Per-sample count files ------------------------------------------------
urls <- unique(pheno_out$supplementary)
urls <- urls[!is.na(urls) & nzchar(urls)]
cat(sprintf("\nCount files to fetch: %d\n", length(urls)))

# GEO returns ftp:// links; use the https mirror which is more reliable.
urls <- sub("^ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov", urls)

n_skip <- 0; n_ok <- 0; failed <- character(0)
for (i in seq_along(urls)) {
  dest <- file.path(count_dir, basename(urls[i]))
  if (file.exists(dest) && file.size(dest) > 1000) { n_skip <- n_skip + 1; next }
  ok <- tryCatch({
    download.file(urls[i], dest, mode = "wb", method = "libcurl", quiet = TRUE)
    file.exists(dest) && file.size(dest) > 1000
  }, error = function(e) FALSE)
  if (ok) n_ok <- n_ok + 1 else failed <- c(failed, basename(urls[i]))
  if (i %% 25 == 0) cat(sprintf("  %d / %d\n", i, length(urls)))
}

cat(sprintf("\nDownloaded: %d, already present: %d, failed: %d\n",
            n_ok, n_skip, length(failed)))
if (length(failed)) {
  cat("Retry these manually:\n  ", paste(failed, collapse = "\n  "), "\n")
}

# ---- Integrity check -------------------------------------------------------
files <- list.files(count_dir, pattern = "\\.counts\\.txt\\.gz$", full.names = TRUE)
sz <- file.info(files)$size
cat(sprintf("\nFiles present: %d\n", length(files)))
cat(sprintf("Size range: %.1f - %.1f kB (median %.1f kB)\n",
            min(sz) / 1024, max(sz) / 1024, median(sz) / 1024))
small <- files[sz < 10000]
if (length(small)) cat("Suspiciously small files:\n  ", paste(basename(small), collapse = "\n  "), "\n")

if (length(files) == 216) {
  cat("\nAll 216 count files are present. Next: 01_discovery/02_count_matrix_and_DEG.R\n")
} else {
  cat("\nFile count is not 216 -- re-run this script to retry missing files.\n")
}
