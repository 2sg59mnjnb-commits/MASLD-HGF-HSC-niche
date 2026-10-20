# ============================================================================
# 02_count_matrix_and_DEG.R
#
# Builds the training expression matrix from the 216 per-sample count files and
# runs the limma differential-expression analysis that defines the fibrosis
# gene set.
#
# Steps
#   1. merge 216 single-sample count files into one gene x sample matrix
#   2. align with the phenotype table (order of rownames, not sample IDs --
#      GEO sample order and matrix column order are identical here)
#   3. filter low-expression genes (count >= 10 in >= 3 samples), log2(x+1)
#   4. limma moderated t-test, High (F2-F4) vs Low (F0-F1)
#   5. DEGs: |log2FC| > 0.5 and BH-adjusted P < 0.05
#
# Output
#   data/derived/GSE135251_aligned_data.rds   counts + phenotype + groups
#   data/derived/GSE135251_DEG_results.rds    deg, deg_sig, log_expr, groups
#   results/01_discovery/GSE135251_DEGs.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({
  library(limma)
  library(GEOquery)
})

count_dir <- file.path(DATA_DIR, "GSE135251_counts")
files <- list.files(count_dir, pattern = "\\.counts\\.txt\\.gz$", full.names = TRUE)
if (length(files) < 200) stop("Run 01_discovery/01_download_GSE135251.R first.")
cat(sprintf("Count files: %d\n", length(files)))

# ---- 1. Merge --------------------------------------------------------------
# Each file has a header row and one gene per row. Reading with row.names = 1
# avoids materialising a data.frame per file.
cat("Reading count files (5-10 min) ...\n")
mat_list <- vector("list", length(files))
gsm_ids  <- sub("_.*", "", basename(files))
ok <- logical(length(files))

for (i in seq_along(files)) {
  ok[i] <- tryCatch({
    df <- read.table(gzfile(files[i]), header = TRUE, row.names = 1,
                     stringsAsFactors = FALSE, check.names = FALSE)
    mat_list[[i]] <- setNames(as.numeric(df[[1]]), rownames(df))
    TRUE
  }, error = function(e) {
    cat(sprintf("  skipping unreadable file: %s\n", basename(files[i]))); FALSE
  })
  if (i %% 50 == 0) cat(sprintf("  %d / %d\n", i, length(files)))
}

mat_list <- mat_list[ok]; gsm_ids <- gsm_ids[ok]
cat(sprintf("Usable samples: %d\n", length(mat_list)))

all_genes <- Reduce(union, lapply(mat_list, names))
count_matrix <- matrix(0, nrow = length(all_genes), ncol = length(mat_list),
                       dimnames = list(all_genes, gsm_ids))
for (i in seq_along(mat_list)) {
  count_matrix[names(mat_list[[i]]), i] <- mat_list[[i]]
}
rm(mat_list); gc(verbose = FALSE)
cat(sprintf("Merged matrix: %d genes x %d samples\n",
            nrow(count_matrix), ncol(count_matrix)))
write.csv(data.frame(gene = rownames(count_matrix)),
          file.path(DATA_DIR, "GSE135251_genes.csv"), row.names = FALSE)

# ---- 2. Phenotype ---------------------------------------------------------
# The matrix column order is the GEO sample order, so align by position after
# checking that the counts match.
pheno <- read.csv(file.path(DATA_DIR, "GSE135251_pheno.csv"),
                  stringsAsFactors = FALSE)
stopifnot(nrow(pheno) == ncol(count_matrix))
if (!identical(pheno$gsm, colnames(count_matrix))) {
  warning("Sample order differs from GEO order; aligning by GSM ID.")
  idx <- match(pheno$gsm, colnames(count_matrix))
  stopifnot(!any(is.na(idx)))
  count_matrix <- count_matrix[, idx]
}

fibro_stage <- pheno$fibrosis_stage
group <- factor(ifelse(fibro_stage >= 2, "High", "Low"), levels = c("Low", "High"))

cat("\nGroup sizes:\n"); print(table(group, useNA = "ifany"))
cat("Disease labels:\n"); print(table(pheno$disease, useNA = "ifany"))

saveRDS(list(count_matrix = count_matrix, pheno = pheno,
             fibro_stage = fibro_stage, group = group),
        file.path(DATA_DIR, "derived", "GSE135251_aligned_data.rds"))

# ---- 3. Filter and transform ----------------------------------------------
keep <- rowSums(count_matrix >= 10) >= 3
count_matrix <- count_matrix[keep, ]
cat(sprintf("\nGenes after low-expression filter: %d\n", nrow(count_matrix)))
log_expr <- log2(count_matrix + 1)

# ---- 4. limma -------------------------------------------------------------
design <- model.matrix(~ group)
fit <- eBayes(lmFit(log_expr, design))
deg <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

# ---- 5. DEGs --------------------------------------------------------------
deg_sig <- deg[abs(deg$logFC) > 0.5 & deg$adj.P.Val < 0.05, ]
cat(sprintf("\nGenes tested: %d\n", nrow(deg)))
cat(sprintf("DEGs (|log2FC| > 0.5, adj.P < 0.05): %d\n", nrow(deg_sig)))
cat(sprintf("  up in high fibrosis:   %d\n", sum(deg_sig$logFC > 0.5)))
cat(sprintf("  down in high fibrosis: %d\n", sum(deg_sig$logFC < -0.5)))

compare("DEG count", nrow(deg_sig), REPORTED$n_deg, digits = 0)

cat("\nTop 15 DEGs:\n")
print(head(deg_sig[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")], 15))

# ---- Save -----------------------------------------------------------------
write.csv(deg, res_path("01_discovery", "GSE135251_all_genes_DE.csv"))
write.csv(deg_sig, res_path("01_discovery", "GSE135251_DEGs.csv"))
saveRDS(list(deg = deg, deg_sig = deg_sig, log_expr = log_expr,
             group = group, fibro_stage = fibro_stage),
        file.path(DATA_DIR, "derived", "GSE135251_DEG_results.rds"))

cat("\nNext: 01_discovery/03_wgcna.R\n")
