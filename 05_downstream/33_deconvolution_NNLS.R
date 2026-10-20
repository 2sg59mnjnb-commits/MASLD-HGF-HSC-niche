# ============================================================================
# 33_deconvolution_NNLS.R
#
# Estimates the hepatic stellate cell fraction in each bulk sample by non-negative
# least squares deconvolution against cell-type mean profiles from the scRNA-seq
# reference, then tests whether that fraction rises with fibrosis stage.
#
# This is an orthogonal check on the bulk data: it asks whether the HSC
# compartment expands with stage, independently of any signature.
#
# Notes
#   - cell-type profiles are computed on the sparse matrix in chunks, so the
#     reference never has to be densified in full
#   - the bulk matrix is log2, the reference is log-normalised; the comparison is
#     therefore on relative profile shape, and proportions are interpreted
#     ordinally, not as absolute fractions
#
# Output
#   results/05_downstream/deconvolution_HSC_proportion.csv
#   data/derived/deconvolution.rds
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(nnls); library(clinfun)
})

# ---- 1. Reference profiles ----------------------------------------------
cat("Building reference profiles ...\n")
sc <- readRDS(file.path(DATA_DIR, "GSE202379_SeuratObject_AllCells.rds.gz"))
DefaultAssay(sc) <- "RNA"

targets <- c("Stellate", "Hepatocyte", "Cholangiocyte", "Endothelial", "Macrophage")
avail <- unique(as.character(sc$cell.annotation))
matched <- unique(unlist(lapply(targets, function(t)
  grep(t, avail, ignore.case = TRUE, value = TRUE)[1])))
matched <- matched[!is.na(matched)]
cat("Cell types used: ", paste(matched, collapse = ", "), "\n")

dat <- GetAssayData(sc, assay = "RNA", layer = "data")
ct  <- as.character(sc$cell.annotation)
rm(sc); gc(verbose = FALSE)

profiles <- list()
for (t in matched) {
  idx <- which(ct == t)
  if (length(idx) > 5000) { set.seed(SEED); idx <- sample(idx, 5000) }
  sub <- dat[, idx, drop = FALSE]

  # Chunked row means keep peak memory bounded
  chunk <- 500
  acc <- rep(0, nrow(sub))
  for (k in seq(1, length(idx), by = chunk)) {
    end <- min(k + chunk - 1, length(idx))
    acc <- acc + Matrix::rowSums(sub[, k:end, drop = FALSE])
  }
  profiles[[t]] <- acc / length(idx)
  cat(sprintf("  %-16s %d cells\n", t, length(idx)))
}
ref <- do.call(cbind, profiles); colnames(ref) <- names(profiles)
rm(dat); gc(verbose = FALSE)
cat(sprintf("Reference matrix: %d genes x %d cell types\n", nrow(ref), ncol(ref)))

# ---- 2. Bulk matrix -----------------------------------------------------
d <- load_bulk()
bulk <- to_symbol(d$log_expr)
stage <- d$fibro_stage
cat(sprintf("Bulk: %d genes x %d samples\n", nrow(bulk), ncol(bulk)))

common <- intersect(rownames(ref), rownames(bulk))
cat(sprintf("Shared genes: %d\n", length(common)))
ref  <- ref[common, , drop = FALSE]
bulk <- bulk[common, , drop = FALSE]

# ---- 3. NNLS per sample -------------------------------------------------
cat("Running NNLS deconvolution ...\n")
prop <- t(vapply(seq_len(ncol(bulk)), function(i) {
  fit <- nnls::nnls(ref, as.numeric(bulk[, i]))
  p <- coef(fit)
  if (sum(p) <= 0) p <- rep(1 / length(p), length(p))
  p / sum(p)
}, numeric(ncol(ref))))
colnames(prop) <- colnames(ref)
rownames(prop) <- colnames(bulk)

hsc_col <- grep("stellate", colnames(prop), ignore.case = TRUE, value = TRUE)[1]
hsc <- prop[, hsc_col]
cat(sprintf("HSC proportion: median %.4f, range %.4f-%.4f\n",
            median(hsc), min(hsc), max(hsc)))

# ---- 4. Trend ------------------------------------------------------------
jt <- jonckheere.test(hsc, as.numeric(factor(stage)), alternative = "increasing")
sp <- cor.test(hsc, stage, method = "spearman")
cat(sprintf("\nHSC proportion across stages:\n"))
cat(sprintf("  Jonckheere-Terpstra P = %.3e\n", jt$p.value))
cat(sprintf("  Spearman rho = %.3f (P = %.3e)\n", sp$estimate, sp$p.value))
cat("  Mean by stage:\n"); print(round(tapply(hsc, stage, mean), 4))

out <- data.frame(
  sample = rownames(prop), stage = stage, hsc_proportion = round(hsc, 4),
  other_proportion = round(1 - hsc, 4)
)
write.csv(out, res_path("05_downstream", "deconvolution_HSC_proportion.csv"),
          row.names = FALSE)

saveRDS(list(proportions = prop, reference = ref, stage = stage,
             JT_P = jt$p.value, rho = as.numeric(sp$estimate),
             celltypes = colnames(ref)),
        file.path(DATA_DIR, "derived", "deconvolution.rds"))

cat("\nNext: 05_downstream/34_multi_cohort_trends.R\n")
