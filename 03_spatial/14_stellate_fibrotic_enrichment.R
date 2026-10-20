# ============================================================================
# 14_stellate_fibrotic_enrichment.R
#
# Tests whether stellate cells are enriched in fibrotic regions.
#
# Circularity is avoided by defining the fibrotic region with a signature that
# does not contain HSC products: an inflammation/injury set (IL6, IL1B, CCL2,
# CXCL9, CXCL10, TNF, MMP2, MMP14). The 75th percentile threshold is computed
# within the non-HSC compartment only.
#
# A collagen-based definition (COL1A1/COL1A2/COL3A1/COL4A1/FN1) is reported as
# a sensitivity analysis: it is expected to give a much larger odds ratio
# because collagen is deposited by the very cells being tested.
#
# Output
#   results/03_spatial/stellate_enrichment.csv
#   results/03_spatial/stellate_enrichment_sensitivity.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

df <- load_derived("cosmx_flat.rds")

# Expression values are needed for the signature; reload the data layer of the
# slide objects and align to the flattened cells.
get_sig_score <- function(genes) {
  scores <- rep(NA_real_, nrow(df))
  for (s in c(1, 2)) {
    obj <- load_cosmx(s)
    dat <- obj[["RNA"]]$data
    g <- intersect(genes, rownames(dat))
    idx <- which(df$slide == paste0("S", s))
    if (length(g) == 0) { rm(obj); next }
    m <- as.matrix(dat[g, df$cell[idx], drop = FALSE])
    scores[idx] <- colMeans(m)
    rm(obj, m); gc(verbose = FALSE)
  }
  scores
}

enrich <- function(score, label, threshold_q = 0.75) {
  is_hsc <- df$celltype == "Stellate.cells"
  thr <- quantile(score[!is_hsc], threshold_q, na.rm = TRUE)
  region <- score >= thr

  tab <- table(factor(ifelse(is_hsc, "Stellate", "Other"),
                      levels = c("Stellate", "Other")),
               factor(ifelse(region, "Fibrotic", "NonFibrotic"),
                      levels = c("Fibrotic", "NonFibrotic")))
  print(tab)

  ft <- fisher.test(tab)
  a <- tab["Stellate", "Fibrotic"]; b <- tab["Stellate", "NonFibrotic"]
  c_ <- tab["Other", "Fibrotic"];    d_ <- tab["Other", "NonFibrotic"]
  or <- (a / b) / (c_ / d_)
  se <- sqrt(1/a + 1/b + 1/c_ + 1/d_)

  out <- data.frame(
    definition   = label,
    threshold    = round(thr, 4),
    stellate_fibrotic     = a,
    stellate_nonfibrotic  = b,
    other_fibrotic        = c_,
    other_nonfibrotic     = d_,
    OR           = round(or, 3),
    CI_low       = round(exp(log(or) - 1.96 * se), 3),
    CI_high      = round(exp(log(or) + 1.96 * se), 3),
    fisher_P     = ft$p.value,
    pct_stellate_fibrotic    = round(100 * a / (a + c_), 2),
    pct_stellate_nonfibrotic = round(100 * b / (b + d_), 2)
  )
  print(out)
  out
}

# ---- Primary definition ---------------------------------------------------
cat("=== Primary: HSC-independent inflammation/injury signature ===\n")
s_indep <- get_sig_score(FIBRO_SIG)
res_primary <- enrich(s_indep, "HSC-independent (8 genes)")
compare("stellate OR (primary)", res_primary$OR, REPORTED$OR_stellate, digits = 2)

# ---- Sensitivity: collagen-based definition -------------------------------
cat("\n=== Sensitivity: collagen-based signature ===\n")
s_coll <- get_sig_score(FIBRO_SIG_COLLAGEN)
res_coll <- enrich(s_coll, "Collagen-based (5 genes)")

out <- rbind(res_primary, res_coll)
write.csv(out, res_path("03_spatial", "stellate_enrichment.csv"), row.names = FALSE)

cat("\nInterpretation: the collagen-based OR is much larger because collagen is\n")
cat("deposited by activated HSCs, which is why it is not used as the primary\n")
cat("definition. The primary estimate is the conservative one.\n")

# ---- Save the scores for downstream/figure use ----------------------------
saveRDS(list(score_indep = s_indep, score_collagen = s_coll,
             results = out, fibrotic_indep = s_indep >= res_primary$threshold),
        file.path(DATA_DIR, "derived", "stellate_enrichment.rds"))

cat("\nNext: 03_spatial/15_HGFpos_niche_scan.R\n")
