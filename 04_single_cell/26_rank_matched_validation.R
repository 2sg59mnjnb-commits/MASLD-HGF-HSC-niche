# ============================================================================
# 26_rank_matched_validation.R
#
# Rank-based depth matching for the 16-gene signature.
#
# Each HGF+ HSC is paired with the HGF- HSC that occupies the same rank in the
# nCount distribution. This equalises ranking but not absolute depth, because
# the two populations do not span the same depth range (HGF+ median nCount ~288
# versus ~114 for HGF-). Rank pairing is adequate for re-scoring a fixed
# signature, but not for distance analyses, which is why the distance result is
# depth-adjusted by covariate/trimming instead (see 20_depth_correction.R).
#
# The 16-gene score of the matched pairs is then compared.
#
# Output
#   results/04_single_cell/rank_matched_validation.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

df <- load_derived("cosmx_flat.rds")
hsc <- df[df$celltype == "Stellate.cells", ]
hsc$HGF_pos <- as.numeric(hsc$HGF > 0)
cat(sprintf("HGF+ %d, HGF- %d\n", sum(hsc$HGF_pos == 1), sum(hsc$HGF_pos == 0)))

# ---- 1. Score the signature on HSCs --------------------------------------
# Only genes present on the panel can contribute; report how many that is.
get_score <- function(genes) {
  out <- rep(NA_real_, nrow(hsc))
  for (s in c(1, 2)) {
    obj <- load_cosmx(s)
    dat <- obj[["RNA"]]$data
    g <- intersect(genes, rownames(dat))
    idx <- which(hsc$slide == paste0("S", s))
    if (length(g) == 0) { rm(obj); next }
    m <- as.matrix(dat[g, hsc$cell[idx], drop = FALSE])
    z <- t(scale(t(m))); z[is.na(z)] <- 0
    out[idx] <- colMeans(z)
    rm(obj, m, z); gc(verbose = FALSE)
  }
  out
}

hsc$sig16 <- get_score(SIG16)
sig_genes_present <- SIG16  # resolved below for reporting
cat(sprintf("16-gene score computed for %d HSCs\n", sum(!is.na(hsc$sig16))))

# ---- 2. Rank matching ----------------------------------------------------
set.seed(SEED)
pos <- hsc[hsc$HGF_pos == 1, ]
neg <- hsc[hsc$HGF_pos == 0, ]
neg <- neg[order(neg$nCount), ]
cat(sprintf("\nDepth: HGF+ median nCount %.0f, HGF- median %.0f\n",
            median(pos$nCount), median(neg$nCount)))

pos$rank_q <- rank(pos$nCount, ties.method = "average") / nrow(pos)
idx <- pmax(1, pmin(nrow(neg), ceiling(pos$rank_q * nrow(neg))))
matched <- neg[idx, ]

cat(sprintf("After rank matching: HGF- median nCount %.0f (target %.0f)\n",
            median(matched$nCount), median(pos$nCount)))
p_depth <- wilcox.test(pos$nCount, matched$nCount)$p.value
cat(sprintf("Residual depth difference between pairs: Wilcoxon P = %.2e\n", p_depth))
cat("  -> rank matching does not equalise absolute depth; this is reported explicitly.\n")

# ---- 3. Signature comparison --------------------------------------------
pair <- data.frame(pos = pos$sig16, neg = matched$sig16)
pair <- pair[complete.cases(pair), ]
w <- wilcox.test(pair$pos, pair$neg, paired = TRUE)
cat(sprintf("\n16-gene score: HGF+ median %.3f, matched HGF- median %.3f\n",
            median(pair$pos), median(pair$neg)))
cat(sprintf("Paired Wilcoxon P = %.2e (n = %d pairs)\n", w$p.value, nrow(pair)))
cat(sprintf("Mean difference = %+.4f\n", mean(pair$pos - pair$neg)))

res <- data.frame(
  comparison = c("HGF+ vs rank-matched HGF-"),
  n_pairs = nrow(pair),
  median_hgfpos = round(median(pair$pos), 4),
  median_hgfneg = round(median(pair$neg), 4),
  mean_difference = round(mean(pair$pos - pair$neg), 4),
  P_paired_wilcox = w$p.value,
  residual_depth_P = p_depth
)
print(res)
write.csv(res, res_path("04_single_cell", "rank_matched_validation.csv"),
          row.names = FALSE)

# ---- 4. Gene-level direction -------------------------------------------
# Per-gene positive rates and their adjusted P values are produced in
# 04_single_cell/24_cosmx_HGFpos_vs_HGFneg_DEG.R and written to
# results/04_single_cell/signature16_positive_rate_S2.csv. They are not
# recomputed here because this script works from the flattened object, which
# carries only a small subset of panel genes.

saveRDS(list(result = res, pairs = pair, median_depth = c(
  hgfpos = median(pos$nCount), matched = median(matched$nCount))),
  file.path(DATA_DIR, "derived", "rank_matched.rds"))

cat("\nNext: 04_single_cell/27_Hep4_identity_scRNA.R\n")
