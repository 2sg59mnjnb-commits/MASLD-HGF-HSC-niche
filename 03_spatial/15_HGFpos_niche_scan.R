# ============================================================================
# 15_HGFpos_niche_scan.R
#
# For every annotated cell type, counts how many neighbours each HSC has within
# 50 um, separately for HGF+ and HGF- HSCs, and normalises by the abundance of
# that cell type so that abundant populations do not dominate the comparison.
#
# Normalisation: neighbours per 1,000 target cells.
# Test: Wilcoxon rank-sum on the per-HSC normalised counts.
# Cross-slide agreement: Pearson correlation of log2 enrichment ratios across
# the cell types that carry a lineage label (the unnamed a-j clusters excluded).
#
# Output
#   results/03_spatial/niche_scan_by_slide.csv
#   results/03_spatial/niche_scan_cross_slide.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

df <- load_derived("cosmx_flat.rds")
hsc  <- df[df$celltype == "Stellate.cells", ]
hsc$HGF_pos <- as.numeric(hsc$HGF > 0)
cat(sprintf("HSC: %d (HGF+ %d, HGF- %d)\n",
            nrow(hsc), sum(hsc$HGF_pos == 1), sum(hsc$HGF_pos == 0)))

cell_types <- sort(setdiff(unique(df$celltype), "Stellate.cells"))

scan_slide <- function(slide_label) {
  d   <- df[df$slide == slide_label, ]
  h   <- hsc[hsc$slide == slide_label, ]
  pos <- h[h$HGF_pos == 1, ]; neg <- h[h$HGF_pos == 0, ]

  rows <- lapply(cell_types, function(ct) {
    tgt <- d[d$celltype == ct, ]
    if (nrow(tgt) < 50) return(NULL)

    c_pos <- count_within(pos, tgt)
    c_neg <- count_within(neg, tgt)
    n_t   <- nrow(tgt)

    # neighbours per 1,000 target cells
    v_pos <- c_pos / n_t * 1000
    v_neg <- c_neg / n_t * 1000
    p <- tryCatch(wilcox.test(v_pos, v_neg)$p.value, error = function(e) NA_real_)

    data.frame(
      slide = slide_label, celltype = ct, n_target = n_t,
      hgfpos_per1000 = mean(v_pos), hgfneg_per1000 = mean(v_neg),
      ratio = mean(v_pos) / (mean(v_neg) + 1e-9),
      log2_ratio = log2((mean(v_pos) + 1e-9) / (mean(v_neg) + 1e-9)),
      P = p, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$FDR <- p.adjust(out$P, method = "BH")
  out[order(out$FDR), ]
}

res1 <- scan_slide("S1")
res2 <- scan_slide("S2")
all_res <- rbind(res1, res2)
all_res <- all_res[order(all_res$slide, all_res$FDR), ]
write.csv(all_res, res_path("03_spatial", "niche_scan_by_slide.csv"), row.names = FALSE)

cat("\n=== Slide1: cell types most enriched around HGF+ HSCs ===\n")
print(head(res1[, c("celltype", "n_target", "ratio", "P", "FDR")], 10))
cat("\n=== Slide1: cell types least enriched ===\n")
print(head(res1[order(res1$ratio), c("celltype", "n_target", "ratio", "P", "FDR")], 6))

hep4_s1 <- res1[res1$celltype == "Hep.4", ]
hep4_s2 <- res2[res2$celltype == "Hep.4", ]
cat(sprintf("\nHep.4 enrichment around HGF+ HSCs:\n"))
cat(sprintf("  Slide1: ratio = %.3f, P = %.2e\n", hep4_s1$ratio, hep4_s1$P))
cat(sprintf("  Slide2: ratio = %.3f, P = %.2e\n", hep4_s2$ratio, hep4_s2$P))

# ---- Cross-slide agreement -----------------------------------------------
merged <- merge(
  res1[res1$n_target >= 200, c("celltype", "log2_ratio")],
  res2[res2$n_target >= 200, c("celltype", "log2_ratio")],
  by = "celltype", suffixes = c("_S1", "_S2")
)
# Restrict to cell types with a lineage label
merged <- merged[!grepl("^[a-j]$", merged$celltype), ]
merged <- merged[complete.cases(merged), ]

ct <- cor.test(merged$log2_ratio_S1, merged$log2_ratio_S2, method = "pearson")
cat(sprintf("\nCross-slide Pearson r = %.3f (P = %.2e) across %d labelled cell types\n",
            ct$estimate, ct$p.value, nrow(merged)))
compare("cross-slide r", ct$estimate, REPORTED$r_cross_slide)

write.csv(merged, res_path("03_spatial", "niche_scan_cross_slide.csv"), row.names = FALSE)

saveRDS(list(slide1 = res1, slide2 = res2,
             cross_slide = merged, cross_slide_r = ct$estimate),
        file.path(DATA_DIR, "derived", "niche_scan.rds"))

cat("\nNext: 03_spatial/16_HGFpos_proximity_mixed_models.R\n")
