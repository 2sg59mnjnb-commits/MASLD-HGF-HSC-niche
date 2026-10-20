# ============================================================================
# 29_TF_activity.R
#
# Transcription-factor activity in HGF+ versus HGF- HSCs, inferred from
# regulon-level target expression with decoupleR's univariate linear model and
# the DoRothEA A/B/C regulons.
#
# A TF is retained when
#   - at least 5 of its targets are detectable in the HSC expression matrix
#   - its activity differs between HGF+ and HGF- HSCs on both slides, with
#     concordant direction, after Benjamini-Hochberg correction
#
# The regulon-level result is hypothesis-generating; 30_TF_HGF_coupling_and_motif.R
# tests whether the inferred regulators could act at the HGF promoter.
#
# Output
#   results/04_single_cell/TF_activity_differential.csv
#   results/04_single_cell/TF_activity_consensus.csv
#   data/derived/tf_activity.rds
# ============================================================================

source("R/config.R")
source("R/helpers.R")

if (!requireNamespace("decoupleR", quietly = TRUE) ||
    !requireNamespace("dorothea", quietly = TRUE)) {
  cat("decoupleR / dorothea not installed; skipping.\n")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(decoupleR); library(dorothea); library(Seurat); library(dplyr)
})

data("dorothea_hs", package = "dorothea")
regulon <- dorothea_hs[dorothea_hs$confidence %in% c("A", "B", "C"), ]

run_slide <- function(slide) {
  obj <- load_cosmx(slide)
  ct  <- as.character(obj@meta.data[[celltype_col(slide)]])
  dat <- obj[["RNA"]]$data

  hsc <- which(ct == "Stellate.cells")
  mat <- as.matrix(dat[, hsc])
  mat <- mat[rowSums(mat > 0) > ncol(mat) * 0.05, , drop = FALSE]
  cat(sprintf("Slide%d: %d genes x %d HSCs\n", slide, nrow(mat), ncol(mat)))

  tf <- decoupleR::run_ulm(mat = mat, net = regulon,
                           .source = "tf", .target = "target",
                           .mor = "mor", minsize = 5)

  score_col <- if ("score" %in% colnames(tf)) "score" else "statistic"
  wide <- tf %>%
    dplyr::select(condition, source, value = dplyr::all_of(score_col)) %>%
    tidyr::pivot_wider(names_from = source, values_from = value) %>%
    as.data.frame()
  rownames(wide) <- wide$condition
  wide$condition <- NULL
  wide <- wide[colnames(mat), , drop = FALSE]

  hgf <- as.numeric(dat["HGF", hsc])
  grp <- ifelse(hgf > 0, "pos", "neg")

  # Keep TFs with at least 5 detectable targets
  n_targets <- regulon %>%
    dplyr::filter(target %in% rownames(mat)) %>%
    dplyr::count(tf, name = "n_targets")
  keep_tf <- n_targets$tf[n_targets$n_targets >= 5]
  wide <- wide[, intersect(keep_tf, colnames(wide)), drop = FALSE]

  diff <- do.call(rbind, lapply(colnames(wide), function(t) {
    a <- wide[[t]][grp == "pos"]; b <- wide[[t]][grp == "neg"]
    data.frame(TF = t,
               delta = mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE),
               P = wilcox.test(a, b)$p.value)
  }))
  diff$FDR <- p.adjust(diff$P, method = "BH")
  diff$slide <- paste0("S", slide)
  diff <- diff[order(diff$FDR), ]

  rm(obj); gc(verbose = FALSE)
  list(matrix = wide, diff = diff, n_targets = n_targets)
}

cat("Running decoupleR ULM ...\n")
r1 <- run_slide(1)
r2 <- run_slide(2)

all_diff <- rbind(r1$diff, r2$diff)
write.csv(all_diff, res_path("04_single_cell", "TF_activity_differential.csv"),
          row.names = FALSE)

cat("\nSlide1 top TFs:\n")
print(head(r1$diff[, c("TF", "delta", "FDR")], 8), row.names = FALSE)

# ---- Consensus across slides --------------------------------------------
m <- merge(r1$diff[, c("TF", "delta", "FDR")],
           r2$diff[, c("TF", "delta", "FDR")],
           by = "TF", suffixes = c("_S1", "_S2"))
consensus <- m[m$FDR_S1 < 0.05 & m$FDR_S2 < 0.05 &
                 sign(m$delta_S1) == sign(m$delta_S2), ]
consensus$abs_sum <- abs(consensus$delta_S1) + abs(consensus$delta_S2)
consensus <- consensus[order(-consensus$abs_sum), ]

cat(sprintf("\nTFs significant on both slides with concordant direction: %d\n",
            nrow(consensus)))
print(head(consensus[, c("TF", "delta_S1", "delta_S2", "FDR_S1", "FDR_S2")], 15),
      row.names = FALSE)

# Reported consensus set uses |delta| > 0.3 on both slides
final <- consensus[abs(consensus$delta_S1) > 0.3 & abs(consensus$delta_S2) > 0.3, ]
cat(sprintf("\nWith |delta| > 0.3 on both slides: %d TFs\n", nrow(final)))
print(final[, c("TF", "delta_S1", "delta_S2")], row.names = FALSE)
cat(sprintf("Overlap with the reported set (%s): %d / %d\n",
            paste(TF_CONSENSUS, collapse = ", "),
            length(intersect(final$TF, TF_CONSENSUS)), length(TF_CONSENSUS)))

write.csv(final, res_path("04_single_cell", "TF_activity_consensus.csv"),
          row.names = FALSE)

# ---- Store per-cell activity for the coupling analysis -------------------
saveRDS(list(slide1 = r1$matrix, slide2 = r2$matrix,
             diff = all_diff, consensus = consensus, final = final,
             n_targets = rbind(cbind(slide = "S1", r1$n_targets),
                               cbind(slide = "S2", r2$n_targets))),
        file.path(DATA_DIR, "derived", "tf_activity.rds"))

cat("\nNext: 04_single_cell/30_TF_HGF_coupling_and_motif.R\n")
