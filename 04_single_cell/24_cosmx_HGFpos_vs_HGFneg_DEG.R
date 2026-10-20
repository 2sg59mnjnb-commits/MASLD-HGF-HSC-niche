# ============================================================================
# 24_cosmx_HGFpos_vs_HGFneg_DEG.R
#
# Differential expression between HGF+ and HGF- HSCs in each CosMx slide.
#
# The per-cell comparison is reported with a strict filter because HGF+ and
# HGF- HSCs differ in transcript count; the strict filter requires a large
# effect, a large positive-rate difference and a stringent adjusted P, so that
# small depth-driven shifts do not pass.
#
# These gene lists are used for
#   - the TFAP2A target-set overlap (the down-regulated list)
#   - the GO enrichment shown in the manuscript
#
# Output
#   results/04_single_cell/cosmx_HGFpos_vs_HGFneg_DEG.csv
#   results/04_single_cell/cosmx_DEG_summary.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(Seurat); library(dplyr)})

run_slide <- function(slide) {
  obj <- load_cosmx(slide)
  ct  <- as.character(obj@meta.data[[celltype_col(slide)]])

  hsc <- subset(obj, cells = colnames(obj)[ct == "Stellate.cells"])
  rm(obj); gc(verbose = FALSE)

  hsc$hgf_status <- ifelse(hsc[["RNA"]]$data["HGF", ] > 0, "HGF_pos", "HGF_neg")
  cat(sprintf("Slide%d: HGF+ %d, HGF- %d\n", slide,
              sum(hsc$hgf_status == "HGF_pos"),
              sum(hsc$hgf_status == "HGF_neg")))

  Idents(hsc) <- hsc$hgf_status
  deg <- FindMarkers(hsc, ident.1 = "HGF_pos", ident.2 = "HGF_neg",
                     logfc.threshold = 0.5, min.pct = 0.1,
                     test.use = "wilcox", verbose = FALSE)
  deg$gene <- rownames(deg)

  # Strict filter: large effect + large positive-rate difference + stringent FDR
  strict <- deg[abs(deg$avg_log2FC) > 1 &
                  deg$p_val_adj < 0.01 &
                  abs(deg$pct.1 - deg$pct.2) > 0.2, ]

  cat(sprintf("  DEGs (nominal): %d | strict: %d (up %d, down %d)\n",
              nrow(deg), nrow(strict),
              sum(strict$avg_log2FC > 0), sum(strict$avg_log2FC < 0)))

  strict$slide <- paste0("S", slide)
  strict[order(-strict$avg_log2FC), ]
}

res <- rbind(run_slide(1), run_slide(2))
write.csv(res, res_path("04_single_cell", "cosmx_HGFpos_vs_HGFneg_DEG.csv"),
          row.names = FALSE)

up   <- res[res$avg_log2FC > 0, ]
down <- res[res$avg_log2FC < 0, ]
cat(sprintf("\nTotal strict DEGs across both slides: %d (up %d, down %d)\n",
            nrow(res), nrow(up), nrow(down)))

summ <- data.frame(
  direction = c("up", "down"),
  n = c(nrow(up), nrow(down)),
  genes = c(paste(head(up$gene[order(-up$avg_log2FC)], 25), collapse = ", "),
            paste(head(down$gene[order(down$avg_log2FC)], 25), collapse = ", "))
)
print(summ)
write.csv(summ, res_path("04_single_cell", "cosmx_DEG_summary.csv"), row.names = FALSE)

saveRDS(list(deg = res, up = up$gene, down = down$gene),
        file.path(DATA_DIR, "derived", "cosmx_deg.rds"))

cat("\nNote: HGF+ HSCs have a markedly higher transcript count than HGF- HSCs in\n")
cat("these sections. The strict filter reduces, but does not remove, depth\n")
cat("confounding; see 20_depth_correction.R and 25 for the depth-adjusted\n")
cat("comparisons that are used for the signature.\n")
cat("\nNext: 04_single_cell/25_pseudobulk_cross_slide.R\n")

# ---- 2. Per-gene positive rates for the 16-gene signature (Slide2) -------
# This is the quantity plotted in Figure 3b: the difference in positive rate
# between HGF+ and HGF- HSCs for each signature gene, computed directly from the
# Slide2 expression matrix rather than from the filtered DEG table, so that no
# gene is excluded by a transcriptome-wide threshold.
cat("\n=== 16-gene signature: positive rates (Slide2) ===\n")
obj  <- load_cosmx(2)
ct2  <- as.character(obj@meta.data[[celltype_col(2)]])
dat2 <- obj[["RNA"]]$data

hsc2 <- which(ct2 == "Stellate.cells")
hgf2 <- as.numeric(dat2["HGF", hsc2])
grp2 <- ifelse(hgf2 > 0, "pos", "neg")
cat(sprintf("HGF+ %d, HGF- %d\n", sum(grp2 == "pos"), sum(grp2 == "neg")))

genes2 <- intersect(SIG16, rownames(dat2))
cat(sprintf("Signature genes on panel: %d / 16\n", length(genes2)))

sig_tab <- do.call(rbind, lapply(genes2, function(g) {
  x <- as.numeric(dat2[g, hsc2])
  a <- x[grp2 == "pos"]; b <- x[grp2 == "neg"]
  data.frame(
    gene = g,
    positive_rate_hgfpos = mean(a > 0),
    positive_rate_hgfneg = mean(b > 0),
    mean_hgfpos = mean(a),
    mean_hgfneg = mean(b),
    P = wilcox.test(a, b)$p.value,
    stringsAsFactors = FALSE
  )
}))
sig_tab$FDR <- p.adjust(sig_tab$P, method = "BH")
sig_tab$delta_positive_rate <- sig_tab$positive_rate_hgfpos - sig_tab$positive_rate_hgfneg
sig_tab <- sig_tab[order(-sig_tab$delta_positive_rate), ]

print(sig_tab, row.names = FALSE)
cat(sprintf("\nGenes with higher positive rate in HGF+ HSCs: %d / %d\n",
            sum(sig_tab$delta_positive_rate > 0), nrow(sig_tab)))
cat(sprintf("Significant at FDR < 0.05: %d / %d\n",
            sum(sig_tab$FDR < 0.05), nrow(sig_tab)))

write.csv(sig_tab, res_path("04_single_cell", "signature16_positive_rate_S2.csv"),
          row.names = FALSE)
rm(obj); gc(verbose = FALSE)

# ---- 3. IL6-family receptor positive rates in HSCs -----------------------
# This is the data behind Figure 4: positive rate of each receptor in HGF+
# versus HGF- HSCs, separately per slide, with a Wilcoxon test per gene.
cat("\n=== IL6-family receptor positive rates (Figure 4) ===\n")
rec_rows <- list()
for (s in c(1, 2)) {
  obj <- load_cosmx(s)
  dat <- obj[["RNA"]]$data
  hsc <- which(as.character(obj@meta.data[[celltype_col(s)]]) == "Stellate.cells")
  hgf <- as.numeric(dat["HGF", hsc])
  grp <- ifelse(hgf > 0, "pos", "neg")

  for (g in intersect(IL6_RECEPTORS, rownames(dat))) {
    x <- as.numeric(dat[g, hsc])
    a <- x[grp == "pos"]; b <- x[grp == "neg"]
    rec_rows[[length(rec_rows) + 1]] <- data.frame(
      slide = paste0("S", s), gene = g,
      n_hgfpos = length(a), n_hgfneg = length(b),
      positive_rate_hgfpos = mean(a > 0),
      positive_rate_hgfneg = mean(b > 0),
      mean_hgfpos = mean(a), mean_hgfneg = mean(b),
      P = wilcox.test(a, b)$p.value,
      stringsAsFactors = FALSE
    )
  }
  cat(sprintf("Slide%d: %d HGF+ / %d HGF- HSCs\n",
              s, sum(grp == "pos"), sum(grp == "neg")))
  rm(obj); gc(verbose = FALSE)
}
rec <- do.call(rbind, rec_rows)
rec$delta_positive_rate <- rec$positive_rate_hgfpos - rec$positive_rate_hgfneg
rec <- rec[order(rec$slide, -abs(rec$delta_positive_rate)), ]
print(rec, row.names = FALSE)
write.csv(rec, res_path("04_single_cell", "IL6_receptor_positive_rates.csv"),
          row.names = FALSE)

cat("\nReceptors higher in HGF+ HSCs: ",
    paste(unique(rec$gene[rec$delta_positive_rate > 0 & rec$P < 0.05]),
          collapse = ", "), "\n")
cat("Receptors lower in HGF+ HSCs: ",
    paste(unique(rec$gene[rec$delta_positive_rate < 0 & rec$P < 0.05]),
          collapse = ", "), "\n")
