# =============================================================================
# Figure 3 -- cell-type distribution of the 16-gene HGF+ HSC signature and of
#             the 11-gene tissue signature (GSE202379)
#
# Input
#   results/signature_HGFpos_up.csv
#   results/GSE202379_seurat_clean.rds            (Part 1 only)
#   results/DEG_pseudobulk_Slide1.csv / _Slide2.csv
#   results/signature_up_single_cell_validation.csv
#   results/sig16_GO_BP.csv
#   results/signature_celltype_scores.rds
# Output
#   results/sig16_vs_sig11_in_scRNA.csv
#   results/signature_celltype_scores.rds
#   figures/final/Figure3_v2.pdf / .png
#   figures/final/FigureS_tissue_sig_celltypes.pdf / .png
#
# Part 1 scores both signatures in GSE202379 (two scoring methods).
# Part 2 builds Figure 3 (panel d = 16-gene HGF+ HSC signature).
# Part 3 builds the supplementary figure for the 11-gene tissue signature.
#
# Usage: cd ~/Desktop/sjfx3 && Rscript Figure3_redraw.R
# =============================================================================

RUN_PART1 <- TRUE    # needs the ~1.3 GB Seurat object; may be set FALSE after one run
RUN_PART2 <- TRUE
RUN_PART3 <- TRUE

setwd("~/Desktop/sjfx3")

# ggrepel places labels with a stochastic repulsion algorithm. Fix the seed so
# that repeated runs produce an identical figure.
set.seed(20240501)

sig11 <- c("EYA2", "FERMT1", "PYROXD2", "EPB41L4A", "CHI3L1",
           "CNKSR1", "ILDR1", "GSN", "TMPRSS3", "ESRRG", "SULF2")
sig16 <- read.csv("results/signature_HGFpos_up.csv")$gene

RED <- "#C73E1D"; BLUE <- "#2E86AB"; LIGHT <- "#F4A582"

# ---------------------------------------------------------------------------
# Saving. These are the exact calls that produced the published figure. NOTE: if
# your R installation has no working cairo device, the stock PDF device may
# silently drop characters it cannot map to its encoding. The PNG is
# authoritative; a faithful vector PDF requires cairo.
# ---------------------------------------------------------------------------
save_fig <- function(plot, stem, w, h) {
  ggsave(paste0(stem, ".pdf"), plot, width = w, height = h)
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 300)
}

# =============================================================================
# Part 1 -- two signatures x two scoring methods
# =============================================================================
if (RUN_PART1) {
  suppressPackageStartupMessages({ library(Seurat); library(dplyr) })

  cat("== Loading scRNA reference ==\n")
  sc <- readRDS("results/GSE202379_seurat_clean.rds")
  DefaultAssay(sc) <- "RNA"
  cat(sprintf("cells = %d | genes = %d\n", ncol(sc), nrow(sc)))
  ct <- sc$cell.annotation

  for (nm in c("sig11", "sig16")) {
    g <- get(nm); p <- g[g %in% rownames(sc)]
    cat(sprintf("%s: %d/%d genes detected", nm, length(p), length(g)))
    if (length(setdiff(g, rownames(sc))))
      cat("  not detected:", paste(setdiff(g, rownames(sc)), collapse = ", "))
    cat("\n")
  }

  cohens_d <- function(x, y) (mean(x) - mean(y)) / sqrt((var(x) + var(y)) / 2)

  ## ---- Method A: AddModuleScore, nbin = 24 (as in the original figure) ----
  sc <- AddModuleScore(sc, features = list(sig11), name = "A11", nbin = 24)
  sc <- AddModuleScore(sc, features = list(sig16), name = "A16", nbin = 24)
  sc$sig11_add <- sc$A111
  sc$sig16_add <- sc$A161

  ## ---- Method B: per-gene z-score, then mean (no control-gene binning) ----
  expr <- GetAssayData(sc, assay = "RNA", layer = "data")
  zmean <- function(genes) {
    g <- genes[genes %in% rownames(sc)]
    m <- as.matrix(expr[g, , drop = FALSE])
    z <- t(scale(t(m))); z[is.na(z)] <- 0
    colMeans(z)
  }
  sc$sig11_z <- zmean(sig11)
  sc$sig16_z <- zmean(sig16)

  ## ---- Summary ----
  tb <- sc[[]] %>%
    group_by(cell.annotation) %>%
    summarise(N = n(),
              sig11_AddModuleScore = mean(sig11_add),
              sig16_AddModuleScore = mean(sig16_add),
              sig11_zscore         = mean(sig11_z),
              sig16_zscore         = mean(sig16_z),
              .groups = "drop") %>%
    arrange(desc(sig16_AddModuleScore))
  cat("\n== Mean score per cell type, both signatures ==\n")
  print(as.data.frame(tb))
  write.csv(tb, "results/sig16_vs_sig11_in_scRNA.csv", row.names = FALSE)

  cat("\n== Stellate vs Endothelial ==\n")
  md <- sc[[]]                     # use meta.data directly to avoid Seurat [[ ]] ambiguity
  for (v in c("sig11_add", "sig11_z", "sig16_add", "sig16_z")) {
    st <- md[[v]][ct == "Stellate"]; en <- md[[v]][ct == "Endothelial"]
    cat(sprintf("%-11s  Stellate %.4f | Endothelial %.4f | Cohen's d = %+.3f | Wilcoxon P = %.2e\n",
                v, mean(st), mean(en), cohens_d(st, en), wilcox.test(st, en)$p.value))
  }

  saveRDS(sc[[]][, c("cell.annotation", "sig11_add", "sig16_add", "sig11_z", "sig16_z")],
          "results/signature_celltype_scores.rds")
  rm(sc, expr); gc()
  cat("\n== Part 1 complete ==\n")
}

# =============================================================================
# Part 2 -- main Figure 3 (panel d = 16-gene HGF+ HSC signature)
# =============================================================================
if (RUN_PART2) {
  suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("ggrepel is required")
  theme_pub <- function(bs = 10) theme_classic(base_size = bs)

  ## panel a: cross-slide pseudobulk DEG (unchanged)
  merged_pb <- inner_join(
    read.csv("results/DEG_pseudobulk_Slide1.csv") %>% dplyr::select(gene, log2FC_1 = logFC),
    read.csv("results/DEG_pseudobulk_Slide2.csv") %>% dplyr::select(gene, log2FC_2 = logFC),
    by = "gene")
  merged_pb$is_sig <- merged_pb$gene %in% sig16
  p3a <- ggplot(merged_pb, aes(log2FC_1, log2FC_2)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(data = subset(merged_pb, !is_sig), alpha = .2, size = .8, color = "grey60") +
    geom_point(data = subset(merged_pb, is_sig), alpha = .9, size = 2, color = RED) +
    ggrepel::geom_text_repel(data = subset(merged_pb, is_sig), aes(label = gene),
                             size = 2.5, max.overlaps = 20) +
    labs(title = "Cross-slide pseudobulk DEG",
         subtitle = "16-gene HGF+ HSC signature (red)",
         x = "log2FC Slide1", y = "log2FC Slide2") +
    coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)) + theme_pub()

  ## panel b: single-cell delta positive rate (unchanged)
  rv <- read.csv("results/signature_up_single_cell_validation.csv")
  rv$gene <- factor(rv$gene, levels = rv$gene[order(rv$delta_pos)])
  p3b <- ggplot(rv, aes(gene, delta_pos, fill = delta_pos)) +
    geom_col(alpha = .9) +
    geom_hline(yintercept = .05, linetype = "dashed", color = "grey40") +
    scale_fill_gradient(low = LIGHT, high = RED) + coord_flip() +
    labs(title = "Single-cell validation (Slide2)",
         subtitle = "\u0394 positive rate (HGF+ minus HGF-)", x = NULL, y = "\u0394 positive rate") +
    theme_pub() + theme(legend.position = "none")

  ## panel c: GO BP (unchanged)
  go <- read.csv("results/sig16_GO_BP.csv") %>% head(10)
  p3c <- ggplot(go, aes(reorder(Description, p.adjust), -log10(p.adjust),
                        fill = -log10(p.adjust))) +
    geom_col(alpha = .9) + scale_fill_gradient(low = LIGHT, high = RED) + coord_flip() +
    labs(title = "GO BP enrichment (16-gene signature)", x = NULL, y = "-log10(adj.P)") +
    theme_pub() + theme(legend.position = "none", axis.text.y = element_text(size = 7))

  ## panel d: corrected to the 16-gene HGF+ HSC signature
  sc16 <- readRDS("results/signature_celltype_scores.rds")
  d16 <- sc16 %>% group_by(cell.annotation) %>%
    summarise(N = n(), Mean = mean(sig16_add), .groups = "drop")
  d16$cell.annotation <- reorder(d16$cell.annotation, d16$Mean)
  p3d <- ggplot(d16, aes(cell.annotation, Mean, fill = Mean > 0)) +
    geom_col(alpha = .9) + scale_fill_manual(values = c("TRUE" = RED, "FALSE" = BLUE)) +
    coord_flip() +
    labs(title = "scRNA validation (GSE202379)",
         subtitle = "16-gene HGF+ HSC signature score", x = NULL, y = "Mean module score") +
    theme_pub() + theme(legend.position = "none")

  fig3 <- (p3a | p3b) / (p3c | p3d) + plot_annotation(tag_levels = "a")
  save_fig(fig3, "figures/final/Figure3_v2", 12, 10)
  cat("\n== Figure3_v2 saved (panel d now shows the 16-gene signature) ==\n")
}

# =============================================================================
# Part 3 -- supplementary figure: cell-type distribution of the 11-gene
#           tissue signature (reported honestly)
# =============================================================================
if (RUN_PART3) {
  suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
  sc11 <- readRDS("results/signature_celltype_scores.rds")
  d11 <- sc11 %>% group_by(cell.annotation) %>%
    summarise(N = n(), Mean = mean(sig11_add), .groups = "drop")
  d11$cell.annotation <- reorder(d11$cell.annotation, d11$Mean)

  p <- ggplot(d11, aes(cell.annotation, Mean, fill = Mean > 0)) +
    geom_col(alpha = .9) + scale_fill_manual(values = c("TRUE" = RED, "FALSE" = BLUE)) +
    coord_flip() +
    labs(title = "Cell-type distribution of the 11-gene tissue signature (GSE202379)",
         subtitle = "Tissue-level readout; not restricted to stellate cells",
         x = NULL, y = "Mean module score") +
    theme_classic(base_size = 10) + theme(legend.position = "none")

  save_fig(p, "figures/final/FigureS_tissue_sig_celltypes", 6, 4.5)
  cat("== Supplementary figure FigureS_tissue_sig_celltypes saved ==\n")
}
