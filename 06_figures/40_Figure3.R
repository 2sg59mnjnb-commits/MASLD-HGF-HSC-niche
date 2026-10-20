# ============================================================================
# 40_Figure3.R
#
# Figure 3. Development and validation of the HGF+ HSC signature.
#   (a) cross-slide pseudobulk log2FC with the 16-gene signature highlighted
#   (b) per-gene positive-rate difference in Slide2
#   (c) GO biological process enrichment of the 16-gene signature
#   (d) signature module score across cell types in GSE202379
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(ggrepel)})

# ---- (a) Pseudobulk scatter ---------------------------------------------
pb <- load_derived("pseudobulk.rds")
m <- pb$merged
m$in_signature <- m$gene %in% SIG16
m <- m[!is.na(m$logFC_S1) & !is.na(m$logFC_S2), ]

p_a <- ggplot(m, aes(logFC_S1, logFC_S2)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(data = m[!m$in_signature, ], alpha = 0.2, size = 0.8,
             colour = "grey60") +
  geom_point(data = m[m$in_signature, ], alpha = 0.9, size = 2,
             colour = COL_HGFPOS) +
  ggrepel::geom_text_repel(data = m[m$in_signature, ], aes(label = gene),
                           size = 2.2, max.overlaps = 20) +
  coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)) +
  labs(title = "Paired FOV pseudobulk",
       subtitle = sprintf("16-gene signature highlighted (%d genes)", sum(m$in_signature)),
       x = "log2FC, Slide1", y = "log2FC, Slide2") +
  theme_pub(9)

# ---- (b) Positive-rate difference ---------------------------------------
pr <- read.csv(res_path("04_single_cell", "signature16_positive_rate_S2.csv"))
pr <- pr[order(pr$delta_positive_rate), ]
pr$gene <- factor(pr$gene, levels = pr$gene)

p_b <- ggplot(pr, aes(delta_positive_rate, gene, fill = delta_positive_rate)) +
  geom_col(alpha = 0.9) +
  geom_vline(xintercept = 0, colour = "grey40") +
  scale_fill_gradient(low = "#F4A582", high = COL_HGFPOS, guide = "none") +
  labs(title = "Single-cell validation (Slide2)",
       subtitle = "Delta positive rate, HGF+ minus HGF-",
       x = "Delta positive rate", y = NULL) +
  theme_pub(9) + theme(axis.text.y = element_text(size = 7))

# ---- (c) GO enrichment ---------------------------------------------------
go_file <- res_path("04_single_cell", "signature16_GO_BP.csv")
if (file.exists(go_file)) {
  go <- read.csv(go_file)
  go <- head(go[order(go$p.adjust), ], 10)
  go$Description <- factor(go$Description, levels = rev(go$Description))
  p_c <- ggplot(go, aes(-log10(p.adjust), Description, fill = -log10(p.adjust))) +
    geom_col(alpha = 0.9) +
    scale_fill_gradient(low = "#F4A582", high = COL_HGFPOS, guide = "none") +
    labs(title = "GO biological process enrichment",
         subtitle = "16-gene HGF+ HSC signature",
         x = "-log10(adjusted P)", y = NULL) +
    theme_pub(9) + theme(axis.text.y = element_text(size = 7))
} else {
  p_c <- ggplot() + theme_void() +
    annotate("text", 0, 0, label = "GO output not available")
}

# ---- (d) Module score by cell type --------------------------------------
spec <- load_derived("signature_specificity.rds")
bt <- spec$by_type
bt <- bt[order(bt$sig16), ]
bt$cell.annotation <- factor(bt$cell.annotation, levels = bt$cell.annotation)

p_d <- ggplot(bt, aes(sig16, cell.annotation, fill = sig16 > 0)) +
  geom_col(alpha = 0.9) +
  scale_fill_manual(values = c(`TRUE` = COL_HGFPOS, `FALSE` = COL_HGFNEG),
                    guide = "none") +
  labs(title = "16-gene score across cell types",
       subtitle = "Independent scRNA-seq reference (GSE202379)",
       x = "Mean module score", y = NULL) +
  theme_pub(9) + theme(axis.text.y = element_text(size = 7))

fig3 <- (p_a | p_b) / (p_c | p_d) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 11))

safe_save(fig3, "Figure3", width = 12, height = 10)
cat("\nFigure 3 written. Next: 06_figures/41_Figure4_Figure5.R\n")
