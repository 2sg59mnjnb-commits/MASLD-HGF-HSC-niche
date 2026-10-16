# =============================================================================
# Figure 3 修正 + 11 基因组织签名的诚实报告
#
# -----------------------------------------------------------------------------
# 问题全貌（三套互不相同的基因集被混在了同一段文字/同一张图里）：
#
#   ① results/HGF_HSC_sig_in_scRNA.csv  —— 原 Figure 3d 实际画的东西
#        来源：HGFpos_HSC_DEGs_strict.csv 中 avg_log2FC > 0 的 63 个基因
#        结果：Stellate 0.293 > Endothelial 0.084（HSC 特异，d ≈ 0.92）
#
#   ② Cohen's d = 0.92 —— 正文 Results §10 / Discussion §2 引用的数字
#        来源：signature_HGFpos_up.csv 的 16 个基因
#
#   ③ 11-gene tissue signature —— 正文与图注真正声称的签名
#        来源：EYA2/FERMT1/PYROXD2/EPB41L4A/CHI3L1/CNKSR1/ILDR1/GSN/TMPRSS3/ESRRG/SULF2
#        与 ① 和 ② 的交集均为 **0 个基因**（完全不相交）
#
#   实测（本机 GSE202379，99,687 细胞）：
#        11 基因签名在 Cholangiocytes 最高（0.194），Stellate 仅 0.053，
#        Endothelial 0.038 → Cohen's d (Stellate vs Endothelial) = 0.102，不是 0.92。
#
#   结论：§10 用一套完全不相交的基因集去"验证"11 基因签名，逻辑上不成立；
#         而真正测 11 基因签名时，它在 scRNA 层面**不是** HSC 特异的。
# -----------------------------------------------------------------------------
#
# 本脚本做的事：
#   Part 1  在同一对象上同时打分：16 基因 HGF+ HSC 签名 + 11 基因组织签名
#           并给出两种打分方法（AddModuleScore 与逐基因 z-score）的交叉验证
#   Part 2  重画 Figure 3：panel d 改为 **16 基因 HGF+ HSC 签名**（与 panel a/b/c
#           同源，图内自洽，且该结果确实成立）
#   Part 3  另出一张补充图：11 基因组织签名在各类细胞中的分布（如实报告）
#
# 用法：  cd ~/Desktop/sjfx3
#         Rscript Figure3_redraw.R
#
# 输出：  results/sig16_vs_sig11_in_scRNA.csv
#         results/signature_celltype_scores.rds      （轻量，便于日后重绘）
#         figures/final/Figure3_v2.pdf / .png        （主图，panel d 已修正）
#         figures/final/FigureS_tissue_sig_celltypes.pdf / .png （补充图，如实报告）
# =============================================================================

RUN_PART1 <- TRUE    # 需要 ~1.3 GB Seurat 对象；跑过一次后可关掉
RUN_PART2 <- TRUE
RUN_PART3 <- TRUE

setwd("~/Desktop/sjfx3")

sig11 <- c("EYA2", "FERMT1", "PYROXD2", "EPB41L4A", "CHI3L1",
           "CNKSR1", "ILDR1", "GSN", "TMPRSS3", "ESRRG", "SULF2")
sig16 <- read.csv("results/signature_HGFpos_up.csv")$gene

RED <- "#C73E1D"; BLUE <- "#2E86AB"; LIGHT <- "#F4A582"

# =============================================================================
# Part 1 —— 两套签名 × 两种打分方法
# =============================================================================
if (RUN_PART1) {
  suppressPackageStartupMessages({ library(Seurat); library(dplyr) })

  cat("== 载入 scRNA 参考集 ==\n")
  sc <- readRDS("results/GSE202379_seurat_clean.rds")
  DefaultAssay(sc) <- "RNA"
  cat(sprintf("cells = %d | genes = %d\n", ncol(sc), nrow(sc)))
  ct <- sc$cell.annotation

  for (nm in c("sig11", "sig16")) {
    g <- get(nm); p <- g[g %in% rownames(sc)]
    cat(sprintf("%s: %d/%d 基因可检出", nm, length(p), length(g)))
    if (length(setdiff(g, rownames(sc))))
      cat("  未检出:", paste(setdiff(g, rownames(sc)), collapse = ", "))
    cat("\n")
  }

  cohens_d <- function(x, y) (mean(x) - mean(y)) / sqrt((var(x) + var(y)) / 2)

  ## ---- 方法 A：AddModuleScore，nbin = 24（与原图/原稿一致）----
  sc <- AddModuleScore(sc, features = list(sig11), name = "A11", nbin = 24)
  sc <- AddModuleScore(sc, features = list(sig16), name = "A16", nbin = 24)
  sc$sig11_add <- sc$A111
  sc$sig16_add <- sc$A161

  ## ---- 方法 B：逐基因 z-score 后取均值（不依赖对照基因分箱）----
  expr <- GetAssayData(sc, assay = "RNA", layer = "data")
  zmean <- function(genes) {
    g <- genes[genes %in% rownames(sc)]
    m <- as.matrix(expr[g, , drop = FALSE])
    z <- t(scale(t(m))); z[is.na(z)] <- 0
    colMeans(z)
  }
  sc$sig11_z <- zmean(sig11)
  sc$sig16_z <- zmean(sig16)

  ## ---- 汇总 ----
  tb <- sc[[]] %>%
    group_by(cell.annotation) %>%
    summarise(N = n(),
              sig11_AddModuleScore = mean(sig11_add),
              sig16_AddModuleScore = mean(sig16_add),
              sig11_zscore         = mean(sig11_z),
              sig16_zscore         = mean(sig16_z),
              .groups = "drop") %>%
    arrange(desc(sig16_AddModuleScore))
  cat("\n== 两套签名在各类细胞中的均值 ==\n")
  print(as.data.frame(tb))
  write.csv(tb, "results/sig16_vs_sig11_in_scRNA.csv", row.names = FALSE)

  cat("\n== Stellate vs Endothelial ==\n")
  md <- sc[[]]                     # 取 meta.data，避免 Seurat [[ ]] 的歧义
  for (v in c("sig11_add", "sig11_z", "sig16_add", "sig16_z")) {
    st <- md[[v]][ct == "Stellate"]; en <- md[[v]][ct == "Endothelial"]
    cat(sprintf("%-11s  Stellate %.4f | Endothelial %.4f | Cohen's d = %+.3f | Wilcoxon P = %.2e\n",
                v, mean(st), mean(en), cohens_d(st, en), wilcox.test(st, en)$p.value))
  }

  saveRDS(sc[[]][, c("cell.annotation", "sig11_add", "sig16_add", "sig11_z", "sig16_z")],
          "results/signature_celltype_scores.rds")
  rm(sc, expr); gc()
  cat("\n== Part 1 完成 ==\n")
}

# =============================================================================
# Part 2 —— 主图 Figure 3（panel d = 16 基因 HGF+ HSC 签名）
# =============================================================================
if (RUN_PART2) {
  suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("需要 ggrepel")
  theme_pub <- function(bs = 10) theme_classic(base_size = bs)

  ## panel a：跨切片 pseudobulk DEG（不变）
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

  ## panel b：单细胞 Δ positive rate（不变）
  rv <- read.csv("results/signature_up_single_cell_validation.csv")
  rv$gene <- factor(rv$gene, levels = rv$gene[order(rv$delta_pos)])
  p3b <- ggplot(rv, aes(gene, delta_pos, fill = delta_pos)) +
    geom_col(alpha = .9) +
    geom_hline(yintercept = .05, linetype = "dashed", color = "grey40") +
    scale_fill_gradient(low = LIGHT, high = RED) + coord_flip() +
    labs(title = "Single-cell validation (Slide2)",
         subtitle = "Δ positive rate (HGF+ minus HGF-)", x = NULL, y = "Δ positive rate") +
    theme_pub() + theme(legend.position = "none")

  ## panel c：GO BP（不变）
  go <- read.csv("results/sig16_GO_BP.csv") %>% head(10)
  p3c <- ggplot(go, aes(reorder(Description, p.adjust), -log10(p.adjust),
                        fill = -log10(p.adjust))) +
    geom_col(alpha = .9) + scale_fill_gradient(low = LIGHT, high = RED) + coord_flip() +
    labs(title = "GO BP enrichment (16-gene signature)", x = NULL, y = "-log10(adj.P)") +
    theme_pub() + theme(legend.position = "none", axis.text.y = element_text(size = 7))

  ## panel d：修正为 16 基因 HGF+ HSC 签名
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
  ggsave("figures/final/Figure3_v2.pdf", fig3, width = 12, height = 10)
  ggsave("figures/final/Figure3_v2.png", fig3, width = 12, height = 10, dpi = 300)
  cat("\n== Figure3_v2 已保存（panel d 已改为 16 基因签名）==\n")
}

# =============================================================================
# Part 3 —— 补充图：11 基因组织签名的细胞类型分布（如实报告）
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

  ggsave("figures/final/FigureS_tissue_sig_celltypes.pdf", p, width = 6, height = 4.5)
  ggsave("figures/final/FigureS_tissue_sig_celltypes.png", p, width = 6, height = 4.5, dpi = 300)
  cat("== 补充图 FigureS_tissue_sig_celltypes 已保存 ==\n")
}
