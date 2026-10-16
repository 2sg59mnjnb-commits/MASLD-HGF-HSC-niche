# =============================================================================
# Figure 2 重建脚本 —— 修正两处「图与图注不符」
#
# ── 问题诊断（本次逐图核对得出）─────────────────────────────────────────────
# 现状：正文与图注描述的是 **9 个面板 (a)–(i)**，其中
#         (d) = HGF 阳性阈值（HGF>0/>1/>2）下 Hep.4 与 Periportal LSEC 的富集比
#         (i) = 25–150 µm 距离梯度分析（正文 Results 引用了 Figure 2i）
#
# 但现有 R 代码拼装的是 **8 个面板**：
#         fig2 <- (p2a|p2b|p2c) / (p2d|p2e|p2f) / (p2g|p2h|plot_spacer())
#   且该代码里的 p2d 用的是「距离阈值 vs 中位距离」（错误内容），
#   而 p2i 在整份代码中根本不存在（只有渲染好的 PNG）。
#
# 因此稿件里嵌入的 Figure 2 面板 (d) 内容错误、且缺少面板 (i)。
#
# ── 本脚本做什么 ─────────────────────────────────────────────────────────
#   1. 把 p2d 改回 HGF 阳性阈值版本（数据源 results/HGF_threshold_sensitivity.csv）
#   2. 补入面板 (i)：直接导入已有的 Figure2i_distance_gradient.png
#   3. 按 3×3 拼装成完整 9 面板图，输出 figures/final/Figure2_v2.pdf / .png
#
# 用法：  cd ~/Desktop/sjfx3 && Rscript Figure2_rebuild.R
# 输出：  figures/final/Figure2_v2.pdf / .png   （确认后改名覆盖 Figure2）
#
# ⚠️ 面板 (i) 的生成代码**不在你的代码仓库里**（analysis_code.R 中搜不到），
#    所以这里只能导入现成的 PNG。建议你找回那段代码并一并补进仓库，
#    否则 Figure 2i 无法被复现——审稿人可能会问。
# =============================================================================

setwd("~/Desktop/sjfx3")
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
  library(png); library(grid); library(ggrepel)
})

OUT <- "figures/final"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
theme_pub <- function(bs = 9) theme_classic(base_size = bs)

# ---------------------------------------------------------------- 2a 生态位扫描
res1 <- read.csv("results/HGF_HSC_niche_scan.csv");        res1$Slide <- "Slide1"
res2 <- read.csv("results/HGF_HSC_niche_scan_Slide2.csv");  res2$Slide <- "Slide2"
niche_df <- rbind(res1, res2) %>% filter(N_total >= 200) %>% mutate(Score = log2(Ratio + 1e-6))

p2a <- ggplot(niche_df, aes(reorder(CellType, Score), Score, fill = Slide)) +
  geom_col(position = position_dodge(0.8), alpha = 0.9) +
  scale_fill_manual(values = c("Slide1" = "#C73E1D", "Slide2" = "#2E86AB")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  coord_flip() +
  labs(title = "HGF+ HSC niche scan (50 µm)",
       subtitle = "log2(enrichment ratio) across two slides",
       x = NULL, y = "log2 ratio") +
  theme_pub() + theme(legend.position = "top", axis.text.y = element_text(size = 7))

# ---------------------------------------------------------------- 2b FOV 一致性
fov_res <- readRDS("results/Reverse_FOV_characterization.rds")
n_cons <- nrow(fov_res$consistent_fovs); n_rev <- nrow(fov_res$reverse_fovs)
fov_df <- data.frame(Direction = factor(c("HGF+ closer", "Other"), levels = c("HGF+ closer", "Other")),
                     N = c(n_cons, n_rev))
fov_df$Pct <- fov_df$N / sum(fov_df$N) * 100

p2b <- ggplot(fov_df, aes(Direction, N, fill = Direction)) +
  geom_col(alpha = 0.9, width = 0.6) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", N, Pct)), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("HGF+ closer" = "#C73E1D", "Other" = "grey70")) +
  labs(title = "FOV-level direction consistency",
       subtitle = sprintf("Binomial P = %.2e; n = %d", fov_res$binom_p, n_cons + n_rev),
       x = NULL, y = "Number of FOVs") +
  theme_pub() + theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, max(fov_df$N) * 1.15))

# ------------------------------------------------------------ 2c 跨切片一致性
merged <- read.csv("results/niche_scan_cross_slide.csv") %>%
  filter(grepl("[A-Z]", CellType),
         grepl("\\.", CellType) | grepl("[A-Z][a-z]", CellType),
         !CellType %in% c("h","j","i","e","c","a","b","d","f","g"))
ct_2c <- cor.test(merged$log2_R1, merged$log2_R2)

p2c <- ggplot(merged, aes(log2_R1, log2_R2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_abline(slope = 1, intercept = 0, color = "grey80", linetype = "dotted") +
  geom_point(size = 2.5, alpha = 0.85, color = "#C73E1D") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey80") +
  ggrepel::geom_text_repel(aes(label = CellType), size = 2.3, max.overlaps = 20) +
  labs(title = "Cross-slide niche consistency",
       subtitle = sprintf("Pearson r = %.3f, P = %.2e", ct_2c$estimate, ct_2c$p.value),
       x = "log2 ratio Slide1", y = "log2 ratio Slide2") +
  coord_fixed() + theme_pub()

# ============================ 2d 【关键修正】HGF 阳性阈值敏感性 ==============
# 旧代码这里读的是 results/revision/A4_threshold_sensitivity.csv，
# 画的是「距离阈值 vs 中位距离」——与图注和正文都不符。改回下面这版。
sens <- read.csv("results/HGF_threshold_sensitivity.csv")
sens_df <- sens %>%
  filter(CellType %in% c("Hep.4", "Periportal.LSECs")) %>%
  mutate(Threshold = factor(Threshold, levels = c("HGF>0", "HGF>1", "HGF>2")))

p2d <- ggplot(sens_df, aes(Threshold, Ratio, fill = CellType)) +
  geom_col(position = position_dodge(0.8), alpha = 0.9, width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Hep.4" = "#C73E1D", "Periportal.LSECs" = "#2E86AB")) +
  labs(title = "Threshold sensitivity (dose-response)",
       subtitle = "HGF+ vs HGF- HSC; enrichment ratio",
       x = "HGF threshold", y = "Enrichment ratio", fill = NULL) +
  theme_pub() + theme(legend.position = "top")

# -------------------------------------------------------- 2e HGF 转录本计数
# 首选从 CosMx 对象实时计算；若对象未落地（iCloud）则回退到已渲染图对应的取值。
counts_tb <- NULL
if (file.exists("data/GSM8855042_Slide1_SeuratObject.RDS")) {
  counts_tb <- tryCatch({
    cosmx1 <- readRDS("data/GSM8855042_Slide1_SeuratObject.RDS")
    ct_col <- grep("clusters$", colnames(cosmx1@meta.data), value = TRUE)[1]
    idx <- which(cosmx1@meta.data[[ct_col]] == "Stellate.cells")
    hc <- as.numeric(cosmx1[["RNA"]]$counts["HGF", idx])
    tb <- data.frame(count = as.numeric(names(table(hc))), n = as.numeric(table(hc)))
    tb$pct <- tb$n / sum(tb$n) * 100
    tb <- tb[tb$count <= 12, ]
    rm(cosmx1); gc(); tb
  }, error = function(e) NULL)
}
if (is.null(counts_tb)) {
  message("⚠️ 未找到 CosMx 对象，改用已渲染图对应的计数值（n = 3,559）")
  counts_tb <- data.frame(count = 0:12,
                          pct = c(58.9, 19.0, 9.8, 4.9, 2.9, 2.1, 1.2, 0.5, 0.4, 0.3, 0.1, 0.0, 0.0))
}
p2e <- ggplot(counts_tb, aes(factor(count), pct)) +
  geom_col(fill = "#C73E1D", alpha = 0.85, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", pct)), vjust = -0.3, size = 2.3, color = "grey30") +
  labs(title = "HGF transcript count in HSCs",
       subtitle = "Zero-inflated long-tail (n = 3,559)",
       x = "HGF count per cell", y = "% of HSCs") +
  coord_cartesian(ylim = c(0, 65)) + theme_pub()

# --------------------------------------------------------------------- 2f SPI
spi <- read.csv("results/SPI_summary_all_HSC_subsets.csv") %>%
  mutate(Subset = factor(Subset, levels = c("HGF+", "HGF-", "PDGFRB+", "PDGFRB-")))
p2f <- ggplot(spi, aes(Subset, Median_SPI, fill = Subset)) +
  geom_col(alpha = 0.9, width = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("HGF+" = "#C73E1D", "HGF-" = "#2E86AB",
                               "PDGFRB+" = "#F18F01", "PDGFRB-" = "#3B8EA5")) +
  labs(title = "Selective Proximity Index", subtitle = "Median SPI per HSC subset",
       x = NULL, y = "Median SPI") +
  theme_pub() + theme(legend.position = "none")

# --------------------------------------------------- 2g/2h 跨患者 + 伪重复校正
a3b <- readRDS("results/revision/A3b_mixed_effects.rds")
pa <- a3b$patient_agg
pa$Sample <- factor(pa$patient_uid, levels = pa$patient_uid[order(pa$delta)])
p2g <- ggplot(pa, aes(delta, Sample)) +
  geom_col(fill = "#C73E1D", alpha = 0.9, width = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_text(aes(label = sprintf("%.1f \u00b5m", delta)), hjust = 1.2, size = 3.5, color = "grey20") +
  labs(title = "Cross-patient consistency (4 independent samples)",
       subtitle = sprintf("Mixed-model P = %.1e (patient RE); %.1e (FOV/patient)", a3b$m_pat$p, a3b$m_nest$p),
       x = "\u0394 median distance (\u00b5m); negative = HGF+ closer", y = NULL) +
  theme_pub() + coord_cartesian(xlim = c(-35, 5), clip = "off")

model_df <- data.frame(
  Model = factor(c("Naive cell-level t-test", "Patient random effect", "FOV nested in Patient"),
                 levels = c("Naive cell-level t-test", "Patient random effect", "FOV nested in Patient")),
  P = c(a3b$naive_t_p, a3b$m_pat$p, a3b$m_nest$p))
model_df$LogP <- -log10(model_df$P)
p2h <- ggplot(model_df, aes(LogP, Model, fill = Model)) +
  geom_col(alpha = 0.9, width = 0.6) +
  geom_text(aes(label = formatC(P, format = "e", digits = 2)), hjust = -0.1, size = 3.2, color = "grey20") +
  scale_fill_manual(values = c("Naive cell-level t-test" = "#F18F01",
                               "Patient random effect" = "#C73E1D",
                               "FOV nested in Patient" = "#2E86AB")) +
  labs(title = "Pseudoreplication correction", subtitle = "Linear mixed-effects models",
       x = "-log10(P)", y = NULL) +
  theme_pub() + theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, max(model_df$LogP) * 1.25))

# ============================ 2i 距离梯度（导入现成 PNG）====================
FIG2I <- "Final/Figure3/Figure2i_distance_gradient.png"
if (!file.exists(FIG2I)) stop("找不到面板 (i) 的图：", FIG2I)
p2i <- wrap_elements(rasterGrob(readPNG(FIG2I), interpolate = TRUE))
cat("已导入面板 (i):", FIG2I, "\n")

# ------------------------------------------------------------------ 拼装 3×3
fig2 <- (p2a | p2b | p2c) / (p2d | p2e | p2f) / (p2g | p2h | p2i) +
  plot_layout(heights = c(1, 1, 0.85)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 11))

ggsave(file.path(OUT, "Figure2_v2.pdf"), fig2, width = 15, height = 14)
ggsave(file.path(OUT, "Figure2_v2.png"), fig2, width = 15, height = 14, dpi = 300)
cat("\n== 已输出 figures/final/Figure2_v2.pdf / .png ==\n")
cat("请打开确认：应为 9 个面板 a–i，(d) 为 HGF 阈值剂量反应，(i) 为距离梯度。\n")
cat("确认无误后改名为 Figure2，替换稿件中的插图与我一并更新。\n")
