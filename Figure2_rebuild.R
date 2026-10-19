# =============================================================================
# Figure 2 -- HGF+ HSC spatial niche (nine-panel display figure)
#
# Input
#   results/HGF_HSC_niche_scan.csv, results/HGF_HSC_niche_scan_Slide2.csv
#   results/Reverse_FOV_characterization.rds, results/niche_scan_cross_slide.csv
#   results/HGF_threshold_sensitivity.csv
#   results/SPI_summary_all_HSC_subsets.csv, results/revision/A3b_mixed_effects.rds
#   Final/Figure3/Figure2i_distance_gradient_HGFpos.csv / _HGFneg.csv
#   data/GSM8855042_Slide1_SeuratObject.RDS (optional; only for panel e)
# Output
#   figures/final/Figure2_v2.pdf / .png
#
# Panels: a niche scan; b FOV-level direction consistency; c cross-slide
# consistency; d HGF-positivity threshold sensitivity; e HGF transcript counts;
# f Selective Proximity Index; g cross-patient consistency; h pseudoreplication
# correction; i distance-gradient enrichment.
#
# Requires the intermediate result tables produced by analysis_code.R.
# Usage: cd ~/Desktop/sjfx3 && Rscript Figure2_rebuild.R
# =============================================================================

setwd("~/Desktop/sjfx3")

# ggrepel places labels with a stochastic repulsion algorithm. Fix the seed so
# that repeated runs produce an identical figure.
set.seed(20240501)
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
  library(ggrepel)
})

OUT <- "figures/final"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
theme_pub <- function(bs = 9) theme_classic(base_size = bs)

# ---------------------------------------------------------------------------
# Saving. The figures published in the manuscript were produced by exactly these
# calls, so running this script reproduces them. NOTE: if your R installation has
# no working cairo device, the stock PDF device silently drops characters it
# cannot map to its encoding (on the original analysis machine the Greek Delta in
# panel g and in Figure 3 panel b was dropped from the PDF while the PNG rendered
# it correctly). The PNG is authoritative; a faithful vector PDF requires cairo
# (e.g. install XQuartz, or the Cairo R package).
# ---------------------------------------------------------------------------
save_fig <- function(plot, stem, w, h) {
  ggsave(paste0(stem, ".pdf"), plot, width = w, height = h)
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 300)
}

# ------------------------------------------------- 2a  niche scan
res1 <- read.csv("results/HGF_HSC_niche_scan.csv");        res1$Slide <- "Slide1"
res2 <- read.csv("results/HGF_HSC_niche_scan_Slide2.csv");  res2$Slide <- "Slide2"
niche_df <- rbind(res1, res2) %>% filter(N_total >= 200) %>% mutate(Score = log2(Ratio + 1e-6))

p2a <- ggplot(niche_df, aes(reorder(CellType, Score), Score, fill = Slide)) +
  geom_col(position = position_dodge(0.8), alpha = 0.9) +
  scale_fill_manual(values = c("Slide1" = "#C73E1D", "Slide2" = "#2E86AB")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  coord_flip() +
  labs(title = "HGF+ HSC niche scan (50 \u00b5m)",
       subtitle = "log2(enrichment ratio) across two slides",
       x = NULL, y = "log2 ratio") +
  theme_pub() + theme(legend.position = "top", axis.text.y = element_text(size = 7))

# ------------------------------------------------- 2b  FOV-level consistency
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

# ------------------------------------------------- 2c  cross-slide consistency
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

# ============ 2d  [KEY CORRECTION] HGF-positivity threshold sensitivity ======
# The old code read results/revision/A4_threshold_sensitivity.csv here and
# plotted "distance threshold vs median distance", which matched neither the
# legend nor the main text. Restored to the correct version below.
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

# ------------------------------------------------- 2e  HGF transcript counts
# Preferred: compute from the CosMx object. If that object is not available
# locally, fall back to the counts underlying the published rendering.
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
  message("NOTE: CosMx object not found; using the counts underlying the published rendering (n = 3,559)")
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

# ------------------------------------------------- 2f  Selective Proximity Index
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

# ------------------------------------------- 2g/2h  cross-patient + pseudoreplication
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

# ============ 2i  distance-gradient enrichment ==============================
# Rebuilt from the underlying result tables rather than pasted in as a rendered
# image. This matches the geometry of the other panels, renders as vector in the
# PDF, and removes the previous dependency on a pre-rendered asset.
dgp <- read.csv("Final/Figure3/Figure2i_distance_gradient_HGFpos.csv")
dgn <- read.csv("Final/Figure3/Figure2i_distance_gradient_HGFneg.csv")
dgi <- rbind(dgp, dgn)
dgi$group  <- factor(dgi$group, levels = c("HGF- Stellate", "HGF+ Stellate"))
dgi$radius <- factor(dgi$radius, levels = c(25, 50, 75, 100, 150))
lab <- subset(dgi, group == "HGF+ Stellate")
lab$txt <- ifelse(lab$p_value <= 0.001, "P=0.001", sprintf("P=%.3f", lab$p_value))

p2i <- ggplot(dgi, aes(radius, enrichment, colour = group, group = group)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.9) +
  geom_text(data = lab, aes(label = txt), vjust = -1.0, size = 2.2, colour = "grey25") +
  scale_colour_manual(values = c("HGF- Stellate" = "#2E86AB", "HGF+ Stellate" = "#C73E1D")) +
  labs(title = "Distance-gradient enrichment of Hep.4 around HSC subsets",
       x = "Radius (\u00b5m)", y = "Hep.4 enrichment ratio (observed / expected)",
       colour = NULL) +
  theme_pub(9) +
  theme(legend.position = "top", legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(size = 8.5))

# ------------------------------------------------- assemble 3x3
fig2 <- (p2a | p2b | p2c) / (p2d | p2e | p2f) / (p2g | p2h | p2i) +
  plot_layout(heights = c(1, 1, 0.85)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 11))

save_fig(fig2, file.path(OUT, "Figure2_v2"), 15, 14)
cat("\n== Wrote figures/final/Figure2_v2.pdf / .png ==\n")
cat("Please verify: nine panels a-i; (d) = HGF threshold dose-response; (i) = distance gradient.\n")
