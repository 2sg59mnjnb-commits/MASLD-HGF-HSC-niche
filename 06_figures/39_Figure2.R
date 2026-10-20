# ============================================================================
# 39_Figure2.R
#
# Figure 2. HGF+ HSC spatial niche and cross-slide consistency (panels a-h).
#   a  niche scan, log2 enrichment ratio, both slides
#   b  FOV-level direction consistency
#   c  cross-slide agreement of niche enrichment
#   d  distance-threshold sensitivity of the median distance
#   e  HGF transcript counts per HSC
#   f  Selective Proximity Index by HSC subset
#   g  per-patient median distance difference
#   h  mixed-model estimates of the HGF+ vs HGF- effect
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(ggrepel)})

# ---- (a) Niche scan ------------------------------------------------------
niche <- read.csv(res_path("03_spatial", "niche_scan_by_slide.csv"))
niche <- niche[niche$n_target >= 200, ]
niche$log2_ratio[!is.finite(niche$log2_ratio)] <- NA
order_ct <- niche %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(m = mean(log2_ratio, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(m)
niche$celltype <- factor(niche$celltype, levels = order_ct$celltype)

p_a <- ggplot(niche, aes(celltype, log2_ratio, fill = slide)) +
  geom_col(position = position_dodge(0.8), alpha = 0.9) +
  scale_fill_manual(values = c(S1 = COL_HGFPOS, S2 = COL_HGFNEG)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  coord_flip() +
  labs(title = "HGF+ HSC niche scan (50 um)",
       subtitle = "log2 enrichment ratio across annotated cell types",
       x = NULL, y = "log2 ratio") +
  theme_pub(9) + theme(legend.position = "top",
                       axis.text.y = element_text(size = 6))

# ---- (b) FOV consistency -------------------------------------------------
fov <- read.csv(res_path("03_spatial", "fov_direction_consistency.csv"))
fov_df <- data.frame(
  direction = factor(c("HGF+ closer", "Other"),
                     levels = c("HGF+ closer", "Other")),
  n = c(sum(fov$direction == "HGF+ closer"), sum(fov$direction != "HGF+ closer"))
)
fov_df$pct <- 100 * fov_df$n / sum(fov_df$n)

p_b <- ggplot(fov_df, aes(direction, n, fill = direction)) +
  geom_col(alpha = 0.9, width = 0.6) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n, pct)), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("HGF+ closer" = COL_HGFPOS, Other = COL_NEUTRAL)) +
  labs(title = "FOV-level direction consistency",
       subtitle = sprintf("Binomial P = %.1e; n = %d FOVs",
                          binom.test(fov_df$n[1], sum(fov_df$n),
                                     alternative = "greater")$p.value,
                          sum(fov_df$n)),
       x = NULL, y = "Number of FOVs") +
  theme_pub(9) + theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, max(fov_df$n) * 1.15))

# ---- (c) Cross-slide -----------------------------------------------------
cs <- read.csv(res_path("03_spatial", "niche_scan_cross_slide.csv"))
ct_cs <- cor.test(cs$log2_ratio_S1, cs$log2_ratio_S2)

p_c <- ggplot(cs, aes(log2_ratio_S1, log2_ratio_S2)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, colour = "grey80", linetype = "dotted") +
  geom_point(size = 2.2, alpha = 0.85, colour = COL_HGFPOS) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", fill = "grey80") +
  ggrepel::geom_text_repel(aes(label = celltype), size = 2, max.overlaps = 20) +
  labs(title = "Cross-slide niche consistency",
       subtitle = sprintf("Pearson r = %.3f, P = %.1e",
                          ct_cs$estimate, ct_cs$p.value),
       x = "log2 ratio, Slide1", y = "log2 ratio, Slide2") +
  coord_fixed() + theme_pub(9)

# ---- (d) Distance-threshold sensitivity ----------------------------------
thr <- read.csv(res_path("03_spatial", "distance_threshold_sensitivity.csv"))
thr_long <- rbind(
  data.frame(threshold_um = thr$threshold_um, group = "HGF+",
             median_um = thr$median_hgfpos),
  data.frame(threshold_um = thr$threshold_um, group = "HGF-",
             median_um = thr$median_hgfneg)
)

p_d <- ggplot(thr_long, aes(threshold_um, median_um, colour = group)) +
  geom_line(linewidth = 1.1) + geom_point(size = 3) +
  scale_colour_manual(values = c("HGF+" = COL_HGFPOS, "HGF-" = COL_HGFNEG)) +
  scale_x_continuous(breaks = thr$threshold_um) +
  labs(title = "Threshold sensitivity of the niche",
       subtitle = "Median distance to nearest Hep.4 within each distance cap",
       x = "Distance threshold (um)", y = "Median distance (um)", colour = NULL) +
  theme_pub(9) + theme(legend.position = "top")

# ---- (e) HGF transcript counts ------------------------------------------
# Raw counts are used here (not the normalised data layer) because the panel
# shows the per-cell transcript count distribution.
obj1 <- load_cosmx(1)
ct1  <- as.character(obj1@meta.data[[celltype_col(1)]])
hsc1_counts <- as.numeric(GetAssayData(obj1, assay = "RNA", layer = "counts")[
  "HGF", ct1 == "Stellate.cells"])
rm(obj1); gc(verbose = FALSE)

tab <- as.data.frame(table(counts = pmin(hsc1_counts, 12)))
tab$pct <- 100 * tab$Freq / sum(tab$Freq)

p_e <- ggplot(tab, aes(factor(counts), pct)) +
  geom_col(fill = COL_HGFPOS, alpha = 0.85, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", pct)), vjust = -0.3, size = 2.2,
            colour = "grey30") +
  labs(title = "HGF transcript counts in HSCs",
       subtitle = sprintf("Zero-inflated long tail (n = %d)", length(hsc1_counts)),
       x = "HGF count per cell (capped at 12)", y = "% of HSCs") +
  coord_cartesian(ylim = c(0, max(tab$pct) * 1.15)) +
  theme_pub(9)

# ---- (f) SPI by subset ---------------------------------------------------
spi <- read.csv(res_path("03_spatial", "spi_by_subset.csv"))
spi$subset <- factor(spi$subset, levels = c("HGF+", "HGF-", "PDGFRB+", "PDGFRB-"))

p_f <- ggplot(spi, aes(subset, median_spi, fill = subset)) +
  geom_col(alpha = 0.9, width = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("HGF+" = COL_HGFPOS, "HGF-" = COL_HGFNEG,
                               "PDGFRB+" = COL_ACCENT, "PDGFRB-" = "#3B8EA5")) +
  labs(title = "Selective Proximity Index",
       subtitle = "log2(d_Hep5 / d_Hep4); higher = closer to Hep.4",
       x = NULL, y = "Median SPI") +
  theme_pub(9) + theme(legend.position = "none")

# ---- (g) Per-patient difference -----------------------------------------
pp <- read.csv(res_path("03_spatial", "per_patient_delta.csv"))
pp$sample <- factor(pp$patient_uid, levels = pp$patient_uid[order(pp$delta)])

p_g <- ggplot(pp, aes(delta, sample)) +
  geom_col(fill = COL_HGFPOS, alpha = 0.9, width = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = sprintf("%.1f", delta)), hjust = 1.2, size = 3.2,
            colour = "grey20") +
  labs(title = "Per-patient median distance difference",
       subtitle = "Negative = HGF+ HSC closer to Hep.4",
       x = "Delta median distance (um)", y = NULL) +
  theme_pub(9) + coord_cartesian(xlim = c(min(pp$delta) * 1.25, 3), clip = "off")

# ---- (h) Mixed models ----------------------------------------------------
mm <- read.csv(res_path("03_spatial", "hgfpos_proximity_models.csv"))
mm <- mm[!is.na(mm$P), ]
mm$neglog10P <- -log10(mm$P)
mm$Model <- factor(mm$Model, levels = mm$Model[order(mm$neglog10P)])

p_h <- ggplot(mm, aes(neglog10P, Model, fill = Model)) +
  geom_col(alpha = 0.9, width = 0.6) +
  geom_text(aes(label = formatC(P, format = "e", digits = 2)),
            hjust = -0.1, size = 2.6, colour = "grey20") +
  scale_fill_manual(values = c("Naive cell-level t-test (biased)" = COL_ACCENT,
                               "Patient random intercept" = COL_HGFPOS,
                               "FOV random intercept" = "#A23B72",
                               "FOV nested in patient" = COL_HGFNEG,
                               "FOV cluster-robust SE" = "#3B8EA5")) +
  labs(title = "Pseudoreplication correction",
       subtitle = "Linear mixed-effects models and robust SE",
       x = "-log10(P)", y = NULL) +
  theme_pub(9) + theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, max(mm$neglog10P) * 1.25))

fig2 <- (p_a | p_b | p_c) / (p_d | p_e | p_f) / (p_g | p_h | plot_spacer()) +
  plot_layout(heights = c(1, 1, 0.8)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 11))

safe_save(fig2, "Figure2", width = 15, height = 14)
cat("\nFigure 2 written. Next: 06_figures/40_Figure3.R\n")
