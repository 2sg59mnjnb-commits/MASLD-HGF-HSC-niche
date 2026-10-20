# ============================================================================
# 41_Figure4_Figure5.R
#
# Figure 4. IL6-family receptor positive rates in HGF+ versus HGF- HSCs,
#           separately for the two CosMx slides.
# Figure 5. Diagnostic performance of the 11-gene signature:
#           (a) training cross-validation ROC
#           (b) external validation ROC (GSE49541)
#           (c) decision curve analysis
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(pROC)})

# ============================ Figure 4 ====================================
rec <- read.csv(res_path("04_single_cell", "IL6_receptor_positive_rates.csv"))
rec_long <- rbind(
  data.frame(slide = rec$slide, gene = rec$gene, group = "HGF+",
             value = rec$positive_rate_hgfpos),
  data.frame(slide = rec$slide, gene = rec$gene, group = "HGF-",
             value = rec$positive_rate_hgfneg)
)
rec_long$group <- factor(rec_long$group, levels = c("HGF-", "HGF+"))
rec_long$slide <- factor(rec_long$slide, levels = c("S1", "S2"),
                         labels = c("Slide1", "Slide2"))

p4 <- ggplot(rec_long, aes(gene, value, fill = group)) +
  geom_col(position = position_dodge(0.8), alpha = 0.9, width = 0.75) +
  facet_wrap(~ slide) +
  scale_fill_manual(values = c("HGF+" = COL_HGFPOS, "HGF-" = COL_HGFNEG)) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +
  labs(title = "IL6-family receptor positive rates in HSC subsets",
       subtitle = "CosMx Slide1 and Slide2; positive rate = fraction of cells with transcript > 0",
       x = NULL, y = "Positive rate", fill = NULL) +
  theme_pub(11) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold", size = 11))

safe_save(p4, "Figure4", width = 8, height = 4.5)

# ============================ Figure 5 ====================================
cv <- load_derived("fixed_cv.rds")
va <- load_derived("GSE49541_validation.rds")
dca <- read.csv(res_path("02_model", "DCA_GSE49541.csv"))

roc_cv <- cv$roc
roc_va <- va$roc

df_cv <- data.frame(FPR = 1 - roc_cv$specificities, TPR = roc_cv$sensitivities)
df_va <- data.frame(FPR = 1 - roc_va$specificities, TPR = roc_va$sensitivities)

p5a <- ggplot(df_cv, aes(FPR, TPR)) +
  geom_line(colour = COL_HGFPOS, linewidth = 1.2) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", colour = "grey60") +
  annotate("text", x = 0.62, y = 0.12, hjust = 0,
           label = sprintf("CV AUC = %.3f\n95%% CI %.3f-%.3f",
                           cv$auc, cv$ci[1], cv$ci[3]),
           size = 3.2, colour = COL_HGFPOS) +
  labs(title = "Training performance",
       subtitle = "GSE135251, 10-fold cross-validation",
       x = "1 - Specificity", y = "Sensitivity") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_pub(10)

p5b <- ggplot(df_va, aes(FPR, TPR)) +
  geom_line(colour = COL_HGFPOS, linewidth = 1.2) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", colour = "grey60") +
  annotate("text", x = 0.62, y = 0.12, hjust = 0,
           label = sprintf("AUC = %.3f\n95%% CI %.3f-%.3f",
                           va$auc, va$ci[1], va$ci[3]),
           size = 3.2, colour = COL_HGFPOS) +
  labs(title = "External validation",
       subtitle = sprintf("GSE49541, n = %d (%d events)",
                          sum(!is.na(va$group)), sum(va$group == 1, na.rm = TRUE)),
       x = "1 - Specificity", y = "Sensitivity") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_pub(10)

dca_mod <- dca[dca$variable %in% c("pred_prob", "all", "none"), ]
dca_mod$variable <- factor(dca_mod$variable,
                           levels = c("pred_prob", "all", "none"),
                           labels = c("11-gene signature", "Treat all", "Treat none"))

p5c <- ggplot(dca_mod, aes(threshold, net_benefit, colour = variable)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c("11-gene signature" = COL_HGFPOS,
                                 "Treat all" = "grey40",
                                 "Treat none" = "black")) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  coord_cartesian(xlim = c(0, 0.6)) +
  labs(title = "Decision curve analysis",
       subtitle = "External validation cohort",
       x = "Threshold probability", y = "Net benefit", colour = NULL) +
  theme_pub(10) +
  theme(legend.position = c(0.72, 0.85),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 8))

fig5 <- (p5a | p5b | p5c) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 11))

safe_save(fig5, "Figure5", width = 13, height = 4.5)

cat("\nFigures 4 and 5 written. Next: 06_figures/42_supplementary_figures.R\n")
