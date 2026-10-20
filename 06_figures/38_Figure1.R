# ============================================================================
# 38_Figure1.R
#
# Figure 1. Spatial landscape of stellate cells in MASLD liver.
#   (a) cell-type composition of CosMx Slide1
#   (b) spatial distribution of stellate cells over the profiled section
#   (c) odds ratio for stellate enrichment in the fibrotic region
#
# Reads the outputs of 03_spatial/13 and 03_spatial/14.
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

comp <- read.csv(res_path("03_spatial", "celltype_composition.csv"))
enr  <- read.csv(res_path("03_spatial", "stellate_enrichment.csv"))
df   <- load_derived("cosmx_flat.rds")
sc   <- load_derived("stellate_enrichment.rds")

# ---- (a) Composition -----------------------------------------------------
comp1 <- comp[comp$slide == "S1", ]
comp1 <- comp1[!grepl("^[a-j]$", comp1$celltype), ]
comp1 <- head(comp1[order(-comp1$n), ], 15)
comp1$highlight <- ifelse(comp1$celltype == "Stellate.cells",
                          "Stellate", "Other")
comp1$celltype <- factor(comp1$celltype, levels = rev(comp1$celltype))

p_a <- ggplot(comp1, aes(celltype, n, fill = highlight)) +
  geom_col(alpha = 0.9) +
  scale_fill_manual(values = c(Stellate = COL_HGFPOS, Other = "grey60")) +
  coord_flip() +
  labs(title = "Cell type composition",
       subtitle = sprintf("CosMx Slide1, n = %s cells",
                          format(sum(comp$n[comp$slide == "S1"]), big.mark = ",")),
       x = NULL, y = "Cell count") +
  theme_pub(9) + theme(legend.position = "none",
                       axis.text.y = element_text(size = 7))

# ---- (b) Spatial distribution -------------------------------------------
plot_df <- df[df$slide == "S1", ]
plot_df$is_stellate <- plot_df$celltype == "Stellate.cells"
plot_df$is_fibrotic <- sc$fibrotic_indep[df$slide == "S1"]

p_b <- ggplot() +
  geom_point(data = plot_df[!plot_df$is_fibrotic, ],
             aes(x, y), size = 0.1, alpha = 0.15, colour = "grey80") +
  geom_point(data = plot_df[plot_df$is_fibrotic, ],
             aes(x, y), size = 0.15, alpha = 0.15, colour = COL_ACCENT) +
  geom_point(data = plot_df[plot_df$is_stellate, ],
             aes(x, y), size = 0.15, alpha = 0.55, colour = COL_HGFPOS) +
  coord_fixed() +
  labs(title = "Stellate cells in fibrotic regions",
       subtitle = sprintf("Orange = fibrotic region; red = stellate cells (n = %s, %.1f%%)",
                          format(sum(plot_df$is_stellate), big.mark = ","),
                          100 * mean(plot_df$is_stellate)),
       x = "x (mm)", y = "y (mm)") +
  theme_pub(9)

# ---- (c) Enrichment OR --------------------------------------------------
prim <- enr[grepl("HSC-independent", enr$definition), ][1, ]

p_c <- ggplot() +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_point(aes(x = prim$OR, y = 1), size = 5, colour = COL_HGFPOS) +
  geom_errorbarh(aes(xmin = prim$CI_low, xmax = prim$CI_high, y = 1),
                 height = 0.15, colour = COL_HGFPOS, linewidth = 1) +
  scale_x_log10() +
  labs(title = "Stellate enrichment in the fibrotic region",
       subtitle = sprintf("OR = %.2f (95%% CI %.2f-%.2f)\nHSC-independent signature",
                          prim$OR, prim$CI_low, prim$CI_high),
       x = "Odds ratio (log scale)", y = NULL) +
  theme_pub(9) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

fig1 <- (p_a | p_b | p_c) + plot_annotation(tag_levels = "a")
safe_save(fig1, "Figure1", width = 15, height = 5)
cat("\nFigure 1 written. Next: 06_figures/39_Figure2.R\n")
