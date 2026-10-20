# ============================================================================
# 42_supplementary_figures.R
#
# Builds the supplementary figures that are generated from tabular outputs of
# the analysis scripts. Figures that require the raw spatial or single-cell
# objects are produced alongside the analysis that computes them; this script
# skips any panel whose input file is missing and prints which ones were skipped,
# so it can be run at any point in the pipeline.
#
# Mapping to the manuscript
#   S1  WGCNA module-trait correlation          -> 01_discovery/03
#   S2  four-algorithm selection                -> 01_discovery/05
#   S3  calibration plot                        -> 02_model/09
#   S4  fixed vs nested CV                      -> 02_model/06, 07
#   S7  inter-sample heterogeneity              -> 03_spatial/16
#   S8  MET pathway sensitivity                 -> results/03_spatial (if run)
#   S9  CHI3L1 model comparison                 -> 02_model/10
#   S11 AUC bootstrap                           -> 02_model/09
#   S12 F0-F4 signature trend                   -> 05_downstream/31
#   S15 Firth vs MLE                            -> 02_model/09
#   S16 calibration bootstrap                   -> 02_model/09
#   S18 F2 subgroup AUC                         -> 02_model/09 (if run)
#   S26 HGF threshold sensitivity               -> 03_spatial/17
#   S27 CHI3L1 ROC overlay                      -> 02_model/10
#   S28 11-gene score in the HGF+ neighbourhood -> results/03_spatial (if run)
#   S30 signature specificity                   -> 04_single_cell/28
#   S31 IL6R model sensitivity                  -> 03_spatial/19
#   S32 depth-corrected distance                -> 03_spatial/20
#   S33 TF activity                             -> 04_single_cell/29
#   S34 motif enrichment                        -> 04_single_cell/30
#   S35 drug repurposing                        -> 05_downstream/36
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(pROC)})

SUPP <- file.path(FIG_DIR, "supplementary")
dir.create(SUPP, showWarnings = FALSE, recursive = TRUE)
skipped <- character(0)

need <- function(...) {
  p <- res_path(...)
  if (!file.exists(p)) { skipped <<- c(skipped, p); return(NULL) }
  p
}

# ---- S1 WGCNA module-trait ----------------------------------------------
f <- need("01_discovery", "WGCNA_module_trait.csv")
if (!is.null(f)) {
  mt <- read.csv(f)
  mt <- mt[mt$module != "ME0", ]
  mt$sig <- ifelse(mt$P < 0.05, "P < 0.05", "n.s.")
  mt <- mt[order(mt$r), ]
  mt$module <- factor(mt$module, levels = mt$module)
  p <- ggplot(mt, aes(r, module, fill = sig)) +
    geom_col(alpha = 0.9) +
    scale_fill_manual(values = c("P < 0.05" = COL_HGFPOS, "n.s." = COL_NEUTRAL)) +
    geom_vline(xintercept = 0, colour = "grey40") +
    labs(title = "WGCNA module-trait correlation",
         subtitle = "Module eigengenes vs fibrosis group",
         x = "Pearson r", y = NULL, fill = NULL) +
    theme_pub(9) + theme(legend.position = "top")
  safe_save(p, "FigureS1_WGCNA", 5, 5, dir = SUPP)
}

# ---- S2 four-algorithm selection ----------------------------------------
f <- need("01_discovery", "feature_selection_genes.csv")
if (!is.null(f)) {
  g <- read.csv(f)
  counts <- data.frame(
    method = c("LASSO", "Random forest", "SVM-RFE", "XGBoost",
               "Intersection (signature)"),
    n = c(sum(g$in_LASSO), sum(g$in_RF), sum(g$in_SVM), sum(g$in_XGB),
          sum(g$n_algorithms == 4))
  )
  counts$type <- ifelse(grepl("Intersection", counts$method),
                        "Consensus", "Individual")
  counts$method <- factor(counts$method, levels = counts$method)
  p <- ggplot(counts, aes(n, method, fill = type)) +
    geom_col(alpha = 0.9) +
    geom_text(aes(label = n), hjust = -0.3, size = 3.2) +
    scale_fill_manual(values = c(Individual = "#92C5DE", Consensus = COL_HGFPOS)) +
    labs(title = "Consensus feature selection",
         subtitle = "Genes selected by each algorithm and by their intersection",
         x = "Number of genes", y = NULL, fill = NULL) +
    theme_pub(9) + theme(legend.position = "top") +
    coord_cartesian(xlim = c(0, max(counts$n) * 1.2))
  safe_save(p, "FigureS2_feature_selection", 6, 4, dir = SUPP)
}

# ---- S3 calibration -----------------------------------------------------
f <- need("02_model", "calibration_GSE49541.csv")
if (!is.null(f)) {
  cal <- read.csv(f)
  perf <- read.csv(res_path("02_model", "external_validation_GSE49541.csv"))
  brier <- perf$value[perf$metric == "brier"]
  p <- ggplot(cal, aes(mean_pred, observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60") +
    geom_point(aes(size = n), colour = COL_HGFPOS, alpha = 0.85) +
    geom_line(colour = COL_HGFPOS) +
    labs(title = "Calibration in GSE49541",
         subtitle = sprintf("Brier score = %.3f", brier),
         x = "Predicted probability", y = "Observed rate", size = "n") +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_pub(10)
  safe_save(p, "FigureS3_calibration", 5, 5, dir = SUPP)
}

# ---- S4 fixed vs nested CV ---------------------------------------------
f1 <- file.path(DATA_DIR, "derived", "fixed_cv.rds")
f2 <- file.path(DATA_DIR, "derived", "nested_cv.rds")
if (file.exists(f1) && file.exists(f2)) {
  a <- readRDS(f1); b <- readRDS(f2)
  d <- rbind(
    data.frame(FPR = 1 - a$roc$specificities, TPR = a$roc$sensitivities,
               method = sprintf("Fixed signature (AUC %.3f)", a$auc)),
    data.frame(FPR = 1 - b$roc$specificities, TPR = b$roc$sensitivities,
               method = sprintf("Fully nested (AUC %.3f)", b$auc))
  )
  p <- ggplot(d, aes(FPR, TPR, colour = method)) +
    geom_line(linewidth = 1.1) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", colour = "grey60") +
    scale_colour_manual(values = c(COL_HGFPOS, COL_HGFNEG)) +
    labs(title = "Fixed versus fully nested cross-validation",
         x = "1 - Specificity", y = "Sensitivity", colour = NULL) +
    coord_fixed() + theme_pub(10) + theme(legend.position = "top")
  safe_save(p, "FigureS4_CV_comparison", 5, 5, dir = SUPP)
} else skipped <- c(skipped, "fixed_cv.rds / nested_cv.rds")

# ---- S7 inter-sample heterogeneity --------------------------------------
f <- need("03_spatial", "per_patient_delta.csv")
if (!is.null(f)) {
  pp <- read.csv(f)
  pp$patient_uid <- factor(pp$patient_uid,
                           levels = pp$patient_uid[order(pp$delta)])
  p <- ggplot(pp, aes(delta, patient_uid)) +
    geom_col(fill = COL_HGFPOS, alpha = 0.9, width = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_text(aes(label = sprintf("%.1f um", delta)), hjust = 1.15, size = 3) +
    labs(title = "Inter-sample heterogeneity",
         subtitle = "Delta median distance (HGF+ minus HGF-) per patient",
         x = "Delta (um)", y = NULL) +
    theme_pub(10) +
    coord_cartesian(xlim = c(min(pp$delta) * 1.3, 3), clip = "off")
  safe_save(p, "FigureS7_inter_sample", 7, 3.5, dir = SUPP)
}

# ---- S9 CHI3L1 comparison ----------------------------------------------
f <- need("02_model", "CHI3L1_comparison.csv")
if (!is.null(f)) {
  cmp <- read.csv(f)
  cmp$model <- factor(cmp$model, levels = rev(cmp$model))
  p <- ggplot(cmp, aes(auc, model)) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2,
                   colour = COL_HGFPOS, linewidth = 0.9) +
    geom_point(size = 3.5, colour = COL_HGFPOS) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey60") +
    geom_text(aes(label = sprintf("%.3f", auc)), vjust = -1.1, size = 3) +
    coord_cartesian(xlim = c(0.5, 1)) +
    labs(title = "CHI3L1 versus multi-gene models",
         subtitle = "External validation (GSE49541)",
         x = "AUC", y = NULL) +
    theme_pub(10)
  safe_save(p, "FigureS9_CHI3L1_comparison", 7, 3.5, dir = SUPP)
}

# ---- S11 AUC bootstrap --------------------------------------------------
f <- need("02_model", "calibration_bootstrap.csv")
if (!is.null(f)) {
  bt <- read.csv(f)
  # AUC bootstrap draws are stored in the rds, not the csv; use the summary here
  p <- ggplot(bt[bt$metric != "auc", ], aes(metric, point)) +
    geom_col(fill = COL_HGFPOS, alpha = 0.9) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15,
                  colour = "grey30") +
    labs(title = "Bootstrap uncertainty of validation metrics",
         subtitle = "2,000 resamples; bars show point estimate and 95% interval",
         x = NULL, y = "Value") +
    theme_pub(10)
  safe_save(p, "FigureS11_metric_bootstrap", 6, 4, dir = SUPP)
}

# ---- S12 F0-F4 trend ---------------------------------------------------
f <- need("05_downstream", "stage_trends_signature.csv")
if (!is.null(f)) {
  st <- load_derived("stage_trends.rds")
  d <- data.frame(stage = factor(st$stage), score = st$score)
  p <- ggplot(d, aes(stage, score, fill = stage)) +
    geom_boxplot(alpha = 0.85, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.4, alpha = 0.3) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                colour = COL_HGFPOS, fill = COL_HGFPOS, alpha = 0.15) +
    scale_fill_brewer(palette = "OrRd") +
    labs(title = "11-gene signature across F0-F4",
         subtitle = "GSE135251 (n = 216)",
         x = "METAVIR stage", y = "Signature score") +
    theme_pub(10) + theme(legend.position = "none")
  safe_save(p, "FigureS12_stage_trend", 7, 4.5, dir = SUPP)
} else skipped <- c(skipped, "stage_trends.csv")

# ---- S15 Firth vs MLE ---------------------------------------------------
f <- need("02_model", "Firth_vs_MLE.csv")
if (!is.null(f)) {
  fm <- read.csv(f)
  fm$model <- factor(fm$model, levels = fm$model)
  p <- ggplot(fm, aes(auc, model, colour = model)) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2,
                   linewidth = 0.9) +
    geom_point(size = 4) +
    scale_colour_manual(values = c(COL_HGFPOS, COL_HGFNEG), guide = "none") +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey60") +
    geom_text(aes(label = sprintf("%.3f", auc)), vjust = -1.1, size = 3) +
    coord_cartesian(xlim = c(0.5, 1)) +
    labs(title = "Firth penalised-likelihood sensitivity analysis",
         subtitle = "External validation, small-sample correction",
         x = "AUC", y = NULL) +
    theme_pub(10)
  safe_save(p, "FigureS15_Firth", 6, 3, dir = SUPP)
}

# ---- S26 HGF threshold effects -----------------------------------------
f <- need("03_spatial", "hgf_threshold_effects.csv")
if (!is.null(f)) {
  th <- read.csv(f)
  p <- ggplot(th, aes(beta, threshold)) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18,
                   colour = COL_HGFPOS, linewidth = 0.9) +
    geom_point(size = 3.5, colour = COL_HGFPOS) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_text(aes(label = sprintf("beta = %.3f", beta)), hjust = 1.15,
              vjust = -0.9, size = 2.8) +
    labs(title = "Sensitivity to the HGF positivity threshold",
         subtitle = "Mixed model with patient random intercept",
         x = "beta (scaled distance)", y = NULL) +
    theme_pub(10) + coord_cartesian(clip = "off")
  safe_save(p, "FigureS26_HGF_thresholds", 7, 3, dir = SUPP)
}

# ---- S27 CHI3L1 ROC overlay --------------------------------------------
f <- file.path(DATA_DIR, "derived", "CHI3L1_comparison.rds")
if (file.exists(f)) {
  cc <- readRDS(f)
  r11 <- cc$rocs[["11-gene signature"]]
  rch <- cc$rocs[["CHI3L1 only"]]
  d <- rbind(
    data.frame(FPR = 1 - r11$specificities, TPR = r11$sensitivities,
               model = sprintf("11-gene (AUC %.3f)", as.numeric(auc(r11)))),
    data.frame(FPR = 1 - rch$specificities, TPR = rch$sensitivities,
               model = sprintf("CHI3L1 (AUC %.3f)", as.numeric(auc(rch))))
  )
  p <- ggplot(d, aes(FPR, TPR, colour = model)) +
    geom_line(linewidth = 1.1) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", colour = "grey60") +
    scale_colour_manual(values = c(COL_HGFPOS, COL_HGFNEG)) +
    labs(title = "11-gene signature versus CHI3L1 alone",
         subtitle = sprintf("DeLong P = %.3f", cc$delong$delong_P),
         x = "1 - Specificity", y = "Sensitivity", colour = NULL) +
    coord_fixed() + theme_pub(10) + theme(legend.position = "top")
  safe_save(p, "FigureS27_CHI3L1_ROC", 5.5, 5.5, dir = SUPP)
} else skipped <- c(skipped, "CHI3L1_comparison.rds")

# ---- S30 signature specificity -----------------------------------------
f <- need("04_single_cell", "signature_score_by_celltype.csv")
if (!is.null(f)) {
  bt <- read.csv(f)
  long <- rbind(
    data.frame(cell = bt$cell.annotation, score = bt$sig11, sig = "11-gene"),
    data.frame(cell = bt$cell.annotation, score = bt$sig16, sig = "16-gene")
  )
  ord <- long[long$sig == "16-gene", ]
  ord <- ord$cell[order(ord$score)]
  long$cell <- factor(long$cell, levels = ord)
  p <- ggplot(long, aes(score, cell, fill = sig)) +
    geom_col(position = position_dodge(0.8), alpha = 0.9, width = 0.7) +
    scale_fill_manual(values = c("11-gene" = COL_HGFNEG, "16-gene" = COL_HGFPOS)) +
    labs(title = "Cell-type specificity of the two signatures",
         subtitle = "Independent scRNA-seq reference (GSE202379)",
         x = "Mean module score", y = NULL, fill = NULL) +
    theme_pub(10) + theme(legend.position = "top",
                          axis.text.y = element_text(size = 7))
  safe_save(p, "FigureS30_signature_specificity", 8, 5, dir = SUPP)
}

# ---- S31 IL6R model sensitivity ----------------------------------------
f <- need("03_spatial", "hep4_il6r_distance_models.csv")
if (!is.null(f)) {
  m <- read.csv(f)
  m <- m[!is.na(m$P), ]
  m$nlog <- -log10(m$P)
  m$Model <- factor(m$Model, levels = m$Model[order(m$nlog)])
  p <- ggplot(m, aes(nlog, Model)) +
    geom_col(fill = COL_HGFPOS, alpha = 0.9, width = 0.6) +
    geom_text(aes(label = formatC(P, format = "e", digits = 1)),
              hjust = -0.1, size = 2.8) +
    labs(title = "Sensitivity of the IL6R proximity result to model choice",
         subtitle = "Only the nested-FOV and cluster-robust models reach P < 0.05",
         x = "-log10(P)", y = NULL) +
    theme_pub(9) + coord_cartesian(xlim = c(0, max(m$nlog) * 1.2))
  safe_save(p, "FigureS31_IL6R_models", 7.5, 3, dir = SUPP)
}

# ---- S32 depth-corrected distance --------------------------------------
f <- need("03_spatial", "depth_corrected_distance.csv")
if (!is.null(f)) {
  d <- read.csv(f)
  d$strategy <- factor(d$strategy, levels = d$strategy[order(d$beta)])
  p <- ggplot(d, aes(beta, strategy)) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2,
                   colour = COL_HGFPOS, linewidth = 0.9) +
    geom_point(size = 3.5, colour = COL_HGFPOS) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_text(aes(label = sprintf("P = %.1e", P)), hjust = 1.15, vjust = -0.9,
              size = 2.8) +
    labs(title = "Depth correction of the HGF+ versus HGF- distance result",
         subtitle = "Covariate adjustment and trimming rather than matching",
         x = "beta (scaled distance)", y = NULL) +
    theme_pub(9) + coord_cartesian(clip = "off")
  safe_save(p, "FigureS32_depth_corrected", 8, 3, dir = SUPP)
}

# ---- S33 TF activity ---------------------------------------------------
f <- need("04_single_cell", "TF_activity_consensus.csv")
if (!is.null(f)) {
  tf <- read.csv(f)
  long <- rbind(
    data.frame(TF = tf$TF, delta = tf$delta_S1, slide = "Slide1"),
    data.frame(TF = tf$TF, delta = tf$delta_S2, slide = "Slide2")
  )
  ord <- tf$TF[order(tf$delta_S1)]
  long$TF <- factor(long$TF, levels = ord)
  p <- ggplot(long, aes(delta, TF, fill = delta)) +
    geom_col(position = position_dodge(0.8), alpha = 0.9, width = 0.7) +
    scale_fill_gradient2(low = COL_HGFNEG, mid = "white", high = COL_HGFPOS,
                         midpoint = 0, guide = "none") +
    geom_vline(xintercept = 0, colour = "grey40") +
    labs(title = "TF activity in HGF+ versus HGF- HSCs",
         subtitle = "DoRothEA A/B/C regulons, concordant across both slides",
         x = "Delta activity", y = NULL) +
    theme_pub(10)
  safe_save(p, "FigureS33_TF_activity", 6.5, 4, dir = SUPP)
}

# ---- S34 motif enrichment ----------------------------------------------
f <- need("04_single_cell", "motif_enrichment.csv")
if (!is.null(f)) {
  mo <- read.csv(f)
  mo$threshold <- factor(mo$threshold, levels = c("80%", "75%", "70%", "65%", "60%"))
  p <- ggplot(mo, aes(threshold, TF, fill = -log10(P_empirical))) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%d\n(%.2fx)", HGF_hits, fold_enrichment)),
              size = 2.4) +
    scale_fill_gradient(low = "grey95", high = COL_HGFPOS, name = "-log10(P)") +
    labs(title = "Motif enrichment in the HGF promoter",
         subtitle = "JASPAR2020, hg38 chr7:81,698-81,702 kb; fold vs 100 shuffled sequences",
         x = "Relative-score threshold", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold"))
  safe_save(p, "FigureS34_motif_enrichment", 8, 4, dir = SUPP)
}

# ---- S35 drug repurposing ----------------------------------------------
f <- need("05_downstream", "DGIdb_drug_summary.csv")
if (!is.null(f)) {
  ds <- read.csv(f)
  ds <- head(ds[order(-ds$n_genes, -ds$max_score), ], 15)
  ds$drug <- factor(ds$drug, levels = ds$drug[order(ds$max_score)])
  p <- ggplot(ds, aes(max_score, drug, fill = n_genes)) +
    geom_col(alpha = 0.9) +
    scale_fill_gradient(low = "#F4A582", high = COL_HGFPOS, name = "genes") +
    labs(title = "Approved drugs targeting the 16-gene HGF+ HSC signature",
         subtitle = "DGIdb query, no score threshold applied",
         x = "Maximum interaction score", y = NULL) +
    theme_pub(9) + theme(axis.text.y = element_text(size = 7))
  safe_save(p, "FigureS35_drug_repurposing", 7, 4.5, dir = SUPP)
}

# ---- Report -------------------------------------------------------------
cat(sprintf("\nSupplementary figures written to %s\n", SUPP))
if (length(skipped)) {
  cat("Skipped (input not found):\n")
  cat(paste0("  ", unique(skipped), collapse = "\n"), "\n")
}
cat("\nNext: 07_tables/43_supplementary_tables.R\n")
