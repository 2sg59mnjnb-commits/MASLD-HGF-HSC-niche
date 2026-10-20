# ============================================================================
# 08_external_validation_GSE49541.R
#
# Primary external validation of the 11-gene signature in GSE49541
# (Affymetrix HG-U133 Plus 2.0, n = 72; mild F0-F1 vs advanced F3-F4).
#
# Preprocessing note
# -----------------
# Methods state that expression is z-scored per gene within each cohort. This
# script therefore z-scores the validation cohort internally and applies the
# logistic model fitted on the training cohort. Scaling the validation cohort
# with the training mean/sd instead shifts the AUC slightly (0.773 vs 0.784);
# the within-cohort z-score is used because it removes platform-specific
# location and scale differences between RNA-seq and microarray.
#
# Output
#   results/02_model/external_validation_GSE49541.csv
#   data/derived/GSE49541_validation.rds
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({
  library(GEOquery); library(pROC); library(org.Hs.eg.db)
  library(AnnotationDbi); library(hgu133plus2.db)
})

# ---- 1. Download / load the series matrix ---------------------------------
f <- file.path(DATA_DIR, "GSE49541_series_matrix.txt.gz")
if (!file.exists(f)) {
  dir.create(DATA_DIR, showWarnings = FALSE)
  gse <- getGEO("GSE49541", GSEMatrix = TRUE, getGPL = FALSE, destdir = DATA_DIR)
  eset <- gse[[1]]
} else {
  gse <- getGEO("GSE49541", GSEMatrix = TRUE, getGPL = FALSE, destdir = DATA_DIR)
  eset <- gse[[1]]
}
expr  <- exprs(eset)
pheno <- pData(eset)
cat(sprintf("GSE49541: %d probes x %d samples\n", nrow(expr), ncol(expr)))

# ---- 2. Group from the fibrosis-stage characteristic ----------------------
stage_raw <- pheno[["fibrosis stage:ch1"]]
if (is.null(stage_raw)) {
  cand <- grep("stage|fibro", colnames(pheno), ignore.case = TRUE, value = TRUE)
  stage_raw <- pheno[[cand[1]]]
}
grp <- ifelse(grepl("advanced", stage_raw, ignore.case = TRUE), 1,
              ifelse(grepl("mild", stage_raw, ignore.case = TRUE), 0, NA))
cat("Group distribution:\n"); print(table(grp, use.names = FALSE))
cat(sprintf("  mild (F0-F1): %d, advanced (F3-F4): %d\n",
            sum(grp == 0, na.rm = TRUE), sum(grp == 1, na.rm = TRUE)))

# ---- 3. Probe -> symbol (mean of duplicate probes) ------------------------
sym <- AnnotationDbi::mapIds(hgu133plus2.db, keys = rownames(expr),
                             column = "SYMBOL", keytype = "PROBEID",
                             multiVals = "first")
ok  <- !is.na(sym)
expr <- expr[ok, ]; sym <- sym[ok]

expr_df <- as.data.frame(expr, check.names = FALSE)
expr_df$symbol <- sym
expr_sym <- expr_df %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean), .groups = "drop") %>%
  as.data.frame()
rownames(expr_sym) <- expr_sym$symbol
expr_sym$symbol <- NULL
expr_sym <- as.matrix(expr_sym)
cat(sprintf("After probe collapse: %d genes\n", nrow(expr_sym)))

missing <- setdiff(SIG11, rownames(expr_sym))
cat(sprintf("Signature genes available: %d / 11\n", length(intersect(SIG11, rownames(expr_sym)))))
if (length(missing)) cat("  missing: ", paste(missing, collapse = ", "), "\n")

# ---- 4. Fit on training, predict validation -------------------------------
ml <- load_derived("ML_inputs.rds")
X_tr <- ml$X[, SIG11, drop = FALSE]
y_tr <- ml$y
Xtr_z <- scale(X_tr); Xtr_z[is.na(Xtr_z)] <- 0
fit <- glm(y_tr ~ ., data = as.data.frame(Xtr_z), family = binomial)

# Within-validation-cohort z-score (see note at the top of the script)
X_va <- t(expr_sym[SIG11, , drop = FALSE])
Xva_z <- scale(X_va); Xva_z[is.na(Xva_z)] <- 0
pred <- predict(fit, newdata = as.data.frame(Xva_z), type = "response")

# ---- 5. Discrimination ----------------------------------------------------
keep <- !is.na(grp)
roc_val <- roc(grp[keep], pred[keep], quiet = TRUE)
auc_val <- as.numeric(auc(roc_val))
ci_val  <- as.numeric(ci.auc(roc_val))
cat(sprintf("\nExternal validation AUC = %.3f (95%% CI %.3f-%.3f)\n",
            auc_val, ci_val[1], ci_val[3]))
compare("GSE49541 AUC", auc_val, REPORTED$auc_GSE49541)

# ---- 6. Calibration and Brier --------------------------------------------
brier <- mean((pred[keep] - grp[keep])^2)
cal <- glm(grp[keep] ~ qlogis(pmin(pmax(pred[keep], 1e-10), 1 - 1e-10)),
           family = binomial)
cat(sprintf("Brier score        = %.3f\n", brier))
cat(sprintf("Calibration slope  = %.3f\n", coef(cal)[2]))
cat(sprintf("Calibration intercept = %.3f\n", coef(cal)[1]))

# Sensitivity analysis: scale with training parameters instead
Xva_tr <- scale(X_va, center = attr(Xtr_z, "scaled:center"),
                scale = attr(Xtr_z, "scaled:scale"))
Xva_tr[is.na(Xva_tr)] <- 0
pred_tr <- predict(fit, newdata = as.data.frame(Xva_tr), type = "response")
auc_tr  <- as.numeric(auc(roc(grp[keep], pred_tr[keep], quiet = TRUE)))
cat(sprintf("\nSensitivity: training-parameter scaling gives AUC = %.3f\n", auc_tr))
cat("  (reported value uses within-cohort scaling)\n")

# ---- 7. Save --------------------------------------------------------------
perf <- data.frame(
  metric = c("auc", "ci_low", "ci_high", "brier",
             "calibration_slope", "calibration_intercept",
             "auc_training_scaling", "n", "events"),
  value  = c(auc_val, ci_val[1], ci_val[3], brier,
             coef(cal)[2], coef(cal)[1], auc_tr,
             sum(keep), sum(grp[keep] == 1))
)
write.csv(perf, res_path("02_model", "external_validation_GSE49541.csv"),
          row.names = FALSE)

saveRDS(list(roc = roc_val, auc = auc_val, ci = ci_val, pred = pred,
             group = grp, brier = brier, cal_slope = coef(cal)[2],
             expr = expr_sym, fit = fit),
        file.path(DATA_DIR, "derived", "GSE49541_validation.rds"))

cat("\nNext: 02_model/09_calibration_DCA_Firth.R\n")
