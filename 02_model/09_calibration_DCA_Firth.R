# ============================================================================
# 09_calibration_DCA_Firth.R
#
# Secondary performance analyses for the external validation cohort:
#   (1) calibration (Brier score, slope, intercept, binned observed vs expected)
#   (2) decision curve analysis (net benefit across threshold probabilities)
#   (3) Firth penalised-likelihood refit, as a small-sample sensitivity analysis
#   (4) bootstrap confidence intervals for AUC, Brier and calibration slope
#
# With 32 events and 11 predictors the events-per-variable ratio is 2.9, so the
# maximum-likelihood estimates are unstable. Firth's penalty is reported
# alongside, not instead of, the primary model.
#
# Output
#   results/02_model/calibration_GSE49541.csv
#   results/02_model/DCA_GSE49541.csv
#   results/02_model/Firth_vs_MLE.csv
#   results/02_model/calibration_bootstrap.csv
#   data/derived/secondary_performance.rds
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({
  library(pROC); library(dcurves); library(logistf)
})

val <- load_derived("GSE49541_validation.rds")
ml  <- load_derived("ML_inputs.rds")

keep <- !is.na(val$group)
obs  <- val$group[keep]
pred <- val$pred[keep]
cat(sprintf("Validation set: n = %d, events = %d\n", length(obs), sum(obs)))

# ---- 1. Calibration -------------------------------------------------------
brier <- mean((pred - obs)^2)
cal   <- glm(obs ~ qlogis(pmin(pmax(pred, 1e-10), 1 - 1e-10)), family = binomial)
slope <- unname(coef(cal)[2]); intercept <- unname(coef(cal)[1])
cat(sprintf("\nBrier = %.3f | calibration slope = %.3f | intercept = %.3f\n",
            brier, slope, intercept))

bins <- cut(pred, breaks = seq(0, 1, 0.2), include.lowest = TRUE)
cal_df <- data.frame(pred = pred, obs = obs, bin = bins) %>%
  dplyr::group_by(bin) %>%
  dplyr::summarise(mean_pred = mean(pred), observed = mean(obs),
                   n = dplyr::n(), .groups = "drop")
print(as.data.frame(cal_df))
write.csv(cal_df, res_path("02_model", "calibration_GSE49541.csv"), row.names = FALSE)

# ---- 2. Decision curve analysis -------------------------------------------
dca_obj <- dcurves::dca(outcome ~ pred_prob,
                        data = data.frame(outcome = obs, pred_prob = pred))
dca_tab <- as.data.frame(dca_obj$dca)
dca_tab <- dca_tab[round(dca_tab$threshold, 2) %in% seq(0, 0.6, 0.05), ]
write.csv(dca_tab, res_path("02_model", "DCA_GSE49541.csv"), row.names = FALSE)

nb_sel <- dca_tab[dca_tab$variable == "pred_prob" &
                    round(dca_tab$threshold, 2) %in% c(0.2, 0.3, 0.4), ]
cat("\nNet benefit of the signature:\n")
print(nb_sel[, c("threshold", "net_benefit")])

# ---- 3. Firth penalised logistic ------------------------------------------
X_tr <- ml$X[, SIG11, drop = FALSE]
Xtr_z <- scale(X_tr); Xtr_z[is.na(Xtr_z)] <- 0
df_tr <- as.data.frame(Xtr_z); df_tr$y <- ml$y

fit_firth <- logistf(y ~ ., data = df_tr)

expr_val <- t(val$expr[SIG11, , drop = FALSE])
Xva_z <- scale(expr_val); Xva_z[is.na(Xva_z)] <- 0
pred_firth <- as.numeric(predict(fit_firth, newdata = as.data.frame(Xva_z),
                                 type = "response"))

roc_firth <- roc(obs, pred_firth, quiet = TRUE)
auc_firth <- as.numeric(auc(roc_firth))
ci_firth  <- as.numeric(ci.auc(roc_firth))
brier_firth <- mean((pred_firth - obs)^2)
cal_firth <- glm(obs ~ qlogis(pmin(pmax(pred_firth, 1e-10), 1 - 1e-10)),
                 family = binomial)

cat(sprintf("\nFirth: AUC = %.3f (95%% CI %.3f-%.3f), Brier = %.3f, slope = %.3f\n",
            auc_firth, ci_firth[1], ci_firth[3], brier_firth, coef(cal_firth)[2]))

firth_tab <- data.frame(
  model = c("MLE", "Firth"),
  auc   = c(as.numeric(val$auc), auc_firth),
  ci_low  = c(val$ci[1], ci_firth[1]),
  ci_high = c(val$ci[3], ci_firth[3]),
  brier = c(brier, brier_firth),
  cal_slope = c(slope, unname(coef(cal_firth)[2])),
  cal_intercept = c(intercept, unname(coef(cal_firth)[1]))
)
print(firth_tab)
write.csv(firth_tab, res_path("02_model", "Firth_vs_MLE.csv"), row.names = FALSE)

# ---- 4. Bootstrap ---------------------------------------------------------
set.seed(SEED)
n_boot <- 2000
boot <- matrix(NA_real_, n_boot, 4,
               dimnames = list(NULL, c("auc", "brier", "slope", "intercept")))

for (i in seq_len(n_boot)) {
  idx <- sample(seq_along(obs), replace = TRUE)
  if (length(unique(obs[idx])) < 2) next
  d <- data.frame(obs = obs[idx], pred = pred[idx])
  boot[i, "auc"]   <- as.numeric(auc(roc(d$obs, d$pred, quiet = TRUE)))
  boot[i, "brier"] <- mean((d$pred - d$obs)^2)
  cf <- tryCatch(glm(obs ~ qlogis(pmin(pmax(pred, 1e-10), 1 - 1e-10)),
                     data = d, family = binomial), error = function(e) NULL)
  if (!is.null(cf)) {
    boot[i, "slope"]     <- coef(cf)[2]
    boot[i, "intercept"] <- coef(cf)[1]
  }
}

boot_tab <- data.frame(
  metric = colnames(boot),
  point  = c(as.numeric(val$auc), brier, slope, intercept),
  ci_low  = apply(boot, 2, quantile, 0.025, na.rm = TRUE),
  ci_high = apply(boot, 2, quantile, 0.975, na.rm = TRUE),
  row.names = NULL
)
print(boot_tab)
write.csv(boot_tab, res_path("02_model", "calibration_bootstrap.csv"),
          row.names = FALSE)

saveRDS(list(calibration = cal_df, dca = dca_tab, firth = firth_tab,
             boot = boot, boot_summary = boot_tab,
             pred_firth = pred_firth, roc_firth = roc_firth),
        file.path(DATA_DIR, "derived", "secondary_performance.rds"))

cat("\nNext: 02_model/10_CHI3L1_comparison.R\n")
