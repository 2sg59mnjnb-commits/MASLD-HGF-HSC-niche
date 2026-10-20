# ============================================================================
# 06_fixed_signature_cv.R
#
# 10-fold cross-validation of the fixed 11-gene signature in the training
# cohort. All preprocessing that depends on the data (gene-wise centring and
# scaling) is fitted inside each training fold and applied to the held-out
# fold, so no information leaks across folds.
#
# Output
#   results/02_model/fixed_signature_cv.csv
#   data/derived/fixed_cv.rds   (per-sample CV predictions, ROC, AUC, CI)
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({
  library(caret)
  library(pROC)
})

ml <- load_derived("ML_inputs.rds")
d  <- load_bulk()

X <- ml$X[, SIG11, drop = FALSE]
y <- as.numeric(d$group == "High")

cat(sprintf("Samples: %d, events: %d, genes: %d\n", nrow(X), sum(y), ncol(X)))

# ---- Cross-validation ------------------------------------------------------
set.seed(SEED)
folds <- createFolds(y, k = 10, list = TRUE)
cv_pred <- rep(NA_real_, length(y))

for (i in seq_along(folds)) {
  test  <- folds[[i]]
  train <- setdiff(seq_along(y), test)

  # Within-fold scaling
  Xtr <- scale(X[train, , drop = FALSE])
  mu  <- attr(Xtr, "scaled:center")
  sdv <- attr(Xtr, "scaled:scale")
  Xtr[is.na(Xtr)] <- 0

  Xte <- scale(X[test, , drop = FALSE], center = mu, scale = sdv)
  Xte[is.na(Xte)] <- 0

  df_tr <- as.data.frame(Xtr); df_tr$y <- y[train]
  fit <- glm(y ~ ., data = df_tr, family = binomial)
  cv_pred[test] <- predict(fit, newdata = as.data.frame(Xte), type = "response")
}

roc_cv <- roc(y, cv_pred, quiet = TRUE)
auc_cv <- as.numeric(auc(roc_cv))
ci_cv  <- as.numeric(ci.auc(roc_cv))

cat(sprintf("\nFixed-signature CV AUC = %.3f (95%% CI %.3f-%.3f)\n",
            auc_cv, ci_cv[1], ci_cv[3]))
compare("CV AUC (fixed signature)", auc_cv, REPORTED$auc_cv_fixed)

# Full-data model, reported for completeness but expected to be optimistic.
X_all <- scale(X); X_all[is.na(X_all)] <- 0
df_all <- as.data.frame(X_all); df_all$y <- y
fit_all <- glm(y ~ ., data = df_all, family = binomial)
auc_full <- as.numeric(auc(roc(y, predict(fit_all, type = "response"), quiet = TRUE)))
cat(sprintf("Full-data AUC (optimistic, not reported in the manuscript) = %.3f\n", auc_full))

# ---- Coefficients of the deployed model -----------------------------------
coef_tab <- summary(fit_all)$coefficients
coef_df <- data.frame(
  term        = rownames(coef_tab),
  coefficient = round(coef_tab[, "Estimate"], 4),
  std_error   = round(coef_tab[, "Std. Error"], 4),
  z           = round(coef_tab[, "z value"], 3),
  P           = coef_tab[, "Pr(>|z|)"],
  OR          = round(exp(coef_tab[, "Estimate"]), 3),
  row.names   = NULL
)
cat("\nModel coefficients (scaled features):\n"); print(coef_df)
write.csv(coef_df, res_path("02_model", "model_coefficients.csv"), row.names = FALSE)

# Training-set scaling parameters, needed to apply the model elsewhere.
write.csv(data.frame(gene = SIG11,
                     mean = as.numeric(attr(X_all, "scaled:center")),
                     sd   = as.numeric(attr(X_all, "scaled:scale"))),
          res_path("02_model", "training_scaling_parameters.csv"), row.names = FALSE)

# ---- Save -----------------------------------------------------------------
perf <- data.frame(
  metric = c("cv_auc", "cv_ci_low", "cv_ci_high",
             "full_auc", "n", "events", "n_genes"),
  value  = c(auc_cv, ci_cv[1], ci_cv[3], auc_full, length(y), sum(y), length(SIG11))
)
write.csv(perf, res_path("02_model", "fixed_signature_cv.csv"), row.names = FALSE)

saveRDS(list(cv_pred = cv_pred, y = y, roc = roc_cv, auc = auc_cv, ci = ci_cv,
             fit_full = fit_all, folds = folds),
        file.path(DATA_DIR, "derived", "fixed_cv.rds"))

cat("\nNext: 02_model/07_fully_nested_cv.R\n")
