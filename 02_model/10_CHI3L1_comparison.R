# ============================================================================
# 10_CHI3L1_comparison.R
#
# Head-to-head comparison of the 11-gene signature against CHI3L1 alone, and
# against two two-gene models, in the external validation cohort.
#
# The point of this analysis is negative: CHI3L1 alone performs at least as
# well as the full signature, so the signature is reported as a distributed
# fibrosis proxy rather than as a CHI3L1-independent biomarker.
#
# All four models use identical preprocessing (within-cohort z-score, logistic
# regression fitted on the training cohort) so the AUCs are comparable.
#
# Output
#   results/02_model/CHI3L1_comparison.csv
#   results/02_model/CHI3L1_DeLong.csv
#   data/derived/CHI3L1_comparison.rds
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(pROC)})

ml  <- load_derived("ML_inputs.rds")
val <- load_derived("GSE49541_validation.rds")

obs  <- val$group
keep <- !is.na(obs)
X_tr <- ml$X
y_tr <- ml$y
expr_val <- val$expr

models <- list(
  "CHI3L1 only"        = c("CHI3L1"),
  "CHI3L1 + ILDR1"     = c("CHI3L1", "ILDR1"),
  "CHI3L1 + SULF2"     = c("CHI3L1", "SULF2"),
  "11-gene signature"  = SIG11
)

# Fit on training, predict validation, with within-cohort scaling on both sides.
fit_predict <- function(genes) {
  genes <- genes[genes %in% colnames(X_tr) & genes %in% rownames(expr_val)]
  if (length(genes) == 0) return(rep(NA_real_, ncol(expr_val)))

  Xtr <- scale(X_tr[, genes, drop = FALSE]); Xtr[is.na(Xtr)] <- 0
  df  <- as.data.frame(Xtr); df$y <- y_tr
  fit <- glm(y ~ ., data = df, family = binomial)

  Xva <- scale(t(expr_val[genes, , drop = FALSE])); Xva[is.na(Xva)] <- 0
  as.numeric(predict(fit, newdata = as.data.frame(Xva), type = "response"))
}

res <- list(); rocs <- list()
for (nm in names(models)) {
  p <- fit_predict(models[[nm]])
  r <- roc(obs[keep], p[keep], quiet = TRUE)
  res[[nm]]  <- c(auc = as.numeric(auc(r)), ci.auc(r))
  rocs[[nm]] <- r
  cat(sprintf("%-20s AUC = %.3f (95%% CI %.3f-%.3f)\n", nm,
              res[[nm]][1], res[[nm]][2], res[[nm]][3]))
}

# ---- DeLong and bootstrap delta AUC for 11-gene vs CHI3L1 ------------------
dl <- roc.test(rocs[["11-gene signature"]], rocs[["CHI3L1 only"]],
               method = "delong")
cat(sprintf("\nDeLong P (11-gene vs CHI3L1) = %.3f\n", dl$p.value))

set.seed(SEED)
n_boot <- 2000
delta <- numeric(n_boot)
p11  <- fit_predict(SIG11)
pchi <- fit_predict("CHI3L1")
for (i in seq_len(n_boot)) {
  idx <- sample(which(keep), replace = TRUE)
  if (length(unique(obs[idx])) < 2) { delta[i] <- NA; next }
  a1 <- as.numeric(auc(roc(obs[idx], p11[idx],  quiet = TRUE)))
  a2 <- as.numeric(auc(roc(obs[idx], pchi[idx], quiet = TRUE)))
  delta[i] <- a1 - a2
}
ci_delta <- quantile(delta, c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("Delta AUC (11-gene - CHI3L1) = %.3f (95%% CI %.3f to %.3f)\n",
            mean(delta, na.rm = TRUE), ci_delta[1], ci_delta[2]))

if (ci_delta[1] <= 0 && ci_delta[2] >= 0) {
  cat("  -> CI includes 0: the signature is not superior to CHI3L1 alone.\n")
}

# ---- Tables ---------------------------------------------------------------
cmp <- data.frame(
  model = names(models),
  n_genes = vapply(models, length, integer(1)),
  auc = vapply(res, `[`, numeric(1), 1),
  ci_low = vapply(res, `[`, numeric(1), 2),
  ci_high = vapply(res, `[`, numeric(1), 3),
  row.names = NULL
)
print(cmp)
write.csv(cmp, res_path("02_model", "CHI3L1_comparison.csv"), row.names = FALSE)

dl_tab <- data.frame(
  comparison = "11-gene vs CHI3L1",
  delta_auc = cmp$auc[cmp$model == "11-gene signature"] -
              cmp$auc[cmp$model == "CHI3L1 only"],
  delta_ci_low = ci_delta[1],
  delta_ci_high = ci_delta[2],
  delong_P = dl$p.value
)
print(dl_tab)
write.csv(dl_tab, res_path("02_model", "CHI3L1_DeLong.csv"), row.names = FALSE)

compare("CHI3L1 AUC", cmp$auc[cmp$model == "CHI3L1 only"], REPORTED$auc_CHI3L1)
compare("DeLong P", dl$p.value, REPORTED$p_delong_CHI3L1)

saveRDS(list(comparison = cmp, delong = dl_tab, rocs = rocs,
             pred_11 = p11, pred_chi = pchi, group = obs),
        file.path(DATA_DIR, "derived", "CHI3L1_comparison.rds"))

cat("\nNext: 02_model/11_GSE130970_exploratory.R\n")
