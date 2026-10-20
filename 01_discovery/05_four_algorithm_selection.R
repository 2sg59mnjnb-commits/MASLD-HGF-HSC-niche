# ============================================================================
# 05_four_algorithm_selection.R
#
# Derives the 11-gene tissue signature. Four feature-selection algorithms are
# applied to the 38-gene candidate set and their selected gene sets are
# intersected:
#
#   LASSO       cv.glmnet, alpha = 1, 10-fold, lambda.min
#   Random for. randomForest, 500 trees, top 30 by MeanDecreaseGini
#   SVM-RFE     recursive feature elimination down to 20 genes
#   XGBoost     max_depth 3, eta 0.1, 100 rounds, top 30 by gain
#
# The intersection of the four sets is the reported 11-gene signature.
#
# Output
#   results/01_discovery/feature_selection_genes.csv
#   data/derived/ML_inputs.rds   (X, y, per-algorithm gene lists, signature)
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({
  library(glmnet)
  library(randomForest)
  library(e1071)
  library(xgboost)
})

d <- load_bulk()
cand <- read.csv(res_path("01_discovery", "candidate_genes_filtered.csv"),
                 stringsAsFactors = FALSE)

# Expression matrix restricted to the candidates: samples x genes.
X <- t(d$log_expr[cand$Ensembl, , drop = FALSE])
colnames(X) <- cand$Symbol
y <- as.numeric(d$group == "High")

cat(sprintf("X: %d samples x %d candidate genes\n", nrow(X), ncol(X)))
cat(sprintf("Events (High fibrosis): %d\n", sum(y)))

# ---- LASSO ----------------------------------------------------------------
# X is centred/scaled internally by glmnet; do not pre-scale here, otherwise
# the coefficient scale differs from the reported model.
set.seed(SEED)
cv_lasso <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 10)
lasso_coef <- coef(cv_lasso, s = "lambda.min")
g_lasso <- setdiff(rownames(lasso_coef)[lasso_coef[, 1] != 0], "(Intercept)")
cat(sprintf("\nLASSO: %d genes\n", length(g_lasso)))
print(g_lasso)

# ---- Random forest --------------------------------------------------------
set.seed(SEED)
rf_fit <- randomForest(x = X, y = as.factor(y), ntree = 500, importance = TRUE)
rf_imp <- importance(rf_fit)
g_rf <- rownames(rf_imp)[order(rf_imp[, "MeanDecreaseGini"], decreasing = TRUE)]
g_rf <- g_rf[seq_len(min(30, ncol(X)))]
cat(sprintf("\nRandom forest: %d genes (top 30 by MeanDecreaseGini)\n", length(g_rf)))

# ---- SVM-RFE --------------------------------------------------------------
# Recursive elimination using an SVM fit. With a radial kernel there are no
# per-feature weights, so the elimination step falls back to dropping the
# feature with the smallest variance -- this is the implementation that was
# used to derive the reported signature and is kept for reproducibility.
# (The nested cross-validation script uses a linear-kernel variant, which does
# provide weights; see 02_model/07_fully_nested_cv.R.)
svm_rfe <- function(X, y, n_target = 20, kernel = "radial", verbose = FALSE) {
  X <- as.data.frame(X)
  feats <- colnames(X)
  yf <- as.factor(y)
  while (length(feats) > n_target) {
    fit <- svm(x = X[, feats, drop = FALSE], y = yf,
               kernel = kernel, scale = TRUE)
    w <- abs(as.numeric(t(fit$coefs) %*% fit$SV))
    if (length(w) != length(feats)) {
      v <- apply(X[, feats, drop = FALSE], 2, var)
      drop <- which.min(v)
    } else {
      names(w) <- feats
      drop <- which.min(w)
    }
    feats <- feats[-drop]
    if (verbose) cat(sprintf("  remaining: %d\n", length(feats)))
  }
  feats
}

set.seed(SEED)
g_svm <- svm_rfe(X, y, n_target = 20, kernel = "radial")
cat(sprintf("\nSVM-RFE: %d genes\n", length(g_svm)))

# ---- XGBoost --------------------------------------------------------------
set.seed(SEED)
dtrain <- xgb.DMatrix(data = X, label = y)
xgb_fit <- xgb.train(
  params = list(objective = "binary:logistic", eval_metric = "auc",
                max_depth = 3, eta = 0.1),
  data = dtrain, nrounds = 100, verbose = 0
)
xgb_imp <- xgb.importance(model = xgb_fit)
g_xgb <- xgb_imp$Feature[seq_len(min(30, nrow(xgb_imp)))]
cat(sprintf("\nXGBoost: %d genes (top 30 by gain)\n", length(g_xgb)))

# ---- Intersection ---------------------------------------------------------
signature <- Reduce(intersect, list(g_lasso, g_rf, g_svm, g_xgb))
cat(sprintf("\nIntersection of the four algorithms: %d genes\n", length(signature)))
print(signature)

cat(sprintf("\nOverlap with the signature reported in the manuscript: %d / 11\n",
            length(intersect(signature, SIG11))))
if (length(setdiff(SIG11, signature))) {
  cat("Reported genes not recovered here: ",
      paste(setdiff(SIG11, signature), collapse = ", "), "\n")
}
if (length(setdiff(signature, SIG11))) {
  cat("Additional genes selected here:     ",
      paste(setdiff(signature, SIG11), collapse = ", "), "\n")
}

# ---- Save -----------------------------------------------------------------
gene_tab <- data.frame(
  gene = union(union(g_lasso, g_rf), union(g_svm, g_xgb)),
  in_LASSO = as.integer(union(union(g_lasso, g_rf), union(g_svm, g_xgb)) %in% g_lasso),
  in_RF    = as.integer(union(union(g_lasso, g_rf), union(g_svm, g_xgb)) %in% g_rf),
  in_SVM   = as.integer(union(union(g_lasso, g_rf), union(g_svm, g_xgb)) %in% g_svm),
  in_XGB   = as.integer(union(union(g_lasso, g_rf), union(g_svm, g_xgb)) %in% g_xgb)
)
gene_tab$n_algorithms <- rowSums(gene_tab[, -1])
gene_tab <- gene_tab[order(-gene_tab$n_algorithms, gene_tab$gene), ]
write.csv(gene_tab, res_path("01_discovery", "feature_selection_genes.csv"),
          row.names = FALSE)

saveRDS(list(X = X, y = y,
             lasso = g_lasso, rf = g_rf, svm = g_svm, xgb = g_xgb,
             signature = signature,
             cv_lasso = cv_lasso, rf_fit = rf_fit, xgb_fit = xgb_fit),
        file.path(DATA_DIR, "derived", "ML_inputs.rds"))

cat("\nNext: 02_model/06_fixed_signature_cv.R\n")
