# ============================================================================
# 07_fully_nested_cv.R
#
# Fully nested cross-validation: the entire discovery pipeline is repeated
# inside every training fold, so that feature selection cannot see the held-out
# samples. Per fold:
#
#   limma DEG  ->  WGCNA (top 5,000 variable genes)  ->  module/DEG candidates
#   ->  protein-coding filter  ->  LASSO / RF / SVM-RFE / XGBoost
#   ->  four-way intersection  ->  logistic model  ->  predict held-out fold
#
# Soft-threshold power: the default (`power_mode = "full"`) reuses the power
# selected on the full cohort, which is what was run for the manuscript. This
# is a mild optimism because the power is chosen once on all samples. Set
# `power_mode = "in_fold"` to re-select it inside every fold (slower, stricter).
#
# Run time: 1-3 h (10 folds x WGCNA).
#
# Output
#   results/02_model/nested_cv_summary.csv
#   results/02_model/nested_cv_stability.csv
#   data/derived/nested_cv.rds
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({
  library(WGCNA); library(limma); library(caret); library(pROC)
  library(glmnet); library(randomForest); library(e1071); library(xgboost)
  library(org.Hs.eg.db)
})
enableWGCNAThreads()
cor <- WGCNA::cor

POWER_MODE <- "full"   # "full" (as published) or "in_fold"

d <- load_bulk()
log_expr  <- d$log_expr
y         <- as.numeric(d$group == "High")
power_ref <- load_derived("WGCNA_results.rds")$power
cat(sprintf("Expression: %d genes x %d samples; reference power = %d\n",
            nrow(log_expr), ncol(log_expr), power_ref))

# Linear-kernel SVM-RFE: unlike the radial kernel it yields per-feature weights,
# so elimination is driven by the model rather than by variance.
svm_rfe_linear <- function(X, y, n_target = 20) {
  X <- as.data.frame(X); feats <- colnames(X); yf <- as.factor(y)
  while (length(feats) > n_target) {
    fit <- svm(x = X[, feats, drop = FALSE], y = yf, kernel = "linear")
    w <- abs(as.numeric(t(fit$coefs) %*% fit$SV))
    if (length(w) != length(feats)) break
    names(w) <- feats
    feats <- feats[-which.min(w)]
  }
  feats
}

# Run the discovery pipeline on one training subset; return selected genes.
discover <- function(log_tr, y_tr, power_mode, power_fixed) {
  design <- model.matrix(~ y_tr)
  fit <- eBayes(lmFit(log_tr, design))
  deg <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  deg_sig <- deg[abs(deg$logFC) > 0.5 & deg$adj.P.Val < 0.05, ]

  # WGCNA on the most variable genes to keep the fold tractable
  v <- apply(log_tr, 1, var)
  keep <- names(sort(v, decreasing = TRUE))[seq_len(min(5000, length(v)))]
  datExpr <- t(log_tr[keep, , drop = FALSE])
  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

  power <- power_fixed
  if (power_mode == "in_fold") {
    sft <- pickSoftThreshold(datExpr, powerVector = c(1:10, 12, 14, 16),
                             networkType = "signed", verbose = 0)
    r2 <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
    cand <- sft$fitIndices$Power[which(r2 > 0.85)[1]]
    if (!is.na(cand)) power <- cand
  }

  net <- blockwiseModules(datExpr, power = power,
                          networkType = "signed", TOMType = "signed",
                          minModuleSize = 30, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = FALSE, verbose = 0)
  MEs <- net$MEs
  r <- WGCNA::cor(MEs, y_tr)
  top <- colnames(MEs)[which.max(abs(r))]
  mod_genes <- colnames(datExpr)[net$colors == as.numeric(gsub("ME", "", top))]

  cand_ens <- intersect(rownames(deg_sig), mod_genes)
  if (length(cand_ens) < 20) {
    cand_ens <- unique(c(cand_ens, rownames(deg_sig)[seq_len(min(50, nrow(deg_sig)))]))
    cand_ens <- cand_ens[cand_ens %in% rownames(log_tr)]
  }

  # Ensembl -> symbol, then restrict to protein-coding genes
  sym_map <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = cand_ens,
                                   column = "SYMBOL", keytype = "ENSEMBL",
                                   multiVals = "first")
  keep <- !is.na(sym_map)
  ens  <- cand_ens[keep]
  sy   <- as.character(sym_map[keep])

  info <- AnnotationDbi::select(org.Hs.eg.db, keys = sy,
                                columns = c("SYMBOL", "GENETYPE"),
                                keytype = "SYMBOL")
  info <- info[!duplicated(info$SYMBOL), ]
  coding <- info$SYMBOL[info$GENETYPE == "protein-coding"]

  keep2 <- sy %in% coding
  ens <- ens[keep2]; sy <- sy[keep2]

  X <- t(log_tr[ens, , drop = FALSE])
  colnames(X) <- sy
  sdv <- apply(X, 2, sd)
  X <- X[, !is.na(sdv) & sdv > 1e-8, drop = FALSE]
  Xz <- scale(X); Xz[is.na(Xz)] <- 0

  g_lasso <- local({
    set.seed(SEED); cv <- cv.glmnet(Xz, y_tr, family = "binomial",
                                    alpha = 1, nfolds = 5)
    cf <- coef(cv, s = "lambda.min")
    setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")
  })
  g_rf <- local({
    set.seed(SEED); rf <- randomForest(x = Xz, y = as.factor(y_tr),
                                       ntree = 500, importance = TRUE)
    im <- importance(rf)
    rownames(im)[order(im[, "MeanDecreaseGini"], decreasing = TRUE)][
      seq_len(min(30, ncol(Xz)))]
  })
  g_svm <- local({
    set.seed(SEED)
    tryCatch(svm_rfe_linear(Xz, y_tr, 20), error = function(e) g_rf)
  })
  g_xgb <- local({
    set.seed(SEED)
    dt <- xgb.DMatrix(data = Xz, label = y_tr)
    m <- xgb.train(params = list(objective = "binary:logistic",
                                 eval_metric = "auc", max_depth = 3, eta = 0.1),
                   data = dt, nrounds = 100, verbose = 0)
    im <- xgb.importance(model = m)
    im$Feature[seq_len(min(30, nrow(im)))]
  })

  core <- Reduce(intersect, list(g_lasso, g_rf, g_svm, g_xgb))
  fallback <- length(core) < 3
  if (fallback) {
    all_g <- unique(c(g_lasso, g_rf, g_svm, g_xgb))
    core <- all_g[seq_len(min(15, length(all_g)))]
  }

  list(genes = core, candidates = colnames(X), fallback = fallback,
       power = power, n_deg = nrow(deg_sig))
}

# ---- Fold loop -------------------------------------------------------------
set.seed(SEED)
folds <- createFolds(y, k = 10, list = TRUE)
cv_pred <- rep(NA_real_, length(y))
fold_tab <- list()

for (i in seq_along(folds)) {
  cat(sprintf("\n=== Fold %d / 10 ===\n", i))
  test  <- folds[[i]]
  train <- setdiff(seq_along(y), test)

  log_tr <- log_expr[, train, drop = FALSE]
  y_tr   <- y[train]

  disc <- discover(log_tr, y_tr, POWER_MODE, power_ref)
  cat(sprintf("  in-fold DEGs: %d | signature: %d genes%s\n",
              disc$n_deg, length(disc$genes),
              ifelse(disc$fallback, " (fallback)", "")))

  # Gene-wise scaling parameters learned on the training fold only
  Xtr_raw <- t(log_tr[disc$candidates, , drop = FALSE])
  mu  <- apply(Xtr_raw, 2, mean)
  sdv <- apply(Xtr_raw, 2, sd)
  sdv[!is.finite(sdv) | sdv < 1e-8] <- 1

  Xtr <- scale(Xtr_raw, center = mu, scale = sdv)
  Xtr[is.na(Xtr)] <- 0
  df_tr <- as.data.frame(Xtr[, disc$genes, drop = FALSE]); df_tr$y <- y_tr
  fit <- glm(y ~ ., data = df_tr, family = binomial)

  Xte_raw <- t(log_expr[disc$candidates, test, drop = FALSE])
  Xte <- scale(Xte_raw, center = mu, scale = sdv)
  Xte[is.na(Xte)] <- 0
  Xte <- as.data.frame(Xte)
  for (g in disc$genes) if (!g %in% colnames(Xte)) Xte[[g]] <- 0

  cv_pred[test] <- predict(fit, newdata = Xte[, disc$genes, drop = FALSE],
                           type = "response")

  fold_tab[[i]] <- data.frame(fold = i, n_genes = length(disc$genes),
                              fallback = disc$fallback, power = disc$power,
                              genes = paste(disc$genes, collapse = ";"))
}

roc_nested <- roc(y, cv_pred, quiet = TRUE)
auc_nested <- as.numeric(auc(roc_nested))
ci_nested  <- as.numeric(ci.auc(roc_nested))

cat(sprintf("\nFully nested CV AUC = %.3f (95%% CI %.3f-%.3f)\n",
            auc_nested, ci_nested[1], ci_nested[3]))
compare("nested CV AUC", auc_nested, REPORTED$auc_cv_nested)

fold_df <- do.call(rbind, fold_tab)
cat(sprintf("Signature size per fold: median %d, range %d-%d\n",
            median(fold_df$n_genes), min(fold_df$n_genes), max(fold_df$n_genes)))
cat(sprintf("Folds that fell back to the union: %d / 10\n", sum(fold_df$fallback)))

# ---- Stability table -------------------------------------------------------
all_genes <- sort(unique(unlist(strsplit(fold_df$genes, ";"))))
stab <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
for (i in seq_len(nrow(fold_df))) {
  stab[[paste0("fold_", fold_df$fold[i])]] <-
    as.integer(all_genes %in% strsplit(fold_df$genes[i], ";")[[1]])
}
stab$n_folds <- rowSums(stab[, grep("^fold_", colnames(stab))])
stab <- stab[order(-stab$n_folds, stab$gene), ]

cat("\nGenes selected in >= 50% of folds:\n")
print(head(stab[stab$n_folds >= 5, c("gene", "n_folds")], 20))
cat(sprintf("\nReported signature recovered in each fold:\n"))
for (g in SIG11) {
  cat(sprintf("  %-9s %d / 10 folds\n", g, sum(sapply(strsplit(fold_df$genes, ";"),
                                                     function(s) g %in% s))))
}

write.csv(fold_df, res_path("02_model", "nested_cv_folds.csv"), row.names = FALSE)
write.csv(stab,    res_path("02_model", "nested_cv_stability.csv"), row.names = FALSE)
write.csv(data.frame(metric = c("nested_auc", "ci_low", "ci_high",
                                "median_signature_size", "n_fallback_folds"),
                     value = c(auc_nested, ci_nested[1], ci_nested[3],
                               median(fold_df$n_genes), sum(fold_df$fallback))),
          res_path("02_model", "nested_cv_summary.csv"), row.names = FALSE)

saveRDS(list(cv_pred = cv_pred, y = y, roc = roc_nested, auc = auc_nested,
             ci = ci_nested, folds = fold_df, stability = stab,
             power_mode = POWER_MODE),
        file.path(DATA_DIR, "derived", "nested_cv.rds"))

cat("\nNext: 02_model/08_external_validation_GSE49541.R\n")
