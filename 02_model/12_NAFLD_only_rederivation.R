# ============================================================================
# 12_NAFLD_only_rederivation.R
#
# Sensitivity analysis: rebuild the whole discovery pipeline after excluding the
# ten healthy controls, so that the comparison is NAFLD-vs-NAFLD (F0-F1 vs
# F2-F4) rather than disease-vs-control.
#
# This is the analysis behind the statement that the re-derived signature
# shares no genes with the reported 11-gene set. It is reported as evidence
# that gene membership is partly a control-versus-disease contrast, while
# discriminative performance is preserved.
#
# Output
#   results/02_model/NAFLD_only_rederivation.csv
#   results/02_model/NAFLD_only_module_correlation.csv
#   data/derived/NAFLD_only.rds
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({
  library(WGCNA); library(limma); library(caret); library(pROC)
  library(glmnet); library(randomForest); library(e1071); library(xgboost)
  library(org.Hs.eg.db); library(AnnotationDbi)
})
enableWGCNAThreads()
cor <- WGCNA::cor

aligned <- load_derived("GSE135251_aligned_data.rds")
pheno   <- aligned$pheno
counts  <- aligned$count_matrix

# ---- 1. Exclude controls --------------------------------------------------
ctrl <- pheno$gsm[grepl("control", pheno$disease, ignore.case = TRUE)]
cat(sprintf("Controls excluded: %d\n", length(ctrl)))
keep <- setdiff(colnames(counts), ctrl)
counts <- counts[, keep]
stage  <- pheno$fibrosis_stage[match(keep, pheno$gsm)]
grp    <- factor(ifelse(stage >= 2, "High", "Low"), levels = c("Low", "High"))
cat(sprintf("NAFLD-only cohort: %d samples (Low %d, High %d)\n",
            length(keep), sum(grp == "Low"), sum(grp == "High")))

# ---- 2. DEG ---------------------------------------------------------------
keep_g <- rowSums(counts >= 10) >= 3
log_expr <- log2(counts[keep_g, ] + 1)
fit <- eBayes(lmFit(log_expr, model.matrix(~ grp)))
deg <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
deg_sig <- deg[abs(deg$logFC) > 0.5 & deg$adj.P.Val < 0.05, ]
cat(sprintf("NAFLD-only DEGs: %d (full cohort: %d)\n",
            nrow(deg_sig), REPORTED$n_deg))

# ---- 3. WGCNA -------------------------------------------------------------
datExpr <- t(log_expr)
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

sft <- pickSoftThreshold(datExpr, powerVector = c(1:10, seq(12, 20, 2)),
                         networkType = "signed", verbose = 0)
r2 <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
power <- sft$fitIndices$Power[which(r2 > 0.85)[1]]
if (is.na(power)) power <- 9
cat(sprintf("Soft-threshold power: %d\n", power))

net <- blockwiseModules(datExpr, power = power,
                        networkType = "signed", TOMType = "signed",
                        minModuleSize = 30, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = FALSE, verbose = 0)

MEs <- net$MEs
r <- WGCNA::cor(MEs, as.numeric(grp == "High"), use = "p")
p <- corPvalueStudent(r, nrow(datExpr))
mod_tab <- data.frame(module = rownames(r), r = as.numeric(r[, 1]),
                      P = as.numeric(p[, 1]))
mod_tab <- mod_tab[order(-mod_tab$r), ]
print(head(mod_tab, 6))
write.csv(mod_tab, res_path("02_model", "NAFLD_only_module_correlation.csv"),
          row.names = FALSE)

top_pos <- mod_tab$module[which.max(mod_tab$r)]
top_neg <- mod_tab$module[which.min(mod_tab$r)]
cat(sprintf("Top positive module: %s (r = %.3f); top negative: %s (r = %.3f)\n",
            top_pos, mod_tab$r[mod_tab$module == top_pos],
            top_neg, mod_tab$r[mod_tab$module == top_neg]))
cat("(full-cohort reference: ME13 r = 0.484, ME14 r = -0.395)\n")

genes_of <- function(lbl) colnames(datExpr)[net$colors == as.numeric(gsub("ME", "", lbl))]
cand_ens <- unique(c(intersect(rownames(deg_sig), genes_of(top_pos)),
                     intersect(rownames(deg_sig), genes_of(top_neg))))

# ---- 4. Protein-coding filter --------------------------------------------
sy <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = cand_ens, column = "SYMBOL",
                            keytype = "ENSEMBL", multiVals = "first")
keep_s <- !is.na(sy)
ens <- cand_ens[keep_s]; sy <- as.character(sy[keep_s])
info <- AnnotationDbi::select(org.Hs.eg.db, keys = sy,
                              columns = c("SYMBOL", "GENETYPE"), keytype = "SYMBOL")
info <- info[!duplicated(info$SYMBOL), ]
coding <- info$SYMBOL[info$GENETYPE == "protein-coding"]
keep_c <- sy %in% coding
ens <- ens[keep_c]; sy <- sy[keep_c]
cat(sprintf("Candidates: %d (full cohort: %d)\n", length(ens), REPORTED$n_candidates))

# ---- 5. Four algorithms ---------------------------------------------------
X <- t(log_expr[ens, , drop = FALSE]); colnames(X) <- sy
sdv <- apply(X, 2, sd); X <- X[, is.finite(sdv) & sdv > 1e-8, drop = FALSE]
Xz <- scale(X); Xz[is.na(Xz)] <- 0
y  <- as.numeric(grp == "High")

set.seed(SEED)
cv <- cv.glmnet(Xz, y, family = "binomial", alpha = 1, nfolds = 10)
cf <- coef(cv, s = "lambda.min")
g_lasso <- setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")

set.seed(SEED)
rf <- randomForest(x = Xz, y = as.factor(y), ntree = 500, importance = TRUE)
im <- importance(rf)
g_rf <- rownames(im)[order(im[, "MeanDecreaseGini"], decreasing = TRUE)][
  seq_len(min(30, ncol(Xz)))]

svm_rfe <- function(X, y, n_target = 20) {
  X <- as.data.frame(X); f <- colnames(X); yf <- as.factor(y)
  while (length(f) > n_target) {
    m <- svm(x = X[, f, drop = FALSE], y = yf, kernel = "linear")
    w <- abs(as.numeric(t(m$coefs) %*% m$SV))
    if (length(w) != length(f)) break
    names(w) <- f; f <- f[-which.min(w)]
  }
  f
}
set.seed(SEED); g_svm <- svm_rfe(Xz, y, 20)

set.seed(SEED)
dt <- xgb.DMatrix(data = Xz, label = y)
xm <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc",
                              max_depth = 3, eta = 0.1),
                data = dt, nrounds = 100, verbose = 0)
imx <- xgb.importance(model = xm)
g_xgb <- imx$Feature[seq_len(min(30, nrow(imx)))]

core <- Reduce(intersect, list(g_lasso, g_rf, g_svm, g_xgb))
fallback <- length(core) < 3
if (fallback) {
  all_g <- unique(c(g_lasso, g_rf, g_svm, g_xgb))
  rk <- setNames(im[, "MeanDecreaseGini"], rownames(im))
  core <- all_g[order(-rk[all_g], na.last = TRUE)][seq_len(min(15, length(all_g)))]
}
cat(sprintf("\nNAFLD-only signature: %d genes%s\n", length(core),
            ifelse(fallback, " (fallback to union)", "")))
print(core)

overlap <- intersect(core, SIG11)
cat(sprintf("\nOverlap with the reported 11-gene signature: %d / 11\n", length(overlap)))
if (length(overlap)) print(overlap)
cat("Genes unique to the NAFLD-only set: ",
    paste(setdiff(core, SIG11), collapse = ", "), "\n")

# ---- 6. CV performance of the NAFLD-only signature -----------------------
set.seed(SEED)
folds <- createFolds(y, k = 10, list = TRUE)
pred <- rep(NA_real_, length(y))
for (i in seq_along(folds)) {
  te <- folds[[i]]; tr <- setdiff(seq_along(y), te)
  Xtr <- scale(X[tr, core, drop = FALSE]); Xtr[is.na(Xtr)] <- 0
  mu <- attr(Xtr, "scaled:center"); sdd <- attr(Xtr, "scaled:scale")
  Xte <- scale(X[te, core, drop = FALSE], center = mu, scale = sdd)
  Xte[is.na(Xte)] <- 0
  df <- as.data.frame(Xtr); df$y <- y[tr]
  m <- glm(y ~ ., data = df, family = binomial)
  pred[te] <- predict(m, newdata = as.data.frame(Xte), type = "response")
}
auc_nafld <- as.numeric(auc(roc(y, pred, quiet = TRUE)))
cat(sprintf("\nNAFLD-only signature CV AUC within the NAFLD-only cohort = %.3f\n",
            auc_nafld))
cat("(the reported 11-gene signature scores 0.821 when locked and applied to the same samples)\n")

# ---- 7. Save --------------------------------------------------------------
out <- data.frame(
  metric = c("n_samples", "n_low", "n_high", "n_deg", "power",
             "r_top_positive", "r_top_negative", "n_candidates",
             "n_signature", "fallback", "overlap_with_SIG11", "cv_auc"),
  value = c(length(y), sum(y == 0), sum(y == 1), nrow(deg_sig), power,
            mod_tab$r[mod_tab$module == top_pos],
            mod_tab$r[mod_tab$module == top_neg],
            length(ens), length(core), as.integer(fallback),
            length(overlap), auc_nafld)
)
write.csv(out, res_path("02_model", "NAFLD_only_rederivation.csv"), row.names = FALSE)

saveRDS(list(signature = core, overlap = overlap, cv_auc = auc_nafld,
             pred = pred, y = y, deg = deg_sig, modules = mod_tab,
             fallback = fallback),
        file.path(DATA_DIR, "derived", "NAFLD_only.rds"))

cat("\nNext: 03_spatial/13_load_CosMx_and_panel.R\n")
