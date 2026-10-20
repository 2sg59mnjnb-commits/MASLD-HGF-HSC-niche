# ============================================================================
# 11_GSE130970_exploratory.R
#
# Exploratory validation in GSE130970 (n = 78, METAVIR F0-F4, bariatric-surgery
# cohort). This cohort meets the inclusion criteria but differs clinically from
# the primary validation cohort, so it is labelled exploratory.
#
# Output
#   results/02_model/GSE130970_validation.csv
#   data/derived/GSE130970_validation.rds
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({
  library(GEOquery); library(pROC); library(org.Hs.eg.db); library(AnnotationDbi)
})

tpm_file <- file.path(DATA_DIR, "GSE130970_TPM.csv.gz")
if (!file.exists(tpm_file)) {
  url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130970/",
                "suppl/GSE130970_all_sample_salmon_tximport_TPM_entrez_gene_ID.csv.gz")
  cat("Downloading GSE130970 TPM matrix ...\n")
  dir.create(DATA_DIR, showWarnings = FALSE)
  download.file(url, tpm_file, mode = "wb", method = "libcurl")
}

# ---- 1. Expression: Entrez -> symbol --------------------------------------
raw <- read.csv(gzfile(tpm_file), row.names = 1, check.names = FALSE)
cat(sprintf("GSE130970 matrix: %d genes x %d samples\n", nrow(raw), ncol(raw)))

sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = rownames(raw),
                             column = "SYMBOL", keytype = "ENTREZID",
                             multiVals = "first")
ok <- !is.na(sym)
expr <- as.matrix(raw[ok, , drop = FALSE]); sym <- sym[ok]

# Collapse duplicate symbols by keeping the highest-mean row
ord <- order(-rowMeans(expr))
expr <- expr[ord, ]; sym <- sym[ord]
keep_first <- !duplicated(sym)
expr <- expr[keep_first, , drop = FALSE]
rownames(expr) <- sym[keep_first]
cat(sprintf("After symbol collapse: %d genes\n", nrow(expr)))

# ---- 2. Phenotype ---------------------------------------------------------
gse <- getGEO("GSE130970", GSEMatrix = TRUE, getGPL = FALSE)
pheno <- pData(gse[[1]])
fs_col <- grep("fibrosis stage", colnames(pheno), ignore.case = TRUE, value = TRUE)
if (length(fs_col) == 0) stop("Fibrosis stage column not found in GSE130970.")
stage <- as.numeric(gsub("[^0-9]", "", as.character(pheno[[fs_col[1]]])))
cat("Fibrosis stage distribution:\n"); print(table(stage, useNA = "ifany"))

# Align expression columns to phenotype rows by sample title
pheno$sample_id <- if (!is.null(pheno$title)) pheno$title else rownames(pheno)
common <- intersect(colnames(expr), pheno$sample_id)
cat(sprintf("Samples matched by title: %d / %d\n", length(common), ncol(expr)))
expr <- expr[, common, drop = FALSE]
stage <- stage[match(common, pheno$sample_id)]

grp <- ifelse(stage >= 2, 1, ifelse(stage <= 1, 0, NA))

# ---- 3. Predict -----------------------------------------------------------
ml <- load_derived("ML_inputs.rds")
X_tr <- ml$X[, SIG11, drop = FALSE]
Xtr_z <- scale(X_tr); Xtr_z[is.na(Xtr_z)] <- 0
fit <- glm(ml$y ~ ., data = as.data.frame(Xtr_z), family = binomial)

genes <- intersect(SIG11, rownames(expr))
cat(sprintf("Signature genes available: %d / 11\n", length(genes)))
if (length(genes) < 11) cat("  missing: ", paste(setdiff(SIG11, genes), collapse = ", "), "\n")

Xva <- t(expr[genes, , drop = FALSE])
Xva_z <- scale(Xva); Xva_z[is.na(Xva_z)] <- 0
p <- predict(fit, newdata = as.data.frame(Xva_z), type = "response")

keep <- !is.na(grp) & !is.na(p)
roc_obj <- roc(grp[keep], p[keep], quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))
ci_val  <- as.numeric(ci.auc(roc_obj))
brier   <- mean((p[keep] - grp[keep])^2)

cat(sprintf("\nGSE130970 AUC = %.3f (95%% CI %.3f-%.3f)\n",
            auc_val, ci_val[1], ci_val[3]))
cat(sprintf("n = %d, events = %d, Brier = %.3f\n",
            sum(keep), sum(grp[keep] == 1), brier))
compare("GSE130970 AUC", auc_val, REPORTED$auc_GSE130970)

# Calibration, reported to document the cross-cohort shift
cal <- glm(grp[keep] ~ qlogis(pmin(pmax(p[keep], 1e-10), 1 - 1e-10)),
           family = binomial)
cat(sprintf("Calibration slope = %.3f, intercept = %.3f\n",
            coef(cal)[2], coef(cal)[1]))

# ---- 4. Save --------------------------------------------------------------
perf <- data.frame(
  metric = c("auc", "ci_low", "ci_high", "brier", "cal_slope", "cal_intercept",
             "n", "events", "event_rate"),
  value  = c(auc_val, ci_val[1], ci_val[3], brier, coef(cal)[2], coef(cal)[1],
             sum(keep), sum(grp[keep] == 1), mean(grp[keep]))
)
write.csv(perf, res_path("02_model", "GSE130970_validation.csv"), row.names = FALSE)

saveRDS(list(roc = roc_obj, auc = auc_val, ci = ci_val, pred = p,
             group = grp, brier = brier),
        file.path(DATA_DIR, "derived", "GSE130970_validation.rds"))

cat("\nNext: 02_model/12_NAFLD_only_rederivation.R\n")
