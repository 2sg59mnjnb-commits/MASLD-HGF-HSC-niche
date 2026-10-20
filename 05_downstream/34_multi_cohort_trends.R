# ============================================================================
# 34_multi_cohort_trends.R
#
# Re-derives the signature score in two additional cohorts that have full
# fibrosis staging but were not used for validation of discrimination:
#   GSE89632  Illumina HumanHT-12 v4, n ~ 55 with stage, NAFLD activity score
#   GSE84044  Affymetrix HG-U133 Plus 2.0, n = 124, Scheuer score
#
# The question is whether the signature tracks stage continuously, not whether
# it discriminates a binary endpoint. Both cohorts are therefore analysed with a
# trend test rather than an AUC.
#
# Output
#   results/05_downstream/multi_cohort_trends.csv
#   results/05_downstream/multi_cohort_gene_trends.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({
  library(GEOquery); library(clinfun); library(AnnotationDbi); library(org.Hs.eg.db)
})

dir.create(file.path(DATA_DIR, "public_cohorts"), showWarnings = FALSE,
           recursive = TRUE)

fetch <- function(acc) {
  f <- file.path(DATA_DIR, "public_cohorts", paste0(acc, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  dir.create(file.path(DATA_DIR, "public_cohorts", acc),
             showWarnings = FALSE, recursive = TRUE)
  gse <- getGEO(acc, GSEMatrix = TRUE, getGPL = FALSE,
                destdir = file.path(DATA_DIR, "public_cohorts", acc))
  out <- list(expr = exprs(gse[[1]]), pheno = pData(gse[[1]]))
  saveRDS(out, f)
  out
}

# Collapse probes to symbols by mean
collapse_probes <- function(expr, symbol) {
  ok <- !is.na(symbol)
  expr <- expr[ok, , drop = FALSE]; symbol <- symbol[ok]
  df <- as.data.frame(expr, check.names = FALSE)
  df$symbol <- symbol
  out <- df %>%
    dplyr::group_by(symbol) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean),
                     .groups = "drop") %>%
    as.data.frame()
  rownames(out) <- out$symbol
  out$symbol <- NULL
  as.matrix(out)
}

trend_row <- function(score, stage, cohort) {
  jt <- jonckheere.test(score, as.numeric(factor(stage)), alternative = "increasing")
  sp <- cor.test(score, stage, method = "spearman")
  data.frame(cohort = cohort, n = length(score),
             n_stages = length(unique(stage)),
             JT_P = jt$p.value,
             spearman_rho = round(as.numeric(sp$estimate), 3),
             spearman_P = sp$p.value)
}

results <- list()
gene_results <- list()

# ---- GSE89632 -------------------------------------------------------------
cat("=== GSE89632 ===\n")
g <- fetch("GSE89632")
suppressPackageStartupMessages(library(illuminaHumanv4.db))
sym <- AnnotationDbi::mapIds(illuminaHumanv4.db, keys = rownames(g$expr),
                             column = "SYMBOL", keytype = "PROBEID",
                             multiVals = "first")
expr <- collapse_probes(g$expr, sym)
cat(sprintf("Expression: %d genes x %d samples\n", nrow(expr), ncol(expr)))

stage <- suppressWarnings(as.numeric(gsub("[^0-9]", "",
          as.character(g$pheno[["fibrosis (stage):ch1"]]))))
keep <- !is.na(stage) & stage %in% 0:4
cat(sprintf("Samples with stage: %d\n", sum(keep)))

if (sum(keep) >= 10) {
  expr_k <- expr[, keep, drop = FALSE]
  score <- sig_score(expr_k, SIG11)
  results[["GSE89632"]] <- trend_row(score, stage[keep], "GSE89632")

  for (gene in intersect(SIG11, rownames(expr_k))) {
    sp <- cor.test(expr_k[gene, ], stage[keep], method = "spearman")
    gene_results[["GSE89632"]] <- rbind(
      gene_results[["GSE89632"]],
      data.frame(cohort = "GSE89632", gene = gene,
                 rho = round(as.numeric(sp$estimate), 3), P = sp$p.value))
  }
}

# ---- GSE84044 -------------------------------------------------------------
cat("\n=== GSE84044 ===\n")
g <- fetch("GSE84044")
suppressPackageStartupMessages(library(hgu133plus2.db))
sym <- AnnotationDbi::mapIds(hgu133plus2.db, keys = rownames(g$expr),
                             column = "SYMBOL", keytype = "PROBEID",
                             multiVals = "first")
expr <- collapse_probes(g$expr, sym)
cat(sprintf("Expression: %d genes x %d samples\n", nrow(expr), ncol(expr)))

stage <- suppressWarnings(as.numeric(
  as.character(g$pheno[["scheuer score s:ch1"]])))
keep <- !is.na(stage)
cat(sprintf("Samples with Scheuer score: %d\n", sum(keep)))

if (sum(keep) >= 10) {
  expr_k <- expr[, keep, drop = FALSE]
  score <- sig_score(expr_k, SIG11)
  results[["GSE84044"]] <- trend_row(score, stage[keep], "GSE84044")

  for (gene in intersect(SIG11, rownames(expr_k))) {
    sp <- cor.test(expr_k[gene, ], stage[keep], method = "spearman")
    gene_results[["GSE84044"]] <- rbind(
      gene_results[["GSE84044"]],
      data.frame(cohort = "GSE84044", gene = gene,
                 rho = round(as.numeric(sp$estimate), 3), P = sp$p.value))
  }
}

# ---- Summarise ------------------------------------------------------------
if (length(results)) {
  tab <- do.call(rbind, results)
  print(tab, row.names = FALSE)
  write.csv(tab, res_path("05_downstream", "multi_cohort_trends.csv"),
            row.names = FALSE)
}
if (length(gene_results)) {
  gt <- do.call(rbind, gene_results)
  write.csv(gt, res_path("05_downstream", "multi_cohort_gene_trends.csv"),
            row.names = FALSE)
}

cat("\nInterpretation: both cohorts test the continuous association with stage,\n")
cat("not discrimination of a binary endpoint. A positive rho is consistent with\n")
cat("the signature tracking fibrosis severity; it does not establish clinical\n")
cat("utility.\n")

saveRDS(list(summary = if (length(results)) do.call(rbind, results) else NULL,
             genes = if (length(gene_results)) do.call(rbind, gene_results) else NULL),
        file.path(DATA_DIR, "derived", "multi_cohort_trends.rds"))

cat("\nNext: 05_downstream/35_incremental_value.R\n")
