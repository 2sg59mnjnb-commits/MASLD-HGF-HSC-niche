# ============================================================================
# 31_stage_trends.R
#
# Do the signature and the IL6-family axis track fibrosis stage across the full
# F0-F4 range, rather than only across the binary F0-F1 / F2-F4 split used for
# model training? The spatial data cover F2 only, so this bulk analysis is what
# addresses the "does it generalise across stage" question.
#
# Tests per feature
#   Jonckheere-Terpstra   ordered monotonic trend across F0..F4
#   Spearman rho          continuous association with stage
#   Kruskal-Wallis        any difference between stages
#
# Output
#   results/05_downstream/stage_trends_signature.csv
#   results/05_downstream/stage_trends_genes.csv
#   results/05_downstream/stage_trends_pathways.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(clinfun)})

d <- load_bulk()
expr <- to_symbol(d$log_expr)
stage <- d$fibro_stage
cat(sprintf("Cohort: %d genes x %d samples\n", nrow(expr), ncol(expr)))
cat("Stage distribution:\n"); print(table(stage))

trend_test <- function(x) {
  x <- as.numeric(x)
  jt <- jonckheere.test(x, as.numeric(factor(stage)), alternative = "increasing")
  sp <- cor.test(x, stage, method = "spearman")
  kw <- kruskal.test(x ~ factor(stage))
  list(JT_P = jt$p.value, rho = as.numeric(sp$estimate), sp_P = sp$p.value,
       KW_P = kw$p.value, means = tapply(x, stage, mean))
}

# ---- 1. Signature score --------------------------------------------------
sig <- sig_score(expr, SIG11)
res <- trend_test(sig)

cat(sprintf("\n11-gene signature across F0-F4:\n"))
cat(sprintf("  Jonckheere-Terpstra P = %.3e\n", res$JT_P))
cat(sprintf("  Spearman rho = %.3f (P = %.3e)\n", res$rho, res$sp_P))
cat("  Mean score by stage:\n"); print(round(res$means, 3))

sig_tab <- data.frame(
  feature = "11-gene signature",
  JT_P = res$JT_P, spearman_rho = round(res$rho, 3), spearman_P = res$sp_P,
  KW_P = res$KW_P,
  mean_F0 = round(res$means["0"], 3), mean_F1 = round(res$means["1"], 3),
  mean_F2 = round(res$means["2"], 3), mean_F3 = round(res$means["3"], 3),
  mean_F4 = round(res$means["4"], 3)
)
write.csv(sig_tab, res_path("05_downstream", "stage_trends_signature.csv"),
          row.names = FALSE)

# ---- 2. Individual signature genes --------------------------------------
genes_present <- intersect(SIG11, rownames(expr))
cat(sprintf("\nSignature genes available: %d / 11\n", length(genes_present)))

gene_tab <- do.call(rbind, lapply(genes_present, function(g) {
  r <- trend_test(expr[g, ])
  data.frame(gene = g, JT_P = r$JT_P, spearman_rho = round(r$rho, 3),
             spearman_P = r$sp_P,
             mean_F0 = round(r$means["0"], 3), mean_F4 = round(r$means["4"], 3))
}))
gene_tab <- gene_tab[order(-gene_tab$spearman_rho), ]
print(gene_tab, row.names = FALSE)
write.csv(gene_tab, res_path("05_downstream", "stage_trends_genes.csv"),
          row.names = FALSE)

# Is CHI3L1 present and rising?
if ("CHI3L1" %in% rownames(expr)) {
  r <- trend_test(expr["CHI3L1", ])
  cat(sprintf("\nCHI3L1: rho = %.3f (P = %.2e), JT P = %.2e\n",
              r$rho, r$sp_P, r$JT_P))
}

# ---- 3. IL6-family axis --------------------------------------------------
il6_all <- c(IL6_LIGANDS, IL6_RECEPTORS, "HGF", "MET")
il6_in <- intersect(il6_all, rownames(expr))
cat(sprintf("\nIL6-family / HGF axis available: %d / %d\n",
            length(il6_in), length(il6_all)))

il6_tab <- do.call(rbind, lapply(il6_in, function(g) {
  r <- trend_test(expr[g, ])
  data.frame(gene = g,
             type = ifelse(g %in% IL6_LIGANDS, "ligand",
                           ifelse(g %in% IL6_RECEPTORS, "receptor", "axis")),
             JT_P = r$JT_P, spearman_rho = round(r$rho, 3),
             spearman_P = r$sp_P,
             mean_F0 = round(r$means["0"], 3), mean_F4 = round(r$means["4"], 3))
}))
il6_tab <- il6_tab[order(-il6_tab$spearman_rho), ]
print(il6_tab, row.names = FALSE)
write.csv(il6_tab, res_path("05_downstream", "stage_trends_IL6_axis.csv"),
          row.names = FALSE)

# ---- 4. Pathway scores ---------------------------------------------------
path_tab <- do.call(rbind, lapply(
  list(JAK_STAT = JAK_STAT_GENES, NFKB = NFKB_GENES), function(gs) {
    g <- intersect(gs, rownames(expr))
    if (length(g) < 3) return(NULL)
    s <- colMeans(zscore_rows(expr[g, , drop = FALSE]))
    r <- trend_test(s)
    data.frame(pathway = paste(g, collapse = "/"),
               n_genes = length(g),
               JT_P = r$JT_P, spearman_rho = round(r$rho, 3),
               spearman_P = r$sp_P,
               mean_F0 = round(r$means["0"], 3), mean_F4 = round(r$means["4"], 3))
  }))
print(path_tab, row.names = FALSE)
write.csv(path_tab, res_path("05_downstream", "stage_trends_pathways.csv"),
          row.names = FALSE)

saveRDS(list(signature = sig_tab, genes = gene_tab, il6 = il6_tab,
             pathways = path_tab, score = sig, stage = stage),
        file.path(DATA_DIR, "derived", "stage_trends.rds"))

cat("\nNext: 05_downstream/32_GSVA.R\n")
