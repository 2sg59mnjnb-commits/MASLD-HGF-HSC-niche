# ============================================================================
# 32_GSVA.R
#
# Pathway-level activity across fibrosis stages using GSVA with curated gene
# sets. MSigDB hallmark sets are used when msigdbr is available; otherwise a
# small curated collection is used, which keeps the analysis reproducible
# offline and avoids a very large dependency.
#
# Output
#   results/05_downstream/GSVA_pathway_trends.csv
#   data/derived/gsva.rds
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(GSVA); library(clinfun)})

d <- load_bulk()
expr <- to_symbol(d$log_expr)
stage <- d$fibro_stage

# ---- 1. Gene sets --------------------------------------------------------
curated <- list(
  IL6_JAK_STAT3 = c("IL6", "IL6R", "IL6ST", "JAK1", "JAK2", "JAK3", "TYK2",
                    "STAT1", "STAT3", "STAT5A", "STAT5B", "SOCS1", "SOCS3",
                    "PIM1", "MYC", "BCL2", "MCL1", "CCND1", "VEGFA", "HIF1A",
                    "IL1B", "TNF", "CXCL8", "CXCL9", "CXCL10", "CCL2",
                    "ICAM1", "VCAM1"),
  TNFA_NFKB = c("TNF", "TNFRSF1A", "TNFRSF1B", "NFKB1", "NFKB2", "RELA",
                "RELB", "NFKBIA", "NFKBIB", "NFKBIE", "IKBKB", "IKBKG", "CHUK",
                "TNFAIP3", "BIRC2", "BIRC3", "TRAF1", "TRAF2", "CXCL1",
                "CXCL2", "CXCL3", "CCL2", "CCL5", "ICAM1", "VCAM1", "SELE",
                "MMP9", "IL1B", "IL6"),
  TGF_BETA = c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "SMAD2",
               "SMAD3", "SMAD4", "SMAD7", "ACVR1", "ACVR2A", "ACVR2B",
               "BMPR1A", "BMPR2", "SERPINE1", "CTGF", "CCN2", "JUN", "FOS",
               "JUNB", "COL1A1", "COL1A2", "COL3A1", "FN1", "TIMP1"),
  EMT = c("CDH1", "CDH2", "VIM", "FN1", "SNAI1", "SNAI2", "TWIST1", "TWIST2",
          "ZEB1", "ZEB2", "TCF3", "TCF4", "COL1A1", "COL1A2", "COL3A1",
          "COL5A1", "COL5A2", "SPARC", "TAGLN", "ACTA2", "TGFB1", "TGFBR1",
          "MMP2", "MMP3", "MMP9", "TIMP1"),
  INFLAMMATORY = c("IL1A", "IL1B", "IL6", "IL6R", "TNF", "TNFRSF1A", "CXCL1",
                   "CXCL2", "CXCL3", "CXCL8", "CXCL9", "CXCL10", "CCL2",
                   "CCL3", "CCL4", "CCL5", "CCL7", "ICAM1", "VCAM1", "SELE",
                   "SELP", "NFKB1", "RELA", "NFKBIA", "TLR2", "TLR4", "MYD88",
                   "PTGS2", "NOS2", "STAT1"),
  IFN_GAMMA = c("IFNG", "IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1", "STAT2",
                "IRF1", "IRF2", "IRF9", "CXCL9", "CXCL10", "CXCL11", "GBP1",
                "GBP2", "GBP4", "GBP5", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA",
                "HLA-DRB1", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9")
)

if (requireNamespace("msigdbr", quietly = TRUE)) {
  cat("msigdbr available: hallmark sets will be added.\n")
  hm <- tryCatch({
    x <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    split(x$gene_symbol, x$gs_name)
  }, error = function(e) NULL)
  if (!is.null(hm)) {
    hm <- lapply(hm, function(g) intersect(g, rownames(expr)))
    hm <- hm[vapply(hm, length, integer(1)) >= 10]
    curated <- c(hm, curated)
    cat(sprintf("Total gene sets: %d\n", length(curated)))
  }
}

curated <- lapply(curated, function(g) intersect(g, rownames(expr)))
curated <- curated[vapply(curated, length, integer(1)) >= 5]
cat(sprintf("Gene sets used: %d\n", length(curated)))

# ---- 2. GSVA -------------------------------------------------------------
param <- gsvaParam(exprData = expr, geneSets = curated, kcdf = "Gaussian")
gsva_res <- gsva(param, verbose = FALSE)
cat(sprintf("GSVA matrix: %d sets x %d samples\n", nrow(gsva_res), ncol(gsva_res)))

# ---- 3. Trend per pathway -----------------------------------------------
trend <- do.call(rbind, lapply(rownames(gsva_res), function(pw) {
  x <- as.numeric(gsva_res[pw, ])
  jt <- jonckheere.test(x, as.numeric(factor(stage)), alternative = "increasing")
  sp <- cor.test(x, stage, method = "spearman")
  data.frame(pathway = pw,
             n_genes = length(curated[[pw]]),
             JT_P = jt$p.value,
             spearman_rho = round(as.numeric(sp$estimate), 3),
             spearman_P = sp$p.value,
             mean_F0 = round(mean(x[stage == 0]), 3),
             mean_F4 = round(mean(x[stage == 4]), 3))
}))
trend <- trend[order(trend$JT_P), ]
print(head(trend, 15), row.names = FALSE)
write.csv(trend, res_path("05_downstream", "GSVA_pathway_trends.csv"),
          row.names = FALSE)

cat(sprintf("\nPathways with JT P < 0.05: %d / %d\n",
            sum(trend$JT_P < 0.05), nrow(trend)))

saveRDS(list(gsva = gsva_res, trends = trend, gene_sets = curated),
        file.path(DATA_DIR, "derived", "gsva.rds"))

cat("\nNext: 05_downstream/33_deconvolution_NNLS.R\n")
