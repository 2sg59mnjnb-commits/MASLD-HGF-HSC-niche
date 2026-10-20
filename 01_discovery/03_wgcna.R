# ============================================================================
# 03_wgcna.R
#
# Weighted gene co-expression network analysis of the training cohort.
# Identifies co-expression modules and correlates module eigengenes with the
# binary fibrosis group (High F2-F4 vs Low F0-F1).
#
# Parameters follow the analysis reported in the manuscript:
#   signed network, signed TOM, minModuleSize = 30, mergeCutHeight = 0.25,
#   soft-threshold power chosen as the smallest power with scale-free R^2 > 0.85
#
# Output
#   results/01_discovery/WGCNA_module_trait.csv
#   results/01_discovery/WGCNA_soft_threshold.csv
#   data/derived/WGCNA_results.rds   (net, MEs, correlations, members)
#
# Run time: 20-40 min.
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(WGCNA)})
enableWGCNAThreads()
# WGCNA masks stats::cor; be explicit from here on.
cor <- WGCNA::cor

d <- load_bulk()
log_expr <- d$log_expr
group    <- d$group

cat(sprintf("Input: %d genes x %d samples\n", nrow(log_expr), ncol(log_expr)))

# ---- 1. Transpose and drop bad samples/genes ------------------------------
datExpr <- t(log_expr)
gsg <- goodSamplesGenes(datExpr, verbose = 1)
cat(sprintf("Samples with missing entries: %d, genes: %d\n",
            sum(!gsg$goodSamples), sum(!gsg$goodGenes)))
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
cat(sprintf("datExpr: %d samples x %d genes\n", nrow(datExpr), ncol(datExpr)))

# ---- 2. Soft-threshold power ---------------------------------------------
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         networkType = "signed", verbose = 0)

sft_df <- data.frame(
  power            = sft$fitIndices$Power,
  scale_free_R2    = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
  mean_connectivity = sft$fitIndices$mean.k.
)
write.csv(sft_df, res_path("01_discovery", "WGCNA_soft_threshold.csv"),
          row.names = FALSE)
print(sft_df)

power <- sft_df$power[which(sft_df$scale_free_R2 > 0.85)[1]]
if (is.na(power)) {
  warning("No power reached R^2 > 0.85; falling back to power 9.")
  power <- 9
}
cat(sprintf("\nSelected soft-threshold power: %d\n", power))

# ---- 3. Build the network -------------------------------------------------
cat("Building network (20-40 min) ...\n")
net <- blockwiseModules(
  datExpr,
  power             = power,
  networkType       = "signed",
  TOMType           = "signed",
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 2
)

module_colors <- labels2colors(net$colors)
cat("\nModule sizes:\n"); print(table(module_colors))

# ---- 4. Module-trait correlation -----------------------------------------
MEs <- net$MEs
fibro_binary <- as.numeric(group == "High")

module_trait_r <- WGCNA::cor(MEs, fibro_binary, use = "p")
module_trait_p <- corPvalueStudent(module_trait_r, nrow(datExpr))

mt <- data.frame(
  module = rownames(module_trait_r),
  r      = as.numeric(module_trait_r[, 1]),
  P      = as.numeric(module_trait_p[, 1]),
  n_genes = as.integer(table(net$colors)[
    as.character(as.numeric(gsub("ME", "", rownames(module_trait_r))))])
)
mt <- mt[order(-mt$r), ]
write.csv(mt, res_path("01_discovery", "WGCNA_module_trait.csv"), row.names = FALSE)

cat("\nModules ranked by correlation with fibrosis:\n")
print(head(mt, 10))
cat(sprintf("\nNumber of modules (excluding the grey/unassigned module): %d\n",
            sum(mt$module != "ME0")))

# ---- 5. Top positive and top negative module ------------------------------
top_pos <- mt$module[which.max(mt$r)]
top_neg <- mt$module[which.min(mt$r)]
cat(sprintf("Top positive module: %s (r = %.3f, P = %.2e, n = %d genes)\n",
            top_pos, mt$r[mt$module == top_pos], mt$P[mt$module == top_pos],
            mt$n_genes[mt$module == top_pos]))
cat(sprintf("Top negative module: %s (r = %.3f, P = %.2e, n = %d genes)\n",
            top_neg, mt$r[mt$module == top_neg], mt$P[mt$module == top_neg],
            mt$n_genes[mt$module == top_neg]))

module_genes <- function(label) {
  k <- as.numeric(gsub("ME", "", label))
  colnames(datExpr)[net$colors == k]
}

saveRDS(list(
  net               = net,
  module_colors     = module_colors,
  MEs               = MEs,
  module_trait_r    = module_trait_r,
  module_trait_p    = module_trait_p,
  module_table      = mt,
  top_positive      = top_pos,
  top_negative      = top_neg,
  genes_positive    = module_genes(top_pos),
  genes_negative    = module_genes(top_neg),
  power             = power,
  datExpr           = datExpr
), file.path(DATA_DIR, "derived", "WGCNA_results.rds"))

cat("\nNext: 01_discovery/04_candidate_genes.R\n")
