# ============================================================================
# 37_MR_exploratory.R
#
# Exploratory cis-Mendelian randomisation, included for transparency.
#
# STATUS: this analysis is NOT reported in the manuscript. It was run, it did not
# support a causal interpretation, and it is retained here so that the negative
# result is on record rather than silently dropped.
#
# What happens and why it fails
#   - cis-eQTLs for CHI3L1 from eQTLGen are in strong linkage disequilibrium.
#   - After distance-based pruning at 1 Mb a single variant remains, so the
#     analysis reduces to a single-instrument Wald ratio.
#   - That single-instrument estimate is not significant, so no causal claim can
#     be made and the analysis was not carried forward.
#
# Method notes
#   - eQTLGen distributes Z-scores rather than beta/SE. Beta and SE are
#     back-computed assuming an allele frequency of 0.25, which is an
#     approximation; it does not affect the sign of the estimate but does affect
#     its standard error. This is another reason the result is exploratory only.
#   - Harmonisation is done manually (allele matching, sign flips) instead of via
#     a package, because the packaged harmoniser drops these variants.
#
# Inputs (download separately; not redistributed)
#   data/MR/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
#   data/MR/GCST90091033_buildGRCh37.tsv.gz
#
# Output
#   results/05_downstream/MR_CHI3L1.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(data.table)})

eqtl_file <- file.path(DATA_DIR, "MR",
  "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
gwas_file <- file.path(DATA_DIR, "MR", "GCST90091033_buildGRCh37.tsv.gz")

if (!file.exists(eqtl_file) || !file.exists(gwas_file)) {
  cat("MR input files not found; skipping.\n")
  cat("See the header of this script for the expected filenames.\n")
  quit(save = "no", status = 0)
}

# ---- 1. Exposure: cis-eQTLs for the gene of interest --------------------
eqtl <- fread(eqtl_file, select = c("SNP", "SNPChr", "SNPPos", "GeneSymbol",
                                    "AssessedAllele", "OtherAllele",
                                    "Zscore", "Pvalue", "NrSamples"))
exp_dat <- eqtl[GeneSymbol == "CHI3L1" & Pvalue < 5e-8]
cat(sprintf("cis-eQTLs for CHI3L1 at P < 5e-8: %d\n", nrow(exp_dat)))
if (nrow(exp_dat) == 0) {
  cat("No instruments; stopping.\n"); quit(save = "no", status = 0)
}

# Back-compute beta / SE from the Z-score (see method note above)
MAF <- 0.25
exp_dat[, denom := sqrt(2 * MAF * (1 - MAF) * (NrSamples + Zscore^2))]
exp_dat[, beta := Zscore / denom]
exp_dat[, se   := 1 / denom]

# ---- 2. Distance-based pruning (1 Mb window) ----------------------------
# Cluster variants that sit within 1 Mb of each other on the same chromosome and
# keep the most significant variant per cluster. Genotype data are not shipped
# with the study, so pruning is positional rather than LD-based; this is noted as
# a limitation of the exploratory analysis.
setDT(exp_dat)
setorder(exp_dat, SNPChr, SNPPos)
exp_dat[, cluster := cumsum(
  c(TRUE, diff(SNPPos) > 1e6 | diff(as.integer(factor(SNPChr))) != 0)
)]
exp_dat <- exp_dat[, .SD[which.min(Pvalue)], by = cluster]
exp_dat[, cluster := NULL]
cat(sprintf("After 1 Mb pruning: %d variant(s)\n", nrow(exp_dat)))

# ---- 3. Outcome ---------------------------------------------------------
gwas <- fread(gwas_file)
need <- c("variant_id", "beta", "standard_error", "effect_allele",
          "other_allele", "p_value")
if (!all(need %in% colnames(gwas))) {
  cat("Unexpected GWAS columns; expected: ", paste(need, collapse = ", "), "\n")
  quit(save = "no", status = 0)
}

m <- merge(exp_dat, gwas, by.x = "SNP", by.y = "variant_id")
cat(sprintf("Variants present in the outcome: %d\n", nrow(m)))
if (nrow(m) == 0) { cat("No overlap; stopping.\n"); quit(save = "no", status = 0) }

# ---- 4. Manual harmonisation -------------------------------------------
same <- toupper(m$AssessedAllele) == toupper(m$effect_allele) &
        toupper(m$OtherAllele)    == toupper(m$other_allele)
flip <- toupper(m$AssessedAllele) == toupper(m$other_allele) &
        toupper(m$OtherAllele)    == toupper(m$effect_allele)

m$beta_out_aligned <- ifelse(same, m$beta,
                             ifelse(flip, -m$beta, NA_real_))
m <- m[!is.na(m$beta_out_aligned)]
cat(sprintf("Variants after allele alignment: %d\n", nrow(m)))
if (nrow(m) == 0) { cat("No alignable variants; stopping.\n"); quit(save = "no", status = 0) }

# ---- 5. Estimate --------------------------------------------------------
if (nrow(m) >= 3) {
  w <- 1 / m$se^2
  beta_ivw <- sum(w * m$beta * m$beta_out_aligned) / sum(w * m$beta^2)
  se_ivw <- sqrt(1 / sum(w * m$beta^2))
  method <- "IVW"
} else {
  beta_ivw <- m$beta_out_aligned[1] / m$beta[1]
  se_ivw <- abs(m$se[1] / m$beta[1])
  method <- "Wald ratio (single instrument)"
}
p <- 2 * pnorm(-abs(beta_ivw / se_ivw))

cat(sprintf("\nMethod: %s\n", method))
cat(sprintf("beta = %.4f, SE = %.4f, P = %.3f\n", beta_ivw, se_ivw, p))
cat(sprintf("OR per SD of expression = %.3f (95%% CI %.3f-%.3f)\n",
            exp(beta_ivw), exp(beta_ivw - 1.96 * se_ivw),
            exp(beta_ivw + 1.96 * se_ivw)))

if (p >= 0.05) {
  cat("\nResult not significant. With a single instrument this analysis cannot\n")
  cat("distinguish a null effect from insufficient power, and it is reported as\n")
  cat("exploratory only (not included in the manuscript).\n")
}

out <- data.frame(
  exposure = "CHI3L1", outcome = "MASLD (GCST90091033)",
  method = method, n_instruments = nrow(m),
  beta = round(beta_ivw, 4), se = round(se_ivw, 4), P = p,
  OR = round(exp(beta_ivw), 3),
  OR_low = round(exp(beta_ivw - 1.96 * se_ivw), 3),
  OR_high = round(exp(beta_ivw + 1.96 * se_ivw), 3)
)
write.csv(out, res_path("05_downstream", "MR_CHI3L1.csv"), row.names = FALSE)

saveRDS(list(instruments = m, estimate = out),
        file.path(DATA_DIR, "derived", "mr.rds"))

cat("\nNext: 06_figures/38_Figure1.R\n")
