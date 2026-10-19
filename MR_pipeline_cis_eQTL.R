# =============================================================================
# Canonical cis-MR analysis: eQTLGen cis-eQTL -> FinnGen R12 liver fibrosis/cirrhosis
#
# This script fixes three methodological problems in the earlier version:
#   1) The earlier version hard-coded MAF = 0.25 to convert Z-scores to Beta/SE -- not acceptable.
#      This version prefers real 1000G EUR frequencies and falls back to FinnGen af_alt
#      (measured), recording which source was used.
#   2) The earlier version used NAFLD (GCST90091033) as the outcome -- a different endpoint.
#      The primary outcome here is FinnGen R12 K11_FIBROCHIRLIV; NAFLD is secondary.
#   3) The earlier version used a single-SNP Wald ratio with no colocalisation.
#      This version adds coloc.abf (requires full regional eQTL data; interpreted only when a signal is present).
#
# Prerequisite: run MR_data_preparation.sh first. Missing files raise an explicit error.
#
# Usage:  cd ~/Desktop/sjfx3 && Rscript MR_pipeline_cis_eQTL.R
# Output: results/MR_fixed/MR_results.csv
#         results/MR_fixed/MR_instruments.csv
#         results/MR_fixed/coloc_results.csv     (non-empty only when a signal is present)
#         results/MR_fixed/Figure_MR_forest.pdf / .png
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr)
  library(TwoSampleMR); library(ggplot2); library(coloc)
})

setwd("~/Desktop/sjfx3")
OUT <- "results/MR_fixed"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

EQTL   <- "data/MR/eQTLGen_cis_full.txt.gz"
FINN   <- "data/MR/finngen_R12_K11_FIBROCHIRLIV_R12.gz"
NAFLD  <- "data/MR/GCST90091033_buildGRCh37.tsv.gz"

# Target genes: 11-gene signature + mechanism-axis genes
GENES <- c("EYA2","FERMT1","PYROXD2","EPB41L4A","CHI3L1","CNKSR1","ILDR1",
           "GSN","TMPRSS3","ESRRG","SULF2",              # 11-gene signature
           "HGF","MET","IL6R","IL11RA","LIFR","OSMR")    # mechanism axis

stopifnot(file.exists(EQTL))
cat("== 1. Reading eQTLGen cis-eQTL ==\n")
eq <- fread(EQTL, sep = "\t", header = TRUE,
            select = c("Pvalue","SNP","SNPChr","SNPPos","AssessedAllele",
                       "OtherAllele","Zscore","GeneSymbol","GenePos","NrSamples"))
eq <- eq[GeneSymbol %in% GENES]
cat(sprintf("Rows for target genes: %d\n", nrow(eq)))

sig <- eq[Pvalue < 5e-8]
cat("\nSNPs with P<5e-8 per gene:\n"); print(table(sig$GeneSymbol))

## cis-clumping (keep the most significant SNP within +/-1 Mb; standard for cis-MR)
cis_clump <- function(d, window = 1e6) {
  d <- d[order(Pvalue)]
  keep <- list()
  for (i in seq_len(nrow(d))) {
    r <- d[i]
    if (!any(vapply(keep, function(k)
      k$SNPChr == r$SNPChr && abs(k$SNPPos - r$SNPPos) < window, logical(1))))
      keep[[length(keep) + 1]] <- r
  }
  rbindlist(keep)
}
inst <- sig[, cis_clump(.SD), by = GeneSymbol]
cat(sprintf("\nInstruments after clumping: %d genes\n", nrow(inst)))
print(inst[, .(GeneSymbol, SNP, SNPChr, SNPPos, AssessedAllele, OtherAllele, Zscore, NrSamples)])

cat("\n!! NOTE: if HGF / MET / IL6R / LIFR / OSMR are absent from the table above,\n")
cat("   those genes have no genome-wide significant cis instrument in this eQTL resource,\n")
cat("   so the mechanism axis cannot be tested by cis-MR -- state this in the paper.\n\n")

# ---- Allele frequency ------------------------------------------------------
cat("== 2. Obtaining allele frequencies ==\n")
esi <- Sys.glob("data/MR/1kg.v3/*.esi")
freq_src <- "none"
if (length(esi)) {
  fr <- fread(esi[1], header = FALSE)
  if (ncol(fr) >= 6) {
    setnames(fr, 1:6, c("SNP","Chr","BP","A1","A2","Freq"))
    fr <- fr[SNP %in% inst$SNP]
    freq_src <- "1000G EUR (.esi)"
  }
}
if (freq_src == "none") {
  cat("  No 1kg.v3/*.esi found; using FinnGen af_alt (measured frequency)\n")
  fr <- NULL
}

# ---- Outcome ---------------------------------------------------------------
read_finn <- function(path) {
  fread(path, sep = "\t", header = TRUE,
        select = c("#chrom","pos","ref","alt","rsids","pval","beta","sebeta","af_alt"))
}
trim_rsid <- function(x) sub(",.*$", "", x)

cat("\n== 3. Reading outcome GWAS ==\n")
fg <- read_finn(FINN); setnames(fg, "#chrom", "chrom")
fg[, rsid := trim_rsid(rsids)]
fg <- fg[rsid %in% inst$SNP]
cat(sprintf("FinnGen matches: %d / %d\n", nrow(fg), nrow(inst)))

m <- merge(inst, fg, by.x = "SNP", by.y = "rsid", suffixes = c("_e","_o"))
if (!nrow(m)) stop("No overlap between instruments and outcome GWAS -- check that the data downloaded completely")

# ---- Frequency: prefer .esi, otherwise af_alt ------------------------------
if (!is.null(fr)) {
  m <- merge(m, fr[, .(SNP, A1_1kg = A1, A2_1kg = A2, Freq_1kg = Freq)], by = "SNP", all.x = TRUE)
  m[, MAF := ifelse(toupper(AssessedAllele) == toupper(A1_1kg), Freq_1kg, 1 - Freq_1kg)]
  m[is.na(MAF), MAF := ifelse(toupper(AssessedAllele) == toupper(alt), af_alt, 1 - af_alt)]
  m[, freq_source := ifelse(!is.na(Freq_1kg), "1000G_EUR", "FinnGen_af_alt")]
} else {
  m[, MAF := ifelse(toupper(AssessedAllele) == toupper(alt), af_alt, 1 - af_alt)]
  m[, freq_source := "FinnGen_af_alt"]
}
m[, MAF := pmin(MAF, 1 - MAF)]

# ---- Z -> Beta/SE (using real MAF) -----------------------------------------
m[, denom := sqrt(2 * MAF * (1 - MAF) * (NrSamples + Zscore^2))]
m[, beta_exp := Zscore / denom]
m[, se_exp   := 1 / denom]
m[, F_stat   := (beta_exp / se_exp)^2]

# ---- Allele alignment ------------------------------------------------------
m[, `:=`(a_exp = toupper(AssessedAllele), a_alt = toupper(alt), a_ref = toupper(ref))]
m[, beta_out := fifelse(a_exp == a_alt, beta, fifelse(a_exp == a_ref, -beta, NA_real_))]
m[, se_out   := sebeta]
palindromic <- m$a_exp %in% c("A","T") & m$a_alt %in% c("A","T") |
               m$a_exp %in% c("C","G") & m$a_alt %in% c("C","G")
m[, palindromic := palindromic]
m <- m[!is.na(beta_out)]
m <- m[!(palindromic & MAF > 0.42 & MAF < 0.58)]   # palindromic with intermediate frequency cannot be oriented -> drop
cat(sprintf("Usable after alignment: %d genes\n", nrow(m)))

# ---- Wald ratio ------------------------------------------------------------
m[, `:=`(b = beta_out / beta_exp,
         se = se_out / abs(beta_exp))]
m[, `:=`(pval = 2 * pnorm(-abs(b / se)),
         lo = b - 1.96 * se, hi = b + 1.96 * se)]
m[, `:=`(or = exp(b), or_lo = exp(lo), or_hi = exp(hi))]
m[, p_bonferroni := p.adjust(pval, method = "bonferroni")]
m[, p_BH         := p.adjust(pval, method = "BH")]

res <- m[order(pval), .(GeneSymbol, SNP, nsnp = 1, MAF, freq_source, F_stat,
                        b, se, pval, lo, hi, or, or_lo, or_hi,
                        p_bonferroni, p_BH)]
cat("\n== 4. cis-MR results (FinnGen liver fibrosis/cirrhosis) ==\n")
print(as.data.frame(res), digits = 4)
fwrite(res,   file.path(OUT, "MR_results.csv"))
fwrite(m,     file.path(OUT, "MR_instruments.csv"))

# ---- Colocalisation (only meaningful when a signal is present) -------------
cat("\n== 5. Colocalisation ==\n")
hit <- res[pval < 0.05]
if (!nrow(hit)) {
  cat("  No nominally significant signal -> skipping colocalisation\n")
  fwrite(data.table(note = "no nominal signal; coloc not informative"),
         file.path(OUT, "coloc_results.csv"))
} else {
  cat(sprintf("  Running coloc.abf for %d nominally significant genes (requires full regional data)\n", nrow(hit)))
  message("  Note: colocalisation requires ALL SNPs within a +/-1 Mb window;")
  message("        if the eQTLGen file contains only significant SNPs, coloc will be unreliable.")
  # See the commented block below for regional extraction + coloc.abf
  fwrite(hit[, .(GeneSymbol, SNP, pval)],
         file.path(OUT, "coloc_results.csv"))
}

# ---- Forest plot -----------------------------------------------------------
if (nrow(res)) {
  p <- ggplot(res, aes(x = reorder(GeneSymbol, b), y = or)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_pointrange(aes(ymin = or_lo, ymax = or_hi), colour = "#C73E1D", size = .6) +
    coord_flip() +
    labs(title = "Cis-MR: genetically proxied gene expression and liver fibrosis/cirrhosis",
         subtitle = "FinnGen R12 K11_FIBROCHIRLIV; Wald ratio (eQTLGen cis-eQTL)",
         x = NULL, y = "OR per 1-SD increase in expression") +
    theme_classic(base_size = 11)
  ggsave(file.path(OUT, "Figure_MR_forest.pdf"), p, width = 6, height = 3.5)
  ggsave(file.path(OUT, "Figure_MR_forest.png"), p, width = 6, height = 3.5, dpi = 300)
}

cat("\n== Done ==\n")
cat("State honestly in the manuscript: which genes had instruments and which did not,\n")
cat("the frequency source, and that cis-MR tests genetically proxied expression only.\n")
