# ============================================================================
# 30_TF_HGF_coupling_and_motif.R
#
# Two orthogonal checks on the TFAP2A / TP53 -> HGF regulatory model:
#
#   Part A. Coupling within single cells. Is TF activity correlated with HGF
#           expression across individual HSCs?
#
#   Part B. Promoter-level support. Are TFAP2A / TP53 motifs present in the HGF
#           promoter more often than expected by chance? JASPAR2020 matrices are
#           scanned against hg38 chr7:81,698,000-81,702,000 (TSS +/- 2 kb) and
#           compared with 100 shuffled sequences.
#
# A regulon-level correlation does not establish direct binding. If Part B shows
# depletion rather than enrichment, the regulatory model is reported as
# hypothesis-generating only, which is what the manuscript does.
#
# Output
#   results/04_single_cell/TF_HGF_coupling.csv
#   results/04_single_cell/motif_enrichment.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

# ============================ Part A: coupling =============================
cat("=== Part A: TF activity vs HGF expression ===\n")
tf <- load_derived("tf_activity.rds")

coupling <- do.call(rbind, lapply(c("slide1", "slide2"), function(s) {
  mat <- tf[[s]]
  slide_no <- ifelse(s == "slide1", 1, 2)

  obj <- load_cosmx(slide_no)
  ct  <- as.character(obj@meta.data[[celltype_col(slide_no)]])
  hgf <- as.numeric(obj[["RNA"]]$data["HGF", which(ct == "Stellate.cells")])
  rm(obj); gc(verbose = FALSE)

  if (nrow(mat) != length(hgf)) {
    cat(sprintf("  %s: cell count mismatch (%d vs %d), skipping\n",
                s, nrow(mat), length(hgf)))
    return(NULL)
  }

  do.call(rbind, lapply(intersect(c("TFAP2A", "TP53", "STAT1", "AR", "SP3"),
                                  colnames(mat)), function(t) {
    ct_test <- cor.test(mat[[t]], hgf, method = "spearman")
    data.frame(slide = s, TF = t,
               rho = round(as.numeric(ct_test$estimate), 4),
               P = ct_test$p.value)
  }))
}))
print(coupling, row.names = FALSE)
write.csv(coupling, res_path("04_single_cell", "TF_HGF_coupling.csv"),
          row.names = FALSE)

# ============================ Part B: motif scan ===========================
cat("\n=== Part B: JASPAR motif scan of the HGF promoter ===\n")

needed <- c("JASPAR2020", "TFBSTools", "Biostrings", "httr", "jsonlite")
if (!all(vapply(needed, requireNamespace, logical(1), quietly = TRUE))) {
  cat("Motif packages missing; skipping Part B.\n")
  cat("Needed: ", paste(needed, collapse = ", "), "\n")
  quit(save = "no", status = 0)
}
suppressPackageStartupMessages({
  library(JASPAR2020); library(TFBSTools); library(Biostrings)
  library(httr); library(jsonlite)
})

# ---- Fetch the promoter from the UCSC API --------------------------------
prom_url <- paste0("https://api.genome.ucsc.edu/getData/sequence?",
                   "genome=hg38&chrom=chr7&start=81698000&end=81702000")
res <- httr::GET(prom_url, httr::timeout(30))
if (httr::status_code(res) != 200) stop("UCSC sequence API unavailable.")
dna <- toupper(jsonlite::fromJSON(httr::content(res, "text", encoding = "UTF-8"))$dna)
cat(sprintf("HGF promoter: %d bp\n", nchar(dna)))

fwd <- DNAString(dna)
rev <- reverseComplement(fwd)

# ---- Build numeric PWMs ---------------------------------------------------
build_pwm <- function(tf) {
  pf <- tryCatch(
    getMatrixSet(JASPAR2020, opts = list(species = 9606, name = tf,
                                         all_versions = FALSE))[[1]],
    error = function(e) NULL)
  if (is.null(pf)) return(NULL)
  m <- as.matrix(pf)
  bases <- c("A", "C", "G", "T")
  full <- matrix(0.25, 4, ncol(m), dimnames = list(bases, colnames(m)))
  for (b in intersect(bases, rownames(m))) full[b, ] <- m[b, ]
  full <- full + 0.25                      # pseudocount
  sweep(full, 2, colSums(full), "/")
}

# ---- Background: 100 shuffled sequences ----------------------------------
set.seed(SEED)
n_bg <- 100
bg <- lapply(seq_len(n_bg), function(i)
  DNAString(paste(sample(c("A", "C", "G", "T"), nchar(dna), replace = TRUE),
                  collapse = "")))

tfs <- c("TFAP2A", "TP53", "TFAP2C", "TFAP2B", "SP3", "STAT1")
thresholds <- c("80%", "75%", "70%", "65%", "60%")

count_hits <- function(pwm, seq, thr) {
  tryCatch(length(Biostrings::matchPWM(pwm, seq, min.score = thr,
                                       with.score = FALSE)),
           error = function(e) 0)
}

motif_tab <- do.call(rbind, lapply(tfs, function(tf) {
  pwm <- build_pwm(tf)
  if (is.null(pwm)) return(NULL)
  do.call(rbind, lapply(thresholds, function(thr) {
    obs <- count_hits(pwm, fwd, thr) + count_hits(pwm, rev, thr)
    bg_hits <- vapply(bg, function(s)
      count_hits(pwm, s, thr) + count_hits(pwm, reverseComplement(s), thr),
      numeric(1))
    data.frame(
      TF = tf, threshold = thr,
      HGF_hits = obs,
      background_mean = round(mean(bg_hits), 2),
      fold_enrichment = round(obs / max(mean(bg_hits), 0.01), 2),
      P_empirical = round((sum(bg_hits >= obs) + 1) / (n_bg + 1), 4)
    )
  }))
}))

print(motif_tab, row.names = FALSE)
write.csv(motif_tab, res_path("04_single_cell", "motif_enrichment.csv"),
          row.names = FALSE)

cat("\nSummary at the most permissive threshold (60%):\n")
print(motif_tab[motif_tab$threshold == "60%",
                c("TF", "HGF_hits", "fold_enrichment", "P_empirical")],
      row.names = FALSE)

cat("\nReading: depletion of TFAP2A and TP53 motifs relative to the shuffled\n")
cat("background means the regulon-level association is not supported at the\n")
cat("promoter. The manuscript therefore labels this analysis hypothesis-generating.\n")

saveRDS(list(coupling = coupling, motif = motif_tab,
             promoter = dna, background = bg),
        file.path(DATA_DIR, "derived", "tf_coupling_motif.rds"))

cat("\nNext: 05_downstream/31_stage_trends.R\n")
