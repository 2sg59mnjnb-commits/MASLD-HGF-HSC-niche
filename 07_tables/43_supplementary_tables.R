# ============================================================================
# 43_supplementary_tables.R
#
# Assembles the supplementary tables as CSV files from the analysis outputs.
# Each table is written only if its input exists; the index at the end records
# what was produced and what was skipped.
#
# Output
#   results/tables/TableS*.csv
#   results/tables/_INDEX.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

OUT <- file.path(RESULTS_DIR, "tables")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

index <- data.frame(table = character(), status = character(),
                    n_rows = integer(), n_cols = integer(),
                    stringsAsFactors = FALSE)

put <- function(df, name) {
  write.csv(df, file.path(OUT, paste0(name, ".csv")), row.names = FALSE)
  index <<- rbind(index, data.frame(table = name, status = "written",
                                    n_rows = nrow(df), n_cols = ncol(df)))
  cat(sprintf("  %-42s %d rows\n", name, nrow(df)))
}

skip <- function(name, why) {
  index <<- rbind(index, data.frame(table = name, status = paste("skipped:", why),
                                    n_rows = NA, n_cols = NA))
  cat(sprintf("  %-42s skipped (%s)\n", name, why))
}

get <- function(...) {
  p <- res_path(...)
  if (!file.exists(p)) return(NULL)
  read.csv(p, stringsAsFactors = FALSE)
}

cat("=== Supplementary tables ===\n")

# ---- S1 Cohort screening (study design table, fixed content) -------------
put(data.frame(
  Dataset = c("GSE135251", "GSE49541", "GSE130970", "GSE84044",
              "GSE126848", "GSE33814"),
  N = c(216, 72, 78, 124, 57, 40),
  Platform = c("RNA-seq", "Affymetrix HG-U133 Plus 2.0", "Illumina RNA-seq",
               "Affymetrix HG-U133 Plus 2.0", "RNA-seq", "Microarray"),
  Signature_coverage = c("11/11", "11/11", "11/11", "partial", "n/a", "n/a"),
  Fibrosis_staging = c("METAVIR F0-F4", "METAVIR F0-F1 vs F3-F4",
                       "METAVIR F0-F4", "Scheuer s0-s4", "none",
                       "incomplete"),
  Role = c("Training / discovery", "Primary external validation",
           "Exploratory validation", "Trend analysis only",
           "Excluded", "Excluded"),
  Reason = c("-", "-", "Different clinical context",
             "Incomplete signature coverage", "No fibrosis staging",
             "Incomplete staging"),
  stringsAsFactors = FALSE
), "TableS1_cohort_screening")

# ---- S2 11-gene signature ----------------------------------------------
put(data.frame(gene = SIG11,
               selection = "intersection of LASSO, random forest, SVM-RFE and XGBoost",
               stringsAsFactors = FALSE), "TableS2_signature_11genes")

# ---- S3 Model coefficients ---------------------------------------------
co <- get("02_model", "model_coefficients.csv")
if (!is.null(co)) put(co, "TableS3_model_coefficients") else skip("TableS3_model_coefficients", "run 02_model/06")

# ---- S4 16-gene signature ----------------------------------------------
sp <- get("04_single_cell", "signature_16genes.csv")
if (!is.null(sp)) {
  put(sp[, c("gene", "logFC_S1", "logFC_S2")], "TableS4_signature_16genes")
} else skip("TableS4_signature_16genes", "run 04_single_cell/25")

# ---- S5 Niche scan -----------------------------------------------------
ns <- get("03_spatial", "niche_scan_by_slide.csv")
if (!is.null(ns)) put(ns, "TableS5_niche_scan") else skip("TableS5_niche_scan", "run 03_spatial/15")

# ---- S6 IL6-family receptor positive rates -----------------------------
rp <- get("04_single_cell", "IL6_receptor_positive_rates.csv")
if (!is.null(rp)) put(rp, "TableS6_IL6_receptor_positive_rates") else skip("TableS6", "run 04_single_cell/24")

# ---- S7 HGF threshold sensitivity --------------------------------------
th <- get("03_spatial", "hgf_threshold_effects.csv")
if (!is.null(th)) put(th, "TableS7_HGF_threshold_sensitivity") else skip("TableS7", "run 03_spatial/17")

# ---- S8 Decision curve --------------------------------------------------
dc <- get("02_model", "DCA_GSE49541.csv")
if (!is.null(dc)) put(dc, "TableS8_decision_curve") else skip("TableS8", "run 02_model/09")

# ---- S9 TF activity -----------------------------------------------------
tf <- get("04_single_cell", "TF_activity_consensus.csv")
if (!is.null(tf)) put(tf, "TableS9_TF_activity") else skip("TableS9", "run 04_single_cell/29")

# ---- S10 Drug prediction -----------------------------------------------
dr <- get("05_downstream", "DGIdb_drug_summary.csv")
if (!is.null(dr)) put(dr, "TableS10_drug_prediction") else skip("TableS10", "run 05_downstream/36")

# ---- S11 GSE130970 validation ------------------------------------------
g130 <- get("02_model", "GSE130970_validation.csv")
if (!is.null(g130)) put(g130, "TableS11_GSE130970_validation") else skip("TableS11", "run 02_model/11")

# ---- S12 Nested CV stability -------------------------------------------
st <- get("02_model", "nested_cv_stability.csv")
if (!is.null(st)) put(st, "TableS12_nested_CV_stability") else skip("TableS12", "run 02_model/07")

# ---- S13 TMM sensitivity ------------------------------------------------
skip("TableS13_TMM_sensitivity",
     "normalisation sensitivity was not run in this release; see README")

# ---- S14 MET pathway sensitivity ---------------------------------------
skip("TableS14_MET_pathway_sensitivity",
     "MET pathway analysis excluded from this release; see README")

# ---- S15 CHI3L1 comparison ---------------------------------------------
cmp <- get("02_model", "CHI3L1_comparison.csv")
if (!is.null(cmp)) put(cmp, "TableS15_CHI3L1_comparison") else skip("TableS15", "run 02_model/10")

# ---- S16 Incremental value ---------------------------------------------
iv <- get("05_downstream", "incremental_value_GSE89632.csv")
if (!is.null(iv)) put(iv, "TableS16_incremental_value") else skip("TableS16", "run 05_downstream/35")

# ---- S17 Depth correction ----------------------------------------------
dc17 <- get("03_spatial", "depth_corrected_distance.csv")
if (!is.null(dc17)) put(dc17, "TableS17_depth_correction") else skip("TableS17", "run 03_spatial/20")

# ---- S18 CHI3L1 DeLong --------------------------------------------------
dl <- get("02_model", "CHI3L1_DeLong.csv")
if (!is.null(dl)) put(dl, "TableS18_CHI3L1_DeLong") else skip("TableS18", "run 02_model/10")

# ---- Index --------------------------------------------------------------
write.csv(index, file.path(OUT, "_INDEX.csv"), row.names = FALSE)
cat(sprintf("\n%d tables written, %d skipped\n",
            sum(index$status == "written"), sum(index$status != "written")))
cat(sprintf("Output: %s\n", OUT))

# ---- Session information -------------------------------------------------
si <- capture.output(sessionInfo())
writeLines(c(paste("R", getRversion()), "", si), "sessionInfo.txt")
cat("sessionInfo.txt written\n")
