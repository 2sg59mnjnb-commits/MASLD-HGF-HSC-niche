# ============================================================================
# helpers.R -- shared functions
#
# Source after R/config.R. Requires: Seurat, dplyr, ggplot2.
# All functions are written to work with Seurat v5 objects.
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# ---- Small utilities -------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# Row-wise z-score; rows with zero variance or NAs become 0.
zscore_rows <- function(m) {
  z <- t(scale(t(as.matrix(m))))
  z[is.na(z)] <- 0
  z
}

# Mean of row-wise z-scores = the signature score used throughout.
sig_score <- function(expr, genes) {
  g <- intersect(genes, rownames(expr))
  if (length(g) == 0) return(rep(NA_real_, ncol(expr)))
  colMeans(zscore_rows(expr[g, , drop = FALSE]))
}

# Drop the unnamed one-letter clusters (a-j) that the original CosMx
# annotation contains; they have no lineage label.
drop_unnamed <- function(x) {
  x[!grepl("^[a-j]$", x)]
}

# ---- Plotting --------------------------------------------------------------

theme_pub <- function(base = 10) {
  theme_classic(base_size = base) +
    theme(plot.title    = element_text(face = "bold", size = base + 1),
          plot.subtitle = element_text(size = base - 1, colour = "grey30"),
          plot.tag      = element_text(face = "bold", size = base + 2))
}

# Save a ggplot as both PDF and PNG. Uses the default pdf() device (no cairo
# dependency) so it works on a bare Linux install.
safe_save <- function(p, name, width, height, dpi = 600,
                      dir = FIG_DIR, formats = c("pdf", "png")) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  if ("pdf" %in% formats) {
    f <- file.path(dir, paste0(name, ".pdf"))
    tryCatch({ ggsave(f, p, width = width, height = height, device = "pdf")
               cat(sprintf("  saved %s\n", basename(f))) },
             error = function(e) cat(sprintf("  FAILED %s: %s\n", basename(f), conditionMessage(e))))
  }
  if ("png" %in% formats) {
    f <- file.path(dir, paste0(name, ".png"))
    tryCatch({ ggsave(f, p, width = width, height = height, dpi = dpi,
                      device = "png")
               cat(sprintf("  saved %s\n", basename(f))) },
             error = function(e) cat(sprintf("  FAILED %s: %s\n", basename(f), conditionMessage(e))))
  }
}

# ---- Data loading ----------------------------------------------------------

# CosMx slide 1 or 2, with the patient label recovered from FOV x-coordinates.
# GSE292268 ships two slides, each holding two patients; the metadata carry no
# patient column, so it is reconstructed from the clear spatial gap at x = 5 mm.
load_cosmx <- function(slide = 1, add_patient_label = TRUE) {
  f <- file.path(DATA_DIR, sprintf("GSM%s_Slide%d_SeuratObject.RDS",
                                   ifelse(slide == 1, "8855042", "8855043"), slide))
  if (!file.exists(f)) {
    stop(sprintf("Missing %s\nDownload it from GSE292268 (see README).", f))
  }
  obj <- readRDS(f)
  if (add_patient_label) obj <- add_patient(obj)
  obj
}

add_patient <- function(obj, x_cutoff = PATIENT_X_CUTOFF) {
  fx <- obj@meta.data %>%
    dplyr::group_by(fov) %>%
    dplyr::summarise(xc = mean(x_slide_mm), .groups = "drop") %>%
    dplyr::mutate(patient = ifelse(xc < x_cutoff, "A", "B"))
  obj@meta.data$patient <- fx$patient[match(obj@meta.data$fov, fx$fov)]
  obj
}

celltype_col <- function(slide) {
  if (slide == 1) CELLTYPE_COL_S1 else CELLTYPE_COL_S2
}

# Flatten a CosMx object into a plain data.frame of coordinates + selected genes.
# Keeps memory low: the expression matrix is subset before conversion.
flatten_cosmx <- function(obj, slide, genes = c("HGF", "IL6R", "MET")) {
  md  <- obj@meta.data
  dat <- obj[["RNA"]]$data
  genes <- intersect(genes, rownames(dat))
  out <- data.frame(
    cell    = colnames(obj),
    x       = md$x_slide_mm,
    y       = md$y_slide_mm,
    fov     = md$fov,
    patient = md$patient,
    slide   = paste0("S", slide),
    celltype = as.character(md[[celltype_col(slide)]]),
    nCount  = md$nCount_RNA,
    nFeature = md$nFeature_RNA,
    stringsAsFactors = FALSE
  )
  for (g in genes) out[[g]] <- as.numeric(dat[g, ])
  out
}

# Merge both slides into one data.frame (genes must be shared).
load_cosmx_flat <- function(genes = c("HGF", "IL6R", "MET", "STAT3"),
                            slides = c(1, 2)) {
  parts <- lapply(slides, function(s) {
    obj <- load_cosmx(s)
    d <- flatten_cosmx(obj, s, genes)
    rm(obj); gc(verbose = FALSE)
    d
  })
  df <- dplyr::bind_rows(parts)
  df$patient_uid <- paste(df$slide, df$patient, sep = "_")
  df
}

# Bulk training cohort, saved by 01_discovery/02_count_matrix_and_DEG.R
load_bulk <- function() {
  f <- file.path(DATA_DIR, "derived", "GSE135251_DEG_results.rds")
  if (!file.exists(f)) stop("Run 01_discovery/02_count_matrix_and_DEG.R first.")
  readRDS(f)
}

# Map Ensembl rownames to symbols, dropping unmapped and duplicated entries.
to_symbol <- function(expr) {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = rownames(expr),
                               column = "SYMBOL", keytype = "ENSEMBL",
                               multiVals = "first")
  keep <- !is.na(sym) & !duplicated(sym)
  out <- as.matrix(expr[keep, , drop = FALSE])
  rownames(out) <- sym[keep]
  out
}

# ---- Spatial geometry ------------------------------------------------------

# FOV labels are NOT unique across slides: both CosMx slides number their fields
# 1..64. Any grouping by field of view must therefore combine slide and fov,
# otherwise fields from the two slides are merged into one group.
fov_key <- function(d) {
  if ("slide" %in% names(d)) paste(d$slide, d$fov, sep = ":") else as.character(d$fov)
}

# Distance (um) from every cell in `src` to the nearest cell in `tgt`,
# computed within the same FOV only. Both inputs need x, y (mm), fov columns.
compute_nn <- function(src, tgt) {
  out <- rep(NA_real_, nrow(src))
  ks <- fov_key(src)
  kt <- fov_key(tgt)
  for (k in unique(ks)) {
    si <- which(ks == k)
    ti <- which(kt == k)
    if (length(si) == 0 || length(ti) == 0) next
    s <- src[si, c("x", "y")]
    t <- tgt[ti, c("x", "y")]
    d <- sqrt(outer(s$x, t$x, "-")^2 + outer(s$y, t$y, "-")^2)
    out[si] <- apply(d, 1, min) * 1000   # mm -> um
  }
  out
}

# Number of target cells within `radius_um` of every source cell, same FOV.
count_within <- function(src, tgt, radius_um = NICHE_RADIUS_UM) {
  r <- radius_um / 1000
  out <- rep(0L, nrow(src))
  ks <- fov_key(src)
  kt <- fov_key(tgt)
  for (k in unique(ks)) {
    si <- which(ks == k)
    ti <- which(kt == k)
    if (length(si) == 0 || length(ti) == 0) next
    s <- src[si, c("x", "y")]
    t <- tgt[ti, c("x", "y")]
    d2 <- outer(s$x, t$x, "-")^2 + outer(s$y, t$y, "-")^2
    out[si] <- rowSums(d2 <= r^2)
  }
  out
}

# FOV area in mm^2, from the coordinate ranges of the cells it contains.
fov_area <- function(df) {
  if (!"slide" %in% names(df)) df$slide <- "S1"
  df %>% dplyr::group_by(slide, fov) %>%
    dplyr::summarise(area = diff(range(x)) * diff(range(y)),
                     n = dplyr::n(), .groups = "drop")
}

# ---- Mixed-model reporting -------------------------------------------------

# Fit lmer and return a one-row data.frame with beta, Wald 95% CI and P.
# P-value comes from lmerTest (Satterthwaite df).
lmm_effect <- function(formula, data, term, label) {
  suppressPackageStartupMessages({library(lme4); library(lmerTest)})
  m <- lmerTest::lmer(formula, data = data, REML = FALSE)
  ct <- summary(m)$coefficients
  b  <- ct[term, "Estimate"]; se <- ct[term, "Std. Error"]
  data.frame(
    Model   = label,
    Term    = term,
    Beta    = round(b, 4),
    CI_low  = round(b - 1.96 * se, 4),
    CI_high = round(b + 1.96 * se, 4),
    t_value = round(ct[term, "t value"], 3),
    df      = round(ct[term, "df"], 1),
    P       = ct[term, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}

# Intraclass correlation for a random-intercept model.
icc_lmm <- function(formula, data) {
  suppressPackageStartupMessages({library(lme4); library(lmerTest)})
  m  <- lmerTest::lmer(formula, data = data, REML = TRUE)
  vc <- as.data.frame(lme4::VarCorr(m))
  vp <- vc$vcov[vc$grp != "Residual"][1]
  vr <- vc$vcov[vc$grp == "Residual"]
  c(icc = vp / (vp + vr), var_random = vp, var_resid = vr)
}

# ---- Reporting -------------------------------------------------------------

# Print an observed vs reported comparison; never feeds back into statistics.
compare <- function(label, observed, reported, digits = 3) {
  cat(sprintf("  %-38s observed %-10s reported %s\n",
              label,
              formatC(observed, format = "f", digits = digits),
              formatC(reported, format = "f", digits = digits)))
}

# ---- Loading derived objects --------------------------------------------

# Which script produces each derived object, so a missing file produces an
# actionable message instead of a bare "cannot open connection".
DERIVED_HINTS <- c(
  "GSE135251_aligned_data.rds" = "01_discovery/02_count_matrix_and_DEG.R",
  "GSE135251_DEG_results.rds"  = "01_discovery/02_count_matrix_and_DEG.R",
  "WGCNA_results.rds"          = "01_discovery/03_wgcna.R",
  "ML_inputs.rds"              = "01_discovery/05_four_algorithm_selection.R",
  "fixed_cv.rds"               = "02_model/06_fixed_signature_cv.R",
  "nested_cv.rds"              = "02_model/07_fully_nested_cv.R",
  "GSE49541_validation.rds"    = "02_model/08_external_validation_GSE49541.R",
  "secondary_performance.rds"  = "02_model/09_calibration_DCA_Firth.R",
  "CHI3L1_comparison.rds"      = "02_model/10_CHI3L1_comparison.R",
  "GSE130970_validation.rds"   = "02_model/11_GSE130970_exploratory.R",
  "cosmx_flat.rds"             = "03_spatial/13_load_CosMx_and_panel.R",
  "stellate_enrichment.rds"    = "03_spatial/14_stellate_fibrotic_enrichment.R",
  "niche_scan.rds"             = "03_spatial/15_HGFpos_niche_scan.R",
  "proximity_models.rds"       = "03_spatial/16_HGFpos_proximity_mixed_models.R",
  "zonation.rds"               = "03_spatial/18_zonation_portal_vs_central.R",
  "hep4_il6r.rds"              = "03_spatial/19_Hep4_IL6R_and_pathways.R",
  "depth_correction.rds"       = "03_spatial/20_depth_correction.R",
  "pseudobulk.rds"             = "04_single_cell/25_pseudobulk_16gene_signature.R",
  "hep4_identity.rds"          = "04_single_cell/27_Hep4_identity_scRNA.R",
  "signature_specificity.rds"  = "04_single_cell/28_signature_specificity_scRNA.R",
  "tf_activity.rds"            = "04_single_cell/29_TF_activity.R",
  "stage_trends.rds"           = "05_downstream/31_stage_trends.R"
)

load_derived <- function(name, hint = NULL) {
  f <- file.path(DATA_DIR, "derived", name)
  if (!file.exists(f)) {
    hint <- hint %||% DERIVED_HINTS[[name]] %||% "the earlier pipeline stage"
    stop(sprintf("Missing derived object: %s\nRun %s first.", f, hint),
         call. = FALSE)
  }
  readRDS(f)
}
