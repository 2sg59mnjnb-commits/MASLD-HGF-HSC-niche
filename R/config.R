# ============================================================================
# config.R -- central configuration for the MASLD HGF+ HSC niche study
#
# Every analysis script starts with
#     source("R/config.R"); source("R/helpers.R")
# and is run with the repository root as the working directory.
#
# All paths below are relative to the repository root. Nothing is written
# outside data/, results/ and figures/.
# ============================================================================

stopifnot(file.exists("R/helpers.R"))

# ---- Directories -----------------------------------------------------------
DATA_DIR     <- "data"
RESULTS_DIR  <- "results"
FIG_DIR      <- "figures"

for (d in c(DATA_DIR, RESULTS_DIR, FIG_DIR,
            file.path(DATA_DIR, "GSE135251_counts"),
            file.path(DATA_DIR, "derived"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Helper to build a results path, creating the sub-directory on demand.
res_path <- function(...) {
  p <- file.path(RESULTS_DIR, ...)
  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  p
}

fig_path <- function(...) {
  p <- file.path(FIG_DIR, ...)
  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  p
}

# ---- Globals ---------------------------------------------------------------
SEED <- 123
set.seed(SEED)
options(stringsAsFactors = FALSE, timeout = 3600)

# ---- Gene sets -------------------------------------------------------------

# 11-gene tissue signature (LASSO n RF n SVM-RFE n XGBoost on 38 candidates)
SIG11 <- c("EYA2", "FERMT1", "PYROXD2", "EPB41L4A", "CHI3L1",
           "CNKSR1", "ILDR1", "GSN", "TMPRSS3", "ESRRG", "SULF2")

# 16-gene HGF+ HSC signature: genes selected by paired FOV pseudobulk in both
# CosMx slides, concordant direction, HGF itself excluded (circularity).
SIG16 <- c("MEG3", "IGFBP3", "ITGA9", "LIFR", "CXCL12", "FGFR2", "BMP1",
           "COL27A1", "CXCL14", "RAMP1", "RXRB", "ITGA1", "RXRA", "RSPO3",
           "DUSP6", "NCAM1")

# Five consensus transcription factors (both slides significant, same direction)
TF_CONSENSUS <- c("TFAP2A", "STAT1", "TP53", "AR", "SP3")

# IL6-family ligands present on the CosMx 1,000-plex panel
IL6_LIGANDS   <- c("IL6", "IL11", "LIF", "OSM", "CTF1", "CNTF", "CLCF1")
IL6_RECEPTORS <- c("IL6R", "IL11RA", "LIFR", "OSMR", "IL6ST", "CNTFR")

# JAK/STAT and NF-kB core genes used for pathway scoring.
# NF-kB is restricted to its three panel-available core genes. Adding the
# upstream cytokines TNF/IL1B/IL6 inflates the score in inflamed hepatocytes and
# makes the pathway difference survive depth correction only under that wider
# definition, so the core set is used throughout.
JAK_STAT_GENES <- c("JAK1", "JAK2", "TYK2", "STAT1", "STAT3", "STAT5A", "STAT5B")
NFKB_GENES     <- c("NFKB1", "RELA", "NFKBIA")

# HSC-independent fibrotic-region signature (avoids circularity with HSCs)
FIBRO_SIG <- c("IL6", "IL1B", "CCL2", "CXCL9", "CXCL10", "TNF", "MMP2", "MMP14")

# Collagen-based fibrotic-region signature (sensitivity analysis only)
FIBRO_SIG_COLLAGEN <- c("COL1A1", "COL1A2", "COL3A1", "COL4A1", "FN1")

# Zonation marker panels used for the Hep.4 identity check (scRNA reference)
ZONE1_GENES <- c("PCK1", "SDS", "CPS1", "ASS1", "G6PC", "FBP1", "HAL")
ZONE2_GENES <- c("HNF4A", "ALB", "APOA1", "APOB", "TTR")
ZONE3_GENES <- c("GLUL", "CYP1A2", "OAT", "CYP2E1", "AXIN2")

# Negative-control genes (not expected in liver) for ambient-RNA assessment
NEG_CONTROL_GENES <- c("KRT5", "KRT14", "INS", "GCG", "SNAP25", "MYH2")

# ---- Spatial / model constants --------------------------------------------
CELLTYPE_COL_S1  <- "RNA_nbclust_8b0a3ae1.8754.4ba9.8269.45e2cf2e855a_1_clusters"
CELLTYPE_COL_S2  <- "RNA_nbclust_ee1320c5.9f48.455c.b54e.1fbfd34347ab_1_clusters"
PATIENT_X_CUTOFF <- 5        # mm; FOVs with x-centre < cutoff belong to patient A
NICHE_RADIUS_UM  <- 50       # HGF+ HSC neighbourhood radius
NEAR_UM          <- 50       # distance-gradient "near" threshold
FAR_UM           <- 150      # distance-gradient "far" threshold
COSMX_UM_PER_PX  <- 0.18     # CellChat spatial conversion factor

# ---- Display ---------------------------------------------------------------
COL_HGFPOS  <- "#C73E1D"
COL_HGFNEG  <- "#2E86AB"
COL_ACCENT  <- "#F18F01"
COL_NEUTRAL <- "grey70"

# ---- Reported values (used only for the consistency tables) ---------------
# These are the numbers that appear in the manuscript. The corresponding
# scripts print an observed-vs-reported comparison rather than hard-coding them
# into any statistical test.
REPORTED <- list(
  n_training          = 216,
  n_training_nafld    = 206,
  n_low_fibrosis      = 94,
  n_high_fibrosis     = 122,
  n_deg               = 986,
  n_wgcna_modules     = 17,
  r_ME13              = 0.484,
  n_candidates        = 38,
  auc_cv_fixed        = 0.841,
  auc_cv_nested       = 0.844,
  auc_GSE49541        = 0.784,
  auc_GSE130970       = 0.661,
  auc_CHI3L1          = 0.802,
  p_delong_CHI3L1     = 0.654,
  stellate_pct        = 6.5,
  OR_stellate         = 1.49,
  hgfpos_pct_slide1   = 41.1,
  hgfpos_pct_slide2   = 37.6,
  r_cross_slide       = 0.897,
  spi_hgfpos          = -0.353,
  spi_hgfneg          = -1.033,
  fov_consistent      = "97/122",
  n_hepg4             = 3028,
  n_hepg4_zonation    = 2763,
  d_central_um        = 30.7,
  d_portal_um         = 126.7,
  r_central           = 1.02,
  r_portal            = 2.21,
  zonation_log2       = -1.13,
  zonation_P          = 1.13e-4,
  n_il6r_pos          = 1921,
  n_il6r_neg          = 1107,
  n_il6r_il11ra       = 814,
  d_il6r_pos          = 49.3,
  d_il6r_neg          = 56.7,
  p_il6r_wilcox       = 1.8e-7,
  # primary estimate is the patient-level random-effects model (S31);
  # the FOV-level and cluster-robust models are sensitivity analyses
  beta_il6r_patient   = -0.030,
  p_il6r_patient      = 0.429,
  jak_pos             = 0.323,
  jak_neg             = 0.000,
  nfkb_pos            = 0.490,
  nfkb_neg            = 0.083,
  beta_jak_depth      = 0.091,
  p_jak_depth         = 0.007,
  beta_nfkb_depth     = 0.068,
  p_nfkb_depth        = 0.058
)
