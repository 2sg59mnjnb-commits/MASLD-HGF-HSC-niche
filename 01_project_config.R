
# Central project configuration
PROJ_ROOT <- normalizePath(".")
DATA_DIR   <- file.path(PROJ_ROOT, "data")
RESULT_DIR <- file.path(PROJ_ROOT, "results")
FIG_DIR    <- file.path(PROJ_ROOT, "figures")

for (d in c(DATA_DIR, RESULT_DIR, FIG_DIR,
            file.path(RESULT_DIR, "supplementary"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

SEED <- 123
N_CV_FOLDS <- 10
N_BOOTSTRAP <- 2000
N_PERMUTATION <- 1000

COSMX_SLIDE1_CELLTYPE_COL <-
  "RNA_nbclust_8b0a3ae1.8754.4ba9.8269.45e2cf2e855a_1_clusters"
COSMX_SLIDE2_CELLTYPE_COL <-
  "RNA_nbclust_ee1320c5.9f48.455c.b54e.1fbfd34347ab_1_clusters"

COL_PRIMARY <- "#C73E1D"
COL_SECOND  <- "#2E86AB"
COL_ACCENT  <- "#F18F01"
COL_NEUTRAL <- "grey70"

`%||%` <- function(a, b) if (is.null(a)) b else a

cat("✅ Project config loaded.\n")

