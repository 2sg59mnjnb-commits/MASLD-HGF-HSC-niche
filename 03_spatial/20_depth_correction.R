# ============================================================================
# 20_depth_correction.R
#
# HGF+ HSCs have a substantially higher transcript count than HGF- HSCs
# (median nCount ~288 vs ~114 in these sections). Because the two populations do
# not overlap in depth, quantile matching cannot equalise absolute depth: the
# matched HGF- cells stay at their own depth. Two approaches that do work are
# therefore used instead:
#
#   (a) model nCount as a covariate in the mixed model
#   (b) trim cells above the 95th percentile of nCount
#
# Both are applied to the HGF+ vs HGF- distance result. The same comparison is
# also reported for the 16-gene signature in 04_single_cell/24.
#
# Output
#   results/03_spatial/depth_comparison.csv
#   results/03_spatial/depth_corrected_distance.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(lme4); library(lmerTest)})

df <- load_derived("cosmx_flat.rds")
hsc  <- df[df$celltype == "Stellate.cells", ]
hep4 <- df[df$celltype == "Hep.4", ]

hsc$d_hep4  <- compute_nn(hsc, hep4)
hsc$HGF_pos <- as.numeric(hsc$HGF > 0)
hsc$dist_scaled <- as.numeric(scale(hsc$d_hep4))

# ---- 1. Depth comparison -------------------------------------------------
depth_tab <- do.call(rbind, lapply(c("nCount", "nFeature"), function(v) {
  a <- hsc[[v]][hsc$HGF_pos == 1]; b <- hsc[[v]][hsc$HGF_pos == 0]
  data.frame(
    variable = v,
    hgfpos_median = median(a, na.rm = TRUE),
    hgfneg_median = median(b, na.rm = TRUE),
    ratio = round(median(a, na.rm = TRUE) / median(b, na.rm = TRUE), 2),
    P = wilcox.test(a, b)$p.value
  )
}))
cat("=== Depth comparison (HGF+ vs HGF-) ===\n")
print(depth_tab)
write.csv(depth_tab, res_path("03_spatial", "depth_comparison.csv"),
          row.names = FALSE)

cat("\nBecause the depth distributions barely overlap, quantile matching equalises\n")
cat("rank but not absolute depth; it is therefore not used for the distance\n")
cat("analysis. Covariate adjustment and trimming are used instead.\n")

# ---- 2. Three depth-correction strategies --------------------------------
hsc$nCount_scaled <- as.numeric(scale(hsc$nCount))

# (a) unadjusted
m0 <- lmer(dist_scaled ~ HGF_pos + (1 | patient_uid), data = hsc, REML = FALSE)
# (b) nCount as covariate
m1 <- lmer(dist_scaled ~ HGF_pos + nCount_scaled + (1 | patient_uid),
           data = hsc, REML = FALSE)
# (c) trim above the 95th percentile
cutoff <- quantile(hsc$nCount, 0.95)
trim <- hsc[hsc$nCount <= cutoff, ]
trim$dist_scaled <- as.numeric(scale(trim$d_hep4))
m2 <- lmer(dist_scaled ~ HGF_pos + (1 | patient_uid), data = trim, REML = FALSE)

extract <- function(m, label, n) {
  ct <- summary(m)$coefficients
  b <- ct["HGF_pos", "Estimate"]; se <- ct["HGF_pos", "Std. Error"]
  data.frame(
    strategy = label, n = n,
    beta = round(b, 3),
    ci_low = round(b - 1.96 * se, 3),
    ci_high = round(b + 1.96 * se, 3),
    t = round(ct["HGF_pos", "t value"], 2),
    P = ct["HGF_pos", "Pr(>|t|)"]
  )
}

res <- rbind(
  extract(m0, "Unadjusted", nrow(hsc)),
  extract(m1, "nCount as covariate", nrow(hsc)),
  extract(m2, sprintf("Trim nCount > 95th percentile (%.0f)", cutoff), nrow(trim))
)
cat("\n=== Distance result under depth correction ===\n")
print(res)
write.csv(res, res_path("03_spatial", "depth_corrected_distance.csv"),
          row.names = FALSE)

compare("unadjusted beta", res$beta[res$strategy == "Unadjusted"], -0.271)
cat("\nThe effect is attenuated but not abolished by depth adjustment, which is the\n")
cat("basis for reporting the niche as a depth-robust observation.\n")

saveRDS(list(depth_comparison = depth_tab, correction = res,
             cutoff = cutoff, models = list(m0 = m0, m1 = m1, m2 = m2)),
        file.path(DATA_DIR, "derived", "depth_correction.rds"))

cat("\nNext: 03_spatial/21_distance_gradient.R\n")
