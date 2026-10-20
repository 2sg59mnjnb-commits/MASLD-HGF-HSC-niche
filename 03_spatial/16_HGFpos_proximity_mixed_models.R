# ============================================================================
# 16_HGFpos_proximity_mixed_models.R
#
# The central spatial result: are HGF+ HSCs closer to Hep.4 hepatocytes than
# HGF- HSCs, and does that survive modelling of within-patient pseudoreplication?
#
# Analyses
#   1. nearest Hep.4 distance per HSC, HGF+ vs HGF-
#   2. Selective Proximity Index: log2(d_Hep5 / d_Hep4), which normalises for
#      local hepatocyte density using a same-compartment reference subtype
#   3. linear mixed-effects models with patient, FOV, or FOV-within-patient as
#      the random intercept, plus FOV cluster-robust standard errors
#   4. ICC and a within-patient permutation test
#   5. per-patient median differences
#
# Output
#   results/03_spatial/hgfpos_proximity_models.csv
#   results/03_spatial/per_patient_delta.csv
#   data/derived/proximity_models.rds
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({
  library(lme4); library(lmerTest); library(sandwich); library(lmtest)
})

df <- load_derived("cosmx_flat.rds")
hsc  <- df[df$celltype == "Stellate.cells", ]
hep4 <- df[df$celltype == "Hep.4", ]
hep5 <- df[df$celltype == "Hep.5", ]

hsc$d_hep4 <- compute_nn(hsc, hep4)
hsc$d_hep5 <- compute_nn(hsc, hep5)
hsc$HGF_pos <- as.numeric(hsc$HGF > 0)
hsc$dist_scaled <- as.numeric(scale(hsc$d_hep4))
hsc$spi <- log2((hsc$d_hep5 + 1) / (hsc$d_hep4 + 1))

cat(sprintf("HSC: %d | Hep.4: %d | Hep.5: %d\n", nrow(hsc), nrow(hep4), nrow(hep5)))
cat(sprintf("HGF+ %d, HGF- %d\n", sum(hsc$HGF_pos == 1), sum(hsc$HGF_pos == 0)))

# ---- 1. Unadjusted distance ----------------------------------------------
w <- wilcox.test(d_hep4 ~ HGF_pos, data = hsc)
cat(sprintf("\nMedian distance to Hep.4: HGF+ %.1f um, HGF- %.1f um\n",
            median(hsc$d_hep4[hsc$HGF_pos == 1], na.rm = TRUE),
            median(hsc$d_hep4[hsc$HGF_pos == 0], na.rm = TRUE)))
cat(sprintf("Wilcoxon P = %.2e\n", w$p.value))

# ---- 2. Selective Proximity Index ----------------------------------------
spi_pos <- hsc$spi[hsc$HGF_pos == 1]; spi_neg <- hsc$spi[hsc$HGF_pos == 0]
w_spi <- wilcox.test(spi_pos, spi_neg)
cat(sprintf("\nSPI (log2 d_Hep5/d_Hep4): HGF+ median %.3f, HGF- median %.3f\n",
            median(spi_pos, na.rm = TRUE), median(spi_neg, na.rm = TRUE)))
cat(sprintf("SPI Wilcoxon P = %.2e\n", w_spi$p.value))
compare("SPI HGF+", median(spi_pos, na.rm = TRUE), REPORTED$spi_hgfpos)
compare("SPI HGF-", median(spi_neg, na.rm = TRUE), REPORTED$spi_hgfneg)

# SPI for other HSC subsets, as a specificity control
if ("MET" %in% colnames(hsc)) {
  met_pos <- as.numeric(hsc$MET > 0)
  if (sum(met_pos) > 10) {
    cat(sprintf("\nSPI (MET+ HSC, n = %d): median %.3f\n",
                sum(met_pos), median(hsc$spi[met_pos == 1], na.rm = TRUE)))
  }
}

# ---- 3. Mixed models -----------------------------------------------------
models <- rbind(
  lmm_effect(dist_scaled ~ HGF_pos + (1 | patient_uid), hsc,
             "HGF_pos", "Patient random intercept"),
  lmm_effect(dist_scaled ~ HGF_pos + (1 | fov), hsc,
             "HGF_pos", "FOV random intercept"),
  lmm_effect(dist_scaled ~ HGF_pos + (1 | patient_uid/fov), hsc,
             "HGF_pos", "FOV nested in patient")
)

lm_fit <- lm(dist_scaled ~ HGF_pos, data = hsc)
cl <- coeftest(lm_fit, vcov = vcovCL(lm_fit, cluster = hsc$fov))
models <- rbind(models, data.frame(
  Model = "FOV cluster-robust SE", Term = "HGF_pos",
  Beta = round(cl["HGF_pos", "Estimate"], 4),
  CI_low = NA, CI_high = NA,
  t_value = round(cl["HGF_pos", "t value"], 3), df = NA,
  P = cl["HGF_pos", "Pr(>|t|)"]
))

naive <- t.test(d_hep4 ~ HGF_pos, data = hsc)
models <- rbind(models, data.frame(
  Model = "Naive cell-level t-test (biased)", Term = "HGF_pos",
  Beta = NA, CI_low = NA, CI_high = NA, t_value = NA, df = NA,
  P = naive$p.value
))

cat("\n=== Mixed-effects models ===\n")
print(models)
write.csv(models, res_path("03_spatial", "hgfpos_proximity_models.csv"),
          row.names = FALSE)

compare("patient-RE beta",
        models$Beta[models$Model == "Patient random intercept"],
        -0.271)

# ---- 4. ICC and permutation test ----------------------------------------
icc <- icc_lmm(dist_scaled ~ HGF_pos + (1 | patient_uid), hsc)
cat(sprintf("\nICC (patient) = %.4f  [between-patient variance %.4f, residual %.4f]\n",
            icc["icc"], icc["var_random"], icc["var_resid"]))

cat("Within-patient label permutation (500 iterations) ...\n")
set.seed(SEED)
obs_beta <- fixef(lmer(dist_scaled ~ HGF_pos + (1 | patient_uid),
                       data = hsc, REML = FALSE))["HGF_pos"]
groups <- split(seq_len(nrow(hsc)), hsc$patient_uid)
n_perm <- 500
null_beta <- numeric(n_perm)
for (i in seq_len(n_perm)) {
  d <- hsc
  for (g in names(groups)) {
    idx <- groups[[g]]
    d$HGF_pos[idx] <- sample(d$HGF_pos[idx], replace = FALSE)
  }
  m <- tryCatch(lmer(dist_scaled ~ HGF_pos + (1 | patient_uid),
                     data = d, REML = FALSE), error = function(e) NULL)
  null_beta[i] <- if (is.null(m)) NA else fixef(m)["HGF_pos"]
}
null_beta <- null_beta[!is.na(null_beta)]
p_perm <- (sum(abs(null_beta) >= abs(obs_beta)) + 1) / (length(null_beta) + 1)
cat(sprintf("Observed beta = %.3f; permutation P = %.4f (%d valid permutations)\n",
            obs_beta, p_perm, length(null_beta)))

# ---- 5. Per-patient consistency -----------------------------------------
per_patient <- hsc %>%
  dplyr::group_by(patient_uid, HGF_pos) %>%
  dplyr::summarise(median_d = median(d_hep4, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = HGF_pos, values_from = median_d,
                     names_prefix = "hgf") %>%
  dplyr::mutate(delta = hgf1 - hgf0)
print(as.data.frame(per_patient))
write.csv(per_patient, res_path("03_spatial", "per_patient_delta.csv"),
          row.names = FALSE)

n_neg <- sum(per_patient$delta < 0, na.rm = TRUE)
bin_p <- binom.test(n_neg, nrow(per_patient), 0.5, alternative = "greater")$p.value
cat(sprintf("\nHGF+ closer in %d / %d patients (binomial P = %.4f)\n",
            n_neg, nrow(per_patient), bin_p))
cat(sprintf("Delta range: %.1f to %.1f um\n",
            min(per_patient$delta, na.rm = TRUE),
            max(per_patient$delta, na.rm = TRUE)))
cat("With n = 4 the smallest attainable one-sided binomial P is 0.0625, so the\n")
cat("per-patient consistency is the informative quantity, not its P value.\n")

saveRDS(list(hsc = hsc, models = models, icc = icc,
             permutation = list(obs = obs_beta, null = null_beta, P = p_perm),
             per_patient = per_patient, spi_test = w_spi, dist_test = w),
        file.path(DATA_DIR, "derived", "proximity_models.rds"))

cat("\nNext: 03_spatial/17_FOV_consistency_and_thresholds.R\n")

# ---- 6. SPI for all HSC subsets (Figure 2f) ------------------------------
# PDGFRB is not carried in the flattened object, so it is read from each slide
# and matched back by cell name.
cat("\n=== Selective Proximity Index by HSC subset ===\n")
pdgfrb <- rep(NA_real_, nrow(hsc))
for (s in c(1, 2)) {
  obj <- load_cosmx(s)
  dat <- obj[["RNA"]]$data
  if (!"PDGFRB" %in% rownames(dat)) { rm(obj); next }
  idx <- which(hsc$slide == paste0("S", s))
  pdgfrb[idx] <- as.numeric(dat["PDGFRB", hsc$cell[idx]])
  rm(obj); gc(verbose = FALSE)
}
hsc$PDGFRB_pos <- as.numeric(pdgfrb > 0)

spi_summary <- data.frame(
  subset = c("HGF+", "HGF-", "PDGFRB+", "PDGFRB-"),
  n = c(sum(hsc$HGF_pos == 1), sum(hsc$HGF_pos == 0),
        sum(hsc$PDGFRB_pos == 1, na.rm = TRUE),
        sum(hsc$PDGFRB_pos == 0, na.rm = TRUE)),
  median_spi = c(median(hsc$spi[hsc$HGF_pos == 1], na.rm = TRUE),
                 median(hsc$spi[hsc$HGF_pos == 0], na.rm = TRUE),
                 median(hsc$spi[hsc$PDGFRB_pos == 1], na.rm = TRUE),
                 median(hsc$spi[hsc$PDGFRB_pos == 0], na.rm = TRUE))
)
print(spi_summary, row.names = FALSE)
write.csv(spi_summary, res_path("03_spatial", "spi_by_subset.csv"), row.names = FALSE)
saveRDS(spi_summary, file.path(DATA_DIR, "derived", "spi_subsets.rds"))
