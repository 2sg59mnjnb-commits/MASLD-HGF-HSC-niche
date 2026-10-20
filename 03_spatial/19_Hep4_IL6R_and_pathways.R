# ============================================================================
# 19_Hep4_IL6R_and_pathways.R
#
# Hep.4 hepatocytes stratified by IL6R expression:
#   - proportions (IL6R+, IL6R-, IL6R+IL11RA+)
#   - distance to the nearest HGF+ HSC, under several models
#   - JAK/STAT and NF-kB pathway activity, unadjusted and depth-corrected
#   - STAT3 / IL6R correlation within Hep.4 hepatocytes
#
# The proximity result is deliberately reported as model-dependent: it is
# significant at the cell level and when FOV is nested in patient, but not in
# the patient-level random-effects model, which is the conservative choice.
#
# Output
#   results/03_spatial/hep4_il6r_stratification.csv
#   results/03_spatial/hep4_il6r_distance_models.csv
#   results/03_spatial/hep4_pathways.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({
  library(lme4); library(lmerTest); library(sandwich); library(lmtest)
})

df <- load_derived("cosmx_flat.rds")
hep4 <- df[df$celltype == "Hep.4", ]
hsc_pos <- df[df$celltype == "Stellate.cells" & df$HGF > 0, ]

hep4$IL6R_pos   <- as.numeric(hep4$IL6R > 0)
hep4$double_pos <- as.numeric(hep4$IL6R > 0 & hep4$IL11RA > 0)

# ---- 1. Proportions ------------------------------------------------------
n <- nrow(hep4)
strat <- data.frame(
  group = c("IL6R+", "IL6R-", "IL6R+IL11RA+ (double)"),
  n = c(sum(hep4$IL6R_pos == 1), sum(hep4$IL6R_pos == 0), sum(hep4$double_pos == 1))
)
strat$pct <- round(100 * strat$n / n, 1)
print(strat)
write.csv(strat, res_path("03_spatial", "hep4_il6r_stratification.csv"),
          row.names = FALSE)
cat(sprintf("\nHep.4 total: %d; HGF+ HSC: %d\n", n, nrow(hsc_pos)))

# ---- 2. Distance to nearest HGF+ HSC ------------------------------------
hep4$d_hgfpos <- compute_nn(hep4, hsc_pos)
hep4$dist_scaled <- as.numeric(scale(hep4$d_hgfpos))
d <- hep4[!is.na(hep4$d_hgfpos), ]

cat(sprintf("\nMedian distance to nearest HGF+ HSC: IL6R+ %.1f um, IL6R- %.1f um\n",
            median(d$d_hgfpos[d$IL6R_pos == 1]),
            median(d$d_hgfpos[d$IL6R_pos == 0])))

w <- wilcox.test(d_hgfpos ~ IL6R_pos, data = d)
cat(sprintf("Wilcoxon P = %.2e\n", w$p.value))

# Permutation test on the median difference
set.seed(SEED)
obs_delta <- median(d$d_hgfpos[d$IL6R_pos == 1]) - median(d$d_hgfpos[d$IL6R_pos == 0])
null <- replicate(1000, {
  s <- sample(d$IL6R_pos)
  median(d$d_hgfpos[s == 1]) - median(d$d_hgfpos[s == 0])
})
p_perm <- (sum(abs(null) >= abs(obs_delta)) + 1) / 1001
cat(sprintf("Median difference %.1f um; permutation P = %.3f\n", obs_delta, p_perm))

# ---- 3. Models -----------------------------------------------------------
models <- rbind(
  lmm_effect(dist_scaled ~ IL6R_pos + (1 | patient_uid), d,
             "IL6R_pos", "Patient random intercept"),
  lmm_effect(dist_scaled ~ IL6R_pos + (1 | fov), d,
             "IL6R_pos", "FOV random intercept"),
  lmm_effect(dist_scaled ~ IL6R_pos + (1 | patient_uid/fov), d,
             "IL6R_pos", "FOV nested in patient")
)
lm_fit <- lm(dist_scaled ~ IL6R_pos, data = d)
cl <- coeftest(lm_fit, vcov = vcovCL(lm_fit, cluster = d$fov))
models <- rbind(models, data.frame(
  Model = "FOV cluster-robust SE", Term = "IL6R_pos",
  Beta = round(cl["IL6R_pos", "Estimate"], 4),
  CI_low = NA, CI_high = NA,
  t_value = round(cl["IL6R_pos", "t value"], 3), df = NA,
  P = cl["IL6R_pos", "Pr(>|t|)"]
))
models <- rbind(models, data.frame(
  Model = "Naive cell-level t-test", Term = "IL6R_pos",
  Beta = NA, CI_low = NA, CI_high = NA, t_value = NA, df = NA,
  P = t.test(d_hgfpos ~ IL6R_pos, data = d)$p.value
))
cat("\n=== IL6R+ vs IL6R- Hep.4: distance to HGF+ HSC ===\n")
print(models)
write.csv(models, res_path("03_spatial", "hep4_il6r_distance_models.csv"),
          row.names = FALSE)

pat_re <- models[models$Model == "Patient random intercept", ]
compare("IL6R patient-RE P", pat_re$P, REPORTED$p_il6r_patient)
cat("\nThe patient-level model is reported as the primary result; the nested-FOV\n")
cat("and cluster-robust models are sensitivity analyses only.\n")

# ---- 4. Pathway activity -------------------------------------------------
# Pathway scores are z-scored gene-wise across ALL cells of a slide and then
# averaged within Hep.4. Scoring within Hep.4 alone would centre every gene on
# the Hep.4 mean and remove the between-compartment contrast that the score is
# meant to capture.
score_pathway <- function(genes) {
  out <- rep(NA_real_, nrow(hep4))
  for (sl in c(1, 2)) {
    obj <- load_cosmx(sl)
    dat <- obj[["RNA"]]$data
    g <- intersect(genes, rownames(dat))
    idx <- which(hep4$slide == paste0("S", sl))
    if (length(g) == 0) { rm(obj); next }
    m <- as.matrix(dat[g, , drop = FALSE])          # all cells
    z <- t(scale(t(m))); z[is.na(z)] <- 0
    sc <- colMeans(z)
    out[idx] <- sc[match(hep4$cell[idx], colnames(dat))]
    rm(obj, dat, m, z); gc(verbose = FALSE)
  }
  out
}

hep4$jak_stat <- score_pathway(JAK_STAT_GENES)
hep4$nfkb     <- score_pathway(NFKB_GENES)

cat("\n=== Pathway activity (z-score within Hep.4) ===\n")
path_tab <- do.call(rbind, lapply(c("jak_stat", "nfkb"), function(v) {
  a <- hep4[[v]][hep4$IL6R_pos == 1]; b <- hep4[[v]][hep4$IL6R_pos == 0]
  data.frame(pathway = v,
             mean_IL6Rpos = round(mean(a, na.rm = TRUE), 3),
             mean_IL6Rneg = round(mean(b, na.rm = TRUE), 3),
             P_wilcox = wilcox.test(a, b)$p.value)
}))
print(path_tab)

# Depth-corrected: nCount as covariate, and trimming above the 95th percentile
hep4$nCount_scaled <- as.numeric(scale(hep4$nCount))
depth_tab <- do.call(rbind, lapply(c("jak_stat", "nfkb"), function(v) {
  d2 <- hep4[!is.na(hep4[[v]]), ]
  d2$y <- as.numeric(scale(d2[[v]]))
  m <- lmer(y ~ IL6R_pos + nCount_scaled + (1 | patient_uid),
            data = d2, REML = FALSE)
  ct <- summary(m)$coefficients

  cutoff <- quantile(d2$nCount, 0.95)
  dt <- d2[d2$nCount <= cutoff, ]
  dt$y <- as.numeric(scale(dt[[v]]))
  mt <- lmer(y ~ IL6R_pos + (1 | patient_uid), data = dt, REML = FALSE)
  ctt <- summary(mt)$coefficients

  data.frame(
    pathway = v,
    beta_depth_covariate = round(ct["IL6R_pos", 1], 3),
    P_depth_covariate = ct["IL6R_pos", "Pr(>|t|)"],
    beta_trim95 = round(ctt["IL6R_pos", 1], 3),
    P_trim95 = ctt["IL6R_pos", "Pr(>|t|)"],
    n_trimmed = nrow(dt)
  )
}))
print(depth_tab)
write.csv(rbind(path_tab, depth_tab),
          res_path("03_spatial", "hep4_pathways.csv"), row.names = FALSE)

compare("JAK/STAT beta (depth covariate)",
        depth_tab$beta_depth_covariate[depth_tab$pathway == "jak_stat"], 0.088)
compare("NF-kB beta (depth covariate)",
        depth_tab$beta_depth_covariate[depth_tab$pathway == "nfkb"], 0.095)

# ---- 5. STAT3 / IL6R correlation ----------------------------------------
if (!all(is.na(hep4$STAT3))) {
  ct <- cor.test(hep4$STAT3, hep4$IL6R, method = "spearman")
  cat(sprintf("\nSTAT3 vs IL6R within Hep.4: rho = %.3f (P = %.2e)\n",
              ct$estimate, ct$p.value))
}

saveRDS(list(hep4 = hep4, stratification = strat, models = models,
             pathways = path_tab, depth = depth_tab,
             wilcox_P = w$p.value, perm_P = p_perm),
        file.path(DATA_DIR, "derived", "hep4_il6r.rds"))

cat("\nNext: 03_spatial/20_depth_correction.R\n")
