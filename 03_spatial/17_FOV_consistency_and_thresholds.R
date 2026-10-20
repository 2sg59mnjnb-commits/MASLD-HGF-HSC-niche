# ============================================================================
# 17_FOV_consistency_and_thresholds.R
#
# Two robustness checks for the HGF+ HSC proximity result:
#
#   (a) FOV-level direction consistency. Within each field of view, is the
#       median distance to Hep.4 shorter for HGF+ than for HGF- HSCs? FOVs with
#       fewer than three HSCs of either class are excluded because a median
#       cannot be estimated.
#
#   (b) HGF positivity threshold. HGF+ is redefined at HGF > 0, > 1 and > 2 and
#       the mixed-model effect is re-estimated. If the effect strengthens with a
#       stricter threshold, the signal is unlikely to come from ambient or
#       low-level transcripts.
#
# Output
#   results/03_spatial/fov_direction_consistency.csv
#   results/03_spatial/hgf_threshold_effects.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(lme4); library(lmerTest)})

df <- load_derived("cosmx_flat.rds")
hsc  <- df[df$celltype == "Stellate.cells", ]
hep4 <- df[df$celltype == "Hep.4", ]
hsc$d_hep4  <- compute_nn(hsc, hep4)

# ---- (a) FOV-level direction ---------------------------------------------
fov_res <- do.call(rbind, lapply(split(hsc, fov_key(hsc)), function(d) {
  pos <- d$d_hep4[d$HGF > 0]; neg <- d$d_hep4[d$HGF == 0]
  pos <- pos[!is.na(pos)]; neg <- neg[!is.na(neg)]
  if (length(pos) < 3 || length(neg) < 3) return(NULL)
  data.frame(
    fov = d$fov[1], slide = d$slide[1], patient = d$patient_uid[1],
    n_pos = length(pos), n_neg = length(neg),
    median_pos = median(pos), median_neg = median(neg),
    delta = median(pos) - median(neg),
    direction = ifelse(median(pos) < median(neg), "HGF+ closer", "HGF- closer"),
    stringsAsFactors = FALSE
  )
}))

n_use <- nrow(fov_res)
n_cons <- sum(fov_res$direction == "HGF+ closer")
bt <- binom.test(n_cons, n_use, 0.5, alternative = "greater")
cat(sprintf("FOVs with >= 3 HSCs of each class: %d (slide-specific fields)\n", n_use))
cat(sprintf("HGF+ closer in %d / %d (%.1f%%), binomial P = %.2e\n",
            n_cons, n_use, 100 * n_cons / n_use, bt$p.value))
compare("FOV consistent", n_cons, 48, digits = 0)
cat(sprintf("(denominator note: 64 FOVs were profiled; %d contained enough HSCs of\n", n_use))
cat(" both classes for a direction test)\n")

write.csv(fov_res, res_path("03_spatial", "fov_direction_consistency.csv"),
          row.names = FALSE)

# ---- (b) HGF threshold ----------------------------------------------------
cat("\n=== HGF positivity threshold ===\n")
thr_tab <- do.call(rbind, lapply(c(0, 1, 2), function(thr) {
  d <- hsc
  d$pos <- as.numeric(d$HGF > thr)
  if (sum(d$pos) < 50 || sum(1 - d$pos) < 50) return(NULL)
  d$dist_scaled <- as.numeric(scale(d$d_hep4))

  m <- lmer(dist_scaled ~ pos + (1 | patient_uid), data = d, REML = FALSE)
  ct <- summary(m)$coefficients
  beta <- ct["pos", "Estimate"]; se <- ct["pos", "Std. Error"]

  per_pat <- d %>%
    dplyr::group_by(patient_uid, pos) %>%
    dplyr::summarise(md = median(d_hep4, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = pos, values_from = md, names_prefix = "p") %>%
    dplyr::filter(!is.na(p1), !is.na(p0)) %>%
    dplyr::mutate(delta = p1 - p0)

  data.frame(
    threshold = paste0("HGF > ", thr),
    n_pos = sum(d$pos), n_neg = sum(1 - d$pos),
    beta = round(beta, 4),
    ci_low = round(beta - 1.96 * se, 4),
    ci_high = round(beta + 1.96 * se, 4),
    P = ct["pos", "Pr(>|t|)"],
    patients_HGFpos_closer = sprintf("%d/%d",
                                     sum(per_pat$delta < 0, na.rm = TRUE),
                                     nrow(per_pat)),
    median_delta_um = round(median(per_pat$delta, na.rm = TRUE), 1)
  )
}))
print(thr_tab)
write.csv(thr_tab, res_path("03_spatial", "hgf_threshold_effects.csv"),
          row.names = FALSE)

cat("\nA monotonic strengthening of beta with the threshold argues against an\n")
cat("artefact of low-level or ambient transcript detection.\n")

# ---- Ambient RNA check ----------------------------------------------------
cat("\n=== Ambient RNA negative controls ===\n")
for (g in intersect(c("KRT5", "KRT14", "ALB"), colnames(df))) {
  cat(sprintf("  %-6s positive in %.1f%% of all cells; %.1f%% of HSCs\n", g,
              100 * mean(df[[g]] > 0, na.rm = TRUE),
              100 * mean(hsc[[g]] > 0, na.rm = TRUE)))
}

saveRDS(list(fov = fov_res, binomial_P = bt$p.value,
             n_consistent = n_cons, n_fov = n_use, thresholds = thr_tab),
        file.path(DATA_DIR, "derived", "fov_thresholds.rds"))

cat("\nNext: 03_spatial/18_zonation_portal_vs_central.R\n")

# ---- (c) Distance-threshold sensitivity ----------------------------------
# Truncating the HSC population by distance to Hep.4 and re-comparing the two
# groups shows whether the direction of effect holds across spatial scales. This
# is the panel plotted as Figure 2d.
cat("\n=== Distance-threshold sensitivity ===\n")
dist_tab <- do.call(rbind, lapply(c(20, 30, 40, 50, 75, 100), function(cut) {
  sub <- hsc[hsc$d_hep4 <= cut, ]
  if (sum(sub$HGF > 0) < 20 || sum(sub$HGF == 0) < 20) return(NULL)
  data.frame(
    threshold_um = cut,
    n_hgfpos = sum(sub$HGF > 0),
    n_hgfneg = sum(sub$HGF == 0),
    median_hgfpos = median(sub$d_hep4[sub$HGF > 0]),
    median_hgfneg = median(sub$d_hep4[sub$HGF == 0]),
    P = wilcox.test(d_hep4 ~ I(HGF > 0), data = sub)$p.value
  )
}))
dist_tab$delta_um <- dist_tab$median_hgfpos - dist_tab$median_hgfneg
print(dist_tab, row.names = FALSE)
write.csv(dist_tab, res_path("03_spatial", "distance_threshold_sensitivity.csv"),
          row.names = FALSE)

saveRDS(list(distance_thresholds = dist_tab),
        file.path(DATA_DIR, "derived", "distance_thresholds.rds"))
