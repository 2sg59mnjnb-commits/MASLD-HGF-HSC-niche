# ============================================================================
# 18_zonation_portal_vs_central.R
#
# Places the niche anatomically: are Hep.4 hepatocytes nearer central-venous or
# nearer portal structures?
#
# Reference populations
#   portal  = Periportal.LSECs + Cholangiocytes
#   central = Central.venous.LSECs
#
# Nearest-neighbour distances are confounded by the abundance of the reference
# population: central-venous endothelial cells are roughly twice as common in
# these sections, which shortens nearest-neighbour distances mechanically. Two
# independent corrections are therefore applied:
#
#   (1) Poisson expectation normalisation. For a random (Poisson) point process
#       of density lambda the expected nearest-neighbour distance is
#       0.5 / sqrt(lambda). Each observed distance is expressed as a ratio to
#       the expectation for that FOV, so the ratio is abundance-free.
#   (2) Same-FOV permutation. Reference sets of equal size are drawn at random
#       from all cells in the FOV and the distance ratio is recomputed.
#
# Output
#   results/03_spatial/zonation_hepg4.csv
#   results/03_spatial/zonation_poisson.csv
#   results/03_spatial/zonation_permutation.csv
#   results/03_spatial/zonation_hsc.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(lme4); library(lmerTest)})

df <- load_derived("cosmx_flat.rds")

PORTAL  <- c("Periportal.LSECs", "Cholangiocytes")
CENTRAL <- "Central.venous.LSECs"

df$zone <- ifelse(df$celltype %in% PORTAL, "portal",
                  ifelse(df$celltype %in% CENTRAL, "central", "other"))
cat("Reference cells:\n"); print(table(df$zone))
if (sum(df$zone == "portal") == 0 || sum(df$zone == "central") == 0) {
  stop("Portal or central reference population missing; check cell-type labels.")
}

# ---- 1. Hep.4 distances --------------------------------------------------
hep4 <- df[df$celltype == "Hep.4", ]
tgt_portal  <- df[df$zone == "portal", ]
tgt_central <- df[df$zone == "central", ]

hep4$d_portal  <- compute_nn(hep4, tgt_portal)
hep4$d_central <- compute_nn(hep4, tgt_central)
hep4 <- hep4[!is.na(hep4$d_portal) & !is.na(hep4$d_central), ]

cat(sprintf("\nHep.4 cells with both distances: %d\n", nrow(hep4)))
cat(sprintf("Median distance to central-venous: %.1f um\n", median(hep4$d_central)))
cat(sprintf("Median distance to portal:         %.1f um\n", median(hep4$d_portal)))
compare("d central", median(hep4$d_central), REPORTED$d_central_um, digits = 1)
compare("d portal",  median(hep4$d_portal),  REPORTED$d_portal_um,  digits = 1)

w <- wilcox.test(hep4$d_portal, hep4$d_central, paired = TRUE)
cat(sprintf("Paired Wilcoxon P = %.2e\n", w$p.value))

# log2(central/portal) > 0 means closer to the portal side
hep4$log_ratio <- log2((hep4$d_central + 1) / (hep4$d_portal + 1))
m_raw <- lmer(log_ratio ~ 1 + (1 | patient_uid), data = hep4, REML = FALSE)
ct <- summary(m_raw)$coefficients
cat(sprintf("Uncorrected log2(central/portal) = %.3f (t = %.1f, P = %.2e)\n",
            ct[1, 1], ct[1, 4], ct[1, 5]))

write.csv(hep4[, c("fov", "patient_uid", "d_portal", "d_central", "log_ratio")],
          res_path("03_spatial", "zonation_hepg4.csv"), row.names = FALSE)

# ---- 2. Poisson density normalisation -----------------------------------
area <- fov_area(df)
area <- area[area$area > 1e-6, ]

dens <- df %>%
  dplyr::group_by(slide, fov) %>%
  dplyr::summarise(n_portal = sum(zone == "portal"),
                   n_central = sum(zone == "central"), .groups = "drop") %>%
  dplyr::left_join(area[, c("slide", "fov", "area")], by = c("slide", "fov")) %>%
  dplyr::filter(!is.na(area), area > 1e-6) %>%
  dplyr::mutate(lambda_portal = n_portal / area,
                lambda_central = n_central / area,
                density_ratio = lambda_central / lambda_portal)
cat(sprintf("\nFOVs with density estimates: %d\n", nrow(dens)))
cat(sprintf("Median density ratio central/portal = %.2f\n",
            median(dens$density_ratio, na.rm = TRUE)))

# Expected NN distance in um for a Poisson process of density lambda (cells/mm^2)
expected_nn <- function(lambda_per_mm2) 0.5 / sqrt(lambda_per_mm2) * 1000

hep4 <- merge(hep4, dens[, c("slide", "fov", "lambda_portal", "lambda_central")],
              by = c("slide", "fov"), all.x = TRUE)
hep4$e_portal  <- expected_nn(hep4$lambda_portal)
hep4$e_central <- expected_nn(hep4$lambda_central)
hep4$r_portal  <- hep4$d_portal  / hep4$e_portal
hep4$r_central <- hep4$d_central / hep4$e_central
hep4$log_norm_ratio <- log2(hep4$r_central / hep4$r_portal)

cat(sprintf("Density-normalised r: portal %.3f, central %.3f\n",
            median(hep4$r_portal, na.rm = TRUE),
            median(hep4$r_central, na.rm = TRUE)))

w_norm <- wilcox.test(hep4$r_portal, hep4$r_central, paired = TRUE)
m_norm <- lmer(log_norm_ratio ~ 1 + (1 | patient_uid), data = hep4, REML = FALSE)
ctn <- summary(m_norm)$coefficients
cat(sprintf("Normalised log2(r_central/r_portal) = %.3f (t = %.2f, P = %.2e)\n",
            ctn[1, 1], ctn[1, 4], ctn[1, 5]))
cat(sprintf("Paired Wilcoxon on normalised distances: P = %.2e\n", w_norm$p.value))
cat(sprintf("=> after density correction Hep.4 is %s\n",
            ifelse(ctn[1, 1] > 0, "closer to PORTAL structures", "closer to CENTRAL structures")))

write.csv(hep4[, c("fov", "patient_uid", "d_portal", "d_central",
                   "r_portal", "r_central", "log_norm_ratio")],
          res_path("03_spatial", "zonation_poisson.csv"), row.names = FALSE)

# ---- 3. Same-FOV permutation --------------------------------------------
cat("\nSame-FOV permutation (500 iterations per FOV) ...\n")
set.seed(SEED)
perm <- do.call(rbind, lapply(split(hep4, fov_key(hep4)), function(d) {
  f <- d$fov[1]
  sld <- d$slide[1]
  pool <- df[df$slide == sld & df$fov == f, ]
  p_ref <- pool[pool$zone == "portal", ]
  c_ref <- pool[pool$zone == "central", ]
  if (nrow(p_ref) < 1 || nrow(c_ref) < 1) return(NULL)
  obs <- median(d$d_portal) / median(d$d_central)

  nb <- nrow(p_ref); nc <- nrow(c_ref)
  nul <- replicate(500, {
    idx <- sample(nrow(pool), nb + nc)
    a <- pool[idx[seq_len(nb)], ]
    b <- pool[idx[(nb + 1):(nb + nc)], ]
    median(compute_nn(d, a), na.rm = TRUE) /
      median(compute_nn(d, b), na.rm = TRUE)
  })
  data.frame(fov = f, obs_ratio = obs,
             null_median = median(nul, na.rm = TRUE),
             # ratio = d_portal/d_central; larger = closer to central, so the
             # upper tail is the hypothesis of interest
             p_perm = (sum(nul >= obs) + 1) / (length(nul) + 1))
}))
perm <- perm[complete.cases(perm), ]
cat(sprintf("FOVs tested: %d\n", nrow(perm)))
cat(sprintf("Observed median d_portal/d_central = %.3f (1.0 = no preference)\n",
            median(perm$obs_ratio)))
cat(sprintf("Null median                        = %.3f\n", median(perm$null_median)))
cat(sprintf("FOVs favouring central with P < 0.05: %d\n",
            sum(perm$p_perm < 0.05 & perm$obs_ratio > 1)))
write.csv(perm, res_path("03_spatial", "zonation_permutation.csv"), row.names = FALSE)

# ---- 4. HGF+ vs HGF- HSC zonation preference ----------------------------
hsc <- df[df$celltype == "Stellate.cells", ]
hsc$d_portal  <- compute_nn(hsc, tgt_portal)
hsc$d_central <- compute_nn(hsc, tgt_central)
hsc <- hsc[!is.na(hsc$d_portal) & !is.na(hsc$d_central), ]
hsc$HGF_pos <- as.numeric(hsc$HGF > 0)
hsc$log_ratio <- log2((hsc$d_central + 1) / (hsc$d_portal + 1))

m_hsc <- lmer(log_ratio ~ HGF_pos + (1 | patient_uid), data = hsc, REML = FALSE)
cth <- summary(m_hsc)$coefficients
cat(sprintf("\nHSC zonation: HGF+ vs HGF- beta = %.3f (t = %.1f, P = %.2e)\n",
            cth["HGF_pos", 1], cth["HGF_pos", 4], cth["HGF_pos", 5]))
cat(sprintf("  HGF+ HSC median d_portal = %.1f, d_central = %.1f\n",
            median(hsc$d_portal[hsc$HGF_pos == 1]), median(hsc$d_central[hsc$HGF_pos == 1])))
write.csv(hsc[, c("fov", "patient_uid", "HGF_pos", "d_portal", "d_central", "log_ratio")],
          res_path("03_spatial", "zonation_hsc.csv"), row.names = FALSE)

saveRDS(list(hepg4 = hep4, permutation = perm, density = dens,
             model_uncorrected = ct, model_normalised = ctn,
             hsc_model = cth, portal = PORTAL, central = CENTRAL),
        file.path(DATA_DIR, "derived", "zonation.rds"))

cat("\nNext: 03_spatial/19_Hep4_IL6R_and_pathways.R\n")
