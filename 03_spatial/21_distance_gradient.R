# ============================================================================
# 21_distance_gradient.R
#
# Distance-gradient analysis: how does Hep.4 enrichment around HGF+ HSCs change
# with radius? A genuine niche should show enrichment at short range that decays
# towards 1.0 (no enrichment) at larger radii.
#
# Enrichment ratio at radius r
#   observed = mean number of Hep.4 cells within r of an HGF+ HSC
#   expected = the same mean for an equal number of non-HSC cells drawn at
#              random within the same FOV (1,000 permutations)
#   ratio    = observed / expected
#
# HGF- HSCs are processed identically as the comparator.
#
# Output
#   results/03_spatial/distance_gradient.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

df <- load_derived("cosmx_flat.rds")
hep4 <- df[df$celltype == "Hep.4", ]
hsc  <- df[df$celltype == "Stellate.cells", ]
anchors <- list(`HGF+ HSC` = hsc[hsc$HGF > 0, ],
                `HGF- HSC` = hsc[hsc$HGF == 0, ])

radii <- c(25, 50, 75, 100, 150)
cat(sprintf("Cells: %d | Hep.4: %d | HGF+ HSC: %d | HGF- HSC: %d\n",
            nrow(df), nrow(hep4), nrow(anchors[[1]]), nrow(anchors[[2]])))

# Count Hep.4 neighbours once per radius for every cell, then reuse the vector
# for the observed statistic and all permutations. This keeps the 1,000-iteration
# permutation test cheap.
cat("Counting Hep.4 neighbours at each radius ...\n")
for (r in radii) df[[paste0("n_", r)]] <- count_within(df, hep4, r)

# The anchor sets are rebuilt here because they were taken from df before the
# per-radius count columns existed; without this they would not carry n_<r>.
hsc <- df[df$celltype == "Stellate.cells", ]
anchors <- list(`HGF+ HSC` = hsc[hsc$HGF > 0, ],
                `HGF- HSC` = hsc[hsc$HGF == 0, ])

# Non-HSC pool used to build the null distribution
pool <- df[df$celltype != "Stellate.cells", ]
pool_fov <- split(seq_len(nrow(pool)), fov_key(pool))

set.seed(SEED)
n_perm <- 1000
results <- list()

for (nm in names(anchors)) {
  a <- anchors[[nm]]
  a_fov <- split(seq_len(nrow(a)), fov_key(a))
  cat(sprintf("\n=== %s (n = %d) ===\n", nm, nrow(a)))

  for (r in radii) {
    col <- paste0("n_", r)
    obs_vals <- a[[col]]

    # Pooled observed mean
    obs <- mean(obs_vals, na.rm = TRUE)

    # Null: for each FOV, draw as many non-HSC cells as there are anchors
    null <- numeric(n_perm)
    for (p in seq_len(n_perm)) {
      vals <- numeric(0)
      for (f in names(a_fov)) {
        n_f <- length(a_fov[[f]])
        cand <- pool_fov[[f]]
        if (is.null(cand) || length(cand) < n_f) next
        vals <- c(vals, pool[[col]][sample(cand, n_f)])
      }
      null[p] <- mean(vals, na.rm = TRUE)
    }

    exp_mean <- mean(null)
    results[[length(results) + 1]] <- data.frame(
      anchor = nm, radius_um = r,
      n_anchor = nrow(a),
      observed = round(obs, 4),
      expected = round(exp_mean, 4),
      enrichment = round(obs / exp_mean, 3),
      P_perm = (sum(null >= obs) + 1) / (n_perm + 1)
    )
  }
}

grad <- do.call(rbind, results)
print(grad, row.names = FALSE)
write.csv(grad, res_path("03_spatial", "distance_gradient.csv"), row.names = FALSE)

cat("\nInterpretation: enrichment > 1 at short range decaying towards 1.0 at\n")
cat("larger radii indicates a local niche rather than a field-wide gradient.\n")

saveRDS(list(gradient = grad, radii = radii, n_perm = n_perm),
        file.path(DATA_DIR, "derived", "distance_gradient.rds"))

cat("\nNext: 03_spatial/22_selectivity_multiple_testing.R\n")
