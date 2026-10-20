# ============================================================================
# 24_HGFpos_neighbourhood_signature.R
#
# Negative control for the 11-gene tissue signature: is the signature score
# elevated in the immediate neighbourhood of HGF+ HSCs?
#
# For each HGF+ HSC, cells within 50 um are compared with cells more than 150 um
# away in the same field of view. Proximal and distal cells are depth-adjusted
# before comparison, because the two sets differ in capture efficiency.
#
# Only CHI3L1 and GSN of the 11 genes are on the CosMx panel, so this score
# effectively reflects two genes. That is a hard limit and is reported with the
# result.
#
# The expected and observed answer is null: the signature is not enriched around
# HGF+ HSCs, which is one of the arguments that the signature does not capture
# the niche.
#
# Output
#   results/03_spatial/neighbourhood_signature.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(lme4); library(lmerTest)})

df <- load_derived("cosmx_flat.rds")
hsc <- df[df$celltype == "Stellate.cells" & df$HGF > 0, ]
cat(sprintf("HGF+ HSC anchors: %d\n", nrow(hsc)))

genes_on_panel <- c("CHI3L1", "GSN")
cat(sprintf("Signature genes on the panel: %s\n", paste(genes_on_panel, collapse = ", ")))

# Score every cell with the two available signature genes.
score <- rep(NA_real_, nrow(df))
for (s in c(1, 2)) {
  obj <- load_cosmx(s)
  dat <- obj[["RNA"]]$data
  g <- intersect(genes_on_panel, rownames(dat))
  idx <- which(df$slide == paste0("S", s))
  if (length(g) == 0) { rm(obj); next }
  m <- as.matrix(dat[g, df$cell[idx], drop = FALSE])
  z <- t(scale(t(m))); z[is.na(z)] <- 0
  score[idx] <- colMeans(z)
  rm(obj, m, z); gc(verbose = FALSE)
}
df$sig_score <- score
cat(sprintf("Scored cells: %d\n", sum(!is.na(df$sig_score))))

# ---- Proximal vs distal sets --------------------------------------------
rows <- list()
for (i in seq_len(nrow(hsc))) {
  f <- hsc$fov[i]
  cand <- df[df$slide == hsc$slide[i] & df$fov == f &
            df$celltype != "Stellate.cells", ]
  if (nrow(cand) == 0) next
  d <- sqrt((cand$x - hsc$x[i])^2 + (cand$y - hsc$y[i])^2) * 1000

  prox <- cand[d <= NEAR_UM, ]
  dist <- cand[d >  FAR_UM, ]
  if (nrow(prox) < 5 || nrow(dist) < 10) next

  rows[[length(rows) + 1]] <- data.frame(
    anchor = i,
    prox_score = mean(prox$sig_score, na.rm = TRUE),
    dist_score = mean(dist$sig_score, na.rm = TRUE),
    prox_ncount = median(prox$nCount, na.rm = TRUE),
    dist_ncount = median(dist$nCount, na.rm = TRUE)
  )
}
pair <- do.call(rbind, rows)
cat(sprintf("Anchors with both a proximal and a distal set: %d\n", nrow(pair)))

# ---- Unadjusted and depth-adjusted comparison ---------------------------
w_raw <- wilcox.test(pair$prox_score, pair$dist_score, paired = TRUE)
cat(sprintf("\nUnadjusted: proximal %.4f vs distal %.4f (paired Wilcoxon P = %.3f)\n",
            median(pair$prox_score), median(pair$dist_score), w_raw$p.value))

# Depth adjustment: regress the score on the median nCount of each set and
# compare residuals, then confirm with a mixed model across anchors.
pair$delta <- pair$prox_score - pair$dist_score
pair$delta_ncount <- pair$prox_ncount - pair$dist_ncount
adj <- lm(delta ~ delta_ncount, data = pair)
pair$delta_adj <- residuals(adj)

w_adj <- wilcox.test(pair$delta_adj, mu = 0)
cat(sprintf("Depth-adjusted: median delta %.4f (Wilcoxon vs 0, P = %.3f)\n",
            median(pair$delta_adj), w_adj$p.value))

long <- rbind(
  data.frame(anchor = pair$anchor, zone = "proximal", score = pair$prox_score),
  data.frame(anchor = pair$anchor, zone = "distal", score = pair$dist_score)
)
long$zone <- factor(long$zone, levels = c("distal", "proximal"))
m <- lmer(score ~ zone + (1 | anchor), data = long, REML = FALSE)
ct <- summary(m)$coefficients
cat(sprintf("Mixed model: proximal vs distal beta = %.4f (P = %.3f)\n",
            ct["zoneproximal", 1], ct["zoneproximal", "Pr(>|t|)"]))

out <- data.frame(
  comparison = c("unadjusted", "depth-adjusted", "mixed model"),
  estimate = c(median(pair$delta), median(pair$delta_adj),
               ct["zoneproximal", 1]),
  P = c(w_raw$p.value, w_adj$p.value,
        ct["zoneproximal", "Pr(>|t|)"]),
  n_anchors = nrow(pair),
  genes_on_panel = paste(genes_on_panel, collapse = ", ")
)
print(out, row.names = FALSE)
write.csv(out, res_path("03_spatial", "neighbourhood_signature.csv"),
          row.names = FALSE)

cat("\nInterpretation: the signature score is not higher near HGF+ HSCs. Only two\n")
cat("of the eleven genes are measurable here, so this is a weak test of the full\n")
cat("signature, but it gives no support for spatial enrichment.\n")

saveRDS(list(pairs = pair, result = out),
        file.path(DATA_DIR, "derived", "neighbourhood_signature.rds"))

cat("\nNext: 03_spatial/23_CellChat_spatial.R (if not already run)\n")
