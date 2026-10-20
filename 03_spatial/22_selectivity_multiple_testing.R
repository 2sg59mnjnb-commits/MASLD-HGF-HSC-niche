# ============================================================================
# 22_selectivity_multiple_testing.R
#
# Is the HGF+ HSC niche specific to Hep.4 hepatocytes, or are HGF+ HSCs simply
# close to everything? For each HGF+ HSC the distance to the nearest Hep.4 is
# compared with its distance to the nearest cell of every other annotated type,
# using paired Wilcoxon tests with Benjamini-Hochberg correction.
#
# Reporting note: the manuscript states that HGF+ HSCs lay closer to Hep.4 than
# to 17 of 24 other annotated cell types. The number of testable types depends on
# the minimum-abundance filter below, so the script prints both the raw and the
# filtered counts.
#
# Output
#   results/03_spatial/selectivity_paired.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

df <- load_derived("cosmx_flat.rds")
hep4 <- df[df$celltype == "Hep.4", ]
hsc  <- df[df$celltype == "Stellate.cells", ]
pos  <- hsc[hsc$HGF > 0, ]
cat(sprintf("HGF+ HSC: %d\n", nrow(pos)))

MIN_CELLS <- 30   # a cell type needs at least this many cells to give a distance

d_hep4 <- compute_nn(pos, hep4)
types <- setdiff(sort(unique(df$celltype)), "Stellate.cells")

rows <- lapply(types, function(ct) {
  tgt <- df[df$celltype == ct, ]
  if (nrow(tgt) < MIN_CELLS) return(NULL)
  d_ct <- compute_nn(pos, tgt)
  ok <- !is.na(d_hep4) & !is.na(d_ct)
  if (sum(ok) < 30) return(NULL)

  w <- wilcox.test(d_hep4[ok], d_ct[ok], paired = TRUE, alternative = "less")
  data.frame(
    celltype = ct,
    n_target = nrow(tgt),
    median_d_hep4 = round(median(d_hep4[ok]), 1),
    median_d_other = round(median(d_ct[ok]), 1),
    # positive delta = Hep.4 is closer than this cell type
    delta_um = round(median(d_ct[ok]) - median(d_hep4[ok]), 1),
    P = w$p.value,
    stringsAsFactors = FALSE
  )
})
res <- do.call(rbind, rows)
res$FDR <- p.adjust(res$P, method = "BH")
res <- res[order(-res$delta_um), ]
print(res, row.names = FALSE)
write.csv(res, res_path("03_spatial", "selectivity_paired.csv"), row.names = FALSE)

n_tested <- nrow(res)
n_closer <- sum(res$delta_um > 0)
n_sig    <- sum(res$delta_um > 0 & res$FDR < 0.05)
cat(sprintf("\nCell types tested: %d\n", n_tested))
cat(sprintf("Hep.4 closer than the other type (delta > 0): %d / %d\n", n_closer, n_tested))
cat(sprintf("Significant after BH correction:              %d / %d\n", n_sig, n_tested))
cat(sprintf("Hep.4 farther than the other type:            %d / %d\n",
            sum(res$delta_um < 0), n_tested))

compare("types where Hep.4 is closer", n_closer, 17, digits = 0)

cat("\nTypes that are closer to HGF+ HSCs than Hep.4 is (delta < 0):\n")
print(head(res[res$delta_um < 0, c("celltype", "delta_um", "FDR")], 10))

saveRDS(list(selectivity = res, n_tested = n_tested,
             n_closer = n_closer, n_sig = n_sig),
        file.path(DATA_DIR, "derived", "selectivity.rds"))

cat("\nNext: 03_spatial/23_CellChat_spatial.R\n")
