# ============================================================================
# 13_load_CosMx_and_panel.R
#
# Loads both CosMx slides (GSE292268), reconstructs the patient label from the
# FOV x-coordinates, summarises cell-type composition and checks which genes of
# interest are present on the 1,000-plex panel.
#
# The panel check matters because absence of a ligand from the panel is a
# detection limit, not evidence of absent signalling.
#
# Output
#   results/03_spatial/celltype_composition.csv
#   results/03_spatial/panel_coverage.csv
#   data/derived/cosmx_flat.rds
# ============================================================================

source("R/config.R")
source("R/helpers.R")

dir.create(file.path(DATA_DIR, "derived"), showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load both slides --------------------------------------------------
genes <- c("HGF", "IL6R", "IL11RA", "MET", "STAT3", "ALB", "KRT5", "KRT14")
parts <- list()
for (s in c(1, 2)) {
  obj <- load_cosmx(s)
  cat(sprintf("Slide%d: %d cells, %d genes on panel\n", s, ncol(obj), nrow(obj)))
  parts[[s]] <- flatten_cosmx(obj, s, genes)
  rm(obj); gc(verbose = FALSE)
}
df <- dplyr::bind_rows(parts)
rm(parts); gc(verbose = FALSE)
df$patient_uid <- paste(df$slide, df$patient, sep = "_")

cat(sprintf("\nTotal cells: %d\n", nrow(df)))
cat("Cells per slide and patient:\n")
print(table(df$slide, df$patient))
cat("Cells per FOV:\n")
print(tapply(df$fov, df$slide, function(x) length(unique(x))))

# ---- 2. Cell-type composition --------------------------------------------
comp <- df %>%
  dplyr::group_by(slide, celltype) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(slide) %>%
  dplyr::mutate(pct = round(100 * n / sum(n), 2)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(slide, dplyr::desc(n))
write.csv(comp, res_path("03_spatial", "celltype_composition.csv"), row.names = FALSE)

for (s in c("S1", "S2")) {
  sub <- comp[comp$slide == s, ]
  st  <- sub[sub$celltype == "Stellate.cells", ]
  cat(sprintf("\n%s: %d cell types; Stellate = %d (%.1f%%)\n",
              s, nrow(sub), st$n, st$pct))
}

# HGF positivity within HSCs
hsc <- df[df$celltype == "Stellate.cells", ]
cat("\nHGF+ HSC fraction:\n")
for (s in c("S1", "S2")) {
  sub <- hsc[hsc$slide == s, ]
  cat(sprintf("  %s: %d / %d = %.1f%%\n", s,
              sum(sub$HGF > 0), nrow(sub), 100 * mean(sub$HGF > 0)))
}

# ---- 3. Panel coverage ----------------------------------------------------
# Use the panel of one slide; both slides carry the same 1,000-plex panel.
panel_obj <- load_cosmx(1, add_patient_label = FALSE)
panel <- rownames(panel_obj[["RNA"]])
rm(panel_obj); gc(verbose = FALSE)

gene_sets <- list(
  "IL6-family ligands"   = IL6_LIGANDS,
  "IL6-family receptors" = IL6_RECEPTORS,
  "HGF-MET axis"         = c("HGF", "MET"),
  "JAK/STAT core"        = JAK_STAT_GENES,
  "NF-kB core"           = NFKB_GENES,
  "11-gene signature"    = SIG11,
  "16-gene signature"    = SIG16,
  "Zonation markers"     = unique(c(ZONE1_GENES, ZONE2_GENES, ZONE3_GENES)),
  "Negative controls"    = NEG_CONTROL_GENES
)

cov <- do.call(rbind, lapply(names(gene_sets), function(nm) {
  g <- gene_sets[[nm]]
  data.frame(
    gene_set = nm,
    n_total = length(g),
    n_on_panel = sum(g %in% panel),
    coverage = sprintf("%d/%d", sum(g %in% panel), length(g)),
    missing = paste(setdiff(g, panel), collapse = ", "),
    stringsAsFactors = FALSE
  )
}))
print(cov)
write.csv(cov, res_path("03_spatial", "panel_coverage.csv"), row.names = FALSE)

cat(sprintf("\nIL6-family ligands on panel: %d / %d\n",
            sum(IL6_LIGANDS %in% panel), length(IL6_LIGANDS)))
cat(sprintf("IL6-family receptors on panel: %d / %d\n",
            sum(IL6_RECEPTORS %in% panel), length(IL6_RECEPTORS)))
cat("Note: missing ligands constrain CellChat, so a negative CellChat result for\n")
cat("      those ligands is a detection limit rather than biological absence.\n")

# ---- 4. Negative-control detection rate ----------------------------------
cat("\nNegative-control positive rates (ambient RNA check):\n")
for (g in intersect(NEG_CONTROL_GENES, colnames(df))) {
  cat(sprintf("  %-7s positive in %.1f%% of cells\n", g, 100 * mean(df[[g]] > 0, na.rm = TRUE)))
}

saveRDS(df, file.path(DATA_DIR, "derived", "cosmx_flat.rds"))

cat("\nNext: 03_spatial/14_stellate_fibrotic_enrichment.R\n")
