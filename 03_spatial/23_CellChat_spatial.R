# ============================================================================
# 23_CellChat_spatial.R
#
# Spatially constrained ligand-receptor inference with CellChat v2.
#
# Cell types are refined so that HSC and Hep.4 identity is resolved by the genes
# of interest:
#   Stellate_HGFpos / Stellate_HGFneg
#   Hep.4_IL6Rpos  / Hep.4_IL6Rneg
#
# Spatial parameters (as reported in the Methods)
#   conversion factor 0.18 um/pixel
#   maximum interaction range 250 um
#   contact range 10 um
#   distance-decay scaling factor 3.0
#   trimean estimator, 10% trimming, 100 bootstrap permutations
#
# Important: the 1,000-plex panel does not contain all IL6-family ligands, so
# the absence of an IL6-family interaction is a detection limit and must not be
# reported as biological absence. See 03_spatial/13 for the coverage table.
#
# Output
#   results/03_spatial/CellChat_all_LR.csv
#   results/03_spatial/CellChat_HGFpos_to_IL6Rpos.csv
#   results/03_spatial/CellChat_HGFpos_pathways.csv
#   data/derived/cellchat.rds
# ============================================================================

source("R/config.R")
source("R/helpers.R")

if (!requireNamespace("CellChat", quietly = TRUE)) {
  cat("CellChat is not installed; skipping.\n")
  cat("Install with: remotes::install_github('jinworks/CellChat')\n")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(CellChat); library(Seurat); library(dplyr)
})

# ---- 1. Build the refined object -----------------------------------------
prepare <- function(slide) {
  obj <- load_cosmx(slide)
  ct  <- as.character(obj@meta.data[[celltype_col(slide)]])
  dat <- obj[["RNA"]]$data
  hgf  <- as.numeric(dat["HGF", ])
  il6r <- as.numeric(dat["IL6R", ])

  refined <- ct
  refined[ct == "Stellate.cells" & hgf  > 0] <- "Stellate_HGFpos"
  refined[ct == "Stellate.cells" & hgf == 0] <- "Stellate_HGFneg"
  refined[ct == "Hep.4" & il6r > 0] <- "Hep.4_IL6Rpos"
  refined[ct == "Hep.4" & il6r == 0] <- "Hep.4_IL6Rneg"

  keep <- c("Stellate_HGFpos", "Stellate_HGFneg",
            "Hep.4_IL6Rpos", "Hep.4_IL6Rneg",
            "Hep.1", "Hep.2", "Hep.3", "Hep.5", "Hep.6",
            "Cholangiocytes", "Inflammatory.macrophages",
            "Non.inflammatory.macrophages", "Endothelial",
            "Periportal.LSECs", "Central.venous.LSECs")
  keep <- keep[keep %in% unique(refined)]

  sub <- subset(obj, cells = colnames(obj)[refined %in% keep])
  sub$celltype_refined <- refined[refined %in% keep]
  sub$slide <- paste0("S", slide)
  rm(obj); gc(verbose = FALSE)
  sub
}

cat("Preparing slides ...\n")
o1 <- prepare(1); o2 <- prepare(2)
obj <- merge(o1, o2); rm(o1, o2); gc(verbose = FALSE)
DefaultAssay(obj) <- "RNA"
if (max(obj[["RNA"]]$data) > 50) obj <- NormalizeData(obj, verbose = FALSE)
cat(sprintf("Merged object: %d cells, %d cell types\n",
            ncol(obj), length(unique(obj$celltype_refined))))
print(table(obj$celltype_refined))

# ---- 2. Spatial coordinates in pixels ------------------------------------
locs <- data.frame(
  x = obj@meta.data$x_slide_mm * 1000 / COSMX_UM_PER_PX,
  y = obj@meta.data$y_slide_mm * 1000 / COSMX_UM_PER_PX
)
rownames(locs) <- colnames(obj)

# Spot tolerance = half the minimum nearest-neighbour distance (CellChat guide)
nn <- RANN::nn2(locs, locs, k = 2)
min_nn <- min(nn$nn.dists[, 2])
spatial_factors <- data.frame(ratio = COSMX_UM_PER_PX, tol = min_nn / 2)
cat(sprintf("Spatial factors: ratio = %.2f, tol = %.3f px (min NN %.2f px = %.2f um)\n",
            spatial_factors$ratio, spatial_factors$tol,
            min_nn, min_nn * COSMX_UM_PER_PX))

# ---- 3. Run CellChat -----------------------------------------------------
meta <- data.frame(labels = obj$celltype_refined,
                   samples = obj$slide,
                   row.names = colnames(obj))

cc <- createCellChat(object = obj[["RNA"]]$data, meta = meta, group.by = "labels",
                     datatype = "spatial", coordinates = locs,
                     spatial.factors = spatial_factors)
cc@DB <- CellChatDB.human
cat(sprintf("Ligand-receptor pairs in database: %d\n", nrow(CellChatDB.human$interaction)))

cc <- subsetData(cc)
cc <- identifyOverExpressedGenes(cc)
cc <- identifyOverExpressedInteractions(cc)

future::plan("multisession", workers = 2)
cc <- computeCommunProb(
  cc,
  type = "triMean", trim = 0.1,
  distance.use = TRUE,
  interaction.range = 250,
  scale.distance = 3.0,
  contact.dependent = TRUE,
  contact.range = 10,
  nboot = 100
)
cc <- filterCommunication(cc, min.cells = 10)

# ---- 4. Extract interactions --------------------------------------------
lr <- subsetCommunication(cc)
cat(sprintf("\nSignificant ligand-receptor pairs: %d\n", nrow(lr)))
write.csv(lr, res_path("03_spatial", "CellChat_all_LR.csv"), row.names = FALSE)

to_target <- lr[lr$source == "Stellate_HGFpos" &
                  lr$target == "Hep.4_IL6Rpos", ]
cat(sprintf("HGF+ HSC -> IL6R+ Hep.4: %d pairs\n", nrow(to_target)))
if (nrow(to_target)) {
  print(to_target[, c("ligand", "receptor", "pathway_name", "prob", "pval")])
}
write.csv(to_target, res_path("03_spatial", "CellChat_HGFpos_to_IL6Rpos.csv"),
          row.names = FALSE)

hgf_met <- lr[lr$ligand == "HGF" & lr$receptor == "MET", ]
cat(sprintf("\nHGF-MET interactions detected: %d\n", nrow(hgf_met)))
if (nrow(hgf_met)) print(hgf_met[, c("source", "target", "prob", "pval")])

il6_fam <- intersect(c(IL6_LIGANDS, IL6_RECEPTORS), colnames(lr))
il6_hits <- lr[lr$ligand %in% IL6_LIGANDS, ]
cat(sprintf("IL6-family ligand interactions detected: %d\n", nrow(il6_hits)))

# ---- 5. Pathway summary --------------------------------------------------
from_pos <- lr[lr$source == "Stellate_HGFpos", ]
path_tab <- from_pos %>%
  dplyr::group_by(pathway_name) %>%
  dplyr::summarise(n_pairs = dplyr::n(),
                   total_prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(total_prob))
cat("\nPathways sent by HGF+ HSCs:\n")
print(head(as.data.frame(path_tab), 15))
write.csv(path_tab, res_path("03_spatial", "CellChat_HGFpos_pathways.csv"),
          row.names = FALSE)

saveRDS(list(cellchat = cc, lr = lr, to_il6rpos = to_target,
             hgf_met = hgf_met, pathways = path_tab),
        file.path(DATA_DIR, "derived", "cellchat.rds"))

cat("\nNote: report the panel-coverage table (03_spatial/13) alongside any negative\n")
cat("CellChat result, because missing ligands cannot be tested.\n")
cat("\nNext: 04_single_cell/24_cosmx_HGFpos_vs_HGFneg_DEG.R\n")
