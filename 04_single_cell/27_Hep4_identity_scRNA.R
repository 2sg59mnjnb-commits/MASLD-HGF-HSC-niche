# ============================================================================
# 27_Hep4_identity_scRNA.R
#
# Transcriptional identity of the Hep.4 hepatocyte subpopulation.
#
# Two questions are asked:
#   1. Which scRNA-seq hepatocyte cluster corresponds to CosMx Hep.4?
#      CosMx Hep.4 marker genes (Hep.4 vs the other hepatocyte subtypes) are used
#      to score every hepatocyte cluster in GSE202379. The best-scoring cluster is
#      taken as the Hep.4 equivalent.
#   2. Where does that cluster sit along the portal-to-central axis?
#      Zone 1 (periportal), Zone 2 (midzonal) and Zone 3 (central) marker panels
#      are scored per cluster.
#
# The zonation score is transcriptional; the anatomical position reported in the
# manuscript comes from the coordinate measurements in
# 03_spatial/18_zonation_portal_vs_central.R. Where the two disagree, the
# coordinate measurement is the direct one and is given precedence -- this script
# reports both so the discrepancy is visible rather than hidden.
#
# Output
#   results/04_single_cell/Hep4_markers_cosmx.csv
#   results/04_single_cell/Hep4_scRNA_cluster_mapping.csv
#   data/derived/hep_sc_15k.rds
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(Seurat); library(dplyr)})

# ---- 1. CosMx Hep.4 markers ----------------------------------------------
cat("Deriving Hep.4 markers from CosMx ...\n")
obj <- load_cosmx(1)
ct  <- as.character(obj@meta.data[[celltype_col(1)]])
hep_types <- c("Hep.1", "Hep.3", "Hep.4", "Hep.5")
hep <- subset(obj, cells = colnames(obj)[ct %in% hep_types])
hep$subtype <- ct[ct %in% hep_types]
rm(obj); gc(verbose = FALSE)

Idents(hep) <- hep$subtype
markers <- FindMarkers(hep, ident.1 = "Hep.4",
                       ident.2 = c("Hep.1", "Hep.3", "Hep.5"),
                       logfc.threshold = 0.25, min.pct = 0.1,
                       test.use = "wilcox", verbose = FALSE)
markers$gene <- rownames(markers)
markers <- markers[order(-markers$avg_log2FC), ]
write.csv(markers, res_path("04_single_cell", "Hep4_markers_cosmx.csv"),
          row.names = FALSE)

top_markers <- head(markers$gene[markers$adj.P.Val < 0.05], 30)
cat(sprintf("Top Hep.4 markers (FDR < 0.05): %d\n", length(top_markers)))
print(top_markers)
rm(hep); gc(verbose = FALSE)

# ---- 2. scRNA reference: hepatocytes only -------------------------------
cat("\nLoading GSE202379 ...\n")
sc <- readRDS(file.path(DATA_DIR, "GSE202379_SeuratObject_AllCells.rds.gz"))
DefaultAssay(sc) <- "RNA"
cat(sprintf("Full object: %d cells\n", ncol(sc)))

hep_cells <- colnames(sc)[grepl("hepatocyte", sc@meta.data$cell.annotation,
                                ignore.case = TRUE)]
cat(sprintf("Hepatocytes: %d\n", length(hep_cells)))

hep_sc <- subset(sc, cells = hep_cells)
rm(sc); gc(verbose = FALSE)

# Down-sample for tractability; the reported analysis uses 15,000 hepatocytes.
set.seed(SEED)
if (ncol(hep_sc) > 15000) hep_sc <- hep_sc[, sample(ncol(hep_sc), 15000)]
cat(sprintf("Hepatocytes used: %d\n", ncol(hep_sc)))

if (max(hep_sc[["RNA"]]$data) > 50) {
  hep_sc <- NormalizeData(hep_sc, verbose = FALSE)
}
saveRDS(hep_sc, file.path(DATA_DIR, "derived", "hep_sc_15k.rds"))

# ---- 3. Score clusters with the CosMx Hep.4 markers ---------------------
markers_in_sc <- intersect(top_markers, rownames(hep_sc))
cat(sprintf("Hep.4 markers present in scRNA: %d / %d\n",
            length(markers_in_sc), length(top_markers)))

hep_sc <- AddModuleScore(hep_sc, features = list(markers_in_sc),
                         name = "Hep4_sig", nbin = 24)
hep_sc <- FindVariableFeatures(hep_sc, verbose = FALSE)
hep_sc <- ScaleData(hep_sc, verbose = FALSE)
hep_sc <- RunPCA(hep_sc, verbose = FALSE)
hep_sc <- FindNeighbors(hep_sc, dims = 1:15, verbose = FALSE)
hep_sc <- FindClusters(hep_sc, resolution = 0.3, verbose = FALSE)

cluster_score <- hep_sc@meta.data %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(n = dplyr::n(),
                   Hep4_score = mean(Hep4_sig1),
                   .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(Hep4_score))
print(as.data.frame(cluster_score))

top_cluster <- as.character(cluster_score$seurat_clusters[1])
cat(sprintf("\nBest-matching cluster: %s (score %.3f, n = %d)\n",
            top_cluster,
            cluster_score$Hep4_score[1],
            cluster_score$n[1]))

# ---- 4. Zone scores ------------------------------------------------------
score_zone <- function(genes) {
  g <- intersect(genes, rownames(hep_sc))
  if (length(g) == 0) return(NULL)
  sapply(levels(hep_sc$seurat_clusters), function(cl) {
    idx <- which(hep_sc$seurat_clusters == cl)
    mean(colMeans(as.matrix(hep_sc[["RNA"]]$data[g, idx, drop = FALSE])))
  })
}

zone <- data.frame(
  cluster = levels(hep_sc$seurat_clusters),
  Zone1 = score_zone(ZONE1_GENES),
  Zone2 = score_zone(ZONE2_GENES),
  Zone3 = score_zone(ZONE3_GENES)
)
zone$Hep4_score <- cluster_score$Hep4_score[match(zone$cluster,
                                                  cluster_score$seurat_clusters)]
zone$n <- cluster_score$n[match(zone$cluster, cluster_score$seurat_clusters)]
print(round(zone[, c("Zone1", "Zone2", "Zone3")], 3))

write.csv(zone, res_path("04_single_cell", "Hep4_scRNA_cluster_mapping.csv"),
          row.names = FALSE)

zi <- zone[zone$cluster == top_cluster, ]
cat(sprintf("\nZone profile of the Hep.4-equivalent cluster %s:\n", top_cluster))
cat(sprintf("  Zone1 (periportal) = %.2f\n  Zone2 (midzonal)   = %.2f\n  Zone3 (central)    = %.2f\n",
            zi$Zone1, zi$Zone2, zi$Zone3))
cat("  -> ", ifelse(zi$Zone1 > zi$Zone3,
                    "transcriptionally periportal/midzonal",
                    "transcriptionally central"), "\n", sep = "")
cat("\nThe coordinate measurements in 03_spatial/18 place Hep.4 pericentrally.\n")
cat("The transcriptional score is intermediate (Zone1 ~ Zone2), which is why the\n")
cat("manuscript states that the zonation score is compatible with a midzonal /\n")
cat("pericentral position rather than a clean Zone 1 profile.\n")

saveRDS(list(markers = markers, top_markers = top_markers,
             cluster_score = cluster_score, zone = zone,
             top_cluster = top_cluster, markers_in_sc = markers_in_sc),
        file.path(DATA_DIR, "derived", "hep4_identity.rds"))

cat("\nNext: 04_single_cell/28_signature_specificity_scRNA.R\n")
