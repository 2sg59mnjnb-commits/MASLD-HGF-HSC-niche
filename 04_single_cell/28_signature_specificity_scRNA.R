# ============================================================================
# 28_signature_specificity_scRNA.R
#
# Cell-type specificity of the two signatures in an independent scRNA-seq
# reference (GSE202379).
#
# The comparison is the point: the 16-gene HGF+ HSC signature is expected to peak
# in stellate cells, whereas the 11-gene tissue signature is expected to be
# distributed and non-specific. The specificity index below quantifies the gap
# between the top-scoring and second-scoring cell type for each signature.
#
# Output
#   results/04_single_cell/signature_score_by_celltype.csv
#   results/04_single_cell/signature_specificity_index.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({library(Seurat); library(dplyr)})

cat("Loading GSE202379 ...\n")
sc <- readRDS(file.path(DATA_DIR, "GSE202379_SeuratObject_AllCells.rds.gz"))
DefaultAssay(sc) <- "RNA"
sc <- subset(sc, subset = cell.annotation != "unknown")
if (max(sc[["RNA"]]$data) > 50) sc <- NormalizeData(sc, verbose = FALSE)
cat(sprintf("Cells: %d, cell types: %d\n",
            ncol(sc), length(unique(sc$cell.annotation))))

sig11_in <- intersect(SIG11, rownames(sc))
sig16_in <- intersect(SIG16, rownames(sc))
cat(sprintf("11-gene available: %d / 11; 16-gene available: %d / 16\n",
            length(sig11_in), length(sig16_in)))

sc <- AddModuleScore(sc, features = list(sig11_in), name = "Sig11", nbin = 24)
sc <- AddModuleScore(sc, features = list(sig16_in), name = "Sig16", nbin = 24)

by_type <- sc@meta.data %>%
  dplyr::group_by(cell.annotation) %>%
  dplyr::summarise(n = dplyr::n(),
                   sig11 = mean(Sig111),
                   sig16 = mean(Sig161),
                   .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(sig16))
print(as.data.frame(by_type))
write.csv(by_type, res_path("04_single_cell", "signature_score_by_celltype.csv"),
          row.names = FALSE)

# ---- Specificity index ---------------------------------------------------
top2 <- function(x) sort(x, decreasing = TRUE)[1:2]
s11 <- top2(by_type$sig11); s16 <- top2(by_type$sig16)
idx11 <- s11[1] - s11[2]
idx16 <- s16[1] - s16[2]

top11 <- by_type$cell.annotation[order(-by_type$sig11)][1:2]
top16 <- by_type$cell.annotation[order(-by_type$sig16)][1:2]

cat(sprintf("\n11-gene: highest in %s (%.3f), second %s (%.3f); index = %.3f\n",
            top11[1], s11[1], top11[2], s11[2], idx11))
cat(sprintf("16-gene: highest in %s (%.3f), second %s (%.3f); index = %.3f\n",
            top16[1], s16[1], top16[2], s16[2], idx16))

spec <- data.frame(
  signature = c("11-gene tissue signature", "16-gene HGF+ HSC signature"),
  top_celltype = c(top11[1], top16[1]),
  top_score = round(c(s11[1], s16[1]), 4),
  second_celltype = c(top11[2], top16[2]),
  second_score = round(c(s11[2], s16[2]), 4),
  specificity_index = round(c(idx11, idx16), 4)
)
print(spec)
write.csv(spec, res_path("04_single_cell", "signature_specificity_index.csv"),
          row.names = FALSE)

cat("\nInterpretation: a large index means the signature is concentrated in one\n")
cat("cell type; an index close to zero means the signal is spread across types and\n")
cat("the signature should not be described as cell-type specific.\n")

# Stellate score relative to endothelial cells, as a direct contrast
stel <- by_type$sig16[grepl("stellate", by_type$cell.annotation, ignore.case = TRUE)]
endo <- by_type$sig16[grepl("endothelial", by_type$cell.annotation, ignore.case = TRUE)]
if (length(stel) && length(endo)) {
  cat(sprintf("\n16-gene score: stellate %.3f vs endothelial %.3f (ratio %.1f)\n",
              stel[1], endo[1], stel[1] / endo[1]))
}

saveRDS(list(by_type = by_type, specificity = spec,
             sig11_in = sig11_in, sig16_in = sig16_in),
        file.path(DATA_DIR, "derived", "signature_specificity.rds"))

cat("\nNext: 04_single_cell/29_TF_activity.R\n")
