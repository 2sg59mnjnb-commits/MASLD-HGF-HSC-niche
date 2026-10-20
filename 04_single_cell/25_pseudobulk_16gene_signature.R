# ============================================================================
# 25_pseudobulk_16gene_signature.R
#
# Derives the 16-gene HGF+ HSC signature from paired FOV pseudobulk.
#
# Why pseudobulk: at the single-cell level HGF+ HSCs have a higher transcript
# count, so a per-cell comparison is confounded by depth. Aggregating all HGF+
# and all HGF- HSCs within a field of view gives two libraries profiled in the
# same physical neighbourhood, and pairing by FOV removes field-to-field
# variation in tissue quality.
#
# Per slide: sum raw counts per FOV and class (>= 3 cells of each class),
# TMM-normalise, voom, limma with ~ group + pair, extract the HGF+ vs HGF- logFC.
# Across slides: keep genes significant in both with concordant direction.
# HGF itself is excluded, because it defines the classes.
#
# Part 2 repeats the comparison with nCount_RNA as a covariate, which is the
# depth-adjusted validation reported in the supplementary figure.
#
# Output
#   results/04_single_cell/pseudobulk_slide1.csv, pseudobulk_slide2.csv
#   results/04_single_cell/signature_16genes.csv
#   results/04_single_cell/pseudobulk_ncount_covariate.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({
  library(Seurat); library(edgeR); library(limma); library(Matrix)
})

# ---- 1. Pseudobulk per slide ---------------------------------------------
pseudobulk <- function(slide) {
  obj <- load_cosmx(slide)
  ct  <- as.character(obj@meta.data[[celltype_col(slide)]])
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  dat    <- obj[["RNA"]]$data
  hgf    <- as.numeric(dat["HGF", ])
  ncount <- obj@meta.data$nCount_RNA
  fov    <- obj@meta.data$fov

  hsc_idx <- which(ct == "Stellate.cells")
  pos <- hsc_idx[hgf[hsc_idx] > 0]
  neg <- hsc_idx[hgf[hsc_idx] == 0]

  fovs <- sort(unique(fov[hsc_idx]))
  keep_fov <- fovs[vapply(fovs, function(f) {
    sum(fov[pos] == f) >= 3 && sum(fov[neg] == f) >= 3
  }, logical(1))]

  mat_pos <- matrix(0, nrow = nrow(counts), ncol = length(keep_fov),
                    dimnames = list(rownames(counts), paste0("F", keep_fov, "_pos")))
  mat_neg <- mat_pos
  nc_pos <- nc_neg <- numeric(length(keep_fov))

  for (i in seq_along(keep_fov)) {
    f <- keep_fov[i]
    p <- pos[fov[pos] == f]; n <- neg[fov[neg] == f]
    mat_pos[, i] <- Matrix::rowSums(counts[, p, drop = FALSE])
    mat_neg[, i] <- Matrix::rowSums(counts[, n, drop = FALSE])
    nc_pos[i] <- mean(ncount[p]); nc_neg[i] <- mean(ncount[n])
  }
  rm(obj); gc(verbose = FALSE)

  cat(sprintf("Slide%d: %d FOV pairs used\n", slide, length(keep_fov)))
  list(pos = mat_pos, neg = mat_neg, nc_pos = nc_pos, nc_neg = nc_neg)
}

pb1 <- pseudobulk(1)
pb2 <- pseudobulk(2)

# ---- 2. Paired limma per slide -------------------------------------------
fit_slide <- function(pb, label, ncount_covariate = FALSE) {
  mat <- cbind(pb$pos, pb$neg)
  grp <- factor(c(rep("pos", ncol(pb$pos)), rep("neg", ncol(pb$neg))),
                levels = c("neg", "pos"))
  pair <- factor(rep(seq_len(ncol(pb$pos)), 2))
  nc <- c(pb$nc_pos, pb$nc_neg)

  keep <- rowSums(mat >= 5) >= 3
  mat <- mat[keep, ]
  dge <- DGEList(counts = mat)
  dge <- calcNormFactors(dge)

  design <- if (ncount_covariate) {
    model.matrix(~ grp + pair + scale(nc))
  } else {
    model.matrix(~ grp + pair)
  }
  # Drop any aliased columns
  qr_ok <- qr(design)$rank == ncol(design)
  if (!qr_ok) design <- design[, qr(design)$pivot[seq_len(qr(design)$rank)],
                               drop = FALSE]

  v <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = "grppos", number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  tt$slide <- label
  tt
}

deg1 <- fit_slide(pb1, "S1")
deg2 <- fit_slide(pb2, "S2")
write.csv(deg1, res_path("04_single_cell", "pseudobulk_slide1.csv"), row.names = FALSE)
write.csv(deg2, res_path("04_single_cell", "pseudobulk_slide2.csv"), row.names = FALSE)

cat(sprintf("\nSlide1 DEGs: %d (FDR<0.05, |logFC|>0.5)\n",
            sum(deg1$adj.P.Val < 0.05 & abs(deg1$logFC) > 0.5)))
cat(sprintf("Slide2 DEGs: %d\n",
            sum(deg2$adj.P.Val < 0.05 & abs(deg2$logFC) > 0.5)))

# ---- 3. Cross-slide intersection -----------------------------------------
merged <- merge(deg1[, c("gene", "logFC", "adj.P.Val")],
                deg2[, c("gene", "logFC", "adj.P.Val")],
                by = "gene", suffixes = c("_S1", "_S2"))

merged$sig_S1 <- merged$adj.P.Val_S1 < 0.05 & abs(merged$logFC_S1) > 0.5
merged$sig_S2 <- merged$adj.P.Val_S2 < 0.05 & abs(merged$logFC_S2) > 0.5

both <- merged[merged$sig_S1 & merged$sig_S2 &
                 sign(merged$logFC_S1) == sign(merged$logFC_S2), ]
cat(sprintf("\nSignificant in both slides with concordant direction: %d\n", nrow(both)))

# HGF defines the classes, so it must not enter the signature
both <- both[both$gene != "HGF", ]
signature <- both[both$logFC_S1 > 0 & both$logFC_S2 > 0, ]
signature <- signature[order(-(signature$logFC_S1 + signature$logFC_S2)), ]

cat(sprintf("Up-regulated in both slides (signature): %d genes\n", nrow(signature)))
print(signature[, c("gene", "logFC_S1", "logFC_S2")], row.names = FALSE)

cat(sprintf("\nOverlap with the reported 16-gene signature: %d / 16\n",
            length(intersect(signature$gene, SIG16))))
if (length(setdiff(SIG16, signature$gene))) {
  cat("Not recovered here: ", paste(setdiff(SIG16, signature$gene), collapse = ", "), "\n")
}
if (length(setdiff(signature$gene, SIG16))) {
  cat("Additional:       ", paste(setdiff(signature$gene, SIG16), collapse = ", "), "\n")
}

write.csv(signature, res_path("04_single_cell", "signature_16genes.csv"),
          row.names = FALSE)

# ---- 4. Depth-adjusted validation ---------------------------------------
cat("\n=== Pseudobulk with nCount covariate ===\n")
deg1_nc <- fit_slide(pb1, "S1", ncount_covariate = TRUE)
deg2_nc <- fit_slide(pb2, "S2", ncount_covariate = TRUE)

check <- data.frame(
  gene = SIG16,
  logFC_S1 = deg1_nc$logFC[match(SIG16, deg1_nc$gene)],
  logFC_S2 = deg2_nc$logFC[match(SIG16, deg2_nc$gene)],
  P_S1 = deg1_nc$adj.P.Val[match(SIG16, deg1_nc$gene)],
  P_S2 = deg2_nc$adj.P.Val[match(SIG16, deg2_nc$gene)]
)
check$concordant <- sign(check$logFC_S1) == sign(check$logFC_S2) & check$logFC_S1 > 0
check$significant <- check$P_S1 < 0.05 & check$P_S2 < 0.05
print(check, row.names = FALSE)

cat(sprintf("\nDirection preserved: %d / 16\n", sum(check$concordant, na.rm = TRUE)))
cat(sprintf("Significant in both slides: %d / 16\n", sum(check$significant, na.rm = TRUE)))
cat("With ~60 pseudobulk pairs per comparison, gene-level significance is limited\n")
cat("by sample size; the direction of effect is the informative readout here.\n")

write.csv(check, res_path("04_single_cell", "pseudobulk_ncount_covariate.csv"),
          row.names = FALSE)

saveRDS(list(slide1 = deg1, slide2 = deg2, merged = merged,
             signature = signature$gene, ncount_check = check,
             n_pairs = c(ncol(pb1$pos), ncol(pb2$pos))),
        file.path(DATA_DIR, "derived", "pseudobulk.rds"))

cat("\nNext: 04_single_cell/26_rank_matched_validation.R\n")

# ---- 5. GO biological process enrichment of the signature ----------------
# Reported as Figure 3c. Enrichment is against the whole genome background
# (clusterProfiler default), not against the CosMx panel, because the signature
# was derived from panel genes but is interpreted as a transcriptional programme.
if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressPackageStartupMessages({
    library(clusterProfiler); library(org.Hs.eg.db)
  })
  genes_sig <- signature$gene
  ids <- tryCatch(
    bitr(genes_sig, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db),
    error = function(e) NULL
  )
  if (!is.null(ids) && nrow(ids) >= 5) {
    ego <- enrichGO(gene = ids$ENTREZID, OrgDb = org.Hs.eg.db,
                    ont = "BP", pAdjustMethod = "BH",
                    pvalueCutoff = 0.5, readable = TRUE)
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      go <- as.data.frame(ego)
      print(head(go[, c("Description", "GeneRatio", "p.adjust")], 10))
      write.csv(go, res_path("04_single_cell", "signature16_GO_BP.csv"),
                row.names = FALSE)
      saveRDS(go, file.path(DATA_DIR, "derived", "signature16_go.rds"))
    } else {
      cat("No enriched GO terms at FDR < 0.5.\n")
    }
  }
} else {
  cat("clusterProfiler / org.Hs.eg.db not installed; skipping GO enrichment.\n")
}
