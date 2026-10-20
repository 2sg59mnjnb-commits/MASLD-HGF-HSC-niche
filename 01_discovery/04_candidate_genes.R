# ============================================================================
# 04_candidate_genes.R
#
# Builds the candidate gene set for feature selection by intersecting the
# fibrosis-associated WGCNA modules with the DEG list, then restricting to
# protein-coding genes with a valid gene symbol.
#
# Input   data/derived/WGCNA_results.rds, data/derived/GSE135251_DEG_results.rds
# Output  results/01_discovery/candidate_genes_filtered.csv  (38 genes)
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

wgcna <- load_derived("WGCNA_results.rds")
deg   <- load_bulk()$deg_sig

module_genes <- unique(c(wgcna$genes_positive, wgcna$genes_negative))
cat(sprintf("%s genes: %d\n", wgcna$top_positive, length(wgcna$genes_positive)))
cat(sprintf("%s genes: %d\n", wgcna$top_negative, length(wgcna$genes_negative)))
cat(sprintf("Union of the two modules: %d\n", length(module_genes)))

deg_genes <- rownames(deg)
cat(sprintf("DEGs: %d\n", length(deg_genes)))

# ---- 1. Intersect with DEGs ------------------------------------------------
candidates <- intersect(module_genes, deg_genes)
cat(sprintf("Module n DEG: %d\n", length(candidates)))

# ---- 2. Ensembl -> symbol --------------------------------------------------
symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = candidates,
                                 column = "SYMBOL", keytype = "ENSEMBL",
                                 multiVals = "first")
map <- data.frame(Ensembl = candidates, Symbol = symbols, stringsAsFactors = FALSE)
map <- map[!is.na(map$Symbol), ]
cat(sprintf("Mapped to a symbol: %d\n", nrow(map)))

# ---- 3. Drop non-coding genes ---------------------------------------------
# The candidate list is restricted to protein-coding genes because the
# downstream models are interpreted at the level of protein-coding transcripts.
info <- AnnotationDbi::select(org.Hs.eg.db, keys = map$Symbol,
                              columns = c("SYMBOL", "GENETYPE"),
                              keytype = "SYMBOL")
info <- info[!duplicated(info$SYMBOL), ]
coding <- info$SYMBOL[info$GENETYPE == "protein-coding"]

filtered <- map[map$Symbol %in% coding, ]
dropped  <- setdiff(map$Symbol, coding)
cat(sprintf("Protein-coding candidates: %d\n", nrow(filtered)))
if (length(dropped)) {
  cat(sprintf("Dropped as non-coding (%d): %s\n",
              length(dropped), paste(dropped, collapse = ", ")))
}

# ---- 4. Save ---------------------------------------------------------------
write.csv(filtered, res_path("01_discovery", "candidate_genes_filtered.csv"),
          row.names = FALSE)
cat(sprintf("\nCandidate set: %d genes\n", nrow(filtered)))
compare("candidate gene count", nrow(filtered), REPORTED$n_candidates, digits = 0)

cat("\nNext: 01_discovery/05_four_algorithm_selection.R\n")
