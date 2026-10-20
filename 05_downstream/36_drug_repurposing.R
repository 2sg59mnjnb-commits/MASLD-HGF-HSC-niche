# ============================================================================
# 36_drug_repurposing.R
#
# Exploratory drug-gene interaction query against DGIdb for the 16-gene HGF+ HSC
# signature.
#
# The 11-gene tissue signature is deliberately not queried: it shares no genes
# with the HSC signature, so a drug list derived from it would not be describing
# the same biology.
#
# No interaction-score threshold is applied. The compounds discussed in the
# manuscript were selected from the returned table by mechanism (retinoid/RXR
# signalling via RXRA/RXRB; FGFR signalling via FGFR2), which is a post hoc
# selection and is labelled as such.
#
# Output
#   results/05_downstream/DGIdb_interactions.csv
#   results/05_downstream/DGIdb_drug_summary.csv
# ============================================================================

source("R/config.R")
source("R/helpers.R")

if (!requireNamespace("httr", quietly = TRUE) ||
    !requireNamespace("jsonlite", quietly = TRUE)) {
  cat("httr / jsonlite not installed; skipping.\n")
  quit(save = "no", status = 0)
}
suppressPackageStartupMessages({library(httr); library(jsonlite)})

genes <- SIG16
cat(sprintf("Querying DGIdb for %d genes\n", length(genes)))

gene_list <- paste(sprintf('"%s"', genes), collapse = ", ")
query <- sprintf('{
  genes(names: [%s]) {
    nodes {
      name
      interactions {
        drug { name approved }
        interactionScore
        interactionTypes { type }
      }
    }
  }
}', gene_list)

res <- tryCatch(
  httr::POST("https://dgidb.org/api/graphql",
             body = list(query = query), encode = "json",
             httr::timeout(60)),
  error = function(e) NULL
)
if (is.null(res) || httr::status_code(res) != 200) {
  cat("DGIdb request failed; no output written.\n")
  quit(save = "no", status = 0)
}

dat <- jsonlite::fromJSON(httr::content(res, "text", encoding = "UTF-8"))
nodes <- dat$data$genes$nodes
if (is.null(nodes) || length(nodes) == 0) {
  cat("No genes returned by DGIdb.\n"); quit(save = "no", status = 0)
}

rows <- list()
for (i in seq_len(nrow(nodes))) {
  g <- nodes$name[i]
  inter <- nodes$interactions[[i]]
  if (is.null(inter) || length(inter) == 0) next
  for (j in seq_len(nrow(inter))) {
    types <- inter$interactionTypes[[j]]
    rows[[length(rows) + 1]] <- data.frame(
      gene = g,
      drug = inter$drug$name[j],
      approved = inter$drug$approved[j],
      score = inter$interactionScore[j],
      interaction_type = if (length(types) && nrow(types)) types$type[1] else NA_character_,
      stringsAsFactors = FALSE
    )
  }
}
interactions <- do.call(rbind, rows)
if (is.null(interactions)) {
  cat("No interactions parsed.\n"); quit(save = "no", status = 0)
}

interactions <- interactions[order(-interactions$score), ]
cat(sprintf("\nInteractions: %d involving %d drugs\n",
            nrow(interactions), length(unique(interactions$drug))))
write.csv(interactions, res_path("05_downstream", "DGIdb_interactions.csv"),
          row.names = FALSE)

approved <- interactions[interactions$approved %in% c(TRUE, "TRUE", "true"), ]
summary_tab <- approved %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(n_genes = dplyr::n(),
                   genes = paste(gene, collapse = ", "),
                   max_score = max(score, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_genes), dplyr::desc(max_score))
print(head(as.data.frame(summary_tab), 20), row.names = FALSE)
write.csv(summary_tab, res_path("05_downstream", "DGIdb_drug_summary.csv"),
          row.names = FALSE)

cat("\nNote: this is a hypothesis-generating query. Mechanism-based selection of\n")
cat("the retinoid and FGFR compounds was made after inspecting the table and is\n")
cat("reported as such; the query itself applies no score threshold.\n")

saveRDS(list(interactions = interactions, summary = summary_tab),
        file.path(DATA_DIR, "derived", "dgidb.rds"))

cat("\nNext: 05_downstream/37_MR_exploratory.R\n")
