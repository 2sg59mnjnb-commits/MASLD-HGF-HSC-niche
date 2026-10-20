# ============================================================================
# 35_incremental_value.R
#
# Does the signature add anything to routine clinical variables? This can only
# be asked in a cohort with paired clinical annotation, which in this study is
# GSE89632 (n = 53 with complete data, 9 events).
#
# Three models are compared on the same samples:
#   clinical only      age, sex, BMI, AST, ALT
#   signature only     the 11-gene score
#   combined           both
#
# Reported: AUC for each model, DeLong for combined vs clinical, likelihood-ratio
# test, continuous NRI and IDI, and the events-per-variable ratio.
#
# The expected answer is that the signature does NOT add measurable value, and
# that is what the manuscript reports. The EPV here is below 1, so the combined
# model is uninterpretable in any strong sense; this is stated in the output.
#
# Output
#   results/05_downstream/incremental_value_GSE89632.csv
#   data/derived/incremental_value.rds
# ============================================================================

source("R/config.R")
source("R/helpers.R")

suppressPackageStartupMessages({
  library(GEOquery); library(pROC); library(illuminaHumanv4.db)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

# ---- 1. Load cohort ------------------------------------------------------
f <- file.path(DATA_DIR, "public_cohorts", "GSE89632.rds")
if (!file.exists(f)) {
  dir.create(file.path(DATA_DIR, "public_cohorts", "GSE89632"),
             showWarnings = FALSE, recursive = TRUE)
  gse <- getGEO("GSE89632", GSEMatrix = TRUE, getGPL = FALSE,
                destdir = file.path(DATA_DIR, "public_cohorts", "GSE89632"))
  g <- list(expr = exprs(gse[[1]]), pheno = pData(gse[[1]]))
  saveRDS(g, f)
} else {
  g <- readRDS(f)
}
pheno <- g$pheno

# ---- 2. Expression and signature ----------------------------------------
sym <- AnnotationDbi::mapIds(illuminaHumanv4.db, keys = rownames(g$expr),
                             column = "SYMBOL", keytype = "PROBEID",
                             multiVals = "first")
ok <- !is.na(sym)
df <- as.data.frame(g$expr[ok, , drop = FALSE], check.names = FALSE)
df$symbol <- sym[ok]
expr <- df %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean), .groups = "drop") %>%
  as.data.frame()
rownames(expr) <- expr$symbol
expr$symbol <- NULL
expr <- as.matrix(expr)

sig <- sig_score(expr, SIG11)

# ---- 3. Parse clinical variables ----------------------------------------
# All characteristics columns are concatenated and parsed with a single
# key:value pattern per variable, which is robust to the column layout.
ch_cols <- grep("^characteristics", colnames(pheno), value = TRUE)
combined_text <- apply(pheno[, ch_cols, drop = FALSE], 1, function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  paste(x, collapse = "; ")
})

# regmatches() on a plain regexpr() drops non-matching rows, which silently
# shortens the vector; positions are preserved explicitly here.
parse_num <- function(text, key) {
  pat <- paste0("(?:^|;\\s*)", key, "[^;:]*:\\s*([-0-9.,]+)")
  m <- regexpr(pat, text, ignore.case = TRUE, perl = TRUE)
  v <- rep(NA_character_, length(text))
  hit <- m > 0
  v[hit] <- sub(pat, "\\1", regmatches(text, m), perl = TRUE)
  v <- gsub(",", "", v)
  suppressWarnings(as.numeric(v))
}

age  <- parse_num(combined_text, "age")
bmi  <- parse_num(combined_text, "body mass index")
ast  <- parse_num(combined_text, "aspartate transaminase")
alt  <- parse_num(combined_text, "alanine transaminase")
fib  <- parse_num(combined_text, "fibrosis")
sex  <- ifelse(grepl("gender\\s*[:=]\\s*male", combined_text, ignore.case = TRUE), 1,
               ifelse(grepl("gender\\s*[:=]\\s*female", combined_text, ignore.case = TRUE), 0, NA))

cat(sprintf("Parsed: age %d, BMI %d, AST %d, ALT %d, sex %d, stage %d\n",
            sum(!is.na(age)), sum(!is.na(bmi)), sum(!is.na(ast)),
            sum(!is.na(alt)), sum(!is.na(sex)), sum(!is.na(fib))))

dataset <- data.frame(
  sample = colnames(expr), signature = sig[colnames(expr)],
  age = age, bmi = bmi, ast = ast, alt = alt, sex = sex,
  stage = fib, stringsAsFactors = FALSE
)
dataset$fibrosis_binary <- ifelse(dataset$stage >= 2, 1,
                                  ifelse(dataset$stage <= 1, 0, NA))

need <- c("signature", "age", "sex", "bmi", "ast", "alt", "fibrosis_binary")
d <- dataset[complete.cases(dataset[, need]), ]
cat(sprintf("\nComplete cases: %d, events: %d\n", nrow(d), sum(d$fibrosis_binary)))

if (nrow(d) < 20 || sum(d$fibrosis_binary) < 5 || sum(1 - d$fibrosis_binary) < 5) {
  cat("Too few complete cases for a multivariable comparison; stopping.\n")
  quit(save = "no", status = 0)
}

for (v in c("age", "bmi", "ast", "alt", "signature")) {
  d[[paste0(v, "_z")]] <- as.numeric(scale(d[[v]]))
}

# ---- 4. Three models -----------------------------------------------------
m_clin <- glm(fibrosis_binary ~ age_z + sex + bmi_z + ast_z + alt_z,
              data = d, family = binomial)
m_sig  <- glm(fibrosis_binary ~ signature_z, data = d, family = binomial)
m_comb <- glm(fibrosis_binary ~ signature_z + age_z + sex + bmi_z + ast_z + alt_z,
              data = d, family = binomial)

roc_clin <- roc(d$fibrosis_binary, predict(m_clin, type = "response"), quiet = TRUE)
roc_sig  <- roc(d$fibrosis_binary, predict(m_sig,  type = "response"), quiet = TRUE)
roc_comb <- roc(d$fibrosis_binary, predict(m_comb, type = "response"), quiet = TRUE)

aucs <- data.frame(
  model = c("Clinical only", "Signature only", "Clinical + signature"),
  auc = c(as.numeric(auc(roc_clin)), as.numeric(auc(roc_sig)), as.numeric(auc(roc_comb))),
  ci_low  = c(ci.auc(roc_clin)[1], ci.auc(roc_sig)[1], ci.auc(roc_comb)[1]),
  ci_high = c(ci.auc(roc_clin)[3], ci.auc(roc_sig)[3], ci.auc(roc_comb)[3])
)
aucs[, 2:4] <- round(aucs[, 2:4], 3)
print(aucs)

dl <- roc.test(roc_comb, roc_clin, method = "delong")
lrt <- anova(m_clin, m_comb, test = "LRT")

cat(sprintf("\nDeLong combined vs clinical: P = %.3f\n", dl$p.value))
cat(sprintf("Likelihood-ratio test: chi2 = %.3f, P = %.3f\n",
            lrt$Deviance[2], lrt$`Pr(>Chi)`[2]))

# NRI (continuous) and IDI
nri <- function(p_old, p_new, y) {
  up <- p_new > p_old; down <- p_new < p_old
  (mean(up[y == 1]) - mean(down[y == 1])) +
    (mean(down[y == 0]) - mean(up[y == 0]))
}
p_clin <- predict(m_clin, type = "response")
p_comb <- predict(m_comb, type = "response")
nri_val <- nri(p_clin, p_comb, d$fibrosis_binary)
idi_val <- (mean(p_comb[d$fibrosis_binary == 1]) - mean(p_comb[d$fibrosis_binary == 0])) -
           (mean(p_clin[d$fibrosis_binary == 1]) - mean(p_clin[d$fibrosis_binary == 0]))
epv <- sum(d$fibrosis_binary) / 11

cat(sprintf("NRI = %.3f, IDI = %.3f, EPV = %.2f\n", nri_val, idi_val, epv))
cat("EPV < 1 means the combined model is not estimable in any strong sense; the\n")
cat("comparison is reported as exploratory only.\n")

# ---- 5. Save -------------------------------------------------------------
out <- data.frame(
  model = aucs$model, auc = aucs$auc, ci_low = aucs$ci_low, ci_high = aucs$ci_high,
  n = nrow(d), events = sum(d$fibrosis_binary),
  delong_combined_vs_clinical = round(dl$p.value, 3),
  lrt_chi2 = round(lrt$Deviance[2], 3),
  lrt_P = round(lrt$`Pr(>Chi)`[2], 3),
  NRI = round(nri_val, 3), IDI = round(idi_val, 3), EPV = round(epv, 2)
)
write.csv(out, res_path("05_downstream", "incremental_value_GSE89632.csv"),
          row.names = FALSE)

saveRDS(list(aucs = aucs, delong = dl, lrt = lrt, nri = nri_val, idi = idi_val,
             epv = epv, data = d),
        file.path(DATA_DIR, "derived", "incremental_value.rds"))

cat("\nNext: 05_downstream/36_drug_repurposing.R\n")
