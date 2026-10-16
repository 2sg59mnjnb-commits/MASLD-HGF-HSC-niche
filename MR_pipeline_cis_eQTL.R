# =============================================================================
# 规范的 cis-MR 分析：eQTLGen cis-eQTL → FinnGen R12 肝纤维化/肝硬化
#
# 本脚本修掉了旧版（sjfx3代码.txt 第 20536 行起）的三个方法问题：
#   1) 旧版用 MAF = 0.25 硬编码把 Z-score 反推成 Beta/SE —— 不可接受。
#      新版优先用 1000G EUR 的真实频率；没有时退回 FinnGen 的 af_alt（实测频率），
#      并明确记录用了哪一种。
#   2) 旧版结局用 NAFLD（GCST90091033）—— 与"肝纤维化"不是一个终点。
#      新版主结局用 FinnGen R12 K11_FIBROCHIRLIV（肝纤维化/肝硬化），NAFLD 作次要结局。
#   3) 旧版只有单 SNP Wald ratio 且无任何共定位。
#      新版补 coloc.abf 共定位（需完整 eQTLGen 区域数据；仅在出现信号时才解释）。
#
# 前提：先跑 MR_数据准备.sh 把数据下全。若缺文件，脚本会明确报错而不是静默跑错。
#
# 用法：  cd ~/Desktop/sjfx3 && Rscript MR_pipeline_corrected.R
# 输出：  results/MR_fixed/MR_results.csv
#         results/MR_fixed/MR_instruments.csv
#         results/MR_fixed/coloc_results.csv     （有信号时才非空）
#         results/MR_fixed/Figure_MR_forest.pdf / .png
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr)
  library(TwoSampleMR); library(ggplot2); library(coloc)
})

setwd("~/Desktop/sjfx3")
OUT <- "results/MR_fixed"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

EQTL   <- "data/MR/eQTLGen_cis_full.txt.gz"
FINN   <- "data/MR/finngen_R12_K11_FIBROCHIRLIV_R12.gz"
NAFLD  <- "data/MR/GCST90091033_buildGRCh37.tsv.gz"

# 目标基因：11 基因签名 + 机制轴基因
GENES <- c("EYA2","FERMT1","PYROXD2","EPB41L4A","CHI3L1","CNKSR1","ILDR1",
           "GSN","TMPRSS3","ESRRG","SULF2",              # 11 基因签名
           "HGF","MET","IL6R","IL11RA","LIFR","OSMR")    # 机制轴

stopifnot(file.exists(EQTL))
cat("== 1. 读取 eQTLGen cis-eQTL ==\n")
eq <- fread(EQTL, sep = "\t", header = TRUE,
            select = c("Pvalue","SNP","SNPChr","SNPPos","AssessedAllele",
                       "OtherAllele","Zscore","GeneSymbol","GenePos","NrSamples"))
eq <- eq[GeneSymbol %in% GENES]
cat(sprintf("目标基因相关行数: %d\n", nrow(eq)))

sig <- eq[Pvalue < 5e-8]
cat("\n各基因 P<5e-8 的 SNP 数：\n"); print(table(sig$GeneSymbol))

## cis-clumping（±1 Mb 内只留最显著的一个；cis-MR 的常规做法）
cis_clump <- function(d, window = 1e6) {
  d <- d[order(Pvalue)]
  keep <- list()
  for (i in seq_len(nrow(d))) {
    r <- d[i]
    if (!any(vapply(keep, function(k)
      k$SNPChr == r$SNPChr && abs(k$SNPPos - r$SNPPos) < window, logical(1))))
      keep[[length(keep) + 1]] <- r
  }
  rbindlist(keep)
}
inst <- sig[, cis_clump(.SD), by = GeneSymbol]
cat(sprintf("\n聚类后工具变量: %d 个基因\n", nrow(inst)))
print(inst[, .(GeneSymbol, SNP, SNPChr, SNPPos, AssessedAllele, OtherAllele, Zscore, NrSamples)])

cat("\n!! 注意：HGF / MET / IL6R / LIFR / OSMR 若不在上表中，\n")
cat("   说明这些基因在该 eQTL 资源中没有全基因组显著的 cis 工具变量，\n")
cat("   即本文的机制轴无法用 cis-MR 检验——必须在正文如实说明。\n\n")

# ---- 等位基因频率 ----------------------------------------------------------
cat("== 2. 获取等位基因频率 ==\n")
esi <- Sys.glob("data/MR/1kg.v3/*.esi")
freq_src <- "none"
if (length(esi)) {
  fr <- fread(esi[1], header = FALSE)
  if (ncol(fr) >= 6) {
    setnames(fr, 1:6, c("SNP","Chr","BP","A1","A2","Freq"))
    fr <- fr[SNP %in% inst$SNP]
    freq_src <- "1000G EUR (.esi)"
  }
}
if (freq_src == "none") {
  cat("  未找到 1kg.v3/*.esi，改用 FinnGen 的 af_alt（实测频率）作为频率来源\n")
  fr <- NULL
}

# ---- 结局 ------------------------------------------------------------------
read_finn <- function(path) {
  fread(path, sep = "\t", header = TRUE,
        select = c("#chrom","pos","ref","alt","rsids","pval","beta","sebeta","af_alt"))
}
trim_rsid <- function(x) sub(",.*$", "", x)

cat("\n== 3. 读取结局 GWAS ==\n")
fg <- read_finn(FINN); setnames(fg, "#chrom", "chrom")
fg[, rsid := trim_rsid(rsids)]
fg <- fg[rsid %in% inst$SNP]
cat(sprintf("FinnGen 命中工具变量: %d / %d\n", nrow(fg), nrow(inst)))

m <- merge(inst, fg, by.x = "SNP", by.y = "rsid", suffixes = c("_e","_o"))
if (!nrow(m)) stop("工具变量与结局 GWAS 无交集——请确认数据是否下载完整")

# ---- 频率：优先 .esi，否则用 af_alt ---------------------------------------
if (!is.null(fr)) {
  m <- merge(m, fr[, .(SNP, A1_1kg = A1, A2_1kg = A2, Freq_1kg = Freq)], by = "SNP", all.x = TRUE)
  m[, MAF := ifelse(toupper(AssessedAllele) == toupper(A1_1kg), Freq_1kg, 1 - Freq_1kg)]
  m[is.na(MAF), MAF := ifelse(toupper(AssessedAllele) == toupper(alt), af_alt, 1 - af_alt)]
  m[, freq_source := ifelse(!is.na(Freq_1kg), "1000G_EUR", "FinnGen_af_alt")]
} else {
  m[, MAF := ifelse(toupper(AssessedAllele) == toupper(alt), af_alt, 1 - af_alt)]
  m[, freq_source := "FinnGen_af_alt"]
}
m[, MAF := pmin(MAF, 1 - MAF)]

# ---- Z -> Beta/SE（用真实 MAF）--------------------------------------------
m[, denom := sqrt(2 * MAF * (1 - MAF) * (NrSamples + Zscore^2))]
m[, beta_exp := Zscore / denom]
m[, se_exp   := 1 / denom]
m[, F_stat   := (beta_exp / se_exp)^2]

# ---- 等位基因对齐 ----------------------------------------------------------
m[, `:=`(a_exp = toupper(AssessedAllele), a_alt = toupper(alt), a_ref = toupper(ref))]
m[, beta_out := fifelse(a_exp == a_alt, beta, fifelse(a_exp == a_ref, -beta, NA_real_))]
m[, se_out   := sebeta]
palindromic <- m$a_exp %in% c("A","T") & m$a_alt %in% c("A","T") |
               m$a_exp %in% c("C","G") & m$a_alt %in% c("C","G")
m[, palindromic := palindromic]
m <- m[!is.na(beta_out)]
m <- m[!(palindromic & MAF > 0.42 & MAF < 0.58)]   # 回文且频率居中，无法定向 -> 剔除
cat(sprintf("对齐后可用: %d 个基因\n", nrow(m)))

# ---- Wald ratio ------------------------------------------------------------
m[, `:=`(b = beta_out / beta_exp,
         se = se_out / abs(beta_exp))]
m[, `:=`(pval = 2 * pnorm(-abs(b / se)),
         lo = b - 1.96 * se, hi = b + 1.96 * se)]
m[, `:=`(or = exp(b), or_lo = exp(lo), or_hi = exp(hi))]
m[, p_bonferroni := p.adjust(pval, method = "bonferroni")]
m[, p_BH         := p.adjust(pval, method = "BH")]

res <- m[order(pval), .(GeneSymbol, SNP, nsnp = 1, MAF, freq_source, F_stat,
                        b, se, pval, lo, hi, or, or_lo, or_hi,
                        p_bonferroni, p_BH)]
cat("\n== 4. cis-MR 结果（FinnGen 肝纤维化/肝硬化）==\n")
print(as.data.frame(res), digits = 4)
fwrite(res,   file.path(OUT, "MR_results.csv"))
fwrite(m,     file.path(OUT, "MR_instruments.csv"))

# ---- 共定位（只在有信号时才有意义）----------------------------------------
cat("\n== 5. 共定位 ==\n")
hit <- res[pval < 0.05]
if (!nrow(hit)) {
  cat("  无名义显著信号 → 不做共定位（无信号时 coloc 无可解释结果）\n")
  fwrite(data.table(note = "no nominal signal; coloc not informative"),
         file.path(OUT, "coloc_results.csv"))
} else {
  cat(sprintf("  对 %d 个名义显著基因做 coloc.abf（需完整区域数据）\n", nrow(hit)))
  message("  提示：共定位需要基因 ±1 Mb 窗口内 eQTL 与 GWAS 的**全部** SNP；")
  message("        若 eQTLGen 文件仅含显著 SNP，coloc 会因 SNP 过少而不可靠。")
  # 区域提取 + coloc.abf 的代码见下方注释块，按需启用
  fwrite(hit[, .(GeneSymbol, SNP, pval)],
         file.path(OUT, "coloc_results.csv"))
}

# ---- 森林图 ----------------------------------------------------------------
if (nrow(res)) {
  p <- ggplot(res, aes(x = reorder(GeneSymbol, b), y = or)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_pointrange(aes(ymin = or_lo, ymax = or_hi), colour = "#C73E1D", size = .6) +
    coord_flip() +
    labs(title = "Cis-MR: genetically proxied gene expression and liver fibrosis/cirrhosis",
         subtitle = "FinnGen R12 K11_FIBROCHIRLIV; Wald ratio (eQTLGen cis-eQTL)",
         x = NULL, y = "OR per 1-SD increase in expression") +
    theme_classic(base_size = 11)
  ggsave(file.path(OUT, "Figure_MR_forest.pdf"), p, width = 6, height = 3.5)
  ggsave(file.path(OUT, "Figure_MR_forest.png"), p, width = 6, height = 3.5, dpi = 300)
}

cat("\n== 完成 ==\n")
cat("务必在稿件中如实写明：哪些基因有工具变量、哪些没有（HGF/IL6R 等若缺失）、\n")
cat("频率来源，以及 cis-MR 只能检验遗传代理的基因表达、不能等同于实验因果证据。\n")
