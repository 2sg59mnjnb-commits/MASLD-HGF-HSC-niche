# A spatially defined pericentral HGF+ HSC niche adjacent to Hep.4 hepatocytes in F2 MASLD fibrosis

Analysis code for the manuscript *"A spatially defined pericentral HGF+ HSC niche adjacent to Hep.4 hepatocytes in F2 MASLD fibrosis"* (Wu, Guo, Long, Yang).

This repository reproduces every analysis and figure reported in the paper, from raw public
data to final panels. It is organised as a linear pipeline: each script reads the outputs of
the previous stage and writes to `results/` or `figures/`.

---

## 1. What the study shows

Two questions are asked and answered separately.

**Spatial question.** In CosMx spatial transcriptomics of liver biopsies from four patients
with F2 MASLD, HGF-expressing hepatic stellate cells sit closer to a specific hepatocyte
subpopulation (Hep.4) than HGF-negative stellate cells do. The effect is present in all four
patients, survives depth adjustment, and strengthens when the HGF positivity threshold is
raised. Coordinate measurements place Hep.4 hepatocytes **pericentrally**, not periportally,
after correcting for the local density of the reference structures.

**Biomarker question.** An 11-gene tissue signature derived from bulk transcriptomics does
**not** capture this niche. It is not cell-type specific, it does not outperform CHI3L1 alone
in external validation (AUC 0.784 vs 0.802, DeLong P = 0.654), its score is not enriched
around HGF+ HSCs, and re-deriving it after excluding the healthy controls yields a completely
different gene set. The signature is reported as a distributed fibrosis-associated proxy, one
solution among several with equivalent performance.

The code is written to make the negative and model-dependent results as easy to reproduce as
the positive ones.

---

## 2. Data

All data are public. Nothing is redistributed here; the scripts download or expect the files
below.

| Accession | Type | Role in the study | Expected local path |
|---|---|---|---|
| GSE135251 | Bulk RNA-seq, 216 samples (F0–F4) | Discovery / training | `data/GSE135251_counts/` (auto-downloaded) |
| GSE49541 | Affymetrix HG-U133 Plus 2.0, 72 samples | Primary external validation | fetched by `getGEO()` |
| GSE130970 | Bulk RNA-seq, 78 samples | Exploratory validation | `data/GSE130970_TPM.csv.gz` (auto-downloaded) |
| GSE89632 | Illumina HumanHT-12 v4 | Continuous stage trend; incremental value | fetched by `getGEO()` |
| GSE84044 | Affymetrix HG-U133 Plus 2.0, 124 samples | Continuous stage trend | fetched by `getGEO()` |
| GSE292268 | CosMx spatial transcriptomics, 2 slides / 4 patients | Spatial niche | `data/GSM8855042_Slide1_SeuratObject.RDS`, `data/GSM8855043_Slide2_SeuratObject.RDS` |
| GSE202379 | scRNA-seq, human liver | Single-cell reference | `data/GSE202379_SeuratObject_AllCells.rds.gz` |

Manual downloads required (the three Seurat objects are large and are not fetched automatically):

- **GSE292268** — `GSM8855042_Slide1_SeuratObject.RDS` and `GSM8855043_Slide2_SeuratObject.RDS`
- **GSE202379** — `GSE202379_SeuratObject_AllCells.rds.gz`

Optional, only for the exploratory Mendelian randomisation script:

- eQTLGen cis-eQTL summary statistics → `data/MR/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz`
- NAFLD GWAS GCST90091033 → `data/MR/GCST90091033_buildGRCh37.tsv.gz`

---

## 3. Repository layout

```
.
├── 00_setup/            package installation
├── R/                   shared configuration and helper functions
├── 01_discovery/        training cohort: matrix, DEG, WGCNA, 38 candidates, 4 algorithms
├── 02_model/            cross-validation, external validation, calibration, CHI3L1
├── 03_spatial/          CosMx: niche, proximity, zonation, depth, CellChat
├── 04_single_cell/      CosMx HGF+/HGF- DEG, pseudobulk signature, scRNA reference, TF
├── 05_downstream/       stage trends, GSVA, deconvolution, incremental value, drugs, MR
├── 06_figures/          Figure 1-5 and the supplementary figures
├── 07_tables/           supplementary tables
├── data/                inputs and derived intermediates (git-ignored)
├── results/             all tabular outputs (git-ignored)
└── figures/             all figure outputs (git-ignored)
```

`data/`, `results/` and `figures/` are created on first run and excluded from version control;
everything in them can be regenerated.

---

## 4. Quick start

```bash
# 1. install packages (once)
Rscript 00_setup/00_install_packages.R

# 2. place the three manually downloaded Seurat objects in data/ (see section 2)

# 3. run from the repository root, in order
Rscript 01_discovery/01_download_GSE135251.R      # 10-30 min, network
Rscript 01_discovery/02_count_matrix_and_DEG.R    # 10 min
Rscript 01_discovery/03_wgcna.R                   # 20-40 min
Rscript 01_discovery/04_candidate_genes.R
Rscript 01_discovery/05_four_algorithm_selection.R

Rscript 02_model/06_fixed_signature_cv.R
Rscript 02_model/07_fully_nested_cv.R             # 1-3 h
Rscript 02_model/08_external_validation_GSE49541.R
Rscript 02_model/09_calibration_DCA_Firth.R
Rscript 02_model/10_CHI3L1_comparison.R
Rscript 02_model/11_GSE130970_exploratory.R
Rscript 02_model/12_NAFLD_only_rederivation.R     # 30-60 min

Rscript 03_spatial/13_load_CosMx_and_panel.R
Rscript 03_spatial/14_stellate_fibrotic_enrichment.R
Rscript 03_spatial/15_HGFpos_niche_scan.R
Rscript 03_spatial/16_HGFpos_proximity_mixed_models.R
Rscript 03_spatial/17_FOV_consistency_and_thresholds.R
Rscript 03_spatial/18_zonation_portal_vs_central.R
Rscript 03_spatial/19_Hep4_IL6R_and_pathways.R
Rscript 03_spatial/20_depth_correction.R
Rscript 03_spatial/21_distance_gradient.R
Rscript 03_spatial/22_selectivity_multiple_testing.R
Rscript 03_spatial/23_CellChat_spatial.R          # optional, 30-60 min
Rscript 03_spatial/24_HGFpos_neighbourhood_signature.R

Rscript 04_single_cell/24_cosmx_HGFpos_vs_HGFneg_DEG.R
Rscript 04_single_cell/25_pseudobulk_16gene_signature.R
Rscript 04_single_cell/26_rank_matched_validation.R
Rscript 04_single_cell/27_Hep4_identity_scRNA.R
Rscript 04_single_cell/28_signature_specificity_scRNA.R
Rscript 04_single_cell/29_TF_activity.R
Rscript 04_single_cell/30_TF_HGF_coupling_and_motif.R

Rscript 05_downstream/31_stage_trends.R
Rscript 05_downstream/32_GSVA.R
Rscript 05_downstream/33_deconvolution_NNLS.R
Rscript 05_downstream/34_multi_cohort_trends.R
Rscript 05_downstream/35_incremental_value.R
Rscript 05_downstream/36_drug_repurposing.R
Rscript 05_downstream/37_MR_exploratory.R         # optional, not in the paper

Rscript 06_figures/38_Figure1.R
Rscript 06_figures/39_Figure2.R
Rscript 06_figures/40_Figure3.R
Rscript 06_figures/41_Figure4_Figure5.R
Rscript 06_figures/42_supplementary_figures.R
Rscript 07_tables/43_supplementary_tables.R
```

Total runtime is roughly 4–6 h on a laptop with 16 GB RAM; the two WGCNA steps and the nested
cross-validation dominate.

Each script prints a short observed-versus-reported comparison for the headline numbers, so a
run can be checked against the paper without opening the outputs.

---

## 5. Figure to script map

| Figure | Panels | Produced by | Output |
|---|---|---|---|
| Figure 1 | cell composition; stellate spatial map; enrichment OR | `03_spatial/13`, `14` then `06_figures/38` | `figures/Figure1.pdf` |
| Figure 2 | niche scan; FOV consistency; cross-slide; threshold; HGF counts; SPI; per-patient; mixed models | `03_spatial/15`, `16`, `17` then `06_figures/39` | `figures/Figure2.pdf` |
| Figure 3 | pseudobulk; positive rates; GO; scRNA score | `04_single_cell/25`, `24`, `28` then `06_figures/40` | `figures/Figure3.pdf` |
| Figure 4 | IL6-family receptor positive rates | `04_single_cell/24` then `06_figures/41` | `figures/Figure4.pdf` |
| Figure 5 | training ROC; external ROC; decision curve | `02_model/06`, `08`, `09` then `06_figures/41` | `figures/Figure5.pdf` |

Supplementary figures are built by `06_figures/42_supplementary_figures.R`, which prints any
panel it had to skip because an upstream output was missing.

| Supplementary | Content | Upstream script |
|---|---|---|
| S1 | WGCNA module–trait correlation | `01_discovery/03` |
| S2 | four-algorithm selection | `01_discovery/05` |
| S3 | calibration in GSE49541 | `02_model/09` |
| S4 | fixed vs fully nested CV | `02_model/06`, `07` |
| S7 | inter-sample heterogeneity | `03_spatial/16` |
| S9 | CHI3L1 versus multi-gene models | `02_model/10` |
| S11 | bootstrap metric uncertainty | `02_model/09` |
| S12 | F0–F4 signature trend | `05_downstream/31` |
| S15 | Firth versus MLE | `02_model/09` |
| S26 | HGF threshold sensitivity | `03_spatial/17` |
| S27 | CHI3L1 ROC overlay | `02_model/10` |
| S30 | signature specificity | `04_single_cell/28` |
| S31 | IL6R model sensitivity | `03_spatial/19` |
| S32 | depth-corrected distance | `03_spatial/20` |
| S33 | TF activity | `04_single_cell/29` |
| S34 | motif enrichment | `04_single_cell/30` |
| S35 | drug repurposing | `05_downstream/36` |

---

## 6. Reproducibility notes

- **R version.** Developed under R 4.6.x. Seurat v5 is required for the spatial and
  single-cell scripts (`layer` argument). `sessionInfo()` after a full run is written to
  `sessionInfo.txt`.
- **Seeds.** All stochastic steps use a single seed (`SEED <- 123`, `R/config.R`) except where
  a different seed is set deliberately inside a fold. Reported values assume this seed.
- **Stochastic steps.** WGCNA module assignment, cross-validation fold composition and the
  permutation tests are all seed-dependent. Small differences from the published numbers are
  expected if the seed or the package versions change; the direction and significance of the
  results are stable.
- **Nested cross-validation optimism.** `02_model/07` reuses the soft-threshold power selected
  on the full cohort (`POWER_MODE <- "full"`), which is what was run for the paper. This is a
  mild optimism because one hyperparameter is chosen once on all samples. Setting
  `POWER_MODE <- "in_fold"` re-selects it inside every fold, which is stricter and slower.
- **Cross-platform scaling.** Validation cohorts are z-scored within cohort, which is the
  harmonisation described in Methods and the one that gives the reported GSE49541 AUC of 0.784.
  `02_model/08_external_validation_GSE49541.R` additionally prints a sensitivity value obtained
  by scaling with the training parameters; that variant is not the one reported.
- **Depth.** HGF+ and HGF- HSCs differ in transcript count (median nCount ≈288 vs ≈114) and
  the distributions barely overlap. Quantile matching equalises rank, not absolute depth, so
  the distance analyses use covariate adjustment and trimming instead
  (`03_spatial/20`). Scripts state which strategy each result depends on.

---

## 7. Scope, and what is deliberately not here

The repository contains the analyses reported in the paper plus two exploratory analyses that
are **not** reported, kept so the negative results are on record rather than silently dropped:

- `05_downstream/37_MR_exploratory.R` — cis-Mendelian randomisation for CHI3L1. The cis-eQTLs
  are in strong linkage disequilibrium; after 1 Mb pruning one instrument remains and the Wald
  ratio is not significant. It cannot distinguish a null effect from insufficient power.
- The scTenifold virtual-knockout analysis run during revision is not included: it did not
  reach a stable result across parameter settings and is not used in the paper.

Two supplementary tables referenced in the manuscript are not regenerated here:

- **Table S13 (TMM normalisation sensitivity)** — not rerun in this release.
- **Table S14 (MET pathway sensitivity)** — the MET pathway analysis was excluded from the
  final manuscript because the apparent signal did not survive depth matching, so the code is
  not carried forward.

Limitations that the code makes visible:

- **Four patients.** With n = 4 the smallest attainable one-sided binomial P is 0.0625. The
  per-patient consistency (4/4, differences of −20 to −27 µm) is the informative quantity, and
  the mixed-model P values are not patient-level evidence.
- **CellChat and the panel.** The 1,000-plex panel carries 5/7 IL6-family ligands and 6/6
  receptors. A negative CellChat result for a missing ligand is a detection limit, not
  biological absence. `03_spatial/13` prints the coverage table that must accompany any such
  claim.
- **Two genes for the neighbourhood control.** Only CHI3L1 and GSN of the 11-gene signature are
  on the panel, so `03_spatial/24` tests a two-gene proxy, not the full signature.
- **Transcriptional versus anatomical zonation.** The Hep.4-equivalent scRNA cluster scores
  highest on Zone 1/Zone 2 markers, while the coordinate measurements place Hep.4 pericentrally.
  `04_single_cell/27` reports both rather than reconciling them; the manuscript gives precedence
  to the direct coordinate measurement.

---

## 8. Citation

If you use this code, please cite the manuscript (see `CITATION.cff`).

## 9. Licence

MIT — see `LICENSE`.

## 10. Contact

Hongju Yang (corresponding author), Department of Geriatrics, Clinical Medicine Research Center,
The First Affiliated Hospital of Kunming Medical University, Kunming, Yunnan, China.
E-mail: hongjuyang105@163.com

---

## 中文快速说明

1. 先运行 `Rscript 00_setup/00_install_packages.R` 安装 R 包。
2. 手动下载三个 Seurat 对象（GSE292268 的两张切片、GSE202379 的单细胞对象），放进 `data/`；
   其余数据由脚本自动从 GEO 下载。
3. 在仓库根目录按第 4 节的顺序依次运行脚本，结果写入 `results/`，图写入 `figures/`。
4. 每个脚本运行时会打印"实际值 vs 论文值"的对照，便于核对，不需要手工翻结果表。
5. 参考基因组/坐标相关说明：CosMx 的空间坐标单位是 mm，脚本内部统一换算成 µm；
   CellChat 用 0.18 µm/pixel 的换算系数。
6. 不需要跑的可选步骤：`03_spatial/23_CellChat_spatial.R`（需另外安装 CellChat）和
   `05_downstream/37_MR_exploratory.R`（论文未采用，仅作为阴性结果存档）。
