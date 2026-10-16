# A spatially defined HGF-expressing HSC niche in MASLD fibrosis

Analysis code for the manuscript by Wu & Yang.

## Repository contents

| File | Contents |
|---|---|
| `analysis_code.R` | The **complete** analysis code for the study, as run. A single consolidated script organised into numbered modules (Module 0 onward) that follow the order of the Methods section. |
| `00_install_packages.R` | Package installation |
| `01_project_config.R` | Central configuration (paths, seed, parameters, colours) |
| `Figure2_rebuild.R` | Rebuilds **main Figure 2** (standalone; see below) |
| `Figure3_redraw.R` | Rebuilds **main Figure 3** and **Additional file 1, Figure S27** (standalone) |
| `MR_data_preparation.sh`, `MR_pipeline_cis_eQTL.R` | cis-Mendelian randomisation analysis — **carried out but not reported**; see `NOTES_exploratory_analyses.md` |
| `NOTES_exploratory_analyses.md` | Why the exploratory analyses above are not part of the published findings |

> Note on structure: the analyses were carried out as one sequentially-executed script
> rather than as independent modules. `analysis_code.R` is therefore provided as a single
> file; module headers inside it correspond to the analysis stages (training cohort, model
> validation, spatial transcriptomics, single-cell, downstream analyses, figures/tables).

## Figure scripts

Two main figures are rebuilt by dedicated scripts, because they were revised after the
consolidated script was last assembled and their panel definitions differ from it:

- `Figure2_rebuild.R` — rebuilds **main Figure 2** (nine panels a–i). It corrects two
  problems in the version inside `analysis_code.R`: panel (d) must be the HGF-positivity
  threshold sensitivity (HGF > 0 / > 1 / > 2, enrichment ratios for Hep.4 and Periportal
  LSECs), and panel (i) — the distance-gradient analysis — is not produced by
  `analysis_code.R` at all. Panel (i) is rebuilt from its result tables
  (`Final/Figure3/Figure2i_distance_gradient_HGFpos.csv` and `..._HGFneg.csv`), so the
  entire figure is reproduced from data.
- `Figure3_redraw.R` — rebuilds **main Figure 3** (panel d = the 16-gene HGF+ HSC
  signature across cell types in GSE202379) and **Figure S27** (cell-type distribution of
  the 11-gene tissue signature).

Both figure scripts fix the random seed before plotting, because `ggrepel` places labels
with a stochastic repulsion algorithm; without a seed, repeated runs produce slightly
different output. With the seed set, running either script twice yields a byte-identical
PNG, and the output matches the corresponding figure in the manuscript.

## Data sources (all publicly available)

| Dataset | Type | Role |
|---|---|---|
| GSE135251 | bulk RNA-seq, n = 216 | development cohort |
| GSE49541  | microarray, n = 72 | primary external validation |
| GSE130970 | bulk RNA-seq, n = 78 | exploratory validation |
| GSE292268 | CosMx spatial transcriptomics | spatial niche analysis |
| GSE202379 | scRNA-seq | single-cell reference |

## How to run

1. `Rscript 00_install_packages.R`
2. Edit the paths at the top of `analysis_code.R` to point at your local copies of the GEO
   datasets (the script uses `setwd("~/Desktop/sjfx3")` throughout; adjust as needed).
3. Run `analysis_code.R` sequentially. The figure scripts require the intermediate result
   tables that this produces, so run it first.

## Scope and limitations

- The code reproduces the analyses reported in the manuscript. Two exploratory analyses
  that are **not** reported in the manuscript are also present in `analysis_code.R`
  (a cis-Mendelian randomisation analysis and a scTenifold virtual-knockout analysis).
  These are documented, with their reasons for non-reporting, in
  `NOTES_exploratory_analyses.md`. They should not be interpreted as part of the published
  findings.
- Absolute paths (`~/Desktop/sjfx3/...`) reflect the working environment of the original
  analysis and must be edited before the code will run elsewhere.
- Transcriptomic data are not redistributed here; download them from GEO using the
  accessions above.

## Licence

See `LICENSE`.

## Contact

Hongju Yang (corresponding author)
