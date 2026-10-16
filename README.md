# A spatially defined HGF-expressing HSC niche in MASLD fibrosis

Analysis code for the manuscript by Wu & Yang.

## What is in this repository

- `analysis_code.R` — the **complete** analysis code for the study, as run. It is a single
  consolidated script organised into numbered modules (Module 0 onward) that follow the
  order of the Methods section.
- `00_install_packages.R` — package installation
- `01_project_config.R` — central configuration (paths, seed, parameters, colours)

> Note on structure: the analyses were carried out as one sequentially-executed script
> rather than as independent modules. `analysis_code.R` is therefore provided as a single
> file; module headers inside it correspond to the analysis stages (training cohort,
> model validation, spatial transcriptomics, single-cell, downstream analyses, figures/tables).

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
2. Edit the paths at the top of `analysis_code.R` to point at your local copies of the
   GEO datasets (the script uses `setwd("~/Desktop/sjfx3")` throughout; adjust as needed).
3. Run `analysis_code.R` sequentially.

## Scope and limitations

- The code reproduces the analyses reported in the manuscript. Some exploratory analyses
  that are **not** reported in the manuscript are also present in `analysis_code.R`
  (notably a cis-Mendelian randomisation analysis and a scTenifold virtual-knockout
  analysis); they are included for completeness and should not be interpreted as part of
  the published findings.
- Absolute paths (`~/Desktop/sjfx3/...`) reflect the working environment of the original
  analysis and must be edited before the code will run elsewhere.
- Transcriptomic data are not redistributed here; download them from GEO using the
  accessions above.

## Licence

See `LICENSE`.

## Contact

Hongju Yang (corresponding author)
