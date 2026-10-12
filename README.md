# Spatially resolved transcriptomics identifies an HGF+ HSC niche in MASLD fibrosis

Analysis code for the manuscript.

## Data sources

- GSE135251 (bulk RNA-seq, training)
- GSE49541 (Affymetrix, external validation)
- GSE130970 (bulk RNA-seq, exploratory validation)
- GSE292268 (CosMx spatial transcriptomics)
- GSE202379 (scRNA-seq reference)

## Pipeline

- 00_setup.R — environment setup
- 01_CosMx_preprocess.R — cell type composition
- 02_CosMx_niche.R — HGF+ HSC niche scan
- 03_CosMx_TF_drug.R — TF activity and drug prediction
- 04_scRNA_HGF_HSC.R — scRNA validation
- 05_GSE135251_WGCNA.R — DEG and WGCNA
- 06_feature_selection_CV.R — four-algorithm selection
- 07_external_validation.R — GSE49541 validation
- 08_downstream_pathways.R — JAK/STAT and NF-kB trends
- 09_figures.R — all figures

## Requirements

R >= 4.3.

## Contact

Corresponding author.