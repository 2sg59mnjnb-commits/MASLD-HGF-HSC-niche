
# A spatially defined HGF-expressing HSC niche in MASLD fibrosis

Analysis code for the manuscript by Wu & Yang.

## Data sources (all from GEO)
- GSE135251 (bulk RNA-seq, n = 216) - training
- GSE49541  (microarray, n = 72) - external validation
- GSE130970 (bulk RNA-seq, n = 78) - exploratory validation
- GSE292268 (CosMx spatial) - spatial niche
- GSE202379 (scRNA-seq) - single-cell reference

## Directory structure
- 00_setup/ - environment and configuration
- 01_training_cohort/ - DEG, WGCNA, feature selection
- 02_model_validation/ - CV, external validation, Firth, calibration
- 03_spatial_transcriptomics/ - CosMx niche analysis
- 04_single_cell/ - scRNA-seq validation, TF activity
- 05_downstream_analyses/ - F0-F4 trends, drug repurposing, MR
- 06_figures_tables/ - final figures and tables

## How to run
1. Run 00_setup/00_install_packages.R
2. Edit paths in 00_setup/01_project_config.R
3. Run scripts in each folder in numerical order

## Contact
Hongju Yang (corresponding author)

