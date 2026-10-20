# Script -> output map

Generated from the `# Output` header block of every script.

Note: two scripts share the prefix `24` (`03_spatial/24_HGFpos_neighbourhood_signature.R`
and `04_single_cell/24_cosmx_HGFpos_vs_HGFneg_DEG.R`); they belong to different stages.

| Script | Purpose | Primary outputs |
|---|---|---|
| `00_setup/00_install_packages.R` | 00_install_packages.R -- install the R packages used in this repository |  |
| `01_discovery/01_download_GSE135251.R` | Downloads the training cohort: 216 per-sample count files from GSE135251. | data/GSE135251_counts/*.counts.txt.gz   216 files; data/GSE135251_pheno.csv                sample metadata; Run time: 10-30 min depending on connection. Re-running skips existing files. |
| `01_discovery/02_count_matrix_and_DEG.R` | Builds the training expression matrix from the 216 per-sample count files and | data/derived/GSE135251_aligned_data.rds   counts + phenotype + groups; data/derived/GSE135251_DEG_results.rds    deg, deg_sig, log_expr, groups; results/01_discovery/GSE135251_DEGs.csv |
| `01_discovery/03_wgcna.R` | Weighted gene co-expression network analysis of the training cohort. | results/01_discovery/WGCNA_module_trait.csv; results/01_discovery/WGCNA_soft_threshold.csv; data/derived/WGCNA_results.rds   (net, MEs, correlations, members) |
| `01_discovery/04_candidate_genes.R` | Builds the candidate gene set for feature selection by intersecting the |  |
| `01_discovery/05_four_algorithm_selection.R` | Derives the 11-gene tissue signature. Four feature-selection algorithms are | results/01_discovery/feature_selection_genes.csv; data/derived/ML_inputs.rds   (X, y, per-algorithm gene lists, signature) |
| `02_model/06_fixed_signature_cv.R` | 10-fold cross-validation of the fixed 11-gene signature in the training | results/02_model/fixed_signature_cv.csv; data/derived/fixed_cv.rds   (per-sample CV predictions, ROC, AUC, CI) |
| `02_model/07_fully_nested_cv.R` | Fully nested cross-validation: the entire discovery pipeline is repeated | results/02_model/nested_cv_summary.csv; results/02_model/nested_cv_stability.csv; data/derived/nested_cv.rds |
| `02_model/08_external_validation_GSE49541.R` | Primary external validation of the 11-gene signature in GSE49541 | results/02_model/external_validation_GSE49541.csv; data/derived/GSE49541_validation.rds |
| `02_model/09_calibration_DCA_Firth.R` | Secondary performance analyses for the external validation cohort: | results/02_model/calibration_GSE49541.csv; results/02_model/DCA_GSE49541.csv; results/02_model/Firth_vs_MLE.csv |
| `02_model/10_CHI3L1_comparison.R` | Head-to-head comparison of the 11-gene signature against CHI3L1 alone, and | results/02_model/CHI3L1_comparison.csv; results/02_model/CHI3L1_DeLong.csv; data/derived/CHI3L1_comparison.rds |
| `02_model/11_GSE130970_exploratory.R` | Exploratory validation in GSE130970 (n = 78, METAVIR F0-F4, bariatric-surgery | results/02_model/GSE130970_validation.csv; data/derived/GSE130970_validation.rds |
| `02_model/12_NAFLD_only_rederivation.R` | Sensitivity analysis: rebuild the whole discovery pipeline after excluding the | results/02_model/NAFLD_only_rederivation.csv; results/02_model/NAFLD_only_module_correlation.csv; data/derived/NAFLD_only.rds |
| `03_spatial/13_load_CosMx_and_panel.R` | Loads both CosMx slides (GSE292268), reconstructs the patient label from the | results/03_spatial/celltype_composition.csv; results/03_spatial/panel_coverage.csv; data/derived/cosmx_flat.rds |
| `03_spatial/14_stellate_fibrotic_enrichment.R` | Tests whether stellate cells are enriched in fibrotic regions. | results/03_spatial/stellate_enrichment.csv; results/03_spatial/stellate_enrichment_sensitivity.csv |
| `03_spatial/15_HGFpos_niche_scan.R` | For every annotated cell type, counts how many neighbours each HSC has within | results/03_spatial/niche_scan_by_slide.csv; results/03_spatial/niche_scan_cross_slide.csv |
| `03_spatial/16_HGFpos_proximity_mixed_models.R` | The central spatial result: are HGF+ HSCs closer to Hep.4 hepatocytes than | results/03_spatial/hgfpos_proximity_models.csv; results/03_spatial/per_patient_delta.csv; data/derived/proximity_models.rds |
| `03_spatial/17_FOV_consistency_and_thresholds.R` | Two robustness checks for the HGF+ HSC proximity result: | results/03_spatial/fov_direction_consistency.csv; results/03_spatial/hgf_threshold_effects.csv |
| `03_spatial/18_zonation_portal_vs_central.R` | Places the niche anatomically: are Hep.4 hepatocytes nearer central-venous or | results/03_spatial/zonation_hepg4.csv; results/03_spatial/zonation_poisson.csv; results/03_spatial/zonation_permutation.csv |
| `03_spatial/19_Hep4_IL6R_and_pathways.R` | Hep.4 hepatocytes stratified by IL6R expression: | results/03_spatial/hep4_il6r_stratification.csv; results/03_spatial/hep4_il6r_distance_models.csv; results/03_spatial/hep4_pathways.csv |
| `03_spatial/20_depth_correction.R` | HGF+ HSCs have a substantially higher transcript count than HGF- HSCs | results/03_spatial/depth_comparison.csv; results/03_spatial/depth_corrected_distance.csv |
| `03_spatial/21_distance_gradient.R` | Distance-gradient analysis: how does Hep.4 enrichment around HGF+ HSCs change | results/03_spatial/distance_gradient.csv |
| `03_spatial/22_selectivity_multiple_testing.R` | Is the HGF+ HSC niche specific to Hep.4 hepatocytes, or are HGF+ HSCs simply | results/03_spatial/selectivity_paired.csv |
| `03_spatial/23_CellChat_spatial.R` | Spatially constrained ligand-receptor inference with CellChat v2. | results/03_spatial/CellChat_all_LR.csv; results/03_spatial/CellChat_HGFpos_to_IL6Rpos.csv; results/03_spatial/CellChat_HGFpos_pathways.csv |
| `03_spatial/24_HGFpos_neighbourhood_signature.R` | Negative control for the 11-gene tissue signature: is the signature score | results/03_spatial/neighbourhood_signature.csv |
| `04_single_cell/24_cosmx_HGFpos_vs_HGFneg_DEG.R` | Differential expression between HGF+ and HGF- HSCs in each CosMx slide. | results/04_single_cell/cosmx_HGFpos_vs_HGFneg_DEG.csv; results/04_single_cell/cosmx_DEG_summary.csv |
| `04_single_cell/25_pseudobulk_16gene_signature.R` | Derives the 16-gene HGF+ HSC signature from paired FOV pseudobulk. | results/04_single_cell/pseudobulk_slide1.csv, pseudobulk_slide2.csv; results/04_single_cell/signature_16genes.csv; results/04_single_cell/pseudobulk_ncount_covariate.csv |
| `04_single_cell/26_rank_matched_validation.R` | Rank-based depth matching for the 16-gene signature. | results/04_single_cell/rank_matched_validation.csv |
| `04_single_cell/27_Hep4_identity_scRNA.R` | Transcriptional identity of the Hep.4 hepatocyte subpopulation. | results/04_single_cell/Hep4_markers_cosmx.csv; results/04_single_cell/Hep4_scRNA_cluster_mapping.csv; data/derived/hep_sc_15k.rds |
| `04_single_cell/28_signature_specificity_scRNA.R` | Cell-type specificity of the two signatures in an independent scRNA-seq | results/04_single_cell/signature_score_by_celltype.csv; results/04_single_cell/signature_specificity_index.csv |
| `04_single_cell/29_TF_activity.R` | Transcription-factor activity in HGF+ versus HGF- HSCs, inferred from | results/04_single_cell/TF_activity_differential.csv; results/04_single_cell/TF_activity_consensus.csv; data/derived/tf_activity.rds |
| `04_single_cell/30_TF_HGF_coupling_and_motif.R` | Two orthogonal checks on the TFAP2A / TP53 -> HGF regulatory model: | results/04_single_cell/TF_HGF_coupling.csv; results/04_single_cell/motif_enrichment.csv |
| `05_downstream/31_stage_trends.R` | Do the signature and the IL6-family axis track fibrosis stage across the full | results/05_downstream/stage_trends_signature.csv; results/05_downstream/stage_trends_genes.csv; results/05_downstream/stage_trends_pathways.csv |
| `05_downstream/32_GSVA.R` | Pathway-level activity across fibrosis stages using GSVA with curated gene | results/05_downstream/GSVA_pathway_trends.csv; data/derived/gsva.rds |
| `05_downstream/33_deconvolution_NNLS.R` | Estimates the hepatic stellate cell fraction in each bulk sample by non-negative | results/05_downstream/deconvolution_HSC_proportion.csv; data/derived/deconvolution.rds |
| `05_downstream/34_multi_cohort_trends.R` | Re-derives the signature score in two additional cohorts that have full | results/05_downstream/multi_cohort_trends.csv; results/05_downstream/multi_cohort_gene_trends.csv |
| `05_downstream/35_incremental_value.R` | Does the signature add anything to routine clinical variables? This can only | results/05_downstream/incremental_value_GSE89632.csv; data/derived/incremental_value.rds |
| `05_downstream/36_drug_repurposing.R` | Exploratory drug-gene interaction query against DGIdb for the 16-gene HGF+ HSC | results/05_downstream/DGIdb_interactions.csv; results/05_downstream/DGIdb_drug_summary.csv |
| `05_downstream/37_MR_exploratory.R` | Exploratory cis-Mendelian randomisation, included for transparency. | results/05_downstream/MR_CHI3L1.csv |
| `06_figures/38_Figure1.R` | Figure 1. Spatial landscape of stellate cells in MASLD liver. |  |
| `06_figures/39_Figure2.R` | Figure 2. HGF+ HSC spatial niche and cross-slide consistency (panels a-h). |  |
| `06_figures/40_Figure3.R` | Figure 3. Development and validation of the HGF+ HSC signature. |  |
| `06_figures/41_Figure4_Figure5.R` | Figure 4. IL6-family receptor positive rates in HGF+ versus HGF- HSCs, |  |
| `06_figures/42_supplementary_figures.R` | Builds the supplementary figures that are generated from tabular outputs of |  |
| `07_tables/43_supplementary_tables.R` | Assembles the supplementary tables as CSV files from the analysis outputs. | results/tables/TableS*.csv; results/tables/_INDEX.csv |
| `R/config.R` | config.R -- central configuration for the MASLD HGF+ HSC niche study |  |
| `R/helpers.R` | helpers.R -- shared functions |  |
