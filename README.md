# Exploratory analyses not reported in the manuscript

The scripts and notes in this folder concern analyses that were carried out during the
study but that are **not** reported in the manuscript. They are deposited for transparency,
so that the corresponding code inside `analysis_code.R` is not left as an unreproducible
fragment. **Their results should not be interpreted as part of the published findings.**

---

## 1. cis-Mendelian randomisation — `MR_data_preparation.sh`, `MR_pipeline_cis_eQTL.R`

A cis-MR analysis testing whether genetically proxied expression of the signature genes
and of the HGF-MET-IL6R axis is associated with liver fibrosis / cirrhosis.

**Status: not completed, not reported.** Three reasons, all documented in the script headers:

1. The full eQTLGen cis-eQTL resource (~3.7 GB) could not be downloaded reliably in the
   analysis environment. The partial local copy contained only **1,204 of ~16,000 genes**.
2. **HGF, MET, IL6R, LIFR and OSMR had no genome-wide significant *cis*-eQTL in the
   available blood eQTL data.** The mechanism axis central to the manuscript therefore
   could not be tested by this method at all. Only four of the eleven signature genes were
   instrumentable: CHI3L1, GSN, IL11RA and SULF2.
3. A GTEx v8 eQTL archive was obtained, but the download was truncated alphabetically at
   `Heart_Atrial_Appendage` and did not include liver - the tissue of interest.

An earlier, methodologically flawed version of this analysis is present in
`analysis_code.R` (module `05_downstream_analyses`). It is superseded by
`MR_pipeline_cis_eQTL.R`, which corrects three problems:

- it uses **measured allele frequencies** instead of an assumed `MAF = 0.25` to convert
  eQTLGen Z-scores to effect sizes;
- it uses a **liver-fibrosis-specific outcome** (FinnGen R12 `K11_FIBROCHIRLIV`) instead
  of a NAFLD GWAS;
- it adds **colocalisation** (`coloc.abf`) rather than relying on a single-SNP Wald ratio.

If this analysis is ever completed, the negative or uninformative result for the mechanism
axis must be reported as such.

---

## 2. scTenifold virtual knockout of HGF in the HSC network

Located in `analysis_code.R` (`05_downstream_analyses`). **Completed but not reported in
the manuscript.**

The in-silico knockout of *HGF* in a gene regulatory network built from the single-cell HSC
compartment (8,327 genes x 3,228 cells, 10 networks of 500 cells each) flagged **3 of 8,327
genes** as significantly differentially regulated at FDR < 0.05:

| gene  | distance | Z    | p.adj |
|-------|----------|------|-------|
| HGF   | 1.74e-04 | 5.53 | 0 |
| RBMS3 | 2.25e-05 | 4.38 | 6.73e-28 |
| PDE4D | 9.71e-06 | 3.93 | 1.34e-03 |

Not used in the manuscript: *HGF* is the knocked-out gene itself, leaving only two
downstream candidates (RBMS3, PDE4D), which was judged too thin a basis for any claim.

> Note on file naming: the outputs retained on disk carry a `_FAILED` suffix from an
> earlier, unsuccessful run of the same script. The values above are from the successful
> run recorded in `hgf_ko.log`. Any future re-run should write to unsuffixed filenames to
> avoid confusion.
