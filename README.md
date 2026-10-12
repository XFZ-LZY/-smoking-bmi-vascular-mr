# Lifetime smoking, BMI and vascular outcomes

Reproducibility materials for **Genetic liability to lifetime smoking and body mass index across coronary, peripheral and cerebral vascular disease: a Mendelian randomization study**.

## What the audit changed

The four stroke inputs were obtained under OpenGWAS IDs ebi-a-GCST005840/41/42/43. Comparison with the original MEGASTROKE files identifies them as **TRANS multiancestry results**, despite European metadata labels. The original inputs are retained as historical records. Use `audit/S14_corrected_stroke_metadata.csv` for reporting, and never use the archived stroke sample-size fields for a sample-size-dependent analysis.

Stroke case counts: AIS 60,341; LAA 6,688; SVS 11,710; CES 9,006. Nominal control counts follow the primary paper's Supplementary Table 1: 454,450 minus deCODE's 255,213 overall controls plus its phenotype-specific controls (AIS 255,213; LAA 146,392; SVS 192,662; CES 204,570). This yields 454,450; 345,629; 391,899; 403,807. These are derived study-level nominal counts, **not per-SNP effective sample sizes or a reconstruction of unique participants**. The current Catalog ancestry records sum to 2,734 more controls in each case. The discrepancy is recorded rather than silently averaged or used for weights. Participant overlap remains unquantified.

A post-hoc instrument audit queries all 188 target variants and 9 additional actually used proxies. The protocol fixes the phenotype categories and significance threshold. The original analysis is retained, and two exclusion sensitivities are reported with all six outcomes, FDR values, conditional F statistics and exact exclusion lists. Hits indicate possible shared pathways, not proven invalid instruments; absence of a curated hit does not prove independence. The screen does not expand to all LD proxies. Cardiometabolic exclusions may remove vertical pathways.

After both exclusion schemes, the most persistent adjusted associations are smoking with CAD/PAD and BMI with CAD. Other associations are instrument-sensitive. This repository does not claim that every publication-readiness criterion is satisfied.

## Run offline

Use R 4.6.1 with TwoSampleMR 0.7.9, MVMR 0.4.8 and MendelianRandomization 0.10.0; Python 3.12. Run from the repository root:

```text
python code/reproduce.py /path/to/Rscript
```

This reconstructs the primary IVW/MVMR estimates, runs the deterministic audit, and reruns the instrument-exclusion sensitivity from the archived query evidence. It verifies the complete 60-row sensitivity output against its frozen reference. Historical bootstrap/mode/PRESSO outputs are retained rather than replaced by a new random realization.

Optional: `python code/plot_independence_audit.py` (NumPy/Matplotlib required). Figure S10 uses the figures4papers scientific-figure-making conventions: Arial, blue smoking, teal BMI, 600-dpi PNG, embedded-font PDF and editable-text SVG.

Refreshing GWAS Catalog evidence is a distinct online action: `python code/repair_three_gaps.py`. Results may change with database updates. Original source verification: `python code/fingerprint_stroke_sources.py --full` (large public downloads). To reproduce the historical query exactly, use the archived JSON, query status and classifications.

## Sources and reuse

MEGASTROKE: Malik et al., Nature Genetics 2018, DOI https://doi.org/10.1038/s41588-018-0058-3 and its Supplementary Table 1. Original source files and readme: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005840/ . Other phenotype accessions and sources are in `audit/S1_corrected_data_sources.csv`.

The MEGASTROKE project received funding from sources specified at http://www.megastroke.org/acknowledgments.html . Original author and consortium credits are in the cited publication. The source readme's data-use conditions continue to apply; in particular, the underlying MEGASTROKE data are not authorized for studies of intelligence, education, social outcomes or sensitive behavioral traits. Here, stroke remains the outcome and external trait annotations are used only to audit instrument validity.

GWAS Catalog API documentation: https://www.ebi.ac.uk/gwas/docs/programmatic-access/rest-api/ . Query date: 2026-10-12. No credentials or individual-level participant records are included. Third-party summary data retain their original terms; this repository does not relicense them.

`reproducibility_bundle.zip` contains the same machine-readable inputs and analysis code for a single-download offline run. `SHA256SUMS.csv` records the unpacked files. The manuscript should cite a fixed commit rather than a moving branch.
