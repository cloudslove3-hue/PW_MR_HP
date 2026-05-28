# PW_MR_HP
Phenome-wide Mendelian randomization of H. pylori antibody-response traits against 17 chronic diseases (Bang, 2026)
# Phenome-wide Mendelian Randomization of *Helicobacter pylori* Antibody-Response Traits

Analysis code for Bang CS. "Phenome-wide Mendelian randomization of seven
*Helicobacter pylori* antibody-response traits: calibration against
multiple-testing correction, HLA-region sensitivity, and eradication-trial
triangulation."

## Requirements
- R 4.6.0
- TwoSampleMR (v0.6+), ieugwasr, MendelianRandomization, MRPRESSO, coloc, MVMR
  (all from the MRC-IEU r-universe)
- An OpenGWAS API token (set `OPENGWAS_JWT` in `.Renviron`; see `.Renviron.example`)

## Reproducibility
All stochastic steps use `set.seed(20260521)`.

## Pipeline order
| Step | Script | Purpose |
|------|--------|---------|
| 1 | `scripts/01_exposure.R` | Extract 7-antigen instruments from Butler-Laporte 2020 |
| … | … | … |
| 12 | `scripts/12_gastric_positive_control.R` | Gastric cancer positive control |
| 14 | `scripts/14_reconciliation.R` | Prior-MR-claim reconciliation (Table 7) |
| 16 | `scripts/16_hla_exclusion.R` | HLA-region exclusion sensitivity (Table 8) |
| 18 | `scripts/18_fig7_triangulation.R` | RCT triangulation figure |

## Outputs
See `results/` for all derived tables (Supplementary Tables S1–S6).

## Data availability
Summary statistics are not redistributed here. Obtain from:
- IEU OpenGWAS (https://api.opengwas.io/)
- EBI GWAS Catalog: GCST90006910–GCST90006916 (exposures), GCST90018849 (gastric cancer)
- FinnGen Release 5 (https://www.finngen.fi/en/access_results)

## Citation
[manuscript citation once published]

## License
[MIT or CC-BY-4.0]
