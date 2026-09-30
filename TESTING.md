# Testing status

Each script asserts its own key results against the values published in the
figures. This records which scripts were executed end to end before deposit,
and which were not, so the reader does not have to infer it.

| Script | Executed end to end | Assertions | Notes |
|---|---|---|---|
| `01_gtex_age_association.py` | yes | pass | Reproduces tibial nerve n = 670, +2.77 per cent per decade, P = 3.39e-8, and finds exactly two nervous-system tissues at q < 0.05 in opposite directions. |
| `02_gtex_reference_genes_and_composition.py` | yes | pass | Reproduces 23 genes at n = 670, 17 at q < 0.05, and the three composition contributions summing to -0.12 per cent per decade. |
| `03_gtex_qtl_lead_variant.py` | yes | pass | Reproduces 23 eQTL tissues with the instrument leading in 12, 35 sQTL tissues with it leading in 34, F = 87.1, and the two opposite-direction blood sQTL junctions. |
| `04_snrna_cellxgene_spinal.py` | **no** | not run | Requires `cellxgene-census` and `tiledbsoma`, which were not installed in the environment used for this deposit. The downstream half of this script, which renders the panel and recomputes every statistic from the deposited donor matrix, WAS executed in a clean directory and reproduced all reported values. The Census query half is syntax-checked only. |
| `05_snrna_geo_spinal.py` | yes | pass | Reproduces the 11 ALS and 9 control donor split and all three per-ACTB ratios, with the direction flip against raw counts per million intact. Observed values below. |
| `06_finngen_mr_and_enrichment.py` | yes | pass | Reproduces 11 of 12 risk-increasing, the rheumatoid arthritis exception at OR 0.699 and P = 0.0053, all five domain denominators, and the 3.81-fold pain enrichment at the permutation floor. |

## Observed output of script 05

Copied from `results/fig5d_reconciliation.csv` as written by the run, not from
recollection. All three assertions on these values passed: every per-ACTB ratio
above 1, every raw counts-per-million ratio below 1, and each per-ACTB ratio
within 15 per cent of the plotted value.

| Dataset | Plotted | Per 1,000 ACTB | Raw CPM | ACTB microglia/astrocyte | n paired |
|---|---|---|---|---|---|
| GSE190442, lumbar, 7 donors | 1.24 | 1.404 | 0.595 | 2.358 | 7 |
| GSE330130 all, 11 ALS + 9 control | 1.30 | 1.351 | 0.607 | 2.227 | 20 |
| GSE330130 control only, 9 donors | 1.68 | 1.756 | 0.906 | 1.937 | 9 |

The paired Wilcoxon reaches significance in none of these three, which is why
the panel must be reported as a consistent direction across datasets anchored
by the one individually significant comparison (the CZ CELLxGENE dataset,
P = 0.003) rather than as independent replications.

## One correction the assertions caught

Writing this repository exposed a real error in an intermediate
reimplementation, which is worth recording because it changes a claim in the
paper.

A single-instrument Wald ratio can be built two ways. Dividing the outcome beta
and its standard error by the same eQTL effect size leaves the z statistic
unchanged, so the resulting p-value is just the raw GWAS p-value and the eQTL
estimate is implicitly treated as known without error. The published analysis
instead uses the delta-method standard error, which propagates uncertainty in
the denominator:

    SE(theta) / theta = sqrt( (SE_Y / beta_Y)^2 + (SE_X / beta_X)^2 )

with the exposure's relative standard error given by the instrument strength,
`1 / sqrt(F)` with F = 87.1.

The first version of `06_finngen_mr_and_enrichment.py` used the naive form. It
reproduced every odds ratio but none of the p-values, and its assertion that
the family p-values are nominal rather than phenome-wide significant failed.
That failure was correct and informative: under the naive standard error the
strongest endpoint reaches 1.6e-5 against a phenome-wide threshold of 2.0e-5
and WOULD be phenome-wide significant, while under the delta method it is
9.0e-5 and is not. The paper's statement is right, and it is right only
because the published analysis used the more conservative standard error.

The relationship was confirmed rather than assumed: the implied exposure
relative standard error, recovered independently from each of the twelve
endpoints, is constant at 0.10716 with a coefficient of variation of 0.09 per
cent, and 1 / 0.10716^2 = 87.1, which is exactly the reported instrument F.
The corrected implementation reproduces the published odds ratios, standard
errors and p-values for all twelve endpoints.
