# ZBTB48 human analyses

Analysis code for the human data in *ZBTB48 maintains neuropathic pain in the aged nervous system*, Muralidharan et al.

Every human analysis in the paper is a reanalysis of public or published data.
There is no primary human data collection anywhere in the study, so each script
here starts from a public accession or API and ends with the numbers that
appear in the figures. Each script also asserts its own key results against the
published values, so a re-run that silently drifts will fail rather than
produce something plausible.

## What is and is not in this repository

**In scope.** The ten public human data sources behind Figure 5 and Extended
Data Figure 7:

| Script | Analysis | Panels | Source |
|---|---|---|---|
| `src/01_gtex_age_association.py` | Age association of *ZBTB48* across GTEx tissues | 5A, 5B | GTEx v10 |
| `src/02_gtex_reference_genes_and_composition.py` | Reference-gene panel and cell-composition sensitivity | ED7 A, B | GTEx v10 |
| `src/03_gtex_qtl_lead_variant.py` | eQTL and sQTL lead-variant status, instrument strength | ED7 D | GTEx v10 |
| `src/04_snrna_cellxgene_spinal.py` | Spinal cord single-nucleus, donor-paired | 5C, ED7 C | CZ CELLxGENE Census |
| `src/05_snrna_geo_spinal.py` | Two spinal single-nucleus replication datasets | 5D | GSE190442, GSE330130 |
| `src/06_finngen_mr_and_enrichment.py` | Mendelian randomisation and phenome-wide pain specificity | 5G, 5H, ED7 E, G, H | FinnGen R12 |
| `figures/` | Panel rendering from the deposited source tables | 5C, 5H | (offline) |

**Not in scope.** The generalised estimating equation models of the mouse
behavioural data reported in the paper were fitted separately and are not part
of this repository. The statistical routines used by the code in this
repository are listed under "Statistics" below.

Also not here: the mouse experiments, the histology and imaging quantification,
the proximity-labelling proteomics, and the osteoporosis replication
meta-analysis of Extended Data Fig. 7F, which uses summary statistics obtained
directly from the GWAS Catalog and BioBank Japan rather than through code in
this repository.

## Reproducing

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python src/01_gtex_age_association.py
python src/02_gtex_reference_genes_and_composition.py
python src/03_gtex_qtl_lead_variant.py
python src/06_finngen_mr_and_enrichment.py
```

Scripts 04 and 05 have extra requirements. Script 04 needs
`cellxgene-census` and `tiledbsoma`, which are heavier than the rest and are
listed separately in `requirements-singlecell.txt`. Script 05 downloads about
0.65 GB from the GEO FTP site into `data/`, cached between runs.

Outputs land in `results/`. Both directories are gitignored apart from a
placeholder.

## Pinned releases

Release identifiers are pinned in `src/common.py` and must stay pinned. Three
of these will drift and silently change the numbers:

| Resource | Pinned value | Why it matters |
|---|---|---|
| GTEx | `gtex_v10` | Gene identifiers are GENCODE-release-specific. Querying v10 with a v26 identifier returns an empty result rather than an error. |
| GTEx GENCODE | `v39` | As above, when resolving gene symbols to identifiers. |
| CZ CELLxGENE Census | `2025-11-08` | `"stable"` advances. Using `"stable"` will eventually return a different set of datasets. |
| FinnGen | R12 | Endpoint definitions and case counts change between releases. |

## Statistics

The routines used, and where:

- Ordinary least squares regression of log2 expression on donor age-bracket
  midpoint (`scipy.stats.linregress`), with Spearman correlation reported
  alongside; scripts 01 and 02.
- Benjamini-Hochberg false discovery rate across all tissues or all panel
  genes (`statsmodels.stats.multitest.multipletests`); scripts 01, 02.
- Wilcoxon signed-rank tests paired within donor, Holm-corrected across the
  family of comparisons; scripts 04 and 05.
- Friedman test as an omnibus across cell populations, on complete cases;
  script 04.
- Exact binomial test of per-donor ranking against a 1/k null, where k is the
  number of populations compared; script 04.
- Single-instrument Wald-ratio Mendelian randomisation, with within-family
  Bonferroni and Benjamini-Hochberg; script 06.
- Effective-sample-size-matched permutation test of enrichment of association,
  20,000 draws across five strata; script 06.
- Recovery of instrument strength as the squared z statistic from the eQTL
  p-value, individual-level data not being available; script 03.

**No generalized estimating equation is fitted anywhere in this repository.**

## Three things that will silently produce a wrong answer

These are the traps that actually bit during the analysis, and they are
documented in the code as well as here.

1. **FinnGen suppresses effect estimates above P = 0.05.** Returned p-values
   cap just under 0.05 and only 292 of 2,470 endpoints carry an estimate. A
   background distribution of effect sizes therefore cannot be built from this
   release, which is why pain specificity is tested as enrichment of
   association rather than by effect direction or magnitude. Treating the
   released subset as if it were the whole phenome inflates every effect-based
   comparison.

2. **A variant being *detected* in a tissue is not the same as it being the
   *lead* variant there.** rs17029626 is a significant *ZBTB48* cis-eQTL in 16
   of 23 tissues but the lead in 12, and a significant cis-sQTL in all 35
   tissues but the lead in 34. Script 03 computes both and asserts both.

3. **`ACTB` is about twice as abundant in microglia as in astrocytes in spinal
   single-nucleus data.** Per-library and per-`ACTB` normalisation therefore
   disagree in *direction*, not merely in magnitude. Never carry an `ACTB`
   ratio from one dataset to another platform or region; script 05 measures it
   per dataset and exports it alongside every ratio.

A fourth, narrower trap: the endpoint-domain classification in script 06 uses
a pattern whose `algia$` alternative matches "Causalgia". Because causalgia is
complex regional pain syndrome type II, this correctly classifies it as pain
and excludes it from the non-pain nervous-system domain, which is why that
domain contains 131 endpoints and not 132. A classification based on the
ICD-10 chapter alone gets this wrong.

## Data availability

| Source | Accession or release | Route |
|---|---|---|
| GTEx | v10 | Portal API v2, `https://gtexportal.org/api/v2/` |
| CZ CELLxGENE Discover Census | 2025-11-08 | `cellxgene-census` Python API |
| Human spinal cord snRNA | GSE190442 | GEO FTP supplementary files |
| Human spinal cord snRNA | GSE330130 | GEO FTP supplementary files |
| FinnGen | R12 | `https://r12.finngen.fi/api/variant/` |
| Ensembl | current | REST API, variant coordinates only |

Analyses not covered by this repository use published human dorsal root
ganglion data (Ray et al. 2023; GSE77968 and GSE78150), the longitudinal whole
blood cohort GSE177034, GWAS Catalog studies GCST90474133 and GCST90038656,
and BioBank Japan summary statistics (PMID 34594039).

## Licence

MIT. See `LICENSE`.
