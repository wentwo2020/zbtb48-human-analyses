"""
Shared helpers for the ZBTB48 human analyses.

Everything in this module is either an API accessor or a classification rule.
The classification rules are the load-bearing part: the endpoint domain
definitions used for the phenome-wide enrichment analysis are regular
expressions, and small changes to them move the denominators. They are
reproduced here verbatim from the analysis that generated the published
figures, with the subtleties documented, so that a re-run reproduces the
published counts rather than something close to them.
"""

import re
import time

import numpy as np
import pandas as pd
import requests

# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

GTEX_API = "https://gtexportal.org/api/v2"
GTEX_DATASET = "gtex_v10"

# GTEx v10 uses GENCODE v39. Gene identifiers are release-specific: querying
# v10 endpoints with a v26 identifier returns an empty result rather than an
# error, so the release must be pinned when resolving symbols to identifiers.
GTEX_GENCODE_VERSION = "v39"
GTEX_GENOME_BUILD = "GRCh38/hg38"

ZBTB48_GENCODE_ID = "ENSG00000204859.13"

FINNGEN_API = "https://r12.finngen.fi/api"
ENSEMBL_API = "https://rest.ensembl.org"
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

# The instrument. G is the alternate allele and raises ZBTB48 expression.
INSTRUMENT_RSID = "rs17029626"
INSTRUMENT_VARIANT_B38 = "chr1_6580056_A_G_b38"
INSTRUMENT_CHROM = "1"
INSTRUMENT_POS_B38 = 6580056

# Whole-blood cis-eQTL normalised effect size for the instrument, from GTEx
# v10. Wald-ratio odds ratios are expressed per standard deviation of
# genetically predicted expression, which means dividing the per-allele log
# odds by this value. Reported in Extended Data Fig. 7D.
NES_WHOLE_BLOOD = 0.2014

# Instrument strength, recovered as the squared z statistic of the whole-blood
# eQTL (individual-level data not being available). For a single instrument
# this gives the relative standard error of the exposure directly, which the
# delta-method Wald-ratio standard error requires. Script 03 recomputes F from
# the live eQTL p-value and will fail if it has drifted from this value.
F_INSTRUMENT = 87.1
REL_SE_EXPOSURE = 1.0 / np.sqrt(F_INSTRUMENT)

# GTEx releases donor age as ten-year brackets, not as a continuous value.
# Regressions use the bracket midpoint. This is a coarsening of the exposure
# and it biases the slope towards zero, so the reported age effects are
# conservative.
AGE_BRACKET_MIDPOINT = {
    "20-29": 24.5, "30-39": 34.5, "40-49": 44.5,
    "50-59": 54.5, "60-69": 64.5, "70-79": 74.5,
}

# ----------------------------------------------------------------------------
# Endpoint classification
# ----------------------------------------------------------------------------
#
# Three distinct regular expressions are used, and they are NOT
# interchangeable. Getting them confused changes the published denominators.
#
#   PAIN_CODES_RE   assigns an endpoint to the coarse "pain" group, which is
#                   what determines the non-pain nervous-system denominator.
#                   It includes `algia$` and `algia,`, and this matters: those
#                   alternatives match "Causalgia" (complex regional pain
#                   syndrome type II), so causalgia is classified as pain and
#                   is therefore EXCLUDED from the non-pain nervous-system
#                   domain. This is why that domain contains 131 endpoints and
#                   not 132. A rule based on the ICD-10 chapter alone would
#                   wrongly count causalgia as non-pain.
#
#   SYMPTOM_RE      defines the narrower "pain symptom" domain used for the
#                   enrichment test. It deliberately omits `algia$`, and
#                   structural spine disease is subtracted from it, giving 28
#                   endpoints.
#
#   STRUCTURAL_RE   defines degenerative spine and joint disease, which is
#                   held out as a power-matched negative control: these
#                   endpoints are larger than the pain-symptom endpoints, so a
#                   null result here cannot be explained by low power.

PAIN_CODES_RE = (
    r"^Pain\b|\bpain\b|algia$|algia,|neuralgia|Migraine|Cluster headache|"
    r"Fibromyalgia|Sciatica|Dorsalgia"
)

SYMPTOM_RE = (
    r"^Pain\b|\bpain\b|neuralgia|Fibromyalgia|Migraine|Cluster headache|"
    r"Sciatica|Dorsalgia"
)

STRUCTURAL_RE = (
    r"Spondylosis|Spondylopath|disc (prolapse|disorder)|prolapse|"
    r"Spinal stenosis|Osteoarthr|arthrosis"
)

# The skeletal sensitivity domain. Deliberately broad: it is a control asking
# whether the instrument looks like a skeletal locus, so a wide net gives the
# skeletal domain more opportunity to appear enriched. It matches 35 endpoints,
# of which 12 are traumatic fractures, 10 are bone neoplasms, 6 are other bone
# or cartilage disorders, 5 are osteoporosis endpoints and 2 are congenital
# malformations. All three associated endpoints are osteoporosis endpoints.
# Note that this domain OVERLAPS the musculoskeletal non-pain domain: 8 of the
# 35 carry an ICD-10 chapter XIII code and are counted in both. Each domain is
# tested against its own matched permutation, so the overlap does not affect
# any p-value, but the domain sizes do not sum to the number of endpoints.
BONE_RE = r"Osteoporo|fracture|bone"


def classify_endpoints(d):
    """Add the domain classification columns to a FinnGen results frame.

    Returns the frame with `grp`, `pain_symptom`, `spine_struct` and `bone`
    added. `grp` is the coarse three-level grouping; the others are the
    domain masks used by the enrichment analysis.
    """
    d = d.copy()
    ps = d.phenostring.fillna("")

    pain_codes = set(d.loc[ps.str.contains(PAIN_CODES_RE, case=False, regex=True), "phenocode"])
    neuro_nonpain = set(
        d.loc[
            d.category.fillna("").str.contains("nervous system")
            & ~d.phenocode.isin(pain_codes),
            "phenocode",
        ]
    )
    d["grp"] = np.where(
        d.phenocode.isin(pain_codes), "pain",
        np.where(d.phenocode.isin(neuro_nonpain), "neuro, non-pain", "non-neuronal / other"),
    )

    d["pain_symptom"] = ps.str.contains(SYMPTOM_RE, case=False, regex=True) & ~ps.str.contains(
        STRUCTURAL_RE, case=False, regex=True
    )
    d["spine_struct"] = ps.str.contains(STRUCTURAL_RE, case=False, regex=True)
    d["bone"] = ps.str.contains(BONE_RE, case=False, regex=True)
    d["msk_nonpain"] = (
        d.category.fillna("").str.contains("musculoskeletal", case=False) & ~d.pain_symptom
    )
    return d


# ----------------------------------------------------------------------------
# Power handling
# ----------------------------------------------------------------------------


def effective_n(n_case, n_control):
    """Case-control effective sample size, 4 / (1/n_case + 1/n_control)."""
    return 4.0 / (1.0 / n_case + 1.0 / n_control)


def quintile_bins(log_neff):
    """Assign endpoints to five strata of effective sample size.

    The enrichment test must be matched on power, because the FinnGen pain
    endpoints are far larger than the average endpoint and an unmatched
    comparison recovers sample size rather than biology.
    """
    bins = np.quantile(log_neff, [0, 0.2, 0.4, 0.6, 0.8, 1.0])
    bins[-1] += 1e-9
    return np.digitize(log_neff, bins[1:-1])


def matched_permutation(dm, mask, nperm=20000, seed=1):
    """Effective-N-matched permutation test for enrichment of association.

    FinnGen releases an effect estimate only for endpoints reaching P < 0.05,
    so a background distribution of effect sizes cannot be built. What IS
    known is the full denominator of endpoints tested per domain, which makes
    enrichment of association the testable quantity rather than effect
    direction or magnitude.

    For each stratum of effective sample size, the same number of endpoints is
    drawn with replacement from non-pain, non-spine endpoints outside the
    target domain, and the number reaching P < 0.05 is counted. Excluding
    spine endpoints from the sampling pool keeps the negative control out of
    its own null.

    Returns (observed, expected, one-sided p). The smallest attainable p is
    1 / (nperm + 1), which at 20,000 draws is 5e-5; a result at that value is
    reported as the permutation floor rather than as an exact p-value.
    """
    rng = np.random.default_rng(seed)
    mask = np.asarray(mask)
    obs = int(dm.loc[mask, "sig05i"].sum())
    pool = dm[~mask & (dm.grp != "pain") & ~dm.spine_struct]
    need = dm.loc[mask, "nb"].value_counts()

    null = np.empty(nperm, dtype=int)
    for i in range(nperm):
        total = 0
        for b, k in need.items():
            cand = pool.loc[pool.nb == b, "sig05i"].values
            if len(cand):
                total += rng.choice(cand, size=k, replace=True).sum()
        null[i] = total
    return obs, float(null.mean()), (1 + int((null >= obs).sum())) / (1 + nperm)


# ----------------------------------------------------------------------------
# API accessors
# ----------------------------------------------------------------------------


def gtex_get(endpoint, params, timeout=300, retries=3):
    """GET one GTEx v2 endpoint, returning the `data` list."""
    params = dict(params)
    params.setdefault("datasetId", GTEX_DATASET)
    for attempt in range(retries):
        try:
            r = requests.get(f"{GTEX_API}/{endpoint}", params=params, timeout=timeout)
            r.raise_for_status()
            return r.json().get("data", [])
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)


def gtex_get_paged(endpoint, params, items_per_page=10000, timeout=300):
    """GET a GTEx v2 endpoint across all pages.

    The association endpoints paginate; a single-page request silently
    truncates, which would understate the number of tissues in which a
    variant is detected.
    """
    rows, page = [], 0
    while True:
        p = dict(params)
        p.setdefault("datasetId", GTEX_DATASET)
        p.update(itemsPerPage=items_per_page, page=page)
        r = requests.get(f"{GTEX_API}/{endpoint}", params=p, timeout=timeout)
        r.raise_for_status()
        j = r.json()
        data = j.get("data", [])
        rows += data
        n_pages = j.get("paging_info", {}).get("numberOfPages", 1)
        if not data or page >= n_pages - 1:
            break
        page += 1
    return rows


def resolve_gencode_ids(symbols):
    """Map gene symbols to GENCODE v39 identifiers via the GTEx gene endpoint."""
    data = gtex_get(
        "reference/gene",
        {
            "geneId": list(symbols),
            "gencodeVersion": GTEX_GENCODE_VERSION,
            "genomeBuild": GTEX_GENOME_BUILD,
            "itemsPerPage": 300,
        },
        timeout=90,
    )
    return {g["geneSymbol"]: g["gencodeId"] for g in data if g["geneSymbol"] in set(symbols)}


def variant_position(rsid=INSTRUMENT_RSID, verify=False):
    """GRCh38 chromosome and position for the instrument.

    The position is a constant on a fixed assembly, so it is pinned rather
    than looked up: making every run depend on a live Ensembl response adds a
    failure mode for no benefit, and Ensembl does intermittently return 500.

    Pass verify=True to cross-check the pinned value against Ensembl. A
    mismatch raises; an unreachable service warns and falls through to the
    pinned value, because an outage is not evidence that the coordinate
    changed.
    """
    if rsid != INSTRUMENT_RSID:
        raise ValueError(
            f"only {INSTRUMENT_RSID} is pinned; look up {rsid} explicitly and add it"
        )
    chrom, pos = INSTRUMENT_CHROM, INSTRUMENT_POS_B38

    if verify:
        try:
            r = requests.get(
                f"{ENSEMBL_API}/variation/human/{rsid}",
                params={"content-type": "application/json"},
                timeout=60,
            )
            r.raise_for_status()
            mappings = [m for m in r.json()["mappings"] if m["assembly_name"] == "GRCh38"]
            got = (mappings[0]["seq_region_name"], mappings[0]["start"])
            if got != (chrom, pos):
                raise RuntimeError(
                    f"Ensembl places {rsid} at {got}, not the pinned {(chrom, pos)}"
                )
        except requests.RequestException as e:
            print(f"  note: could not reach Ensembl to verify the position ({e}); "
                  f"using the pinned {chrom}:{pos}")
    return chrom, pos


def fetch_finngen_variant(rsid=INSTRUMENT_RSID):
    """Fetch the full FinnGen R12 phenome-wide scan for one variant.

    Returns a DataFrame of all endpoints tested. Effect estimates are present
    only for endpoints with P < 0.05; the rest carry a null beta. This is a
    property of the release, not of the query, and it is the reason the
    analysis tests enrichment of association rather than effect direction.
    """
    chrom, pos = variant_position(rsid)
    url = f"{FINNGEN_API}/variant/{chrom}-{pos}-A-G"
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    d = pd.DataFrame(r.json()["results"])
    for c in ("n_case", "n_control", "beta", "sebeta", "pval", "mlogp"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    return d


def wald_ratio_per_sd(beta, sebeta, nes=NES_WHOLE_BLOOD, rel_se_exposure=None):
    """Single-instrument Wald-ratio MR, with the delta-method standard error.

    The point estimate is the outcome log odds divided by the eQTL effect
    size, which rescales the exposure to one standard deviation of genetically
    predicted expression.

    The standard error is the part that is easy to get wrong. Dividing the
    outcome beta and its standard error by the same constant leaves the z
    statistic unchanged, so a Wald ratio built that way inherits the raw GWAS
    p-value and implicitly treats the eQTL estimate as known without error.
    That is anti-conservative. The delta method propagates uncertainty in the
    denominator as well:

        SE(theta) / theta = sqrt( (SE_Y / beta_Y)^2 + (SE_X / beta_X)^2 )

    and the relative standard error of the exposure follows from the
    instrument strength, since for a single instrument F is the squared z
    statistic of the eQTL: SE_X / beta_X = 1 / sqrt(F).

    The difference is not cosmetic here. Under the raw GWAS p-value the
    strongest endpoint in the pain family reaches 1.6e-5, which would be
    phenome-wide significant against a 2.0e-5 threshold. Under the delta
    method it is 9.0e-5, which is nominal. The paper's statement that
    individual associations are nominal rather than phenome-wide significant
    is correct under the correct standard error and would be wrong under the
    naive one.

    Returns (OR, lower, upper, se_log_or, p_two_sided).

    Note separately that the instrument is a considerably stronger splicing
    QTL than expression QTL, so the instrumented exposure may incorporate
    isoform usage as well as total transcript abundance. That is a limitation
    of the exposure definition, not of this calculation.
    """
    if rel_se_exposure is None:
        rel_se_exposure = REL_SE_EXPOSURE
    beta = np.asarray(beta, dtype=float)
    sebeta = np.asarray(sebeta, dtype=float)
    log_or = beta / nes
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.sqrt((sebeta / np.abs(beta)) ** 2 + rel_se_exposure ** 2)
        z = 1.0 / rel
    se_log_or = np.abs(log_or) * rel
    p = 2.0 * _norm_sf(z)
    return (
        np.exp(log_or),
        np.exp(log_or - 1.96 * se_log_or),
        np.exp(log_or + 1.96 * se_log_or),
        se_log_or,
        p,
    )


def _norm_sf(z):
    from scipy import stats as _st

    return _st.norm.sf(z)
