"""
Figure 5G, 5H and Extended Data Fig. 7E, 7G, 7H: Mendelian randomisation and
phenome-wide pain specificity.

Three analyses on one FinnGen R12 phenome-wide scan of rs17029626.

1. Single-instrument Mendelian randomisation on a family of 12 pain endpoints
   chosen a priori. Wald-ratio odds ratios per standard deviation of
   genetically predicted ZBTB48, with within-family Bonferroni and
   Benjamini-Hochberg. Genetically predicted higher ZBTB48 raises risk in 11
   of the 12; the exception is other or unspecified rheumatoid arthritis.
   Four of the 12 are overlapping back-pain or broad-pain definitions that
   share cases, so a sign test over them assumes more independence than the
   endpoints have. This is stated rather than corrected, because the
   informative result is the enrichment analysis below, not the sign count.

2. Pain specificity by ENRICHMENT rather than by direction. This choice is
   forced by the data, not preferred: FinnGen releases an effect estimate only
   for endpoints reaching P < 0.05, so no background distribution of effect
   sizes exists. What is more, 88 per cent of the released estimates at this
   variant are risk-increasing, so a direction-concordance argument would
   carry almost no information. The full denominator of endpoints tested per
   domain IS known, which makes enrichment of association testable.

3. Domain controls. Degenerative spine disease, non-pain musculoskeletal
   disease, skeletal endpoints and non-pain nervous-system disease are each
   tested against their own effective-N-matched permutation. The spine control
   is the important one: those endpoints have MORE power than the pain-symptom
   endpoints, so their null result cannot be explained by low power.

Output: results/zbtb48_finngen_r12_ALL_endpoints.csv
        results/zbtb48_mr_prespecified_family.csv
        results/zbtb48_mr_enrichment_by_domain.csv
Runtime: about one minute; the permutations dominate.
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats as st
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (  # noqa: E402
    NES_WHOLE_BLOOD,
    classify_endpoints,
    effective_n,
    fetch_finngen_variant,
    matched_permutation,
    quintile_bins,
    wald_ratio_per_sd,
)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")

# The pain family, fixed before the phenome scan was inspected. Rheumatoid
# arthritis is retained even though it is the one discordant endpoint: an
# autoimmune endpoint in the family makes the test harder, not easier, and
# dropping it after seeing the result would be selection on the outcome.
PAIN_FAMILY = [
    "PAIN", "M13_DORSALGIA", "M13_LOWBACKPAINORANDSCIATICA", "M13_LOWBACKPAIN",
    "M13_LIMBPAIN", "G6_MIGRAINE", "M13_INTERVERTEB", "M13_SPONDYLOPATHY",
    "G6_MIGRAINE_WITH_AURA", "M13_FIBROMYALGIA", "G6_CLUSTHEADACHE_WIDE",
    "RHEUMA_NOS",
]

# The four endpoints in the family that are overlapping definitions of back or
# broad pain and therefore share cases.
NESTED_IN_FAMILY = [
    "PAIN", "M13_DORSALGIA", "M13_LOWBACKPAIN", "M13_LOWBACKPAINORANDSCIATICA",
]


def main():
    d = fetch_finngen_variant()
    d = classify_endpoints(d)
    d["sig05"] = d.beta.notna()

    n_tested = len(d)
    released = d[d.sig05]
    phenome_threshold = 0.05 / n_tested

    # Effect scales. Both are exported so the reader can see which is which:
    # the per-allele odds ratio is the raw GWAS scale, the per-SD odds ratio is
    # the Wald ratio and is what the figures plot.
    d["OR_per_allele"] = np.exp(d.beta)
    (
        d["OR_per_sd"], d["lo_per_sd"], d["hi_per_sd"],
        d["se_log_or_per_sd"], d["mr_p_two_sided"],
    ) = wald_ratio_per_sd(d.beta, d.sebeta)
    d["estimate_released"] = d.beta.notna()

    # Two different p-values are in play and they must not be conflated.
    #
    #   `pval`             FinnGen's own association p-value. Used ONLY to
    #                      decide whether an endpoint counts as associated in
    #                      the enrichment test, because that threshold is a
    #                      property of the data release: FinnGen releases an
    #                      estimate exactly when pval < 0.05.
    #
    #   `mr_p_two_sided`   the delta-method Wald-ratio p-value, which
    #                      propagates uncertainty in the eQTL denominator.
    #                      Used for every Mendelian randomisation claim. It is
    #                      larger than `pval`, by a factor that depends on the
    #                      outcome's own precision.
    d["neff"] = np.where(
        (d.n_case > 0) & (d.n_control > 0), effective_n(d.n_case, d.n_control), np.nan
    )
    d["excluded_from_models"] = ~((d.n_case > 0) & (d.n_control > 0))
    d["exclusion_reason"] = np.where(
        d.excluded_from_models, "n_control = 0, effective N undefined", ""
    )
    d["group"] = np.where(
        d.pain_symptom, "pain symptom",
        np.where(
            d.spine_struct, "degenerative spine / structural",
            np.where(d.grp == "neuro, non-pain", "nervous system, non-pain", "other"),
        ),
    )

    os.makedirs(OUT, exist_ok=True)
    keep = [
        "phenocode", "phenostring", "category", "group", "n_case", "n_control", "neff",
        "beta", "sebeta", "pval", "mlogp", "OR_per_allele", "OR_per_sd", "lo_per_sd",
        "hi_per_sd", "se_log_or_per_sd", "mr_p_two_sided", "estimate_released",
        "pain_symptom", "spine_struct", "bone", "excluded_from_models", "exclusion_reason",
    ]
    # Endpoints dropped from the models are retained in the export with a
    # reason, so that nothing is silently missing from the released data.
    d[keep].sort_values(["group", "pval"], na_position="last").to_csv(
        os.path.join(OUT, "zbtb48_finngen_r12_ALL_endpoints.csv"), index=False
    )

    # --- 1. Mendelian randomisation on the pre-specified family -------------
    fam = d[d.phenocode.isin(PAIN_FAMILY)].copy()
    if len(fam) != len(PAIN_FAMILY):
        missing = set(PAIN_FAMILY) - set(fam.phenocode)
        raise RuntimeError(f"pain-family endpoints absent from the scan: {sorted(missing)}")
    fam["mr_or_per_SD"] = fam.OR_per_sd
    fam["mr_se"] = fam.se_log_or_per_sd
    k = len(fam)
    fam["bonferroni_threshold_within_family"] = 0.05 / k
    fam["q_fam_bh"] = multipletests(fam.mr_p_two_sided, method="fdr_bh")[1]
    fam["passes_within_family_bonferroni"] = fam.mr_p_two_sided < 0.05 / k
    fam["is_nested_back_or_broad_pain"] = fam.phenocode.isin(NESTED_IN_FAMILY)
    fam["family_definition"] = (
        f"{k} pain endpoints selected a priori; phenome-wide correction does not apply"
    )
    fam[
        [
            "phenocode", "phenostring", "n_case", "n_control", "beta", "sebeta",
            "mr_or_per_SD", "mr_se", "lo_per_sd", "hi_per_sd", "pval", "mr_p_two_sided", "q_fam_bh",
            "bonferroni_threshold_within_family", "passes_within_family_bonferroni",
            "is_nested_back_or_broad_pain", "family_definition",
        ]
    ].sort_values("mr_or_per_SD").to_csv(
        os.path.join(OUT, "zbtb48_mr_prespecified_family.csv"), index=False
    )

    n_up = int((fam.mr_or_per_SD > 1).sum())
    discordant = fam[fam.mr_or_per_SD < 1]

    print("=== 1. Mendelian randomisation, pre-specified pain family ===")
    print(f"endpoints                     {k}")
    print(f"OR per SD > 1                 {n_up} of {k}                 (published 11 of 12)")
    for _, r in discordant.iterrows():
        print(
            f"discordant                    {r.phenostring} "
            f"OR {r.mr_or_per_SD:.3f}, P {r.mr_p_two_sided:.4f}   (published OR 0.699, P 0.0053)"
        )
    print(f"nested back/broad-pain        {int(fam.is_nested_back_or_broad_pain.sum())} of {k}")
    print(f"passing within-family Bonf.   {int(fam.passes_within_family_bonferroni.sum())} of {k}")
    print(f"phenome-wide threshold        {phenome_threshold:.2e}")
    print(
        f"smallest family P (delta)     {fam.mr_p_two_sided.min():.1e}  "
        f"-> nominal, not phenome-wide significant"
    )
    print(
        f"smallest family P (raw GWAS)  {fam.pval.min():.1e}  "
        f"-> would be phenome-wide significant; this is why the delta-method SE matters"
    )

    # --- 2. why enrichment and not direction --------------------------------
    frac_up = (released.beta > 0).mean()
    print("\n=== 2. Why specificity is tested by enrichment, not direction ===")
    print(f"endpoints tested              {n_tested:,}")
    print(f"estimates released (P<0.05)  {len(released)}")
    print(
        f"risk-increasing              {int((released.beta > 0).sum())} of {len(released)} "
        f"({100 * frac_up:.1f} per cent)"
    )
    print(
        "direction carries almost no information at this variant, so the test is "
        "enrichment of association"
    )

    # --- 3. domain enrichment ------------------------------------------------
    dm = d[(d.n_case > 0) & (d.n_control > 0)].copy()
    dm["neff"] = effective_n(dm.n_case, dm.n_control)
    dm["log_neff"] = np.log10(dm.neff)
    dm["nb"] = quintile_bins(dm.log_neff)
    dm["sig05i"] = dm.sig05.astype(int)

    domains = [
        ("Pain symptoms", dm.pain_symptom.values),
        ("Degenerative spine", dm.spine_struct.values),
        ("Musculoskeletal, non-pain", dm.msk_nonpain.values),
        ("Bone / osteoporosis", dm.bone.values),
        ("Nervous system, non-pain", (dm.grp == "neuro, non-pain").values),
    ]
    rows = []
    for label, mask in domains:
        obs, exp, p = matched_permutation(dm, mask)
        rows.append(
            dict(
                Group=label,
                n_endpoints_tested=int(mask.sum()),
                obs=obs,
                exp=exp,
                fold=obs / exp if exp else np.nan,
                p=p,
                permutation_draws=20000,
                strata=5,
            )
        )
    enr = pd.DataFrame(rows)
    enr.to_csv(os.path.join(OUT, "zbtb48_mr_enrichment_by_domain.csv"), index=False)

    print("\n=== 3. Domain enrichment, effective-N matched ===")
    print(enr.round(4).to_string(index=False))
    print(
        "\nNOTE: the skeletal domain OVERLAPS the musculoskeletal domain "
        "(8 endpoints are in both), so the domain sizes do not sum to the number "
        "of endpoints. Each domain is tested against its own matched permutation."
    )

    # --- checks against the published values ---------------------------------
    g = enr.set_index("Group")
    assert n_up == 11, f"expected 11 of 12 risk-increasing, got {n_up}"
    assert list(discordant.phenocode) == ["RHEUMA_NOS"], list(discordant.phenocode)
    assert int(fam.is_nested_back_or_broad_pain.sum()) == 4
    assert fam.mr_p_two_sided.min() > phenome_threshold, (
        "family p-values must be nominal rather than phenome-wide significant"
    )
    assert frac_up > 0.8, frac_up
    # The delta-method p-values must reproduce the published family, which the
    # naive Wald ratio does not: RA is published at P = 0.0053, and the raw
    # GWAS p-value for that endpoint is 0.0035.
    ra = fam[fam.phenocode == "RHEUMA_NOS"].iloc[0]
    assert abs(ra.mr_p_two_sided - 0.0053) < 2e-4, (
        f"RA delta-method P is {ra.mr_p_two_sided:.4f}, published 0.0053; "
        "check F_INSTRUMENT and NES_WHOLE_BLOOD in common.py"
    )
    assert abs(ra.mr_or_per_SD - 0.699) < 2e-3, ra.mr_or_per_SD
    assert int(g.loc["Pain symptoms", "n_endpoints_tested"]) == 28
    assert int(g.loc["Degenerative spine", "n_endpoints_tested"]) == 29
    assert int(g.loc["Musculoskeletal, non-pain", "n_endpoints_tested"]) == 224
    assert int(g.loc["Bone / osteoporosis", "n_endpoints_tested"]) == 35
    assert int(g.loc["Nervous system, non-pain", "n_endpoints_tested"]) == 131, (
        f"expected 131 non-pain nervous-system endpoints, got "
        f"{int(g.loc['Nervous system, non-pain', 'n_endpoints_tested'])}; the "
        "PAIN_CODES_RE pattern must classify causalgia as pain"
    )
    assert int(g.loc["Pain symptoms", "obs"]) == 18
    assert abs(g.loc["Pain symptoms", "fold"] - 3.8) < 0.2, g.loc["Pain symptoms", "fold"]
    assert g.loc["Pain symptoms", "p"] <= 1 / 20001 + 1e-12, g.loc["Pain symptoms", "p"]
    for neg in ["Degenerative spine", "Musculoskeletal, non-pain", "Bone / osteoporosis",
                "Nervous system, non-pain"]:
        assert g.loc[neg, "p"] > 0.05, (neg, g.loc[neg, "p"])
    assert dm.loc[dm.spine_struct, "neff"].median() > dm.loc[dm.pain_symptom, "neff"].median(), (
        "the spine control must have MORE power than the pain endpoints for its "
        "null result to be meaningful"
    )
    print(f"\nwrote three tables to {OUT}")
    return d, fam, enr


if __name__ == "__main__":
    main()
