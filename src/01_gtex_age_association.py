"""
Figure 5A, 5B: age association of ZBTB48 across GTEx v10 tissues.

Regresses sample-level log2(TPM + 1) on donor age-bracket midpoint in every
GTEx v10 tissue with at least four age brackets and 40 samples, and applies
Benjamini-Hochberg correction across all tissues analysed.

Two nervous-system tissues reach q < 0.05 and they differ in direction:
tibial nerve increases with age, and cerebellum decreases. Both are reported;
the figure should mark both.

Output: results/zbtb48_age_gtex_trends.csv
Runtime: about one minute, dominated by one large API response.
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (  # noqa: E402
    AGE_BRACKET_MIDPOINT,
    ZBTB48_GENCODE_ID,
    gtex_get,
)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")


def main():
    records = gtex_get(
        "expression/geneExpression",
        {
            "gencodeId": ZBTB48_GENCODE_ID,
            "attributeSubset": "ageBracket",
            "itemsPerPage": 100000,
        },
    )
    records = [x for x in records if len(x.get("data") or []) > 0]
    if not records:
        raise RuntimeError(
            "GTEx returned no expression data. The most likely cause is a "
            "GENCODE identifier that does not match the pinned dataset release."
        )

    long_rows = []
    for x in records:
        values = np.asarray(x["data"], dtype=float)
        age = AGE_BRACKET_MIDPOINT[x["subsetGroup"]]
        long_rows += [(x["tissueSiteDetailId"], age, t) for t in values]
    lg = pd.DataFrame(long_rows, columns=["tissue", "age_mid", "tpm"])
    lg["l2"] = np.log2(lg.tpm + 1)

    fits = []
    for tissue, d in lg.groupby("tissue"):
        # Four brackets is the minimum for a slope that is not driven by two
        # points; 40 samples keeps single-donor-dominated tissues out.
        if d.age_mid.nunique() < 4 or len(d) < 40:
            continue
        slope, _, _, p, se = stats.linregress(d.age_mid, d.l2)
        rho, p_rho = stats.spearmanr(d.age_mid, d.l2)
        fits.append(
            dict(
                tissue=tissue,
                n=len(d),
                slope_l2_per_decade=slope * 10,
                se_per_decade=se * 10,
                pct_per_decade=(2 ** (slope * 10) - 1) * 100,
                p_linreg=p,
                spearman_rho=rho,
                spearman_p=p_rho,
                median_tpm=d.tpm.median(),
            )
        )
    f = pd.DataFrame(fits).sort_values("slope_l2_per_decade").reset_index(drop=True)

    # Correction is across every tissue analysed, not across the nervous-system
    # subset shown in the figure. Correcting only within the displayed subset
    # would be a smaller and less honest denominator.
    f["q_bh"] = multipletests(f.p_linreg, method="fdr_bh")[1]

    def tissue_class(t):
        tl = t.lower()
        if tl.startswith("brain") or "nerve" in tl or "pituitary" in tl:
            return "nervous system"
        if t in ("Whole_Blood", "Spleen", "Cells_EBV-transformed_lymphocytes"):
            return "blood / immune"
        return "other tissue"

    f["tissue_class"] = f.tissue.map(tissue_class)
    f["tissue_label"] = f.tissue.str.replace("_", " ")

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "zbtb48_age_gtex_trends.csv")
    f.to_csv(path, index=False)

    # --- checks against the published values ---------------------------------
    nerve = f[f.tissue == "Nerve_Tibial"]
    assert len(nerve) == 1, "tibial nerve missing from the fit table"
    nerve = nerve.iloc[0]
    nervous_sig = f[(f.tissue_class == "nervous system") & (f.q_bh < 0.05)]

    print(f"tissues fitted                {len(f)}")
    print(f"tibial nerve n                {nerve.n}                (published 670)")
    print(f"tibial nerve % per decade     {nerve.pct_per_decade:+.2f}            (published +2.77)")
    print(f"tibial nerve P (two-sided)    {nerve.p_linreg:.3e}        (published 3.39e-8)")
    print(f"tibial nerve q (BH)           {nerve.q_bh:.3e}")
    print(f"nervous-system tissues q<0.05 {len(nervous_sig)}")
    for _, r in nervous_sig.iterrows():
        print(f"    {r.tissue_label:36s} {r.pct_per_decade:+6.2f} %/decade  q={r.q_bh:.2e}")
    print(f"\nwrote {path}")

    assert nerve.n == 670, f"expected n=670 for tibial nerve, got {nerve.n}"
    assert abs(nerve.pct_per_decade - 2.77) < 0.05, nerve.pct_per_decade
    assert nerve.p_linreg < 1e-7, nerve.p_linreg
    assert len(nervous_sig) == 2, (
        f"expected 2 nervous-system tissues at q<0.05, got {len(nervous_sig)}; "
        "the figure must mark both, and they differ in direction"
    )
    return f


if __name__ == "__main__":
    main()
