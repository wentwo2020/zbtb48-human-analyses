"""
Extended Data Fig. 7A, 7B: is the tibial-nerve age effect an artefact?

Two sensitivity analyses on the Figure 5B result.

7A. A 23-gene reference panel is regressed on age in tibial nerve alongside
ZBTB48: cell-type markers for the populations that could plausibly shift with
age in a peripheral nerve biopsy, plus four housekeeping genes. The payload of
this panel is usually read as "the housekeeping genes are flat", but that is
only half true and the half that is true is the important half: ACTB and
RPL13A are stable with age, while TBP and GAPDH are not. It is the stability
of ACTB and RPL13A that excludes a global normalisation drift, and ACTB is the
normaliser used for the single-nucleus panels, so its stability is
load-bearing rather than incidental.

7B. A composition test. For each of three cell populations that could dilute
or enrich the nerve biopsy with age, the analysis asks how much of the
observed ZBTB48 age effect could be produced purely by that population's
fraction changing, given the marker's own age slope and the population's
ZBTB48 level in a reference tissue. The three predicted contributions are
small and one is negative; combined they predict -0.12 per cent per decade
against an observed +2.77, so composition shift does not account for the
result and in fact points the other way.

Output: results/zbtb48_nerve_reference_genes.csv
        results/zbtb48_nerve_composition_tests.csv
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
    gtex_get,
    resolve_gencode_ids,
)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")

NERVE_PANEL = {
    "ZBTB48": "target",
    "ADIPOQ": "adipocyte", "PLIN1": "adipocyte", "LEP": "adipocyte",
    "FABP4": "adipocyte", "CFD": "adipocyte",
    "MPZ": "Schwann / myelin", "PMP22": "Schwann / myelin",
    "S100B": "Schwann / myelin", "SOX10": "Schwann / myelin",
    "EGR2": "Schwann / myelin",
    "COL1A1": "fibroblast / ECM", "COL3A1": "fibroblast / ECM",
    "DCN": "fibroblast / ECM",
    "PTPRC": "immune", "CD68": "immune", "C1QA": "immune",
    "PECAM1": "endothelial", "VWF": "endothelial",
    "ACTB": "housekeeping", "GAPDH": "housekeeping",
    "TBP": "housekeeping", "RPL13A": "housekeeping",
}

# Each population is paired with a reference tissue in which it is
# effectively pure, so that its ZBTB48 level can be read off directly.
COMPOSITION_TESTS = [
    ("ADIPOQ", "Adipose_Subcutaneous", "adipocyte"),
    ("PTPRC", "Whole_Blood", "leukocyte"),
    ("COL1A1", "Cells_Cultured_fibroblasts", "fibroblast"),
]


def fit_panel():
    sym2gid = resolve_gencode_ids(list(NERVE_PANEL))
    missing = set(NERVE_PANEL) - set(sym2gid)
    if missing:
        raise RuntimeError(f"could not resolve GENCODE ids for {sorted(missing)}")
    gid2sym = {v: k for k, v in sym2gid.items()}

    data = gtex_get(
        "expression/geneExpression",
        {
            "gencodeId": list(sym2gid.values()),
            "tissueSiteDetailId": "Nerve_Tibial",
            "attributeSubset": "ageBracket",
            "itemsPerPage": 100000,
        },
    )
    rows = []
    for x in data:
        if not x.get("data"):
            continue
        sym = gid2sym[x["gencodeId"]]
        rows += [
            (sym, NERVE_PANEL[sym], AGE_BRACKET_MIDPOINT[x["subsetGroup"]], t)
            for t in x["data"]
        ]
    nb = pd.DataFrame(rows, columns=["gene", "category", "age_mid", "tpm"])
    nb["l2"] = np.log2(nb.tpm + 1)

    out = []
    for gene, d in nb.groupby("gene"):
        slope, _, _, p, se = stats.linregress(d.age_mid, d.l2)
        out.append(
            dict(
                gene=gene,
                category=NERVE_PANEL[gene],
                n=len(d),
                slope_l2_per_decade=slope * 10,
                se_l2_per_decade=se * 10,
                pct_per_decade=(2 ** (slope * 10) - 1) * 100,
                p=p,
                median_tpm=d.tpm.median(),
            )
        )
    nf = pd.DataFrame(out).sort_values("pct_per_decade").reset_index(drop=True)
    nf["q_bh"] = multipletests(nf.p, method="fdr_bh")[1]
    return nf, sym2gid, gid2sym


def median_matrix(sym2gid, gid2sym, genes, tissues):
    data = gtex_get(
        "expression/medianGeneExpression",
        {
            "gencodeId": [sym2gid[g] for g in genes],
            "tissueSiteDetailId": list(tissues),
            "itemsPerPage": 200,
        },
        timeout=90,
    )
    return pd.DataFrame(
        [
            {
                "gene": gid2sym[x["gencodeId"]],
                "tissue": x["tissueSiteDetailId"],
                "median_tpm": x["median"],
            }
            for x in data
        ]
    ).pivot(index="gene", columns="tissue", values="median_tpm")


def main():
    nf, sym2gid, gid2sym = fit_panel()
    os.makedirs(OUT, exist_ok=True)
    nf.to_csv(os.path.join(OUT, "zbtb48_nerve_reference_genes.csv"), index=False)
    z = nf[nf.gene == "ZBTB48"].iloc[0]

    immune = median_matrix(
        sym2gid, gid2sym,
        ["ZBTB48", "PTPRC", "COL1A1"],
        ["Whole_Blood", "Nerve_Tibial", "Cells_Cultured_fibroblasts"],
    )
    adipose = median_matrix(
        sym2gid, gid2sym,
        ["ZBTB48", "ADIPOQ", "FABP4", "PLIN1"],
        ["Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Nerve_Tibial"],
    )
    med = pd.concat([adipose, immune.drop(columns=["Nerve_Tibial"])], axis=1)

    # Use the age-bracket medians from the panel fit for the nerve column, so
    # that the nerve values come from the same samples as the slopes.
    nerve_med = nf.set_index("gene").median_tpm
    for g in ["PTPRC", "COL1A1", "ADIPOQ", "ZBTB48"]:
        med.loc[g, "Nerve_Tibial"] = nerve_med[g]

    def composition_test(marker, source_tissue, label):
        # Implied fraction of the nerve biopsy contributed by this population,
        # from the marker's nerve level relative to its level in a pure source.
        f0 = med.loc[marker, "Nerve_Tibial"] / med.loc[marker, source_tissue]
        # How fast that fraction appears to change with age, from the marker.
        d_rel = nf.loc[nf.gene == marker, "pct_per_decade"].iloc[0] / 100
        df = f0 * d_rel
        a = med.loc["ZBTB48", source_tissue]
        n_obs = med.loc["ZBTB48", "Nerve_Tibial"]
        # ZBTB48 level in the nerve with this population's contribution removed.
        n_pure = (n_obs - f0 * a) / (1 - f0)
        pred = df * (a - n_pure) / n_obs * 100
        need = (z.pct_per_decade / 100) * n_obs / (a - n_pure)
        return dict(
            population=label,
            marker=marker,
            reference_tissue=source_tissue,
            marker_pct_per_decade_in_nerve=nf.loc[nf.gene == marker, "pct_per_decade"].iloc[0],
            implied_fraction_pct=f0 * 100,
            fraction_change_pp_per_decade=df * 100,
            zbtb48_tpm_in_population=a,
            zbtb48_tpm_in_nerve=n_obs,
            predicted_pct_per_decade=pred,
            pct_of_observed_explained=abs(pred / z.pct_per_decade) * 100,
            fraction_change_needed_pp_per_decade=need * 100,
        )

    ct = pd.DataFrame([composition_test(*t) for t in COMPOSITION_TESTS])
    ct.to_csv(os.path.join(OUT, "zbtb48_nerve_composition_tests.csv"), index=False)

    # --- checks against the published values ---------------------------------
    hk = nf[nf.category == "housekeeping"].set_index("gene")
    print(f"panel genes fitted            {len(nf)}                  (published 23)")
    print(f"all at n = 670                {bool((nf.n == 670).all())}")
    print(f"genes at q < 0.05             {int((nf.q_bh < 0.05).sum())} of {len(nf)}       (published 17 of 23)")
    print("\nhousekeeping genes, which are NOT uniformly stable:")
    for g in ["ACTB", "RPL13A", "TBP", "GAPDH"]:
        r = hk.loc[g]
        flag = "stable" if r.q_bh >= 0.05 else "CHANGES with age"
        print(f"    {g:8s} {r.pct_per_decade:+6.2f} %/decade  q={r.q_bh:.3f}   {flag}")
    print("\ncomposition contributions to the observed +2.77 %/decade:")
    for _, r in ct.iterrows():
        print(f"    {r.population:12s} {r.predicted_pct_per_decade:+6.2f} %/decade")
    print(f"    {'combined':12s} {ct.predicted_pct_per_decade.sum():+6.2f} %/decade   (published -0.12)")
    print(f"\nwrote {OUT}/zbtb48_nerve_reference_genes.csv and .../zbtb48_nerve_composition_tests.csv")

    assert len(nf) == 23, len(nf)
    assert (nf.n == 670).all()
    assert hk.loc["ACTB", "q_bh"] >= 0.05, "ACTB must be age-stable for the per-ACTB panels to hold"
    assert abs(ct.predicted_pct_per_decade.sum() + 0.12) < 0.05, ct.predicted_pct_per_decade.sum()
    return nf, ct


if __name__ == "__main__":
    main()
