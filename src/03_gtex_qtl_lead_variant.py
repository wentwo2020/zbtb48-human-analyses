"""
Extended Data Fig. 7D: the instrument, and whether it is the lead variant.

Pulls every ZBTB48 cis-eQTL and cis-sQTL association in GTEx v10 and asks, per
tissue, whether rs17029626 is merely detected or is the lead variant for the
gene.

The distinction matters and is easy to get wrong. rs17029626 is a significant
cis-eQTL in 16 of the 23 tissues with a detected ZBTB48 eQTL, but it is the
LEAD variant in only 12 of them: it is out-ranked in aorta, lung, oesophageal
muscularis and breast. On the splicing side it is significant in all 35
tissues with a detected sQTL and the lead in 34, the exception being bladder.
Reporting "detected" as "lead" overstates the instrument, and a referee can
falsify it in one query.

Note also that the splicing signal comprises two junctions in the same cluster
with opposite effect directions, which is a competing-donor-site pattern
rather than a single sQTL, and is a stronger mechanistic statement than a
single junction would be.

Output: results/zbtb48_eqtl_lead_by_tissue_v10.csv
        results/zbtb48_sqtl_lead_by_tissue_v10.csv
        results/zbtb48_instrument_strength.csv
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats as st

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (  # noqa: E402
    INSTRUMENT_VARIANT_B38,
    ZBTB48_GENCODE_ID,
    gtex_get_paged,
)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")


def lead_status(assoc, label):
    """Per tissue: is the instrument present, and is it the lead variant?

    "Lead" means lowest p-value for this gene in this tissue among all
    associations GTEx returns.
    """
    rows = []
    for tissue, g in assoc.groupby("tissueSiteDetailId"):
        best = g.pValue.min()
        here = g[g.variantId == INSTRUMENT_VARIANT_B38]
        present = len(here) > 0
        p_here = here.pValue.min() if present else np.nan
        rows.append(
            dict(
                tissue=tissue,
                n_assoc=len(g),
                variant_present=present,
                variant_p=p_here,
                best_p=best,
                is_lead=bool(present and p_here == best),
                rank=(int((g.pValue < p_here).sum()) + 1) if present else np.nan,
            )
        )
    r = pd.DataFrame(rows).sort_values("variant_p")
    print(f"\n{label}")
    print(f"  tissues with a detected ZBTB48 signal   {len(r)}")
    print(f"  instrument is a significant association {int(r.variant_present.sum())}")
    print(f"  instrument is the LEAD variant          {int(r.is_lead.sum())}")
    notlead = r[r.variant_present & ~r.is_lead]
    if len(notlead):
        print("  detected but out-ranked:")
        for _, x in notlead.sort_values("rank").iterrows():
            print(
                f"      {x.tissue:26s} rank {int(x['rank']):3d}   "
                f"its P {x.variant_p:.2e}   lead P {x.best_p:.2e}"
            )
    return r


def main():
    eqtl = pd.DataFrame(
        gtex_get_paged("association/singleTissueEqtl", {"gencodeId": ZBTB48_GENCODE_ID})
    )
    sqtl = pd.DataFrame(
        gtex_get_paged("association/singleTissueSqtl", {"gencodeId": ZBTB48_GENCODE_ID})
    )
    if eqtl.empty or sqtl.empty:
        raise RuntimeError(
            "GTEx returned no associations. Check that ZBTB48_GENCODE_ID matches "
            "the GENCODE release used by the pinned dataset."
        )

    re_ = lead_status(eqtl, "cis-eQTL")
    rs_ = lead_status(sqtl, "cis-sQTL")

    os.makedirs(OUT, exist_ok=True)
    re_.to_csv(os.path.join(OUT, "zbtb48_eqtl_lead_by_tissue_v10.csv"), index=False)
    rs_.to_csv(os.path.join(OUT, "zbtb48_sqtl_lead_by_tissue_v10.csv"), index=False)

    # Instrument strength. F is recovered from the eQTL p-value rather than
    # from an F-test on individual-level data, which is not available: for a
    # single instrument F equals the squared z statistic.
    inst = []
    for label, assoc, tissues in [
        ("eQTL", eqtl, ["Whole_Blood", "Nerve_Tibial"]),
        ("sQTL", sqtl, ["Whole_Blood"]),
    ]:
        for t in tissues:
            g = assoc[(assoc.tissueSiteDetailId == t) & (assoc.variantId == INSTRUMENT_VARIANT_B38)]
            if not len(g):
                continue
            for _, x in g.sort_values("pValue").iterrows():
                z = st.norm.isf(x.pValue / 2)
                inst.append(
                    dict(
                        qtl_type=label,
                        tissue=t,
                        phenotype_id=x.get("phenotypeId", ""),
                        nes=x.get("nes", np.nan),
                        pvalue=x.pValue,
                        z=z,
                        F_equals_z_squared=z ** 2,
                    )
                )
    ins = pd.DataFrame(inst)
    ins.to_csv(os.path.join(OUT, "zbtb48_instrument_strength.csv"), index=False)

    blood_e = ins[(ins.qtl_type == "eQTL") & (ins.tissue == "Whole_Blood")].iloc[0]
    print("\ninstrument strength")
    print(
        f"  whole-blood eQTL NES {blood_e.nes:+.4f}  P {blood_e.pvalue:.2e}  "
        f"F = z^2 = {blood_e.F_equals_z_squared:.1f}    (published F approx 87)"
    )
    blood_s = ins[(ins.qtl_type == "sQTL") & (ins.tissue == "Whole_Blood")]
    print(f"  whole-blood sQTL junctions returned: {len(blood_s)}")
    for _, x in blood_s.iterrows():
        print(f"      {x.phenotype_id}  NES {x.nes:+.4f}  P {x.pvalue:.2e}")

    # --- checks against the published values ---------------------------------
    assert len(re_) == 23, f"expected 23 eQTL tissues, got {len(re_)}"
    assert int(re_.variant_present.sum()) == 16, int(re_.variant_present.sum())
    assert int(re_.is_lead.sum()) == 12, (
        f"expected the instrument to lead in 12 eQTL tissues, got {int(re_.is_lead.sum())}; "
        "do not report 'detected' as 'lead'"
    )
    assert len(rs_) == 35, f"expected 35 sQTL tissues, got {len(rs_)}"
    assert int(rs_.variant_present.sum()) == 35, int(rs_.variant_present.sum())
    assert int(rs_.is_lead.sum()) == 34, (
        f"expected the instrument to lead in 34 of 35 sQTL tissues, got {int(rs_.is_lead.sum())}; "
        "it is not the lead in bladder, so 'all 35' is false"
    )
    assert abs(blood_e.F_equals_z_squared - 87) < 3, blood_e.F_equals_z_squared
    assert len(blood_s) >= 2, (
        "expected at least two whole-blood sQTL junctions; the opposite-direction "
        "pair is the competing-donor-site result"
    )
    print(f"\nwrote {OUT}/zbtb48_eqtl_lead_by_tissue_v10.csv and two more")
    return re_, rs_, ins


if __name__ == "__main__":
    main()
