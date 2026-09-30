"""
Figure 5D: astrocyte and microglial ZBTB48 across four independent human
spinal cord single-nucleus datasets.

This script handles the two GEO datasets. The CZ CELLxGENE bar comes from
script 04.

BOTH datasets use the depositors' own cell-type annotations, verbatim. No
re-clustering, no marker-based reassignment, no label transfer, and no
cell-level filtering beyond what the depositors applied. Two annotation quirks
are preserved rather than corrected, because a re-run has to match them:
GSE190442 spells its astrocyte label "Astrcytes", and GSE330130 has a double
space in "Excitatory  Neuron".

WHY THE METRIC MATTERS. ZBTB48 is expressed at low level and detected in only
a few per cent of nuclei, and ACTB is roughly twice as abundant in microglia
as in astrocytes in these datasets. Per-library and per-ACTB normalisation
therefore differ in DIRECTION, not just in magnitude: on raw counts per
million microglia look higher, and on ZBTB48 per 1,000 ACTB astrocytes are
higher. This is not an artefact to be hidden, it is the reason the
normalisation was chosen, and the ACTB level per cell type is exported
alongside every ratio so the choice is auditable. Do not extrapolate an ACTB
ratio measured in one dataset to another platform or region; measure it in
each.

TWO THINGS TO DISCLOSE when reporting this panel. GSE330130 comprises 11
donors with amyotrophic lateral sclerosis and 9 controls, and ALS lumbar
spinal cord carries marked astrogliosis and microgliosis, so the control-only
comparison is the primary one and the two GSE330130 bars are nested rather
than independent. And the within-donor paired test reaches significance in
only one of the four datasets, so the panel is a consistent direction of
effect across four datasets anchored by one individually significant paired
comparison, not four independent replications.

Output: results/gse190442_per1000ACTB_by_celltype.csv
        results/gse190442_per1000ACTB_per_donor.csv
        results/gse330130_per1000ACTB_by_celltype.csv
        results/gse330130_per1000ACTB_per_donor.csv
        results/fig5d_reconciliation.csv
Downloads: about 0.65 GB in total, cached in data/.
Runtime: a few minutes, dominated by the downloads.
"""

import csv
import gzip
import os
import sys
from urllib.request import urlretrieve

import numpy as np
import pandas as pd
from scipy import stats as st

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import GEO_FTP  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
DATA = os.path.join(HERE, "..", "data")

GENES = {"ZBTB48", "ACTB", "TBP", "GAPDH", "RPL13A"}

GSE190442 = {
    "counts": f"{GEO_FTP}/GSE190nnn/GSE190442/suppl/GSE190442_aggregated_counts_postqc.csv.gz",
    "meta": f"{GEO_FTP}/GSE190nnn/GSE190442/suppl/GSE190442_aggregated_metadata_postqc.csv.gz",
}
GSE330130_H5AD = (
    f"{GEO_FTP}/GSE330nnn/GSE330130/suppl/GSE330130_LSC_SNRNASEQ_COUNTS.h5ad"
)

# Depositors' labels, verbatim, including the GSE190442 misspelling.
MAIN_190442 = ["Astrcytes", "Microglia", "Oligodendrocytes", "OPC", "Neurons"]
MAIN_330130 = [
    "Astrocyte", "Microglia", "Oligodendrocytes", "OPCs",
    "Inhibitory Neuron", "Excitatory  Neuron", "Motor Neurons",
]
ASTRO_190442, ASTRO_330130 = "Astrcytes", "Astrocyte"

MIN_NUCLEI_PER_DONOR_CELLTYPE = 50


def cache(url, name):
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, name)
    if not os.path.exists(path):
        print(f"downloading {name} ...")
        urlretrieve(url, path)
    return path


def per_actb(df, populations, label):
    """Pooled per-cell-type summary, with ACTB reported alongside the ratio."""
    return (
        df[df.ct.isin(populations)]
        .groupby("ct")
        .apply(
            lambda g: pd.Series(
                {
                    "cohort": label,
                    "n_cells": len(g),
                    "donors": g.donor.nunique(),
                    "ZBTB48_per_1000_ACTB": 1000 * g.ZBTB48.sum() / g.ACTB.sum(),
                    "ACTB_CPM": 1e6 * g.ACTB.sum() / g.lib.sum(),
                    "ZBTB48_CPM": 1e6 * g.ZBTB48.sum() / g.lib.sum(),
                    "detection_pct": 100 * (g.ZBTB48 > 0).mean(),
                }
            ),
            include_groups=False,
        )
        .sort_values("ZBTB48_per_1000_ACTB", ascending=False)
    )


def per_donor(df, astro, micro="Microglia"):
    """Donor-level per-ACTB values for the two populations being compared."""
    return (
        df[df.ct.isin([astro, micro])]
        .groupby(["donor", "ct"])
        .apply(
            lambda g: (1000 * g.ZBTB48.sum() / g.ACTB.sum())
            if (g.ACTB.sum() > 0 and len(g) >= MIN_NUCLEI_PER_DONOR_CELLTYPE)
            else np.nan,
            include_groups=False,
        )
        .unstack()
    )


def load_gse190442():
    counts = cache(GSE190442["counts"], "GSE190442_counts.csv.gz")
    meta_path = cache(GSE190442["meta"], "GSE190442_metadata.csv.gz")

    kept = {}
    with gzip.open(counts, "rt", newline="") as f:
        rd = csv.reader(f)
        cells = next(rd)[1:]
        for row in rd:
            if row[0] in GENES:
                kept[row[0]] = np.asarray(row[1:], dtype=np.int32)
    sc = pd.DataFrame(kept, index=cells).T
    # The deposited counts CSV has R-mangled barcode headers, in which the
    # hyphen of each 10x suffix became a period. Repair before joining; the
    # join then has to match every barcode, which is the check that the
    # repair was right.
    sc.columns = [c.replace(".1_", "-1_") for c in sc.columns]

    md = pd.read_csv(meta_path, index_col=0, low_memory=False)
    assert set(sc.columns) <= set(md.index), "barcode repair failed: unmatched cells"
    md = md.loc[sc.columns]

    return pd.DataFrame(
        {
            "ct": md.top_level_annotation.values,
            "donor": md["sample"].values,
            "lib": md.nCount_RNA.values.astype(float),
            "ZBTB48": sc.loc["ZBTB48"].values.astype(float),
            "ACTB": sc.loc["ACTB"].values.astype(float),
        }
    )


def load_gse330130():
    import h5py

    path = cache(GSE330130_H5AD, "GSE330130_LSC_SNRNASEQ_COUNTS.h5ad")
    f = h5py.File(path, "r")

    def categorical(name):
        g = f[f"obs/{name}"]
        cats = np.array(
            [s.decode() if isinstance(s, bytes) else s for s in g["categories"][:]]
        )
        return cats[g["codes"][:]]

    # The parent series also contains motor cortex. Verify rather than assume
    # that this file is spinal cord only; that is why no tissue filter is
    # applied downstream.
    tissues = set(categorical("Tissue Type"))
    assert tissues == {"Lumbar Spinal Cord"}, f"unexpected tissues: {tissues}"

    var_names = np.array(
        [s.decode() if isinstance(s, bytes) else s for s in f["var/common_name/categories"][:]]
    )[f["var/common_name/codes"][:]]
    targets = {g: i for i, g in enumerate(var_names) if g in GENES}

    indptr = f["X/indptr"][:]
    nnz = int(indptr[-1])
    n_cell = len(indptr) - 1
    out = {g: np.zeros(n_cell, dtype=np.float32) for g in targets}
    cols = np.array(sorted(targets.values()))
    col2gene = {v: k for k, v in targets.items()}

    # Stream the CSR matrix in chunks; the full nonzero block is too large to
    # hold alongside everything else.
    for start in range(0, nnz, 30_000_000):
        end = min(start + 30_000_000, nnz)
        idx = f["X/indices"][start:end]
        dat = f["X/data"][start:end]
        hit = np.isin(idx, cols)
        if hit.any():
            positions = np.nonzero(hit)[0] + start
            rows = np.searchsorted(indptr, positions, side="right") - 1
            for c, r, v in zip(idx[hit], rows, dat[hit]):
                out[col2gene[c]][r] = v
        del idx, dat, hit

    df = pd.DataFrame(
        {
            "ct": categorical("Cell_Type"),
            "donor": categorical("Subject ID"),
            "cohort": categorical("Level 1"),
            "lib": f["obs/n_counts"][:],
        }
    )
    f.close()
    for k, v in out.items():
        df[k] = v
    return df


def main():
    os.makedirs(OUT, exist_ok=True)

    # --- GSE190442 -----------------------------------------------------------
    a = load_gse190442()
    sa = per_actb(a, MAIN_190442, "all 7 donors")
    wa = per_donor(a, ASTRO_190442)
    sa.reset_index().to_csv(os.path.join(OUT, "gse190442_per1000ACTB_by_celltype.csv"), index=False)
    wa.reset_index().to_csv(os.path.join(OUT, "gse190442_per1000ACTB_per_donor.csv"), index=False)

    # --- GSE330130 -----------------------------------------------------------
    b = load_gse330130()
    sb_all = per_actb(b, MAIN_330130, "all (11 ALS + 9 control)")
    sb_ctl = per_actb(b[b.cohort == "Control"], MAIN_330130, "control only")
    wb = per_donor(b, ASTRO_330130)
    wb["cohort"] = b.groupby("donor").cohort.first()
    pd.concat([sb_all.reset_index(), sb_ctl.reset_index()]).to_csv(
        os.path.join(OUT, "gse330130_per1000ACTB_by_celltype.csv"), index=False
    )
    wb.reset_index().to_csv(os.path.join(OUT, "gse330130_per1000ACTB_per_donor.csv"), index=False)

    # --- reconcile all four bars --------------------------------------------
    def ratio(summary, astro):
        return (
            summary.loc[astro, "ZBTB48_per_1000_ACTB"]
            / summary.loc["Microglia", "ZBTB48_per_1000_ACTB"],
            summary.loc[astro, "ZBTB48_CPM"] / summary.loc["Microglia", "ZBTB48_CPM"],
            summary.loc["Microglia", "ACTB_CPM"] / summary.loc[astro, "ACTB_CPM"],
        )

    def paired(w, astro):
        p = w[[astro, "Microglia"]].dropna()
        s = st.wilcoxon(p[astro], p.Microglia)
        return (
            len(p),
            (p[astro] / p.Microglia).median(),
            int((p[astro] > p.Microglia).sum()),
            s.pvalue,
        )

    rows = []
    for label, published, summary, astro, w in [
        ("GSE190442 (lumbar, 7 donors)", 1.24, sa, ASTRO_190442, wa),
        ("GSE330130 all (11 ALS + 9 control, lumbar)", 1.30, sb_all, ASTRO_330130, wb),
        ("GSE330130 control only (lumbar, 9 donors)", 1.68, sb_ctl, ASTRO_330130,
         wb[wb.cohort == "Control"]),
    ]:
        r_actb, r_cpm, r_actb_level = ratio(summary, astro)
        n, med, higher, pv = paired(w, astro)
        rows.append(
            dict(
                dataset=label,
                figure_5D_value=published,
                recomputed_per_1000_ACTB=r_actb,
                raw_CPM_astro_over_micro=r_cpm,
                ACTB_micro_over_astro=r_actb_level,
                n_paired=n,
                median_paired_ratio=med,
                donors_astro_higher=higher,
                paired_wilcoxon_P=pv,
            )
        )
    fig = pd.DataFrame(rows)
    fig.to_csv(os.path.join(OUT, "fig5d_reconciliation.csv"), index=False)

    print("\n=== Figure 5D: recomputed against the plotted values ===")
    print(fig.round(3).to_string(index=False))
    print(
        f"\nACTB is {fig.ACTB_micro_over_astro.min():.2f} to "
        f"{fig.ACTB_micro_over_astro.max():.2f}-fold higher in microglia than astrocytes, "
        "which is why raw CPM and per-ACTB differ in direction"
    )
    n_sig = int((fig.paired_wilcoxon_P < 0.05).sum())
    print(
        f"within-donor paired test significant in {n_sig} of {len(fig)} GEO comparisons; "
        "report the panel as a consistent direction, not as independent replications"
    )
    print(
        f"GSE330130 cohort: {b[b.cohort == 'ALS'].donor.nunique()} ALS + "
        f"{b[b.cohort == 'Control'].donor.nunique()} control donors -- disclose this in the legend"
    )

    # --- checks against the published values ---------------------------------
    assert len(a) == 55289, f"GSE190442 should have 55,289 nuclei, got {len(a)}"
    assert len(b) == 121881, f"GSE330130 should have 121,881 nuclei, got {len(b)}"
    assert b[b.cohort == "ALS"].donor.nunique() == 11
    assert b[b.cohort == "Control"].donor.nunique() == 9
    assert (fig.recomputed_per_1000_ACTB > 1).all(), (
        "all three GEO bars must be astrocyte-favouring on the per-ACTB scale"
    )
    assert (fig.raw_CPM_astro_over_micro < 1).all(), (
        "and all three must be microglia-favouring on raw CPM; if this fails the "
        "normalisation no longer explains the direction difference"
    )
    rel = ((fig.recomputed_per_1000_ACTB - fig.figure_5D_value).abs() / fig.figure_5D_value)
    assert rel.max() < 0.15, f"recomputed ratios drift from plotted by {rel.max():.1%}"
    print(f"\nwrote five tables to {OUT}")
    return fig


if __name__ == "__main__":
    main()
