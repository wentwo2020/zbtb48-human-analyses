"""
Figure 5C: ZBTB48 across cell populations of the human spinal cord, from the
CZ CELLxGENE Discover Census.

This script does two things, in order.

FIRST it re-derives the full Census pull: the spinal-cord query, the grouped
sums per dataset x donor x cell type, and the eight-row cell-type summary over
all three deposited datasets. That configuration contains 42,186 nuclei from
21 donors.

SECOND it produces the panel that is actually published, which is NOT that
configuration. Pooling the three datasets is not defensible for a cell-type
comparison, because two of the three contribute a single donor each, one of
them has no annotated microglia at all and the other contributes only six
astrocyte nuclei. The published panel is therefore restricted to the one
dataset with genuine donor replication, and every comparison is paired within
donor. Both configurations are written out so the difference is auditable.

One upstream detail is worth repeating because it caused real confusion: the
Census query filters on `tissue_general`, and the original pipeline did not
carry the specific `tissue` label into the group key, which made the spinal
region unrecoverable from the saved table. This script keeps `tissue` in the
group key.

Requires: cellxgene-census==1.18.0, tiledbsoma==2.3.0 (Python 3.11)
Note:     census_version is PINNED. "stable" advances and will drift.

Output: results/fig5c_reproduced_per_donor_celltype.csv
        results/fig5c_reproduced_celltype_summary.csv
        results/fig5c_donor_by_celltype_matrix.csv   (the published panel)
        results/fig5c_paired_statistics.csv
"""
import numpy as np, pandas as pd, tiledbsoma as soma, cellxgene_census, re, os

CENSUS_VERSION = "2025-11-08"          # what "stable" resolved to on the original run
VALUE_FILTER   = "is_primary_data == True and tissue_general == 'spinal cord'"
MARK = ["ZBTB48","P2RY12","CX3CR1","TMEM119","AIF1","GFAP","AQP4","SLC1A2","SNAP25",
        "RBFOX3","MBP","PLP1","PDGFRA","ACTB","GAPDH","TBP","RPL13A"]
CT_MAP = {"microglial cell":"Microglia","central nervous system macrophage":"CNS macrophage",
          "astrocyte":"Astrocyte","oligodendrocyte":"Oligodendrocyte",
          "oligodendrocyte precursor cell":"OPC","neuron":"Neuron",
          "GABAergic neuron":"Inhibitory neuron","glutamatergic neuron":"Excitatory neuron"}

def age_years(stage):
    m = re.search(r"(\d+)", str(stage) or "")
    return float(m.group(1)) if m else np.nan

# If the sandbox routes S3 through a proxy, tiledbsoma needs it passed explicitly.
cfg = {}
px = os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy")
if px:
    hp = px.split("://")[-1].rstrip("/")
    host_, _, port_ = hp.partition(":")
    cfg = {"vfs.s3.proxy_host": host_, "vfs.s3.proxy_port": port_ or "80",
           "vfs.s3.region": "us-west-2", "vfs.s3.no_sign_request": "true"}

census = cellxgene_census.open_soma(census_version=CENSUS_VERSION,
                                    tiledb_config=cfg or None)
h = census["census_data"]["homo_sapiens"]

obs = h.obs.read(column_names=["soma_joinid","tissue","cell_type","development_stage",
                               "donor_id","dataset_id","assay","suspension_type","disease"],
                 value_filter=VALUE_FILTER).concat().to_pandas()
obs["age"] = obs.development_stage.map(age_years)
for c in ["tissue","cell_type","donor_id","dataset_id","assay","suspension_type"]:
    obs[c] = obs[c].astype(str)
print(f"spinal nuclei pulled: {len(obs):,} | datasets: {obs.dataset_id.nunique()} | donors: {obs.donor_id.nunique()}")

vm = h.ms["RNA"].var.read(value_filter=f"feature_name in {MARK}").concat().to_pandas()
jid2sym = dict(zip(vm.soma_joinid.astype(int), vm.feature_name.astype(str)))

# NOTE: the original group key omits `tissue`. Keeping it here so the specific
# spinal region is recoverable, which the published table could not do.
obs["grp"] = obs.dataset_id + "|" + obs.donor_id + "|" + obs.cell_type + "|" + obs.tissue
codes, uniq = pd.factorize(obs.grp, sort=False)
jid2g = pd.Series(codes, index=obs.soma_joinid.values)
gidx = {j: i for i, j in enumerate(sorted(jid2sym))}
S = np.zeros((len(uniq), len(vm))); P = np.zeros((len(uniq), len(vm)), dtype=np.int64)
n_cells = np.bincount(codes, minlength=len(uniq)).astype(np.int64)

q = h.axis_query(measurement_name="RNA",
                 obs_query=soma.AxisQuery(coords=(np.sort(obs.soma_joinid.values),)),
                 var_query=soma.AxisQuery(coords=(sorted(jid2sym),)))
for tbl in q.X("normalized").tables():
    d0 = tbl.column("soma_dim_0").to_numpy(); d1 = tbl.column("soma_dim_1").to_numpy()
    v  = tbl.column("soma_data").to_numpy()
    np.add.at(S, (jid2g.reindex(d0).to_numpy(), np.array([gidx[x] for x in d1])), v)
    np.add.at(P, (jid2g.reindex(d0).to_numpy(), np.array([gidx[x] for x in d1])), 1)
q.close(); census.close()

sym = [jid2sym[j] for j in sorted(jid2sym)]
D = pd.DataFrame({"grp": uniq, "n_cells": n_cells})
for i, g in enumerate(sym):
    D[f"sum_{g}"] = S[:, i]; D[f"pos_{g}"] = P[:, i]
D = D.merge(obs.drop_duplicates("grp")[["grp","dataset_id","donor_id","cell_type",
                                        "tissue","age","assay","suspension_type"]], on="grp")
D["ct"] = D.cell_type.map(CT_MAP)
for g in sym: D[f"m_{g}"] = D[f"sum_{g}"] / D.n_cells
D["zbtb48_per_1000_ACTB"] = 1000 * D.m_ZBTB48 / D.m_ACTB
D.to_csv("fig5c_reproduced_per_donor_celltype.csv", index=False)

F = D[D.ct.notna()]
print(f"\nFig 5c configuration: {int(F.n_cells.sum()):,} nuclei | {F.donor_id.nunique()} donors "
      f"| {F.dataset_id.nunique()} datasets   (published legend: 42,186 / 21)")
out = (F.groupby("ct")
        .apply(lambda d: pd.Series({
            "nuclei": int(d.n_cells.sum()), "donors": d.donor_id.nunique(),
            "datasets": d.dataset_id.nunique(),
            "per_1000_ACTB": 1000*d.sum_ZBTB48.sum()/d.sum_ACTB.sum(),
            "detection": d.pos_ZBTB48.sum()/d.n_cells.sum()}), include_groups=False)
        .sort_values("per_1000_ACTB", ascending=False))
out.to_csv("fig5c_reproduced_celltype_summary.csv")
print(out.round(4).to_string())

# ===========================================================================
# The published panel: donor-paired within one dataset
# ===========================================================================
from scipy import stats as st
from statsmodels.stats.multitest import multipletests

# The only one of the three datasets with genuine donor replication.
# Seeker et al., Acta Neuropathol Commun 2023, doi:10.1186/s40478-023-01568-z.
# Cervical spinal cord WHITE matter, which is a real limitation: the mouse
# result implicates dorsal horn grey matter, so this panel establishes
# astrocyte enrichment in human spinal cord but not in the specific
# compartment the mouse experiments target.
REPLICATED_DATASET = "c05e6940-729c-47bd-a2a6-6ce3730c4919"

MIN_NUCLEI_PER_GROUP = 20   # a donor x cell-type group below this is noise
MIN_DONORS_PER_POP = 8      # a population below this cannot support a paired test

sub = D[(D.dataset_id == REPLICATED_DATASET) & D.ct.notna() & (D.n_cells >= MIN_NUCLEI_PER_GROUP)]
donors_per_pop = sub.groupby("ct").donor_id.nunique().sort_values(ascending=False)
keep = list(donors_per_pop[donors_per_pop >= MIN_DONORS_PER_POP].index)
dropped = [c for c in donors_per_pop.index if c not in keep]

W = sub[sub.ct.isin(keep)].pivot_table(
    index="donor_id", columns="ct", values="zbtb48_per_1000_ACTB")
W.to_csv("fig5c_donor_by_celltype_matrix.csv")

print(f"\n=== published panel: {W.shape[0]} donors x {W.shape[1]} populations ===")
print(f"populations retained: {keep}")
print(f"populations dropped (single-donor or underpowered): {dropped}")

# Paired comparisons against astrocytes, Holm-corrected across the five tests.
rows = []
for c in [x for x in keep if x != "Astrocyte"]:
    p = W[["Astrocyte", c]].dropna()
    stat, pv = st.wilcoxon(p.Astrocyte, p[c])
    rows.append(dict(comparison=f"Astrocyte vs {c}", n_paired_donors=len(p),
                     median_astro=p.Astrocyte.median(), median_other=p[c].median(),
                     median_ratio=p.Astrocyte.median() / p[c].median(), W=stat, p_raw=pv))
S = pd.DataFrame(rows)
S["p_holm"] = multipletests(S.p_raw, method="holm")[1]
S["significant_holm_0.05"] = S.p_holm < 0.05

# Omnibus across the four best-covered populations, on complete cases only.
CORE = ["Astrocyte", "Microglia", "Oligodendrocyte", "Inhibitory neuron"]
Wc = W[[c for c in CORE if c in W.columns]].dropna()
fr = st.friedmanchisquare(*[Wc[c].values for c in Wc.columns])

# Per-donor ranking. The null is 1/k, not 1/2: the question is how often
# astrocytes rank FIRST among k populations, not whether they beat one other.
rk = W[keep].rank(axis=1, ascending=False)
top = int((rk["Astrocyte"] == 1).sum())
tot = int(W.Astrocyte.notna().sum())
binom = st.binomtest(top, tot, 1 / len(keep), alternative="greater")

S.to_csv("fig5c_paired_statistics.csv", index=False)
print("\nwithin-donor paired tests (Wilcoxon signed-rank, Holm-corrected):")
print(S.round(5).to_string(index=False))
print(f"\nFriedman, {Wc.shape[1]} populations x {Wc.shape[0]} complete donors: "
      f"chi2 = {fr.statistic:.2f}, P = {fr.pvalue:.3e}")
print(f"astrocytes ranked first of {len(keep)} in {top} of {tot} donors; "
      f"exact binomial vs 1/{len(keep)}: P = {binom.pvalue:.3e}")

assert "Astrocyte" in keep and len(keep) == 6, keep
assert S.p_holm.max() < 0.05, (
    "astrocytes must exceed every other population after Holm correction")
assert fr.pvalue < 1e-4, fr.pvalue
