"""
Regenerate Figure 5c/d and every reported statistic from the deposited source data.
No network access required. Run from the folder containing the CSVs.

    python make_fig5c_panel.py
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as st
from statsmodels.stats.multitest import multipletests

W = pd.read_csv("fig5c_donor_by_celltype_matrix.csv", index_col=0)
ORDER = list(W.median().sort_values(ascending=False).index)
FOCAL, OTHER, GREY = "#B2182B", "#4D4D4D", "#6E6E6E"

# ---- statistics (recomputed, not read from file) -------------------------------
rows = []
for c in [x for x in ORDER if x != "Astrocyte"]:
    p = W[["Astrocyte", c]].dropna()
    s, pv = st.wilcoxon(p.Astrocyte, p[c])
    rows.append(dict(comparison=c, n=len(p), W=s, p_raw=pv,
                     median_ratio=p.Astrocyte.median()/p[c].median()))
S = pd.DataFrame(rows)
S["p_holm"] = multipletests(S.p_raw, method="holm")[1]
pmap = dict(zip(S.comparison, S.p_holm))

pa = W[["Astrocyte", "Microglia"]].dropna()
wst = st.wilcoxon(pa.Astrocyte, pa.Microglia)
CORE = ["Astrocyte", "Microglia", "Oligodendrocyte", "Inhibitory neuron"]
Wc = W[CORE].dropna()
fr = st.friedmanchisquare(*[Wc[c].values for c in CORE])
rk = W[ORDER].rank(axis=1, ascending=False)
top, tot = int((rk.Astrocyte == 1).sum()), int(W.Astrocyte.notna().sum())
sgn = st.binomtest(top, tot, 1/len(ORDER), alternative="greater")

print(f"Astrocyte vs Microglia: n={len(pa)} medians {pa.Astrocyte.median():.1f} vs "
      f"{pa.Microglia.median():.1f} ratio {pa.Astrocyte.median()/pa.Microglia.median():.2f} "
      f"W={wst.statistic:.0f} P={wst.pvalue:.4f}")
print(f"Friedman ({len(CORE)} types, {Wc.shape[0]} donors): chi2={fr.statistic:.2f} P={fr.pvalue:.3g}")
print(f"Astrocyte ranked 1st in {top}/{tot} donors, binomial vs 1/{len(ORDER)} P={sgn.pvalue:.3g}")
print(S.round(5).to_string(index=False))
S.to_csv("fig5c_statistics_regenerated.csv", index=False)

# ---- figure -------------------------------------------------------------------
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 6, "ytick.labelsize": 6,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 300, "savefig.dpi": 300})
SHORT = {"Inhibitory neuron": "Inhib. neuron", "Excitatory neuron": "Excit. neuron"}
fig = plt.figure(figsize=(7.2, 3.2))
gs = fig.add_gridspec(1, 2, width_ratios=[2.05, 1.0], wspace=0.40)
rng = np.random.default_rng(0)

ax = fig.add_subplot(gs[0])
for i, c in enumerate(ORDER):
    v = W[c].dropna().values
    col = FOCAL if c == "Astrocyte" else OTHER
    ax.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=11, color=col, alpha=0.55,
               linewidths=0, zorder=2)
    ax.plot([i-0.30, i+0.30], [np.median(v)]*2, color=col, lw=1.8, zorder=3,
            solid_capstyle="butt")
ax.set_ylim(bottom=0); ax.margins(x=0.06, y=0.16)
for i, c in enumerate(ORDER):
    if c == "Astrocyte": continue
    p = pmap[c]
    ax.annotate("P < 0.001" if p < 0.001 else f"P = {p:.3f}",
                (i, W[c].dropna().max()), xytext=(0, 5), textcoords="offset points",
                ha="center", fontsize=6, color=GREY)
ax.set_xticks(range(len(ORDER)))
ax.set_xticklabels([SHORT.get(c, c) for c in ORDER], rotation=28, ha="right",
                   rotation_mode="anchor")
ax.set_ylabel("$\\it{ZBTB48}$ per 1,000 $\\it{ACTB}$")
ax.set_title("Astrocytes carry the highest $\\it{ZBTB48}$ of any spinal cell type", loc="left")
ax.annotate(f"{int(W[ORDER].notna().sum().min())} to {int(W[ORDER].notna().sum().max())}"
            " donors per cell type\nWilcoxon vs astrocyte, Holm-corrected",
            (0.985, 0.965), xycoords="axes fraction", ha="right", va="top",
            fontsize=6, color=GREY, linespacing=1.4)
ax.text(-0.14, 1.06, "c", transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")

ax2 = fig.add_subplot(gs[1])
for _, r in pa.iterrows():
    ax2.plot([0, 1], [r.Astrocyte, r.Microglia], color=GREY, lw=0.6, alpha=0.55, zorder=1)
ax2.scatter([0]*len(pa), pa.Astrocyte, s=13, color=FOCAL, zorder=3, linewidths=0)
ax2.scatter([1]*len(pa), pa.Microglia, s=13, color=OTHER, zorder=3, linewidths=0)
for x, v, cl in [(0, pa.Astrocyte.median(), FOCAL), (1, pa.Microglia.median(), OTHER)]:
    ax2.plot([x-0.18, x+0.18], [v]*2, color=cl, lw=2.0, zorder=4, solid_capstyle="butt")
ax2.set_xticks([0, 1])
ax2.set_xticklabels(["Astrocyte", "Microglia"], rotation=28, ha="right", rotation_mode="anchor")
ax2.set_xlim(-0.45, 1.45); ax2.set_ylim(bottom=0); ax2.margins(y=0.18)
ax2.set_ylabel("$\\it{ZBTB48}$ per 1,000 $\\it{ACTB}$")
ax2.set_title(f"Paired within donor ($n$ = {len(pa)})", loc="left")
ax2.annotate(f"P = {wst.pvalue:.3f}", (0.5, pa.values.max()), xytext=(0, 7),
             textcoords="offset points", ha="center", fontsize=6, color=GREY)
ax2.text(-0.28, 1.06, "d", transform=ax2.transAxes, fontsize=11, fontweight="bold", va="top")

fig.savefig("fig5c_regenerated.png", dpi=300, bbox_inches="tight")
print("\nwrote fig5c_regenerated.png and fig5c_statistics_regenerated.csv")
