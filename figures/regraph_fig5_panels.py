"""Regenerate the three ZBTB48 phenome-specificity panels from the saved primary data.
Reads only the CSV files shipped alongside; no network access and no session state required.
    python regraph_fig5_panels.py
"""
import pandas as pd, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt

H  = pd.read_csv("panel_h_source_data.csv")
I  = pd.read_csv("panel_i_source_data.csv")
J  = pd.read_csv("panel_j_source_data.csv")
LB = pd.read_csv("figure_display_labels.csv").set_index("phenocode").label_as_plotted

FOCAL, OTHER, GREY = "#B2182B", "#4D4D4D", "#8C8C8C"
XL = "Odds ratio per s.d. higher genetically predicted $\\it{ZBTB48}$ (95% CI)"
mpl.rcParams.update({"font.size":8,"axes.titlesize":7,"axes.labelsize":7,
                     "xtick.labelsize":6,"ytick.labelsize":6,"axes.spines.top":False,
                     "axes.spines.right":False,"figure.dpi":110})

def ticks(ax, vals):
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(mpl.ticker.FixedLocator(vals))
    ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.xaxis.set_major_formatter(mpl.ticker.FixedFormatter([str(v) for v in vals]))

def forest(ax, df, tv, title, colour, marker):
    d = df.sort_values("beta")
    lab = [LB.get(c, s) for c, s in zip(d.phenocode, d.phenostring)]
    ax.errorbar(d.OR_per_sd, np.arange(len(d)),
                xerr=[d.OR_per_sd-d.lo_per_sd, d.hi_per_sd-d.OR_per_sd],
                fmt=marker, ms=4.2, color=colour, lw=0.9, capsize=1.6)
    ax.axvline(1, color=GREY, lw=0.8, ls="--")
    ax.set_yticks(range(len(d))); ax.set_yticklabels(lab, fontsize=6)
    ticks(ax, tv); ax.set_xlabel(XL); ax.set_title(title, loc="left")
    ax.margins(x=0.05, y=0.02)

fig = plt.figure(figsize=(7.4, 6.9))
gs  = fig.add_gridspec(3, 1, height_ratios=[0.95, 2.1, 0.95], hspace=0.50)

ax = fig.add_subplot(gs[0])
E  = H.iloc[::-1]
cols = [FOCAL if g == "Pain symptoms" else OTHER for g in E.Group]
ax.hlines(range(len(E)), 1, E.fold, color=cols, lw=1.2)
ax.scatter(E.fold, range(len(E)), s=46, color=cols, zorder=3)
ax.axvline(1, color=GREY, lw=0.8, ls="--")
ax.set_yticks(range(len(E))); ax.set_yticklabels(E.Group)
ticks(ax, [0.25, 0.5, 1, 2, 4])
ax.set_xlabel("Nominally associated endpoints, observed / matched expectation")
ax.set_title("Association at rs17029626 is specific to pain, not to nervous-system disease", loc="left")
for y, (f, pv) in enumerate(zip(E.fold, E.p)):
    txt = (f"{f:.1f}x, P=5x10$^{{-5}}$" if pv < 1e-4 else
           (f"{f:.1f}x, P>0.99" if pv > 0.99 else f"{f:.1f}x, P={pv:.2g}"))
    ha, dx = ("left", 7) if f >= 1 else ("right", -7)
    ax.annotate(txt, (f, y), xytext=(dx, 0), textcoords="offset points",
                va="center", ha=ha, fontsize=6, color=cols[y])
ax.set_xlim(0.16, 7.4); ax.margins(y=0.18)

forest(fig.add_subplot(gs[1]), I, [1, 1.5, 2, 3],
       f"Pain-symptom endpoints associated with ZBTB48: {len(I)} of 28 tested", FOCAL, "o")
forest(fig.add_subplot(gs[2]), J, [1, 2, 4, 8],
       f"Nervous-system endpoints without pain: {len(J)} of 131 tested", OTHER, "s")
for ax, l in zip(fig.axes, "hij"):
    ax.annotate(l, (0, 1), xycoords="axes fraction", xytext=(-46, 8),
                textcoords="offset points", fontweight="bold", fontsize=10)
fig.savefig("fig5_pain_specificity_REGRAPHED.png", dpi=300, bbox_inches="tight")
print("wrote fig5_pain_specificity_REGRAPHED.png")
