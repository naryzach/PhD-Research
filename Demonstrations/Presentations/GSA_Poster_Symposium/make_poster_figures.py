"""Poster-scale versions of the figures (larger type, no baked-in footers).

Data sources:
  Figure 2 (matched replicates) and Figure 3 (fold by supplier):
      Local/Aggregate_FCS_Analysis/aggregate_summary.csv (QC-passing trials only)
  Figure 5 (calibration), computed here, not typed:
      Local/Calibration/recipe_scores.csv (AF3 ipTM, loop pLDDT per construct x target),
      Demonstrations/SharedAssets/data/De_Novo_Binder_Generation/experimental_binding.json
      (WT-normalized Pos Med Ratio, pooled across suppliers), and the aggregate above (Expr+ %).
  Figure 4 reuses the four 3Dmol.js panel renders in figures/src/.
"""
import os
import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(HERE, "figures")
plt.rcParams.update({"font.family": "DejaVu Sans", "figure.facecolor": "white", "savefig.facecolor": "white"})

# ---------------- Figure 2: matched (Enzo-only) replicates ----------------
agg = pd.read_csv(os.path.join(ROOT, "Local", "Aggregate_FCS_Analysis", "aggregate_summary.csv"))
agg = agg[(agg["Trial Failed"] == False) & (agg["Low Expression"] == False) & (agg["Low Events"] == False)]
agg = agg.dropna(subset=["Binding Efficiency"])
enzo = agg[agg["Source"] == "Enzo"]
M = "Pos Med Ratio"
MATCHED = ["AB 1", "AB 2", "AB 6", "C 12", "C 15", "TIMP 3"]
LABEL = {"TIMP 3": "TIMP3-WT"}
c2, c9 = "#e15759", "#59a14f"


def vals(c, t):
    s = enzo[(enzo["Construct"] == c) & (enzo["Target"] == t)]
    return s[M].astype(float).values, s["Date"].astype(int).values


fig, ax = plt.subplots(figsize=(10.5, 3.9))
x = np.arange(len(MATCHED))
w = 0.38
rng = np.random.default_rng(0)
tops = []
ax.axvspan(len(MATCHED) - 1.5, len(MATCHED) - 0.5, color="#eef2f7", zorder=0)
for i, c in enumerate(MATCHED):
    for t, off, col in (("MMP2", -w / 2, c2), ("MMP9", w / 2, c9)):
        v, d = vals(c, t)
        ax.bar(i + off, v.mean(), w, color=col, alpha=0.4, edgecolor="white")
        for y, dd, j in zip(v, d, rng.uniform(-0.07, 0.07, len(v))):
            ax.scatter(i + off + j, y, s=70, color=col, edgecolor="black", linewidth=0.8,
                       marker="o" if dd == 20260424 else "s", zorder=4)
    a, _ = vals(c, "MMP9")
    b, _ = vals(c, "MMP2")
    p = stats.ttest_ind(a, b, equal_var=False)[1]
    top = a.max() + 0.035
    tops.append(top)
    ax.plot([i - w / 2, i + w / 2], [top, top], color="#1f77b4", lw=1.8)
    ax.text(i, top + 0.006, f"p={p:.3f}", ha="center", va="bottom", color="#1f77b4", fontsize=13, fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.85), zorder=6)
ax.set_xticks(x)
ax.set_xticklabels([f"{LABEL.get(c, c)}\nn={len(vals(c, 'MMP2')[0])}/{len(vals(c, 'MMP9')[0])}" for c in MATCHED], fontsize=14)
ax.set_ylabel("Pos Med Ratio (APC/FITC)", fontsize=15)
ax.tick_params(axis="y", labelsize=13)
ax.set_ylim(0, max(tops) + 0.07)
ax.grid(axis="y", alpha=0.3)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
leg = [mpatches.Patch(color=c2, alpha=0.5, label="MMP2"), mpatches.Patch(color=c9, alpha=0.5, label="MMP9"),
       plt.Line2D([], [], marker="o", ls="", color="gray", markeredgecolor="black", label="2026-04-24"),
       plt.Line2D([], [], marker="s", ls="", color="gray", markeredgecolor="black", label="2026-05-09")]
ax.legend(handles=leg, loc="lower left", bbox_to_anchor=(0.0, 1.0), fontsize=13, ncol=4, frameon=False)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "poster_fig2_selectivity.png"), dpi=300)
plt.close()

# ---------------- Figure 3: MMP9:MMP2 fold by supplier ----------------
ORDER = ["AB 1", "AB 2", "AB 6", "C 12", "C 15", "TIMP 3", "AB 3", "AB 4", "AB 7", "C 11", "C 13", "C 14", "ABC 22"]


def fold(c, vendor):
    s = agg[(agg["Construct"] == c) & (agg["Source"] == vendor)]
    a, b = s[s["Target"] == "MMP9"][M].values, s[s["Target"] == "MMP2"][M].values
    return (a.mean() / b.mean() if len(a) and len(b) else np.nan), len(a), len(b)


fig, ax = plt.subplots(figsize=(10.5, 3.5))
xs = np.arange(len(ORDER))
ax.axhline(1, color="black", lw=1)
ax.axvspan(-0.5, 4.5, color="#e8f3e8", zorder=0)
ax.axvspan(4.5, 5.5, color="#eef2f7", zorder=0)
for i, c in enumerate(ORDER):
    fe, n9, n2 = fold(c, "Enzo")
    fs, s9, s2 = fold(c, "Sino")
    ax.scatter(i - 0.1, fe, s=110, color="#2b7a2b", edgecolor="black", zorder=4,
               marker="o" if min(n9, n2) >= 2 else "^")
    ax.scatter(i + 0.1, fs, s=90, facecolor="white", edgecolor="#d97706", linewidth=2, zorder=4, marker="D")
ax.set_yscale("log")
ax.set_ylim(0.9, 4.9)
ax.set_yticks([1, 2, 3, 4])
ax.set_yticklabels(["1", "2", "3", "4"], fontsize=13)
ax.minorticks_off()
ax.set_xticks(xs)
ax.set_xticklabels(["WT" if c == "TIMP 3" else c for c in ORDER], fontsize=12)
ax.set_ylabel("MMP9 : MMP2 fold\n(Pos Med Ratio)", fontsize=14)
ax.text(2, 4.35, "designed for MMP9", ha="center", fontsize=13, color="#2b7a2b", fontweight="bold")
ax.text(5, 4.35, "WT", ha="center", fontsize=13, color="#334", fontweight="bold")
ax.text(9.0, 4.35, "other constructs", ha="center", fontsize=13, color="#555", fontweight="bold")
ax.grid(axis="y", alpha=0.3)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
leg = [plt.Line2D([], [], marker="o", ls="", color="#2b7a2b", markeredgecolor="black", markersize=10, label="Enzo, matched (n >= 2 per target)"),
       plt.Line2D([], [], marker="^", ls="", color="#2b7a2b", markeredgecolor="black", markersize=10, label="Enzo, matched (n = 1 per target)"),
       plt.Line2D([], [], marker="D", ls="", markerfacecolor="white", markeredgecolor="#d97706", markeredgewidth=2, markersize=9, label="Sino, not matched")]
ax.legend(handles=leg, loc="lower left", bbox_to_anchor=(0.0, 1.0), fontsize=12, ncol=3, frameon=False,
          columnspacing=1.0, handletextpad=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "poster_fig3_fold_by_supplier.png"), dpi=300)
plt.close()

# ---------------- Figure 5: calibration (computed, not typed) ----------------
rec = pd.read_csv(os.path.join(ROOT, "Local", "Calibration", "recipe_scores.csv"))
with open(os.path.join(ROOT, "Demonstrations", "SharedAssets", "data", "De_Novo_Binder_Generation",
                       "experimental_binding.json")) as f:
    _exp = json.load(f)
_nmr = {(r["Construct"], t): r["Norm Median Ratio"] for t, rows in _exp.items() for r in rows}
rec["NMR"] = [_nmr[(c, t)] for c, t in zip(rec["Construct"], rec["Target"])]
rec = rec.merge(agg.groupby(["Target", "Construct"])["Expr+ %"].mean().reset_index(), on=["Target", "Construct"], how="left")
_g = rec["NMR"].mean()
rec["resid"] = (rec["NMR"] - rec.groupby("Construct")["NMR"].transform("mean")
                - rec.groupby("Target")["NMR"].transform("mean") + _g)
_ssc = (rec.groupby("Construct")["NMR"].mean() - _g).pow(2).sum() * rec["Target"].nunique()
_sst = (rec.groupby("Target")["NMR"].mean() - _g).pow(2).sum() * rec["Construct"].nunique()
_ssi = rec["resid"].pow(2).sum()
_tot = _ssc + _sst + _ssi
vals_a = [round(100 * _ssc / _tot), round(100 * _sst / _tot), round(100 * _ssi / _tot)]
d_expr = [round(stats.spearmanr(rec[m], rec["Expr+ %"])[0], 2) for m in ("LpLDDT", "ipTM")]
d_bind = [round(stats.spearmanr(rec[m], rec["resid"])[0], 2) for m in ("LpLDDT", "ipTM")]
print("calibration (computed):", vals_a, d_expr, d_bind, "constructs", rec["Construct"].nunique(), "pairs", len(rec))
assert vals_a == [64, 7, 29] and d_bind == [0.20, 0.09] and d_expr == [-0.06, 0.37], \
    "calibration values changed; update the p-values and text in poster_data.yaml"

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.6), gridspec_kw={"width_ratios": [1, 1.25]})
labels = ["Construct\n(avidity)", "Target\nbaseline", "Target-specific\ninteraction"]
a1.bar(labels, vals_a, color=["#9a9a9a", "#4e79a7", "#59a14f"])
for i, v in enumerate(vals_a):
    a1.text(i, v + 1.5, f"{v}%", ha="center", fontsize=17, fontweight="bold")
a1.set_ylabel("Share of binding variance (%)", fontsize=14)
a1.set_ylim(0, 78)
a1.tick_params(axis="x", labelsize=13)
a1.tick_params(axis="y", labelsize=12)
a1.set_title("Variance decomposition", fontsize=15)
xx = np.arange(2)
bw = 0.34
a2.bar(xx - bw / 2, d_expr, bw, color="#f28e2b", label="display level (Expr+ %)")
a2.bar(xx + bw / 2, d_bind, bw, color="#4e79a7", label="target-specific binding")
for i in range(2):
    a2.text(i - bw / 2, d_expr[i] + (0.02 if d_expr[i] >= 0 else -0.05), f"{d_expr[i]:.2f}", ha="center", fontsize=13)
    a2.text(i + bw / 2, d_bind[i] + 0.02, f"{d_bind[i]:.2f}", ha="center", fontsize=13)
a2.axhline(0, color="black", lw=1)
a2.set_xticks(xx)
a2.set_xticklabels(["loop pLDDT", "ipTM"], fontsize=14)
a2.set_ylabel("Spearman rho", fontsize=14)
a2.set_ylim(-0.15, 0.55)
a2.tick_params(axis="y", labelsize=12)
a2.set_title("AlphaFold3 metric vs.", fontsize=15)
a2.legend(fontsize=12, loc="upper left")
for sp in ("top", "right"):
    a1.spines[sp].set_visible(False)
    a2.spines[sp].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "poster_fig4_calibration.png"), dpi=300)
plt.close()

# ---------------- Figure 4: interface panels ----------------
src = os.path.join(OUT, "src")
panels = [("pocket_c12_mmp9.png", "C 12 + MMP9"), ("pocket_c12_mmp2.png", "C 12 + MMP2"),
          ("pocket_ab2_mmp9.png", "AB 2 + MMP9"), ("pocket_ab2_mmp2.png", "AB 2 + MMP2")]
fig = plt.figure(figsize=(8.6, 9.4))
gs = fig.add_gridspec(2, 2, hspace=0.10, wspace=0.02, left=0.01, right=0.99, top=0.95, bottom=0.11)
for i, (f, t) in enumerate(panels):
    ax = fig.add_subplot(gs[i // 2, i % 2])
    ax.imshow(Image.open(os.path.join(src, f)))
    ax.axis("off")
    ax.set_title(t, fontsize=22, fontweight="bold", pad=4)
leg = [mpatches.Patch(color="#2c6fbb", label="TIMP3 construct"),
       mpatches.Patch(color="#c0392b", label="within 5 Å of target"),
       mpatches.Patch(color="#d0d0d0", label="MMP9 or MMP2")]
fig.legend(handles=leg, loc="lower center", ncol=3, fontsize=17, frameon=False, bbox_to_anchor=(0.5, 0.0))
plt.savefig(os.path.join(OUT, "poster_fig3_interface.png"), dpi=300)
plt.close()
print("poster figures written")
