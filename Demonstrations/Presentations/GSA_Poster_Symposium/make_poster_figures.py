"""Poster-scale versions of the figures (larger type, no baked-in footers).

Data sources: Local/Aggregate_FCS_Analysis/aggregate_summary.csv (Figure 2) and the
calibration statistics reported in the paper (Figure 4). Figure 3 reuses the four 3Dmol.js
panel renders in figures/src/.
"""
import os
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

# ---------------- Figure 2: replicates ----------------
agg = pd.read_csv(os.path.join(ROOT, "Local", "Aggregate_FCS_Analysis", "aggregate_summary.csv"))
agg = agg[(agg["Trial Failed"] == False) & (agg["Low Expression"] == False) & (agg["Low Events"] == False)]
agg = agg.dropna(subset=["Binding Efficiency"])
enzo = agg[agg["Source"] == "Enzo"]
M = "Pos Med Ratio"
MATCHED = ["AB 1", "AB 2", "AB 6", "C 12", "C 15", "TIMP 3"]
order = ["AB 1", "AB 2", "AB 4", "AB 5", "AB 6", "AB 7", "C 11", "C 12", "C 13", "C 14", "C 15", "ABC 22", "TIMP 3"]


def vals(c, t):
    src = enzo if c in MATCHED else agg
    s = src[(src["Construct"] == c) & (src["Target"] == t)]
    return s[M].astype(float).values, s["Source"].astype(str).values


fig, ax = plt.subplots(figsize=(10.5, 5.1))
x = np.arange(len(order))
w = 0.38
c2, c9 = "#e15759", "#59a14f"
rng = np.random.default_rng(0)
tops = []
for i, c in enumerate(order):
    for t, off, col in (("MMP2", -w / 2, c2), ("MMP9", w / 2, c9)):
        v, src = vals(c, t)
        ax.bar(i + off, v.mean(), w, color=col, alpha=0.4, edgecolor="white")
        for y, s, j in zip(v, src, rng.uniform(-0.08, 0.08, len(v))):
            ax.scatter(i + off + j, y, s=60, color=col, edgecolor="black", linewidth=0.8,
                       marker="o" if s == "Enzo" else "D", zorder=4)
    if c in MATCHED:
        a, _ = vals(c, "MMP9"); b, _ = vals(c, "MMP2")
        p = stats.ttest_ind(a, b, equal_var=False)[1]
        top = a.max() + 0.035
        tops.append(top)
        ax.plot([i - w / 2, i + w / 2], [top, top], color="#1f77b4", lw=1.8)
        ax.text(i, top + 0.006, f"p={p:.3f}", ha="center", va="bottom", color="#1f77b4", fontsize=12, fontweight="bold", bbox=dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.85), zorder=6)
wt2 = vals("TIMP 3", "MMP2")[0].mean(); wt9 = vals("TIMP 3", "MMP9")[0].mean()
ax.axhline(wt2, ls="--", lw=1.3, color=c2, alpha=0.7)
ax.axhline(wt9, ls="--", lw=1.3, color=c9, alpha=0.7)
ax.set_xticks(x)
ax.set_xticklabels([f"{c}\nn={len(vals(c,'MMP2')[0])}/{len(vals(c,'MMP9')[0])}" for c in order], fontsize=12)
ax.set_ylabel("Pos Med Ratio (APC/FITC)", fontsize=15)
ax.tick_params(axis="y", labelsize=13)
ax.set_ylim(0, max(0.72, max(tops) + 0.06))
ax.grid(axis="y", alpha=0.3)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
leg = [mpatches.Patch(color=c2, alpha=0.5, label="MMP2"), mpatches.Patch(color=c9, alpha=0.5, label="MMP9"),
       plt.Line2D([], [], marker="o", ls="", color="gray", markeredgecolor="black", label="Enzo"),
       plt.Line2D([], [], marker="D", ls="", color="gray", markeredgecolor="black", label="other supplier")]
ax.legend(handles=leg, loc="lower left", bbox_to_anchor=(0.0, 1.0), fontsize=13, ncol=4, frameon=False)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "poster_fig2_selectivity.png"), dpi=300)
plt.close()

# ---------------- Figure 4: calibration ----------------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.6), gridspec_kw={"width_ratios": [1, 1.25]})
labels = ["Construct\n(avidity)", "Target\nbaseline", "Target-specific\ninteraction"]
vals_a = [64, 7, 29]
a1.bar(labels, vals_a, color=["#9a9a9a", "#4e79a7", "#59a14f"])
for i, v in enumerate(vals_a):
    a1.text(i, v + 1.5, f"{v}%", ha="center", fontsize=17, fontweight="bold")
a1.set_ylabel("Share of binding variance (%)", fontsize=14)
a1.set_ylim(0, 78)
a1.tick_params(axis="x", labelsize=13); a1.tick_params(axis="y", labelsize=12)
a1.set_title("Variance decomposition", fontsize=15)
xx = np.arange(2)
bw = 0.34
d_expr = [-0.06, 0.37]; d_bind = [0.20, 0.09]
a2.bar(xx - bw / 2, d_expr, bw, color="#f28e2b", label="display level (Expr+ %)")
a2.bar(xx + bw / 2, d_bind, bw, color="#4e79a7", label="target-specific binding")
for i in range(2):
    a2.text(i - bw / 2, d_expr[i] + (0.02 if d_expr[i] >= 0 else -0.05), f"{d_expr[i]:.2f}", ha="center", fontsize=13)
    a2.text(i + bw / 2, d_bind[i] + 0.02, f"{d_bind[i]:.2f}", ha="center", fontsize=13)
a2.axhline(0, color="black", lw=1)
a2.set_xticks(xx); a2.set_xticklabels(["loop pLDDT", "ipTM"], fontsize=14)
a2.set_ylabel("Spearman rho", fontsize=14)
a2.set_ylim(-0.15, 0.55)
a2.tick_params(axis="y", labelsize=12)
a2.set_title("AlphaFold3 metric vs.", fontsize=15)
a2.legend(fontsize=12, loc="upper left")
for sp in ("top", "right"):
    a1.spines[sp].set_visible(False); a2.spines[sp].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "poster_fig4_calibration.png"), dpi=300)
plt.close()

# ---------------- Figure 3: interface panels ----------------
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
