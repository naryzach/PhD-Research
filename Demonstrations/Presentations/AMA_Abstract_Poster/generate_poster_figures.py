"""
generate_poster_figures.py
Generates publication-quality figures for the De Novo Binder Generation poster,
presentation, and paper. Focuses on MMP9 vs MMP2 specificity validated by
one-way ANOVA + Tukey-HSD (Binding Efficiency metric).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import matplotlib.ticker as ticker
import os, sys

# ── shared style ─────────────────────────────────────────────────────────────
DARK_BG   = "#0f172a"
PANEL_BG  = "#1e2d4a"
ACCENT    = "#38bdf8"
GREEN     = "#34d399"
AMBER     = "#fbbf24"
RED       = "#f87171"
VIOLET    = "#a78bfa"
TEXT_PRI  = "#f8fafc"
TEXT_SEC  = "#94a3b8"
GRID_CLR  = "#1e3a5f"

plt.rcParams.update({
    "figure.facecolor":  DARK_BG,
    "axes.facecolor":    PANEL_BG,
    "axes.edgecolor":    "#334155",
    "axes.labelcolor":   TEXT_PRI,
    "xtick.color":       TEXT_SEC,
    "ytick.color":       TEXT_SEC,
    "text.color":        TEXT_PRI,
    "grid.color":        GRID_CLR,
    "grid.linewidth":    0.6,
    "font.family":       "DejaVu Sans",
    "axes.titlesize":    13,
    "axes.labelsize":    11,
    "xtick.labelsize":   9,
    "ytick.labelsize":   9,
    "legend.fontsize":   9,
    "legend.framealpha": 0.3,
    "legend.edgecolor":  "#334155",
})

# Output directory: SharedAssets is the canonical location for all referenced figures.
# Original script location: Demonstrations/Presentations/AMA_Abstract_Poster/
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.normpath(os.path.join(
    _SCRIPT_DIR, "..", "..", "SharedAssets", "figures", "De_Novo_Binder_Generation"
))
os.makedirs(OUTDIR, exist_ok=True)

# ── Data (from selectivity_summary.csv / ANOVA_Results) ──────────────────────
# Binding Efficiency (DP/FITC+), mean ± SEM, per construct per target
# Only MMP2 and MMP9 data for the specificity story; ADAM17 as secondary
# Format: {construct: {target: (mean, sem, count)}}

SELDATA = {
    # ─── MMP9>MMP2 designed / significant ──────────────────────────────────
    "C 12":   {"MMP2": (0.321, 0.129, 2), "MMP9": (0.913, 0.003, 2), "ADAM17": (0.306, 0.0,  1)},
    "C 15":   {"MMP2": (0.338, 0.074, 2), "MMP9": (0.886, 0.027, 2), "ADAM17": (0.238, 0.0,  1)},
    "AB 6":   {"MMP2": (0.238, 0.011, 2), "MMP9": (0.899, 0.006, 2), "ADAM17": (0.189, 0.0,  1)},
    # ─── non-specific / low (designed low, ANOVA non-sig) ──────────────────
    "C 13":   {"MMP2": (0.160, 0.009, 3), "MMP9": (0.454, 0.403, 3), "ADAM17": (0.237, 0.025, 5)},
    "AB 5":   {"MMP2": (0.118, 0.033, 2), "MMP9": (0.376, 0.440, 3), "ADAM17": (0.216, 0.027, 5)},
    # ─── TIMP3 WT reference ─────────────────────────────────────────────────
    "TIMP 3": {"MMP2": (0.219, 0.078, 5), "MMP9": (0.545, 0.132, 7), "ADAM17": (0.366, 0.092, 7)},
}

ANOVA_P = {
    "C 12":   0.0115,
    "C 15":   0.0177,
    "AB 6":   0.00024,
    "C 13":   0.0961,
    "AB 5":   0.6006,
    "TIMP 3": 0.0240,
}

# Tukey significant pairs (group1, group2) for each construct
TUKEY_SIG = {
    "C 12":   [("MMP2","MMP9")],
    "C 15":   [("MMP2","MMP9"), ("TIMP3","MMP9")],   # MMP3 not shown
    "AB 6":   [("MMP2","MMP9")],
    "C 13":   [],
    "AB 5":   [],
    "TIMP 3": [],
}

def sig_label(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "n.s."


# ── Figure 1: MMP9 vs MMP2 Specificity Bar Chart ────────────────────────────
def fig_mmp9_vs_mmp2():
    constructs = ["C 12", "C 15", "AB 6", "C 13", "AB 5", "TIMP 3"]
    targets    = ["MMP2", "MMP9"]
    colors     = {"MMP2": RED, "MMP9": GREEN}
    bar_width  = 0.32
    x          = np.arange(len(constructs))

    fig, ax = plt.subplots(figsize=(11, 6))
    fig.patch.set_facecolor(DARK_BG)
    ax.set_facecolor(PANEL_BG)

    for i, (t, c) in enumerate([(t, colors[t]) for t in targets]):
        means = [SELDATA[cn].get(t, (0,0,0))[0] for cn in constructs]
        sems  = [SELDATA[cn].get(t, (0,0,0))[1] for cn in constructs]
        offset = (i - 0.5) * bar_width
        rects = ax.bar(x + offset, means, bar_width,
                       color=c, alpha=0.85, label=t,
                       error_kw=dict(ecolor="#f8fafc", capsize=4, lw=1.2))
        ax.errorbar(x + offset, means, yerr=sems,
                    fmt="none", ecolor="#f8fafc", capsize=4, lw=1.2)

    # Significance brackets
    bracket_targets = ["C 12", "C 15", "AB 6"]
    for cn in bracket_targets:
        xi = constructs.index(cn)
        m2 = SELDATA[cn]["MMP2"][0]
        m9 = SELDATA[cn]["MMP9"][0]
        y_top = max(m2, m9) + 0.07
        p     = ANOVA_P[cn]
        lbl   = sig_label(p)
        ax.annotate("", xy=(xi + bar_width/2, y_top + 0.025),
                    xytext=(xi - bar_width/2, y_top + 0.025),
                    arrowprops=dict(arrowstyle="-", color=ACCENT, lw=1.5))
        ax.text(xi, y_top + 0.04, lbl, ha="center", va="bottom",
                color=ACCENT, fontsize=11, fontweight="bold")

    # Non-sig labels
    for cn in ["C 13", "AB 5"]:
        xi = constructs.index(cn)
        m9 = SELDATA[cn]["MMP9"][0]
        ax.text(xi, m9 + 0.08, "n.s.", ha="center", va="bottom",
                color=TEXT_SEC, fontsize=9, style="italic")

    # TIMP3 dashed reference lines
    timp3_m2 = SELDATA["TIMP 3"]["MMP2"][0]
    timp3_m9 = SELDATA["TIMP 3"]["MMP9"][0]
    ax.axhline(timp3_m2, color=RED,   linestyle="--", lw=1.0, alpha=0.5)
    ax.axhline(timp3_m9, color=GREEN, linestyle="--", lw=1.0, alpha=0.5)

    # Vertical separator before TIMP 3
    ax.axvline(len(constructs) - 1 - 0.55, color="#334155", lw=1.0, ls="--")

    # Region labels
    ax.text(1.0, 1.02, "Designed MMP9-Selective", ha="center", va="bottom",
            transform=ax.get_xaxis_transform(), fontsize=8,
            color=ACCENT, style="italic")
    ax.text(3.5, 1.02, "Low / Non-Selective Controls", ha="center", va="bottom",
            transform=ax.get_xaxis_transform(), fontsize=8,
            color=TEXT_SEC, style="italic")

    ax.set_xticks(x)
    ax.set_xticklabels(constructs, fontsize=10)
    ax.set_ylabel("Binding Efficiency (DP / FITC⁺)", fontsize=11)
    ax.set_ylim(0, 1.18)
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
    ax.grid(axis="y", alpha=0.4)
    ax.set_title("MMP9 vs MMP2 Binding Efficiency — Specificity Validation",
                 fontsize=13, color=TEXT_PRI, pad=12)

    legend_handles = [
        mpatches.Patch(color=RED,   label="MMP2"),
        mpatches.Patch(color=GREEN, label="MMP9"),
    ]
    ax.legend(handles=legend_handles, loc="upper right")

    plt.tight_layout()
    out = os.path.join(OUTDIR, "fig_mmp9_vs_mmp2.png")
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor=DARK_BG)
    plt.close(fig)
    print(f"  ✓  Saved: {out}")


# ── Figure 2: Specificity Score (MMP9/MMP2 ratio) bar chart ─────────────────
def fig_specificity_ratio():
    """MMP9 binding efficiency / MMP2 binding efficiency ratio per construct."""
    constructs = ["C 12", "C 15", "AB 6", "C 13", "AB 5", "TIMP 3"]
    ratios = []
    for cn in constructs:
        m9 = SELDATA[cn].get("MMP9", (0,0,0))[0]
        m2 = SELDATA[cn].get("MMP2", (1e-9,0,0))[0]
        ratios.append(m9 / max(m2, 1e-9))

    fig, ax = plt.subplots(figsize=(10, 5.5))
    fig.patch.set_facecolor(DARK_BG)
    ax.set_facecolor(PANEL_BG)

    bar_colors = []
    for cn, r in zip(constructs, ratios):
        p = ANOVA_P.get(cn, 1.0)
        if p < 0.05 and r > 2:
            bar_colors.append(GREEN)
        elif cn == "TIMP 3":
            bar_colors.append(ACCENT)
        else:
            bar_colors.append(TEXT_SEC)

    bars = ax.bar(constructs, ratios, color=bar_colors, alpha=0.85, width=0.55)

    # ANOVA annotation
    for i, (cn, r) in enumerate(zip(constructs, ratios)):
        p = ANOVA_P.get(cn, 1.0)
        lbl = sig_label(p)
        color = ACCENT if p < 0.05 else TEXT_SEC
        ax.text(i, r + 0.1, lbl, ha="center", va="bottom",
                fontsize=10, fontweight="bold", color=color)

    ax.axhline(1.0, color=ACCENT, linestyle="--", lw=1.2, alpha=0.6,
               label="No Preference (ratio = 1)")
    ax.set_ylabel("MMP9 / MMP2 Binding Efficiency Ratio", fontsize=11)
    ax.set_title("MMP9 Selectivity Index vs MMP2\n(Binding Efficiency Ratio, ANOVA annotated)",
                 fontsize=12, color=TEXT_PRI)
    ax.grid(axis="y", alpha=0.4)
    ax.legend(loc="upper right")
    plt.tight_layout()
    out = os.path.join(OUTDIR, "fig_specificity_ratio.png")
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor=DARK_BG)
    plt.close(fig)
    print(f"  ✓  Saved: {out}")


# ── Figure 3: Binding Efficiency Heatmap (MMP2, MMP9, ADAM17) ───────────────
def fig_binding_heatmap():
    constructs = ["AB 1","AB 2","AB 3","AB 4","AB 5","AB 6","AB 7",
                  "C 11","C 12","C 13","C 14","C 15","TIMP 3"]

    # Full data from selectivity_summary.csv — binding efficiency means
    full_data = {
        "AB 1":  {"MMP2": 0.300, "MMP9": 0.555, "ADAM17": 0.401},
        "AB 2":  {"MMP2": 0.193, "MMP9": 0.617, "ADAM17": 0.393},
        "AB 3":  {"MMP2": 0.115, "MMP9": 0.554, "ADAM17": 0.235},
        "AB 4":  {"MMP2": 0.097, "MMP9": 0.347, "ADAM17": 0.155},
        "AB 5":  {"MMP2": 0.118, "MMP9": 0.376, "ADAM17": 0.216},
        "AB 6":  {"MMP2": 0.238, "MMP9": 0.899, "ADAM17": 0.189},
        "AB 7":  {"MMP2": 0.119, "MMP9": 0.648, "ADAM17": 0.243},
        "C 11":  {"MMP2": 0.112, "MMP9": 0.440, "ADAM17": 0.190},
        "C 12":  {"MMP2": 0.321, "MMP9": 0.913, "ADAM17": 0.306},
        "C 13":  {"MMP2": 0.160, "MMP9": 0.454, "ADAM17": 0.237},
        "C 14":  {"MMP2": 0.113, "MMP9": 0.454, "ADAM17": 0.324},
        "C 15":  {"MMP2": 0.338, "MMP9": 0.886, "ADAM17": 0.238},
        "TIMP 3":{"MMP2": 0.219, "MMP9": 0.545, "ADAM17": 0.366},
    }

    targets = ["MMP2", "MMP9", "ADAM17"]
    mat = np.array([[full_data[cn].get(t, np.nan) for t in targets]
                    for cn in constructs])

    fig, ax = plt.subplots(figsize=(11, 13))
    fig.patch.set_facecolor(DARK_BG)
    ax.set_facecolor(DARK_BG)

    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        "custom", ["#0f172a", "#1e3a5f", "#38bdf8", "#34d399", "#fbbf24"], N=256)

    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1.0)

    # Annotate cells
    for i in range(len(constructs)):
        for j in range(len(targets)):
            v = mat[i, j]
            if not np.isnan(v):
                txt_color = "black" if v > 0.6 else TEXT_PRI
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=9, color=txt_color, fontweight="bold")

    # Star annotations for significant constructs
    sig_rows = {"C 12": "*", "C 15": "*", "AB 6": "***"}
    for cn, star in sig_rows.items():
        i = constructs.index(cn)
        ax.text(len(targets) - 0.35, i, star, ha="left", va="center",
                fontsize=11, color=ACCENT, fontweight="bold")

    ax.set_xticks(range(len(targets)))
    ax.set_xticklabels(targets, fontsize=11, color=TEXT_PRI)
    ax.set_yticks(range(len(constructs)))
    ax.set_yticklabels(constructs, fontsize=10, color=TEXT_PRI)

    # Horizontal line before TIMP 3
    ax.axhline(len(constructs) - 1.5, color="#334155", lw=1.5, ls="--")

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label("Binding Efficiency (DP/FITC⁺)", color=TEXT_PRI, fontsize=9)
    cbar.ax.yaxis.set_tick_params(color=TEXT_PRI)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=TEXT_PRI)

    ax.set_title("Binding Efficiency Heatmap\n(* = ANOVA p<0.05 across targets)",
                 fontsize=12, color=TEXT_PRI, pad=10)
    plt.tight_layout()
    out = os.path.join(OUTDIR, "fig_binding_heatmap.png")
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor=DARK_BG)
    plt.close(fig)
    print(f"  ✓  Saved: {out}")


# ── Figure 4: Pipeline overview diagram (for paper/slides) ──────────────────
def fig_pipeline():
    """Horizontal pipeline diagram with 6 stages."""
    stages = [
        ("RFdiffusion\n(Backbone Gen.)", "AB / C / EF loops\n6–24 aa expansion", "#3b82f6"),
        ("ProteinMPNN\n(Sequence Design)", "1000 seqs / target\nT=0.2, fixed scaffold", "#8b5cf6"),
        ("AlphaFold Server\n(Co-folding)", "Batch JSON, 5 targets\nipTM / PAE / pLDDT", "#06b6d4"),
        ("Best Binder\nSelection", "T-score > 2.0\nMulti-metric consensus", "#10b981"),
        ("Twist Bioscience\n(DNA Synthesis)", "13 constructs\nCodon-opt., Dec 2025", "#f59e0b"),
        ("Yeast Display\n& Flow Cytometry", "Binding Efficiency\nANOVA + Tukey-HSD", "#ef4444"),
    ]

    fig, ax = plt.subplots(figsize=(14, 3.5))
    fig.patch.set_facecolor(DARK_BG)
    ax.set_facecolor(DARK_BG)
    ax.axis("off")

    box_w, box_h = 1.9, 1.6
    gap = 0.45
    total = len(stages) * box_w + (len(stages) - 1) * gap
    x0 = (14 - total) / 2

    for i, (title, subtitle, color) in enumerate(stages):
        x = x0 + i * (box_w + gap)
        rect = FancyBboxPatch((x, 0.6), box_w, box_h,
                              boxstyle="round,pad=0.08",
                              facecolor=color + "22",
                              edgecolor=color, linewidth=2,
                              transform=ax.transData, clip_on=False)
        ax.add_patch(rect)
        ax.text(x + box_w/2, 0.6 + box_h - 0.25, title,
                ha="center", va="top", fontsize=9, fontweight="bold",
                color=color, transform=ax.transData)
        ax.text(x + box_w/2, 0.6 + 0.25, subtitle,
                ha="center", va="bottom", fontsize=7.5, color=TEXT_SEC,
                transform=ax.transData)
        if i < len(stages) - 1:
            ax.annotate("", xy=(x + box_w + gap, 0.6 + box_h/2),
                        xytext=(x + box_w, 0.6 + box_h/2),
                        arrowprops=dict(arrowstyle="->", color=TEXT_SEC,
                                        lw=1.5, mutation_scale=14))

    ax.set_xlim(0, 14)
    ax.set_ylim(0, 2.8)
    ax.set_title("Computational-to-Experimental Design Pipeline",
                 fontsize=13, color=TEXT_PRI, pad=8)
    plt.tight_layout()
    out = os.path.join(OUTDIR, "fig_pipeline.png")
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor=DARK_BG)
    plt.close(fig)
    print(f"  ✓  Saved: {out}")


# ── main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Generating poster / paper figures …")
    fig_mmp9_vs_mmp2()
    fig_specificity_ratio()
    fig_binding_heatmap()
    fig_pipeline()
    print("Done.")
