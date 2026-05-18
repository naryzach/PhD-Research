"""
generate_poster_figures.py
Generates publication-quality figures for the De Novo Binder Generation poster,
presentation, and paper. All figures are written to SharedAssets/figures/De_Novo_Binder_Generation/.

Primary poster metric: Binding Median of Expression-Positive Events (Bind Med (Expr+))
  — raw median APC-A fluorescence restricted to FITC-positive (expressing) cells.
  NOT normalized — raw values enable direct cross-target selectivity comparison:
  a higher bar for MMP9 than MMP2 on the same construct is direct evidence of preference.
  TIMP3-WT reference lines are drawn at the actual TIMP3-WT MFI per target.

Additional figures (used by paper and presentation):
  - Binding Efficiency (DP/FITC+) bar chart and heatmap
  - Norm Bind Med (Expr+) heatmap (within-target construct ranking, paper/slides)
  - MMP9/MMP2 selectivity ratio
  - 6-stage pipeline diagram

Data source priority:
  1. SharedAssets/data/De_Novo_Binder_Generation/aggregate_summary.csv
  2. Hardcoded fallback constants (Binding Efficiency only)

Run:
    python generate_poster_figures.py
    python generate_poster_figures.py --csv path/to/aggregate_summary.csv
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
import os, sys, argparse, json
import warnings

try:
    from scipy import stats as _scipy_stats
    def _ci95(vals):
        n = len(vals)
        if n < 2:
            return 0.0
        return float(_scipy_stats.t.ppf(0.975, n - 1) * np.std(vals, ddof=1) / np.sqrt(n))
except ImportError:
    def _ci95(vals):
        n = len(vals)
        return float(1.96 * np.std(vals, ddof=1) / np.sqrt(n)) if n > 1 else 0.0

warnings.filterwarnings("ignore")

# ── Shared style ──────────────────────────────────────────────────────────────
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

# ── Output directory ──────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.normpath(os.path.join(
    _SCRIPT_DIR, "..", "..", "SharedAssets", "figures", "De_Novo_Binder_Generation"
))
DATA_DIR = os.path.normpath(os.path.join(
    _SCRIPT_DIR, "..", "..", "SharedAssets", "data", "De_Novo_Binder_Generation"
))
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

# ── Construct order for plots ─────────────────────────────────────────────────
CONSTRUCT_ORDER = [
    "C 12", "C 15", "AB 6",          # MMP9-selective (ANOVA confirmed)
    "AB 1",                           # MMP9-predicted selective (ANOVA n.s. — high variance)
    "AB 2", "AB 3",                   # High affinity, non-significant
    "C 13", "AB 5",                   # Low binder controls
    "AB 4", "AB 7", "C 11", "C 14",  # ADAM17-targeted
    "TIMP 3",                          # WT reference
]

# ── Fallback hardcoded data (Binding Efficiency, DP/FITC+) ───────────────────
# Used when aggregate_summary.csv is not available.
# Format: {construct: {target: (mean, sem, n_trials)}}
_BIND_EFF_FALLBACK = {
    "C 12":   {"MMP2": (0.321, 0.129, 2), "MMP9": (0.913, 0.003, 2), "ADAM17": (0.306, 0.000, 1)},
    "C 15":   {"MMP2": (0.338, 0.074, 2), "MMP9": (0.886, 0.027, 2), "ADAM17": (0.238, 0.000, 1)},
    "AB 6":   {"MMP2": (0.238, 0.011, 2), "MMP9": (0.899, 0.006, 2), "ADAM17": (0.189, 0.000, 1)},
    "AB 1":   {"MMP2": (0.300, 0.000, 2), "MMP9": (0.555, 0.434, 4), "ADAM17": (0.401, 0.000, 1)},
    "AB 2":   {"MMP2": (0.193, 0.150, 2), "MMP9": (0.617, 0.420, 2), "ADAM17": (0.393, 0.000, 1)},
    "AB 3":   {"MMP2": (0.115, 0.120, 2), "MMP9": (0.554, 0.380, 2), "ADAM17": (0.235, 0.000, 1)},
    "C 13":   {"MMP2": (0.160, 0.018, 3), "MMP9": (0.454, 0.403, 3), "ADAM17": (0.237, 0.025, 5)},
    "AB 5":   {"MMP2": (0.118, 0.065, 2), "MMP9": (0.376, 0.440, 3), "ADAM17": (0.216, 0.027, 5)},
    "AB 4":   {"MMP2": (0.097, 0.082, 3), "MMP9": (0.347, 0.320, 3), "ADAM17": (0.155, 0.000, 1)},
    "AB 7":   {"MMP2": (0.119, 0.075, 3), "MMP9": (0.648, 0.410, 3), "ADAM17": (0.243, 0.000, 1)},
    "C 11":   {"MMP2": (0.112, 0.088, 3), "MMP9": (0.440, 0.360, 3), "ADAM17": (0.190, 0.000, 1)},
    "C 14":   {"MMP2": (0.113, 0.069, 3), "MMP9": (0.454, 0.390, 3), "ADAM17": (0.324, 0.000, 1)},
    "TIMP 3": {"MMP2": (0.219, 0.078, 5), "MMP9": (0.545, 0.132, 7), "ADAM17": (0.366, 0.092, 7)},
}

# One-way ANOVA p-values (Binding Efficiency, across targets per construct)
_ANOVA_P_BIND_EFF = {
    "C 12":   0.0115,  "C 15":  0.0177,  "AB 6":  0.00024,
    "AB 1":   0.3500,  # Predicted selective; n.s. due to high trial-to-trial variance
    "AB 2":   0.1200,  "AB 3":  0.1850,  "C 13":  0.0961,  "AB 5":  0.6006,
    "AB 4":   0.2500,  "AB 7":  0.3100,  "C 11":  0.4200,  "C 14":  0.5800,
    "TIMP 3": 0.0240,
}


def sig_label(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "n.s."


# ── Construct label colors — match poster table CSS row classes ───────────────
_CONSTRUCT_LABEL_CLR = {
    # sig-row: ANOVA-confirmed MMP9-selective (green)
    "C 12":   GREEN,
    "C 15":   GREEN,
    "AB 6":   GREEN,
    # default: predicted / non-selective / high-affinity (no special class)
    "AB 1":   TEXT_PRI,
    "AB 2":   TEXT_PRI,
    "AB 3":   TEXT_PRI,
    "ABC 22": TEXT_PRI,
    # ctrl-row: low-binder controls (muted)
    "C 13":   TEXT_SEC,
    "AB 5":   TEXT_SEC,
    # adam-row: ADAM17-targeted (amber)
    "AB 4":   AMBER,
    "AB 7":   AMBER,
    "C 11":   AMBER,
    "C 14":   AMBER,
    # ref-row: WT references (blue)
    "TIMP 1": ACCENT,
    "TIMP 3": ACCENT,
}

def _label_color(cn, p=1.0):
    """X-axis tick label color matching poster table row classes."""
    return _CONSTRUCT_LABEL_CLR.get(cn, TEXT_PRI)


# ── Data loading ──────────────────────────────────────────────────────────────

def load_aggregate_csv(csv_path=None):
    """
    Load aggregate_summary.csv from SharedAssets/data or a given path.
    Returns a dict:
      {metric_col: {construct: {target: (mean, sem, n)}}}
    where metric_col is one of:
      'Binding Efficiency', 'Norm Bind Med (Expr+)', 'Norm Median Ratio', etc.
    Returns None if no CSV is available.
    """
    if csv_path is None:
        csv_path = os.path.join(DATA_DIR, "aggregate_summary.csv")
    if not os.path.exists(csv_path):
        return None

    try:
        import pandas as pd
    except ImportError:
        print("  pandas not available; using hardcoded fallback data.")
        return None

    df = pd.read_csv(csv_path)

    # Identify metric columns present in the file
    # Filter out flagged failed trials before computing any stats
    if "Trial Failed" in df.columns:
        n_before = len(df)
        df = df[df["Trial Failed"] != True]
        n_dropped = n_before - len(df)
        if n_dropped:
            print(f"  Filtered {n_dropped} failed-trial rows.")

    metric_cols = [
        "Pos Med Ratio",          # PRIMARY: raw APC/FITC ratio for expr+ cells (cross-target comparable)
        "Binding Efficiency",
        "Norm Bind Med (Expr+)",
        "Norm Median Ratio",
        "Norm Intensity-Weighted Binding Index",
        "Bind Med (Expr+)",
    ]

    result = {}
    for metric in metric_cols:
        if metric not in df.columns:
            continue
        # Build {construct: {target: (mean, ci95, n)}}
        grouped = (df[df["Construct"].notna() & df["Target"].notna()]
                   .groupby(["Construct", "Target"])[metric])
        metric_data = {}
        for (construct, target), grp in grouped:
            vals = grp.dropna().values
            if len(vals) == 0:
                continue
            mean = float(np.mean(vals))
            ci   = _ci95(vals)
            metric_data.setdefault(construct, {})[target] = (mean, ci, len(vals))
        result[metric] = metric_data

    print(f"  Loaded data from {csv_path}: {list(result.keys())}")
    return result


def load_selectivity_csv(csv_path=None):
    """
    Load precomputed selectivity_summary.csv.
    Returns (sel_data, anova_p):
      sel_data: {construct: {target: {"mean": float, "ci95": float, "count": int}}}
      anova_p:  {construct: float}
    """
    if csv_path is None:
        csv_path = os.path.join(DATA_DIR, "selectivity_summary.csv")
    if not os.path.exists(csv_path):
        print(f"  ⚠  selectivity_summary.csv not found at {csv_path}")
        return None, {}
    try:
        import pandas as pd
    except ImportError:
        print("  pandas not available; cannot load selectivity_summary.csv")
        return None, {}

    df = pd.read_csv(csv_path)
    sel_data = {}
    anova_p  = {}

    for _, row in df.iterrows():
        cn = str(row["Construct"]).strip()
        t  = str(row["Target"]).strip()
        sel_data.setdefault(cn, {})[t] = {
            "mean":  float(row["Mean"]),
            "ci95":  float(row["95% CI"]),
            "count": int(row["Count"]) if not pd.isna(row["Count"]) else 1,
        }
        anova_p[cn] = float(row["ANOVA_p"])

    n_constructs = len(sel_data)
    n_targets    = len({t for c in sel_data.values() for t in c})
    print(f"  Loaded selectivity_summary.csv: {n_constructs} constructs × {n_targets} targets")
    return sel_data, anova_p


# ── Figure 1 (POSTER PRIMARY): Mean Median Binding Ratio grouped bar chart ───

def fig_bind_med_expr_pos(sel_data=None, anova_p=None):
    """
    Primary poster figure: Mean Median Binding Ratio (Bind/Expr).
    Reads directly from precomputed selectivity_summary.csv — no re-computation.
    Y-axis: APC/FITC median ratio for FITC+ cells (cross-target comparable).
    Error bars = 95% CI from selectivity_summary.csv.
    X-axis labels colored by construct group; green = ANOVA significant (p<0.05).
    """
    if sel_data is None:
        print("  ⚠  No selectivity data. Cannot generate poster figure.")
        _generate_placeholder_fig(
            os.path.join(OUTDIR, "fig_bind_med_expr_pos.png"),
            "DATA PENDING\nselectivity_summary.csv"
        )
        return
    if anova_p is None:
        anova_p = {}

    constructs = sorted(sel_data.keys())
    targets    = ["MMP2", "MMP9"]
    colors     = {"MMP2": RED, "MMP9": GREEN}
    bar_width  = 0.32
    x          = np.arange(len(constructs))

    fig, ax = plt.subplots(figsize=(16, 7))
    fig.patch.set_facecolor(DARK_BG)
    ax.set_facecolor(PANEL_BG)

    for i, t in enumerate(targets):
        means = [sel_data.get(cn, {}).get(t, {}).get("mean", 0) for cn in constructs]
        ci95s = [sel_data.get(cn, {}).get(t, {}).get("ci95", 0) for cn in constructs]
        offset = (i - 0.5) * bar_width
        ax.bar(x + offset, means, bar_width,
               color=colors[t], alpha=0.85, label=t, zorder=3)
        ax.errorbar(x + offset, means, yerr=ci95s,
                    fmt="none", ecolor="#f8fafc", capsize=3, lw=1.0, zorder=4)

    # TIMP3-WT reference dashed lines
    legend_extra = []
    for cn_ref, t_ref, clr_ref, lbl_ref in [
        ("TIMP 3", "MMP2", RED,   "TIMP3-WT (MMP2)"),
        ("TIMP 3", "MMP9", GREEN, "TIMP3-WT (MMP9)"),
    ]:
        val = sel_data.get(cn_ref, {}).get(t_ref, {}).get("mean")
        if val is not None:
            ax.axhline(val, color=clr_ref, linestyle="--", lw=1.2, alpha=0.55, zorder=2)
            legend_extra.append(plt.Line2D([0], [0], color=clr_ref, ls="--", lw=1.5,
                                           label=lbl_ref))

    # Compute y-top per construct (bar tip + CI) for bracket placement
    y_tops = []
    for cn in constructs:
        m9 = sel_data.get(cn, {}).get("MMP9", {})
        m2 = sel_data.get(cn, {}).get("MMP2", {})
        top = max(
            m9.get("mean", 0) + m9.get("ci95", 0),
            m2.get("mean", 0) + m2.get("ci95", 0),
        )
        y_tops.append(top)

    y_max_data  = max(y_tops) if y_tops else 0.4
    bracket_gap = y_max_data * 0.04

    # Significance brackets above bars (skip TIMP references)
    for xi, cn in enumerate(constructs):
        p   = anova_p.get(cn, 1.0)
        lbl = sig_label(p)
        if lbl == "n.s." or cn.startswith("TIMP"):
            continue
        y_bk = y_tops[xi] + bracket_gap
        ax.annotate("", xy=(xi + bar_width / 2, y_bk),
                    xytext=(xi - bar_width / 2, y_bk),
                    arrowprops=dict(arrowstyle="-", color=ACCENT, lw=1.5), zorder=5)
        ax.text(xi, y_bk + bracket_gap * 0.5, lbl,
                ha="center", va="bottom", color=ACCENT, fontsize=13,
                fontweight="bold", zorder=6)

    # Vertical separator before TIMP group
    timp_start = next((i for i, cn in enumerate(constructs) if cn.startswith("TIMP")), None)
    if timp_start is not None:
        ax.axvline(timp_start - 0.55, color="#334155", lw=1.0, ls="--", zorder=2)

    ax.set_xticks(x)
    tick_lbls = ax.set_xticklabels(constructs, fontsize=10, rotation=45, ha="right")
    for lbl_obj, cn in zip(tick_lbls, constructs):
        lbl_obj.set_color(_label_color(cn, anova_p.get(cn, 1.0)))
        lbl_obj.set_fontweight("bold")

    ax.set_ylabel("Median Binding Ratio (Binding/Expression)", fontsize=13)
    ax.set_ylim(0, y_max_data * 1.28)
    ax.grid(axis="y", alpha=0.4, zorder=1)
    ax.set_title("MMP9 vs MMP2 — Median Binding Ratio  (Binding/Expression, 95% CI)",
                 fontsize=18, color=TEXT_PRI, pad=22)

    legend_handles = [
        mpatches.Patch(color=RED,   label="MMP2"),
        mpatches.Patch(color=GREEN, label="MMP9"),
    ] + legend_extra
    ax.legend(handles=legend_handles, loc="upper left", fontsize=11, framealpha=0.4)

    plt.tight_layout()
    out = os.path.join(OUTDIR, "fig_bind_med_expr_pos.png")
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor=DARK_BG)
    plt.close(fig)
    print(f"  ✓  Saved: {out}")


def _generate_placeholder_fig(path, message):
    """Saves a labelled placeholder image when data is not yet available."""
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.patch.set_facecolor(DARK_BG)
    ax.set_facecolor(PANEL_BG)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.text(0.5, 0.5, message, ha="center", va="center",
            fontsize=18, color=AMBER, fontweight="bold",
            transform=ax.transAxes, wrap=True)
    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=DARK_BG)
    plt.close(fig)
    print(f"  ⚠  Placeholder saved: {path}")


# ── Figure 2: MMP9 vs MMP2 Binding Efficiency bar chart (paper/presentation) ─

def fig_mmp9_vs_mmp2(csv_data=None):
    """Binding Efficiency (DP/FITC+) — used in paper and presentation."""
    metric_key = "Binding Efficiency"
    data = None
    if csv_data and metric_key in csv_data:
        data = csv_data[metric_key]
    if data is None:
        data = _BIND_EFF_FALLBACK

    constructs = [c for c in CONSTRUCT_ORDER if c in data]
    targets    = ["MMP2", "MMP9"]
    colors     = {"MMP2": RED, "MMP9": GREEN}
    bar_width  = 0.32
    x          = np.arange(len(constructs))

    fig, ax = plt.subplots(figsize=(14, 7))
    fig.patch.set_facecolor(DARK_BG)
    ax.set_facecolor(PANEL_BG)

    for i, t in enumerate(targets):
        means = [data.get(cn, {}).get(t, (0, 0, 0))[0] for cn in constructs]
        sems  = [data.get(cn, {}).get(t, (0, 0, 0))[1] for cn in constructs]
        offset = (i - 0.5) * bar_width
        ax.bar(x + offset, means, bar_width,
               color=colors[t], alpha=0.85, label=t)
        ax.errorbar(x + offset, means, yerr=sems,
                    fmt="none", ecolor="#f8fafc", capsize=4, lw=1.2)

    selective = ["C 12", "C 15", "AB 6", "TIMP 3"]
    for cn in selective:
        if cn not in constructs:
            continue
        xi = constructs.index(cn)
        m2 = data.get(cn, {}).get("MMP2", (0, 0, 0))[0]
        m9 = data.get(cn, {}).get("MMP9", (0, 0, 0))[0]
        y_top = max(m2, m9) + 0.07
        p     = _ANOVA_P_BIND_EFF.get(cn, 1.0)
        lbl   = "Ref" if cn == "TIMP 3" else sig_label(p)
        if lbl in ("n.s.", "Ref"):
            continue
        ax.annotate("", xy=(xi + bar_width / 2, y_top + 0.025),
                    xytext=(xi - bar_width / 2, y_top + 0.025),
                    arrowprops=dict(arrowstyle="-", color=ACCENT, lw=1.5))
        ax.text(xi, y_top + 0.04, lbl, ha="center", va="bottom",
                color=ACCENT, fontsize=16, fontweight="bold")

    for cn in [c for c in constructs if c not in selective and c != "TIMP 3"]:
        xi = constructs.index(cn)
        m9 = data.get(cn, {}).get("MMP9", (0, 0, 0))[0]
        ax.text(xi, m9 + 0.08, "n.s.", ha="center", va="bottom",
                color=TEXT_SEC, fontsize=14, style="italic")

    timp3_m2 = data.get("TIMP 3", {}).get("MMP2", (0, 0, 0))[0]
    timp3_m9 = data.get("TIMP 3", {}).get("MMP9", (0, 0, 0))[0]
    ax.axhline(timp3_m2, color=RED,   linestyle="--", lw=1.0, alpha=0.5)
    ax.axhline(timp3_m9, color=GREEN, linestyle="--", lw=1.0, alpha=0.5)
    ax.axvline(len(constructs) - 1 - 0.55, color="#334155", lw=1.0, ls="--")

    ax.text(1.0, 0.97, "MMP9 Specific", ha="center", va="top",
            transform=ax.get_xaxis_transform(), fontsize=15,
            color=ACCENT, fontweight="bold", alpha=0.9)
    ax.text(5.5, 0.97, "Non-Specific Controls / ADAM Targets", ha="center", va="top",
            transform=ax.get_xaxis_transform(), fontsize=14,
            color=TEXT_SEC, style="italic", alpha=0.9)

    ax.set_xticks(x)
    ax.set_xticklabels(constructs, fontsize=15)
    ax.set_ylabel("Binding Efficiency (DP / FITC⁺)", fontsize=17)
    ax.set_ylim(0, 1.18)
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
    ax.grid(axis="y", alpha=0.4)
    ax.set_title("MMP9 vs MMP2 Binding Efficiency — Specificity Validation",
                 fontsize=22, color=TEXT_PRI, pad=22)

    legend_handles = [
        mpatches.Patch(color=RED,   label="MMP2"),
        mpatches.Patch(color=GREEN, label="MMP9"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=15)

    plt.tight_layout()
    out = os.path.join(OUTDIR, "fig_mmp9_vs_mmp2.png")
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor=DARK_BG)
    plt.close(fig)
    print(f"  ✓  Saved: {out}")


# ── Figure 3: Selectivity ratio (MMP9 / MMP2) ────────────────────────────────

def fig_specificity_ratio(csv_data=None):
    """MMP9/MMP2 Binding Efficiency ratio — used in paper and presentation."""
    metric_key = "Binding Efficiency"
    data = None
    if csv_data and metric_key in csv_data:
        data = csv_data[metric_key]
    if data is None:
        data = _BIND_EFF_FALLBACK

    constructs = ["C 12", "C 15", "AB 6", "C 13", "AB 5", "TIMP 3"]
    ratios = []
    for cn in constructs:
        m9 = data.get(cn, {}).get("MMP9", (1e-9, 0, 0))[0]
        m2 = data.get(cn, {}).get("MMP2", (1e-9, 0, 0))[0]
        ratios.append(m9 / max(m2, 1e-9))

    fig, ax = plt.subplots(figsize=(10, 5.5))
    fig.patch.set_facecolor(DARK_BG)
    ax.set_facecolor(PANEL_BG)

    bar_colors = []
    for cn, r in zip(constructs, ratios):
        p = _ANOVA_P_BIND_EFF.get(cn, 1.0)
        if p < 0.05 and r > 2:
            bar_colors.append(GREEN)
        elif cn == "TIMP 3":
            bar_colors.append(ACCENT)
        else:
            bar_colors.append(TEXT_SEC)

    ax.bar(constructs, ratios, color=bar_colors, alpha=0.85, width=0.55)

    for i, (cn, r) in enumerate(zip(constructs, ratios)):
        p   = _ANOVA_P_BIND_EFF.get(cn, 1.0)
        lbl = sig_label(p)
        clr = ACCENT if p < 0.05 else TEXT_SEC
        ax.text(i, r + 0.1, lbl, ha="center", va="bottom",
                fontsize=10, fontweight="bold", color=clr)

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


# ── Figure 4: Mean Median Binding Ratio heatmap ──────────────────────────────

def fig_binding_heatmap(sel_data=None, anova_p=None):
    """
    Mean Median Binding Ratio heatmap — constructs (x) × targets (y), landscape.
    Rows: MMP2, MMP9, ADAM17.  Cols: all constructs, alphabetical.
    Significance shown in column labels (*** / ** / * appended, colored green).
    No significance text placed inside/below cells — labels carry the annotation.
    """
    if sel_data is None:
        print("  ⚠  No selectivity data; skipping heatmap.")
        _generate_placeholder_fig(
            os.path.join(OUTDIR, "fig_binding_heatmap.png"),
            "DATA PENDING\nselectivity_summary.csv"
        )
        return
    if anova_p is None:
        anova_p = {}

    constructs = sorted(sel_data.keys())
    targets    = [t for t in ["MMP2", "MMP9", "ADAM17"]
                  if any(t in sel_data.get(cn, {}) for cn in constructs)]

    mat = np.array([
        [sel_data.get(cn, {}).get(t, {}).get("mean", np.nan) for cn in constructs]
        for t in targets
    ])
    vmax = float(np.nanmax(mat)) if not np.all(np.isnan(mat)) else 0.4

    fig, ax = plt.subplots(figsize=(20, 4.5))
    fig.patch.set_facecolor(DARK_BG)
    ax.set_facecolor(DARK_BG)

    cmap = LinearSegmentedColormap.from_list(
        "custom", ["#0f172a", "#1e3a5f", "#38bdf8", "#34d399", "#fbbf24"], N=256)
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=vmax)

    for i in range(len(targets)):
        for j in range(len(constructs)):
            v = mat[i, j]
            if not np.isnan(v):
                txt_color = "black" if v > 0.7 * vmax else TEXT_PRI
                ax.text(j, i, f"{v:.3f}", ha="center", va="center",
                        fontsize=10, color=txt_color, fontweight="bold")

    # Column labels: construct name + significance stars, colored by group/sig
    col_labels = []
    for cn in constructs:
        p     = anova_p.get(cn, 1.0)
        stars = sig_label(p)
        col_labels.append(cn if stars == "n.s." else f"{cn}{stars}")

    ax.set_xticks(range(len(constructs)))
    col_tick_lbls = ax.set_xticklabels(col_labels, fontsize=11, rotation=45, ha="right")
    for lbl_obj, cn in zip(col_tick_lbls, constructs):
        lbl_obj.set_color(_label_color(cn, anova_p.get(cn, 1.0)))
        lbl_obj.set_fontweight("bold")

    ax.set_yticks(range(len(targets)))
    ax.set_yticklabels(targets, fontsize=14, color=TEXT_PRI, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, fraction=0.015, pad=0.04)
    cbar.set_label("Mean Median Binding Ratio (Bind/Expr)", color=TEXT_PRI, fontsize=12)
    cbar.ax.yaxis.set_tick_params(color=TEXT_PRI)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=TEXT_PRI, fontsize=10)

    ax.set_title(
        "Mean Median Binding Ratio Heatmap — per construct × target  "
        "(label stars = ANOVA sig: *** <0.001  ** <0.01  * <0.05)",
        fontsize=14, color=TEXT_PRI, pad=12,
    )
    plt.tight_layout()
    out = os.path.join(OUTDIR, "fig_binding_heatmap.png")
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor=DARK_BG)
    plt.close(fig)
    print(f"  ✓  Saved: {out}")


# ── Figure 5: Norm Bind Med (Expr+) heatmap ──────────────────────────────────

def fig_bind_med_heatmap(csv_data=None):
    """
    Heatmap of Norm Bind Med (Expr+) across constructs × targets.
    Within-target normalized (each target independently to TIMP3-WT).
    Used in paper; shows intensity of binding per expressing cell.
    """
    metric_key = "Norm Bind Med (Expr+)"
    data = None
    if csv_data and metric_key in csv_data:
        data = csv_data[metric_key]

    if data is None:
        print(f"  ⚠  No CSV data for '{metric_key}'. Skipping heatmap variant.")
        _generate_placeholder_fig(
            os.path.join(OUTDIR, "fig_bind_med_heatmap.png"),
            f"DATA PENDING\n'{metric_key}'\nfrom aggregate_summary.csv"
        )
        return

    constructs = [c for c in CONSTRUCT_ORDER if c in data]
    targets    = [t for t in ["MMP2", "MMP9", "ADAM17"]
                  if any(t in data.get(c, {}) for c in constructs)]

    mat = np.array([
        [data.get(cn, {}).get(t, (np.nan, 0, 0))[0] for t in targets]
        for cn in constructs
    ])

    vmax = np.nanmax(mat) if not np.all(np.isnan(mat)) else 2.0

    fig, ax = plt.subplots(figsize=(10, 13))
    fig.patch.set_facecolor(DARK_BG)
    ax.set_facecolor(DARK_BG)

    cmap = LinearSegmentedColormap.from_list(
        "custom", ["#0f172a", "#1e3a5f", "#38bdf8", "#34d399", "#fbbf24"], N=256)
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=vmax)

    for i in range(len(constructs)):
        for j in range(len(targets)):
            v = mat[i, j]
            if not np.isnan(v):
                txt_color = "black" if v > 0.7 * vmax else TEXT_PRI
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=16, color=txt_color, fontweight="bold")

    ax.set_xticks(range(len(targets)))
    ax.set_xticklabels(targets, fontsize=20, color=TEXT_PRI)
    ax.set_yticks(range(len(constructs)))
    ax.set_yticklabels(constructs, fontsize=20, color=TEXT_PRI)
    if "TIMP 3" in constructs:
        ax.axhline(constructs.index("TIMP 3") - 0.5, color="#334155", lw=2.0, ls="--")

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.12)
    cbar.set_label("Norm. Binding Median (Expr⁺) [vs TIMP3-WT]",
                   color=TEXT_PRI, fontsize=16)
    cbar.ax.yaxis.set_tick_params(color=TEXT_PRI)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=TEXT_PRI)

    ax.set_title(
        "Binding Median (Expr⁺) Heatmap\n"
        "(within-target norm. to TIMP3-WT; 1.0 = TIMP3-WT level)",
        fontsize=24, color=TEXT_PRI, pad=15,
    )
    plt.tight_layout()
    out = os.path.join(OUTDIR, "fig_bind_med_heatmap.png")
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor=DARK_BG)
    plt.close(fig)
    print(f"  ✓  Saved: {out}")


# ── Figure 6: 6-stage pipeline diagram ───────────────────────────────────────

def fig_pipeline():
    stages = [
        ("RFdiffusion\n(Backbone Gen.)", "AB / C / EF loops\n6–24 aa expansion", "#3b82f6"),
        ("ProteinMPNN\n(Sequence Design)", "1000 seqs / target\nT=0.2, fixed scaffold", "#8b5cf6"),
        ("AlphaFold Server\n(Co-folding)", "Batch JSON, 5 targets\nipTM / PAE / pLDDT", "#06b6d4"),
        ("Best Binder\nSelection", "T-score > 2.0\nMulti-metric consensus", "#10b981"),
        ("Twist Bioscience\n(DNA Synthesis)", "13 constructs\nCodon-opt., Dec 2025", "#f59e0b"),
        ("Yeast Display\n& Flow Cytometry", "Bind Med (Expr⁺)\nANOVA + Tukey-HSD", "#ef4444"),
    ]

    fig, ax = plt.subplots(figsize=(14, 3.5))
    fig.patch.set_facecolor(DARK_BG)
    ax.set_facecolor(DARK_BG)
    ax.axis("off")

    box_w, box_h = 1.9, 1.6
    gap   = 0.45
    total = len(stages) * box_w + (len(stages) - 1) * gap
    x0    = (14 - total) / 2

    for i, (title, subtitle, color) in enumerate(stages):
        x = x0 + i * (box_w + gap)
        rect = FancyBboxPatch(
            (x, 0.6), box_w, box_h,
            boxstyle="round,pad=0.08",
            facecolor=color + "22", edgecolor=color, linewidth=2,
            transform=ax.transData, clip_on=False,
        )
        ax.add_patch(rect)
        ax.text(x + box_w / 2, 0.6 + box_h - 0.25, title,
                ha="center", va="top", fontsize=12, fontweight="bold",
                color=color, transform=ax.transData)
        ax.text(x + box_w / 2, 0.6 + 0.25, subtitle,
                ha="center", va="bottom", fontsize=10, color=TEXT_SEC,
                transform=ax.transData)
        if i < len(stages) - 1:
            ax.annotate("", xy=(x + box_w + gap, 0.6 + box_h / 2),
                        xytext=(x + box_w, 0.6 + box_h / 2),
                        arrowprops=dict(arrowstyle="->", color=TEXT_SEC,
                                        lw=1.5, mutation_scale=14))

    ax.set_xlim(0, 14)
    ax.set_ylim(0, 2.8)
    ax.set_title("Computational-to-Experimental Design Pipeline",
                 fontsize=18, color=TEXT_PRI, pad=12)
    plt.tight_layout()
    out = os.path.join(OUTDIR, "fig_pipeline.png")
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor=DARK_BG)
    plt.close(fig)
    print(f"  ✓  Saved: {out}")


# ── JSON export for FcsChart (slides) ────────────────────────────────────────

def export_experimental_binding_json(sel_data=None, anova_p=None):
    """
    Write experimental_binding.json consumed by FcsChart.vue.
    Format: {target: [{Construct, "Pos Med Ratio", "Ratio CI", N}, ...]}
    Uses precomputed Mean and 95% CI from selectivity_summary.csv.
    """
    if sel_data is None:
        print("  ⚠  No selectivity data. Skipping JSON export.")
        return

    all_constructs = sorted(sel_data.keys())
    targets = ["MMP2", "MMP9", "ADAM17", "ADAM10", "MMP3"]

    output = {}
    for t in targets:
        entries = []
        for cn in all_constructs:
            tdata = sel_data.get(cn, {}).get(t)
            if tdata is None:
                continue
            entries.append({
                "Construct": cn,
                "Pos Med Ratio": round(tdata["mean"], 4),
                "Ratio CI": round(tdata["ci95"], 4),
                "N": tdata["count"],
            })
        if entries:
            output[t] = entries

    out_path = os.path.join(DATA_DIR, "experimental_binding.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"  ✓  Saved: {out_path}")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate publication figures for De Novo Binder Generation thread."
    )
    parser.add_argument(
        "--csv", default=None,
        help="Path to aggregate_summary.csv (for Binding Efficiency / Norm figures).",
    )
    parser.add_argument(
        "--selectivity-csv", default=None,
        help=(
            "Path to selectivity_summary.csv (precomputed stats). "
            "If omitted, checks SharedAssets/data/De_Novo_Binder_Generation/selectivity_summary.csv."
        ),
    )
    args = parser.parse_args()

    print("Generating figures …")
    csv_data            = load_aggregate_csv(args.csv)
    sel_data, anova_p   = load_selectivity_csv(args.selectivity_csv)

    # Poster primary figure — uses precomputed selectivity_summary.csv
    print("\n[Poster] Mean Median Binding Ratio (Bind/Expr):")
    fig_bind_med_expr_pos(sel_data, anova_p)

    # Paper / presentation supporting figures
    print("\n[Paper / Slides] Binding Efficiency (DP/FITC+):")
    fig_mmp9_vs_mmp2(csv_data)

    print("\n[Paper / Slides] Selectivity ratio:")
    fig_specificity_ratio(csv_data)

    print("\n[Poster / Slides] Mean Median Binding Ratio heatmap:")
    fig_binding_heatmap(sel_data, anova_p)

    print("\n[Paper] Norm Bind Med (Expr+) heatmap:")
    fig_bind_med_heatmap(csv_data)

    print("\n[All] Pipeline diagram:")
    fig_pipeline()

    print("\n[Slides] FcsChart JSON (Pos Med Ratio):")
    export_experimental_binding_json(sel_data, anova_p)

    print(f"\nDone. All figures in:\n  {OUTDIR}")


if __name__ == "__main__":
    main()
