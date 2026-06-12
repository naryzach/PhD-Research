#!/usr/bin/env python3
"""
Compare the AlphaFold-driven predictions behind the Dec-2025 Twist order
(Local/Twist_Order_Dec2025) against the experimental FCS results
(Local/Aggregate_FCS_Analysis), and build a comprehensive per-construct x
per-target rollup of every metric in the *_Renamed_Analysis aggregation
(not just the selectivity slice).

Outputs
-------
Local/Construct_Metric_Summary/
    construct_target_metric_summary.csv        (long: Construct,Target,Metric,Mean,Std,SEM,N)
    construct_target_metric_summary_wide.csv   (wide: one row per Construct x Target)
Local/Prediction_vs_Result_Analysis/
    predicted_vs_actual.csv                    (AF metrics + experimental means per construct)
    correlation_summary.csv                    (Pearson/Spearman per target x metric pair)
    prediction_scorecard.csv                   (categorical expected-vs-observed verdicts)
    plot_LpLDDT_vs_binding.png
    plot_overall_affinity_rank.png
"""
import os, glob
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Local")  # repo/Local (gitignored outputs)
AFM  = os.path.join(ROOT, "AlphaFoldMetrics")
AGGF = os.path.join(ROOT, "Aggregate_FCS_Analysis", "aggregate_summary.csv")
OUT_PRED = os.path.join(ROOT, "Prediction_vs_Result_Analysis")
OUT_SUMM = os.path.join(ROOT, "Construct_Metric_Summary")
os.makedirs(OUT_PRED, exist_ok=True)
os.makedirs(OUT_SUMM, exist_ok=True)

# ----------------------------------------------------------------------------
# 1. AlphaFold predicted metrics (per variant x target).  The ordered library
#    is a mix of the original "Exp" designs and ProteinMPNN-redesigned "PMPNN"
#    sequences, so both matrix sets are loaded and the right one is used per
#    construct (flavor column records which).
# ----------------------------------------------------------------------------
AF_METRICS = ["ipTM", "LpLDDT", "PAE", "ApTM", "BpTM", "pTM", "pLDDT"]
AF_TARGETS = ["adam10", "adam17", "mmp2", "mmp9"]

def load_matrices(flavor):
    d = {}
    for p in glob.glob(os.path.join(AFM, f"Comparison_Matrix_*_{flavor}.xlsx")):
        name = os.path.basename(p).split("Comparison_Matrix_")[1].split(f"_{flavor}.xlsx")[0]
        df = pd.read_excel(p)
        df["Loop"] = df["Loop"].ffill()
        df["Variant"] = df["Variant"].astype(str).str.lower()
        d[name] = df
    return d

EXP = load_matrices("Exp")
PMP = load_matrices("PMPNN")

# Construct -> (loop, sequence, AF source flavor, expected-label)
CONSTRUCTS = {
    "AB 1": ("ab", "dgptge",         "PMPNN", "M9>M2"),
    "AB 2": ("ab", "eversghkvke",    "Exp",   "M9+"),
    "AB 3": ("ab", "kgpyge",         "PMPNN", "High"),
    "AB 4": ("ab", "knpdgtlt",       "Exp",   "A17+, A17>A10"),
    "AB 5": ("ab", "patptstrgaggee", "Exp",   "Low"),
    "AB 6": ("ab", "tdtfptanwtgev",  "Exp",   "M9>M2"),
    "AB 7": ("ab", "tlpdgske",       "Exp",   "A17+, A17>A10"),
    "C 11": ("c",  "anpeyc",         "PMPNN", "A17+"),
    "C 12": ("c",  "asgpitvngetiw",  "Exp",   "M9+, M9>M2"),
    "C 13": ("c",  "asveavetgfs",    "Exp",   "Low"),
    "C 14": ("c",  "ggnygsck",       "Exp",   "A17>A10"),
    "C 15": ("c",  "ltqeelpdpnavspc","Exp",   "M9+, M9>M2"),
    "C 16": ("c",  "sveslc",         "PMPNN", "High"),
}

def af_lookup(metric, loop, seq, flavor, target):
    src = EXP if flavor == "Exp" else PMP
    df = src.get(metric)
    if df is None:
        return np.nan
    hit = df[(df["Loop"] == loop) & (df["Variant"] == seq)]
    return hit[target].values[0] if len(hit) else np.nan

af_rows = []
for c, (loop, seq, flav, pred) in CONSTRUCTS.items():
    rec = {"Construct": c, "Loop": loop, "Sequence": seq, "AF_Source": flav, "Expected": pred}
    for met in AF_METRICS:
        for t in AF_TARGETS:
            rec[f"{met}_{t}"] = af_lookup(met, loop, seq, flav, t)
    af_rows.append(rec)
af = pd.DataFrame(af_rows)

# ----------------------------------------------------------------------------
# 2. Experimental rollup from aggregate_summary (valid / non-failed trials).
#    Build BOTH the comprehensive per-construct x per-target summary across
#    every metric, and the per-target binding means used for the comparison.
# ----------------------------------------------------------------------------
agg = pd.read_csv(AGGF)
valid = agg[agg["Trial Failed"] != True].copy()

METRIC_COLS = [
    "Pos Med Ratio", "Norm Median Ratio", "Pos Mean Ratio", "Norm Mean Ratio",
    "Double+ %", "Expr+ %", "Gated Events",
    "Bind Med (Expr+)", "Norm Bind Med (Expr+)", "Bind Mean (Expr+)", "Norm Bind Mean (Expr+)",
    "Expr Med (Bind+)", "Norm Expr Med (Bind+)", "Binding Efficiency",
    "Intensity-Weighted Binding Index", "Norm Intensity-Weighted Binding Index",
]

g = valid.groupby(["Construct", "Target"])
mean_df = g[METRIC_COLS].mean()
std_df  = g[METRIC_COLS].std()
n_ser   = g.size().rename("N")

# ---- comprehensive long-format summary (request: "summarize everything") ----
long_rows = []
for (constr, targ), row in mean_df.iterrows():
    n = int(n_ser.loc[(constr, targ)])
    for m in METRIC_COLS:
        mu = row[m]
        sd = std_df.loc[(constr, targ), m]
        sem = sd / np.sqrt(n) if (n and pd.notna(sd)) else (0.0 if n == 1 else np.nan)
        long_rows.append({"Construct": constr, "Target": targ, "Metric": m,
                          "Mean": mu, "Std": sd, "SEM": sem, "N": n})
long_df = pd.DataFrame(long_rows)
long_df.to_csv(os.path.join(OUT_SUMM, "construct_target_metric_summary.csv"), index=False)

# wide version: one row per construct x target, mean of each metric + N
wide = mean_df.copy()
wide["N_trials"] = n_ser
wide = wide.reset_index()
wide.to_csv(os.path.join(OUT_SUMM, "construct_target_metric_summary_wide.csv"), index=False)

# ----------------------------------------------------------------------------
# 3. Merge predicted + actual into one comparison table (single-loop constructs)
# ----------------------------------------------------------------------------
VERIF = {"ADAM17": "adam17", "MMP2": "mmp2", "MMP9": "mmp9"}   # usable targets
EXP_KEY = "Norm Median Ratio"                                  # primary readout

pva = af.copy()
for T, t in VERIF.items():
    for met in ["Norm Median Ratio", "Norm Intensity-Weighted Binding Index",
                "Binding Efficiency", "Double+ %"]:
        col = f"exp_{met}_{T}"
        pva[col] = [mean_df.loc[(c, T), met] if (c, T) in mean_df.index else np.nan
                    for c in pva["Construct"]]
    pva[f"n_{T}"] = [int(n_ser.loc[(c, T)]) if (c, T) in n_ser.index else 0
                     for c in pva["Construct"]]
pva.to_csv(os.path.join(OUT_PRED, "predicted_vs_actual.csv"), index=False)

# ----------------------------------------------------------------------------
# 4. Numerical correlations: AF metric vs experimental readout, per target,
#    plus pooled across all construct x target pairs.
# ----------------------------------------------------------------------------
corr_rows = []
AF_FOR_CORR = ["ipTM", "LpLDDT", "PAE", "ApTM", "BpTM"]
EXP_FOR_CORR = ["Norm Median Ratio", "Norm Intensity-Weighted Binding Index", "Binding Efficiency"]

def add_corr(label, xs, ys, amet, emet):
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    ok = np.isfinite(xs) & np.isfinite(ys)
    xs, ys = xs[ok], ys[ok]
    if len(xs) >= 4 and np.ptp(xs) > 0 and np.ptp(ys) > 0:
        r, p = stats.pearsonr(xs, ys)
        rho, pp = stats.spearmanr(xs, ys)
        corr_rows.append({"Scope": label, "AF_Metric": amet, "Exp_Metric": emet,
                          "N": len(xs), "Pearson_r": r, "Pearson_p": p,
                          "Spearman_rho": rho, "Spearman_p": pp})

for emet in EXP_FOR_CORR:
    for amet in AF_FOR_CORR:
        pooled_x, pooled_y = [], []
        for T, t in VERIF.items():
            xs = [af_lookup(amet, *CONSTRUCTS[c][:3], t) for c in CONSTRUCTS]
            ys = [mean_df.loc[(c, T), emet] if (c, T) in mean_df.index else np.nan
                  for c in CONSTRUCTS]
            add_corr(f"per-target:{T}", xs, ys, amet, emet)
            pooled_x += xs
            pooled_y += ys
        add_corr("pooled(all targets)", pooled_x, pooled_y, amet, emet)

corr = pd.DataFrame(corr_rows)
corr.to_csv(os.path.join(OUT_PRED, "correlation_summary.csv"), index=False)

# ----------------------------------------------------------------------------
# 5. Categorical scorecard: did each expected label pan out?
# ----------------------------------------------------------------------------
def binding(c, T):
    return mean_df.loc[(c, T), EXP_KEY] if (c, T) in mean_df.index else np.nan

# overall affinity rank (mean Norm Median Ratio across the 3 usable targets)
overall = {c: np.nanmean([binding(c, T) for T in VERIF]) for c in CONSTRUCTS}
ranked = sorted([c for c in CONSTRUCTS if np.isfinite(overall[c])],
                key=lambda c: -overall[c])
rank_of = {c: i + 1 for i, c in enumerate(ranked)}
n_rank = len(ranked)

score_rows = []
for c, (loop, seq, flav, pred) in CONSTRUCTS.items():
    a17, m2, m9 = binding(c, "ADAM17"), binding(c, "MMP2"), binding(c, "MMP9")
    ov = overall[c]
    rk = rank_of.get(c)
    notes, verdict = [], []
    labels = [p.strip() for p in pred.split(",")]
    for lab in labels:
        if lab == "High":
            if rk is None:
                v = "Untestable (no data)"
            elif rk <= max(1, n_rank // 3):
                v = "Hit"
            elif rk <= 2 * n_rank // 3:
                v = "Partial (mid-pack)"
            else:
                v = "Miss (low)"
            notes.append(f"overall rank {rk}/{n_rank} (mean NMR {ov:.2f})")
        elif lab == "Low":
            if rk is None:
                v = "Untestable (no data)"
            elif rk >= 2 * n_rank // 3:
                v = "Hit"
            elif rk > n_rank // 3:
                v = "Partial (mid-pack)"
            else:
                v = "Miss (bound well)"
            notes.append(f"overall rank {rk}/{n_rank} (mean NMR {ov:.2f})")
        elif lab == "A17+":
            # ADAM17 is the construct's strongest of the 3 usable targets?
            vals = {"ADAM17": a17, "MMP2": m2, "MMP9": m9}
            best = max((t for t in vals if np.isfinite(vals[t])),
                       key=lambda t: vals[t], default=None)
            if not np.isfinite(a17):
                v = "Untestable (no data)"
            elif best == "ADAM17":
                v = "Hit"
            else:
                v = "Partial" if a17 >= 1.0 else "Miss"
            notes.append(f"A17 NMR {a17:.2f}; strongest target = {best}")
        elif lab == "A17>A10":
            v = "Untestable (no ADAM10 data)"
            notes.append("ADAM10 has no valid non-control trials")
        elif lab == "M9+":
            if not np.isfinite(m9):
                v = "Untestable (no data)"
            else:
                v = "Untestable (MMP9 ANOVA n.s.)"
            notes.append(f"M9 NMR {m9:.2f}; MMP9 shows no between-construct spread")
        elif lab == "M9>M2":
            if np.isfinite(m9) and np.isfinite(m2):
                v = "Hit (direction)" if m9 > m2 else "Miss (direction)"
                notes.append(f"M9 {m9:.2f} vs M2 {m2:.2f}; (MMP9 spread n.s.)")
            else:
                v = "Untestable (no data)"
        else:
            v = "?"
        verdict.append(f"{lab}: {v}")
    score_rows.append({
        "Construct": c, "Sequence": seq, "AF_Source": flav, "Expected": pred,
        "ADAM17_NMR": a17, "MMP2_NMR": m2, "MMP9_NMR": m9,
        "Overall_mean_NMR": ov, "Overall_rank": rk,
        "Verdict": " | ".join(verdict), "Notes": "; ".join(notes),
    })
score = pd.DataFrame(score_rows)
score.to_csv(os.path.join(OUT_PRED, "prediction_scorecard.csv"), index=False)

# ----------------------------------------------------------------------------
# 6. Plots
# ----------------------------------------------------------------------------
# 6a. LpLDDT (best dynamic-range AF metric) vs experimental NMR, per target
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
for ax, (T, t) in zip(axes, VERIF.items()):
    xs, ys, labs = [], [], []
    for c in CONSTRUCTS:
        x = af_lookup("LpLDDT", *CONSTRUCTS[c][:3], t)
        y = binding(c, T)
        if np.isfinite(x) and np.isfinite(y):
            xs.append(x); ys.append(y); labs.append(c)
    ax.scatter(xs, ys, c="#2b6cb0")
    for x, y, l in zip(xs, ys, labs):
        ax.annotate(l, (x, y), fontsize=7, xytext=(3, 3), textcoords="offset points")
    if len(xs) >= 4:
        r, p = stats.pearsonr(xs, ys)
        m, b = np.polyfit(xs, ys, 1)
        xx = np.linspace(min(xs), max(xs), 50)
        ax.plot(xx, m * xx + b, "--", c="gray", lw=1)
        ax.set_title(f"{T}\nPearson r={r:.2f} (p={p:.2f}), n={len(xs)}")
    ax.axhline(1.0, color="red", lw=0.8, ls=":")
    ax.set_xlabel("AlphaFold interface pLDDT (LpLDDT)")
    ax.set_ylabel("Experimental Norm Median Ratio")
fig.suptitle("Predicted interface confidence vs measured binding (red dashed = positive-control level)")
fig.tight_layout()
fig.savefig(os.path.join(OUT_PRED, "plot_LpLDDT_vs_binding.png"), dpi=130)
plt.close(fig)

# 6b. overall predicted affinity tier vs observed overall binding
tier_color = {"High": "#2f855a", "Low": "#c53030"}
fig, ax = plt.subplots(figsize=(9, 4.6))
xs = list(range(len(ranked)))
ys = [overall[c] for c in ranked]
cols = []
for c in ranked:
    lbl = CONSTRUCTS[c][3]
    if "High" in lbl: cols.append("#2f855a")
    elif "Low" in lbl: cols.append("#c53030")
    elif "A17" in lbl: cols.append("#dd6b20")
    elif "M9" in lbl: cols.append("#3182ce")
    else: cols.append("#718096")
ax.bar(xs, ys, color=cols)
ax.set_xticks(xs)
ax.set_xticklabels([f"{c}\n{CONSTRUCTS[c][3].split(',')[0]}" for c in ranked],
                   fontsize=7, rotation=0)
ax.axhline(1.0, color="black", lw=0.8, ls=":")
ax.set_ylabel("Observed overall binding\n(mean Norm Median Ratio, 3 targets)")
ax.set_title("Constructs ranked by measured overall binding, colored by PREDICTED tier\n"
             "green=High  red=Low  orange=ADAM17  blue=MMP9 (dotted = positive-control level)")
fig.tight_layout()
fig.savefig(os.path.join(OUT_PRED, "plot_overall_affinity_rank.png"), dpi=130)
plt.close(fig)

# 6c. observed binding separated BY TARGET (constructs ranked within each target)
def _tier_color(c):
    lbl = CONSTRUCTS[c][3]
    if "High" in lbl: return "#2f855a"
    if "Low" in lbl:  return "#c53030"
    if "A17" in lbl:  return "#dd6b20"
    if "M9" in lbl:   return "#3182ce"
    return "#718096"

fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
for ax, T in zip(axes, ["ADAM17", "MMP2", "MMP9"]):
    vals = {c: binding(c, T) for c in CONSTRUCTS}
    rk = sorted([c for c in CONSTRUCTS if np.isfinite(vals[c])], key=lambda c: -vals[c])
    xs = list(range(len(rk)))
    ax.bar(xs, [vals[c] for c in rk], color=[_tier_color(c) for c in rk])
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{c}\n{CONSTRUCTS[c][3].split(',')[0]}" for c in rk], fontsize=6, rotation=90)
    ax.axhline(1.0, color="black", lw=0.8, ls=":")
    ax.set_title(f"{T}" + ("  [ANOVA n.s.]" if T == "MMP9" else ""), fontsize=10)
    ax.set_ylabel("Observed Norm Median Ratio")
fig.suptitle("Observed binding per target — constructs ranked WITHIN each target, colored by PREDICTED tier\n"
             "green=High  red=Low  orange=ADAM17  blue=MMP9  (dotted = positive-control level)")
fig.tight_layout()
fig.savefig(os.path.join(OUT_PRED, "plot_overall_affinity_by_target.png"), dpi=130)
plt.close(fig)

# ----------------------------------------------------------------------------
print("=== scorecard ===")
print(score[["Construct", "Expected", "ADAM17_NMR", "MMP2_NMR", "MMP9_NMR",
             "Overall_rank", "Verdict"]].to_string(index=False))
print("\n=== strongest correlations (|Pearson_r|, pooled + per-target) ===")
cc = corr.reindex(corr["Pearson_r"].abs().sort_values(ascending=False).index)
print(cc.head(12).to_string(index=False))
print("\nWrote:")
for f in ["construct_target_metric_summary.csv", "construct_target_metric_summary_wide.csv"]:
    print("  ", os.path.join("Local/Construct_Metric_Summary", f))
for f in ["predicted_vs_actual.csv", "correlation_summary.csv", "prediction_scorecard.csv",
          "plot_LpLDDT_vs_binding.png", "plot_overall_affinity_rank.png"]:
    print("  ", os.path.join("Local/Prediction_vs_Result_Analysis", f))
