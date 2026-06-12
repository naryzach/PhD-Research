#!/usr/bin/env python3
"""
Ground-truth calibration of AlphaFold metrics against experimental FCS binding.

Three stages:
  1. DECOMPOSE the experimental signal (constructs x targets) into
        stickiness (construct main effect) + target baseline + target-specific
        (interaction).  The cross-target correlation r~0.7 means a single
        construct-level avidity factor dominates; we must calibrate against the
        target-SPECIFIC residual, not raw binding, or we just predict avidity.
  2. CALIBRATE each AF metric against (a) the stickiness factor, (b) the
        target-specific residual, (c) off-target selectivity gaps.
  3. EVALUATE a multi-term prediction recipe (binding_recipe.py) against
        single-metric baselines, with leave-one-out honesty.

n = 12 constructs x 3 usable targets (ADAM17, MMP2, MMP9). Small; treat as
directional. Re-run when a new ordered batch's FCS results land
(re-point AGGF / add the constructs to CONSTRUCTS).

Outputs -> Local/Calibration/
"""
import os, glob, re, sys
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))   # repo/Binder-Calibration
REPO = os.path.dirname(HERE)                          # repo root
ROOT = os.path.join(REPO, "Local")                    # repo/Local (gitignored outputs)
GEN  = os.path.join(REPO, "Generation")               # binding_recipe.py lives here (tracked)
sys.path.insert(0, GEN)
from binding_recipe import score_design, score_breakdown, RECIPE_WEIGHTS, normalize  # noqa: E402

import argparse
_ap = argparse.ArgumentParser(description="Calibrate AF-server (default) or ESMFold2 metrics vs FCS binding.")
_ap.add_argument("--esm-scores", default=None,
                 help="CSV of ESMFold2 scores to use INSTEAD of the AF-server matrices. Columns: "
                      "Construct,Target,esm_iptm,esm_ptm,esm_plddt[,esm_lplddt,esm_pae]; "
                      "Target in {ADAM17,MMP2,MMP9}. Runs the SAME Stage 2/3 on ESMFold2 metrics so "
                      "the AF3 and ESMFold2 calibrations are directly comparable.")
ARGS, _ = _ap.parse_known_args()
PRED_SOURCE = "esmfold2" if ARGS.esm_scores else "afserver"
TAG = "_esmfold2" if PRED_SOURCE == "esmfold2" else ""   # suffix outputs so AF3 results aren't overwritten

AFM  = os.path.join(ROOT, "AlphaFoldMetrics")
AGGF = os.path.join(ROOT, "Aggregate_FCS_Analysis", "aggregate_summary.csv")
OUT  = os.path.join(ROOT, "Calibration")             # outputs stay in Local/Calibration
os.makedirs(OUT, exist_ok=True)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

# ----------------------------------------------------------- AlphaFold metrics
AF_METRICS = ["ipTM", "LpLDDT", "PAE", "ApTM", "BpTM", "pTM", "pLDDT"]
EXPECT_SIGN = {"ipTM": +1, "LpLDDT": +1, "PAE": -1, "ApTM": +1, "BpTM": +1, "pTM": +1, "pLDDT": +1}
def load_matrices(flavor):
    d = {}
    for p in glob.glob(os.path.join(AFM, f"Comparison_Matrix_*_{flavor}.xlsx")):
        name = os.path.basename(p).split("Comparison_Matrix_")[1].split(f"_{flavor}.xlsx")[0]
        df = pd.read_excel(p); df["Loop"] = df["Loop"].ffill()
        df["Variant"] = df["Variant"].astype(str).str.lower(); d[name] = df
    return d
EXP_M, PMP_M = load_matrices("Exp"), load_matrices("PMPNN")
CONSTRUCTS = {
    "AB 1": ("ab", "dgptge", "PMPNN"), "AB 2": ("ab", "eversghkvke", "Exp"),
    "AB 3": ("ab", "kgpyge", "PMPNN"), "AB 4": ("ab", "knpdgtlt", "Exp"),
    "AB 5": ("ab", "patptstrgaggee", "Exp"), "AB 6": ("ab", "tdtfptanwtgev", "Exp"),
    "AB 7": ("ab", "tlpdgske", "Exp"), "C 11": ("c", "anpeyc", "PMPNN"),
    "C 12": ("c", "asgpitvngetiw", "Exp"), "C 13": ("c", "asveavetgfs", "Exp"),
    "C 14": ("c", "ggnygsck", "Exp"), "C 15": ("c", "ltqeelpdpnavspc", "Exp"),
    "C 16": ("c", "sveslc", "PMPNN"),
}
TARGET_MAP = {"ADAM17": "adam17", "MMP2": "mmp2", "MMP9": "mmp9"}
ALL_AF_TARGETS = ["adam10", "adam17", "mmp2", "mmp9"]

# ESMFold2 source (optional): {(Construct, target_code): {metric: value}}
ESM_PRED = {}
if PRED_SOURCE == "esmfold2":
    _e = pd.read_csv(ARGS.esm_scores)
    _colmap = {"esm_iptm": "ipTM", "esm_ptm": "pTM", "esm_plddt": "pLDDT",
               "esm_lplddt": "LpLDDT", "esm_pae": "PAE"}
    for _, r in _e.iterrows():
        if r["Target"] not in TARGET_MAP:
            continue
        tc = TARGET_MAP[r["Target"]]
        ESM_PRED[(str(r["Construct"]).strip(), tc)] = {
            m: float(r[col]) for col, m in _colmap.items()
            if col in _e.columns and pd.notna(r[col])}
    print(f"[ESMFold2 mode] loaded {len(ESM_PRED)} construct-target rows from {ARGS.esm_scores}")

def af(metric, c, target_code):
    if PRED_SOURCE == "esmfold2":
        return ESM_PRED.get((c, target_code), {}).get(metric, np.nan)
    loop, seq, flav = CONSTRUCTS[c]
    df = (EXP_M if flav == "Exp" else PMP_M).get(metric)
    if df is None: return np.nan
    hit = df[(df["Loop"] == loop) & (df["Variant"] == seq)]
    return hit[target_code].values[0] if len(hit) else np.nan

def af_metrics_dict(c, target_code):
    return {m: af(m, c, target_code) for m in AF_METRICS}

# ----------------------------------------------------------- Experimental side
agg = pd.read_csv(AGGF)
valid = agg[agg["Trial Failed"] != True].copy()
# enrich with % of Pos Ctrl from per-folder summary_stats (disambiguate by Gated Events)
ss_rows = []
for p in glob.glob(os.path.join(ROOT, "*_Renamed_Analysis", "summary_stats.csv")):
    folder = os.path.basename(os.path.dirname(p)); m = re.match(r"([A-Za-z0-9]+)_(\d{8})", folder)
    if not m: continue
    df = pd.read_csv(p); df["Target"] = m.group(1); df["Date"] = int(m.group(2)); ss_rows.append(df)
ss = pd.concat(ss_rows, ignore_index=True).rename(columns={"Filename": "Raw Name"})
valid = valid.reset_index(drop=True); valid["_id"] = valid.index
mg = valid.merge(ss[["Target", "Date", "Raw Name", "Gated Events", "% of Pos Ctrl"]],
                 on=["Target", "Date", "Raw Name"], how="left", suffixes=("", "_ss"))
mg["_d"] = (mg["Gated Events"] - mg["Gated Events_ss"]).abs()
mg = mg.sort_values("_d").drop_duplicates("_id").sort_values("_id")
valid["% of Pos Ctrl"] = mg.set_index("_id")["% of Pos Ctrl"].reindex(valid["_id"]).values

# consensus experimental binding = mean of within-(c,t) means, z-scored across all
# points, over consistently-oriented binding readouts (a "variety of metrics",
# avoiding the sign-unstable Binding Efficiency / Stain Index).
BINDING_READOUTS = ["Pos Med Ratio", "Norm Median Ratio", "% of Pos Ctrl",
                    "Norm Intensity-Weighted Binding Index"]
g = valid.groupby(["Construct", "Target"])
mean_df = g[BINDING_READOUTS].mean()
n_df = g.size().rename("N")

CN = list(CONSTRUCTS)
TG = list(TARGET_MAP)
def cell(metric, c, T):
    return mean_df.loc[(c, T), metric] if (c, T) in mean_df.index else np.nan

# z-score each readout across all (c,t), then average -> consensus
def zmat(metric):
    M = np.array([[cell(metric, c, T) for T in TG] for c in CN], float)
    flat = M[np.isfinite(M)]
    mu, sd = flat.mean(), flat.std()
    return (M - mu) / sd
consensus = np.nanmean([zmat(m) for m in BINDING_READOUTS], axis=0)   # 13 x 3
cons_df = pd.DataFrame(consensus, index=CN, columns=TG)

# keep only fully-observed constructs for the balanced decomposition
obs = cons_df.dropna()
print(f"Consensus binding built: {obs.shape[0]} constructs x {obs.shape[1]} targets (balanced)")

# ============================================================ 1. DECOMPOSITION
X = obs.values
grand = X.mean()
row = X.mean(axis=1)          # construct / stickiness effect
col = X.mean(axis=0)          # target baseline
inter = X - row[:, None] - col[None, :] + grand
ss_constr = ((row - grand) ** 2).sum() * X.shape[1]
ss_target = ((col - grand) ** 2).sum() * X.shape[0]
ss_inter = (inter ** 2).sum()
ss_tot = ss_constr + ss_target + ss_inter
decomp = pd.DataFrame({
    "Component": ["construct (stickiness)", "target (baseline)", "interaction (target-specific)"],
    "SS": [ss_constr, ss_target, ss_inter],
    "Frac_of_total": [ss_constr / ss_tot, ss_target / ss_tot, ss_inter / ss_tot],
})
decomp.to_csv(os.path.join(OUT, "variance_decomposition.csv"), index=False)
inter_df = pd.DataFrame(inter, index=obs.index, columns=obs.columns)
stick = pd.Series(row - grand, index=obs.index, name="stickiness")

# ====================================================== 2. CALIBRATE AF METRICS
def corr(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y); x, y = x[ok], y[ok]
    if len(x) < 4 or np.ptp(x) == 0 or np.ptp(y) == 0: return (len(x), np.nan, np.nan)
    rho, p = stats.spearmanr(x, y); return (len(x), rho, p)

rows = []
OC = list(obs.index)
# (a) does AF predict the stickiness (construct) factor?  use construct-mean of AF
for amet in AF_METRICS:
    s = EXPECT_SIGN[amet]
    af_constr = [np.nanmean([af(amet, c, TARGET_MAP[T]) for T in TG]) for c in OC]
    n, rho, p = corr(af_constr, stick.values)
    rows.append(["stickiness(construct)", amet, n, (rho * s) if pd.notna(rho) else np.nan, p])
    # (b) does AF predict the target-specific interaction?  double-center AF too
    AFm = np.array([[af(amet, c, TARGET_MAP[T]) for T in TG] for c in OC], float)
    af_inter = AFm - AFm.mean(1)[:, None] - AFm.mean(0)[None, :] + AFm.mean()
    n, rho, p = corr(af_inter.flatten(), inter_df.loc[OC].values.flatten())
    rows.append(["target_specific(interaction)", amet, n, (rho * s) if pd.notna(rho) else np.nan, p])
cal = pd.DataFrame(rows, columns=["Calibration", "AF_Metric", "n", "oriented_rho", "p"])
cal.to_csv(os.path.join(OUT, f"afmetric_calibration{TAG}.csv"), index=False)

# (c) SELECTIVITY: AF on-off gap vs experimental on-off gap, per target pair
sel_rows = []
pairs = [("ADAM17", "MMP2"), ("ADAM17", "MMP9"), ("MMP2", "MMP9"),
         ("MMP2", "ADAM17"), ("MMP9", "ADAM17"), ("MMP9", "MMP2")]
for amet in AF_METRICS:
    s = EXPECT_SIGN[amet]
    for on, off in pairs:
        exp_gap = [cons_df.loc[c, on] - cons_df.loc[c, off] for c in OC]
        af_gap = [s * (af(amet, c, TARGET_MAP[on]) - af(amet, c, TARGET_MAP[off])) for c in OC]
        n, rho, p = corr(af_gap, exp_gap)
        sel_rows.append([amet, f"{on}>{off}", n, rho, p])
sel = pd.DataFrame(sel_rows, columns=["AF_Metric", "Selectivity_pair", "n", "rho", "p"])
sel.to_csv(os.path.join(OUT, f"selectivity_calibration{TAG}.csv"), index=False)

# ====================================================== 3. EVALUATE THE RECIPE
# Build per (construct,target) recipe score using off-targets = the OTHER usable
# proteases (in production use all predicted off-targets).
def recipe_score(c, on_T, weights=RECIPE_WEIGHTS):
    on = af_metrics_dict(c, TARGET_MAP[on_T])
    offs = [af_metrics_dict(c, TARGET_MAP[oT]) for oT in TG if oT != on_T]
    return score_design(on, offs, weights)

rec_rows = []
for c in OC:
    for T in TG:
        rec_rows.append({"Construct": c, "Target": T,
                         "recipe_score": recipe_score(c, T),
                         "consensus_binding": cons_df.loc[c, T],
                         "ipTM": af("ipTM", c, TARGET_MAP[T]),
                         "LpLDDT": af("LpLDDT", c, TARGET_MAP[T])})
rec = pd.DataFrame(rec_rows)
rec.to_csv(os.path.join(OUT, f"recipe_scores{TAG}.csv"), index=False)

# performance vs baselines: on-target binding (pooled within-target z) + selectivity
def pooled_within_target(score_col):
    xs, ys = [], []
    for T in TG:
        sub = rec[rec.Target == T]
        x = sub[score_col].values.astype(float); y = sub["consensus_binding"].values.astype(float)
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() > 2 and np.ptp(x[ok]) > 0:
            xs.append((x[ok] - x[ok].mean()) / (x[ok].std() or 1))
            ys.append((y[ok] - y[ok].mean()) / (y[ok].std() or 1))
    X = np.concatenate(xs); Y = np.concatenate(ys)
    return stats.spearmanr(X, Y)

val_rows = []
for col, label in [("recipe_score", "RECIPE (multi-term)"), ("ipTM", "ipTM only (current pipeline)"),
                   ("LpLDDT", "LpLDDT only")]:
    rho, p = pooled_within_target(col)
    val_rows.append({"Scorer": label, "vs": "on-target binding (pooled within-target)",
                     "Spearman_rho": rho, "p": p})
# selectivity performance: recipe gap vs experimental gap for ADAM17>MMP2 (the pair with spread)
for on, off in [("ADAM17", "MMP2"), ("ADAM17", "MMP9")]:
    exp_gap = [cons_df.loc[c, on] - cons_df.loc[c, off] for c in OC]
    rec_gap = [recipe_score(c, on) - recipe_score(c, off) for c in OC]
    rho, p = stats.spearmanr(rec_gap, exp_gap)
    val_rows.append({"Scorer": "RECIPE (multi-term)", "vs": f"selectivity {on}>{off}",
                     "Spearman_rho": rho, "p": p})
val = pd.DataFrame(val_rows)

# --- leave-one-CONSTRUCT-out cross-check: does a DATA-FITTED blend of the 4 recipe
# terms beat the reasoned fixed weights? (honest test that we are not over-tuning) ---
from binding_recipe import score_breakdown
TERMS = ["term_affinity", "term_selectivity", "term_independent", "term_expressibility"]
feat_rows = []
for c in OC:
    for T in TG:
        on = af_metrics_dict(c, TARGET_MAP[T])
        offs = [af_metrics_dict(c, TARGET_MAP[oT]) for oT in TG if oT != T]
        bd = score_breakdown(on, offs)
        feat_rows.append({"Construct": c, "Target": T, "y": cons_df.loc[c, T],
                          **{k: bd[k] for k in TERMS}})
F = pd.DataFrame(feat_rows).dropna()
# within-target z of y (rank designs within a target)
F["yz"] = F.groupby("Target")["y"].transform(lambda s: (s - s.mean()) / (s.std() or 1))
loo_pred = np.full(len(F), np.nan)
for i, c in enumerate(OC):
    tr = F[F.Construct != c]; te = F[F.Construct == c]
    if len(te) == 0: continue
    A = np.column_stack([tr[TERMS].values, np.ones(len(tr))])
    coef, *_ = np.linalg.lstsq(A, tr["yz"].values, rcond=None)
    Ate = np.column_stack([te[TERMS].values, np.ones(len(te))])
    loo_pred[F.index.get_indexer(te.index)] = Ate @ coef
ok = np.isfinite(loo_pred)
loo_rho, loo_p = stats.spearmanr(loo_pred[ok], F["yz"].values[ok])
val_rows.append({"Scorer": "DATA-FITTED blend (leave-1-construct-out)",
                 "vs": "on-target binding (pooled within-target)", "Spearman_rho": loo_rho, "p": loo_p})

# --- weight sensitivity: a few sensible variants (transparency, NOT test-set tuning) ---
VARIANTS = {
    "v1_reasoned (default)": RECIPE_WEIGHTS,
    "v2_affinity_heavy": {"affinity": 0.55, "selectivity": 0.25, "independent": 0.10, "expressibility": 0.10},
    "v3_selectivity_heavy": {"affinity": 0.30, "selectivity": 0.45, "independent": 0.10, "expressibility": 0.15},
    "v4_no_iptm": {"affinity": 0.50, "selectivity": 0.35, "independent": 0.15, "expressibility": 0.0},
}
sens_rows = []
for name, wts in VARIANTS.items():
    rec[f"_s"] = [score_design(af_metrics_dict(c, TARGET_MAP[T]),
                               [af_metrics_dict(c, TARGET_MAP[o]) for o in TG if o != T], wts)
                  for c, T in zip(rec.Construct, rec.Target)]
    rho, p = pooled_within_target("_s")
    sens_rows.append({"Variant": name, "Spearman_rho_vs_binding": rho, "p": p})
sens = pd.DataFrame(sens_rows); sens.to_csv(os.path.join(OUT, f"weight_sensitivity{TAG}.csv"), index=False)

# --- bootstrap 95% CI on the headline within-target Spearman (small-n humility) ---
def _within_target_rho(df_):
    xs, ys = [], []
    for T in TG:
        gg = df_[df_.Target == T]
        x = gg["recipe_score"].values.astype(float); y = gg["consensus_binding"].values.astype(float)
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() > 2 and np.ptp(x[ok]) > 0 and np.ptp(y[ok]) > 0:
            xs.append((x[ok] - x[ok].mean()) / x[ok].std())
            ys.append((y[ok] - y[ok].mean()) / y[ok].std())
    return stats.spearmanr(np.concatenate(xs), np.concatenate(ys))[0] if xs else np.nan

rng = np.random.default_rng(0)
boot = [r for r in (_within_target_rho(pd.concat([rec[rec.Construct == c] for c in
        rng.choice(OC, len(OC), replace=True)], ignore_index=True)) for _ in range(1000)) if np.isfinite(r)]
ci_lo, ci_hi = np.percentile(boot, [2.5, 97.5])
val_rows.append({"Scorer": "RECIPE 95% bootstrap CI (resample constructs)",
                 "vs": "on-target binding (pooled within-target)",
                 "Spearman_rho": float(np.median(boot)), "p": np.nan,
                 "ci_low": float(ci_lo), "ci_high": float(ci_hi)})
val = pd.DataFrame(val_rows)
val.to_csv(os.path.join(OUT, f"recipe_validation{TAG}.csv"), index=False)
print(f"\n[humility] RECIPE within-target rho: point {_within_target_rho(rec):+.2f}, "
      f"95% bootstrap CI [{ci_lo:+.2f}, {ci_hi:+.2f}] over n={len(OC)} constructs "
      f"(wide CI -> treat as directional, not a fitted estimate)")

# ----------------------------------------------------------------------- plots
fig, ax = plt.subplots(figsize=(7.5, 5))
xs = rec["recipe_score"].values; ys = rec["consensus_binding"].values
cols = {"ADAM17": "#dd6b20", "MMP2": "#3182ce", "MMP9": "#718096"}
for T in TG:
    sub = rec[rec.Target == T]
    ax.scatter(sub["recipe_score"], sub["consensus_binding"], label=T, c=cols[T])
    for _, r in sub.iterrows():
        ax.annotate(r["Construct"], (r["recipe_score"], r["consensus_binding"]),
                    fontsize=6, xytext=(2, 2), textcoords="offset points")
ax.set_xlabel("Recipe score (multi-term composite)")
ax.set_ylabel("Consensus experimental binding (z)")
rho, p = stats.spearmanr(xs[np.isfinite(xs) & np.isfinite(ys)], ys[np.isfinite(xs) & np.isfinite(ys)])
wt_rho, wt_p = pooled_within_target("recipe_score")
ax.set_title(f"Recipe vs experimental binding\nwithin-target ρ={wt_rho:.2f} (p={wt_p:.2f}, the ranking metric)  |  raw pooled ρ={rho:.2f}")
ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(OUT, f"plot_recipe_vs_binding{TAG}.png"), dpi=130); plt.close(fig)

# ---------------------------------------------------------------------- console
pd.set_option("display.width", 200)
print(f"\n########## CALIBRATION SOURCE: {PRED_SOURCE.upper()} "
      f"{'(AF-server == AF3, same model)' if PRED_SOURCE=='afserver' else '(outputs suffixed '+TAG+')'} ##########")
print("\n=== 1. VARIANCE DECOMPOSITION (where the experimental signal lives) ===")
print(decomp.round(3).to_string(index=False))
print("\n=== 2a/2b. AF metric calibration (oriented Spearman) ===")
print(cal.pivot(index="AF_Metric", columns="Calibration", values="oriented_rho").round(3).to_string())
print("\n=== 2c. Selectivity calibration — best AF metric per pair (oriented gap) ===")
best_sel = sel.sort_values("rho", ascending=False).groupby("Selectivity_pair").head(1)
print(best_sel.round(3).to_string(index=False))
print("\n=== 3. RECIPE vs BASELINES (+ data-fitted LOO cross-check) ===")
print(val.round(3).to_string(index=False))
print("\n=== Weight sensitivity (transparency; v1 is the reasoned default, NOT test-tuned) ===")
print(sens.round(3).to_string(index=False))
print(f"\nWrote calibration outputs to Local/Calibration/")
