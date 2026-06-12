#!/usr/bin/env python3
"""
Metric regression sweep: every AlphaFold predicted metric vs every experimental
FCS readout, per target and pooled, to learn WHICH predicted metric best tracks
WHICH experimental quantity. Informs the metric blend used for selection and the
later ground-truth calibration.

Predicted (AlphaFold, per variant x target, Exp + PMPNN matrices):
    ipTM, LpLDDT, PAE, ApTM, BpTM, pTM, pLDDT
    (PAE is "lower = better"; all others "higher = better" — captured in EXPECT_SIGN)

Experimental (per construct x target, mean over valid trials):
    - from aggregate_summary.csv: 16 metrics (binding raw / binding normalized /
      expression / count)
    - enriched from per-folder summary_stats.csv: Bind Stain Index, Bind FC vs Pos
      Ctrl, % of Pos Ctrl, Bind MFI  (joined on Target+Date+Raw Name, disambiguated
      by Gated Events when a date has multiple source/gate folders)

Usable targets: ADAM17, MMP2, MMP9 (ADAM10 has no valid test-construct data; MMP3
was not an AlphaFold target). MMP9 has no significant between-construct spread
(ANOVA p=0.084) so it is reported but down-weighted in the headline ranking.

Outputs -> Local/Prediction_vs_Result_Analysis/Regression_Sweep/
"""
import os, glob, re
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Local")  # repo/Local (gitignored outputs)
AFM  = os.path.join(ROOT, "AlphaFoldMetrics")
AGGF = os.path.join(ROOT, "Aggregate_FCS_Analysis", "aggregate_summary.csv")
OUT  = os.path.join(ROOT, "Prediction_vs_Result_Analysis", "Regression_Sweep")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- AlphaFold side
AF_METRICS = ["ipTM", "LpLDDT", "PAE", "ApTM", "BpTM", "pTM", "pLDDT"]
EXPECT_SIGN = {"ipTM": +1, "LpLDDT": +1, "PAE": -1, "ApTM": +1, "BpTM": +1, "pTM": +1, "pLDDT": +1}

def load_matrices(flavor):
    d = {}
    for p in glob.glob(os.path.join(AFM, f"Comparison_Matrix_*_{flavor}.xlsx")):
        name = os.path.basename(p).split("Comparison_Matrix_")[1].split(f"_{flavor}.xlsx")[0]
        df = pd.read_excel(p); df["Loop"] = df["Loop"].ffill()
        df["Variant"] = df["Variant"].astype(str).str.lower()
        d[name] = df
    return d
EXP_M, PMP_M = load_matrices("Exp"), load_matrices("PMPNN")

CONSTRUCTS = {  # construct -> (loop, sequence, AF source flavor)
    "AB 1": ("ab", "dgptge", "PMPNN"), "AB 2": ("ab", "eversghkvke", "Exp"),
    "AB 3": ("ab", "kgpyge", "PMPNN"), "AB 4": ("ab", "knpdgtlt", "Exp"),
    "AB 5": ("ab", "patptstrgaggee", "Exp"), "AB 6": ("ab", "tdtfptanwtgev", "Exp"),
    "AB 7": ("ab", "tlpdgske", "Exp"),
    "C 11": ("c", "anpeyc", "PMPNN"), "C 12": ("c", "asgpitvngetiw", "Exp"),
    "C 13": ("c", "asveavetgfs", "Exp"), "C 14": ("c", "ggnygsck", "Exp"),
    "C 15": ("c", "ltqeelpdpnavspc", "Exp"), "C 16": ("c", "sveslc", "PMPNN"),
}
TARGET_MAP = {"ADAM17": "adam17", "MMP2": "mmp2", "MMP9": "mmp9"}

def af_lookup(metric, loop, seq, flavor, target):
    df = (EXP_M if flavor == "Exp" else PMP_M).get(metric)
    if df is None: return np.nan
    hit = df[(df["Loop"] == loop) & (df["Variant"] == seq)]
    return hit[target].values[0] if len(hit) else np.nan

# ------------------------------------------------------------- Experimental side
agg = pd.read_csv(AGGF)
valid = agg[agg["Trial Failed"] != True].copy()

AGG_METRICS = {
    "Pos Med Ratio": "binding_raw", "Pos Mean Ratio": "binding_raw",
    "Bind Med (Expr+)": "binding_raw", "Bind Mean (Expr+)": "binding_raw",
    "Intensity-Weighted Binding Index": "binding_raw", "Binding Efficiency": "binding_raw",
    "Double+ %": "binding_raw",
    "Norm Median Ratio": "binding_norm", "Norm Mean Ratio": "binding_norm",
    "Norm Bind Med (Expr+)": "binding_norm", "Norm Bind Mean (Expr+)": "binding_norm",
    "Norm Intensity-Weighted Binding Index": "binding_norm",
    "Expr+ %": "expression", "Expr Med (Bind+)": "expression", "Norm Expr Med (Bind+)": "expression",
    "Gated Events": "count",
}
# ---- enrich with richer binding metrics from per-folder summary_stats.csv ----
EXTRA = {"Bind Stain Index": "binding_quality", "Bind FC vs Pos Ctrl": "binding_norm",
         "% of Pos Ctrl": "binding_norm", "Bind MFI": "binding_raw"}
ss_rows = []
for p in glob.glob(os.path.join(ROOT, "*_Analysis", "summary_stats.csv")):
    folder = os.path.basename(os.path.dirname(p))
    m = re.match(r"([A-Za-z0-9]+)_(\d{8})", folder)
    if m:
        df = pd.read_csv(p)
        df["Target"] = m.group(1)
        df["Date"] = int(m.group(2))
        for col in EXTRA:
            if col not in df.columns:
                df[col] = np.nan
        ss_rows.append(df)
    else:
        m2 = re.match(r"(\d{8})_(\d{6})_Analysis", folder)
        if m2:
            df = pd.read_csv(p)
            df["Date"] = int(m2.group(1))
            for col in EXTRA:
                if col not in df.columns:
                    df[col] = np.nan
            ss_rows.append(df)
ss = pd.concat(ss_rows, ignore_index=True).rename(columns={"Filename": "Raw Name"})
ss_keep = ss[["Target", "Date", "Raw Name", "Gated Events"] + list(EXTRA)].copy()

valid = valid.reset_index(drop=True); valid["_id"] = valid.index
merged = valid.merge(ss_keep, on=["Target", "Date", "Raw Name"], how="left",
                     suffixes=("", "_ss"))
# disambiguate multi-folder dates: keep the summary_stats row whose Gated Events matches
merged["_ge_diff"] = (merged["Gated Events"] - merged["Gated Events_ss"]).abs()
merged = merged.sort_values("_ge_diff").drop_duplicates("_id", keep="first").sort_values("_id")
for col in EXTRA:
    valid[col] = merged.set_index("_id")[col].reindex(valid["_id"]).values

ALL_METRICS = {**AGG_METRICS, **EXTRA}
g = valid.groupby(["Construct", "Target"])
mean_df = g[list(ALL_METRICS)].mean()
n_ser = g.size().rename("N")

# ------------------------------------------------------------------ correlations
def vec_af(metric, target_code):
    return np.array([af_lookup(metric, *CONSTRUCTS[c][:3], target_code) for c in CONSTRUCTS], float)
def vec_exp(metric, T):
    return np.array([mean_df.loc[(c, T), metric] if (c, T) in mean_df.index else np.nan
                     for c in CONSTRUCTS], float)

def corr(xs, ys):
    ok = np.isfinite(xs) & np.isfinite(ys)
    xs, ys = xs[ok], ys[ok]
    if len(xs) < 4 or np.ptp(xs) == 0 or np.ptp(ys) == 0:
        return (len(xs), np.nan, np.nan, np.nan, np.nan)
    r, pr = stats.pearsonr(xs, ys); rho, ps = stats.spearmanr(xs, ys)
    return (len(xs), r, pr, rho, ps)

def zwithin(v):
    s = np.nanstd(v)
    return (v - np.nanmean(v)) / s if s > 0 else v * np.nan

rows = []
for amet in AF_METRICS:
    sign = EXPECT_SIGN[amet]
    for emet, ecat in ALL_METRICS.items():
        per_t_rho = {}
        zx_all, zy_all = [], []
        for T, tc in TARGET_MAP.items():
            xs, ys = vec_af(amet, tc), vec_exp(emet, T)
            n, r, pr, rho, ps = corr(xs, ys)
            per_t_rho[T] = rho
            rows.append({"AF_Metric": amet, "Exp_Metric": emet, "Exp_Category": ecat,
                         "Scope": T, "n": n, "Pearson_r": r, "Pearson_p": pr,
                         "Spearman_rho": rho, "Spearman_p": ps,
                         "Oriented_rho": rho * sign if pd.notna(rho) else np.nan})
            ok = np.isfinite(xs) & np.isfinite(ys)
            zx_all.append(zwithin(np.where(ok, xs, np.nan)))
            zy_all.append(zwithin(np.where(ok, ys, np.nan)))
        zx, zy = np.concatenate(zx_all), np.concatenate(zy_all)
        n, r, pr, rho, ps = corr(zx, zy)
        rows.append({"AF_Metric": amet, "Exp_Metric": emet, "Exp_Category": ecat,
                     "Scope": "pooled_within_target_z", "n": n, "Pearson_r": r,
                     "Pearson_p": pr, "Spearman_rho": rho, "Spearman_p": ps,
                     "Oriented_rho": rho * sign if pd.notna(rho) else np.nan})
        # consistency across ADAM17/MMP2 (targets with real spread)
        rows[-1]["Consistency_A17_M2"] = int(np.nansum([
            1 for T in ("ADAM17", "MMP2")
            if pd.notna(per_t_rho[T]) and np.sign(per_t_rho[T] * sign) > 0]))
sweep = pd.DataFrame(rows)
sweep.to_csv(os.path.join(OUT, "sweep_long.csv"), index=False)

# ---- headline ranking: which AF metric best predicts BINDING -----------------
binding_cats = {"binding_raw", "binding_norm", "binding_quality"}
bind = sweep[(sweep["Exp_Category"].isin(binding_cats))].copy()

# primary score = mean oriented Spearman over ADAM17 & MMP2 (real-spread targets)
prim = bind[bind["Scope"].isin(["ADAM17", "MMP2"])]
af_rank = (prim.groupby("AF_Metric")["Oriented_rho"]
           .agg(mean_oriented_rho="mean", median_oriented_rho="median")
           .sort_values("mean_oriented_rho", ascending=False))
# best single exp readout for each AF metric (by pooled_within_target_z oriented rho)
pz = sweep[sweep["Scope"] == "pooled_within_target_z"]
best_exp = (pz[pz["Exp_Category"].isin(binding_cats)]
            .sort_values("Oriented_rho", ascending=False)
            .groupby("AF_Metric").head(1).set_index("AF_Metric"))
af_rank["pooledz_oriented_rho"] = pz.groupby("AF_Metric")["Oriented_rho"].mean()
af_rank["best_exp_metric"] = best_exp["Exp_Metric"]
af_rank["best_exp_rho"] = best_exp["Oriented_rho"]
af_rank.to_csv(os.path.join(OUT, "af_predictor_ranking.csv"))

# ---- which experimental readout is most predictable (ground-truth candidate) --
exp_rank = (pz.assign(absrho=pz["Oriented_rho"].abs())
            .sort_values("Oriented_rho", ascending=False)
            .groupby(["Exp_Metric", "Exp_Category"]).head(1)
            .loc[:, ["Exp_Metric", "Exp_Category", "AF_Metric", "Oriented_rho", "n"]]
            .rename(columns={"AF_Metric": "best_AF_metric", "Oriented_rho": "best_oriented_rho"})
            .sort_values("best_oriented_rho", ascending=False))
exp_rank.to_csv(os.path.join(OUT, "exp_readout_predictability.csv"), index=False)

# ---- AF-metric redundancy (pooled across all construct x target points) -------
afmat = {}
for amet in AF_METRICS:
    col = []
    for T, tc in TARGET_MAP.items():
        col.append(vec_af(amet, tc))
    afmat[amet] = np.concatenate(col)
afdf = pd.DataFrame(afmat)
af_inter = afdf.corr(method="spearman")
af_inter.to_csv(os.path.join(OUT, "af_metric_intercorrelation.csv"))

# ---- DYNAMIC RANGE diagnostic: a correlation is meaningless if the metric is ---
# flat. Report within-target spread of each AF metric across the constructs.
dr_rows = []
for amet in AF_METRICS:
    stds, rngs = [], []
    for T, tc in TARGET_MAP.items():
        v = vec_af(amet, tc); v = v[np.isfinite(v)]
        if len(v) > 1:
            stds.append(np.std(v)); rngs.append(np.ptp(v))
    dr_rows.append({"AF_Metric": amet, "within_target_std": np.mean(stds),
                    "within_target_range": np.mean(rngs),
                    "cross_construct_min": np.nanmin([np.nanmin(vec_af(amet, tc)) for tc in TARGET_MAP.values()]),
                    "cross_construct_max": np.nanmax([np.nanmax(vec_af(amet, tc)) for tc in TARGET_MAP.values()])})
dyn = pd.DataFrame(dr_rows).set_index("AF_Metric")
dyn["flat_flag"] = np.where(dyn["within_target_range"] < 0.05, "FLAT (rank=noise)", "has range")
dyn.to_csv(os.path.join(OUT, "af_metric_dynamic_range.csv"))
af_rank["within_target_range"] = dyn["within_target_range"]
af_rank["flat_flag"] = dyn["flat_flag"]
af_rank.to_csv(os.path.join(OUT, "af_predictor_ranking.csv"))

# ---- EXPRESSION-confound view: does the AF metric track binding or expression? -
conf = (pz[pz["Exp_Metric"].isin(["Expr+ %", "% of Pos Ctrl", "Pos Med Ratio", "Norm Median Ratio"])]
        .pivot(index="AF_Metric", columns="Exp_Metric", values="Oriented_rho")
        .reindex(AF_METRICS).round(3))
conf.to_csv(os.path.join(OUT, "expression_vs_binding_confound.csv"))

# ------------------------------------------------------------------------ heatmap
def heatmap(scope, fname, title):
    sub = sweep[sweep["Scope"] == scope]
    piv = sub.pivot(index="AF_Metric", columns="Exp_Metric", values="Oriented_rho")
    piv = piv.reindex(index=AF_METRICS)
    order = [m for m in ALL_METRICS]  # category-grouped order
    piv = piv.reindex(columns=order)
    fig, ax = plt.subplots(figsize=(15, 4.6))
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-0.7, vmax=0.7, aspect="auto")
    ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns, rotation=60, ha="right", fontsize=7)
    ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index, fontsize=9)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if pd.notna(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6,
                        color="white" if abs(v) > 0.45 else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.018, pad=0.01); cb.set_label("Oriented Spearman ρ\n(+ = predicts as expected)")
    ax.set_title(title, fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, fname), dpi=130); plt.close(fig)

heatmap("pooled_within_target_z", "heatmap_pooledz.png",
        "AF metric (rows) vs experimental readout (cols) — pooled within-target z, oriented Spearman ρ")
heatmap("ADAM17", "heatmap_ADAM17.png", "ADAM17 only — oriented Spearman ρ (AF metric vs experimental readout)")
heatmap("MMP2", "heatmap_MMP2.png", "MMP2 only — oriented Spearman ρ (AF metric vs experimental readout)")
heatmap("MMP9", "heatmap_MMP9.png",
        "MMP9 only — oriented Spearman ρ (AF metric vs experimental readout)  [NOTE: MMP9 ANOVA n.s. — no real spread]")

# ---- per-pair scatter: every AF metric x every experimental readout, 3 panels --
# (generalizes the original plot_LpLDDT_vs_binding.png to all pairs)
def _san(s):
    out = (s.replace("%", "pct").replace("+", "pos").replace("/", "_")
            .replace("(", "").replace(")", "").replace(" ", "_").replace("-", "_"))
    while "__" in out:
        out = out.replace("__", "_")
    return out.strip("_")

PS_DIR = os.path.join(OUT, "Pairwise_Scatter")
TARGET_COLORS = {"ADAM17": "#dd6b20", "MMP2": "#3182ce", "MMP9": "#718096"}
n_made = 0
for amet in AF_METRICS:
    sign = EXPECT_SIGN[amet]
    for emet, ecat in ALL_METRICS.items():
        d = os.path.join(PS_DIR, ecat); os.makedirs(d, exist_ok=True)
        fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
        for ax, (T, tc) in zip(axes, TARGET_MAP.items()):
            xs, ys = vec_af(amet, tc), vec_exp(emet, T)
            labs = list(CONSTRUCTS)
            ok = np.isfinite(xs) & np.isfinite(ys)
            x, y = xs[ok], ys[ok]
            lb = [l for l, o in zip(labs, ok) if o]
            ax.scatter(x, y, c=TARGET_COLORS[T])
            for xi, yi, li in zip(x, y, lb):
                ax.annotate(li, (xi, yi), fontsize=6, xytext=(2, 2), textcoords="offset points")
            if len(x) >= 4 and np.ptp(x) > 0 and np.ptp(y) > 0:
                r, pr = stats.pearsonr(x, y); rho, ps = stats.spearmanr(x, y)
                m, b = np.polyfit(x, y, 1); xx = np.linspace(x.min(), x.max(), 50)
                ax.plot(xx, m * xx + b, "--", c="gray", lw=1)
                ax.set_title(f"{T}\nPearson r={r:.2f} (p={pr:.2f})  Spearman ρ={rho:.2f}", fontsize=8)
            else:
                ax.set_title(f"{T}\n(insufficient / flat — n={len(x)})", fontsize=8)
            ax.set_xlabel(amet + ("  (lower = better)" if sign < 0 else ""), fontsize=8)
            ax.set_ylabel(emet, fontsize=8)
        fig.suptitle(f"{amet}  vs  {emet}    [{ecat}]", fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(d, f"{_san(emet)}__{amet}.png"), dpi=105)
        plt.close(fig)
        n_made += 1
print(f"Pairwise_Scatter: wrote {n_made} figures across {len(set(ALL_METRICS.values()))} category subfolders")

# ------------------------------------------------------------------------- print
pd.set_option("display.width", 200)
print("=== AF predictor ranking (primary = mean oriented Spearman over ADAM17 & MMP2 binding readouts) ===")
print(af_rank.round(3).to_string())
print("\n=== Most predictable experimental readouts (best AF predictor, pooled within-target z) ===")
print(exp_rank.head(12).round(3).to_string(index=False))
print("\n=== AF-metric redundancy (Spearman) ===")
print(af_inter.round(2).to_string())
print("\n=== AF-metric DYNAMIC RANGE within a target (correlation is noise if flat) ===")
print(dyn.round(3).to_string())
print("\n=== Expression vs binding confound (oriented Spearman, pooled within-target z) ===")
print(conf.to_string())
print(f"\nWrote outputs to {os.path.relpath(OUT, ROOT+'/..')}")
