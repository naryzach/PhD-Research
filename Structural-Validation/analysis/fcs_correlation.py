"""
Structure -> binding correlation layer.

Joins the ~40 in-silico structural metrics (master_complex_metrics.csv, per
construct x target x source) to the lab's aggregate flow-cytometry binding
readouts (Local/Aggregate_FCS_Analysis/aggregate_summary.csv), and asks: does
any structural signal actually track measured binding?

Framing (locked with RG 2026-07-16):
  * HEADLINE = within-target Spearman (across constructs, per target) then
    averaged across targets. This removes the target-identity confound — some
    targets simply bind more — and matches the lab's "target-specific, not raw"
    calibration finding. A pooled-across-all-pairs correlation is also computed
    and shown alongside, flagged for Simpson's-paradox risk.
  * Binding value per pair = median over valid trials (Trial Failed / Low Events /
    Low Expression dropped). Two variants: all valid trials, and purchased-target
    vendors only (a sensitivity check, since MMP3/MMP9 have in-house preps).
  * All ~10 binding readouts are correlated; the report lets the data pick the
    most-predictable readout as the headline.

Outputs (analysis/fcs/, figures/fcs/, reports/):
  binding_by_pair.csv, binding_by_pair_purchased.csv,
  structure_binding_merged.csv, correlations_long.csv, best_predictors.csv,
  figures/fcs/*.png, reports/fcs_correlation_report.html

Run:  python Structural-Validation/analysis/fcs_correlation.py
"""
from __future__ import annotations

import base64
import sys
from io import BytesIO
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C  # noqa: E402

OUT_FCS = C.OUT_ANALYSIS / "fcs"
FIG_FCS = C.OUT_FIG / "fcs"

# Binding readouts to correlate (binding-focused first). Expression-only columns
# (Expr+ %, Expr Med) are kept out of the headline but still in the long table.
READOUTS = [
    "Norm Median Ratio", "Pos Med Ratio", "Norm Bind Med (Expr+)",
    "Bind Med (Expr+)", "Binding Efficiency",
    "Norm Intensity-Weighted Binding Index", "Intensity-Weighted Binding Index",
    "Double+ %", "Norm Mean Ratio", "Pos Mean Ratio",
]
# Higher readout = more binding for all of the above.

# Structural metrics tested as predictors. Direction annotations only affect how
# we describe a correlation, not the maths.
METRICS = [
    "bsa", "f_apolar_bsa", "n_contacts_5A", "n_hbonds", "n_salt_bridges",
    "n_hydrophobic", "contact_density", "sc_shape_complementarity",
    "charge_complementarity", "n_buried_unsat_hbond", "min_ca_ca_zincloop",
    "catalytic_occlusion", "zinc_bsa_buried", "iptm", "interface_plddt",
    "interface_pae", "pdockq", "pdockq2", "lis", "dockq", "fnat", "fnonnat",
    "complex_tm", "iface_composite", "haddock_score",
]
SOURCE_ORDER = ["AF3_cofold", "ESMFold2_cofold", "HADDOCK:AF3xAF3",
                "HADDOCK:AF3xCrystal", "HADDOCK:ESMFold2xCrystal",
                "HADDOCK:ESMFold2xESMFold2"]
SRC_LABEL = {"AF3_cofold": "AF3 co-fold", "ESMFold2_cofold": "ESM co-fold",
             "HADDOCK:AF3xAF3": "HADDOCK AF3×AF3",
             "HADDOCK:AF3xCrystal": "HADDOCK AF3×Xtal",
             "HADDOCK:ESMFold2xCrystal": "HADDOCK ESM×Xtal",
             "HADDOCK:ESMFold2xESMFold2": "HADDOCK ESM×ESM"}
MIN_WITHIN = 4   # min constructs to compute a within-target correlation
MIN_TARGETS = 3  # min targets to report a within-target-then-pooled value


def _b(x):
    return str(x).strip().lower() in ("true", "1", "1.0", "yes")


def load_binding(purchased_only: bool = False) -> pd.DataFrame:
    """Per (construct_id, target_id) median binding over valid trials."""
    df = pd.read_csv(C.FCS_AGG_CSV)
    valid = df[~df["Trial Failed"].apply(_b) & ~df["Low Events"].apply(_b)
               & ~df["Low Expression"].apply(_b)].copy()
    if purchased_only:
        valid = valid[valid["Source"].isin(C.FCS_PURCHASED_VENDORS)]
    valid["construct_id"] = valid["Construct"].map(C.FCS_CONSTRUCT_MAP)
    valid["target_id"] = valid["Target"]
    valid = valid.dropna(subset=["construct_id"])
    for r in READOUTS:
        valid[r] = pd.to_numeric(valid[r], errors="coerce")
    agg = (valid.groupby(["construct_id", "target_id"])[READOUTS]
           .median().reset_index())
    agg["n_trials"] = (valid.groupby(["construct_id", "target_id"]).size()
                       .reset_index(drop=True))
    return agg


def _spear(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < MIN_WITHIN or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return np.nan, np.nan, int(m.sum())
    rho, p = spearmanr(x[m], y[m])
    return float(rho), float(p), int(m.sum())


def _auc(scores, labels):
    """AUROC of `scores` (metric) discriminating binder (label 1) from non-binder,
    via the Mann-Whitney U rank statistic. AUC>0.5: higher metric -> binder;
    <0.5: metric anti-predicts. NaN if a class is empty."""
    from scipy.stats import rankdata
    s = np.asarray(scores, float)
    y = np.asarray(labels, int)
    m = np.isfinite(s)
    s, y = s[m], y[m]
    n1, n0 = int((y == 1).sum()), int((y == 0).sum())
    if n1 == 0 or n0 == 0 or len(s) < MIN_WITHIN:
        return np.nan, n1 + n0
    r = rankdata(s)
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)), n1 + n0


def _within_target_auc(sub, metric, readout):
    """Binder = binding above the per-target median; pool across targets, one AUC.
    This is the binary-classification framing that lines up with the lab's prior
    calibration (best in-silico AUC ~0.68)."""
    scores, labels = [], []
    for _, g in sub.groupby("target_id"):
        yb = pd.to_numeric(g[readout], errors="coerce")
        med = yb.median()
        for xv, yv in zip(pd.to_numeric(g[metric], errors="coerce"), yb):
            if np.isfinite(xv) and np.isfinite(yv):
                scores.append(xv)
                labels.append(1 if yv > med else 0)
    return _auc(scores, labels)


def correlate(merged: pd.DataFrame) -> pd.DataFrame:
    """Long table: source × readout × metric × {within-target ρ, pooled ρ, AUC}."""
    rows = []
    for source in SOURCE_ORDER:
        sub = merged[merged.source == source]
        if sub.empty:
            continue
        for readout in READOUTS:
            for metric in METRICS:
                if metric not in sub or sub[metric].notna().sum() < MIN_WITHIN:
                    continue
                # within-target rho per target, then pool
                trhos, tps, tn = [], [], 0
                for t, g in sub.groupby("target_id"):
                    rho, p, n = _spear(g[metric], g[readout])
                    if not np.isnan(rho):
                        trhos.append(rho); tps.append(p); tn += 1
                within = np.mean(trhos) if len(trhos) >= MIN_TARGETS else np.nan
                n_sig = int(np.sum(np.array(tps) < 0.10)) if trhos else 0
                # pooled across all pairs
                prho, pp, pn = _spear(sub[metric], sub[readout])
                # binder/non-binder AUC (within-target median split), vs prior calibration.
                # auc is signed (dir of the metric); auc_oriented flips anti-predictors
                # to the 0.5-1 scale so it is directly comparable to the prior best ~0.68.
                auc, n_auc = _within_target_auc(sub, metric, readout)
                auc_oriented = max(auc, 1 - auc) if auc == auc else np.nan
                rows.append({
                    "source": source, "readout": readout, "metric": metric,
                    "within_target_rho": round(within, 3) if within == within else np.nan,
                    "n_targets": len(trhos), "n_targets_sig": n_sig,
                    "pooled_rho": round(prho, 3) if prho == prho else np.nan,
                    "pooled_p": round(pp, 4) if pp == pp else np.nan, "pooled_n": pn,
                    "auc": round(auc, 3) if auc == auc else np.nan,
                    "auc_oriented": round(auc_oriented, 3) if auc_oriented == auc_oriented else np.nan,
                    "auc_n": n_auc,
                })
    return pd.DataFrame(rows)


def pick_headline_readout(corr: pd.DataFrame) -> str:
    """Readout whose best within-target |rho| (any source/metric) is largest."""
    c = corr.dropna(subset=["within_target_rho"]).copy()
    c["abs"] = c["within_target_rho"].abs()
    best = c.groupby("readout")["abs"].max().sort_values(ascending=False)
    return best.index[0] if len(best) else READOUTS[0]


def build_merged(binding: pd.DataFrame) -> pd.DataFrame:
    cx = pd.read_csv(C.OUT_ANALYSIS / "master_complex_metrics.csv")
    cx = cx[cx.status == "ok"].copy()
    for m in METRICS:
        if m in cx:
            cx[m] = pd.to_numeric(cx[m], errors="coerce")
    return cx.merge(binding, on=["construct_id", "target_id"], how="inner")


def main():
    OUT_FCS.mkdir(parents=True, exist_ok=True)
    FIG_FCS.mkdir(parents=True, exist_ok=True)

    binding = load_binding(False)
    binding_p = load_binding(True)
    binding.to_csv(OUT_FCS / "binding_by_pair.csv", index=False)
    binding_p.to_csv(OUT_FCS / "binding_by_pair_purchased.csv", index=False)

    merged = build_merged(binding)
    merged_p = build_merged(binding_p)
    merged.to_csv(OUT_FCS / "structure_binding_merged.csv", index=False)

    corr = correlate(merged)
    corr_p = correlate(merged_p)
    corr["variant"], corr_p["variant"] = "all_valid", "purchased_only"
    allcorr = pd.concat([corr, corr_p], ignore_index=True)
    allcorr.to_csv(OUT_FCS / "correlations_long.csv", index=False)

    headline = pick_headline_readout(corr)

    # best predictors for the headline readout (all-valid), ranked by |within rho|
    best = (corr[corr.readout == headline].dropna(subset=["within_target_rho"])
            .assign(abs_rho=lambda d: d.within_target_rho.abs())
            .sort_values("abs_rho", ascending=False))
    best.drop(columns="abs_rho").to_csv(OUT_FCS / "best_predictors.csv", index=False)

    from fcs_report import build_figures_and_report  # noqa: E402
    build_figures_and_report(binding, merged, merged_p, corr, corr_p, headline)

    # best binder/non-binder AUC across all readouts (direction-agnostic, vs prior ~0.68)
    best_auc = (allcorr[allcorr.variant == "all_valid"]
                .dropna(subset=["auc_oriented"])
                .sort_values("auc_oriented", ascending=False))

    print(f"FCS binding pairs: {len(binding)} (all) / {len(binding_p)} (purchased)")
    print(f"Merged rows: {len(merged)}  |  headline readout: {headline!r}")
    print("\nTop 12 predictors (within-target rho vs headline readout):")
    cols = ["source", "metric", "within_target_rho", "n_targets", "n_targets_sig",
            "pooled_rho", "pooled_p", "auc"]
    print(best[cols].head(12).to_string(index=False))
    print("\nBest binder/non-binder AUC, direction-agnostic (any readout; prior best ~0.68):")
    acols = ["source", "metric", "readout", "auc", "auc_oriented", "auc_n"]
    print(best_auc[acols].head(8).to_string(index=False))
    print(f"\n-> {OUT_FCS.relative_to(C.REPO_ROOT)}")


if __name__ == "__main__":
    main()
