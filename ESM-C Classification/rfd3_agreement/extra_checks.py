"""Two extra views on ESM-C vs the RFd3 pipeline (run after score_designs + analyze_agreement).

1. Tier enrichment: are the manufacturing candidates (AF3-validated / shortlist / ordered / hall-of-fame)
   given higher ESM-C binding probability than ordinary pool designs, and than the pipeline's WEAKEST
   designs (bottom sv_pdockq quartile)?
2. Selectivity: does P(ADAM17) - P(MMP9) separate ADAM17-designed from MMP9-designed binders?

    python extra_checks.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
LOCAL = HERE.parents[1] / "Local"
RUN = LOCAL / "esmc_rfd3_agreement"
OUT = RUN / "analysis"


def boot_auc(y, p, n=500, seed=0):
    rng = np.random.default_rng(seed)
    y, p = np.asarray(y), np.asarray(p)
    a = roc_auc_score(y, p)
    bs = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(set(y[i])) == 2:
            bs.append(roc_auc_score(y[i], p[i]))
    return a, *np.percentile(bs, [2.5, 97.5])


def tiers():
    d = pd.read_csv(OUT / "design_scores_long.csv")
    d = d[d.loops_ok & d.on_window & (d["head"] == d["target"]) & d.target.isin(["ADAM17", "MMP9"])].copy()
    q = d.groupby(["variant", "head", "window", "target"])["sv_pdockq"].transform(lambda s: s.quantile(.25))
    # tier per row (a design can be in several sets; pool rows are split by pdockq quartile)
    def tier_rows(g):
        out = {}
        s = g.sets.str.split(",")
        for name in ("af3", "shortlist", "order", "hof"):
            out[name] = g[s.apply(lambda x: name in x)]
        pool = g[s.apply(lambda x: "pool" in x)]
        out["pool_bottom_q"] = pool[pool.sv_pdockq <= pool.sv_pdockq.quantile(.25)]
        out["pool_top_q"] = pool[pool.sv_pdockq >= pool.sv_pdockq.quantile(.75)]
        return out
    thr_cache, rows = {}, []
    for (v, h, w, t), g in d.groupby(["variant", "head", "window", "target"]):
        thr = thr_cache.setdefault(v, json.loads((LOCAL / "esmc_multirun" / v / "model" / "model_meta.json").read_text())["thresholds"])[h]
        tr = tier_rows(g)
        base = tr["pool_bottom_q"]
        for name, s in tr.items():
            if len(s) < 5:
                continue
            row = {"variant": v, "head": h, "window": w, "tier": name, "n": len(s),
                   "median_prob": s.prob.median(), "call_rate": (s.prob >= thr).mean(), "threshold": thr}
            if name not in ("pool_bottom_q",) and len(base) >= 20:
                y = np.r_[np.ones(len(s)), np.zeros(len(base))]
                a, lo, hi = boot_auc(y, np.r_[s.prob.values, base.prob.values])
                row.update({"auc_vs_pool_bottom_q": a, "auc_lo": lo, "auc_hi": hi})
            rows.append(row)
    return pd.DataFrame(rows)


def selectivity():
    meta = pd.read_csv(RUN / "design_meta.csv")
    seqs = pd.read_csv(RUN / "seqs.csv")[["seq_id", "full_seq"]]
    meta = meta.merge(seqs, on="full_seq")
    meta = meta[meta.target.isin(["ADAM17", "MMP9"]) & meta.sets.str.contains("pool")]
    rows = []
    for v in ("cloop_only", "all3_original", "everything_combined"):
        sc = pd.read_csv(RUN / "scores" / f"{v}.csv")
        for w in ("C", "AB"):
            s = sc[sc.window == w][["seq_id", "prob_ADAM17", "prob_MMP9"]].dropna()
            if s.empty:
                continue
            g = meta.merge(s, on="seq_id")
            g = g.assign(sel=g.prob_ADAM17 - g.prob_MMP9, is_a17=(g.target == "ADAM17").astype(int))
            a, lo, hi = boot_auc(g.is_a17, g.sel)
            rows.append({"variant": v, "window": w, "n_ADAM17_designs": int(g.is_a17.sum()),
                         "n_MMP9_designs": int((1 - g.is_a17).sum()),
                         "auc_ADAM17_vs_MMP9_designs": a, "lo": lo, "hi": hi,
                         "mean_sel_on_A17_designs": g.loc[g.is_a17 == 1, "sel"].mean(),
                         "mean_sel_on_M9_designs": g.loc[g.is_a17 == 0, "sel"].mean()})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    pd.set_option("display.width", 220)
    t = tiers(); t.to_csv(OUT / "tier_enrichment.csv", index=False)
    print("Tier enrichment (matched head, on-window):")
    print(t.round(3).to_string(index=False))
    s = selectivity(); s.to_csv(OUT / "selectivity_auc.csv", index=False)
    print("\nSelectivity (0.5 = no discrimination):")
    print(s.round(3).to_string(index=False))
