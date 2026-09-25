"""Do the trained ESM-C models agree with the RFd3-pipeline structural predictions?

Reads Local/esmc_rfd3_agreement/{design_meta,seqs,controls_ids}.csv and scores/<variant>.csv.
Writes tables + figures to Local/esmc_rfd3_agreement/analysis/.

Comparisons are (variant, window, head):
  window = the loop the head was trained on (ADAM17 -> C, MMP3 -> AB, MMP9 -> the variant's
  own loop; the mixed models' MMP9 head saw both, so both windows are reported).
  EF is scored as a null window (no model saw EF variation).

    python analyze_agreement.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
LOCAL = HERE.parents[1] / "Local"
RUN = LOCAL / "esmc_rfd3_agreement"
OUT = RUN / "analysis"
OUT.mkdir(exist_ok=True)

# variant -> {head: [windows the head was trained on]}
TRAINED = {
    "abloop_only":         {"MMP3": ["AB"], "MMP9": ["AB"]},
    "cloop_only":          {"ADAM17": ["C"], "MMP9": ["C"]},
    "mmp9_other":          {"MMP9": ["C"]},
    "all3_original":       {"ADAM17": ["C"], "MMP3": ["AB"], "MMP9": ["AB", "C"]},
    "everything_combined": {"ADAM17": ["C"], "MMP3": ["AB"], "MMP9": ["AB", "C"]},
}
STRUCT = ["sv_pdockq", "af3_iptm", "esm_iptm", "esm_plddt"]
BLUE, GRAY, INK, INK2, GRID = "#2a78d6", "#9a9a94", "#0b0b0b", "#52514e", "#e3e2dd"
RNG = np.random.default_rng(0)

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2,
    "ytick.color": INK2, "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "figure.dpi": 130,
    "axes.axisbelow": True,
})


def spearman_ci(x, y, n_boot=1000):
    m = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[m], np.asarray(y)[m]
    if len(x) < 8 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan, np.nan, np.nan, int(m.sum())
    r = stats.spearmanr(x, y)[0]
    bs = []
    for _ in range(n_boot):
        i = RNG.integers(0, len(x), len(x))
        bs.append(stats.spearmanr(x[i], y[i])[0])
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    return r, lo, hi, int(m.sum())


def partial_spearman(x, y, Z):
    """Spearman of x,y after regressing out the ranks of the covariate columns in Z."""
    Z = np.asarray(Z, float).reshape(len(x), -1)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(Z).all(1)
    if m.sum() < 8 + Z.shape[1]:
        return np.nan
    rx, ry = stats.rankdata(x[m]), stats.rankdata(y[m])
    A = np.c_[np.ones(m.sum()), np.column_stack([stats.rankdata(Z[m, k]) for k in range(Z.shape[1])])]
    resid = lambda a: a - A @ np.linalg.lstsq(A, a, rcond=None)[0]
    return float(np.corrcoef(resid(rx), resid(ry))[0, 1])


def load():
    meta = pd.read_csv(RUN / "design_meta.csv")
    seqs = pd.read_csv(RUN / "seqs.csv")[["seq_id", "full_seq"]]
    meta = meta.merge(seqs, on="full_seq", how="left")
    # drop designs whose flank-located loop lengths disagree with the pipeline's own record
    for w, col in (("AB", "loop_AB_seq"), ("C", "loop_C_seq")):
        has = meta[col].notna()
        meta.loc[has & (meta[f"len_{w}"] != meta[col].str.len()), "loops_ok"] = False
    meta["loops_ok"] = meta["loops_ok"].fillna(True).astype(bool)
    scores = {}
    for f in sorted((RUN / "scores").glob("*.csv")):
        scores[f.stem] = pd.read_csv(f)
    ctrl = pd.read_csv(RUN / "controls_ids.csv")
    return meta, scores, ctrl


def controls_table(scores, ctrl):
    """Held-out lab rows through this scoring path (sanity: should match test_report AUCs)."""
    rows = []
    for v, sc in scores.items():
        c = ctrl[ctrl.variant == v]
        if c.empty:
            continue
        w = np.where(c.loop_start.values == 30, "AB", "C")
        c = c.assign(window=w).merge(sc, on=["seq_id", "window"])
        for t in c.target.unique():
            g = c[c.target == t]
            if g.label.nunique() == 2:
                rows.append({"variant": v, "head": t, "n": len(g), "auc": roc_auc_score(g.label, g[f"prob_{t}"])})
    return pd.DataFrame(rows)


def design_table(meta, scores):
    """Long: one row per (variant, window, head, design) with prob + structural metrics."""
    parts = []
    for v, sc in scores.items():
        for head, wins in TRAINED[v].items():
            for w in wins + ["EF"]:
                s = sc[sc.window == w][["seq_id", f"prob_{head}"]].rename(columns={f"prob_{head}": "prob"})
                d = meta.merge(s, on="seq_id")
                d = d.assign(variant=v, head=head, window=w, on_window=(w != "EF"))
                parts.append(d)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def agreement(d):
    rows = []
    for (v, head, w), g in d.groupby(["variant", "head", "window"]):
        for tgt in ("ADAM17", "MMP9", "MMP2", "ADAM10"):
            gg = g[(g.target == tgt) & g.loops_ok]
            for m in STRUCT:
                if gg[m].notna().sum() < 8:
                    continue
                r, lo, hi, n = spearman_ci(gg["prob"].values, gg[m].values)
                # Structure metrics track loop lengths and Cys count, and so does the ESM-C
                # output; partial those out (ranks) to see whether any loop-specific signal is left.
                Z = gg[["len_AB", "len_C", "len_EF", "n_cys"]].values.astype(float)
                # Native-length slice: only designs whose loop in this window is the 6 aa the
                # models were trained on (in-distribution for the pooled window).
                nat = gg[gg[f"len_{w}"] == 6] if w != "EF" else gg.iloc[0:0]
                rn, lon, hin, nn = spearman_ci(nat["prob"].values, nat[m].values)
                rows.append({"variant": v, "head": head, "window": w, "design_target": tgt, "metric": m,
                             "n": n, "rho": r, "lo": lo, "hi": hi,
                             "rho_partial": partial_spearman(gg["prob"].values, gg[m].values, Z),
                             "n_native6": nn, "rho_native6": rn, "lo_native6": lon, "hi_native6": hin,
                             "matched": head == tgt})
    return pd.DataFrame(rows)


def quartile_auc(d):
    """AUC of ESM-C prob for top-vs-bottom sv_pdockq quartile within each design target."""
    rows = []
    for (v, head, w), g in d.groupby(["variant", "head", "window"]):
        for tgt in ("ADAM17", "MMP9", "MMP2", "ADAM10"):
            gg = g[(g.target == tgt) & g.loops_ok & g.sv_pdockq.notna()]
            if len(gg) < 40:
                continue
            q1, q3 = gg.sv_pdockq.quantile([.25, .75])
            hi, lo = gg[gg.sv_pdockq >= q3], gg[gg.sv_pdockq <= q1]
            y = np.r_[np.ones(len(hi)), np.zeros(len(lo))]
            p = np.r_[hi.prob.values, lo.prob.values]
            rows.append({"variant": v, "head": head, "window": w, "design_target": tgt,
                         "n_hi": len(hi), "n_lo": len(lo), "auc_top_vs_bottom_pdockq": roc_auc_score(y, p)})
    return pd.DataFrame(rows)


def candidates_table(d, scores):
    """Manufacturing candidates: ESM-C call at each model's tuned threshold vs structural evidence."""
    rows = []
    for v in scores:
        thr = json.loads((LOCAL / "esmc_multirun" / v / "model" / "model_meta.json").read_text())["thresholds"]
        g = d[(d.variant == v) & d.on_window & d.sets.str.contains("shortlist|order|af3|hof") & (d["head"] == d["target"])]
        for r in g.itertuples():
            rows.append({"variant": v, "design_id": r.design_id, "target": r.target, "window": r.window,
                         "prob": r.prob, "threshold": thr[r.head], "esmc_call": r.prob >= thr[r.head],
                         "sv_pdockq": r.sv_pdockq, "af3_iptm": r.af3_iptm, "sets": r.sets, "full_seq": r.full_seq})
    return pd.DataFrame(rows)


def fig_scatter(d, variant, head, window, fname):
    """ESM-C prob vs sv_pdockq for designs against the head's own target; EF window as gray null."""
    g = d[(d.variant == variant) & (d["head"] == head) & (d.target == head) & d.loops_ok & d.sv_pdockq.notna()]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), sharey=True, sharex=True)
    for ax, (w, col, ttl) in zip(axes, [(window, BLUE, f"{window} window (trained)"), ("EF", GRAY, "EF window (null control)")]):
        s = g[g.window == w]
        ax.scatter(s.sv_pdockq, s.prob, s=9, color=col, alpha=0.55, linewidths=0)
        r = spearman_ci(s.prob.values, s.sv_pdockq.values, n_boot=200)
        ax.set_title(f"{ttl}\nSpearman {r[0]:+.2f}  (n={r[3]})", fontsize=8.5, color=INK2, loc="left")
        ax.set_xlabel("sv_pdockq (structure)")
    axes[0].set_ylabel(f"ESM-C P({head} binds)")
    fig.suptitle(f"{variant}: {head} designs", fontsize=9.5, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(OUT / fname); plt.close(fig)


def fig_forest(ag, fname):
    """Spearman with CI vs sv_pdockq, matched heads, on-window vs EF null."""
    a = ag[(ag.metric == "sv_pdockq") & ag.matched].copy()
    if a.empty:
        return
    a["label"] = a.variant + " / " + a["head"] + " · " + a.design_target
    labels = sorted(a.label.unique())
    fig, ax = plt.subplots(figsize=(6.8, 0.42 * len(labels) + 1.3))
    for i, lab in enumerate(labels):
        for w, col, dy in (("on", BLUE, -0.12), ("EF", GRAY, 0.12)):
            s = a[(a.label == lab) & ((a.window != "EF") if w == "on" else (a.window == "EF"))]
            if s.empty:
                continue
            s = s.iloc[0] if w == "EF" else s.sort_values("n").iloc[-1]
            ax.plot([s.lo, s.hi], [i + dy] * 2, color=col, lw=1.6)
            ax.plot(s.rho, i + dy, "o", color=col, ms=5, mec="white", mew=1)
    ax.axvline(0, color=INK2, lw=0.8)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Spearman ρ, ESM-C probability vs sv_pdockq (95% bootstrap CI)")
    ax.plot([], [], "o", color=BLUE, label="trained loop window"); ax.plot([], [], "o", color=GRAY, label="EF window (null)")
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    fig.tight_layout(); fig.savefig(OUT / fname); plt.close(fig)


def main():
    meta, scores, ctrl = load()
    print("variants scored:", list(scores))
    ct = controls_table(scores, ctrl)
    ct.to_csv(OUT / "controls_auc.csv", index=False)
    print("\nHeld-out lab controls through this scoring path:\n", ct.round(3).to_string(index=False))

    d = design_table(meta, scores)
    d.to_csv(OUT / "design_scores_long.csv", index=False)
    ag = agreement(d)
    ag.to_csv(OUT / "agreement_spearman.csv", index=False)
    qa = quartile_auc(d)
    qa.to_csv(OUT / "agreement_quartile_auc.csv", index=False)
    cand = candidates_table(d, scores)
    cand.to_csv(OUT / "candidates_calls.csv", index=False)

    pd.set_option("display.width", 200)
    print("\nSpearman vs sv_pdockq, matched head (head == design target):")
    print(ag[(ag.metric == "sv_pdockq") & ag.matched].round(3).drop(columns=["matched"]).to_string(index=False))
    print("\nSpearman vs af3_iptm, matched head:")
    print(ag[(ag.metric == "af3_iptm") & ag.matched].round(3).drop(columns=["matched"]).to_string(index=False))
    print("\nTop-vs-bottom pdockq quartile AUC, matched head:")
    print(qa[qa["head"] == qa.design_target].round(3).to_string(index=False))

    fig_forest(ag, "forest_spearman_pdockq.png")
    for v in scores:
        for head, wins in TRAINED[v].items():
            if head in ("ADAM17", "MMP9"):
                fig_scatter(d, v, head, wins[-1], f"scatter_{v}_{head}.png")
    print("\nWrote tables + figures to", OUT)


if __name__ == "__main__":
    main()
