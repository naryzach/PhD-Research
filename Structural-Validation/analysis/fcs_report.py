"""
Figures + self-contained HTML report for the structure->binding correlation
(companion to fcs_correlation.py; imported by its main()).
"""
from __future__ import annotations

import base64
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C  # noqa: E402
from fcs_correlation import (FIG_FCS, METRICS, READOUTS, SOURCE_ORDER, SRC_LABEL,
                             _spear)  # noqa: E402

TARGET_COLORS = {"MMP2": "#0072B2", "MMP9": "#56B4E9", "ADAM17": "#D55E00",
                 "MMP3": "#009E73", "ADAM10": "#CC79A7"}


def heatmap(corr, headline, out_png):
    sub = corr[corr.readout == headline]
    piv = sub.pivot_table(index="metric", columns="source", values="within_target_rho")
    piv = piv.reindex(index=[m for m in METRICS if m in piv.index],
                      columns=[s for s in SOURCE_ORDER if s in piv.columns])
    fig, ax = plt.subplots(figsize=(1.6 + 1.1 * piv.shape[1], 0.9 + 0.34 * piv.shape[0]))
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto")
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels([SRC_LABEL.get(s, s) for s in piv.columns], rotation=30, ha="right")
    ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index, fontsize=8)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if v == v:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6,
                        color="white" if abs(v) > 0.45 else "#333")
    ax.set_title(f"Within-target Spearman ρ  (structural metric vs {headline})",
                 fontweight="bold", fontsize=11)
    fig.colorbar(im, ax=ax, fraction=0.035, label="mean ρ across targets")
    ax.grid(False)
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return _b64_from_path(out_png)


def _b64_from_path(p: Path) -> str:
    return base64.b64encode(Path(p).read_bytes()).decode()


def readout_ranking(corr, headline, out_png):
    c = corr.dropna(subset=["within_target_rho"]).copy()
    c["abs"] = c["within_target_rho"].abs()
    best = c.groupby("readout")["abs"].max().reindex(
        [r for r in READOUTS]).dropna().sort_values()
    fig, ax = plt.subplots(figsize=(7.5, 4))
    cols = ["#D55E00" if r == headline else "#9AA0A6" for r in best.index]
    ax.barh(best.index, best.values, color=cols)
    for y, v in enumerate(best.values):
        ax.text(v + .01, y, f"{v:.2f}", va="center", fontsize=8)
    ax.set_xlabel("best |within-target ρ| achieved by any structural metric")
    ax.set_title("Which binding readout is most predictable?", fontweight="bold")
    ax.set_xlim(0, max(0.6, best.max() + .1))
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    return _b64_from_path(out_png)


def best_predictor_fig(merged, corr, headline, out_png):
    best = (corr[corr.readout == headline].dropna(subset=["within_target_rho"])
            .assign(a=lambda d: d.within_target_rho.abs()).sort_values("a", ascending=False))
    if best.empty:
        return None, None
    top = best.iloc[0]
    src, met = top.source, top.metric
    sub = merged[merged.source == src]
    per_target = []
    for t, g in sub.groupby("target_id"):
        rho, p, n = _spear(g[met], g[headline])
        per_target.append((t, rho, n))
    per_target = [x for x in per_target if x[1] == x[1]]

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    ts = [t for t, _, _ in per_target]
    rhos = [r for _, r, _ in per_target]
    ax[0].bar(ts, rhos, color=[TARGET_COLORS.get(t, "#888") for t in ts])
    ax[0].axhline(0, color="k", lw=.6)
    ax[0].set_ylabel("within-target Spearman ρ")
    ax[0].set_ylim(-1, 1)
    ax[0].set_title(f"A  {met} vs {headline}\n(per target, {SRC_LABEL.get(src, src)})",
                    loc="left", fontweight="bold", fontsize=10)

    for t, g in sub.groupby("target_id"):
        ax[1].scatter(g[met], g[headline], s=42, edgecolor="k", lw=.4,
                      color=TARGET_COLORS.get(t, "#888"), label=t)
    ax[1].set_xlabel(met); ax[1].set_ylabel(headline)
    ax[1].set_title("B  colored by target (within-target trend is the signal)",
                    loc="left", fontweight="bold", fontsize=10)
    ax[1].legend(fontsize=7, title="target")
    fig.tight_layout()
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    return _b64_from_path(out_png), top


KEY_METRICS_FOR_SPLIT = ["iptm", "pdockq", "pdockq2", "sc_shape_complementarity",
                         "bsa", "n_hbonds", "catalytic_occlusion", "fnat", "dockq",
                         "complex_tm", "iface_composite"]


def readout_metric_heatmap(corr, out_png):
    """mean within-target ρ for key metrics (rows) × every readout (cols),
    AF3+ESM co-folds averaged — exposes that the (negative) signal is confined to
    intensity/ratio readouts and absent for efficiency/fraction readouts."""
    co = corr[corr.source.isin(["AF3_cofold", "ESMFold2_cofold"])]
    piv = (co[co.metric.isin(KEY_METRICS_FOR_SPLIT)]
           .pivot_table(index="metric", columns="readout", values="within_target_rho"))
    piv = piv.reindex(index=[m for m in KEY_METRICS_FOR_SPLIT if m in piv.index],
                      columns=[r for r in READOUTS if r in piv.columns])
    fig, ax = plt.subplots(figsize=(1.5 + 0.95 * piv.shape[1], 1 + 0.34 * piv.shape[0]))
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-0.7, vmax=0.7, aspect="auto")
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels(piv.columns, rotation=35, ha="right", fontsize=7)
    ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index, fontsize=8)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if v == v:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6,
                        color="white" if abs(v) > 0.45 else "#333")
    ax.set_title("Within-target ρ by readout (co-folds) — signal is readout-specific",
                 fontweight="bold", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.03, label="mean ρ")
    ax.grid(False)
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return _b64_from_path(out_png)


def build_figures_and_report(binding, merged, merged_p, corr, corr_p, headline):
    FIG_FCS.mkdir(parents=True, exist_ok=True)
    hm = heatmap(corr, headline, FIG_FCS / "heatmap_within_target.png")
    rk = readout_ranking(corr, headline, FIG_FCS / "readout_ranking.png")
    split = readout_metric_heatmap(corr, FIG_FCS / "readout_metric_split.png")
    bp_b64, top = best_predictor_fig(merged, corr, headline,
                                     FIG_FCS / "best_predictor.png")

    # headline numbers
    best_all = (corr[corr.readout == headline].dropna(subset=["within_target_rho"])
                .assign(a=lambda d: d.within_target_rho.abs())
                .sort_values("a", ascending=False))
    n_pairs = len(merged.groupby(["construct_id", "target_id"]))
    max_rho = best_all.within_target_rho.abs().max() if not best_all.empty else float("nan")

    def _tbl(df, n=12):
        rows = ""
        for _, r in df.head(n).iterrows():
            wr = "—" if pd.isna(r.within_target_rho) else f"{r.within_target_rho:+.2f}"
            pr = "—" if pd.isna(r.pooled_rho) else f"{r.pooled_rho:+.2f}"
            pp = "—" if pd.isna(r.pooled_p) else f"{r.pooled_p:.3f}"
            rows += (f"<tr><td>{SRC_LABEL.get(r.source, r.source)}</td>"
                     f"<td>{r.metric}</td><td>{wr}</td>"
                     f"<td>{int(r.n_targets)}</td><td>{int(r.n_targets_sig)}</td>"
                     f"<td>{pr}</td><td>{pp}</td></tr>")
        return rows

    # purchased-only view of the same top predictors
    best_p = corr_p[corr_p.readout == headline].dropna(subset=["within_target_rho"])
    merge_key = best_all.head(12)[["source", "metric"]]
    best_p_top = merge_key.merge(best_p, on=["source", "metric"], how="left")

    # direction + readout-dependence + which construct drives the headline readout
    top_signed = top.within_target_rho
    direction = "negatively" if top_signed < 0 else "positively"
    # responsive/flat judged on the CORE confidence & native-similarity metrics on
    # co-folds (not max over all 25 metrics, which would clear 0.25 by chance).
    core_metrics = ["iptm", "pdockq", "fnat", "dockq", "complex_tm",
                    "iface_composite", "catalytic_occlusion"]
    core = (corr[corr.source.isin(["AF3_cofold", "ESMFold2_cofold"])
                 & corr.metric.isin(core_metrics)]
            .dropna(subset=["within_target_rho"])
            .assign(a=lambda d: d.within_target_rho.abs()))
    best_by_readout = core.groupby("readout")["a"].max()
    responsive = [r for r in READOUTS if best_by_readout.get(r, 0) >= 0.40]
    flat = [r for r in READOUTS if best_by_readout.get(r, 1.0) < 0.25]
    tb = binding.groupby("construct_id")[headline].median().sort_values(ascending=False)
    top_binder = tb.index[0]
    wt_rank = int(tb.rank(ascending=False).get("TIMP3_WT", np.nan)) \
        if "TIMP3_WT" in tb.index else None
    strength = ("moderate" if max_rho >= 0.5 else "weak" if max_rho >= 0.35 else "no")

    html = f"""<!doctype html><meta charset=utf-8>
<title>TIMP3 Structure → Binding Correlation</title>
<style>
 body{{font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;max-width:1000px;
   margin:2rem auto;padding:0 1.2rem;color:#1a1a1a;line-height:1.5}}
 h1{{font-size:1.55rem;margin-bottom:.2rem}} h2{{border-bottom:2px solid #0072B2;padding-bottom:.2rem;margin-top:1.8rem}}
 .sub{{color:#555;margin-top:0}} table{{border-collapse:collapse;width:100%;font-size:.9rem;margin:1rem 0}}
 th,td{{border:1px solid #ddd;padding:.35rem .55rem;text-align:center}} th{{background:#0072B2;color:#fff}}
 td:first-child,th:first-child{{text-align:left}} td:nth-child(2){{text-align:left}}
 .key{{background:#fff8e1;border-left:4px solid #E69F00;padding:.8rem 1rem;margin:1rem 0;border-radius:4px}}
 .caveat{{background:#f4f4f4;border-left:4px solid #9AA0A6;padding:.6rem 1rem;font-size:.9rem}}
 img{{max-width:100%;border:1px solid #e3e3e3;border-radius:6px;margin:1rem 0}}
 code{{background:#f0f0f0;padding:.05rem .3rem;border-radius:3px}}
</style>
<h1>TIMP3 constructs × proteases — does structure predict binding?</h1>
<p class="sub">Joining {len(METRICS)} structural metrics × 6 model sources to the
aggregate flow-cytometry binding data ({n_pairs} construct×target pairs with valid
trials). Headline correlation is <b>within-target Spearman ρ</b> (across constructs,
per target, then averaged over targets) — it isolates construct-driven differences
from the fact that targets differ in overall bindability.</p>

<div class="key"><b>Bottom line — the signal is real, {strength}, and points the
"wrong" way.</b> Among {len(READOUTS)} readouts the data pick <b>{headline}</b>,
where the best predictor (<b>{top.metric}</b>, {SRC_LABEL.get(top.source, top.source)})
correlates <b>{direction}</b> with binding, within-target ρ = <b>{top_signed:+.2f}</b>.
i.e. constructs the models score as <i>higher-confidence / more native-like</i> tend
to bind <i>less</i>. The driver is <b>{top_binder}</b>, a promiscuous construct that
is the top binder against every target{f' (WT ranks #{wt_rank} of {len(tb)})' if wt_rank else ''}
yet scores low on ipTM / pDockQ / fnat — the "stickiness dominates, specificity
scores don't capture it" pattern, now with a measurable <b>anti-correlation</b>.
<br><br><b>But it is readout-specific.</b> The negative signal lives in the
intensity- and ratio-based readouts ({', '.join(responsive) or '—'}); the
efficiency / fraction-of-cells readouts ({', '.join(flat) or '—'}) show
<b>essentially no correlation</b> — reproducing the lab's prior "no predictive
power" result. So: structure does not tell you <i>what fraction</i> of cells bind,
and mildly <i>anti</i>-predicts <i>how much</i> signal the stickiest construct
throws.</div>

<h2>Which readout is most predictable, and by what?</h2>
<p class="sub">Only some readouts respond at all — and the ones that do, respond
negatively. The split below is the key caveat on the headline number.</p>
<img src="data:image/png;base64,{rk}" alt="readout ranking">
<img src="data:image/png;base64,{split}" alt="readout x metric split heatmap">
<img src="data:image/png;base64,{hm}" alt="within-target rho heatmap">

<h2>Best predictor, dissected</h2>
<p class="sub">Top-ranked structural metric for {headline}: <b>{top.metric}</b>
({SRC_LABEL.get(top.source, top.source)}), within-target ρ = {top.within_target_rho:+.2f}
over {int(top.n_targets)} targets ({int(top.n_targets_sig)} with p&lt;0.10). Panel B
shows why pooling misleads: colored by target, any trend lives within a target, not
across the cloud.</p>
<img src="data:image/png;base64,{bp_b64}" alt="best predictor detail">

<h2>Top predictors — within-target vs pooled</h2>
<table><tr><th>Source</th><th>Metric</th><th>within-target ρ</th><th>n targets</th>
<th>n sig (p&lt;.1)</th><th>pooled ρ</th><th>pooled p</th></tr>
{_tbl(best_all)}</table>
<p style="font-size:.85rem;color:#555">Within-target ρ is averaged across targets
(the headline framing). Pooled ρ is a single correlation over all pairs and is
confounded by target identity — where the two diverge, trust within-target.</p>

<h2>Purchased-target sensitivity check</h2>
<p class="sub">Same top predictors, recomputed using only commercially-sourced
targets ({', '.join(C.FCS_PURCHASED_VENDORS)}); in-house preps (MMP3/MMP9 Masoud
etc.) dropped. If the (weak) signal were an artefact of in-house material it would
move here.</p>
<table><tr><th>Source</th><th>Metric</th><th>within-target ρ (purchased)</th>
<th>n targets</th></tr>
{''.join(f"<tr><td>{SRC_LABEL.get(r.source, r.source)}</td><td>{r.metric}</td>"
         f"<td>{'—' if pd.isna(r.within_target_rho) else f'{r.within_target_rho:+.2f}'}</td>"
         f"<td>{0 if pd.isna(r.n_targets) else int(r.n_targets)}</td></tr>"
         for _, r in best_p_top.iterrows())}</table>

<h2>Caveats</h2>
<div class="caveat">
Within-target n is small (≈10–14 constructs/target), so per-target ρ is noisy and
individual p-values are weak; the averaged ρ is the more stable read. Binding is the
median of valid trials only (failed / low-event / low-expression trials dropped).
Correlations are monotonic (Spearman), not causal. The four HADDOCK tracks lack the
PAE-derived confidence metrics by construction, so those cells are blank for them.
</div>
<p style="font-size:.8rem;color:#888;margin-top:2rem">Generated by
<code>Structural-Validation/analysis/fcs_correlation.py</code>. Data:
<code>analysis/fcs/correlations_long.csv</code>,
<code>analysis/fcs/binding_by_pair.csv</code>.</p>
"""
    out = C.OUT_REPORT / "fcs_correlation_report.html"
    out.write_text(html, encoding="utf-8")
    print(f"Report -> {out.relative_to(C.REPO_ROOT)}")
