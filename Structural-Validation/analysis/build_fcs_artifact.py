"""
Assemble a polished, theme-aware, self-contained HTML artifact of the structure->
binding correlation finding (for sharing on claude.ai). Reuses the figures and the
numbers already produced by fcs_correlation.py; recomputes only the few headline
values so it stays decoupled.

Output: reports/fcs_correlation_artifact.html  (body content only — no <doctype>/
<head>/<body>; the Artifact host wraps it).

Run:  python Structural-Validation/analysis/build_fcs_artifact.py
"""
from __future__ import annotations

import base64
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C  # noqa: E402

FCS = C.OUT_ANALYSIS / "fcs"
FIG = C.OUT_FIG / "fcs"
CORE = ["iptm", "pdockq", "fnat", "dockq", "complex_tm", "iface_composite",
        "catalytic_occlusion"]
SRC = {"AF3_cofold": "AF3 co-fold", "ESMFold2_cofold": "ESMFold2 co-fold"}


def _img(name: str) -> str:
    p = FIG / name
    return base64.b64encode(p.read_bytes()).decode() if p.exists() else ""


def compute():
    corr = pd.read_csv(FCS / "correlations_long.csv")
    binding = pd.read_csv(FCS / "binding_by_pair.csv")
    corr = corr[corr.variant == "all_valid"].copy()   # headline uses the primary variant
    av = corr.dropna(subset=["within_target_rho"]).copy()
    av["abs"] = av["within_target_rho"].abs()
    headline = av.groupby("readout")["abs"].max().sort_values(ascending=False).index[0]
    top = av[av.readout == headline].sort_values("abs", ascending=False).iloc[0]
    auc = corr.dropna(subset=["auc_oriented"]).sort_values("auc_oriented", ascending=False).iloc[0]
    tb = binding.groupby("construct_id")[headline].median().sort_values(ascending=False)
    top_binder = tb.index[0]
    wt_rank = int(tb.rank(ascending=False).get("TIMP3_WT")) if "TIMP3_WT" in tb.index else None
    core = av[av.source.isin(SRC) & av.metric.isin(CORE)]
    best_by_r = core.groupby("readout")["abs"].max()
    responsive = [r for r in best_by_r.index if best_by_r[r] >= 0.40]
    flat = [r for r in best_by_r.index if best_by_r[r] < 0.25]
    n_pairs = len(binding)
    return dict(headline=headline, top=top, auc=auc, top_binder=top_binder,
                wt_rank=wt_rank, responsive=responsive, flat=flat, n_pairs=n_pairs,
                n_constructs=binding.construct_id.nunique())


STYLE = """<style>
  :root{
    --ground:#FBFCFD; --panel:#FFFFFF; --tile:#F1F6F8; --ink:#17212B;
    --muted:#586772; --faint:#7C8994; --rule:#E4EAEE; --accent:#0E7C86;
    --neg:#B5432B; --pos:#1F7A52; --shadow:0 1px 2px rgba(20,40,55,.06),0 6px 20px rgba(20,40,55,.05);
  }
  @media (prefers-color-scheme:dark){:root{
    --ground:#0E1519; --panel:#151E25; --tile:#18242C; --ink:#E7EEF3;
    --muted:#98A8B4; --faint:#6E7F8B; --rule:#243139; --accent:#45BCC6;
    --neg:#E17A5B; --pos:#59B98C; --shadow:0 1px 2px rgba(0,0,0,.3),0 8px 24px rgba(0,0,0,.28);
  }}
  :root[data-theme="dark"]{
    --ground:#0E1519; --panel:#151E25; --tile:#18242C; --ink:#E7EEF3;
    --muted:#98A8B4; --faint:#6E7F8B; --rule:#243139; --accent:#45BCC6;
    --neg:#E17A5B; --pos:#59B98C; --shadow:0 1px 2px rgba(0,0,0,.3),0 8px 24px rgba(0,0,0,.28);
  }
  :root[data-theme="light"]{
    --ground:#FBFCFD; --panel:#FFFFFF; --tile:#F1F6F8; --ink:#17212B;
    --muted:#586772; --faint:#7C8994; --rule:#E4EAEE; --accent:#0E7C86;
    --neg:#B5432B; --pos:#1F7A52;
  }
  .fcs{--sans:system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
    --mono:ui-monospace,"SF Mono","Cascadia Code",Menlo,Consolas,monospace;
    background:var(--ground); color:var(--ink); font-family:var(--sans);
    line-height:1.6; -webkit-font-smoothing:antialiased;
    padding:clamp(1.2rem,4vw,3rem) clamp(1rem,4vw,2rem); }
  .fcs *{box-sizing:border-box;}
  .wrap{max-width:1000px; margin:0 auto;}
  .prose{max-width:68ch;}
  .eyebrow{font-family:var(--mono); font-size:.72rem; letter-spacing:.16em;
    text-transform:uppercase; color:var(--accent); margin:0 0 .7rem;}
  h1{font-size:clamp(1.7rem,4vw,2.5rem); line-height:1.12; letter-spacing:-.02em;
    text-wrap:balance; margin:0 0 .6rem; font-weight:680;}
  h2{font-size:1.28rem; letter-spacing:-.01em; text-wrap:balance; margin:2.6rem 0 .4rem;
    padding-bottom:.35rem; border-bottom:1px solid var(--rule);}
  .stand{font-size:1.08rem; color:var(--muted); max-width:64ch; margin:.2rem 0 0;}
  p{margin:.7rem 0;} a{color:var(--accent);}
  code{font-family:var(--mono); font-size:.86em; background:var(--tile);
    padding:.08em .38em; border-radius:4px;}
  .tiles{display:grid; grid-template-columns:repeat(auto-fit,minmax(210px,1fr));
    gap:1rem; margin:2rem 0;}
  .tile{background:var(--panel); border:1px solid var(--rule); border-radius:12px;
    padding:1.1rem 1.2rem; box-shadow:var(--shadow);}
  .tile .lab{font-family:var(--mono); font-size:.68rem; letter-spacing:.12em;
    text-transform:uppercase; color:var(--faint); margin:0 0 .5rem;}
  .tile .val{font-size:2rem; font-weight:680; letter-spacing:-.02em;
    font-variant-numeric:tabular-nums; line-height:1;}
  .tile .sub{font-size:.82rem; color:var(--muted); margin:.45rem 0 0;}
  .neg{color:var(--neg);} .pos{color:var(--pos);}
  .callout{background:var(--panel); border:1px solid var(--rule);
    border-left:4px solid var(--accent); border-radius:10px; padding:1rem 1.2rem;
    margin:1.3rem 0; box-shadow:var(--shadow);}
  .callout.warn{border-left-color:var(--neg);}
  .callout .lab{font-family:var(--mono); font-size:.7rem; letter-spacing:.12em;
    text-transform:uppercase; color:var(--faint); margin:0 0 .4rem;}
  figure{margin:1.4rem 0;}
  .fig{background:#fff; border:1px solid var(--rule); border-radius:12px;
    padding:.9rem; overflow-x:auto; box-shadow:var(--shadow);}
  .fig img{display:block; width:100%; min-width:640px; height:auto; border-radius:4px;}
  figcaption{font-size:.85rem; color:var(--muted); margin-top:.55rem; max-width:70ch;}
  .foot{margin-top:2.6rem; padding-top:1rem; border-top:1px solid var(--rule);
    font-size:.8rem; color:var(--faint);}
  :focus-visible{outline:2px solid var(--accent); outline-offset:2px;}
</style>"""


def build():
    d = compute()
    top, auc = d["top"], d["auc"]
    rho = top.within_target_rho
    dir_word = "anti-correlates with" if rho < 0 else "tracks"
    rho_cls = "neg" if rho < 0 else "pos"
    resp = ", ".join(d["responsive"]) or "—"
    flat = ", ".join(d["flat"]) or "—"
    wt = f" (wild-type ranks #{d['wt_rank']} of {d['n_constructs']})" if d["wt_rank"] else ""
    auc_dir = "anti-directional" if auc.auc < 0.5 else "positive"

    html = f"""<title>TIMP3 · Does Structure Predict Binding?</title>
{STYLE}
<div class="fcs"><div class="wrap">
  <p class="eyebrow">TIMP3 protease inhibitors · in-silico ↔ wet-lab</p>
  <h1>Does predicted structure predict measured binding?</h1>
  <p class="stand prose">A ~40-metric structural battery (AlphaFold3 &amp; ESMFold2
  co-folds, HADDOCK docking) scored against the lab's aggregate flow-cytometry
  binding data across {d['n_pairs']} construct×target pairs. Correlations are
  <b>within-target</b> Spearman — comparing constructs against each other for a
  fixed target, so target-to-target differences in bindability can't masquerade
  as a predictive signal.</p>

  <div class="tiles">
    <div class="tile">
      <p class="lab">Best structural predictor</p>
      <div class="val {rho_cls}">{rho:+.2f}</div>
      <p class="sub">Spearman ρ · <code>{top.metric}</code>, {SRC.get(top.source, top.source)}
      vs {d['headline']}</p>
    </div>
    <div class="tile">
      <p class="lab">Binder / non-binder AUC</p>
      <div class="val">{auc.auc_oriented:.2f}</div>
      <p class="sub">direction-agnostic ({auc_dir}); prior calibration best ≈ 0.68</p>
    </div>
    <div class="tile">
      <p class="lab">Who actually binds most</p>
      <div class="val">{d['top_binder']}</div>
      <p class="sub">top binder vs every target{wt} — yet scores low in-silico</p>
    </div>
  </div>

  <div class="callout warn">
    <p class="lab">Bottom line</p>
    <p style="margin:.2rem 0"><b>No usable structural predictor of binding — and the
    strongest signal points the wrong way.</b> For intensity/ratio readouts the
    confidence &amp; native-similarity metrics <span class="{rho_cls}">{dir_word}</span>
    binding (best ρ = {rho:+.2f}, <code>{top.metric}</code>), because a promiscuous,
    sticky construct (<b>{d['top_binder']}</b>) dominates the assay while the models
    rate it poorly. For efficiency / fraction-of-cells readouts there is essentially
    <b>no</b> correlation — reproducing the lab's earlier "no predictive power" result.</p>
  </div>

  <h2>The signal is readout-specific</h2>
  <p class="prose">Only some binding readouts respond at all, and the ones that do
  respond <i>negatively</i>. Intensity- and ratio-based readouts ({resp}) carry the
  anti-correlation; efficiency / fraction readouts ({flat}) are flat. This split is
  the essential caveat on any single headline number.</p>
  <figure><div class="fig"><img alt="within-target rho by readout and metric"
    src="data:image/png;base64,{_img('readout_metric_split.png')}"></div>
    <figcaption>Within-target ρ (co-folds) for each metric × readout. Blue = the
    metric anti-predicts binding; white = no relationship. The negative block is
    confined to the intensity/ratio readouts on the left-of-center columns.</figcaption>
  </figure>

  <h2>The best predictor, dissected</h2>
  <p class="prose">The top-ranked metric, per target and pooled. Colouring the
  scatter by target shows the trend lives <i>within</i> each target — pooling all
  points together would blur it.</p>
  <figure><div class="fig"><img alt="best predictor per target and scatter"
    src="data:image/png;base64,{_img('best_predictor.png')}"></div>
    <figcaption>Per-target Spearman ρ (left) and the underlying scatter coloured by
    target (right). The stickiest, lowest-scoring constructs bind hardest.</figcaption>
  </figure>

  <div class="callout">
    <p class="lab">Is the anti-correlation causal? Probably not</p>
    <p style="margin:.2rem 0">The lab's exact-sequence calibration
    (<code>calibrated_scoring.py</code>) saw the same negative slope and attributed
    it to <b>range restriction / collider bias</b>: only constructs that were made,
    and that fold and express, are ever observed — the non-folding non-binders that
    would anchor the low-confidence end are missing. So "low confidence → binds
    more" must <b>not</b> become a design rule; rewarding low ipTM / pDockQ on new
    designs would be optimizing a within-sample artifact. Structural models stay a
    foldability <i>filter</i>, never a directional ranker.</p>
  </div>

  <h2>As a binder / non-binder classifier</h2>
  <p class="prose">Framed as classification (binder = above the per-target median),
  the best structural metric reaches an oriented AUC of
  <b>{auc.auc_oriented:.2f}</b> ({auc.metric}, {SRC.get(auc.source, auc.source)}) —
  marginally above the prior best of ~0.68 but the same weak ballpark, and
  {auc_dir}. The richer battery does not turn structure into a binder ranker.</p>

  <p class="foot">Within-target n ≈ 10–14 constructs/target, so per-target ρ is
  noisy — the averaged value is the stable read. Binding = median of valid trials
  (failed / low-event / low-expression dropped). Monotonic (Spearman) associations,
  not causal. Generated from <code>analysis/fcs/correlations_long.csv</code>; full
  report, per-target detail and purchased-only sensitivity in
  <code>reports/fcs_correlation_report.html</code>.</p>
</div></div>"""
    out = C.OUT_REPORT / "fcs_correlation_artifact.html"
    out.write_text(html, encoding="utf-8")
    print(f"Artifact HTML -> {out}")
    print(f"headline={d['headline']!r} top={top.metric}/{top.source} rho={rho:+.2f} "
          f"auc={auc.auc_oriented:.2f} top_binder={d['top_binder']}")
    return out


if __name__ == "__main__":
    build()
