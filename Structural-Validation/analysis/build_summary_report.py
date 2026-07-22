"""
Build a polished, PI-facing summary of the TIMP3 structural-validation study:
one clean multi-panel figure + a self-contained HTML report.

Reads the master tables in analysis/ and writes:
  figures/summary_overview.png
  reports/summary_report.html   (figure embedded as base64 -> single shareable file)

Run:  python Structural-Validation/analysis/build_summary_report.py
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C  # noqa: E402

A = C.OUT_ANALYSIS
CB = {"AF3": "#0072B2", "ESMFold2": "#E69F00", "HADDOCK": "#009E73",
      "accent": "#D55E00", "grey": "#9AA0A6"}
CONSTRUCT_ORDER = ["TIMP3_WT"] + [f"AB_{i}" for i in range(1, 8)] + \
                  [f"C_{i}" for i in range(11, 17)]
SRC_LABEL = {
    "AF3_cofold": "AF3 co-fold", "ESMFold2_cofold": "ESMFold2 co-fold",
    "HADDOCK:AF3xAF3": "HADDOCK AF3×AF3",
    "HADDOCK:ESMFold2xESMFold2": "HADDOCK ESM×ESM",
    "HADDOCK:AF3xCrystal": "HADDOCK AF3×Xtal",
    "HADDOCK:ESMFold2xCrystal": "HADDOCK ESM×Xtal",
}
SRC_ORDER = ["AF3_cofold", "ESMFold2_cofold", "HADDOCK:AF3xAF3",
             "HADDOCK:AF3xCrystal", "HADDOCK:ESMFold2xCrystal",
             "HADDOCK:ESMFold2xESMFold2"]


def _num(df, col):
    return pd.to_numeric(df[col], errors="coerce")


def load():
    cx = pd.read_csv(A / "master_complex_metrics.csv")
    mono = pd.read_csv(A / "master_monomer_metrics.csv")
    agree = pd.read_csv(A / "af3_vs_esmfold_monomer.csv")
    return cx, mono, agree


def per_source_table(cx):
    ok = cx[cx.status == "ok"].copy()
    numcols = ["bsa", "n_hbonds", "min_ca_ca_zincloop", "haddock_score", "dockq",
               "fnonnat", "complex_tm", "pdockq", "pdockq2", "lis", "interface_pae",
               "sc_shape_complementarity", "catalytic_occlusion",
               "n_buried_unsat_hbond"]
    for c in numcols:
        ok[c] = _num(ok, c)

    def m(sub, col):
        return round(sub[col].mean(), 3) if sub[col].notna().any() else np.nan

    rows = []
    for s in SRC_ORDER:
        sub = ok[ok.source == s]
        if sub.empty:
            continue
        a17 = sub[sub.target_id == "ADAM17"]          # native-benchmark rows only
        rows.append({
            "source": SRC_LABEL[s], "n": len(sub),
            "dockq_native_ADAM17": m(a17, "dockq") if len(a17) else np.nan,
            "complex_tm_ADAM17": m(a17, "complex_tm") if len(a17) else np.nan,
            "fnonnat_ADAM17": m(a17, "fnonnat") if len(a17) else np.nan,
            "pdockq2": m(sub, "pdockq2"), "lis": m(sub, "lis"),
            "interface_pae": m(sub, "interface_pae"),
            "sc": m(sub, "sc_shape_complementarity"),
            "cat_occl": m(sub, "catalytic_occlusion"),
            "bsa": round(sub.bsa.mean(), 0),
            "unsat": round(sub.n_buried_unsat_hbond.mean(), 1)
            if sub.n_buried_unsat_hbond.notna().any() else np.nan,
            "haddock_score": round(sub.haddock_score.mean(), 1)
            if sub.haddock_score.notna().any() else np.nan,
        })
    return pd.DataFrame(rows)


def make_figure(cx, mono, agree, out_png):
    ok = cx[cx.status == "ok"].copy()
    for c in ["dockq", "bsa", "n_hbonds"]:
        ok[c] = _num(ok, c)
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle("TIMP3 construct × protease — structural validation summary",
                 fontsize=15, fontweight="bold")

    # (A) DockQ vs native ADAM17 by source
    a17 = ok[ok.target_id == "ADAM17"]
    order = [s for s in SRC_ORDER if s in a17.source.unique()]
    vals = [a17[a17.source == s].dockq.mean() for s in order]
    cols = [CB["AF3"] if "AF3_co" in s else CB["ESMFold2"] if "ESMFold2_co" in s
            else CB["HADDOCK"] for s in order]
    b = ax[0, 0].barh([SRC_LABEL[s] for s in order], vals, color=cols)
    ax[0, 0].axvline(0.23, ls="--", lw=.8, c=CB["grey"])
    ax[0, 0].axvline(0.49, ls="--", lw=.8, c=CB["grey"])
    ax[0, 0].text(0.23, -0.6, "acceptable", fontsize=7, c=CB["grey"])
    ax[0, 0].text(0.49, -0.6, "medium", fontsize=7, c=CB["grey"])
    for r, v in zip(b, vals):
        ax[0, 0].text(v + .01, r.get_y() + r.get_height()/2, f"{v:.2f}", va="center", fontsize=8)
    ax[0, 0].set_xlim(0, 0.9); ax[0, 0].invert_yaxis()
    ax[0, 0].set_xlabel("DockQ vs native TIMP3:ADAM17 crystal")
    ax[0, 0].set_title("A  Native-mode recovery: co-folds win", loc="left", fontweight="bold")

    # (B) monomer fold accuracy AF3 vs ESMFold2 (TM-score vs crystal)
    m = mono[mono.status == "ok"].copy(); m["tm_score"] = _num(m, "tm_score")
    piv = m.pivot_table(index="entity_id", columns="folder", values="tm_score")
    piv = piv.reindex([e for e in CONSTRUCT_ORDER if e in piv.index] +
                      [e for e in piv.index if e not in CONSTRUCT_ORDER])
    x = np.arange(len(piv))
    ax[0, 1].bar(x - .2, piv.get("AF3", pd.Series()), .4, label="AF3", color=CB["AF3"])
    ax[0, 1].bar(x + .2, piv.get("ESMFold2", pd.Series()), .4, label="ESMFold2", color=CB["ESMFold2"])
    ax[0, 1].set_xticks(x); ax[0, 1].set_xticklabels(piv.index, rotation=90, fontsize=6)
    ax[0, 1].set_ylim(0.6, 1.0); ax[0, 1].set_ylabel("TM-score vs crystal")
    ax[0, 1].legend(fontsize=8, loc="lower left")
    ax[0, 1].set_title("B  Monomer fold accuracy (both modelers strong)", loc="left", fontweight="bold")

    # (C) interface size: co-fold vs HADDOCK
    grp = ok.groupby("source").agg(bsa=("bsa", "mean"), hb=("n_hbonds", "mean"))
    order2 = [s for s in SRC_ORDER if s in grp.index]
    xx = np.arange(len(order2))
    ax[1, 0].bar(xx, [grp.loc[s, "bsa"] for s in order2],
                 color=[CB["AF3"] if "AF3_co" in s else CB["ESMFold2"] if "ESMFold2_co" in s
                        else CB["HADDOCK"] for s in order2])
    ax[1, 0].set_xticks(xx); ax[1, 0].set_xticklabels([SRC_LABEL[s] for s in order2],
                                                      rotation=30, ha="right", fontsize=7)
    ax[1, 0].set_ylabel("Buried surface area (Å²)")
    ax[1, 0].set_title("C  Co-folds build larger, native-like interfaces", loc="left", fontweight="bold")

    # (D) AF3 vs ESMFold2 co-fold DockQ heatmap (construct × target)
    co = ok[ok.source == "AF3_cofold"]
    mat = co.pivot_table(index="construct_id", columns="target_id", values="dockq")
    mat = mat.reindex([c for c in CONSTRUCT_ORDER if c in mat.index])
    mat = mat.reindex(columns=[t for t in C.TARGET_ORDER if t in mat.columns])
    im = ax[1, 1].imshow(mat.values, cmap="viridis", aspect="auto", vmin=0, vmax=0.9)
    ax[1, 1].set_xticks(range(mat.shape[1])); ax[1, 1].set_xticklabels(mat.columns, rotation=45)
    ax[1, 1].set_yticks(range(mat.shape[0])); ax[1, 1].set_yticklabels(mat.index, fontsize=7)
    ax[1, 1].set_title("D  AF3 co-fold DockQ (ADAM17=native; others approx.)",
                       loc="left", fontweight="bold")
    fig.colorbar(im, ax=ax[1, 1], fraction=.046, label="DockQ")
    ax[1, 1].grid(False)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def make_extended_figure(cx, out_png):
    """The expanded metric panel: confidence vs correctness + per-source batteries."""
    ok = cx[cx.status == "ok"].copy()
    for c in ["dockq", "pdockq2", "pdockq", "sc_shape_complementarity",
              "catalytic_occlusion", "interface_pae", "fnonnat"]:
        ok[c] = _num(ok, c)
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle("Extended interface-metric panel", fontsize=15, fontweight="bold")

    # (A) confidence vs correctness — does pDockQ2 track native DockQ? (AF3 co-fold)
    co = ok[(ok.source == "AF3_cofold") & ok.pdockq2.notna() & ok.dockq.notna()]
    ax[0, 0].scatter(co.pdockq2, co.dockq, c=CB["AF3"], edgecolor="k", lw=.4, s=36)
    if len(co) > 2:
        r = np.corrcoef(co.pdockq2, co.dockq)[0, 1]
        ax[0, 0].text(.05, .92, f"Pearson r = {r:.2f}", transform=ax[0, 0].transAxes,
                      fontsize=9)
    ax[0, 0].set_xlabel("pDockQ2 (AF3 interface confidence)")
    ax[0, 0].set_ylabel("DockQ vs native")
    ax[0, 0].set_title("A  Does confidence predict native recovery?",
                       loc="left", fontweight="bold")

    # (B) catalytic-cleft occlusion by source
    order = [s for s in SRC_ORDER if s in ok.source.unique()]
    def barsrc(axx, col, ylab, title):
        vals = [ok[ok.source == s][col].mean() for s in order]
        cols = [CB["AF3"] if "AF3_co" in s else CB["ESMFold2"] if "ESMFold2_co" in s
                else CB["HADDOCK"] for s in order]
        axx.bar(range(len(order)), vals, color=cols)
        axx.set_xticks(range(len(order)))
        axx.set_xticklabels([SRC_LABEL[s] for s in order], rotation=30, ha="right",
                            fontsize=7)
        axx.set_ylabel(ylab)
        axx.set_title(title, loc="left", fontweight="bold")
    barsrc(ax[0, 1], "catalytic_occlusion",
           "fraction of zinc loop buried", "B  Catalytic-cleft occlusion")
    barsrc(ax[1, 0], "sc_shape_complementarity",
           "Sc (surface-dot approx.)", "C  Shape complementarity")
    barsrc(ax[1, 1], "fnonnat",
           "fraction non-native contacts", "D  Non-native contact fraction (ADAM17)")

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def _b64(png: Path) -> str:
    return base64.b64encode(png.read_bytes()).decode()


def build_haddock_section() -> str:
    """'Can HADDOCK be improved?' — legacy-generation comparison, the zinc trial,
    and the input-conformation effect. Reads analysis/legacy_haddock_eval.csv."""
    p = A / "legacy_haddock_eval.csv"
    if not p.exists():
        return ""
    df = pd.read_csv(p)
    df = df[(df.get("status") == "ok") & (df.target == "ADAM17")].copy()
    if df.empty:
        return ""
    df["dockq"] = _num(df, "dockq")

    def agg(name):
        s = df[df["set"] == name]["dockq"].dropna()
        if s.empty:
            return None
        return (s.mean(), s.min(), s.max(), len(s))

    rows, order = "", [
        ("legacy_HADDOCK2", "HADDOCK2 (web server, first run)", "unbound monomers"),
        ("legacy_HADDOCK3", "HADDOCK3 (the improvement)", "unbound monomers"),
        ("indep_noZn", "HADDOCK3 reproduction", "unbound monomers"),
        ("indep_Zn", "HADDOCK3 + catalytic Zn²⁺", "unbound monomers"),
        ("indep_ens", "HADDOCK3 + ensemble &amp; flexible edge", "unbound, 3×3 conformers"),
        ("cofoldsplit_noZn", "HADDOCK3 seeded with co-fold", "<b>bound</b> conformations"),
    ]
    for key, label, inp in order:
        a = agg(key)
        if a is None:
            continue
        mean, lo, hi, n = a
        rng = f"[{lo:.3f}–{hi:.3f}]" if n > 1 else "—"
        hl = ' style="background:#fdecea"' if key == "cofoldsplit_noZn" else ""
        rows += (f"<tr{hl}><td>{label}</td><td>{inp}</td><td>{mean:.3f}</td>"
                 f"<td>{rng}</td><td>{n}</td></tr>")
    return f"""
<h2>Can HADDOCK be improved? (three things we tried)</h2>
<p class="sub">All numbers are DockQ against the native TIMP3:ADAM17 co-crystal.
For scale, the AF3/ESMFold2 co-folds reach <b>0.76</b> and CAPRI “medium”; DockQ
0.23 is the “acceptable” floor. Every HADDOCK row below is CAPRI <b>incorrect</b>.</p>
<table><tr><th>Run</th><th>Input conformations</th><th>DockQ mean</th>
<th>range over seeds</th><th>n</th></tr>{rows}</table>
<div class="key"><b>What actually moves the needle — and what doesn't.</b>
<b>(1) HADDOCK2 → HADDOCK3</b> was a real improvement (DockQ 0.030 → 0.085, AIR
energy 166 → 95, violations 5 → 3) — and it is already what the current pipeline
runs, so there is nothing left to recover there.
<b>(2) Adding the catalytic Zn²⁺</b> — absent from every previous run — nudged the
mean (0.046 → 0.078) but the seed ranges overlap and at n=3 vs 3 the comparison
cannot reach significance (Mann-Whitney floor p=0.1). Tellingly HADDOCK's own score
got <i>worse</i> with the metal. Zinc is not the missing ingredient.
<b>(3) Ensemble + interface flexibility also fails</b> (0.034 vs 0.046, i.e. no
better and nominally worse). Note AF3's own five seeds are useless as an ensemble —
they differ by only 0.2–0.8 Å — so this used genuine cross-method conformers
(AF3 + ESMFold2, differing up to 2.6 Å at the reactive edge) plus explicit
semi-flexible segments. Its AIR energies were the best of any arm (50–161) while its
DockQ was the worst: restraint satisfaction improved, realism did not.
<b>(4) Input conformation dominates everything.</b> Seeding HADDOCK with the
co-fold's <i>bound</i> conformations gives 0.78; independently-folded <i>unbound</i>
monomers give 0.046 — a ~15× gap that dwarfs every other factor. HADDOCK's failure
here is an induced-fit problem: rigid-body docking cannot rebuild the bound
interface from unbound monomers, and no amount of metal, sampling or segment
flexibility substitutes for it.</div>
<p style="font-size:.85rem;color:#555">Caveat: the co-fold-seeded row is
<i>semi-circular</i> — it needs the answer (the co-fold) as input, so it is a
legitimate modelling route but not a validation. Also, HADDOCK parameterises the ion
at charge +0.96 rather than a formal +2, so the null zinc result shows “zinc as
HADDOCK models it doesn’t rescue this dock”, not that the metal is irrelevant.
Full data: <code>analysis/legacy_haddock_eval.csv</code>.</p>
"""


def build_fcs_section() -> str:
    """A 'does structure predict binding?' block folded in from the FCS correlation
    layer (fcs_correlation.py). Returns '' if that analysis hasn't been run."""
    corr_csv = C.OUT_ANALYSIS / "fcs" / "correlations_long.csv"
    split_png = C.OUT_FIG / "fcs" / "readout_metric_split.png"
    if not corr_csv.exists():
        return ""
    corr = pd.read_csv(corr_csv)
    av = corr[corr.variant == "all_valid"].dropna(subset=["within_target_rho"]).copy()
    if av.empty:
        return ""
    av["abs"] = av["within_target_rho"].abs()
    headline = av.groupby("readout")["abs"].max().sort_values(ascending=False).index[0]
    top = av[av.readout == headline].sort_values("abs", ascending=False).iloc[0]
    best_auc = corr.dropna(subset=["auc_oriented"]).sort_values(
        "auc_oriented", ascending=False)
    auc_txt = (f" As a binder/non-binder classifier the best structural metric reaches "
               f"oriented AUC ≈ {best_auc.iloc[0].auc_oriented:.2f} (prior best ≈ 0.68) — "
               f"still not a usable ranker." if not best_auc.empty else "")
    img = (f'<img src="data:image/png;base64,{_b64(split_png)}" '
           f'alt="readout x metric split">' if split_png.exists() else "")
    return f"""
<h2>Does structure predict binding? (flow-cytometry correlation)</h2>
<div class="key"><b>No usable predictor — and a caution.</b> Joining the metric
battery to the aggregate FCS data (within-target Spearman), the signal is
<b>readout-specific</b>: for intensity/ratio readouts, confidence &amp; native-
similarity metrics <b>anti</b>-correlate with binding (best |ρ| ≈
{top['abs']:.2f}, e.g. <code>{top.metric}</code>), driven by a promiscuous sticky
construct that scores low in-silico; for efficiency/fraction readouts there is
<b>no correlation</b>.{auc_txt} The negative slope is most likely a
range-restriction / collider artifact (only folding, expressing constructs are
observed) — it must <b>not</b> be used to reward low confidence in design.</div>
{img}
<p style="font-size:.85rem;color:#555">Full analysis, per-target breakdowns and the
purchased-only sensitivity check: <code>reports/fcs_correlation_report.html</code>.</p>
"""


def build_html(tbl: pd.DataFrame, png: Path, ext_png: Path, contact_png: Path,
               cx, out_html, fcs_html: str = "", haddock_html: str = ""):
    n_dock = (cx.source.str.startswith("HADDOCK") & (cx.status == "ok")).sum()

    def cell(v, fmt="{:.2f}"):
        return "—" if pd.isna(v) else fmt.format(v)

    def img(path: Path, alt: str) -> str:
        return (f'<img src="data:image/png;base64,{_b64(path)}" alt="{alt}">'
                if path.exists() else
                f'<p class="caveat">[{alt} not generated]</p>')

    cols = ["Source", "n", "DockQ*", "complex-TM*", "fnonnat*", "pDockQ2†",
            "LIS†", "iPAE† Å", "Sc", "Cat.occl", "BSA Å²", "unsat"]
    thead = "".join(f"<th>{c}</th>" for c in cols)
    trows = ""
    for _, r in tbl.iterrows():
        hl = ' style="background:#eaf6ef;font-weight:600"' if "co-fold" in r.source else ""
        trows += (
            f"<tr{hl}><td>{r.source}</td><td>{r.n}</td>"
            f"<td>{cell(r.dockq_native_ADAM17)}</td>"
            f"<td>{cell(r.complex_tm_ADAM17)}</td>"
            f"<td>{cell(r.fnonnat_ADAM17)}</td>"
            f"<td>{cell(r.pdockq2, '{:.3f}')}</td>"
            f"<td>{cell(r.lis, '{:.3f}')}</td>"
            f"<td>{cell(r.interface_pae, '{:.1f}')}</td>"
            f"<td>{cell(r.sc, '{:.3f}')}</td>"
            f"<td>{cell(r.cat_occl, '{:.2f}')}</td>"
            f"<td>{cell(r.bsa, '{:.0f}')}</td>"
            f"<td>{cell(r.unsat, '{:.1f}')}</td></tr>")
    html = f"""<!doctype html><meta charset=utf-8>
<title>TIMP3 Structural Validation — Summary</title>
<style>
 body{{font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;max-width:960px;
   margin:2rem auto;padding:0 1.2rem;color:#1a1a1a;line-height:1.5}}
 h1{{font-size:1.6rem;margin-bottom:.2rem}} h2{{border-bottom:2px solid #0072B2;padding-bottom:.2rem;margin-top:1.8rem}}
 .sub{{color:#555;margin-top:0}} table{{border-collapse:collapse;width:100%;font-size:.92rem;margin:1rem 0}}
 th,td{{border:1px solid #ddd;padding:.4rem .6rem;text-align:center}} th{{background:#0072B2;color:#fff}}
 td:first-child,th:first-child{{text-align:left}}
 .key{{background:#fff8e1;border-left:4px solid #E69F00;padding:.8rem 1rem;margin:1rem 0;border-radius:4px}}
 .caveat{{background:#f4f4f4;border-left:4px solid #9AA0A6;padding:.6rem 1rem;font-size:.9rem}}
 img{{max-width:100%;border:1px solid #e3e3e3;border-radius:6px;margin:1rem 0}}
 code{{background:#f0f0f0;padding:.05rem .3rem;border-radius:3px}}
</style>
<h1>TIMP3 constructs × proteases — in-silico structural validation</h1>
<p class="sub">14 constructs (WT + AB 1–7, C 11–16) × 5 targets (MMP2, MMP9, ADAM17,
MMP3, ADAM10). Two fold modelers (AF3, ESMFold2), four HADDOCK docking tracks,
and native/approximate DockQ benchmarking. {n_dock} HADDOCK complexes scored.</p>

<div class="key"><b>Headline.</b> Co-folding — with either AlphaFold3 or the free,
local ESMFold2 — recovers the experimentally-observed TIMP3:ADAM17 binding mode
(<b>DockQ 0.76</b>, CAPRI “medium”). Restraint-based HADDOCK docking does not
(DockQ 0.05–0.08), <b>even when docking against the real crystal target</b> — so
the gap is the docking protocol, not the input structure. ESMFold2 co-fold matches
AF3 co-fold on interface recovery.</div>

<h2>Cross-method comparison</h2>
<table><tr>{thead}</tr>{trows}</table>
<p style="font-size:.85rem;color:#555">
<b>*</b> DockQ / complex-TM / fnonnat need a native complex, so they are averaged
over the <b>ADAM17</b> rows only (the sole native co-crystal). <b>†</b> pDockQ2 /
LIS / interface-PAE need the AF3 PAE matrix, so they are populated for AF3 co-folds
(ESMFold2 gets pLDDT-based pDockQ; HADDOCK poses carry neither). <b>Sc</b> =
Lawrence–Colman shape complementarity (surface-dot approximation); <b>Cat.occl</b> =
fraction of the catalytic zinc loop buried by the construct; <b>unsat</b> = buried
unsatisfied interface H-bond donors/acceptors. The full ~40-metric battery per pair
is in <code>master_complex_metrics.csv</code>.</p>

<h2>Overview figure</h2>
{img(png, "summary figure")}

<h2>Extended interface-metric panel</h2>
<p class="sub">Confidence-vs-correctness and the per-source metric batteries that the
richer panel adds on top of DockQ.</p>
{img(ext_png, "extended metrics figure")}

<h2>Reactive-edge × zinc-motif contact matrices</h2>
<p class="sub">Minimum heavy-atom distance (Å) between the TIMP3 reactive N-terminal
edge and each target's catalytic zinc motif, for the best construct per target
(AF3 co-fold). Dark = close contact. This is the notebook's residue-interaction
matrix, generalised across sources; raw matrices for every pair are in
<code>analysis/contact_matrices/</code>.</p>
{img(contact_png, "contact matrix gallery")}
{haddock_html}
{fcs_html}
<h2>What this supports</h2>
<ul>
<li><b>Use co-folds for binding-mode / interface questions.</b> AF3 and ESMFold2
co-folds reproduce the native interface (larger BSA, more H-bonds) and correct
orientation; ESMFold2 is nearly identical to AF3 and is free/local.</li>
<li><b>HADDOCK’s value is different.</b> With loose zinc-loop restraints it finds
compact, favorably-scored poses (score ≈ −97) that satisfy zinc proximity but are
not the crystallographic orientation. It provides physics-based energies, not
native-mode recovery.</li>
<li><b>Both modelers fold the monomers well</b> (TM-score ≈ 0.85–0.99 vs crystal),
so differences at the complex level are about interface prediction, not fold.</li>
</ul>

<h2>Caveats</h2>
<div class="caveat">
Native DockQ exists only for ADAM17; the other targets use approximate homologous
references (1UEA / TIMP3:ADAM17) and should be read as ballpark. All results here
are structural geometry — correlation with the lab’s flow-cytometry binding data is
the next step and is intentionally not included yet.
</div>
<p style="font-size:.8rem;color:#888;margin-top:2rem">Generated by
<code>Structural-Validation/analysis/build_summary_report.py</code>. Full per-pair
data: <code>analysis/master_complex_metrics.csv</code>.</p>
"""
    out_html.write_text(html, encoding="utf-8")


def main():
    C.OUT_FIG.mkdir(parents=True, exist_ok=True)
    C.OUT_REPORT.mkdir(parents=True, exist_ok=True)
    cx, mono, agree = load()
    tbl = per_source_table(cx)
    png = C.OUT_FIG / "summary_overview.png"
    ext_png = C.OUT_FIG / "extended_metrics.png"
    make_figure(cx, mono, agree, png)
    make_extended_figure(cx, ext_png)
    contact_png = C.OUT_FIG / "contact_matrices" / "contact_matrix_gallery_AF3_cofold.png"
    fcs_html = build_fcs_section()
    haddock_html = build_haddock_section()
    html = C.OUT_REPORT / "summary_report.html"
    build_html(tbl, png, ext_png, contact_png, cx, html, fcs_html, haddock_html)
    print("Summary table:\n", tbl.to_string(index=False))
    print(f"\nFigure   -> {png.relative_to(C.REPO_ROOT)}")
    print(f"Extended -> {ext_png.relative_to(C.REPO_ROOT)}")
    print(f"Report   -> {html.relative_to(C.REPO_ROOT)}")


if __name__ == "__main__":
    main()
