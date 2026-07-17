"""
Visualize the validation results and assemble an HTML dashboard.

Always produces a pipeline-status overview. As real outputs arrive it adds:
  * model-vs-crystal RMSD / TM-score bars (per entity, per fold modeler)
  * AF3 vs ESMFold2 structural-agreement scatter
  * construct x target heatmaps for each key complex metric (per method)
  * complex-metric Spearman correlation heatmap
  * DockQ / CAPRI summary for targets with a native complex reference (ADAM17)

Figures -> figures/ ; dashboard -> reports/dashboard.html

Run:  python Structural-Validation/analysis/visualize.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C  # noqa: E402

A, F = C.OUT_ANALYSIS, C.OUT_FIG
plt.rcParams.update({"figure.dpi": 120, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.3, "savefig.bbox": "tight"})
# colorblind-safe
CB = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9"]
CONSTRUCT_ORDER = ["TIMP3_WT"] + [f"AB_{i}" for i in range(1, 8)] + \
                  [f"C_{i}" for i in range(11, 17)]


def _read(name):
    p = A / name
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


def _num(df, col):
    return pd.to_numeric(df[col], errors="coerce") if col in df else pd.Series(dtype=float)


def fig_status(figs: list[str]) -> None:
    mono, cplx = _read("master_monomer_metrics.csv"), _read("master_complex_metrics.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, df, title in [(axes[0], mono, "Monomer models"),
                          (axes[1], cplx, "Complexes")]:
        if df.empty or "status" not in df:
            ax.text(0.5, 0.5, "no data", ha="center"); ax.set_title(title); continue
        vc = df["status"].value_counts()
        colors = ["#009E73" if s == "ok" else "#BBBBBB" if "pending" in s
                  else "#D55E00" for s in vc.index]
        ax.barh(vc.index.astype(str), vc.values, color=colors)
        ax.set_title(f"{title} (n={len(df)})")
        ax.set_xlabel("rows")
    fig.suptitle("Pipeline completion status", fontweight="bold")
    p = F / "00_pipeline_status.png"
    fig.savefig(p); plt.close(fig); figs.append(p.name)


def fig_model_vs_crystal(figs: list[str]) -> None:
    mono = _read("master_monomer_metrics.csv")
    ok = mono[mono["status"] == "ok"] if "status" in mono else pd.DataFrame()
    if ok.empty:
        return
    for metric, better in [("rmsd_ca", "lower"), ("tm_score", "higher"),
                           ("lddt", "higher")]:
        piv = ok.pivot_table(index="entity_id", columns="folder",
                             values=metric, aggfunc="first")
        piv = piv.reindex([e for e in CONSTRUCT_ORDER if e in piv.index] +
                          [e for e in piv.index if e not in CONSTRUCT_ORDER])
        ax = piv.plot(kind="bar", figsize=(12, 4), color=CB[:piv.shape[1]])
        ax.set_title(f"Model vs crystal — {metric} ({better} is better)")
        ax.set_ylabel(metric); ax.set_xlabel("")
        p = F / f"10_model_vs_crystal_{metric}.png"
        ax.get_figure().savefig(p); plt.close(ax.get_figure()); figs.append(p.name)


def fig_af3_vs_esm(figs: list[str]) -> None:
    ag = _read("af3_vs_esmfold_monomer.csv")
    ok = ag[ag["status"] == "ok"] if "status" in ag else pd.DataFrame()
    if ok.empty:
        return
    fig, ax = plt.subplots(figsize=(6, 5))
    x, y = _num(ok, "tm_af3_esm"), _num(ok, "rmsd_af3_esm")
    colors = ["#0072B2" if k == "construct" else "#E69F00" for k in ok["kind"]]
    ax.scatter(x, y, c=colors, s=60, edgecolor="k", linewidth=0.5)
    for _, r in ok.iterrows():
        ax.annotate(r["entity_id"], (r["tm_af3_esm"], r["rmsd_af3_esm"]),
                    fontsize=6, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("TM-score (AF3 vs ESMFold2)")
    ax.set_ylabel("CA RMSD Å (AF3 vs ESMFold2)")
    ax.set_title("AF3 vs ESMFold2 structural agreement")
    p = F / "20_af3_vs_esmfold_agreement.png"
    fig.savefig(p); plt.close(fig); figs.append(p.name)


def fig_complex_heatmaps(figs: list[str]) -> None:
    cplx = _read("master_complex_metrics.csv")
    if cplx.empty:
        return
    if "source" not in cplx:
        return
    metrics = ["bsa", "n_hbonds", "contact_density", "sc_shape_complementarity",
               "catalytic_occlusion", "min_ca_ca_zincloop", "haddock_score",
               "iptm", "interface_pae", "pdockq2", "dockq", "complex_tm"]
    for source in sorted(cplx["source"].dropna().unique()):
        sub = cplx[cplx["source"] == source]
        present = [m for m in metrics if m in sub and _num(sub, m).notna().any()]
        if not present:
            continue
        n = len(present)
        fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 4), squeeze=False)
        for ax, m in zip(axes[0], present):
            mat = sub.assign(v=_num(sub, m)).pivot_table(
                index="construct_id", columns="target_id", values="v", aggfunc="first")
            mat = mat.reindex([c for c in CONSTRUCT_ORDER if c in mat.index])
            mat = mat.reindex(columns=[t for t in C.TARGET_ORDER if t in mat.columns])
            im = ax.imshow(mat.values, aspect="auto", cmap="viridis")
            ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=45)
            ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index, fontsize=6)
            ax.set_title(m); fig.colorbar(im, ax=ax, fraction=0.046)
            ax.grid(False)
        fig.suptitle(f"Complex metrics — {source}", fontweight="bold")
        p = F / f"30_complex_heatmaps_{source.replace(':','_').replace('x','X')}.png"
        fig.savefig(p); plt.close(fig); figs.append(p.name)


def fig_correlations(figs: list[str]) -> None:
    for corr in sorted(A.glob("metric_correlations_*.csv")):
        df = pd.read_csv(corr, index_col=0)
        if df.empty:
            continue
        fig, ax = plt.subplots(figsize=(1 + 0.6 * len(df), 1 + 0.6 * len(df)))
        im = ax.imshow(df.values, cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(len(df))); ax.set_xticklabels(df.columns, rotation=90)
        ax.set_yticks(range(len(df))); ax.set_yticklabels(df.index)
        for i in range(len(df)):
            for j in range(len(df)):
                ax.text(j, i, f"{df.values[i, j]:.2f}", ha="center", va="center",
                        fontsize=6)
        ax.set_title(f"Spearman — {corr.stem.replace('metric_correlations_', '')}")
        fig.colorbar(im, ax=ax, fraction=0.046); ax.grid(False)
        p = F / f"40_{corr.stem}.png"
        fig.savefig(p); plt.close(fig); figs.append(p.name)


def fig_dockq(figs: list[str]) -> None:
    cplx = _read("master_complex_metrics.csv")
    if cplx.empty or "dockq" not in cplx:
        return
    dq = cplx.assign(dockq=_num(cplx, "dockq")).dropna(subset=["dockq"])
    if dq.empty:
        return
    fig, ax = plt.subplots(figsize=(11, 4))
    sources = sorted(dq["source"].unique())
    width = 0.8 / max(len(sources), 1)
    for i, source in enumerate(sources):
        s = dq[dq["source"] == source]
        ax.bar(np.arange(len(s)) + i * width, s["dockq"], width=width,
               label=source, color=CB[i % len(CB)])
        ax.set_xticks(np.arange(len(s)))
        ax.set_xticklabels(s["construct_id"] + "→" + s["target_id"], rotation=90, fontsize=6)
    for thr, lab in [(0.23, "acceptable"), (0.49, "medium"), (0.80, "high")]:
        ax.axhline(thr, ls="--", lw=0.7, color="grey")
        ax.text(0, thr, lab, fontsize=6, va="bottom")
    ax.set_ylabel("DockQ"); ax.set_title("DockQ vs native complex reference")
    ax.legend()
    p = F / "50_dockq.png"
    fig.savefig(p); plt.close(fig); figs.append(p.name)


def build_dashboard(figs: list[str]) -> None:
    C.OUT_REPORT.mkdir(parents=True, exist_ok=True)
    rel_fig = Path("..") / "figures"
    csvs = sorted(A.glob("*.csv"))
    html = ["<!doctype html><meta charset=utf-8>",
            "<title>TIMP3 Structural Validation</title>",
            "<style>body{font-family:system-ui;margin:2rem;max-width:1100px}"
            "img{max-width:100%;border:1px solid #ddd;margin:.5rem 0}"
            "h2{border-bottom:2px solid #0072B2;padding-bottom:.2rem}"
            "code{background:#f4f4f4;padding:.1rem .3rem}</style>",
            "<h1>TIMP3 Construct × Target Structural Validation</h1>",
            f"<p>Panel: 14 constructs × 5 targets · folders: {', '.join(C.FOLDERS)} · "
            "docking: HADDOCK3 · generated by <code>analysis/visualize.py</code>.</p>"]
    for fn in figs:
        html.append(f"<h2>{fn}</h2><img src='{rel_fig / fn}'>")
    html.append("<h2>Data tables</h2><ul>")
    for c in csvs:
        html.append(f"<li><code>analysis/{c.name}</code></li>")
    html.append("</ul>")
    (C.OUT_REPORT / "dashboard.html").write_text("\n".join(html))


def main() -> None:
    F.mkdir(parents=True, exist_ok=True)
    for old in F.glob("*.png"):   # avoid stale figures from a prior data state
        old.unlink()
    figs: list[str] = []
    fig_status(figs)
    fig_model_vs_crystal(figs)
    fig_af3_vs_esm(figs)
    fig_complex_heatmaps(figs)
    fig_correlations(figs)
    fig_dockq(figs)
    build_dashboard(figs)
    print(f"Generated {len(figs)} figure(s) -> {F.relative_to(C.REPO_ROOT)}")
    for fn in figs:
        print(f"  {fn}")
    print(f"Dashboard -> {(C.OUT_REPORT / 'dashboard.html').relative_to(C.REPO_ROOT)}")


if __name__ == "__main__":
    main()
