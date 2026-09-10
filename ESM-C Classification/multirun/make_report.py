"""Generate a standalone HTML report for one multirun variant.

Reads whatever this variant's run_all.py steps already produced (data
manifest, test_report, SHAP summaries, UMAP pngs, and -- if enumeration was
run -- the enumeration analysis + any bench_candidates_novel.csv sitting next
to it) and writes a single self-contained report.html into the variant's
output_dir, next to those files. Nothing here re-runs the pipeline; it only
narrates what's already on disk, so it's safe to re-run any time (e.g. after
re-running a step with --force) to refresh the report.

Usage:
    python multirun/make_report.py --config multirun/configs_small/mmp9_other.yaml
    python multirun/make_report.py --config multirun/configs/mmp9_other.yaml --out /tmp/x.html
"""
from __future__ import annotations

import argparse
import html
import json
import sys
from pathlib import Path

MULTIRUN_DIR = Path(__file__).resolve().parent
ESMC_DIR = MULTIRUN_DIR.parent
sys.path.insert(0, str(ESMC_DIR))
from esmc_utils import load_config, resolve_path  # noqa: E402

LOOP_TAGS = ["abloop", "cloop"]
TAG_LABEL = {"abloop": "AB-loop", "cloop": "C-loop"}


def read_json(path: Path):
    if not path.exists():
        return None
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def fmt(x, nd=3):
    if x is None:
        return "&mdash;"
    if isinstance(x, float):
        if x != x:  # NaN
            return "&mdash;"
        return f"{x:.{nd}f}"
    return html.escape(str(x))


def metrics_table(per_target: dict, thresholds: dict) -> str:
    rows = []
    for target, m in per_target.items():
        rows.append(f"""
        <tr>
          <td class="tgt">{html.escape(target)}</td>
          <td>{m.get('n', '&mdash;')}</td>
          <td>{fmt(m.get('pos_rate'), 3)}</td>
          <td class="hi">{fmt(m.get('pr_auc'))}</td>
          <td class="hi">{fmt(m.get('roc_auc'))}</td>
          <td>{fmt(m.get('mcc'))}</td>
          <td>{fmt(m.get('f1'))}</td>
          <td>{fmt(m.get('fbeta'))}</td>
          <td>{fmt((thresholds or {}).get(target))}</td>
        </tr>""")
    return "".join(rows)


def shap_block(tag: str, summary: dict) -> str:
    targets = summary.get("targets", [])
    parts = [f'<h3>{TAG_LABEL.get(tag, tag)} window <span class="dim">(n_explain={summary.get("n_explain", "?")})</span></h3>']
    for t in targets:
        s = summary.get(t, {})
        rank = s.get("position_importance_rank", [])
        aa = s.get("top_aa_per_position", {})
        pos_rows = "".join(
            f"<tr><td>{i+1}</td><td class='tgt'>{p}</td>"
            f"<td class='hi'>{html.escape(aa.get(p, {}).get('most_positive', '?'))}</td>"
            f"<td class='lo'>{html.escape(aa.get(p, {}).get('most_negative', '?'))}</td></tr>"
            for i, p in enumerate(rank)
        )
        parts.append(f"""
        <div class="shap-target">
          <div class="shap-target-head">
            <span class="tgt">{html.escape(t)}</span>
            <span class="dim">mean predicted P(bind) over explained set = {fmt(s.get('mean_pred_prob'))}</span>
          </div>
          <table>
            <thead><tr><th>Rank</th><th>Position</th><th>Favors binding</th><th>Disfavors binding</th></tr></thead>
            <tbody>{pos_rows}</tbody>
          </table>
        </div>""")
    return "".join(parts)


def img_tag(path: Path, rel_base: Path, alt: str) -> str:
    if not path.exists():
        return f'<div class="missing">missing: {html.escape(path.name)}</div>'
    rel = Path(path).resolve().relative_to(rel_base.resolve()) if _is_relative(path, rel_base) else path
    return f'<img src="{html.escape(str(rel).replace(chr(92), "/"))}" alt="{html.escape(alt)}">'


def _is_relative(path: Path, base: Path) -> bool:
    try:
        Path(path).resolve().relative_to(base.resolve())
        return True
    except ValueError:
        return False


def build_nav(out_dir: Path, variant_name: str) -> str:
    """Link to sibling variants' report.html, if present (built or not yet built)."""
    siblings_root = out_dir.parent  # .../esmc_multirun[_small]/
    if not siblings_root.is_dir():
        return ""
    items = []
    for d in sorted(siblings_root.iterdir()):
        if not d.is_dir():
            continue
        label = d.name
        target = d / "report.html"
        cls = "current" if d.name == variant_name else ""
        if target.exists() or d.name == variant_name:
            rel = f"../{d.name}/report.html"
            items.append(f'<a class="navlink {cls}" href="{html.escape(rel)}">{html.escape(label)}</a>')
    if not items:
        return ""
    return f'<nav class="topnav">{"".join(items)}</nav>'


def build_summary(variant_name, targets, n_seq, test_report, smoke_note, shap_summaries, enum_sections) -> str:
    """A short prose paragraph synthesizing the numbers below, so the report reads
    top-to-bottom rather than requiring the reader to assemble the takeaway themselves."""
    sentences = []
    n_t = len(targets)
    sentences.append(
        f"This variant fine-tunes on {n_seq:,} labeled loop sequences across "
        f"{n_t} target{'s' if n_t != 1 else ''} ({', '.join(targets)})."
        if n_seq else f"This variant covers {n_t} target{'s' if n_t != 1 else ''} ({', '.join(targets)})."
    )
    if smoke_note:
        sentences.append("Only a tiny smoke-test model has been trained so far, so no real "
                          "performance conclusions can be drawn yet.")
    elif test_report:
        per_t = test_report.get("per_target", {})
        ranked = sorted(
            ((t, m.get("roc_auc")) for t, m in per_t.items() if m.get("roc_auc") == m.get("roc_auc")),
            key=lambda x: x[1], reverse=True,
        )
        if ranked:
            best_t, best_auc = ranked[0]
            sentences.append(
                f"On held-out test data, <b>{html.escape(best_t)}</b> is discriminated best "
                f"(ROC-AUC={best_auc:.3f}, PR-AUC={per_t[best_t].get('pr_auc', float('nan')):.3f})."
            )
            if len(ranked) > 1:
                worst_t, worst_auc = ranked[-1]
                if worst_auc < 0.65:
                    sentences.append(
                        f"<b>{html.escape(worst_t)}</b> discriminates poorly (ROC-AUC={worst_auc:.3f}), "
                        f"likely reflecting a small or heavily imbalanced test slice rather than a model deficiency."
                    )
                elif worst_t != best_t:
                    sentences.append(f"<b>{html.escape(worst_t)}</b> trails at ROC-AUC={worst_auc:.3f}.")
    if shap_summaries:
        first_tag = next(iter(shap_summaries))
        s = shap_summaries[first_tag]
        first_target = next(iter(s.get("targets", [])), None)
        if first_target:
            rank = s.get(first_target, {}).get("position_importance_rank", [])
            if rank:
                sentences.append(
                    f"SHAP attribution over the {TAG_LABEL.get(first_tag, first_tag)} window ranks "
                    f"position {rank[0]} as most influential for {html.escape(first_target)}."
                )
    if enum_sections:
        names = ", ".join(html.escape(n) for n, _, _ in enum_sections)
        sentences.append(f"A full enumeration sweep ({names}) has also been completed for this variant.")
    return " ".join(sentences)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, help="path to a multirun variant yaml (relative to ESM-C Classification/, or absolute)")
    ap.add_argument("--out", default=None, help="output html path (default: <output_dir>/report.html)")
    args = ap.parse_args()

    cfg_path = Path(args.config)
    if not cfg_path.is_absolute():
        cfg_path = ESMC_DIR / cfg_path
    cfg = load_config(cfg_path)
    variant_name = cfg_path.stem
    out_dir = resolve_path(cfg, cfg["output_dir"])
    out_path = Path(args.out) if args.out else out_dir / "report.html"

    manifest = read_json(out_dir / "data" / "manifest.json") or {}
    test_report = read_json(out_dir / "model" / "test_report.json")
    if test_report is None:
        test_report = read_json(out_dir / "model_smoke" / "test_report.json")
        smoke_note = True
    else:
        smoke_note = False

    shap_summaries = {}
    for tag in LOOP_TAGS:
        s = read_json(out_dir / "shap" / tag / "summary.json")
        if s:
            shap_summaries[tag] = s

    enum_root = out_dir / "enumeration"
    enum_dirs = sorted(enum_root.glob("*")) if enum_root.is_dir() else []
    enum_sections = []
    for ed in enum_dirs:
        if not ed.is_dir():
            continue
        analysis = read_json(ed / "analysis" / "summary.json")
        bench_csv = ed / "analysis" / "bench_candidates_novel.csv"
        if analysis is None and not bench_csv.exists():
            continue
        enum_sections.append((ed.name, analysis, bench_csv if bench_csv.exists() else None))

    targets = cfg.get("targets", [])
    model_id = cfg.get("model", {}).get("model_id", "?")
    data_cfg = cfg.get("data", {})

    n_seq = manifest.get("n_sequences")
    splits = manifest.get("splits", {})
    train_n = splits.get("train", {}).get("n_sequences")
    val_n = splits.get("val", {}).get("n_sequences")
    test_n = splits.get("test", {}).get("n_sequences")

    dataset_rows = []
    for split_name in ("train", "val", "test"):
        sp = splits.get(split_name, {})
        per_t = sp.get("per_target", {})
        for t in targets:
            m = per_t.get(t)
            if not m:
                continue
            dataset_rows.append(
                f"<tr><td>{split_name}</td><td class='tgt'>{html.escape(t)}</td>"
                f"<td>{m.get('n_assayed', '&mdash;')}</td><td>{m.get('n_pos', '&mdash;')}</td>"
                f"<td>{fmt(m.get('pos_rate'), 3)}</td></tr>"
            )
    dataset_table = "".join(dataset_rows)

    hamming_exposure = manifest.get("hamming1_exposure", {})

    summary_text = build_summary(variant_name, targets, n_seq, test_report, smoke_note, shap_summaries, enum_sections)

    smoke_banner = ""
    if smoke_note:
        smoke_banner = ('<div class="callout warn"><b>Note:</b> only a <code>--smoke</code> '
                         '(tiny sanity-pass) model exists for this variant &mdash; the numbers '
                         'below are not a real trained model\'s performance.</div>')

    metrics_html = ""
    if test_report:
        metrics_html = f"""
        <table class="metrics">
          <thead><tr><th>Target</th><th>N (test)</th><th>Pos. rate</th><th>PR-AUC</th><th>ROC-AUC</th>
          <th>MCC</th><th>F1</th><th>F&beta;(0.5)</th><th>Threshold</th></tr></thead>
          <tbody>{metrics_table(test_report.get('per_target', {}), test_report.get('thresholds', {}))}</tbody>
        </table>"""
    else:
        metrics_html = '<div class="missing">no test_report.json found</div>'

    shap_html = "".join(shap_block(tag, s) for tag, s in shap_summaries.items()) or '<div class="missing">no SHAP summaries found</div>'

    viz_dir = out_dir / "visualizations"
    viz_html = f"""
    <div class="viz-grid">
      <figure>{img_tag(viz_dir / 'target.png', out_dir, 'UMAP colored by assayed target')}<figcaption>colored by assayed target</figcaption></figure>
      <figure>{img_tag(viz_dir / 'binding.png', out_dir, 'UMAP colored by binding label')}<figcaption>colored by true binding label</figcaption></figure>
    </div>"""

    enum_html = ""
    for name, analysis, bench_csv in enum_sections:
        consensus = ""
        novelty = ""
        if analysis:
            consensus = "".join(
                f"<li><b>{html.escape(t)}</b>: <code>{html.escape(m)}</code></li>"
                for t, m in analysis.get("top_consensus", {}).items()
            )
            novelty = "".join(
                f"<li><b>{html.escape(t)}</b>: {fmt(n.get('frac_in_train'), 4)} of top-K already in "
                f"training data; {fmt(n.get('frac_within1_edit_first5k'), 4)} of top-5K within 1 edit of a "
                f"measured loop</li>"
                for t, n in analysis.get("novelty_topK", {}).items()
            )
        bench_note = ""
        if bench_csv:
            rel = bench_csv.relative_to(out_dir) if _is_relative(bench_csv, out_dir) else bench_csv
            bench_note = (f'<div class="callout"><b>Bench shortlist:</b> curated novel candidates at '
                           f'<code>{html.escape(str(rel).replace(chr(92), "/"))}</code></div>')
        enum_html += f"""
        <div class="enum-block">
          <h3>Enumeration sweep: {html.escape(name)}</h3>
          <ul class="dim">{consensus}</ul>
          <ul class="dim">{novelty}</ul>
          {bench_note}
        </div>"""

    nav = build_nav(out_dir, variant_name)

    sections = [
        ("Dataset", f"""
        <div class="dim">source: <code>{html.escape(str(data_cfg.get('csv_path', '?')))}</code></div>
        <div class="dim">{n_seq} total sequences &middot; train/val/test = {train_n}/{val_n}/{test_n} (loop-group split, leakage-safe)</div>
        <table>
          <thead><tr><th>Split</th><th>Target</th><th>N assayed</th><th>N positive</th><th>Pos. rate</th></tr></thead>
          <tbody>{dataset_table}</tbody>
        </table>
        {f'<div class="dim" style="margin-top:8px;">Hamming-1 exposure (fraction of val/test within 1 edit of a train loop): val={fmt(hamming_exposure.get("val"),3)}, test={fmt(hamming_exposure.get("test"),3)}</div>' if hamming_exposure else ''}
        """),
        ("Held-out test performance", f"{smoke_banner}{metrics_html}"),
        ("SHAP hotspots", shap_html),
        ("UMAP embedding (test split)", viz_html),
    ]
    if enum_html:
        sections.append(("Enumeration &amp; bench candidates", enum_html))

    section_html = "".join(
        f'<section><h2><span class="secnum">{i+1}</span>{title}</h2>{body}</section>'
        for i, (title, body) in enumerate(sections)
    )

    html_doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{html.escape(variant_name)} &mdash; ESM-C multirun report</title>
<style>
  :root {{
    --bg: #0f1117; --panel: #171a23; --border: #262b38; --text: #e5e9f0;
    --dim: #8a93a6; --blue: #6fa8ff; --amber: #f0b93d; --green: #4fd19a; --red: #e07a7a;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    background: var(--bg); color: var(--text); margin: 0; padding: 0;
    font-family: Georgia, "Times New Roman", serif;
    line-height: 1.6;
  }}
  .wrap {{ max-width: 860px; margin: 0 auto; padding: 40px 28px 90px; }}
  .topnav {{ display: flex; gap: 8px; flex-wrap: wrap; margin-bottom: 30px; font-family: -apple-system, "Segoe UI", sans-serif; }}
  .navlink {{
    font-size: 12px; padding: 5px 10px; border-radius: 6px; border: 1px solid var(--border);
    color: var(--dim); text-decoration: none;
  }}
  .navlink.current {{ color: var(--blue); border-color: var(--blue); }}
  .navlink:hover {{ color: var(--text); }}
  header.title {{ margin-bottom: 8px; border-bottom: 2px solid var(--border); padding-bottom: 20px; }}
  header.title h1 {{ font-size: 28px; margin: 0 0 8px; }}
  header.title .sub {{ color: var(--dim); font-size: 13px; font-family: -apple-system, "Segoe UI", sans-serif; }}
  .summary {{ font-size: 15px; margin: 22px 0 6px; padding: 16px 18px; background: var(--panel);
              border-left: 3px solid var(--blue); border-radius: 4px; }}
  .summary-label {{ font-family: -apple-system, "Segoe UI", sans-serif; font-size: 10.5px; color: var(--blue);
                     text-transform: uppercase; letter-spacing: .06em; margin-bottom: 8px; display: block; }}
  section {{ margin: 34px 0; }}
  h2 {{ font-size: 16px; color: var(--text); margin: 0 0 14px; font-weight: 700;
        font-family: -apple-system, "Segoe UI", sans-serif; border-bottom: 1px solid var(--border); padding-bottom: 8px; }}
  .secnum {{ color: var(--blue); margin-right: 8px; }}
  h3 {{ font-size: 13.5px; margin: 16px 0 8px; font-family: -apple-system, "Segoe UI", sans-serif; color: var(--dim); font-weight: 600; }}
  .dim {{ color: var(--dim); font-size: 12.5px; font-family: -apple-system, "Segoe UI", sans-serif; }}
  code {{ background: #1b2130; padding: 1px 5px; border-radius: 4px; font-size: 12px; font-family: Consolas, monospace; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 13px; margin-top: 6px;
           font-family: -apple-system, "Segoe UI", sans-serif; }}
  th, td {{ padding: 6px 10px; border-bottom: 1px solid var(--border); text-align: right; }}
  th:first-child, td:first-child {{ text-align: left; }}
  td.tgt {{ text-align: left; font-weight: 600; }}
  th {{ color: var(--dim); font-weight: 500; text-transform: uppercase; font-size: 10.5px; letter-spacing: .03em; }}
  td.hi {{ color: var(--green); font-weight: 600; }}
  td.lo {{ color: var(--red); }}
  .callout {{
    background: rgba(240, 185, 61, 0.08); border: 1px solid rgba(240, 185, 61, 0.25);
    border-radius: 6px; padding: 10px 14px; font-size: 12.5px; margin-top: 12px;
    font-family: -apple-system, "Segoe UI", sans-serif;
  }}
  .callout.warn {{ background: rgba(255, 90, 90, 0.08); border-color: rgba(255, 90, 90, 0.3); }}
  .missing {{ color: var(--dim); font-size: 12.5px; font-style: italic; font-family: -apple-system, "Segoe UI", sans-serif; }}
  .shap-target {{ margin-bottom: 20px; }}
  .shap-target-head {{ display: flex; justify-content: space-between; align-items: baseline; margin-bottom: 6px;
                        font-family: -apple-system, "Segoe UI", sans-serif; }}
  .viz-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 18px; }}
  .viz-grid figure {{ margin: 0; }}
  .viz-grid img {{ width: 100%; border-radius: 6px; border: 1px solid var(--border); }}
  .viz-grid figcaption {{ font-size: 11px; color: var(--dim); text-align: center; margin-top: 4px;
                           font-family: -apple-system, "Segoe UI", sans-serif; }}
  .enum-block {{ margin-top: 10px; }}
  .enum-block ul {{ margin: 4px 0; padding-left: 18px; font-family: -apple-system, "Segoe UI", sans-serif; }}
  footer {{ margin-top: 50px; color: var(--dim); font-size: 11px; font-family: -apple-system, "Segoe UI", sans-serif; }}
  @media (max-width: 700px) {{ .viz-grid {{ grid-template-columns: 1fr; }} }}
</style>
</head>
<body>
<div class="wrap">
  {nav}
  <header class="title">
    <h1>{html.escape(variant_name)}</h1>
    <div class="sub">ESM-C multirun variant report &middot; model: <code>{html.escape(model_id)}</code> &middot; targets: {", ".join(html.escape(t) for t in targets)}</div>
  </header>

  <div class="summary"><span class="summary-label">Summary</span>{summary_text}</div>

  {section_html}

  <footer>Generated by multirun/make_report.py from files under this variant's output_dir. Re-run it any time to refresh after re-running a step.</footer>
</div>
</body>
</html>
"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html_doc, encoding="utf-8")
    print(f"[ok] wrote {out_path}")


if __name__ == "__main__":
    main()
