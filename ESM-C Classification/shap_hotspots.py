"""SHAP attribution over the grafted loop: which positions/amino acids drive
the model's binding calls ("hot spots").

Treats the loop as L categorical features (one per position, 20 amino-acid
levels each) and runs KernelSHAP against the trained multi-task model, using
the same fixed-scaffold graft mechanics as enumerate_cloop.py (loop_window.py)
so results are exactly consistent with the exhaustive-enumeration outputs.

With only L=6 positions KernelExplainer enumerates coalitions essentially
exactly (2^6 = 64) -- this is cheap per explained sequence, not a 64M sweep.
Background + explained sequences are drawn from the lab's own measured data
for the given loop subtype, so the attributions explain real observed binding
signal, not the unconstrained sequence space.

Each explained loop still needs ~(2^L) x background_size forward passes
through the model, so at 300 explained loops this is genuinely GPU-hours, not
seconds. To avoid a multi-hour silent black box: the explain loop is chunked
(--shap-chunk-size, default 10), progress prints after every chunk, and
shap_values.csv is written incrementally (flushed + fsynced per chunk) so a
kill mid-run only loses the current chunk. Re-running the same command
auto-resumes from shap_values.csv (same explain-set, since it's derived
deterministically from --seed) instead of starting over; pass --no-resume to
force a clean restart. Chunking does not change the SHAP values themselves --
each explained loop is computed independently either way.

Usage
-----
    python shap_hotspots.py --config config.yaml --loop-subtype C-loop
    python shap_hotspots.py --config config.yaml --loop-subtype AB-loop \
        --template-csv ../Data/TIMP3_Binding_Results/TIMP_binder_all.csv
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from esmc_utils import load_config, resolve_path
from inference import load_trained_model
from loop_window import AAS, build_scorer, derive_loop_template

AA_IDX = {a: i for i, a in enumerate(AAS)}


def loops_to_matrix(loops):
    return np.array([[AA_IDX[c] for c in s] for s in loops], dtype=int)


def matrix_to_loops(X):
    return ["".join(AAS[int(v)] for v in row) for row in X]


def load_measured_loops(cfg, loop_subtype, template_csv):
    d = cfg["data"]
    path = template_csv or resolve_path(cfg, d["csv_path"])
    df = pd.read_csv(path)
    sub = "Ligand_Subtype"
    if loop_subtype and sub in df.columns:
        df = df[df[sub] == loop_subtype]
    loops = df[d["loop_col"]].dropna().astype(str)
    L = None
    if len(loops):
        L = len(loops.iloc[0])
        loops = loops[loops.str.len() == L]
    return loops.tolist()


def shap_cols(targets, L):
    return [f"shap_{t}_pos{p+1}" for t in targets for p in range(L)]


def load_resume_state(csv_path, run_meta_path, expected_meta, explain_loops, targets, L):
    """Return (n_done, shap_vals_partial) if a valid, matching partial run exists."""
    if not (csv_path.exists() and run_meta_path.exists()):
        return 0, None
    try:
        prev_meta = json.loads(run_meta_path.read_text())
    except (json.JSONDecodeError, OSError):
        return 0, None
    if any(prev_meta.get(k) != v for k, v in expected_meta.items()):
        print("[warn] existing run_meta.json doesn't match this run's parameters; starting fresh")
        return 0, None
    try:
        prev = pd.read_csv(csv_path)
    except (pd.errors.EmptyDataError, OSError):
        return 0, None
    n_explain = len(explain_loops)
    n_prev = min(len(prev), n_explain)
    if n_prev == 0:
        return 0, None
    if list(prev["loop"].iloc[:n_prev]) != list(explain_loops[:n_prev]):
        print("[warn] existing shap_values.csv doesn't match this run's explain-loop order "
              "(measured data may have changed); starting fresh")
        return 0, None
    shap_vals = np.full((n_explain, L, len(targets)), np.nan)
    for ti, t in enumerate(targets):
        for p in range(L):
            shap_vals[:n_prev, p, ti] = prev[f"shap_{t}_pos{p+1}"].iloc[:n_prev].to_numpy()
    return n_prev, shap_vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--model-dir", default=None)
    ap.add_argument("--loop-subtype", default="C-loop",
                    help="Ligand_Subtype value defining the window to explain "
                         "(e.g. 'C-loop' or 'AB-loop'); pass '' for unfiltered")
    ap.add_argument("--template-csv", default=None,
                    help="CSV to derive the loop window + measured loops from "
                         "(default: cfg data.csv_path)")
    ap.add_argument("--n-explain", type=int, default=300,
                    help="number of measured loops to explain (random sample)")
    ap.add_argument("--background-size", type=int, default=50,
                    help="KernelSHAP background sample size")
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default=None, help="output subfolder name (default: loop_subtype tag)")
    ap.add_argument("--shap-chunk-size", type=int, default=10,
                    help="explain this many loops per KernelSHAP call, checkpointing "
                         "shap_values.csv after each (progress visibility + crash safety)")
    ap.add_argument("--no-resume", action="store_true",
                    help="ignore any existing shap_values.csv and start this run from scratch")
    args = ap.parse_args()

    import shap  # deferred: heavy import, only needed here

    cfg = load_config(args.config)
    model_dir = args.model_dir or (resolve_path(cfg, cfg["output_dir"]) / "model")
    model, tokenizer, meta, device = load_trained_model(model_dir)
    targets = meta["targets"]
    if device.type == "cuda":
        import torch
        cc = torch.cuda.get_device_capability(0)[0]
        model._amp_dtype = torch.bfloat16 if cc >= 8 else torch.float16

    subtype = args.loop_subtype or None
    templ_masked, start_char, L = derive_loop_template(cfg, loop_subtype=subtype,
                                                        csv_path=args.template_csv)
    score = build_scorer(model, tokenizer, meta, templ_masked, start_char, L,
                         batch_size=args.batch_size)
    print(f"window: char {start_char}..{start_char+L-1} | {len(targets)} targets: {targets}")

    measured = load_measured_loops(cfg, subtype, args.template_csv)
    measured = sorted(set(l for l in measured if len(l) == L))
    if len(measured) < args.background_size + 5:
        raise SystemExit(f"Only {len(measured)} unique measured loops for subtype={subtype!r} "
                          f"-- too few for background={args.background_size}")
    print(f"{len(measured)} unique measured loops (subtype={subtype})")

    rng = np.random.default_rng(args.seed)
    bg_loops = rng.choice(measured, size=args.background_size, replace=False)
    remaining = [l for l in measured if l not in set(bg_loops)]
    n_explain = min(args.n_explain, len(remaining))
    explain_loops = rng.choice(remaining, size=n_explain, replace=False).tolist()

    background = loops_to_matrix(bg_loops)

    def f(X: np.ndarray) -> np.ndarray:
        return score(matrix_to_loops(X))

    print(f"Running KernelSHAP: background={len(bg_loops)}, explain={n_explain}, "
          f"positions={L} (=> {2**L} coalitions/instance), chunk_size={args.shap_chunk_size}")
    explainer = shap.KernelExplainer(f, background)

    out_dir = (resolve_path(cfg, cfg["output_dir"]) / "shap" /
              (args.out or (subtype or "all").replace("-", "").lower()))
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "shap_values.csv"
    run_meta_path = out_dir / "run_meta.json"
    run_meta = {"loop_subtype": subtype, "start_char": start_char, "loop_len": L,
                "targets": targets, "n_explain": n_explain, "background_size": args.background_size,
                "seed": args.seed}

    n_done, shap_vals = (0, None)
    if not args.no_resume:
        n_done, shap_vals = load_resume_state(csv_path, run_meta_path, run_meta,
                                              explain_loops, targets, L)
    if shap_vals is None:
        shap_vals = np.full((n_explain, L, len(targets)), np.nan)
    if n_done:
        print(f"[resume] {n_done}/{n_explain} loops already explained, continuing")
    else:
        if csv_path.exists():
            csv_path.unlink()
        run_meta_path.write_text(json.dumps(run_meta, indent=2))

    cols = ["loop"] + shap_cols(targets, L)
    t0 = time.time()
    idx = n_done
    while idx < n_explain:
        j = min(idx + args.shap_chunk_size, n_explain)
        chunk_loops = explain_loops[idx:j]
        chunk_X = loops_to_matrix(chunk_loops)
        sv_chunk = explainer.shap_values(chunk_X, nsamples="auto", silent=True)
        if isinstance(sv_chunk, list):
            sv_chunk = np.stack(sv_chunk, axis=-1)  # (n_chunk, L, n_targets)
        elif sv_chunk.ndim == 2:
            sv_chunk = sv_chunk[:, :, None]
        shap_vals[idx:j] = sv_chunk

        rows = []
        for k, loop in enumerate(chunk_loops):
            row = {"loop": loop}
            for ti, t in enumerate(targets):
                for p in range(L):
                    row[f"shap_{t}_pos{p+1}"] = float(sv_chunk[k, p, ti])
            rows.append(row)
        write_header = idx == 0 and not csv_path.exists()
        with open(csv_path, "a", newline="", encoding="utf-8") as fh:
            pd.DataFrame(rows, columns=cols).to_csv(fh, index=False, header=write_header)
            fh.flush()
            os.fsync(fh.fileno())

        idx = j
        el = time.time() - t0
        done_now = idx - n_done
        rate = done_now / max(el, 1e-9)
        eta_min = (n_explain - idx) / max(rate, 1e-9) / 60
        print(f"  {idx}/{n_explain} explained  {rate*60:.2f} loops/min  "
              f"ETA {eta_min:.1f} min", flush=True)

    if n_done >= n_explain:
        print(f"[resume] all {n_explain} loops already explained; regenerating plots/summary only")

    # --- raw values (already fully written incrementally above; re-derive probs) ---
    probs = score(list(explain_loops))
    summary = {"loop_subtype": subtype, "targets": targets, "n_explain": n_explain}

    # --- per target: position importance + per-position-AA hotspot heatmap ---
    for ti, t in enumerate(targets):
        sv = shap_vals[:, :, ti]  # (n_explain, L)

        pos_importance = np.abs(sv).mean(axis=0)  # (L,)
        plt.figure(figsize=(1.2 * L + 3, 4))
        plt.bar([f"P{p+1}" for p in range(L)], pos_importance)
        plt.ylabel("mean |SHAP|"); plt.title(f"{t}: loop-position importance")
        plt.tight_layout(); plt.savefig(out_dir / f"position_importance_{t}.png", dpi=150)
        plt.close()

        # mean signed SHAP per (position, amino acid actually present at that position)
        hot = np.full((L, 20), np.nan)
        counts = np.zeros((L, 20), dtype=int)
        for i, loop in enumerate(explain_loops):
            for p in range(L):
                a = AA_IDX[loop[p]]
                prev = hot[p, a] if counts[p, a] else 0.0
                hot[p, a] = (prev * counts[p, a] + sv[i, p]) / (counts[p, a] + 1)
                counts[p, a] += 1

        plt.figure(figsize=(1.1 * L + 4, 9))
        sns.heatmap(hot.T, yticklabels=AAS, xticklabels=[f"P{p+1}" for p in range(L)],
                    cmap="RdBu_r", center=0.0, linewidths=0.3,
                    cbar_kws={"label": f"mean SHAP contribution to P({t} bind)"})
        plt.title(f"{t}: per-position amino-acid hotspots (n={counts.sum(axis=1).max()} obs/pos max)")
        plt.xlabel("loop position"); plt.ylabel("amino acid")
        plt.tight_layout(); plt.savefig(out_dir / f"aa_hotspot_{t}.png", dpi=150)
        plt.close()

        # --- advanced SHAP plots ---
        base_val = float(explainer.expected_value[ti]) if isinstance(explainer.expected_value, (list, np.ndarray)) else float(explainer.expected_value)
        features_str = np.array([list(loop) for loop in explain_loops])
        feature_names = [f"P{p+1}" for p in range(L)]
        exp = shap.Explanation(values=sv,
                               base_values=np.full(sv.shape[0], base_val),
                               data=features_str,
                               feature_names=feature_names)

        # 1. Waterfall Plot (first sample)
        plt.figure()
        shap.plots.waterfall(exp[0], show=False)
        plt.savefig(out_dir / f"waterfall_sample_0_{t}.png", bbox_inches='tight', dpi=150)
        plt.close()

        # 2. Individual Force Plot (first sample)
        shap.plots.force(exp[0], show=False, matplotlib=True)
        plt.savefig(out_dir / f"force_sample_0_{t}.png", bbox_inches='tight', dpi=150)
        plt.close()

        # 3. Stacked/Global Force Plot (Interactive HTML)
        try:
            force_html = shap.plots.force(base_val, sv, features_str, feature_names=feature_names, show=False)
            shap.save_html(str(out_dir / f"stacked_force_plot_{t}.html"), force_html)
        except Exception as e:
            print(f"Could not generate stacked force HTML for {t}: {e}")

        # 4. Beeswarm Plot
        try:
            plt.figure()
            # Bypassing color mapping issues with strings by converting them to numeric indices for coloring
            numeric_features = np.vectorize(AA_IDX.get)(features_str)
            exp_numeric = shap.Explanation(values=sv,
                                           base_values=np.full(sv.shape[0], base_val),
                                           data=numeric_features,
                                           feature_names=feature_names)
            shap.plots.beeswarm(exp_numeric, show=False, color_bar=False)
            plt.title(f'{t}: SHAP Beeswarm Plot (Color=AA index)')
            plt.savefig(out_dir / f"beeswarm_{t}.png", bbox_inches='tight', dpi=150)
            plt.close()
        except Exception as e:
            print(f"Could not generate beeswarm for {t}: {e}")

        pd.DataFrame(hot.T, index=list(AAS), columns=[f"P{p+1}" for p in range(L)]
                    ).to_csv(out_dir / f"aa_hotspot_{t}.csv")

        rank = np.argsort(pos_importance)[::-1]
        top_aa_per_pos = {}
        for p in range(L):
            valid = np.where(counts[p] > 0)[0]
            if len(valid):
                best = valid[np.nanargmax(hot[p, valid])]
                worst = valid[np.nanargmin(hot[p, valid])]
                top_aa_per_pos[f"P{p+1}"] = {"most_positive": AAS[best], "most_negative": AAS[worst]}
        summary[t] = {
            "position_importance_rank": [f"P{p+1}" for p in rank],
            "mean_pred_prob": float(probs[:, ti].mean()),
            "top_aa_per_position": top_aa_per_pos,
        }

    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    lines = ["SHAP HOTSPOT SUMMARY", "=" * 40, f"loop_subtype: {subtype}", f"n_explain: {n_explain}", ""]
    for t in targets:
        s = summary[t]
        lines.append(f"{t}: position importance (most->least) = {s['position_importance_rank']}")
        for p, aa in s["top_aa_per_position"].items():
            lines.append(f"    {p}: pushes-toward-bind={aa['most_positive']} "
                         f"pushes-away={aa['most_negative']}")
    (out_dir / "summary.txt").write_text("\n".join(lines))
    print("\n" + "\n".join(lines))
    print(f"\nFigures + tables written to {out_dir}")


if __name__ == "__main__":
    main()
