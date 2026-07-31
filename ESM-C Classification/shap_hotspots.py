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

Usage
-----
    python shap_hotspots.py --config config.yaml --loop-subtype C-loop
    python shap_hotspots.py --config config.yaml --loop-subtype AB-loop \
        --template-csv ../Data/TIMP3_Binding_Results/TIMP_binder_all.csv
"""
from __future__ import annotations

import argparse
import json
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
    explain_loops = rng.choice(remaining, size=n_explain, replace=False)

    background = loops_to_matrix(bg_loops)

    def f(X: np.ndarray) -> np.ndarray:
        return score(matrix_to_loops(X))

    print(f"Running KernelSHAP: background={len(bg_loops)}, explain={n_explain}, "
          f"positions={L} (=> {2**L} coalitions/instance)")
    explainer = shap.KernelExplainer(f, background)
    explain_X = loops_to_matrix(explain_loops)
    shap_vals = explainer.shap_values(explain_X, nsamples="auto", silent=True)
    # shap_vals: list[n_targets] of (n_explain, L), or a single (n_explain, L, n_targets) array
    if isinstance(shap_vals, list):
        shap_vals = np.stack(shap_vals, axis=-1)  # (n_explain, L, n_targets)
    elif shap_vals.ndim == 2:
        shap_vals = shap_vals[:, :, None]

    out_dir = (resolve_path(cfg, cfg["output_dir"]) / "shap" /
              (args.out or (subtype or "all").replace("-", "").lower()))
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "run_meta.json").write_text(json.dumps({
        "loop_subtype": subtype, "start_char": start_char, "loop_len": L,
        "targets": targets, "n_explain": n_explain, "background_size": args.background_size,
        "seed": args.seed,
    }, indent=2))

    # --- raw values -----------------------------------------------------
    rows = []
    for i, loop in enumerate(explain_loops):
        row = {"loop": loop}
        for ti, t in enumerate(targets):
            for p in range(L):
                row[f"shap_{t}_pos{p+1}"] = float(shap_vals[i, p, ti])
        rows.append(row)
    pd.DataFrame(rows).to_csv(out_dir / "shap_values.csv", index=False)

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
