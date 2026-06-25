"""Does the model think binders are non-selective ('binds one -> binds all')?

The lab's wet-lab data shows little target selectivity (binding is dominated by
a construct-level 'stickiness' factor). This script checks whether the fine-tuned
model reproduces that: it predicts per-target binding probabilities for a set of
sequences and reports

  * the correlation matrix of predicted probabilities across targets, and
  * the conditional co-binding rate  P(call B | call A)  using stored thresholds.

High off-diagonal correlation / co-binding => the model also sees binders as
largely non-selective (a finding consistent with the assay), not a bug.

Usage
-----
    python selectivity.py --config config.yaml --split test
    python selectivity.py --config config.yaml --input seqs.csv --seq-col "Full Seq" --loop-col "Residues"
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from esmc_utils import load_config, resolve_path
from inference import load_trained_model, predict_sequences


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--split", default="test", choices=["train", "val", "test", "all"])
    ap.add_argument("--input", default=None, help="CSV of sequences (overrides --split)")
    ap.add_argument("--seq-col", default=None)
    ap.add_argument("--loop-col", default=None)
    ap.add_argument("--model-dir", default=None)
    ap.add_argument("--batch-size", type=int, default=16)
    args = ap.parse_args()

    cfg = load_config(args.config)
    model_dir = args.model_dir or (resolve_path(cfg, cfg["output_dir"]) / "model")
    model, tokenizer, meta, device = load_trained_model(model_dir)
    targets = meta["targets"]

    # --- gather sequences --------------------------------------------------
    if args.input:
        in_path = (resolve_path(cfg, args.input)
                   if not Path(args.input).is_absolute() else Path(args.input))
        df = pd.read_csv(in_path)
        seq_col = args.seq_col or cfg["data"]["seq_col"]
        loop_col = args.loop_col or cfg["data"]["loop_col"]
        df = df.dropna(subset=[seq_col])
        sequences = df[seq_col].astype(str).str.strip().tolist()
        loops = df[loop_col].tolist() if loop_col in df.columns else None
        source = str(in_path.name)
    else:
        data_dir = resolve_path(cfg, cfg["output_dir"]) / "data"
        parts = ["train", "val", "test"] if args.split == "all" else [args.split]
        df = pd.concat([pd.read_parquet(data_dir / f"{p}.parquet") for p in parts],
                       ignore_index=True)
        sequences, loops = df["sequence"].tolist(), df["loop"].tolist()
        source = f"split={args.split}"

    print(f"--- Predicting {len(sequences)} sequences ({source}) ---")
    preds = predict_sequences(model, tokenizer, meta, sequences, loops=loops,
                              device=device, batch_size=args.batch_size)

    prob_cols = [f"prob_{t}" for t in targets]
    probs = preds[prob_cols]
    pearson = probs.corr(method="pearson")
    spearman = probs.corr(method="spearman")

    # conditional co-binding P(call B | call A) using stored thresholds
    calls = preds[[f"call_{t}" for t in targets]].to_numpy()
    cobind = np.full((len(targets), len(targets)), np.nan)
    for i in range(len(targets)):
        a = calls[:, i] == 1
        na = a.sum()
        for j in range(len(targets)):
            cobind[i, j] = (calls[a, j] == 1).mean() if na else np.nan
    cobind = pd.DataFrame(cobind, index=targets, columns=targets)

    def offdiag_mean(m):
        a = m.to_numpy().astype(float)
        return float(np.nanmean(a[~np.eye(len(a), dtype=bool)]))

    r_mean = offdiag_mean(pearson)

    print("\n=== Predicted-probability correlation (Pearson) ===")
    print(pearson.round(3).to_string())
    print("\n=== Predicted-probability correlation (Spearman) ===")
    print(spearman.round(3).to_string())
    print("\n=== Co-binding rate  P(call=binder for column | call=binder for row) ===")
    print(cobind.round(3).to_string())
    verdict = ("HIGH cross-target correlation -> model sees binders as largely "
               "non-selective (consistent with the assay's 'stickiness')"
               if r_mean >= 0.5 else
               "MODERATE/LOW cross-target correlation -> model retains some "
               "target specificity")
    print(f"\nMean off-diagonal Pearson r = {r_mean:.3f}  ->  {verdict}")

    # --- save --------------------------------------------------------------
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = resolve_path(cfg, cfg["output_dir"]) / "selectivity"
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {"source": source, "n": len(sequences),
               "pearson": pearson.round(4).to_dict(),
               "spearman": spearman.round(4).to_dict(),
               "cobinding_PB_given_A": cobind.round(4).to_dict(),
               "mean_offdiagonal_pearson": round(r_mean, 4)}
    (out_dir / f"selectivity_{ts}.json").write_text(json.dumps(payload, indent=2))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        sns.heatmap(pearson, annot=True, vmin=-1, vmax=1, cmap="RdBu_r", ax=ax[0])
        ax[0].set_title("Predicted-prob correlation (Pearson)")
        sns.heatmap(cobind, annot=True, vmin=0, vmax=1, cmap="Reds", ax=ax[1])
        ax[1].set_title("Co-binding  P(call col | call row)")
        plt.tight_layout()
        png = out_dir / f"selectivity_{ts}.png"
        plt.savefig(png, dpi=150)
        plt.close()
        print(f"Saved heatmaps to {png}")
    except Exception as e:
        print(f"[info] skipped heatmap ({e})")
    print(f"Saved selectivity report to {out_dir / f'selectivity_{ts}.json'}")


if __name__ == "__main__":
    main()
