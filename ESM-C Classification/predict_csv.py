"""Predict per-target binding probabilities for a CSV of sequences.

Usage
-----
    python predict_csv.py --config config.yaml --input my_seqs.csv
    python predict_csv.py --config config.yaml --input my_seqs.csv \
        --seq-col "Full Seq" --loop-col "Residues" --out results.csv

The input CSV needs a sequence column (default "Full Seq"). If it also has a
loop-motif column (default "Residues"), it is used for accurate loop-region
pooling; otherwise the script falls back to mean pooling with a warning.
"""
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd

from esmc_utils import load_config, resolve_path
from inference import load_trained_model, predict_sequences


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--input", required=True, help="input CSV of sequences")
    ap.add_argument("--seq-col", default=None, help="sequence column (default from config)")
    ap.add_argument("--loop-col", default=None, help="loop-motif column (default from config)")
    ap.add_argument("--out", default=None, help="output CSV path")
    ap.add_argument("--model-dir", default=None)
    ap.add_argument("--batch-size", type=int, default=16)
    args = ap.parse_args()

    cfg = load_config(args.config)
    seq_col = args.seq_col or cfg["data"]["seq_col"]
    loop_col = args.loop_col or cfg["data"]["loop_col"]
    model_dir = args.model_dir or (resolve_path(cfg, cfg["output_dir"]) / "model")

    in_path = resolve_path(cfg, args.input) if not Path(args.input).is_absolute() else Path(args.input)
    df_in = pd.read_csv(in_path)
    if seq_col not in df_in.columns:
        raise SystemExit(f"[error] sequence column '{seq_col}' not in {list(df_in.columns)}")
    df_in = df_in.dropna(subset=[seq_col]).copy()
    df_in[seq_col] = df_in[seq_col].astype(str).str.strip()
    df_in = df_in[df_in[seq_col] != ""]

    sequences = df_in[seq_col].tolist()
    loops = df_in[loop_col].tolist() if loop_col in df_in.columns else None
    if loops is None:
        print(f"[info] no loop column '{loop_col}' found; mean pooling will be used.")

    print(f"--- Loading model from {model_dir} ---")
    model, tokenizer, meta, device = load_trained_model(model_dir)
    print(f"--- Predicting {len(sequences)} sequences (targets: {meta['targets']}) ---")
    preds = predict_sequences(model, tokenizer, meta, sequences, loops=loops,
                              device=device, batch_size=args.batch_size)

    # merge prediction columns back onto the original rows (keep their metadata)
    pred_cols = [c for c in preds.columns if c not in ("sequence", "loop")]
    out = pd.concat([df_in.reset_index(drop=True), preds[pred_cols].reset_index(drop=True)], axis=1)

    if args.out:
        out_path = Path(args.out)
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = resolve_path(cfg, cfg["output_dir"]) / "predictions"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"predictions_{ts}.csv"
    out.to_csv(out_path, index=False)
    print(f"\nSaved {len(out)} predictions to {out_path}")

    # quick summary
    for t in meta["targets"]:
        print(f"  {t:<7} predicted binders: {int(out[f'call_{t}'].sum())}/{len(out)}")


if __name__ == "__main__":
    main()
