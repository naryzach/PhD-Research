"""Predict binding for a CSV of TIMP3-variant sequences.

Three modes:
  (default)      score every row against every KNOWN target (wide output,
                 prob_<target>/call_<target> columns) -- same shape as the
                 ESM-C pipeline's predict_csv.py.
  --target-col   a column names a KNOWN target per row (looked up via
                 model_meta.json's target_sequences) -> one prob_bind column,
                 scored against just that row's target.
  --seqB-col     a column holds an EXPLICIT target sequence per row -- works
                 for a target the model never saw a name for.

Usage
-----
    python predict_csv.py --config config.yaml --input seqs.csv --seq-col "Full Seq"
    python predict_csv.py --config config.yaml --input seqs.csv --seq-col "Full Seq" --target-col Target
    python predict_csv.py --config config.yaml --input seqs.csv --seq-col "Full Seq" --seqB-col "Target Seq"
"""
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from inference import load_trained_model, predict_against_known_targets, predict_pairs
from mint_utils import load_config, resolve_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--model-dir", default=None)
    ap.add_argument("--input", required=True)
    ap.add_argument("--seq-col", default=None, help="default: cfg data.seqA_col")
    ap.add_argument("--target-col", default=None, help="column naming a KNOWN target per row")
    ap.add_argument("--seqB-col", default=None, help="column holding an explicit target sequence per row")
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    model_dir = args.model_dir or (resolve_path(cfg, cfg["output_dir"]) / "model")
    model, collate, meta, device = load_trained_model(model_dir)

    in_path = (resolve_path(cfg, args.input) if not Path(args.input).is_absolute()
              else Path(args.input))
    df = pd.read_csv(in_path)
    seq_col = args.seq_col or cfg["data"]["seqA_col"]
    df = df.dropna(subset=[seq_col]).reset_index(drop=True)
    sequences = df[seq_col].astype(str).str.strip().tolist()

    if args.seqB_col:
        seqB = df[args.seqB_col].astype(str).tolist()
        probs = predict_pairs(model, collate, sequences, seqB, device=device, batch_size=args.batch_size)
        out_df = df.copy()
        out_df["prob_bind"] = probs
        out_df["call_bind"] = (probs >= 0.5).astype(int)
    elif args.target_col:
        tnames = df[args.target_col].astype(str).tolist()
        tseqs = meta.get("target_sequences", {})
        missing = sorted(set(t for t in tnames if t not in tseqs))
        if missing:
            raise SystemExit(f"[error] unknown target name(s) in --target-col: {missing}; "
                             f"known targets: {sorted(tseqs)}. Use --seqB-col for novel targets.")
        seqB = [tseqs[t] for t in tnames]
        probs = predict_pairs(model, collate, sequences, seqB, device=device, batch_size=args.batch_size)
        thresholds = meta.get("thresholds", {})
        thr = np.array([float(thresholds.get(t, 0.5)) for t in tnames])
        out_df = df.copy()
        out_df["prob_bind"] = probs
        out_df["call_bind"] = (probs >= thr).astype(int)
    else:
        wide = predict_against_known_targets(model, collate, meta, sequences, device=device,
                                             batch_size=args.batch_size)
        out_df = pd.concat([df.reset_index(drop=True), wide.drop(columns=["sequence"])], axis=1)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = resolve_path(cfg, cfg["output_dir"]) / "predictions"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = Path(args.out) if args.out else out_dir / f"predictions_{ts}.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Saved {len(out_df)} predictions to {out_path}")


if __name__ == "__main__":
    main()
