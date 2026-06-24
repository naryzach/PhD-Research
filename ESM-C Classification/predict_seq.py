"""Predict per-target binding probabilities for a SINGLE sequence.

Usage
-----
    python predict_seq.py --config config.yaml --seq MRQS...DP
    python predict_seq.py --config config.yaml --seq MRQS...DP --loop AADSIY
    python predict_seq.py --config config.yaml --seq MRQS...DP --json

``--loop`` is the 6-aa grafted loop. Provide it whenever you know it: the model
was trained with loop-region pooling, so giving the loop yields the accurate
prediction. Without it the script falls back to whole-sequence mean pooling.
"""
from __future__ import annotations

import argparse
import json

from esmc_utils import load_config, resolve_path
from inference import load_trained_model, predict_sequences


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--seq", required=True, help="protein sequence")
    ap.add_argument("--loop", default=None, help="grafted loop motif (e.g. AADSIY)")
    ap.add_argument("--model-dir", default=None,
                    help="override model dir (default: <output_dir>/model)")
    ap.add_argument("--json", action="store_true", help="emit JSON")
    args = ap.parse_args()

    cfg = load_config(args.config)
    model_dir = args.model_dir or (resolve_path(cfg, cfg["output_dir"]) / "model")

    model, tokenizer, meta, device = load_trained_model(model_dir)
    df = predict_sequences(model, tokenizer, meta, [args.seq],
                           loops=[args.loop] if args.loop else None, device=device)

    row = df.iloc[0]
    result = {t: {"prob": round(float(row[f"prob_{t}"]), 4),
                  "call": int(row[f"call_{t}"])} for t in meta["targets"]}
    if args.json:
        print(json.dumps({"top_target": row["top_target"],
                          "max_prob": round(float(row["max_prob"]), 4),
                          "per_target": result}, indent=2))
    else:
        print(f"\nSequence length: {len(args.seq)}  loop: {args.loop or '(unknown)'}")
        print(f"{'target':<8} {'P(bind)':>8} {'call':>5}")
        for t in meta["targets"]:
            print(f"{t:<8} {result[t]['prob']:>8.4f} {result[t]['call']:>5}")
        print(f"\nStrongest target: {row['top_target']} (P={row['max_prob']:.4f})")


if __name__ == "__main__":
    main()
