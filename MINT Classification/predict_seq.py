"""Predict binding of a single TIMP3-variant sequence against a target.

Usage
-----
    python predict_seq.py --config config.yaml --seq <SEQ>                    # vs every known target
    python predict_seq.py --config config.yaml --seq <SEQ> --target MMP9      # vs one known target
    python predict_seq.py --config config.yaml --seq <SEQ> --seqB <TARGET_SEQ> # arbitrary target sequence
"""
from __future__ import annotations

import argparse

from inference import load_trained_model, predict_against_known_targets, predict_pairs
from mint_utils import load_config, resolve_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--model-dir", default=None)
    ap.add_argument("--seq", required=True, help="TIMP3-variant full sequence (chain A)")
    ap.add_argument("--target", default=None,
                    help="known target name (e.g. MMP9); default: score against every known target")
    ap.add_argument("--seqB", default=None,
                    help="explicit target sequence (chain B) -- overrides --target; works for a "
                         "target the model never saw a name for")
    args = ap.parse_args()

    cfg = load_config(args.config)
    model_dir = args.model_dir or (resolve_path(cfg, cfg["output_dir"]) / "model")
    model, collate, meta, device = load_trained_model(model_dir)

    if args.seqB:
        prob = predict_pairs(model, collate, [args.seq], [args.seqB], device=device)[0]
        print(f"P(bind) vs supplied target sequence: {prob:.4f}")
        return

    targets = [args.target] if args.target else None
    df = predict_against_known_targets(model, collate, meta, [args.seq], targets=targets, device=device)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
