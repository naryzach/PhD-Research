"""Compare pooling strategies (loop vs mean vs cls) under identical conditions.

Trains a fresh model for each pooling mode on the SAME prepared splits and
reports held-out-test metrics side by side, so you can confirm that loop-region
pooling actually helps (the hypothesis: the 6-aa loop is the whole signal, so
pooling only the loop beats mean/cls).

Run `data_prep.py` first. Each mode is saved to <output_dir>/model_<mode>/.

Usage
-----
    python compare_pooling.py --config config.yaml                 # full
    python compare_pooling.py --config config.yaml --smoke          # fast
    python compare_pooling.py --config config.yaml --modes loop mean
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime

import pandas as pd

from esmc_utils import load_config, mean_metric, resolve_path
from train import train_and_evaluate


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--modes", nargs="+", default=["loop", "mean", "cls"],
                    choices=["loop", "mean", "cls"])
    args = ap.parse_args()

    cfg = load_config(args.config)
    targets = cfg["targets"]
    sel_metric = cfg["train"].get("metric", "pr_auc")

    rows = []
    results = {}
    for mode in args.modes:
        print(f"\n{'#'*70}\n# POOLING = {mode}\n{'#'*70}")
        res = train_and_evaluate(cfg, smoke=args.smoke, pooling=mode,
                                 out_name=f"model_{mode}")
        results[mode] = res
        per = res["per_test"]
        row = {"pooling": mode, "best_val": round(res["best_val"], 4),
               f"test_mean_{sel_metric}": round(mean_metric(per, sel_metric), 4),
               "test_mean_mcc": round(mean_metric(per, "mcc"), 4)}
        for t in targets:
            row[f"{t}_pr_auc"] = round(per.get(t, {}).get("pr_auc", float("nan")), 4)
        rows.append(row)

    table = pd.DataFrame(rows)
    print(f"\n{'='*70}\nPOOLING COMPARISON (held-out test)\n{'='*70}")
    print(table.to_string(index=False))

    best = max(results, key=lambda m: mean_metric(results[m]["per_test"], sel_metric))
    print(f"\nBest pooling by mean test {sel_metric}: {best}")

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = resolve_path(cfg, cfg["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"pooling_comparison_{ts}.csv"
    table.to_csv(csv_path, index=False)
    (out_dir / f"pooling_comparison_{ts}.json").write_text(
        json.dumps({m: results[m]["per_test"] for m in results}, indent=2))
    print(f"Saved comparison to {csv_path}")


if __name__ == "__main__":
    main()
