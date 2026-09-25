"""Strict novel-loop evaluation for the trained multirun variants.

``train.py --strict-test`` reports metrics on the test rows that have no training
neighbour within one mutation (``near_train_h1`` False), but the five multirun
variants were trained without that flag, so ``test_report.json`` has
``per_target_novel = null`` (not run; the subset itself is not empty).

This script scores each variant's saved model on its own ``data/test.parquet`` once and
reports metrics for three row sets: all test rows (a check against the stored report),
novel rows (no training neighbour within one mutation), and near rows. Thresholds are the
ones saved in ``model_meta.json`` (chosen on validation); the threshold-free metrics
(ROC-AUC, PR-AUC) do not depend on them.

    python strict_novel_eval.py --variants all
    python strict_novel_eval.py --variants all3_original --limit 200    # quick check
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader

HERE = Path(__file__).resolve().parent
ESMC_DIR = HERE.parent
sys.path.insert(0, str(ESMC_DIR))

from esmc_utils import compute_target_metrics, get_device, get_tokenizer  # noqa: E402
from model import Collator, MultiTaskESMC, SeqDataset  # noqa: E402
from train import evaluate  # noqa: E402

LOCAL = HERE.parents[1] / "Local"
RUNS = LOCAL / "esmc_multirun"
OUT = RUNS / "_cross_sweep_analysis"
VARIANTS = ["all3_original", "abloop_only", "cloop_only", "mmp9_other", "everything_combined"]


def metrics_rows(variant, subset, labels, probs, targets, thresholds, beta):
    rows = []
    for ti, t in enumerate(targets):
        m = labels[:, ti] != -100
        if m.sum() == 0:
            continue
        r = compute_target_metrics(labels[m, ti], probs[m, ti], threshold=thresholds.get(t, 0.5), beta=beta)
        rows.append(dict(variant=variant, subset=subset, target=t, n=int(r["n"]), n_pos=int(r["n_pos"]),
                         pos_rate=r["pos_rate"], roc_auc=r["roc_auc"], pr_auc=r["pr_auc"], mcc=r["mcc"],
                         f1=r["f1"]))
    return rows


def run_variant(variant, limit, batch_size, device):
    mdir = RUNS / variant / "model"
    meta = json.loads((mdir / "model_meta.json").read_text())
    targets = meta["targets"]
    model = MultiTaskESMC(model_id=meta["model_id"], targets=targets, pooling=meta["pooling"],
                          dropout=meta.get("dropout", 0.1), pos_weight=None)
    model.load_state_dict(torch.load(mdir / "model_state.pt", map_location="cpu"))
    model.to(device).eval()
    tok = get_tokenizer(meta["model_id"])
    collate = Collator(tok, bos_offset=meta["bos_offset"], max_length=meta["max_length"])

    ds = SeqDataset(RUNS / variant / "data" / "test.parquet", targets, count_weighting="log")
    df = ds.df.reset_index(drop=True)
    if limit:
        df = df.sample(n=min(limit, len(df)), random_state=0).reset_index(drop=True)
    near = df["near_train_h1"].astype(bool).to_numpy()
    loader = DataLoader(SeqDataset(df, targets, count_weighting="log"), batch_size=batch_size,
                        shuffle=False, collate_fn=collate)
    t0 = time.time()
    _, probs, labels = evaluate(model, loader, device, targets, thresholds=meta["thresholds"])
    print(f"{variant}: scored {len(df)} test rows in {time.time() - t0:.0f}s", flush=True)
    beta = 0.5
    rows = []
    for name, sel in (("all", np.ones(len(df), bool)), ("novel", ~near), ("near", near)):
        rows += metrics_rows(variant, name, labels[sel], probs[sel], targets, meta["thresholds"], beta)
    del model
    torch.cuda.empty_cache()
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", default="all")
    ap.add_argument("--limit", type=int, default=0, help="score a random subset of test rows (quick check)")
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--out", default=str(OUT / "strict_novel_eval_2026-09-24.csv"))
    args = ap.parse_args()
    variants = VARIANTS if args.variants == "all" else args.variants.split(",")
    device = get_device()
    print("device:", device, flush=True)
    allrows = []
    for v in variants:
        allrows += run_variant(v, args.limit, args.batch_size, device)
        pd.DataFrame(allrows).to_csv(args.out, index=False)      # incremental save
    df = pd.DataFrame(allrows)
    pd.set_option("display.width", 200)
    print(df.round(3).to_string(index=False))
    print("saved", args.out)


if __name__ == "__main__":
    main()
