"""UMAP visualization of the fine-tuned ESM-C embedding space.

Embeds sequences with the trained backbone (using the model's own pooling, so
the picture reflects what the classifier actually sees), reduces to 2D with
**UMAP** (falls back to PCA if ``umap-learn`` is not installed), and saves a
scatter plot + the 2D coordinates.

Two input modes:
  --split {train,val,test,all}   use the prepared data (color by assayed target
                                 or true binding label)
  --input seqs.csv               arbitrary sequences (color by predicted target)

Usage
-----
    python visualize.py --config config.yaml --split test
    python visualize.py --config config.yaml --split test --color-by binding
    python visualize.py --config config.yaml --input gen.csv --seq-col "Full Seq" \
        --loop-col "Residues"
"""
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from esmc_utils import load_config, locate_loop, resolve_path
from inference import load_trained_model
from model import Collator


# --------------------------------------------------------------------------- #
def reduce_2d(emb, reducer, n_neighbors, min_dist, seed):
    if reducer == "umap":
        try:
            import umap  # umap-learn
            xy = umap.UMAP(n_components=2, n_neighbors=n_neighbors,
                           min_dist=min_dist, random_state=seed).fit_transform(emb)
            return xy, "UMAP"
        except ImportError:
            print("[warn] umap-learn not installed (`pip install umap-learn`); "
                  "falling back to PCA. (See README for the recommended Python env.)")
    from sklearn.decomposition import PCA
    return PCA(n_components=2).fit_transform(emb), "PCA"


@torch.no_grad()
def embed_and_predict(model, tokenizer, meta, sequences, loops, device, batch_size=16):
    """Return (embeddings [N,H], probs [N,T]) in a single pass."""
    items = []
    for seq, loop in zip(sequences, loops if loops is not None else [None] * len(sequences)):
        start = locate_loop(seq, loop) if loop else -1
        items.append({"sequence": str(seq), "labels": [-100] * len(meta["targets"]),
                      "count_weight": [0.0] * len(meta["targets"]),
                      "loop_start": start,
                      "loop_len": len(loop) if isinstance(loop, str) and loop else 0})
    collate = Collator(tokenizer, bos_offset=meta["bos_offset"], max_length=meta["max_length"])
    loader = DataLoader(items, batch_size=batch_size, shuffle=False, collate_fn=collate)
    use_amp = device.type == "cuda"
    embs, logits = [], []
    for batch in loader:
        batch = {k: (v.to(device) if torch.is_tensor(v) else v) for k, v in batch.items()}
        with torch.amp.autocast("cuda", enabled=use_amp):
            pooled = model.embed(batch["input_ids"], batch["attention_mask"],
                                 batch["loop_start"], batch["loop_len"])
            lg = model.head(pooled)
        embs.append(pooled.float().cpu().numpy())
        logits.append(lg.float().cpu().numpy())
    emb = np.vstack(embs)
    probs = 1.0 / (1.0 + np.exp(-np.vstack(logits)))
    return emb, probs


def load_split(cfg, split):
    data_dir = resolve_path(cfg, cfg["output_dir"]) / "data"
    parts = ["train", "val", "test"] if split == "all" else [split]
    df = pd.concat([pd.read_parquet(data_dir / f"{p}.parquet") for p in parts], ignore_index=True)
    return df


def assayed_fields(df, targets):
    a_target, a_binding = [], []
    for _, r in df.iterrows():
        chosen, binding = "none", -1
        for t in targets:
            if int(r[f"mask_{t}"]) == 1:
                chosen, binding = t, int(r[f"label_{t}"])
                break
        a_target.append(chosen)
        a_binding.append("binder" if binding == 1 else ("non-binder" if binding == 0 else "n/a"))
    return a_target, a_binding


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--split", default=None, choices=["train", "val", "test", "all"])
    ap.add_argument("--input", default=None, help="CSV of sequences (alternative to --split)")
    ap.add_argument("--seq-col", default=None)
    ap.add_argument("--loop-col", default=None)
    ap.add_argument("--color-by", default=None,
                    choices=["target", "binding", "top_target"],
                    help="default: 'target' for --split, 'top_target' for --input")
    ap.add_argument("--model-dir", default=None)
    ap.add_argument("--reducer", default=None, choices=["umap", "pca"])
    ap.add_argument("--n-neighbors", type=int, default=None)
    ap.add_argument("--min-dist", type=float, default=None)
    ap.add_argument("--max-points", type=int, default=None,
                    help="subsample for speed/legibility on large sets")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if not args.split and not args.input:
        raise SystemExit("[error] provide --split or --input")

    cfg = load_config(args.config)
    seed = int(cfg["split"].get("seed", 42))
    # CLI overrides config 'visualize' section, which overrides built-in defaults.
    vz = cfg.get("visualize", {}) or {}
    reducer = args.reducer or vz.get("reducer", "umap")
    n_neighbors = args.n_neighbors or int(vz.get("n_neighbors", 15))
    min_dist = args.min_dist if args.min_dist is not None else float(vz.get("min_dist", 0.1))
    max_points = args.max_points or int(vz.get("max_points", 5000))
    targets = cfg["targets"]
    model_dir = args.model_dir or (resolve_path(cfg, cfg["output_dir"]) / "model")
    model, tokenizer, meta, device = load_trained_model(model_dir)

    # --- gather sequences + metadata ---------------------------------------
    if args.split:
        df = load_split(cfg, args.split)
        seq_col, loop_col = "sequence", "loop"
        a_target, a_binding = assayed_fields(df, targets)
        df = df.assign(assayed_target=a_target, binding=a_binding)
        default_color = "target"
    else:
        in_path = (resolve_path(cfg, args.input)
                   if not Path(args.input).is_absolute() else Path(args.input))
        df = pd.read_csv(in_path)
        seq_col = args.seq_col or cfg["data"]["seq_col"]
        loop_col = args.loop_col or cfg["data"]["loop_col"]
        df = df.dropna(subset=[seq_col]).reset_index(drop=True)
        default_color = "top_target"

    if len(df) > max_points:
        df = df.sample(n=max_points, random_state=seed).reset_index(drop=True)
        print(f"[info] subsampled to {max_points} points")

    sequences = df[seq_col].astype(str).str.strip().tolist()
    loops = df[loop_col].tolist() if loop_col in df.columns else None
    if loops is None and meta["pooling"] == "loop":
        print(f"[warn] no loop column '{loop_col}'; loop pooling falls back to mean.")

    # --- embed + predict ---------------------------------------------------
    print(f"--- Embedding {len(sequences)} sequences (pooling={meta['pooling']}) ---")
    emb, probs = embed_and_predict(model, tokenizer, meta, sequences, loops, device)
    df["top_target"] = [targets[i] for i in probs.argmax(axis=1)]
    df["max_prob"] = probs.max(axis=1)

    # --- reduce ------------------------------------------------------------
    print(f"--- Reducing to 2D with {reducer.upper()} ---")
    xy, method = reduce_2d(emb, reducer, n_neighbors, min_dist, seed)
    df["dim1"], df["dim2"] = xy[:, 0], xy[:, 1]

    color_by = args.color_by or default_color
    if color_by == "target":
        hue = "assayed_target"
    elif color_by == "binding":
        hue = "binding"
    else:
        hue = "top_target"
    if hue not in df.columns:
        print(f"[warn] '{hue}' unavailable for this input; coloring by top_target.")
        hue = "top_target"

    # --- plot --------------------------------------------------------------
    plt.figure(figsize=(11, 9))
    sns.scatterplot(data=df, x="dim1", y="dim2", hue=hue, s=28, alpha=0.75,
                    edgecolor="none")
    plt.title(f"{method} of fine-tuned ESM-C embeddings (color = {hue}, pooling = {meta['pooling']})")
    plt.xlabel(f"{method}-1"); plt.ylabel(f"{method}-2")
    plt.legend(title=hue, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = resolve_path(cfg, cfg["output_dir"]) / "visualizations"
    out_dir.mkdir(parents=True, exist_ok=True)
    png = Path(args.out) if args.out else out_dir / f"{method.lower()}_{color_by}_{ts}.png"
    plt.savefig(png, dpi=150)
    plt.close()
    coords = png.with_suffix(".csv")
    keep = [c for c in (seq_col, "assayed_target", "binding", "top_target", "max_prob",
                        "dim1", "dim2") if c in df.columns]
    df[keep].to_csv(coords, index=False)
    print(f"\nSaved plot to {png}\nSaved coordinates to {coords}")


if __name__ == "__main__":
    main()
