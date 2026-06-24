"""Prepare the lab binding table for multi-task ESM-C fine-tuning.

What it does
------------
1. Load the long table (one row per sequence x target assay).
2. Resolve conflicting labels (same seq+target observed as both 0 and 1) by
   count-weighted majority vote; drop exact ties.
3. Pivot to ONE row per unique sequence, with a binary label + assay mask +
   read-count for each target. A sequence assayed against only MMP3 has a real
   label for MMP3 and "not assayed" (mask 0) for ADAM17/MMP9.
4. Locate the grafted loop inside each sequence (for loop-region pooling).
5. Split train/val/test by loop *group* so near-identical loops never straddle
   splits (optionally merging edit-distance<=1 loops into one group).
6. Save parquet splits + a manifest (class balance, pos_weights, diagnostics).

Usage
-----
    python data_prep.py --config config.yaml
    python data_prep.py --config config.yaml --limit 1500      # quick smoke data
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from esmc_utils import load_config, resolve_path, set_seed, locate_loop


# --------------------------------------------------------------------------- #
# Conflict resolution + pivot
# --------------------------------------------------------------------------- #
def resolve_and_pivot(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Collapse the long table to one row per sequence with per-target fields."""
    d = cfg["data"]
    seq_col, tgt_col, enc_col = d["seq_col"], d["target_col"], d["encoding_col"]
    loop_col, count_col = d["loop_col"], d["count_col"]
    targets = cfg["targets"]

    df = df.dropna(subset=[seq_col, tgt_col, enc_col]).copy()
    df = df[df[tgt_col].isin(targets)]
    df[enc_col] = df[enc_col].astype(int)
    if count_col in df.columns:
        df[count_col] = pd.to_numeric(df[count_col], errors="coerce").fillna(1).clip(lower=1)
    else:
        df[count_col] = 1

    # Count-weighted vote per (seq, target).
    df["_pos"] = np.where(df[enc_col] == 1, df[count_col], 0.0)
    df["_neg"] = np.where(df[enc_col] == 0, df[count_col], 0.0)
    agg = (
        df.groupby([seq_col, tgt_col])
        .agg(pos=("_pos", "sum"), neg=("_neg", "sum"),
             loop=(loop_col, "first"))
        .reset_index()
    )
    n_conflict = int(((agg["pos"] > 0) & (agg["neg"] > 0)).sum())
    n_tie = int((agg["pos"] == agg["neg"]).sum())
    agg["label"] = np.where(agg["pos"] > agg["neg"], 1,
                            np.where(agg["neg"] > agg["pos"], 0, -1))
    agg["count"] = agg["pos"] + agg["neg"]
    agg = agg[agg["label"] >= 0]  # drop ties

    min_count = int(cfg["preprocess"].get("min_count", 1))
    if min_count > 1:
        agg = agg[agg["count"] >= min_count]

    # Pivot to per-sequence wide form.
    rows = []
    for seq, g in agg.groupby(seq_col):
        rec = {"sequence": seq, "loop": g["loop"].iloc[0]}
        for t in targets:
            sub = g[g[tgt_col] == t]
            if len(sub):
                rec[f"label_{t}"] = int(sub["label"].iloc[0])
                rec[f"mask_{t}"] = 1
                rec[f"count_{t}"] = float(sub["count"].iloc[0])
            else:
                rec[f"label_{t}"] = -100      # sentinel = not assayed
                rec[f"mask_{t}"] = 0
                rec[f"count_{t}"] = 0.0
        rows.append(rec)
    wide = pd.DataFrame(rows)
    wide.attrs["n_conflict"] = n_conflict
    wide.attrs["n_tie"] = n_tie
    return wide


# --------------------------------------------------------------------------- #
# Loop location
# --------------------------------------------------------------------------- #
def add_loop_offsets(wide: pd.DataFrame) -> pd.DataFrame:
    starts, lens, ok = [], [], []
    for seq, loop in zip(wide["sequence"], wide["loop"]):
        s = locate_loop(seq, loop)
        starts.append(s)
        lens.append(len(loop) if isinstance(loop, str) else 0)
        ok.append(s >= 0)
    wide = wide.copy()
    wide["loop_start"] = starts
    wide["loop_len"] = lens
    n_missing = int((~np.array(ok)).sum())
    if n_missing:
        print(f"  [warn] {n_missing} sequences: loop motif not found in sequence "
              f"(loop pooling will fall back to mean for those).")
    return wide


# --------------------------------------------------------------------------- #
# Leakage-safe grouping
# --------------------------------------------------------------------------- #
def assign_groups(loops: pd.Series, cluster_hamming1: bool) -> np.ndarray:
    """Group id per row. Either the exact loop, or an edit-distance<=1 cluster.

    The edit<=1 clustering uses deletion-neighborhood keys: two strings within
    one substitution/insertion/deletion share at least one key, so union-find
    over shared keys merges them. O(N * L), not O(N^2).
    """
    loops = loops.fillna("").astype(str)
    if not cluster_hamming1:
        uniq = {v: i for i, v in enumerate(sorted(loops.unique()))}
        return loops.map(uniq).to_numpy()

    parent: dict = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    key_to_node: dict = {}
    uniq_loops = sorted(set(loops))
    for s in uniq_loops:
        find(s)
        keys = [s] + [s[:i] + s[i + 1:] for i in range(len(s))]  # self + deletions
        for k in keys:
            k = ("d", k)
            if k in key_to_node:
                union(key_to_node[k], s)
            else:
                key_to_node[k] = s
    root_to_gid: dict = {}
    gids = []
    for s in loops:
        r = find(s)
        if r not in root_to_gid:
            root_to_gid[r] = len(root_to_gid)
        gids.append(root_to_gid[r])
    return np.asarray(gids)


def grouped_split(wide: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    from sklearn.model_selection import GroupShuffleSplit

    sp = cfg["split"]
    seed = int(sp.get("seed", 42))
    groups = assign_groups(wide["loop"], bool(sp.get("cluster_hamming1", True)))
    wide = wide.copy()
    wide["group"] = groups

    idx = np.arange(len(wide))
    test_frac = float(sp["test_frac"])
    val_frac = float(sp["val_frac"])

    gss1 = GroupShuffleSplit(n_splits=1, test_size=test_frac, random_state=seed)
    trainval_idx, test_idx = next(gss1.split(idx, groups=groups))

    rel_val = val_frac / (1.0 - test_frac)
    gss2 = GroupShuffleSplit(n_splits=1, test_size=rel_val, random_state=seed)
    tv_groups = groups[trainval_idx]
    tr_rel, val_rel = next(gss2.split(trainval_idx, groups=tv_groups))
    train_idx = trainval_idx[tr_rel]
    val_idx = trainval_idx[val_rel]

    split = np.empty(len(wide), dtype=object)
    split[train_idx] = "train"
    split[val_idx] = "val"
    split[test_idx] = "test"
    wide["split"] = split

    # Leakage assertion: no group shared across splits.
    overlap = (
        wide.groupby("group")["split"].nunique().gt(1).sum()
    )
    assert overlap == 0, f"{overlap} loop-groups leaked across splits!"
    return wide


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
def build_manifest(wide: pd.DataFrame, cfg: dict) -> dict:
    targets = cfg["targets"]
    manifest = {"targets": targets, "n_sequences": int(len(wide)),
                "n_conflicts_resolved": int(wide.attrs.get("n_conflict", 0)),
                "n_ties_dropped": int(wide.attrs.get("n_tie", 0)),
                "splits": {}, "pos_weight": {}}
    for sp in ("train", "val", "test"):
        s = wide[wide["split"] == sp]
        per = {}
        for t in targets:
            m = s[f"mask_{t}"] == 1
            n = int(m.sum())
            npos = int((s.loc[m, f"label_{t}"] == 1).sum())
            per[t] = {"n_assayed": n, "n_pos": npos,
                      "pos_rate": round(npos / n, 4) if n else None}
        manifest["splits"][sp] = {"n_sequences": int(len(s)), "per_target": per}
    # pos_weight from TRAIN only = n_neg / n_pos (for BCE pos_weight).
    tr = wide[wide["split"] == "train"]
    for t in targets:
        m = tr[f"mask_{t}"] == 1
        npos = int((tr.loc[m, f"label_{t}"] == 1).sum())
        nneg = int((tr.loc[m, f"label_{t}"] == 0).sum())
        manifest["pos_weight"][t] = round(nneg / npos, 4) if npos else 1.0
    return manifest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--limit", type=int, default=None,
                    help="Cap number of unique sequences (quick smoke data).")
    args = ap.parse_args()

    cfg = load_config(args.config)
    set_seed(int(cfg["split"].get("seed", 42)))

    csv_path = resolve_path(cfg, cfg["data"]["csv_path"])
    print(f"--- Loading {csv_path} ---")
    df = pd.read_csv(csv_path)
    print(f"  raw rows: {len(df)}")

    wide = resolve_and_pivot(df, cfg)
    print(f"  unique sequences after pivot: {len(wide)} "
          f"(resolved {wide.attrs['n_conflict']} conflicts, "
          f"dropped {wide.attrs['n_tie']} ties)")

    wide = add_loop_offsets(wide)

    if args.limit:
        wide = wide.sample(n=min(args.limit, len(wide)),
                           random_state=int(cfg["split"].get("seed", 42))).reset_index(drop=True)
        print(f"  [smoke] limited to {len(wide)} sequences")

    wide = grouped_split(wide, cfg)

    out_dir = resolve_path(cfg, cfg["output_dir"]) / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    for sp in ("train", "val", "test"):
        s = wide[wide["split"] == sp].reset_index(drop=True)
        s.to_parquet(out_dir / f"{sp}.parquet", index=False)

    manifest = build_manifest(wide, cfg)
    with open(out_dir / "manifest.json", "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)

    print("\n--- Split summary ---")
    for sp in ("train", "val", "test"):
        info = manifest["splits"][sp]
        print(f"{sp:>5}: {info['n_sequences']} seqs")
        for t in cfg["targets"]:
            pt = info["per_target"][t]
            print(f"        {t:<7} assayed={pt['n_assayed']:<6} "
                  f"pos={pt['n_pos']:<6} pos_rate={pt['pos_rate']}")
    print("\npos_weight (train n_neg/n_pos):", manifest["pos_weight"])
    print(f"\nSaved splits + manifest to {out_dir}")


if __name__ == "__main__":
    main()
