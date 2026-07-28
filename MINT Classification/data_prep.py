"""Prepare the lab binding table for MINT pair fine-tuning.

Unlike the ESM-C pipeline, there is no pivot-to-multi-task step here: MINT
takes a PAIR (TIMP3-variant sequence, target sequence) per example, and the
target is just a column value that comes along with each row -- not a
separate output head. So one assayed (sequence, target) row IS one training
example. What's left to do:

1. Load the long table (one row per sequence x target assay).
2. Resolve conflicting labels (same seq+target observed as both 0 and 1) by
   count-weighted majority vote; drop exact ties. (Identical logic to the
   ESM-C pipeline, just not followed by a pivot.)
3. Split train/val/test by loop *group* so near-identical loops never
   straddle splits (same leakage-safe grouping as the ESM-C pipeline).
4. Save parquet splits + a manifest (per-target-name class balance, overall
   pos_weight, diagnostics).

Usage
-----
    python data_prep.py --config config.yaml
    python data_prep.py --config config.yaml --limit 300      # quick smoke data
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from mint_utils import load_config, resolve_path, set_seed


# --------------------------------------------------------------------------- #
# Conflict resolution (long format -- one row per (seqA, target) pair)
# --------------------------------------------------------------------------- #
def resolve_conflicts(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    d = cfg["data"]
    a_col, b_col = d["seqA_col"], d["seqB_col"]
    tgt_col, enc_col = d["target_name_col"], d["encoding_col"]
    loop_col, count_col = d["loop_col"], d["count_col"]

    df = df.dropna(subset=[a_col, b_col, tgt_col, enc_col]).copy()
    df[enc_col] = df[enc_col].astype(int)
    if count_col in df.columns:
        df[count_col] = pd.to_numeric(df[count_col], errors="coerce").fillna(1).clip(lower=1)
    else:
        df[count_col] = 1

    df["_pos"] = np.where(df[enc_col] == 1, df[count_col], 0.0)
    df["_neg"] = np.where(df[enc_col] == 0, df[count_col], 0.0)
    agg = (
        df.groupby([a_col, tgt_col])
        .agg(pos=("_pos", "sum"), neg=("_neg", "sum"),
             seqB=(b_col, "first"), loop=(loop_col, "first"))
        .reset_index()
        .rename(columns={a_col: "seqA", tgt_col: "target_name"})
    )
    n_conflict = int(((agg["pos"] > 0) & (agg["neg"] > 0)).sum())
    n_tie = int((agg["pos"] == agg["neg"]).sum())
    agg["label"] = np.where(agg["pos"] > agg["neg"], 1,
                            np.where(agg["neg"] > agg["pos"], 0, -1))
    agg["count"] = agg["pos"] + agg["neg"]
    agg = agg[agg["label"] >= 0].reset_index(drop=True)  # drop ties

    min_count = int(cfg["preprocess"].get("min_count", 1))
    if min_count > 1:
        agg = agg[agg["count"] >= min_count].reset_index(drop=True)

    agg.attrs["n_conflict"] = n_conflict
    agg.attrs["n_tie"] = n_tie
    return agg[["seqA", "seqB", "target_name", "loop", "label", "count"]]


# --------------------------------------------------------------------------- #
# Leakage-safe grouping (identical to the ESM-C pipeline: group by exact loop
# motif, or an edit-distance<=1 cluster, so near-duplicate grafts don't
# straddle train/val/test)
# --------------------------------------------------------------------------- #
def assign_groups(loops: pd.Series, cluster_hamming1: bool) -> np.ndarray:
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
    for s in sorted(set(loops)):
        find(s)
        for k in [s] + [s[:i] + s[i + 1:] for i in range(len(s))]:
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


def grouped_split(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    sp = cfg["split"]
    seed = int(sp.get("seed", 42))
    cluster = bool(sp.get("cluster_hamming1", False))
    groups = assign_groups(df["loop"], cluster)
    df = df.copy()
    df["group"] = groups

    n = len(df)
    sizes = df.groupby("group").size()
    biggest = int(sizes.max())
    if biggest > 0.5 * n:
        print(f"  [warn] largest loop-group is {biggest} rows ({biggest/n:.1%} of data) -- "
              f"the edit-distance<=1 loop graph has percolated; set split.cluster_hamming1: false.")

    test_target, val_target = float(sp["test_frac"]) * n, float(sp["val_frac"]) * n
    gids = sizes.index.to_numpy().copy()
    np.random.default_rng(seed).shuffle(gids)

    assign, c_test, c_val = {}, 0, 0
    for g in gids:
        s = int(sizes[g])
        if c_test < test_target:
            assign[g], c_test = "test", c_test + s
        elif c_val < val_target:
            assign[g], c_val = "val", c_val + s
        else:
            assign[g] = "train"
    df["split"] = df["group"].map(assign)

    overlap = df.groupby("group")["split"].nunique().gt(1).sum()
    assert overlap == 0, f"{overlap} loop-groups leaked across splits!"
    return df


def hamming1_flags(df: pd.DataFrame) -> pd.Series:
    train_keys = set()
    for s in df.loc[df["split"] == "train", "loop"].fillna("").astype(str):
        train_keys.add(("s", s))
        for i in range(len(s)):
            train_keys.add(("d", s[:i] + s[i + 1:]))

    def exposed(s: str) -> bool:
        if ("s", s) in train_keys:
            return True
        return any(("d", s[:i] + s[i + 1:]) in train_keys for i in range(len(s)))

    flags = []
    for split, loop in zip(df["split"], df["loop"].fillna("").astype(str)):
        flags.append(False if split == "train" else exposed(loop))
    return pd.Series(flags, index=df.index)


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
def build_manifest(df: pd.DataFrame, n_conflict: int, n_tie: int) -> dict:
    targets = sorted(df["target_name"].unique())
    manifest = {"targets": targets, "n_rows": int(len(df)),
                "n_conflicts_resolved": n_conflict, "n_ties_dropped": n_tie,
                "splits": {}, "pos_weight": {}}
    for sp in ("train", "val", "test"):
        s = df[df["split"] == sp]
        per = {}
        for t in targets:
            st = s[s["target_name"] == t]
            n, npos = len(st), int((st["label"] == 1).sum())
            per[t] = {"n": n, "n_pos": npos, "pos_rate": round(npos / n, 4) if n else None}
        manifest["splits"][sp] = {"n_rows": int(len(s)), "per_target": per}
    tr = df[df["split"] == "train"]
    npos, nneg = int((tr["label"] == 1).sum()), int((tr["label"] == 0).sum())
    manifest["pos_weight"]["overall"] = round(nneg / npos, 4) if npos else 1.0
    for t in targets:
        st = tr[tr["target_name"] == t]
        npos_t, nneg_t = int((st["label"] == 1).sum()), int((st["label"] == 0).sum())
        manifest["pos_weight"][t] = round(nneg_t / npos_t, 4) if npos_t else 1.0
    # target_name -> canonical target sequence, so predict_seq.py/predict_csv.py can
    # take --target-name instead of the caller having to paste the full protease sequence.
    manifest["target_sequences"] = {t: str(df.loc[df["target_name"] == t, "seqB"].iloc[0])
                                    for t in targets}
    return manifest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--limit", type=int, default=None,
                    help="Cap number of rows (quick smoke data).")
    args = ap.parse_args()

    cfg = load_config(args.config)
    set_seed(int(cfg["split"].get("seed", 42)))

    csv_path = resolve_path(cfg, cfg["data"]["csv_path"])
    print(f"--- Loading {csv_path} ---")
    raw = pd.read_csv(csv_path)
    print(f"  raw rows: {len(raw)}")

    df = resolve_conflicts(raw, cfg)
    n_conflict, n_tie = df.attrs.get("n_conflict", 0), df.attrs.get("n_tie", 0)
    print(f"  (seqA,target) pairs after conflict resolution: {len(df)} "
          f"(resolved {n_conflict} conflicts, dropped {n_tie} ties)")

    if args.limit:
        df = df.sample(n=min(args.limit, len(df)),
                       random_state=int(cfg["split"].get("seed", 42))).reset_index(drop=True)
        print(f"  [smoke] limited to {len(df)} rows")

    df = grouped_split(df, cfg)
    df["near_train_h1"] = hamming1_flags(df)

    out_dir = resolve_path(cfg, cfg["output_dir"]) / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    for sp in ("train", "val", "test"):
        s = df[df["split"] == sp].reset_index(drop=True)
        s.to_parquet(out_dir / f"{sp}.parquet", index=False)

    manifest = build_manifest(df, n_conflict, n_tie)
    manifest["hamming1_exposure"] = {
        sp: round(float(df.loc[df["split"] == sp, "near_train_h1"].mean()), 4)
        for sp in ("val", "test")
    }
    with open(out_dir / "manifest.json", "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)

    print("\n--- Split summary ---")
    for sp in ("train", "val", "test"):
        info = manifest["splits"][sp]
        print(f"{sp:>5}: {info['n_rows']} rows")
        for t in manifest["targets"]:
            pt = info["per_target"][t]
            print(f"        {t:<7} n={pt['n']:<6} pos={pt['n_pos']:<6} pos_rate={pt['pos_rate']}")
    print("\npos_weight (train n_neg/n_pos):", manifest["pos_weight"])
    exp = manifest["hamming1_exposure"]
    print(f"Hamming-1 exposure (val/test loops within 1 edit of a train loop): {exp}")
    print(f"\nSaved splits + manifest to {out_dir}")


if __name__ == "__main__":
    main()
