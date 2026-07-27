#!/usr/bin/env python3
"""
rescore_pool.py — recover ESMFold2 scores for designs that were GENERATED but
left UNSCORED (e.g. by a GPU-contention CUDA-OOM during a run). Re-folds ONLY the
unscored designs with ESMFold2 — the cheap stage — WITHOUT re-running RFd3/
LigandMPNN, computes the calibrated interface features, fills the affected
round_summary.csv rows in place (each backed up to <name>.bak), and rebuilds the
per-target HOF from the now-complete pool.

Use this after a run where round summaries show many rows with esm_iptm = NaN
(composite_score = 0). Typical case: ADAM10/ADAM17 starved by OOM while the MMPs
scored fine.

RUN IN THE GENERATION ENV (foundry + ESMFOLD2_PYTHON), GPUs FREE:
    export ESMFOLD2_PYTHON=$(command -v python)          # if ESMFold2 is in foundry
    python Generation/rescore_pool.py                     # recover all targets
    python Generation/rescore_pool.py --targets ADAM10 ADAM17   # just the starved ones
    python Generation/rescore_pool.py --esmfold2-gpus auto

It reuses the pipeline's own ESMFold2 scorer and composite, so recovered rows are
identical in format to a clean run. Idempotent: already-scored rows are skipped.
"""
import argparse
import shutil
from pathlib import Path

import pandas as pd
import numpy as np

import iterative_refinement as ir
from iterative_refinement import (
    IterativeRefiner, OUT_BASE, _esm_iface_feats, calc_composite,
)


def main():
    ap = argparse.ArgumentParser(description="Recover ESMFold2 scores for unscored designs.")
    ap.add_argument("--targets", nargs="+", default=["MMP2", "MMP9", "ADAM10", "ADAM17"])
    ap.add_argument("--loops", nargs="+", default=["AB", "C", "EF"])
    ap.add_argument("--esmfold2-gpus", default="auto", metavar="N|auto")
    ap.add_argument("--dry-run", action="store_true", help="Report the unscored counts and exit.")
    args = ap.parse_args()
    ir.ESMFOLD2_GPUS = args.esmfold2_gpus

    refiner = IterativeRefiner(active_targets=args.targets, active_loops=args.loops)

    csvs = sorted(OUT_BASE.glob("it_*/round_summary.csv"))
    if not csvs:
        print(f"No round_summary.csv under {OUT_BASE}; nothing to do.")
        return

    frames = {c: pd.read_csv(c) for c in csvs}

    # Collect unscored rows (esm_iptm NaN) for the requested targets.
    to_score = {}   # design_id -> row dict for the scorer
    per_target = {}
    for c, df in frames.items():
        if "esm_iptm" in df.columns:
            need = df["esm_iptm"].isna()
        else:
            need = pd.Series(True, index=df.index)
        need &= df["target_name"].isin(args.targets)
        for _, r in df[need].iterrows():
            tseq = refiner.target_seqs.get(r["target_name"], "")
            did = r.get("design_id")
            if isinstance(r.get("full_seq"), str) and tseq and did and did not in to_score:
                to_score[did] = {"design_id": did, "target_name": r["target_name"],
                                 "full_seq": r["full_seq"], "target_seq": tseq}
                per_target[r["target_name"]] = per_target.get(r["target_name"], 0) + 1

    print(f"Unscored designs to recover: {len(to_score)}  {per_target}")
    if args.dry_run or not to_score:
        return

    # Re-fold with ESMFold2 (reuses the pipeline's GPU-sharded scorer).
    out_dir = OUT_BASE / "rescore"
    out_dir.mkdir(parents=True, exist_ok=True)
    scores = refiner._score_sequences_esmfold2(list(to_score.values()), out_dir)
    print(f"ESMFold2 returned {len(scores)} scores.")
    if not scores:
        print("No scores produced (check ESMFOLD2_PYTHON / GPU). Round summaries untouched.")
        return

    # Add the calibrated interface features from each predicted complex CIF.
    feats = {did: {**m, **_esm_iface_feats(m.get("esm_cif"))} for did, m in scores.items()}

    # Fill the affected round_summary rows in place (with .bak backups).
    n_filled = 0
    for c, df in frames.items():
        changed = False
        for i, r in df.iterrows():
            did = r.get("design_id")
            already = ("esm_iptm" in df.columns) and pd.notna(r.get("esm_iptm"))
            if did in feats and not already:
                m = feats[did]
                merged = {**r.to_dict(), **m, "source": "ESMFold2"}
                for k, v in m.items():
                    df.at[i, k] = v
                df.at[i, "source"] = "ESMFold2"
                df.at[i, "promising"] = bool(ir._safe(m.get("esm_iptm")) >= ir.IPTM_PROMISING)
                df.at[i, "composite_score"] = calc_composite(merged)
                changed = True
                n_filled += 1
        if changed:
            shutil.copy2(c, str(c) + ".bak")
            df.to_csv(c, index=False)
    print(f"Filled {n_filled} rows across {sum(1 for c in frames)} round summaries (backups: *.bak).")

    # Rebuild the per-target HOF from the now-complete scored pool.
    all_scored = []
    for c in csvs:
        df = pd.read_csv(c)
        df = df[pd.to_numeric(df.get("esm_iptm"), errors="coerce").notna()]
        all_scored.extend(df.to_dict("records"))
    for t in refiner.active_targets:
        refiner.state["hof"][t] = []
    refiner.update_hof(all_scored)
    refiner._save_state()
    hof_counts = {t: len(refiner.state["hof"].get(t, [])) for t in refiner.active_targets}
    print(f"HOF rebuilt from {len(all_scored)} scored designs. HOF sizes: {hof_counts}")
    print("Done. Re-run select_binders / stratified export to use the recovered pool.")


if __name__ == "__main__":
    main()
