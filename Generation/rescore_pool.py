#!/usr/bin/env python3
"""
rescore_pool.py — recover ESMFold2 scores for designs that were GENERATED but
left UNSCORED (e.g. by a GPU-contention CUDA-OOM during a run). Re-folds ONLY the
unscored designs with ESMFold2 — the cheap stage — WITHOUT re-running RFd3/
LigandMPNN, computes the calibrated interface features, fills the affected
round_summary.csv rows in place (each backed up to <name>.bak), and rebuilds the
per-target HOF from the now-complete pool.

Use this after a run where round summaries show many rows with esm_plddt = NaN
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
import subprocess
import sys
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
    ap.add_argument("--chunk-size", type=int, default=300, metavar="N",
                    help="Fold this many designs per scorer invocation, writing the recovered "
                         "rows back to the round summaries after each chunk. A crash then "
                         "costs at most one chunk, and re-running picks up where it stopped.")
    ap.add_argument("--timeout", type=int, default=None, metavar="SEC",
                    help="Per-chunk scorer timeout (default: 20 s/design, min 3600).")
    ap.add_argument("--skip-preflight", action="store_true",
                    help="Skip the ESMFold2 fold test. Not recommended — the test costs "
                         "~2 min against a multi-day job.")
    args = ap.parse_args()
    ir.ESMFOLD2_GPUS = args.esmfold2_gpus
    # --resume lets a shard killed mid-chunk restart from its partial CSV.
    ir.ESMFOLD2_EXTRA_ARGS = ["--resume"]
    ir.ESMFOLD2_TIMEOUT_S = args.timeout or max(3600, 20 * args.chunk_size)

    refiner = IterativeRefiner(active_targets=args.targets, active_loops=args.loops)

    csvs = sorted(OUT_BASE.glob("it_*/round_summary.csv"))
    if not csvs:
        print(f"No round_summary.csv under {OUT_BASE}; nothing to do.")
        return

    frames = {c: pd.read_csv(c) for c in csvs}

    # Collect unscored rows (no esm_plddt) for the requested targets.
    to_score = {}   # design_id -> row dict for the scorer
    per_target = {}
    for c, df in frames.items():
        # Key on esm_plddt, NOT esm_iptm: esm_plddt (with the interface features) is
        # what calc_composite's ESMFold2 branch actually reads, and recover_from_cifs.py
        # fills it from saved structures while leaving esm_iptm blank by design. Keying
        # on esm_iptm would re-fold every CIF-recovered design for nothing.
        if "esm_plddt" in df.columns:
            need = df["esm_plddt"].isna()
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

    # Prove ESMFold2 can fold BEFORE committing a multi-day job to it. `import esm`
    # is not evidence — the foundry env passes that and still fails at fold time.
    if not args.skip_preflight:
        print("Preflight: folding a test complex through the real scorer path (~1-2 min)...",
              flush=True)
        rc = subprocess.call([sys.executable, str(Path(ir.__file__).parent / "esmfold2_smoketest.py")])
        if rc != 0:
            sys.exit("\nABORTING: ESMFold2 cannot fold in this environment (see output above).\n"
                     "  Fix the env — activate the esmfold2 env, or set ESMFOLD2_PYTHON, or\n"
                     "  export DISABLE_CUEQUIVARIANCE=1 — then re-run. Nothing was modified.")
        print("Preflight OK.\n", flush=True)

    # Re-fold with ESMFold2 (reuses the pipeline's GPU-sharded scorer), in chunks so
    # that a hard kill costs one chunk instead of the whole pool — the failure mode
    # that silently wiped the Aug-3 run (structures on disk, scores never written).
    items = list(to_score.values())
    chunks = [items[i:i + args.chunk_size] for i in range(0, len(items), args.chunk_size)]
    print(f"Scoring in {len(chunks)} chunk(s) of up to {args.chunk_size} "
          f"(timeout {ir.ESMFOLD2_TIMEOUT_S}s/chunk).")

    n_filled, n_empty_streak = 0, 0
    for ci, chunk in enumerate(chunks, 1):
        out_dir = OUT_BASE / "rescore" / f"chunk_{ci:04d}"
        out_dir.mkdir(parents=True, exist_ok=True)
        scores = refiner._score_sequences_esmfold2(chunk, out_dir)
        print(f"[chunk {ci}/{len(chunks)}] ESMFold2 returned {len(scores)}/{len(chunk)} scores.",
              flush=True)
        if not scores:
            # Don't grind through 45 chunks of nothing the way the Aug-3 run ground
            # through 37 iterations of nothing.
            n_empty_streak += 1
            if n_empty_streak >= 2:
                sys.exit(f"\nABORTING: {n_empty_streak} consecutive chunks scored 0 designs — "
                         "the ranker is dead, not flaky.\n"
                         f"  Check:  python {Path(ir.__file__).parent / 'esmfold2_smoketest.py'}\n"
                         "  Rows filled so far are saved; re-running skips them.")
            print(f"[chunk {ci}/{len(chunks)}] no scores — retrying with the next chunk.")
            continue
        n_empty_streak = 0

        # Add the calibrated interface features from each predicted complex CIF.
        feats = {did: {**m, **_esm_iface_feats(m.get("esm_cif"))} for did, m in scores.items()}

        # Fill the affected round_summary rows in place (with .bak backups).
        for c, df in frames.items():
            changed = False
            for i, r in df.iterrows():
                did = r.get("design_id")
                already = ("esm_plddt" in df.columns) and pd.notna(r.get("esm_plddt"))
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
                if not Path(str(c) + ".bak").exists():
                    shutil.copy2(c, str(c) + ".bak")   # one backup of the pre-rescore state
                df.to_csv(c, index=False)
        print(f"[chunk {ci}/{len(chunks)}] round summaries updated "
              f"({n_filled} rows filled so far).", flush=True)

    if not n_filled:
        print("No rows filled — round summaries unchanged. Check ESMFOLD2_PYTHON / GPU.")
        return
    print(f"Filled {n_filled} rows across {sum(1 for c in frames)} round summaries (backups: *.bak).")

    # Rebuild the per-target HOF from the now-complete scored pool.
    all_scored = []
    for c in csvs:
        df = pd.read_csv(c)
        df = df[pd.to_numeric(df.get("esm_plddt"), errors="coerce").notna()]
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
