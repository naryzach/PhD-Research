#!/usr/bin/env python3
"""
partial_t_ladder.py — find the noise level at which a parent's interface survives.

The it_54-63 conformation campaign seeded ~75% of every round from HOF parents
averaging sv_pdockq 0.46-0.55, yet the round means did not move off the
sequence-only baseline (0.34-0.48) on any of the four targets. The beam was
healthy (10 seeds/target/round on disk) and partial diffusion ran, so the
failure is not plumbing: partial_t=8.0 A nominal (~4 A of actual loop movement,
per partial_diffusion_smoketest) re-rolls the interface faster than selection
can exploit it. The children are random draws with extra steps.

This sweeps partial_t over the SAME seeds and scores the children through the
real ranker, so "did the child inherit the parent" becomes a measured number:

    transmission = (mean child pdockq - population mean) /
                   (mean parent pdockq - population mean)

  ~1.0  children hold the parent's interface   (t too low -> no exploration;
                                                check the loop RMSD column)
  ~0.0  children are random draws              (t too high -> what it_54-63 did)

Pick the largest t that still transmits; that is the most exploration the
parent's information can pay for.

    python Generation/partial_t_ladder.py --target MMP2
    python Generation/partial_t_ladder.py --target ADAM17 --t 1,2,3,5 --seeds 4

Needs a free GPU (it competes with a live generation job) and runs ~1 h at the
defaults: 5 noise levels x 3 seeds x 4 children = 60 ESMFold2 folds.
"""
import argparse
import glob
import os
import re
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
os.environ.setdefault("DISABLE_CUEQUIVARIANCE", "1")

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))


def population_mean(out_base, target, metric, window=10):
    """
    Baseline a child must beat to have inherited anything: what the pipeline
    produces right now, i.e. the last `window` iterations rather than the whole
    campaign. Averaging all of history drags the baseline down (MMP2: 0.286
    all-time vs 0.336 over it_46+) and inflates transmission for free.
    """
    per = {}
    for f in glob.glob(str(out_base / "it_*" / "round_summary.csv")):
        d = pd.read_csv(f)
        if metric not in d.columns or "esm_plddt" not in d.columns:
            continue
        d = d[(d.esm_plddt.notna()) & (d.target_name == target)]
        v = [x for x in d[metric].tolist() if x == x]
        if v:
            per[int(re.search(r"it_(\d+)", f).group(1))] = v
    if not per:
        return float("nan")
    keep = sorted(per)[-window:] if window else sorted(per)
    vals = [x for i in keep for x in per[i]]
    return float(np.mean(vals))


def main() -> int:
    import conformation_seed as cseed
    import iterative_refinement as ir
    from biotite.structure import superimpose

    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="MMP2")
    ap.add_argument("--t", default="1,2,3,5,8",
                    help="comma-separated partial_t values (Angstroms)")
    ap.add_argument("--seeds", type=int, default=3, help="beam seeds to test")
    ap.add_argument("--children", type=int, default=4, help="children per seed per t")
    ap.add_argument("--metric", default="sv_pdockq")
    ap.add_argument("--pop-window", type=int, default=10,
                    help="iterations of history defining the baseline (0 = all)")
    ap.add_argument("--no-sv", action="store_true",
                    help="skip the SV battery (sv_* metrics unavailable)")
    ap.add_argument("--out-base", default=None)
    ap.add_argument("--work", default=None, help="scratch dir (default: <out_base>/ladder)")
    args = ap.parse_args()

    tvals = [float(x) for x in args.t.split(",")]
    # SV_BATTERY is a module global that iterative_refinement.main() flips from its
    # --sv-battery flag. Importing the module leaves it False, so sv_pdockq is never
    # computed and every child scores as "nothing" -- turn it on the same way.
    if not args.no_sv:
        ir.SV_BATTERY = True
    out_base = Path(args.out_base) if args.out_base else ir.OUT_BASE
    work = Path(args.work) if args.work else out_base / "ladder" / args.target
    work.mkdir(parents=True, exist_ok=True)
    loops = ["AB", "C", "EF"]
    B = ir.DESIGN_BINDER_CHAIN

    pop = population_mean(out_base, args.target, args.metric, args.pop_window)
    print(f"target {args.target} | population mean {args.metric} = {pop:.3f} "
          f"(last {args.pop_window} iterations)")

    # Same beam the live run uses, so the seeds are the ones conformation mode
    # would actually have picked.
    refiner = ir.IterativeRefiner(active_targets=[args.target], active_loops=loops)
    seeds = refiner._beam_seeds(args.target, args.seeds)
    if not seeds:
        print("FAIL: beam returned no seeds.")
        return 1

    look = {}
    for f in glob.glob(str(out_base / "it_*" / "round_summary.csv")):
        d = pd.read_csv(f)
        if args.metric in d.columns:
            look.update(dict(zip(d.design_id, d[args.metric])))
    parents = [look.get(sid, float("nan")) for sid, _ in seeds]
    pmean = float(np.nanmean(parents))
    print(f"seeds: {len(seeds)} | parent mean {args.metric} = {pmean:.3f}")
    for sid, _ in seeds:
        print(f"    {sid}  {look.get(sid, float('nan')):.3f}")

    rows = []
    for t in tvals:
        kids, loop_rmsds, rmsd_err = [], [], 0
        for sid, seed in seeds:
            const_ids, _ = cseed._split_const_loop(seed, B, loops)
            sl = cseed.loop_residue_ids(seed, B, loops)
            lids = [int(i) for v in sl.values() for i in v]
            try:
                arrs = refiner.run_rfd3_partial(
                    args.target, seed, args.children, t,
                    seed_path=work / f"t{t:g}" / f"{sid}.pdb")
            except Exception as exc:
                print(f"  t={t:g} {sid}: partial diffusion failed: {exc}")
                continue
            for arr in arrs:
                # How far did the loops actually travel from this parent? Pin on the
                # constant scaffold first, so this is loop movement and not a rigid
                # shift. _ca_of returns an AtomArray, so take .coord -- subtracting the
                # arrays themselves is not elementwise and raises.
                try:
                    a = cseed._ca_of(seed, B, const_ids)
                    b = cseed._ca_of(arr, B, const_ids)
                    k = min(a.array_length(), b.array_length())
                    _, tr = superimpose(a[:k], b[:k])
                    la = cseed._ca_of(seed, B, lids).coord
                    lb = tr.apply(cseed._ca_of(arr, B, lids)).coord
                    k = min(len(la), len(lb))
                    loop_rmsds.append(float(np.sqrt(
                        ((la[:k] - lb[:k]) ** 2).sum(axis=1).mean())))
                except Exception as exc:
                    # Report, never swallow: a silently-NaN diagnostic column is how
                    # a broken measurement survives a whole campaign unnoticed.
                    if not rmsd_err:
                        print(f"    loop RMSD unavailable: {type(exc).__name__}: {exc}")
                    rmsd_err += 1
            kids += arrs

        if not kids:
            continue
        d = work / f"t{t:g}"
        cands = refiner.run_lmpnn(args.target, kids, d / "lmpnn",
                                  refiner.state.get("temperature", 0.15),
                                  origins=[("perturbed", "ladder")] * len(kids))
        scored = refiner.run_esmfold2(args.target, cands, d / "esmfold2")
        vals = [s.get(args.metric) for s in scored
                if s.get(args.metric) is not None and s.get(args.metric) == s.get(args.metric)]
        if not vals:
            keys = sorted(scored[0].keys()) if scored else []
            print(f"  t={t:g}: {len(scored)} records but none carry '{args.metric}'.")
            if keys:
                print(f"    available: {', '.join(k for k in keys if not k.startswith('_'))}")
            continue
        cmean = float(np.mean(vals))
        trans = (cmean - pop) / (pmean - pop) if pmean != pop else float("nan")
        rows.append({"partial_t": t, "n": len(vals), "loop_rmsd": float(np.mean(loop_rmsds))
                     if loop_rmsds else float("nan"), "child_mean": cmean,
                     "child_best": float(np.max(vals)), "transmission": trans})
        print(f"  t={t:g}: n={len(vals)} loopRMSD={rows[-1]['loop_rmsd']:.2f}A "
              f"child_mean={cmean:.3f} transmission={trans:+.2f}")

    if not rows:
        print("FAIL: no noise level produced scored children.")
        return 1
    df = pd.DataFrame(rows)
    df.to_csv(work / "ladder.csv", index=False)
    print(f"\npopulation {pop:.3f}  ->  parents {pmean:.3f}\n")
    print(df.round(3).to_string(index=False))
    good = df[df.transmission >= 0.5]
    if len(good):
        print(f"\nLargest t still transmitting >=0.5: {good.partial_t.max():g} A"
              f"  -> set PARTIAL_T_INIT to this.")
    else:
        print("\nNo tested t transmits >=0.5; extend the ladder downward "
              "(--t 0.5,1,1.5,2).")
    print(f"wrote {work / 'ladder.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
