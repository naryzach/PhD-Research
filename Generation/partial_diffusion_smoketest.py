#!/usr/bin/env python3
"""
partial_diffusion_smoketest.py — verify the RFd3 partial-diffusion call in isolation.

Exercises exactly the path conformation mode uses, without running an iteration:
  1. pick the best scored design for a target from the pool
  2. graft its loops onto the canonical template (conformation_seed)
  3. call RFd3 with atom_array_input + partial_t + select_fixed_atoms
  4. measure how far the loops actually moved, and confirm the pinned scaffold
     and target did NOT move

Step 4 is the real test. A call that "succeeds" but perturbs the scaffold, or moves
the loops by ~0 A, is a silent failure that would waste days in a live run.

    python Generation/partial_diffusion_smoketest.py                 # MMP2, t=8
    python Generation/partial_diffusion_smoketest.py --target MMP9 --partial-t 5 -n 2

Run with the GPUs free — it competes with any live generation job.
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="MMP2")
    ap.add_argument("--partial-t", type=float, default=8.0)
    ap.add_argument("-n", "--n-designs", type=int, default=2)
    ap.add_argument("--out-base", default=None)
    args = ap.parse_args()

    import conformation_seed as cseed
    import iterative_refinement as ir

    out_base = Path(args.out_base) if args.out_base else ir.OUT_BASE
    loops = ["AB", "C", "EF"]

    # 1. best scored design for this target
    rows = []
    for f in glob.glob(str(out_base / "it_*" / "round_summary.csv")):
        it = int(re.search(r"it_(\d+)", f).group(1))
        d = pd.read_csv(f)
        if "esm_plddt" not in d.columns:
            continue
        d = d[(d.esm_plddt.notna()) & (d.target_name == args.target)]
        if len(d):
            rows.append(d.assign(it=it))
    if not rows:
        print(f"FAIL: no scored {args.target} designs under {out_base}")
        return 1
    best = pd.concat(rows, ignore_index=True).nlargest(1, "composite_score").iloc[0]
    cif = out_base / f"it_{best.it}" / "esmfold2" / args.target / "structures" / f"{best.design_id}.cif"
    tpl = str(ir.DATA_DIR / ir.TARGETS[args.target]["pdb"])
    print(f"seed design : {best.design_id}  composite={best.composite_score:.3f} "
          f"pdockq={best.get('sv_pdockq', float('nan')):.3f}")
    if not cif.exists():
        print(f"FAIL: predicted complex missing: {cif}")
        return 1

    # 2. graft
    scaf, tgt = cseed.drift_rmsd(str(cif), tpl, loops)
    seed = cseed.build_seed_template(str(cif), tpl, loops)
    print(f"raw drift   : scaffold {scaf:.2f} A | target {tgt:.2f} A")
    s2, t2 = 0.0, 0.0
    print(f"grafted seed: {seed.array_length()} atoms, scaffold/target reset to canonical")

    # 3. partial diffusion through the same method the pipeline uses
    refiner = ir.IterativeRefiner(active_targets=[args.target], active_loops=loops)
    print(f"\ncalling RFd3: partial_t={args.partial_t} A, n={args.n_designs} ...")
    try:
        outs = refiner.run_rfd3_partial(args.target, seed, args.n_designs, args.partial_t)
    except Exception as exc:
        import traceback
        print(f"FAIL: partial diffusion raised {type(exc).__name__}: {exc}")
        traceback.print_exc()
        return 1
    if not outs:
        print("FAIL: partial diffusion returned no backbones.")
        return 1
    print(f"returned {len(outs)} backbone(s)")

    # 4. did the right atoms move?
    from biotite.structure import superimpose, rmsd as bio_rmsd
    B, T = ir.DESIGN_BINDER_CHAIN, ir.DESIGN_TARGET_CHAIN
    sl = cseed.loop_residue_ids(seed, B, loops)
    loop_ids = np.concatenate([v for v in sl.values() if len(v)])
    const_ids, _ = cseed._split_const_loop(seed, B, loops)

    print(f"\n{'backbone':10} {'loop RMSD':>10} {'scaffold':>9} {'target':>8}")
    ok = True
    for i, arr in enumerate(outs):
        try:
            sc_s = cseed._ca_of(seed, B, const_ids)
            sc_o = cseed._ca_of(arr, B, const_ids)
            k = min(len(sc_s), len(sc_o))
            _, tr = superimpose(sc_s[:k], sc_o[:k])
            arr_f = tr.apply(arr)
            r_sc = float(bio_rmsd(sc_s[:k].coord, cseed._ca_of(arr_f, B, const_ids)[:k].coord))
            lo_s, lo_o = cseed._ca_of(seed, B, loop_ids), cseed._ca_of(arr_f, B, loop_ids)
            m = min(len(lo_s), len(lo_o))
            r_lo = float(bio_rmsd(lo_s[:m].coord, lo_o[:m].coord))
            tg_s, tg_o = cseed._ca_of(seed, T, np.unique(seed.res_id[seed.chain_id == T])), None
            to = arr_f[(arr_f.chain_id == T) & (arr_f.atom_name == "CA")]
            j = min(len(tg_s), len(to))
            r_tg = float(bio_rmsd(tg_s[:j].coord, to[:j].coord))
        except Exception as exc:
            print(f"  backbone {i}: could not measure ({exc})")
            ok = False
            continue
        print(f"{i:10} {r_lo:10.2f} {r_sc:9.2f} {r_tg:8.2f}")
        if r_lo < 0.2:
            print("   ^ loops barely moved — partial_t may not be taking effect")
            ok = False
        if r_sc > 1.0 or r_tg > 1.0:
            print("   ^ pinned region moved — select_fixed_atoms may not be honoured")
            ok = False

    print("\nPASS: loops perturbed, scaffold and target held." if ok else
          "\nCHECK THE NUMBERS ABOVE before running a live campaign.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
