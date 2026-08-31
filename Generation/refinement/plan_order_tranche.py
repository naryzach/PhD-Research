#!/usr/bin/env python3
"""
plan_order_tranche.py - plan a construct order that spans the sv_pdockq range.

Why a spread and not a top-N. Two independent lines of evidence say ranking by
in-silico score does not by itself justify what you synthesise:

  * Local/TIMP3_Redesign_2026-07/docs/09_exact_calibration_final.md - on exact-matched,
    expression-controlled sequences, NO in-silico metric reliably predicted binding
    (best AUC ~0.68). Its recommendation 3 is to deliberately span predicted-strong to
    predicted-weak so the next calibration has variance to fit.
  * The 2026-08-26 AF3 tranche - above sv_pdockq ~0.5, AF3 ipTM saturates (upper-half
    IQR 0.023, 11/12 at >=0.80), so AF3 can no longer rank the top designs at all.

A top-N order answers "are our best designs good?". A spread order answers the question
the calibration loop was built for: does sv_pdockq track real binding at the top of its
range, or is everything above ~0.5 equivalent at the bench?

Emits two lists:
  ORDER   - the proposed constructs, banded, preferring designs already AF3-validated
  AF3     - designs that need folding first, because ordering gates on AF3 (and the
            top band is currently under-sampled, MMP9 worst)

    python Generation/refinement/plan_order_tranche.py --af3-joined <joined.csv>
    python Generation/refinement/plan_order_tranche.py --per-band 2 --min-iteration 69
"""
import argparse
import glob
import re
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parents[1]

TARGETS = ["MMP2", "MMP9", "ADAM10", "ADAM17"]
# Bands are fixed on the sv_pdockq scale, not quantiles: the point is to span the
# PREDICTION axis evenly so the FCS regression has leverage across it. Quantile bands
# would follow the pool's own distribution and pile up wherever the search converged.
BANDS = [("weak", 0.15, 0.30), ("mid", 0.30, 0.45),
         ("strong", 0.45, 0.58), ("frontier", 0.58, 1.01)]
PTM_MIN, IPTM_MIN = 0.70, 0.50


def load_pool(out_base, min_iteration):
    rows = []
    for f in glob.glob(str(out_base / "it_*" / "round_summary.csv")):
        it = int(re.search(r"it_(\d+)", f).group(1))
        try:
            d = pd.read_csv(f, on_bad_lines="skip")
        except Exception:
            continue
        if "sv_pdockq" not in d.columns or "esm_plddt" not in d.columns:
            continue
        d = d[d.esm_plddt.notna() & d.sv_pdockq.notna()]
        if len(d):
            rows.append(d.assign(it=it))
    if not rows:
        return pd.DataFrame()
    a = pd.concat(rows, ignore_index=True)
    a = a[a.it >= min_iteration] if min_iteration else a
    # Cysteines in the DESIGNED loops. TIMP3 already carries 6 structural disulfides,
    # so a free loop Cys risks scrambling on yeast display. It is also a confound:
    # measured on the it_0-80 pool, mean loop-Cys rises monotonically with sv_pdockq
    # (0.68 / 0.78 / 0.95 / 1.20 / 1.39 across the bands below), and the search has
    # converged on a motif -- EF loops ending TICD / SRCD / RLCD / SLCD / SICD recur
    # across independent designs and different targets. Spreading the order on
    # sv_pdockq WITHOUT holding Cys count down would vary the two together, so a
    # binding correlation could not be attributed to either.
    loops = [c for c in ("loop_AB_seq", "loop_C_seq", "loop_EF_seq") if c in a.columns]
    a["loop_cys"] = a[loops].fillna("").agg("".join, axis=1).str.upper().str.count("C")
    # One row per SEQUENCE: the same binder can recur across rounds, and ordering the
    # same construct twice wastes a synthesis slot.
    return a.sort_values("sv_pdockq", ascending=False).drop_duplicates("full_seq")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-base", default=str(_ROOT / "Local" / "iterative_refinement"))
    ap.add_argument("--af3-joined", default=None,
                    help="CSV from analyze_af3_tranche.py --out (all zips)")
    ap.add_argument("--per-band", type=int, default=1,
                    help="constructs per target per band")
    ap.add_argument("--min-iteration", type=int, default=0)
    ap.add_argument("--af3-budget", type=int, default=30, help="AF3 jobs/day")
    ap.add_argument("--max-loop-cys", type=int, default=1,
                    help="max cysteines across the designed loops (see load_pool); "
                         "-1 disables the filter")
    ap.add_argument("--emit-af3", default=None, help="write the AF3 list to this CSV")
    args = ap.parse_args()

    pool = load_pool(Path(args.out_base), args.min_iteration)
    if pool.empty:
        print("no scored designs found")
        return 1
    print(f"pool: {len(pool)} unique sequences, sv_pdockq "
          f"{pool.sv_pdockq.min():.3f}-{pool.sv_pdockq.max():.3f}")
    pool["band_"] = pd.cut(pool.sv_pdockq, [b[1] for b in BANDS] + [BANDS[-1][2]],
                           labels=[b[0] for b in BANDS])
    cys_by_band = pool.groupby("band_").loop_cys.mean().round(2)
    print(f"  mean loop-Cys by band: "
          + "  ".join(f"{k}={v}" for k, v in cys_by_band.items()))
    if args.max_loop_cys >= 0:
        before = len(pool)
        pool = pool[pool.loop_cys <= args.max_loop_cys]
        after = pool.groupby("band_").loop_cys.mean().round(2)
        print(f"  loop_cys <= {args.max_loop_cys}: kept {len(pool)} of {before}")
        print(f"  after filter          : "
              + "  ".join(f"{k}={v}" for k, v in after.items()))
    print()

    af3 = pd.DataFrame()
    if args.af3_joined and Path(args.af3_joined).exists():
        af3 = pd.read_csv(args.af3_joined)
        keep = [c for c in ("design_id", "af3_iptm", "af3_ptm") if c in af3.columns]
        af3 = af3[keep].drop_duplicates("design_id")
        pool = pool.merge(af3, on="design_id", how="left")
    for c in ("af3_iptm", "af3_ptm"):
        if c not in pool.columns:
            pool[c] = np.nan
    pool["af3_ok"] = (pool.af3_ptm >= PTM_MIN) & (pool.af3_iptm >= IPTM_MIN)
    pool["has_af3"] = pool.af3_iptm.notna()

    order, need = [], []
    for t in TARGETS:
        for band, lo, hi in BANDS:
            g = pool[(pool.target_name == t) & (pool.sv_pdockq >= lo) & (pool.sv_pdockq < hi)]
            if g.empty:
                print(f"  !! {t:7} {band:9} EMPTY - no design in {lo:.2f}-{hi:.2f}")
                continue
            # Prefer an AF3-validated design; fall back to the best-folded unvalidated
            # one, which then has to be folded before it can be ordered.
            ok = g[g.af3_ok == True]                                   # noqa: E712
            pick = (ok.nlargest(args.per_band, "esm_plddt") if len(ok)
                    else g.nlargest(args.per_band, "esm_plddt"))
            for _, r in pick.iterrows():
                rec = {"target": t, "band": band, "design_id": r.design_id,
                       "sv_pdockq": r.sv_pdockq, "esm_plddt": r.esm_plddt,
                       "it": int(r.it), "af3_iptm": r.af3_iptm,
                       "loop_cys": int(r.loop_cys),
                       "validated": bool(r.af3_ok is True or r.af3_ok == True),
                       "full_seq": r.full_seq}
                order.append(rec)
                if not rec["validated"]:
                    need.append(rec)

    o = pd.DataFrame(order)
    print("=== PROPOSED ORDER ===")
    print(f"{'target':8} {'band':9} {'design_id':26} {'pdockq':>7} {'plddt':>6} "
          f"{'it':>4} {'cys':>4} {'af3_iptm':>9}  status")
    for _, r in o.iterrows():
        a = "-" if pd.isna(r.af3_iptm) else f"{r.af3_iptm:.3f}"
        print(f"{r.target:8} {r.band:9} {r.design_id:26} {r.sv_pdockq:7.3f} "
              f"{r.esm_plddt:6.1f} {r.it:4} {r.loop_cys:4} {a:>9}  "
              f"{'AF3-validated' if r.validated else 'NEEDS AF3'}")

    print(f"\n{len(o)} constructs | {int(o.validated.sum())} already AF3-validated | "
          f"{len(need)} need folding")

    if need:
        n = pd.DataFrame(need)
        print(f"\n=== RUN THESE IN AF3 ({len(n)} jobs, budget {args.af3_budget}/day) ===")
        for _, r in n.iterrows():
            print(f"\n>{r.target}_{r.band}_{r.design_id}   "
                  f"[sv_pdockq {r.sv_pdockq:.3f}, it_{r.it}]")
            print(r.full_seq)
        if args.emit_af3:
            n.to_csv(args.emit_af3, index=False)
            print(f"\nwrote {args.emit_af3}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
