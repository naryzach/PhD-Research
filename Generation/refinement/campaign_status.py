#!/usr/bin/env python3
"""
campaign_status.py - is the conformation campaign healthy, and is the frontier moving?

Separates the two questions that got conflated for ten rounds at it_54-63:

  MECHANISM   are perturbed children beating fresh ones in the same round?
              That is the loop working. It shows up within 2-3 rounds and is the
              only thing worth restarting the run over.
  FRONTIER    is the per-target best sv_pdockq climbing?
              That is the loop compounding. Round-to-round scatter on per-round
              best is ~0.02, so this needs 10+ rounds before it means anything.

Run it daily. Intervene ONLY on a mechanism flag; a flat frontier this early is
not a result, and a flat frontier after the full window is the answer (ceiling).

    python Generation/refinement/campaign_status.py
    python Generation/refinement/campaign_status.py --since 69      # rounds after the it_68 fix
"""
import argparse
import glob
import json
import re
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd

TARGETS = ["MMP2", "MMP9", "ADAM10", "ADAM17"]
METRIC = "sv_pdockq"


def load(base):
    rows = []
    for f in glob.glob(str(base / "it_*" / "round_summary.csv")):
        try:
            d = pd.read_csv(f, on_bad_lines="skip")
        except Exception:
            continue
        if METRIC not in d.columns or "esm_plddt" not in d.columns:
            continue
        d = d[d.esm_plddt.notna()]
        if len(d):
            rows.append(d.assign(it=int(re.search(r"it_(\d+)", f).group(1))))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base",
                    default=str(Path(__file__).resolve().parents[2]
                                / "Local" / "iterative_refinement"))
    ap.add_argument("--since", type=int, default=69,
                    help="first round of the current configuration")
    ap.add_argument("--window", type=int, default=10,
                    help="rounds to run before the frontier verdict is meaningful")
    args = ap.parse_args()

    base = Path(args.base)
    a = load(base)
    if a.empty:
        print("no scored rounds found")
        return 1

    cur = a[a.it >= args.since]
    prior = a[a.it < args.since]
    n_rounds = cur.it.nunique()
    print(f"rounds since it_{args.since}: {n_rounds} "
          f"(frontier verdict at {args.window})\n")

    st = base / "refinement_state.json"
    if st.exists():
        try:
            s = json.loads(st.read_text())
            cs = s.get("conformation_start")
            if cs is not None and len(cur):
                k = int(cur.it.max()) - int(cs)
                pt = max(1.0, 3.0 * 0.90 ** k)
                print(f"iteration {s.get('iteration')} | partial_t ~ {pt:.2f} A "
                      f"(k={k} rounds since conformation_start={cs})\n")
        except Exception:
            pass

    print("MECHANISM - perturbed vs fresh, same rounds")
    print(f"{'target':8} {'fresh':>8} {'perturbed':>10} {'delta':>8} {'n_pert':>7}  flag")
    flags, measured = [], 0
    for t in TARGETS:
        g = cur[cur.target_name == t]
        if "bb_origin" not in g.columns or g.bb_origin.isna().all():
            print(f"{t:8}  (no provenance tag yet)")
            continue
        fr = g[g.bb_origin == "fresh"][METRIC].dropna()
        pe = g[g.bb_origin == "perturbed"][METRIC].dropna()
        if not len(fr) or not len(pe):
            print(f"{t:8}  (waiting for both groups)")
            continue
        measured += 1
        delta = pe.mean() - fr.mean()
        flag = ""
        if delta < 0.03:
            flag = "** transmission lost - investigate"
            flags.append(t)
        elif len(pe) < 50:
            flag = "(thin)"
        print(f"{t:8} {fr.mean():8.3f} {pe.mean():10.3f} {delta:+8.3f} "
              f"{len(pe):7} {flag}")

    print("\nFRONTIER - best sv_pdockq, before vs since")
    print(f"{'target':8} {'prior best':>11} {'best since':>11} {'gain':>8}")
    gains = {}
    for t in TARGETS:
        pb = prior[prior.target_name == t][METRIC].max() if len(prior) else np.nan
        cb = cur[cur.target_name == t][METRIC].max() if len(cur) else np.nan
        gains[t] = cb - pb
        fmt = lambda v, w: f"{v:{w}.3f}" if v == v else "-".rjust(w)
        g = f"{cb - pb:+8.3f}" if (cb == cb and pb == pb) else "-".rjust(8)
        print(f"{t:8} {fmt(pb, 11)} {fmt(cb, 11)} {g}")

    print()
    if not measured:
        print(f"WAITING: no completed round under this configuration yet "
              f"(nothing since it_{args.since}). A round takes several hours across "
              "four targets; the summary is written when it finishes. Nothing has "
              "been measured, so nothing is known - check back later.")
    elif flags:
        print(f"ACTION: transmission has collapsed on {', '.join(flags)}. "
              "The loop is not working there - worth a look.")
    elif n_rounds < args.window:
        print(f"HOLD: mechanism healthy on {measured}/{len(TARGETS)} targets, "
              f"{args.window - n_rounds} more rounds before the frontier number "
              "means anything. Change nothing.")
    elif not [v for v in gains.values() if v == v]:
        print("WAITING: rounds exist but none carry sv_pdockq yet.")
    elif max(v for v in gains.values() if v == v) >= 0.02:
        best = max((v, k) for k, v in gains.items() if v == v)
        print(f"CLIMBING: frontier gained {best[0]:+.3f} on {best[1]}. Keep running.")
    else:
        print("CEILING: mechanism worked and the frontier did not move over the full "
              "window. sv_pdockq ~0.6 is a real limit for this scaffold and these "
              "three loops - stop and order constructs spanning the range.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
