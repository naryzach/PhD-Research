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
    ap.add_argument("--marginal-window", type=int, default=5,
                    help="rounds over which to measure whether the frontier is STILL "
                         "moving (the number that should drive the decision)")
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

    # CUMULATIVE gain answers "did this configuration ever help", which stays true
    # forever once it is true. MARGINAL gain answers "is it still helping", which is
    # the question that decides whether to keep burning GPU. Reporting only the
    # cumulative number made the tool say CLIMBING through eleven flat rounds at
    # it_74-85, where the running max moved +0.005/+0.001/0.000/0.000.
    print("\nFRONTIER - best sv_pdockq")
    print(f"{'target':8} {'prior best':>11} {'best since':>11} {'cumulative':>11} "
          f"{'last ' + str(args.marginal_window):>9}")
    gains, marginal = {}, {}
    cutoff = (cur.it.max() - args.marginal_window) if len(cur) else np.nan
    for t in TARGETS:
        g_prior = prior[prior.target_name == t][METRIC]
        g_cur = cur[cur.target_name == t]
        pb = g_prior.max() if len(g_prior) else np.nan
        cb = g_cur[METRIC].max() if len(g_cur) else np.nan
        # what the max was `marginal_window` rounds ago, vs what it is now
        earlier = g_cur[g_cur.it <= cutoff][METRIC].max() if len(g_cur) else np.nan
        gains[t] = cb - pb
        marginal[t] = (cb - earlier) if earlier == earlier else np.nan
        fmt = lambda v, w: f"{v:{w}.3f}" if v == v else "-".rjust(w)
        g = f"{cb - pb:+11.3f}" if (cb == cb and pb == pb) else "-".rjust(11)
        m = f"{marginal[t]:+9.3f}" if marginal[t] == marginal[t] else "-".rjust(9)
        print(f"{t:8} {fmt(pb, 11)} {fmt(cb, 11)} {g} {m}")

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
    elif [v for v in marginal.values() if v == v] and \
            max(v for v in marginal.values() if v == v) >= 0.02:
        best = max((v, k) for k, v in marginal.items() if v == v)
        print(f"CLIMBING: frontier gained {best[0]:+.3f} on {best[1]} in the last "
              f"{args.marginal_window} rounds. Keep running.")
    else:
        mx = max((v for v in marginal.values() if v == v), default=float("nan"))
        cum = max((v for v in gains.values() if v == v), default=float("nan"))
        print(f"PLATEAU: the last {args.marginal_window} rounds moved the frontier "
              f"by at most {mx:+.3f} (cumulative gain for this configuration is still "
              f"{cum:+.3f}, but that is history, not headroom).")
        print("  The configuration worked and has now converged. More rounds buy "
              "designs no in-silico instrument you have can tell apart - stop and "
              "order constructs spanning the range.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
