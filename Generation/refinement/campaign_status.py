#!/usr/bin/env python3
"""
campaign_status.py - is a refinement campaign healthy, and is it still making progress?

Works on BOTH pools and auto-detects which one it is pointed at:
  generation   Local/iterative_refinement    metric: sv_pdockq
  specificity  Local/specificity_refinement  metric: composite_score (on-quality +
                                             selectivity, which is what the beam breeds on)

Four questions, deliberately separated, because conflating them has cost this project
weeks:

  YIELD       what fraction of designs actually got scored?
              A near-constant count each round is the signature of the whole-batch
              ESMFold2 timeout, not a model failure. The specificity run sat at
              35/150 and 60/150 for a week and every log line looked routine.
  MECHANISM   are conformation-seeded designs beating freshly generated ones in the
              SAME round? That is the loop working. Visible in 2-3 rounds, and the
              only thing worth restarting a run over.
  FRONTIER    is the best score still climbing? Reported as MARGINAL gain over the
              last N rounds, not cumulative -- cumulative stays true forever once
              true, and kept reporting CLIMBING through eleven flat rounds.
  SELECTIVITY (specificity pools only) is the on/off gap moving off zero?

    python Generation/refinement/campaign_status.py
    python Generation/refinement/campaign_status.py --base Local/specificity_refinement
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

_ROOT = Path(__file__).resolve().parents[2]


def load(base):
    """Every row, scored or not -- yield needs the unscored ones."""
    rows = []
    for f in glob.glob(str(base / "it_*" / "round_summary.csv")):
        try:
            d = pd.read_csv(f, on_bad_lines="skip")
        except Exception:
            continue
        if len(d):
            rows.append(d.assign(it=int(re.search(r"it_(\d+)", f).group(1))))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def fmt(v, w, prec=3):
    return f"{v:{w}.{prec}f}" if v == v else "-".rjust(w)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=str(_ROOT / "Local" / "iterative_refinement"))
    ap.add_argument("--since", type=int, default=None,
                    help="first round of the current configuration (default: all)")
    ap.add_argument("--metric", default=None, help="override the ranking metric")
    ap.add_argument("--marginal-window", type=int, default=5)
    ap.add_argument("--window", type=int, default=10,
                    help="rounds before the frontier verdict is meaningful")
    args = ap.parse_args()

    base = Path(args.base)
    a = load(base)
    if a.empty:
        print(f"no round summaries under {base}")
        return 1

    spec = "on_target" in a.columns or "selectivity_cd" in a.columns
    metric = args.metric or ("composite_score" if spec else "sv_pdockq")
    if metric not in a.columns:
        print(f"metric '{metric}' not in these summaries; columns include "
              f"{[c for c in a.columns if c.startswith(('sv_','composite','selectivity'))][:6]}")
        return 1
    targets = [t for t in a.target_name.dropna().unique()]
    since = args.since if args.since is not None else int(a.it.min())
    cur, prior = a[a.it >= since], a[a.it < since]
    n_rounds = cur.it.nunique()

    print(f"{'SPECIFICITY' if spec else 'GENERATION'} pool: {base}")
    print(f"rounds since it_{since}: {n_rounds}  |  ranking metric: {metric}\n")

    st = next((p for p in (base / "specificity_state.json",
                           base / "refinement_state.json") if p.exists()), None)
    if st:
        try:
            s = json.loads(st.read_text())
            cs_ = s.get("conformation_start")
            line = f"iteration {s.get('iteration')} | T={s.get('temperature', float('nan')):.3f}"
            if cs_ is not None and len(cur):
                k = int(cur.it.max()) - int(cs_)
                line += f" | partial_t ~ {max(1.0, 3.0 * 0.90 ** k):.2f} A (k={k})"
            print(line + "\n")
        except Exception:
            pass

    # ── YIELD ────────────────────────────────────────────────────────────────
    scored_col = "on_target" if spec else "esm_plddt"
    print("YIELD - designs actually scored")
    print(f"{'target':8} {'rounds':>7} {'mean':>7} {'min':>6} {'spread':>7}  flag")
    yield_flags = []
    for t in targets:
        g = cur[cur.target_name == t]
        per = g.groupby("it").apply(
            lambda d: 100.0 * d[scored_col].notna().sum() / max(1, len(d)))
        if per.empty:
            continue
        spread = per.max() - per.min()
        flag = ""
        # A LOW and near-CONSTANT count is the timeout signature: a model fault varies,
        # a wall-clock cut-off does not.
        if per.mean() < 70 and spread < 12:
            flag = "** batch timeout signature (low + near-constant)"
            yield_flags.append(t)
        elif per.mean() < 70:
            flag = "** low yield, varying - check for OOM"
            yield_flags.append(t)
        print(f"{t:8} {len(per):7} {per.mean():6.0f}% {per.min():5.0f}% "
              f"{spread:6.0f}%  {flag}")

    # ── MECHANISM ────────────────────────────────────────────────────────────
    print("\nMECHANISM - perturbed vs fresh, same rounds")
    print(f"{'target':8} {'fresh':>8} {'perturbed':>10} {'delta':>8} {'n_pert':>7}  flag")
    mech_flags, measured = [], 0
    for t in targets:
        g = cur[cur.target_name == t]
        if "bb_origin" not in g.columns or g.bb_origin.isna().all():
            print(f"{t:8}  (no provenance tag yet)")
            continue
        fr = g[g.bb_origin == "fresh"][metric].dropna()
        pe = g[g.bb_origin == "perturbed"][metric].dropna()
        if not len(fr) or not len(pe):
            print(f"{t:8}  (waiting for both groups)")
            continue
        measured += 1
        delta = pe.mean() - fr.mean()
        flag = "** transmission lost" if delta < 0.03 else ("(thin)" if len(pe) < 50 else "")
        if delta < 0.03:
            mech_flags.append(t)
        print(f"{t:8} {fr.mean():8.3f} {pe.mean():10.3f} {delta:+8.3f} {len(pe):7}  {flag}")

    # ── FRONTIER ─────────────────────────────────────────────────────────────
    print(f"\nFRONTIER - best {metric}")
    print(f"{'target':8} {'prior':>9} {'best now':>10} {'cumulative':>11} "
          f"{'last ' + str(args.marginal_window):>9}")
    gains, marginal = {}, {}
    cutoff = (cur.it.max() - args.marginal_window) if len(cur) else np.nan
    for t in targets:
        gp = prior[prior.target_name == t][metric]
        gc = cur[cur.target_name == t]
        pb = gp.max() if len(gp) else np.nan
        cb = gc[metric].max() if len(gc) else np.nan
        earlier = gc[gc.it <= cutoff][metric].max() if len(gc) else np.nan
        gains[t] = cb - pb
        marginal[t] = (cb - earlier) if earlier == earlier else np.nan
        g_ = f"{gains[t]:+11.3f}" if (cb == cb and pb == pb) else "-".rjust(11)
        m_ = f"{marginal[t]:+9.3f}" if marginal[t] == marginal[t] else "-".rjust(9)
        print(f"{t:8} {fmt(pb, 9)} {fmt(cb, 10)} {g_} {m_}")

    # ── SELECTIVITY (specificity only) ───────────────────────────────────────
    if spec:
        print("\nSELECTIVITY - on-target minus off-target (positive = wanted direction)")
        print(f"{'target':8} {'metric':20} {'n':>6} {'mean':>8} {'best':>8} {'frac>0':>7}")
        for t in targets:
            g = cur[cur.target_name == t]
            for col in ("selectivity_pdockq", "selectivity_cd"):
                if col not in g.columns:
                    continue
                s = g[col].dropna()
                if not len(s):
                    continue
                print(f"{t:8} {col:20} {len(s):6} {s.mean():+8.3f} {s.max():+8.3f} "
                      f"{(s > 0).mean():6.0%}")

    # ── VERDICT ──────────────────────────────────────────────────────────────
    print()
    if yield_flags:
        print(f"ACTION: yield is low on {', '.join(yield_flags)}. Fix this FIRST - every "
              "other number here is computed on a truncated, order-biased sample.")
    elif not measured:
        print(f"WAITING: no completed round with provenance under this configuration "
              f"yet. Nothing measured, so nothing known.")
    elif mech_flags:
        print(f"ACTION: transmission has collapsed on {', '.join(mech_flags)}. "
              "The conformation loop is not working there.")
    elif n_rounds < args.window:
        print(f"HOLD: mechanism healthy on {measured}/{len(targets)} targets, "
              f"{args.window - n_rounds} more rounds before the frontier number means "
              "anything. Change nothing.")
    elif [v for v in marginal.values() if v == v] and \
            max(v for v in marginal.values() if v == v) >= 0.02:
        best = max((v, k) for k, v in marginal.items() if v == v)
        print(f"CLIMBING: frontier gained {best[0]:+.3f} on {best[1]} in the last "
              f"{args.marginal_window} rounds. Keep running.")
    else:
        mx = max((v for v in marginal.values() if v == v), default=float("nan"))
        print(f"PLATEAU: the last {args.marginal_window} rounds moved the frontier by at "
              f"most {mx:+.3f}. The configuration worked and has converged - stop and "
              "take what it produced.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
