"""
loop_probe_sweep.py

Wrapper around `loop_probe.run_probe` that sweeps a RANGE of loop lengths for
each loop and iterates over every target.

Default (marginal) mode: for each selected loop, vary THAT loop's length across
a range while the other selected loops stay at their native length.  This
isolates each loop's length effect on LigandMPNN's per-position amino-acid
preferences, for every target.  (A full cartesian `--grid` mode is also
available but scales as the product of the ranges — use with care.)

For each (target, loop) it additionally writes two cross-length summaries via
the GPU-free `loop_probe_analysis` helpers:
  * length_montage_<loop>.png    — the raw 20-AA heatmap at every swept length,
                                    side by side
  * length_trend_<loop>_<scheme>.png / .csv — group frequency (position-averaged)
                                    vs length, per biochemical scheme

Runs in the GPU `foundry` conda env (it drives RFd3 + LigandMPNN through
loop_probe).

Examples
--------
    # native-centered marginal sweep of AB, C, EF across all six targets
    python Generation/Loop-Probe/loop_probe_sweep.py

    # sweep only AB over 4..12, on MMP2 and ADAM17, warmer sampling
    python Generation/Loop-Probe/loop_probe_sweep.py --targets MMP2 ADAM17 --loops AB \
        --length-range AB=4-12 --temperature 0.4
"""

from __future__ import annotations

import sys
import json
import logging
import argparse
import itertools
from pathlib import Path
from datetime import datetime

import pandas as pd

_HERE = Path(__file__).parent.resolve()
# own dir for sibling modules (loop_probe, loop_probe_analysis); parent
# Generation/ for iterative_refinement.
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent))

from iterative_refinement import LOOP_CONFIGS  # noqa: E402
import loop_probe as lp  # noqa: E402
import loop_probe_analysis as lpa  # noqa: E402

logger = logging.getLogger("loop_probe_sweep")


# ── Length-range parsing ────────────────────────────────────────────────────────
def default_range(loop: str) -> range:
    """Sensible per-loop default: native-2 .. native+4, clamped to [1, max]."""
    lc = LOOP_CONFIGS[loop]
    lo = max(1, lc["normal"] - 2)
    hi = min(lc["max"], lc["normal"] + 4)
    return range(lo, hi + 1)


def parse_ranges(spec: str | None, loops: list[str], full: bool) -> dict[str, list[int]]:
    """
    Build {loop: [lengths]}.  `spec` like 'AB=4-12,C=5-8' overrides per loop;
    `full` uses native-2..max; otherwise default_range().
    """
    overrides: dict[str, tuple[int, int]] = {}
    if spec:
        for tok in spec.replace(" ", "").split(","):
            if not tok:
                continue
            name, _, rng = tok.partition("=")
            lo, _, hi = rng.partition("-")
            overrides[name] = (int(lo), int(hi))
    out = {}
    for lp_name in loops:
        if lp_name in overrides:
            lo, hi = overrides[lp_name]
        elif full:
            lo, hi = max(1, LOOP_CONFIGS[lp_name]["normal"] - 2), LOOP_CONFIGS[lp_name]["max"]
        else:
            r = default_range(lp_name)
            lo, hi = r.start, r.stop - 1
        out[lp_name] = list(range(lo, hi + 1))
    return out


# ── Cross-length summaries for one (target, loop) ───────────────────────────────
def _write_length_summaries(target: str, loop: str,
                            counts_by_length: dict[int, pd.DataFrame],
                            out_dir: Path) -> None:
    if not counts_by_length:
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    freq_by_length = {L: lpa.to_frequency(c) for L, c in counts_by_length.items()}
    lpa.length_montage(
        freq_by_length, f"{target} | {loop}: AA frequency vs loop length",
        out_dir / f"length_montage_{loop}.png")
    for scheme in lpa.PROPERTY_SCHEMES:
        trend = lpa.length_trend_matrix(counts_by_length, scheme)
        trend.to_csv(out_dir / f"length_trend_{loop}_{scheme}.csv")
        lpa.heatmap(
            trend, f"{target} | {loop}: {scheme} vs length",
            out_dir / f"length_trend_{loop}_{scheme}.png",
            cbar_label="mean group frequency", annotate=True, vmax=1.0)


# ── Sweep drivers ───────────────────────────────────────────────────────────────
def marginal_sweep(target: str, loops: list[str], ranges: dict[str, list[int]],
                   base_out: Path, gen_kwargs: dict) -> list[dict]:
    """Vary one loop at a time; others held at native length."""
    native = lp.resolve_lengths(loops, {})
    summaries = []
    for swept in loops:
        counts_by_length: dict[int, pd.DataFrame] = {}
        for L in ranges[swept]:
            lengths = dict(native)
            lengths[swept] = L
            out_dir = base_out / target / swept / f"L{L:02d}"
            summary = lp.run_probe(target=target, active_loops=loops, lengths=lengths,
                                   out_dir=out_dir, **gen_kwargs)
            summaries.append(summary)
            # reload the swept loop's counts for the cross-length summary
            cpath = out_dir / f"position_counts_{swept}.csv"
            if cpath.exists():
                c = pd.read_csv(cpath, index_col=0).reindex(lpa.AA20).fillna(0).astype(int)
                c.attrs["n_sequences"] = summary["per_loop"][swept]["n_usable"]
                counts_by_length[L] = c
        _write_length_summaries(target, swept, counts_by_length, base_out / target / swept)
    return summaries


def grid_sweep(target: str, loops: list[str], ranges: dict[str, list[int]],
               base_out: Path, gen_kwargs: dict) -> list[dict]:
    """Full cartesian product of all loop lengths (can be large)."""
    combos = list(itertools.product(*(ranges[l] for l in loops)))
    logger.warning(f"[{target}] grid sweep = {len(combos)} length combinations "
                   f"x RFd3/LigandMPNN each")
    summaries = []
    for combo in combos:
        lengths = dict(zip(loops, combo))
        sig = "_".join(f"{l}{v}" for l, v in lengths.items())
        out_dir = base_out / target / "grid" / sig
        summaries.append(lp.run_probe(target=target, active_loops=loops,
                                      lengths=lengths, out_dir=out_dir, **gen_kwargs))
    return summaries


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--targets", nargs="+", default=list(lp.PROBE_TARGETS),
                    choices=list(lp.PROBE_TARGETS),
                    help="targets to iterate (default: all six)")
    ap.add_argument("--loops", nargs="+", default=["AB", "C", "EF"],
                    choices=list(LOOP_CONFIGS), help="loops to design (default AB C EF)")
    ap.add_argument("--length-range", default=None,
                    help="override sweep ranges, e.g. 'AB=4-12,C=5-8'")
    ap.add_argument("--full-range", action="store_true",
                    help="sweep native-2..max for each loop (wide)")
    ap.add_argument("--grid", action="store_true",
                    help="full cartesian product instead of marginal (expensive)")
    ap.add_argument("--n-backbones", type=int, default=40)
    ap.add_argument("--seqs-per-backbone", type=int, default=3)
    ap.add_argument("--temperature", type=float, default=0.3)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-dir", default=None,
                    help="sweep root (default Local/loop_probe/sweep_<timestamp>)")
    ap.add_argument("--no-plots", action="store_true",
                    help="per-run CSVs only; cross-length montages/trends still drawn")
    args = ap.parse_args()

    lp.setup_env()
    ranges = parse_ranges(args.length_range, args.loops, args.full_range)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base_out = Path(args.out_dir) if args.out_dir else lp.OUT_BASE / f"sweep_{stamp}"
    base_out.mkdir(parents=True, exist_ok=True)

    gen_kwargs = dict(
        n_backbones=args.n_backbones, seqs_per_backbone=args.seqs_per_backbone,
        temperature=args.temperature, seed=args.seed, make_plots=not args.no_plots,
    )
    manifest = {"created": stamp, "targets": args.targets, "loops": args.loops,
                "ranges": ranges, "mode": "grid" if args.grid else "marginal",
                "gen_kwargs": gen_kwargs}
    with open(base_out / "sweep_manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    logger.info(f"Sweep -> {base_out}  targets={args.targets} ranges={ranges}")

    all_summaries = []
    for target in args.targets:
        try:
            fn = grid_sweep if args.grid else marginal_sweep
            all_summaries.extend(fn(target, args.loops, ranges, base_out, gen_kwargs))
        except Exception as exc:  # isolate per-target failures for long unattended runs
            logger.error(f"[{target}] sweep failed: {exc}", exc_info=True)

    if all_summaries:
        rows = [{"target": s["target"], "lengths": json.dumps(s["lengths"]),
                 "temperature": s["temperature"], "n_designs": s["n_designs_total"],
                 **{f"{k}_usable": v["n_usable"] for k, v in s["per_loop"].items()}}
                for s in all_summaries]
        pd.DataFrame(rows).to_csv(base_out / "sweep_summary.csv", index=False)
    logger.info(f"Sweep complete: {len(all_summaries)} runs -> {base_out}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    main()
