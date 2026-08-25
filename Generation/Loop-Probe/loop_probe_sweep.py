"""
loop_probe_sweep.py

Wrapper around `loop_probe.run_probe` that sweeps loop LENGTHS and iterates over
targets.

Each `--loops` token is a UNIT that designs ONLY its own loops — every other
loop keeps the native template sequence, fully fixed.  So the units are distinct
experiments, and the single-loop units are the baselines for the joint one:

  * `--loops AB C EF`   -> three independent MARGINAL sweeps: AB alone (C, EF
                           fixed), C alone, EF alone; each varies only its own
                           loop's length.
  * `--loops AB_C EF`   -> AB and C designed TOGETHER with their lengths swept
                           jointly (every AB x C combination) — captures their
                           interaction — plus EF alone.  Compare AB_C against the
                           standalone AB and C units to isolate the interaction.
  * `--grid`            -> one joint unit: design ALL selected loops together,
                           full cartesian product of their lengths.

(A single `loop_probe.py` run still co-designs whatever loops you pass it; the
per-unit isolation is what this sweep wrapper adds.)

Per-loop length ranges come from `--length-range 'AB=4-12,C=5-9'`, a config
file, or a sensible default (native-2 … native+4, clamped to LOOP_CONFIGS max).

For each sweep unit it also writes cross-sweep summaries via the GPU-free
`loop_probe_analysis` helpers:
  * `length_montage_<unit>.png`             — 20-AA heatmap at every sweep point
  * `length_trend_<unit>_<scheme>.png/.csv` — group frequency vs sweep point

Runs in the GPU `foundry` conda env (it drives RFd3 + LigandMPNN via loop_probe).

Examples
--------
    # marginal sweep of AB, C, EF across all six targets (defaults: 100 bb x 5 seq, T=0.5)
    python Generation/Loop-Probe/loop_probe_sweep.py

    # AB and C swept jointly, EF marginally, HADDOCK templates
    python Generation/Loop-Probe/loop_probe_sweep.py --loops AB_C EF --template-set haddock

    # everything from a config file
    python Generation/Loop-Probe/loop_probe_sweep.py --config Generation/Loop-Probe/sweep_config.example.yaml
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
sys.path.insert(0, str(_HERE.parent / "refinement"))  # pipeline moved 2026-08-25

from iterative_refinement import LOOP_CONFIGS  # noqa: E402
import loop_probe as lp  # noqa: E402
import loop_probe_analysis as lpa  # noqa: E402

logger = logging.getLogger("loop_probe_sweep")


# ── Length-range parsing ───────────────────────────────────────────────────────
def default_range(loop: str) -> tuple[int, int]:
    """Sensible per-loop default: native-2 … native+4, clamped to [1, max]."""
    lc = LOOP_CONFIGS[loop]
    return max(1, lc["normal"] - 2), min(lc["max"], lc["normal"] + 4)


def parse_ranges(spec, loops: list[str], full: bool) -> dict[str, list[int]]:
    """
    Build {loop: [lengths]}.

    `spec` may be a string ('AB=4-12,C=5-9') or a mapping from a config file
    ({'AB': [4, 12]} or {'AB': '4-12'} or {'AB': [4,5,6]} for an explicit list).
    `full` widens unspecified loops to native-2 … LOOP_CONFIGS max.
    """
    overrides: dict[str, list[int]] = {}
    if isinstance(spec, dict):
        for name, val in spec.items():
            if isinstance(val, str):
                lo, _, hi = val.partition("-")
                overrides[name] = list(range(int(lo), int(hi) + 1))
            elif isinstance(val, (list, tuple)):
                # [lo, hi] is a range; longer lists are explicit length sets
                vals = [int(v) for v in val]
                overrides[name] = (list(range(vals[0], vals[1] + 1))
                                   if len(vals) == 2 else vals)
            else:
                overrides[name] = [int(val)]
    elif spec:
        for tok in str(spec).replace(" ", "").split(","):
            if not tok:
                continue
            name, _, rng = tok.partition("=")
            lo, _, hi = rng.partition("-")
            overrides[name] = list(range(int(lo), int(hi or lo) + 1))

    out = {}
    for name in loops:
        if name in overrides:
            out[name] = overrides[name]
        elif full:
            out[name] = list(range(max(1, LOOP_CONFIGS[name]["normal"] - 2),
                                   LOOP_CONFIGS[name]["max"] + 1))
        else:
            lo, hi = default_range(name)
            out[name] = list(range(lo, hi + 1))
        if not out[name]:
            raise ValueError(f"empty length range for loop {name}")
    return out


# ── Cross-sweep summaries for one sweep unit ───────────────────────────────────
def _write_sweep_summaries(target: str, unit_label: str, loops_in_unit: list[str],
                           counts_by_key: dict[str, dict], out_dir: Path) -> None:
    """
    counts_by_key: {loop_name: {sweep_key: counts_df}} for the loops in this unit.
    Writes a montage + per-scheme trend figure for each loop of the unit.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    for loop in loops_in_unit:
        by_key = counts_by_key.get(loop, {})
        if not by_key:
            continue
        freq_by_key = {k: lpa.to_frequency(c) for k, c in by_key.items()}
        lpa.length_montage(
            freq_by_key, f"{target} | {loop} (unit {unit_label}): AA frequency vs length",
            out_dir / f"length_montage_{loop}.png")
        for scheme in lpa.PROPERTY_SCHEMES:
            trend = lpa.length_trend_matrix(by_key, scheme)
            trend.to_csv(out_dir / f"length_trend_{loop}_{scheme}.csv")
            lpa.heatmap(
                trend, f"{target} | {loop} (unit {unit_label}): {scheme} vs length",
                out_dir / f"length_trend_{loop}_{scheme}.png",
                cbar_label="mean group frequency", annotate=True, vmax=1.0)


def _load_counts(run_dir: Path, loop: str, n_usable: int):
    """Reload a run's per-position counts for cross-sweep aggregation."""
    p = run_dir / f"position_counts_{loop}.csv"
    if not p.exists():
        return None
    c = pd.read_csv(p, index_col=0).reindex(lpa.AA20).fillna(0).astype(int)
    c.attrs["n_sequences"] = n_usable
    return c


# ── Sweep drivers ──────────────────────────────────────────────────────────────
def run_unit(target: str, unit: list[str], ranges: dict[str, list[int]],
             base_out: Path, gen_kwargs: dict, resume: bool = True) -> list[dict]:
    """
    Run one sweep unit for one target.

    A unit designs ONLY its own loops (every other loop keeps the native template
    sequence, fully fixed):
      * a single loop 'AB'   -> redesign AB alone, sweep its length (marginal);
      * a group 'AB_C'       -> redesign AB and C together, joint cartesian sweep
                                of their lengths (captures loop-loop interaction).
    So 'AB', 'C', and 'AB_C' are distinct experiments — the single-loop units are
    the baselines you compare the joint unit against.

    With `resume=True`, any configuration whose summary.json already exists is
    skipped (its results are reloaded), so a restarted multi-day run continues
    where it left off.
    """
    unit_label = "_".join(unit)
    combos = list(itertools.product(*(ranges[l] for l in unit)))
    logger.info(f"[{target}] unit {unit_label} (designs {unit} only): "
                f"{len(combos)} length combination(s) x "
                f"({gen_kwargs['n_backbones']} backbones "
                f"x {gen_kwargs['seqs_per_backbone']} seqs)")

    summaries: list[dict] = []
    counts_by_key: dict[str, dict] = {l: {} for l in unit}
    for combo in combos:
        # Only this unit's loops are designed, at the swept lengths.
        lengths = dict(zip(unit, combo))
        # single-loop units key by the integer length so figures read "L06";
        # joint units key by the combination label, e.g. "AB6_C7".
        key = combo[0] if len(unit) == 1 else "_".join(
            f"{l}{v}" for l, v in zip(unit, combo))
        run_dir = (base_out / target / unit_label /
                   (f"L{combo[0]:02d}" if len(unit) == 1 else key))

        done = run_dir / "summary.json"
        if resume and done.exists():
            with open(done) as fh:
                summary = json.load(fh)
            logger.info(f"[{target}] {unit_label} {key}: already complete — skipping")
        else:
            summary = lp.run_probe(target=target, active_loops=unit,
                                   lengths=lengths, out_dir=run_dir, **gen_kwargs)
        summary["sweep_unit"] = unit_label
        summary["sweep_key"] = str(key)
        summaries.append(summary)

        for loop in unit:
            n_us = summary.get("per_loop", {}).get(loop, {}).get("n_usable", 0)
            c = _load_counts(run_dir, loop, n_us)
            if c is not None:
                counts_by_key[loop][key] = c

    _write_sweep_summaries(target, unit_label, unit, counts_by_key,
                           base_out / target / unit_label)
    return summaries


def order_jobs(targets: list[str], units: list[list[str]], order: str
               ) -> list[tuple[str, list[str]]]:
    """
    Flatten (target, unit) work into an execution order.
      * 'unit'   -> all targets of unit 1, then all of unit 2, …  (units run in
                    the order listed, so a cheap baseline unit listed first
                    finishes across every target before an expensive later unit).
      * 'target' -> finish one target completely before the next.
    """
    if order == "unit":
        return [(t, u) for u in units for t in targets]
    return [(t, u) for t in targets for u in units]


def _write_summary_csv(summaries: list[dict], base_out: Path) -> None:
    """(Re)write sweep_summary.csv from all summaries collected so far."""
    if not summaries:
        return
    rows = []
    for s in summaries:
        row = {"target": s["target"], "unit": s.get("sweep_unit"),
               "key": s.get("sweep_key"), "lengths": json.dumps(s["lengths"]),
               "temperature": s["temperature"], "n_designs": s["n_designs_total"]}
        for k, v in s.get("per_loop", {}).items():
            row[f"{k}_usable"] = v.get("n_usable")
            row[f"{k}_unique"] = v.get("n_unique")
        rows.append(row)
    pd.DataFrame(rows).to_csv(base_out / "sweep_summary.csv", index=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--targets", nargs="+", default=None, choices=lp.TARGET_NAMES,
                    help="targets to iterate (default: all six)")
    ap.add_argument("--length-range", default=None,
                    help="sweep ranges, e.g. 'AB=4-12,C=5-9'")
    ap.add_argument("--full-range", action="store_true",
                    help="sweep native-2..max for each loop (wide)")
    ap.add_argument("--grid", action="store_true",
                    help="one joint sweep over ALL selected loops (expensive)")
    ap.add_argument("--order", default=None, choices=["unit", "target"],
                    help="execution order: 'unit' (default) runs each unit across "
                         "all targets before the next unit — so a cheap unit "
                         "listed first finishes everywhere before an expensive "
                         "later unit; 'target' finishes one target at a time")
    ap.add_argument("--fresh", action="store_true",
                    help="recompute everything (default resumes, skipping configs "
                         "whose summary.json already exists)")
    lp.add_common_args(ap)
    args = ap.parse_args()

    cfg = lp.load_config(args.config)
    targets = lp.cfg_get(args.targets, cfg, "targets", lp.TARGET_NAMES)
    if isinstance(targets, str):
        targets = targets.split()
    loop_tokens = lp.cfg_get(args.loops, cfg, "loops", ["AB", "C", "EF"])
    if isinstance(loop_tokens, str):
        loop_tokens = loop_tokens.split()
    groups, design_loops = lp.parse_loop_tokens(loop_tokens)

    ranges = parse_ranges(
        lp.cfg_get(args.length_range, cfg, "length_range", None) or cfg.get("ranges"),
        design_loops,
        args.full_range or bool(cfg.get("full_range", False)))

    lp.setup_env()
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = lp.cfg_get(args.out_dir, cfg, "out_dir", None)
    base_out = Path(out_dir) if out_dir else lp.OUT_BASE / f"sweep_{stamp}"
    base_out.mkdir(parents=True, exist_ok=True)

    gen_kwargs = dict(
        n_backbones=lp.cfg_get(args.n_backbones, cfg, "n_backbones", lp.DEFAULT_N_BACKBONES),
        seqs_per_backbone=lp.cfg_get(args.seqs_per_backbone, cfg, "seqs_per_backbone",
                                     lp.DEFAULT_SEQS_PER_BACKBONE),
        temperature=lp.cfg_get(args.temperature, cfg, "temperature", lp.DEFAULT_TEMPERATURE),
        template_set=lp.cfg_get(args.template_set, cfg, "template_set", lp.DEFAULT_TEMPLATE_SET),
        template_dir=lp.cfg_get(args.template_dir, cfg, "template_dir", None),
        template_map=cfg.get("template_map"),
        binder_chain=lp.cfg_get(args.binder_chain, cfg, "binder_chain", None),
        target_chain=lp.cfg_get(args.target_chain, cfg, "target_chain", None),
        scaffold_len=lp.cfg_get(args.scaffold_len, cfg, "scaffold_len", None),
        seed=lp.cfg_get(args.seed, cfg, "seed", lp.DEFAULT_SEED),
        make_plots=not (args.no_plots or cfg.get("no_plots", False)),
    )
    grid = args.grid or bool(cfg.get("grid", False))
    order = lp.cfg_get(args.order, cfg, "order", "unit")
    resume = not (args.fresh or bool(cfg.get("fresh", False)))

    units = [list(design_loops)] if grid else groups
    jobs = order_jobs(targets, units, order)
    n_runs = sum(len(list(itertools.product(*(ranges[l] for l in u))))
                 for _t, u in jobs)
    manifest = {"created": stamp, "targets": targets, "loops": loop_tokens,
                "groups": groups, "design_loops": design_loops, "ranges": ranges,
                "mode": "grid" if grid else "units", "order": order,
                "resume": resume, "planned_runs": n_runs,
                "execution_order": [{"target": t, "unit": "_".join(u)} for t, u in jobs],
                "gen_kwargs": {k: (str(v) if isinstance(v, Path) else v)
                               for k, v in gen_kwargs.items()}}
    with open(base_out / "sweep_manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    logger.info(f"Sweep -> {base_out}")
    logger.info(f"  targets={targets} units={units} order={order} resume={resume}")
    logger.info(f"  {n_runs} runs x {gen_kwargs['n_backbones']}"
                f"x{gen_kwargs['seqs_per_backbone']} seqs, T={gen_kwargs['temperature']}")

    all_summaries: list[dict] = []
    for target, unit in jobs:  # ordered per `order`
        try:
            all_summaries.extend(
                run_unit(target, unit, ranges, base_out, gen_kwargs, resume=resume))
            _write_summary_csv(all_summaries, base_out)  # visible progress mid-run
        except Exception as exc:  # isolate per-(target,unit) failures on long runs
            logger.error(f"[{target}] unit {'_'.join(unit)} failed: {exc}", exc_info=True)

    logger.info(f"Sweep complete: {len(all_summaries)} runs -> {base_out}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    main()
