"""
loop_probe_audit.py

Audit a loop-probe sweep directory for completeness and silent failures, and
(optionally) clean out incomplete/suspect runs so a resume recomputes exactly
those.  GPU-free — run it anywhere after (or during) a sweep.

A run directory is judged from its own `summary.json` + CSVs against the sweep's
`sweep_manifest.json` (which records the planned run set and n_backbones ×
seqs_per_backbone).  Classification per planned run:

  OK        summary.json present, full design count, healthy per-loop usable
  SUSPECT   summary.json present but something is off (unreadable JSON, design
            shortfall, a loop with zero/low usable seqs, parse-fail spike, or a
            missing position_counts CSV) — likely a silent error
  MISSING   never ran (no directory / no summary.json)

Usage:
    python loop_probe_audit.py <sweep_dir>                 # report only
    python loop_probe_audit.py <sweep_dir> --clean         # + delete SUSPECT and
                                                            #   partial (no-summary) dirs
    python loop_probe_audit.py <sweep_dir> --min-usable-frac 0.9

After `--clean`, re-run the SAME sweep command with the SAME `--out-dir <sweep_dir>`
— resume skips the OK runs and recomputes the cleaned/MISSING ones.
"""

from __future__ import annotations

import json
import shutil
import argparse
import itertools
from pathlib import Path


def planned_runs(manifest: dict) -> list[tuple[str, str, str]]:
    """Reconstruct the full planned set of (target, unit_label, sub) from a manifest."""
    targets = manifest["targets"]
    ranges = manifest["ranges"]
    grid = manifest.get("mode") == "grid"
    units = [manifest["design_loops"]] if grid else manifest["groups"]
    runs = []
    for t in targets:
        for unit in units:
            combos = itertools.product(*(ranges[l] for l in unit))
            label = "_".join(unit)
            for c in combos:
                sub = f"L{c[0]:02d}" if len(unit) == 1 else \
                    "_".join(f"{l}{v}" for l, v in zip(unit, c))
                runs.append((t, label, sub))
    return runs


def audit_run(run_dir: Path, expect_designs: int,
              min_usable_frac: float) -> list[str]:
    """Return a list of issue strings for one run dir ([] == OK, ['MISSING'] == absent)."""
    summ = run_dir / "summary.json"
    if not summ.exists():
        return ["MISSING"] if not run_dir.exists() else ["NO_SUMMARY"]
    try:
        d = json.loads(summ.read_text())
    except Exception as e:
        return [f"BAD_JSON:{type(e).__name__}"]

    issues = []
    nd = d.get("n_designs_total", 0)
    if expect_designs and nd < expect_designs:
        issues.append(f"designs={nd}/{expect_designs}")
    for lp, v in d.get("per_loop", {}).items():
        nu = v.get("n_usable", 0)
        npf = v.get("n_parse_fail", 0)
        if nu == 0:
            issues.append(f"{lp}:ZERO_usable")
        elif expect_designs and nu < expect_designs * min_usable_frac:
            issues.append(f"{lp}:low_usable={nu}")
        if expect_designs and npf > expect_designs * 0.1:
            issues.append(f"{lp}:parsefail={npf}")
        if not (run_dir / f"position_counts_{lp}.csv").exists():
            issues.append(f"{lp}:no_counts_csv")
    return issues


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweep_dir", help="a Local/loop_probe/sweep_* directory")
    ap.add_argument("--clean", action="store_true",
                    help="delete SUSPECT and partial (no-summary) run dirs so a "
                         "resume recomputes them")
    ap.add_argument("--min-usable-frac", type=float, default=0.9,
                    help="flag a loop whose usable seqs fall below this fraction "
                         "of expected designs (default 0.9)")
    args = ap.parse_args()

    sweep = Path(args.sweep_dir)
    man_path = sweep / "sweep_manifest.json"
    if not man_path.exists():
        ap.error(f"no sweep_manifest.json in {sweep}")
    man = json.loads(man_path.read_text())
    gk = man.get("gen_kwargs", {})
    expect = (gk.get("n_backbones") or 0) * (gk.get("seqs_per_backbone") or 0)

    planned = planned_runs(man)
    ok, suspect, missing, to_clean = [], [], [], []
    for t, label, sub in planned:
        rd = sweep / t / label / sub
        iss = audit_run(rd, expect, args.min_usable_frac)
        rel = f"{t}/{label}/{sub}"
        if iss == ["MISSING"]:
            missing.append(rel)
        elif not iss:
            ok.append(rel)
        else:
            suspect.append((rel, iss))
            if rd.exists():
                to_clean.append(rd)

    print(f"Sweep: {sweep}")
    print(f"  expected designs/run: {expect}  (n_backbones x seqs_per_backbone)")
    print(f"  planned: {len(planned)}")
    print(f"  OK:      {len(ok)}")
    print(f"  SUSPECT: {len(suspect)}")
    print(f"  MISSING: {len(missing)}")
    for rel, iss in suspect[:40]:
        print(f"    SUSPECT {rel}: {', '.join(iss)}")
    if len(suspect) > 40:
        print(f"    ... +{len(suspect) - 40} more")
    if missing:
        from collections import Counter
        by_unit = Counter(r.split("/", 2)[1] for r in missing)
        print(f"  MISSING by unit: {dict(by_unit)}")

    remaining = len(missing) + len(suspect)
    if args.clean:
        for rd in to_clean:
            shutil.rmtree(rd, ignore_errors=True)
        print(f"\nCleaned {len(to_clean)} suspect/partial run dir(s).")
    print(f"\nTo finish: rerun the sweep with --out-dir {sweep} "
          f"(resume recomputes {remaining} run(s); {len(ok)} are kept).")
    if suspect and not args.clean:
        print("Re-run with --clean first to delete the SUSPECT dirs so they recompute.")


if __name__ == "__main__":
    main()
