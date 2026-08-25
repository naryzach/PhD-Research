#!/usr/bin/env python3
"""
analyze_af3_tranche.py — score AF3 result ZIPs against the local pool, WITHOUT
needing the stratified manifest.

Why manifest-free: the manifest is a sidecar that a re-run of the stratified export
can overwrite, and a careless sync can clobber. The AF3 ZIP, by contrast, carries
everything needed to identify each job:

    <job>/fold_<job>_job_request.json        -> the submitted sequences (ground truth)
    <job>/fold_<job>_summary_confidences_0.json -> ipTM / pTM / chain pTM / ranking

Chain A's sequence is the binder, which joins straight back to `full_seq` in the
round summaries for design_id, target, composite_score and the SV battery. Band
labels are recomputed as terciles of the local score within each target, which is
what the original banding meant anyway.

    python Generation/refinement/analyze_af3_tranche.py Local/iterative_refinement/folds_*.zip
    python Generation/refinement/analyze_af3_tranche.py *.zip --out af3_results.csv

Prints per-metric Spearman correlations (raw and target-centred) between the local
composite and each AF3 metric, plus any SV metrics present, and writes a tidy CSV.
"""
import argparse
import glob as _glob
import json
import os
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parents[1]


def jobs_from_zip(zpath: Path) -> pd.DataFrame:
    """One row per AF3 job: submitted binder/target sequence + confidence metrics."""
    z = zipfile.ZipFile(zpath)
    names = z.namelist()
    out = []
    for job in sorted({n.split("/")[0] for n in names if "/" in n}):
        req = f"{job}/fold_{job}_job_request.json"
        summ = f"{job}/fold_{job}_summary_confidences_0.json"   # model_0 = top ranked
        if req not in names or summ not in names:
            continue
        try:
            j = json.loads(z.read(req))
            j = j[0] if isinstance(j, list) else j
            chains = [s["proteinChain"]["sequence"] for s in j.get("sequences", [])
                      if "proteinChain" in s]
            s = json.loads(z.read(summ))
        except Exception as exc:
            print(f"  ! {job}: {exc}")
            continue
        if not chains:
            continue
        cp = s.get("chain_ptm") or [np.nan]
        out.append({
            "zip": zpath.name, "job": job,
            "binder_seq": chains[0],
            "target_len": len(chains[1]) if len(chains) > 1 else np.nan,
            "af3_iptm": s.get("iptm"), "af3_ptm": s.get("ptm"),
            "af3_aptm": cp[0] if len(cp) else np.nan,
            "af3_bptm": cp[1] if len(cp) > 1 else np.nan,
            "af3_ranking": s.get("ranking_score"),
        })
    return pd.DataFrame(out)


def load_pool(out_base: Path) -> pd.DataFrame:
    frames = []
    for csv in sorted(out_base.glob("it_*/round_summary.csv")):
        try:
            frames.append(pd.read_csv(csv))
        except Exception:
            pass
    if not frames:
        sys.exit(f"No round_summary.csv under {out_base}")
    pool = pd.concat(frames, ignore_index=True)
    return pool.drop_duplicates("full_seq")


def main():
    ap = argparse.ArgumentParser(description="Analyze AF3 tranche ZIPs (manifest-free).")
    ap.add_argument("zips", nargs="+", help="AF3 result ZIP(s); globs allowed.")
    ap.add_argument("--out-base", default=None,
                    help="Run directory holding it_*/ (default $REFINE_OUT_BASE or "
                         "Local/iterative_refinement).")
    ap.add_argument("--out", default=None, help="Write the joined table here (CSV).")
    args = ap.parse_args()

    out_base = Path(args.out_base) if args.out_base else Path(
        os.environ.get("REFINE_OUT_BASE") or (_ROOT / "Local" / "iterative_refinement"))
    out_base = out_base.resolve()

    paths = []
    for pat in args.zips:
        paths.extend(Path(p) for p in _glob.glob(pat))
    paths = [p for p in paths if p.suffix.lower() == ".zip"]
    if not paths:
        sys.exit("No ZIPs matched.")

    af = pd.concat([jobs_from_zip(p) for p in sorted(paths)], ignore_index=True)
    print(f"AF3 jobs parsed: {len(af)} from {len(paths)} zip(s)")

    pool = load_pool(out_base)
    m = af.merge(pool, left_on="binder_seq", right_on="full_seq", how="left")
    matched = m["design_id"].notna().sum()
    print(f"Joined to the pool by sequence: {matched}/{len(m)}"
          + ("" if matched == len(m) else "  (unmatched = designed in a run not in this pool)"))

    u = m.dropna(subset=["design_id", "composite_score"]).drop_duplicates("binder_seq").copy()
    n_dup = len(m.dropna(subset=["design_id"])) - len(u)
    print(f"Unique designs: {len(u)}" + (f"  ({n_dup} re-submitted across tranches)" if n_dup else ""))
    if u.empty:
        sys.exit("Nothing joined — is --out-base the run these designs came from?")

    # Terciles of the local score within each target: what "band" originally meant.
    u["band"] = (u.groupby("target_name")["composite_score"]
                  .transform(lambda s: pd.qcut(s.rank(method="first"), 3,
                                               labels=["LO", "MID", "HI"])))

    from scipy.stats import spearmanr
    af3_cols = [c for c in u.columns if c.startswith("af3_")]
    sv_cols = [c for c in u.columns
               if c.startswith("sv_") and pd.api.types.is_numeric_dtype(u[c])
               and u[c].notna().sum() >= 10]

    def corr(xcol, ycol):
        g = u.dropna(subset=[xcol, ycol])
        if len(g) < 8:
            return None
        r, p = spearmanr(g[xcol], g[ycol])
        x = g[xcol] - g.groupby("target_name")[xcol].transform("mean")
        y = g[ycol] - g.groupby("target_name")[ycol].transform("mean")
        rc, pc = spearmanr(x, y)
        return len(g), r, p, rc, pc

    print(f"\n=== local composite vs AF3 (n={len(u)}) ===")
    print(f"{'metric':14} {'n':>4} {'raw rho':>9} {'p':>7} {'centred rho':>12} {'p':>7}")
    for c in af3_cols:
        res = corr("composite_score", c)
        if res:
            n, r, p, rc, pc = res
            print(f"{c:14} {n:4} {r:+9.3f} {p:7.3f} {rc:+12.3f} {pc:7.3f}"
                  + ("  *" if min(p, pc) < 0.05 else ""))

    if sv_cols:
        print(f"\n=== SV structural metrics vs AF3 ipTM (the untested hypothesis) ===")
        rows = []
        for c in sv_cols:
            res = corr(c, "af3_iptm")
            if res:
                rows.append((c,) + res)
        rows.sort(key=lambda t: -abs(t[4]))
        print(f"{'sv metric':30} {'n':>4} {'centred rho':>12} {'p':>7}")
        for c, n, r, p, rc, pc in rows:
            print(f"{c:30} {n:4} {rc:+12.3f} {pc:7.3f}" + ("  *" if pc < 0.05 else ""))
    else:
        print("\n(no SV metrics on these designs — they predate --sv-battery; "
              "use --strat-min-iteration to draw a tranche from SV-covered rounds)")

    print(f"\n=== band means ===")
    print(u.groupby("band", observed=True).agg(
        n=("af3_iptm", "size"), local=("composite_score", "mean"),
        iptm=("af3_iptm", "mean"), aptm=("af3_aptm", "mean")).round(3).to_string())

    if args.out:
        keep = (["zip", "job", "design_id", "target_name", "band", "composite_score",
                 "binder_seq"] + af3_cols + sv_cols)
        u[[c for c in keep if c in u.columns]].to_csv(args.out, index=False)
        print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
