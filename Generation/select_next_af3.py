"""
select_next_af3.py

Build the next AF3 Server submission from a full-pool ESMFold2 scoring
(esmfold2_scores.csv produced by `score_with_esmfold2.py --all`).

Selection policy:
  - Rank each target's designs by the ESMFold2 composite (0.5*ipTM + 0.5*pLDDT/100,
    matching COMPOSITE_ESMFOLD2 in iterative_refinement.py).
  - Per-target top-K (equal quota, so no target monopolizes the daily 30 slots).
  - Exclude designs already sent to AF3 (default: the stratified_manifest.json
    sequences) so we don't spend slots re-testing.

Output:
  Local/iterative_refinement/af3_submission_esmfold2.json   <- upload to AF3
  Local/iterative_refinement/esmfold2_selection.csv          <- record of picks

When the results come back, `validate_boltz_filter.py <zip>` will automatically
join them to esmfold2_scores.csv by sequence, growing the ESMFold2-vs-AF3
calibration set (currently n=24) toward significance.

Usage (from repo root):
    python Generation/select_next_af3.py                 # 7 per target (28 total)
    python Generation/select_next_af3.py --per-target 6  # 24 total
"""

import json
import argparse
from pathlib import Path

import pandas as pd

OUT_BASE = (Path(__file__).parent / ".." / "Local" / "iterative_refinement").resolve()


def load_target_seqs() -> dict:
    """Map target_name -> target_seq from the round summaries (one per target)."""
    frames = []
    for it_dir in sorted(OUT_BASE.glob("it_*")):
        csv = it_dir / "round_summary.csv"
        if csv.exists():
            frames.append(pd.read_csv(csv))
    if not frames:
        raise SystemExit(f"No round_summary.csv under {OUT_BASE}")
    df = pd.concat(frames, ignore_index=True)
    return {t: g["target_seq"].dropna().iloc[0] for t, g in df.groupby("target_name")}


def load_excluded() -> set:
    """Sequences already sent to AF3 (default: the stratified manifest)."""
    excluded = set()
    man = OUT_BASE / "stratified_manifest.json"
    if man.exists():
        excluded |= {e["full_seq"] for e in json.loads(man.read_text())}
    return excluded


def main():
    ap = argparse.ArgumentParser(description="Select next AF3 batch by ESMFold2 composite.")
    ap.add_argument("--per-target", type=int, default=7, help="Designs per target (cap ~30/4).")
    ap.add_argument("--scores", default=str(OUT_BASE / "esmfold2_scores.csv"))
    ap.add_argument("--out", default=str(OUT_BASE / "af3_submission_esmfold2.json"))
    args = ap.parse_args()

    esm = pd.read_csv(args.scores)
    esm["esm_composite"] = 0.5 * esm["esm_iptm"] + 0.5 * (esm["esm_plddt"] / 100.0)

    target_seqs = load_target_seqs()
    excluded    = load_excluded()
    esm = esm[~esm["full_seq"].isin(excluded)]
    print(f"Pool after excluding {len(excluded)} already-tested designs: {len(esm)}")

    jobs, picks = [], []
    for tname in sorted(esm["target_name"].unique()):
        tseq = target_seqs.get(tname)
        if not tseq:
            print(f"  [{tname}] no target_seq found; skipping.")
            continue
        top = esm[esm.target_name == tname].nlargest(args.per_target, "esm_composite")
        for _, r in top.iterrows():
            i = len(jobs)
            name = f"esmsel_{tname}_{i:02d}"
            jobs.append({
                "name": name,
                "modelSeeds": [42],
                "sequences": [
                    {"proteinChain": {"sequence": r["full_seq"], "count": 1}},
                    {"proteinChain": {"sequence": tseq,          "count": 1}},
                ],
            })
            picks.append({"job_name": name, "design_id": r["design_id"], "target": tname,
                          "esm_iptm": r["esm_iptm"], "esm_plddt": r["esm_plddt"],
                          "esm_composite": r["esm_composite"]})
        print(f"  [{tname}] selected {len(top)} "
              f"(composite {top['esm_composite'].min():.3f}-{top['esm_composite'].max():.3f})")

    Path(args.out).write_text(json.dumps(jobs, indent=2))
    sel_csv = OUT_BASE / "esmfold2_selection.csv"
    pd.DataFrame(picks).to_csv(sel_csv, index=False)
    print(f"\nWrote {len(jobs)} jobs -> {Path(args.out).name}")
    print(f"Record -> {sel_csv.name}")
    print("Upload the JSON to AF3; afterwards run validate_boltz_filter.py on the results.")


if __name__ == "__main__":
    main()
