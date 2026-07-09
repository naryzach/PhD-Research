"""
Generate AlphaFold3-server batch JSON from the sequence registry.

Emits two kinds of jobs (AlphaFold Server schema, same as AlphaFold/AF_batch_gen.py):
  1. MONOMERS   - every construct and target folded alone (19 jobs -> one file).
                  These are the models compared to crystal structures and fed
                  to HADDOCK.
  2. COMPLEXES  - every construct x target co-fold (14 x 5 = 70 jobs). Optional
                  but recommended: gives an independent complex model plus AF3
                  ipTM/PAE. Chunked into files of <= MAX_JOBS_PER_DAY because the
                  public server caps submissions per day.

Run:  python Structural-Validation/folding/make_af3_json.py
Then: upload each JSON at https://alphafoldserver.com (one file per day for
      the complex batches), download results into af3/results/.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C  # noqa: E402

MAX_JOBS_PER_DAY = 30


def load_registry() -> list[dict]:
    with open(C.OUT_SEQ / "registry.csv", newline="") as fh:
        return list(csv.DictReader(fh))


def _job(name: str, chains: list[str]) -> dict:
    return {
        "name": name,
        "sequences": [
            {"proteinChain": {"sequence": s, "count": 1}} for s in chains
        ],
        "modelSeeds": [],
    }


def build_monomers(rows: list[dict]) -> list[dict]:
    return [_job(r["id"], [r["sequence"]]) for r in rows]


def build_complexes(rows: list[dict]) -> list[dict]:
    constructs = [r for r in rows if r["kind"] == "construct"]
    targets = [r for r in rows if r["kind"] == "target"]
    jobs = []
    for c in constructs:
        for t in targets:
            jobs.append(_job(f"{c['id']}__{t['id']}", [c["sequence"], t["sequence"]]))
    return jobs


def write_chunked(jobs: list[dict], stem: str) -> list[Path]:
    C.OUT_AF3.mkdir(parents=True, exist_ok=True)
    paths = []
    if len(jobs) <= MAX_JOBS_PER_DAY:
        p = C.OUT_AF3 / f"{stem}.json"
        p.write_text(json.dumps(jobs, indent=2))
        paths.append(p)
    else:
        for i in range(0, len(jobs), MAX_JOBS_PER_DAY):
            chunk = jobs[i:i + MAX_JOBS_PER_DAY]
            n = i // MAX_JOBS_PER_DAY + 1
            p = C.OUT_AF3 / f"{stem}_batch{n:02d}.json"
            p.write_text(json.dumps(chunk, indent=2))
            paths.append(p)
    return paths


def main() -> None:
    rows = load_registry()
    (C.OUT_AF3 / "results").mkdir(parents=True, exist_ok=True)

    mono = build_monomers(rows)
    cplx = build_complexes(rows)

    mono_paths = write_chunked(mono, "af3_monomers")
    cplx_paths = write_chunked(cplx, "af3_complexes")

    print(f"Monomer jobs : {len(mono)}")
    for p in mono_paths:
        print(f"  {p.relative_to(C.REPO_ROOT)}")
    print(f"\nComplex jobs : {len(cplx)}  "
          f"({len(cplx_paths)} daily batches of <= {MAX_JOBS_PER_DAY})")
    for p in cplx_paths:
        print(f"  {p.relative_to(C.REPO_ROOT)}")
    print("\nSubmit at https://alphafoldserver.com ; save downloads under "
          f"{(C.OUT_AF3 / 'results').relative_to(C.REPO_ROOT)}/")


if __name__ == "__main__":
    main()
