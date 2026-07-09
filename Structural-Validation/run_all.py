"""
Orchestrator for the TIMP3 structural-validation pipeline.

Stages
  prep      (runs anywhere) build sequences, AF3 JSON, ESMFold2 inputs,
            crystal catalog, and HADDOCK dry-run inputs.
  compute   (YOUR side)     submit AF3 jobs, sbatch ESMFold2, run HADDOCK.
            This stage only prints instructions — it does not launch compute.
  analysis  (runs anywhere) score monomers vs crystal, score complexes,
            aggregate, visualise. Safe to run before compute is done — it
            reports what is still pending.

Usage:
  python Structural-Validation/run_all.py prep
  python Structural-Validation/run_all.py analysis
  python Structural-Validation/run_all.py all        # prep + analysis
  python Structural-Validation/run_all.py compute     # print run instructions
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable


def run(script: str, *args: str) -> None:
    path = HERE / script
    print(f"\n{'='*70}\n>> {script} {' '.join(args)}\n{'='*70}")
    subprocess.run([PY, str(path), *args], check=True)


def prep() -> None:
    run("sequences/build_sequences.py")
    run("folding/make_af3_json.py")
    run("folding/make_esmfold_inputs.py")
    run("crystal/catalog_structures.py")
    run("crystal/extract_targets.py")
    run("docking/run_haddock_all.py", "--dry-run")


def analysis() -> None:
    run("analysis/model_vs_crystal.py")
    run("analysis/complex_metrics.py")
    run("analysis/aggregate.py")
    run("analysis/visualize.py")


COMPUTE_MSG = """
COMPUTE STAGE — run these on your side, then re-run `analysis`:

1. AF3 (https://alphafoldserver.com), ~30 jobs/day cap:
   - upload af3/af3_monomers.json               (19 monomer jobs)
   - upload af3/af3_complexes_batch01..03.json  (optional co-folds, 1 file/day)
   - unzip each result set under af3/results/

2. ESMFold2 (GPU cluster, `esmfold2` conda env) — monomers AND complexes:
   sbatch esmfold2/run_esmfold2.slurm
   (writes esmfold2/results/<id>.pdb, esmfold2/results_complex/<c>__<t>.cif,
    esmfold2/esmfold2_metrics.csv, esmfold2/esmfold2_complex_metrics.csv)

3. HADDOCK3 (`haddock` conda env; needs the monomer models from steps 1-2 and
   the extracted crystal targets from prep):
   python Structural-Validation/docking/run_haddock_all.py --matrix full --jobs 8
   Tracks (construct×target inputs): AF3×AF3, ESMFold2×ESMFold2, AF3×Crystal,
   ESMFold2×Crystal  -> docking/<track>/best_models/<c>__<t>_HADDOCK.pdb
   (subset with --matrix compare|crystal, or --construct-source/--target-source)
"""


def main() -> None:
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    if stage == "prep":
        prep()
    elif stage == "analysis":
        analysis()
    elif stage == "compute":
        print(COMPUTE_MSG)
    elif stage == "all":
        prep()
        print(COMPUTE_MSG)
        analysis()
    else:
        print(__doc__)
        sys.exit(1)


if __name__ == "__main__":
    main()
