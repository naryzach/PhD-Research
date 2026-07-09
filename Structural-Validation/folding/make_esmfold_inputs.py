"""
Build ESMFold2 inputs from the SAME registry the AF3 JSON is built from, so the
two folders model an identical set of sequences.

Emits (into Local/TIMP3_Structural_Validation_2026-07/esmfold2/):
  * esmfold2_input.fasta   every construct + target, one record each (19)
  * run_esmfold2.py        the folding runner (copied from folding/run_esmfold2.py)
  * run_esmfold2.slurm     GPU batch script matching the lab's SLURM conventions

The runner (folding/run_esmfold2.py, biohub/ESMFold2 API) writes:
  * results/<id>.pdb|cif + esmfold2_metrics.csv      (monomers)
  * results_complex/<c>__<t>.cif|pdb + esmfold2_complex_metrics.csv
    with esm_iptm/esm_ptm/esm_plddt/esm_lplddt/esm_pae (construct×target co-folds)

Run:  python Structural-Validation/folding/make_esmfold_inputs.py
"""
from __future__ import annotations

import csv
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C  # noqa: E402

HERE = Path(__file__).resolve().parent

SLURM_TEMPLATE = """#!/bin/bash
#SBATCH --job-name=TIMP3_valid_ESMFold2
#SBATCH --time=2-00:00:00
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=esmfold2_%j.log
#SBATCH --error=esmfold2_%j.err

# ESMFold2 runs in its own env (Python >= 3.12); see Generation/score_with_esmfold2.py
source ~/miniconda3/etc/profile.d/conda.sh
conda activate esmfold2

cd "{outdir}"
python run_esmfold2.py \\
    --mode both \\
    --fasta esmfold2_input.fasta \\
    --pairs esmfold2_complex_pairs.csv \\
    --outdir results --complex-outdir results_complex \\
    --metrics esmfold2_metrics.csv \\
    --complex-metrics esmfold2_complex_metrics.csv

conda deactivate
"""


def load_registry() -> list[dict]:
    with open(C.OUT_SEQ / "registry.csv", newline="") as fh:
        return list(csv.DictReader(fh))


def main() -> None:
    C.OUT_ESM.mkdir(parents=True, exist_ok=True)
    (C.OUT_ESM / "results").mkdir(parents=True, exist_ok=True)
    (C.OUT_ESM / "results_complex").mkdir(parents=True, exist_ok=True)

    rows = load_registry()
    constructs = [r for r in rows if r["kind"] == "construct"]
    targets = [r for r in rows if r["kind"] == "target"]

    # monomer FASTA (all 19)
    fasta = C.OUT_ESM / "esmfold2_input.fasta"
    with open(fasta, "w") as fh:
        for r in rows:
            fh.write(f">{r['id']}\n{r['sequence']}\n")

    # complex pairs (14 x 5) with the grafted-loop string for loop-pLDDT
    pairs = C.OUT_ESM / "esmfold2_complex_pairs.csv"
    with open(pairs, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["construct_id", "target_id",
                                           "construct_seq", "target_seq", "loop"])
        w.writeheader()
        for c in constructs:
            for t in targets:
                w.writerow({"construct_id": c["id"], "target_id": t["id"],
                            "construct_seq": c["sequence"], "target_seq": t["sequence"],
                            "loop": c["loop"]})

    shutil.copy(HERE / "run_esmfold2.py", C.OUT_ESM / "run_esmfold2.py")
    slurm = C.OUT_ESM / "run_esmfold2.slurm"
    slurm.write_text(SLURM_TEMPLATE.format(outdir=str(C.OUT_ESM)))

    print(f"Monomers: {len(rows)} -> {fasta.relative_to(C.REPO_ROOT)}")
    print(f"Complexes: {len(constructs)*len(targets)} -> {pairs.relative_to(C.REPO_ROOT)}")
    print(f"Runner -> {(C.OUT_ESM / 'run_esmfold2.py').relative_to(C.REPO_ROOT)}")
    print(f"SLURM  -> {slurm.relative_to(C.REPO_ROOT)}")
    print(f"\nSubmit on the cluster:  sbatch {slurm.relative_to(C.REPO_ROOT)}")


if __name__ == "__main__":
    main()
