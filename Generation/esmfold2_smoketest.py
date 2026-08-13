#!/usr/bin/env python3
"""
esmfold2_smoketest.py — prove ESMFold2 can actually FOLD before a long run starts.

Why this exists: `python -c "import esm"` succeeds in environments where the fold
itself fails, so the import check is not evidence. This folds a real two-chain
complex through the EXACT path the pipeline uses (subprocess -> $ESMFOLD2_PYTHON
-> score_with_esmfold2.py) and checks a finite esm_iptm came back. ~1-2 min.

It does NOT install, patch, or shim anything. It also distinguishes the two very
different reasons a fold can fail, because telling someone to "fix your
environment" when the real problem is a busy GPU sends them chasing a bug they
do not have:

    exit 0  folded — safe to launch
    exit 1  ENVIRONMENT fault: cannot load/run the model at all
    exit 2  RESOURCE fault: model loaded fine, GPU was out of memory

    python Generation/esmfold2_smoketest.py
    ESMFOLD2_PYTHON=/path/to/env/bin/python python Generation/esmfold2_smoketest.py
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

_HERE = Path(__file__).resolve().parent
SCORER = _HERE / "score_with_esmfold2.py"
TEMPLATE = _HERE / ".." / "Data" / "TIMP_Complexes" / "AF3_Templates" / "MMP2_TIMP3_AF3.pdb"

# Same default as iterative_refinement.ESMFOLD2_PYTHON — keep them in step.
ESMFOLD2_PYTHON = os.environ.get(
    "ESMFOLD2_PYTHON",
    str(Path.home() / "miniconda3" / "envs" / "esmfold2" / "bin" / "python"),
)

_AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


def chain_seqs(pdb: Path) -> dict:
    """{chain_id: one-letter sequence} from CA records. Dependency-free on purpose
    so this runs in the ESMFold2 env without pulling in atomworks/foundry."""
    seqs, seen = {}, set()
    for line in pdb.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        ch, resid = line[21], line[22:27]
        if (ch, resid) in seen:
            continue
        seen.add((ch, resid))
        seqs.setdefault(ch, []).append(_AA3.get(line[17:20].strip().upper(), "X"))
    return {c: "".join(v) for c, v in seqs.items()}


def _print_gpu_state() -> None:
    """Show who is holding GPU memory, so the real culprit is named rather than guessed."""
    for title, query in (
        ("  GPUs:", ["--query-gpu=index,name,memory.used,memory.total", "--format=csv,noheader"]),
        ("  Processes on GPU:",
         ["--query-compute-apps=pid,used_memory,process_name", "--format=csv,noheader"]),
    ):
        try:
            out = subprocess.run(["nvidia-smi"] + query, capture_output=True,
                                 text=True, timeout=20).stdout.strip()
        except Exception as exc:
            out = f"(nvidia-smi unavailable: {exc})"
        print(title)
        for ln in (out.splitlines() or ["(none)"]):
            print(f"    {ln}")


def main() -> int:
    if not SCORER.exists():
        print(f"FAIL: scorer not found at {SCORER}")
        return 1
    tpl = TEMPLATE.resolve()
    if not tpl.exists():
        print(f"FAIL: AF3 template not found at {tpl}")
        return 1

    seqs = chain_seqs(tpl)
    binder, target = seqs.get("A", ""), seqs.get("B", "")
    if len(binder) < 30 or len(target) < 30:
        print(f"FAIL: could not read two chains from {tpl.name} "
              f"(A={len(binder)} aa, B={len(target)} aa)")
        return 1

    tmp = Path(tempfile.gettempdir())
    in_csv, out_csv = tmp / "esm_smoke.csv", tmp / "esm_smoke_out.csv"
    out_csv.unlink(missing_ok=True)          # never PASS on a stale result
    in_csv.write_text(
        "design_id,target_name,full_seq,target_seq\n"
        f"SMOKE_TEST,MMP2,{binder},{target}\n"
    )
    print(f"Test complex: binder(chainA)={len(binder)} aa, "
          f"target(chainB)={len(target)} aa -> {in_csv}")

    if not os.access(ESMFOLD2_PYTHON, os.X_OK):
        print(f"FAIL: ESMFOLD2_PYTHON is not an executable python: {ESMFOLD2_PYTHON}\n"
              "  Launch from your ESMFold2 env, or export ESMFOLD2_PYTHON=/path/to/env/bin/python")
        return 1
    print(f"Folding with: {ESMFOLD2_PYTHON}")

    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    env.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
    try:
        proc = subprocess.run(
            [ESMFOLD2_PYTHON, str(SCORER), "--input", str(in_csv), "--out", str(out_csv)],
            capture_output=True, text=True, env=env,
        )
    except OSError as exc:
        print(f"FAIL: could not launch the scorer under {ESMFOLD2_PYTHON}: {exc}")
        return 1

    iptm = plddt = None
    if out_csv.exists():
        try:
            import csv
            with out_csv.open() as fh:
                for row in csv.DictReader(fh):
                    iptm, plddt = float(row["esm_iptm"]), float(row["esm_plddt"])
                    break
        except Exception as exc:
            print(f"  (could not parse {out_csv}: {exc})")

    if iptm is not None and iptm == iptm:      # not NaN
        print(f"PASS: ESMFold2 folded — esm_iptm={iptm:.3f}, esm_plddt={plddt:.1f}")
        print("Safe to launch the long run.")
        return 0

    combined = (proc.stdout or "") + (proc.stderr or "")

    # RESOURCE fault, not an environment fault. The model loaded and CUDA
    # initialised; the card is simply occupied. Never tell the user to fix an
    # environment that demonstrably works.
    if "out of memory" in combined.lower():
        print("FAIL: ESMFold2 ran OUT OF GPU MEMORY.")
        print("      This is NOT an environment fault — the model loaded and CUDA "
              "initialised fine.")
        for ln in combined.splitlines():
            if "out of memory" in ln.lower():
                print(f"  {ln.strip()[:320]}")
                break
        _print_gpu_state()
        print("\n  Free the RESOURCE, do not change the env:\n"
              "    * a leftover process is the usual culprit — killing a shell job does NOT\n"
              "      kill the scorer child it spawned. Find the PID above and kill it.\n"
              "    * or pin this job to a free card:  CUDA_VISIBLE_DEVICES=<n> <command>\n"
              "    * export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True (reduces fragmentation)\n"
              "    * never co-schedule two folding jobs on one card")
        return 2

    print("FAIL: ESMFold2 did not produce a score, and it was not an out-of-memory error.")
    print(f"  scorer exit code: {proc.returncode}")
    if proc.stdout.strip():
        print("  --- scorer stdout (tail) ---")
        print("  " + "\n  ".join(proc.stdout.strip().splitlines()[-15:]))
    if proc.stderr.strip():
        print("  --- scorer stderr (tail) ---")
        print("  " + "\n  ".join(proc.stderr.strip().splitlines()[-25:]))
    print("\n  Read the error above before assuming the environment is at fault. If it IS\n"
          "  an install/runtime problem, common causes are: the wrong env (imports `esm`\n"
          "  but cannot fold), a cuBLAS/cuequivariance mismatch (export\n"
          "  DISABLE_CUEQUIVARIANCE=1), or torch built for a different CUDA than the driver.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
