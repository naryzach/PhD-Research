"""
Zinc-aware HADDOCK trial — runner (executes inside WSL, stdlib only).

Reads the manifest written by prep_zn_trial.py and runs the legacy HADDOCK_Outputs
protocol for six TIMP3:protease pairs in two arms:

  noZn : target as-is                       (faithful reproduction / control)
  Zn   : target + its catalytic Zn(2+)      (the hypothesis under test)

Everything else is held constant: same AF3-derived monomers, same AIRs, same
HADDOCK3 workflow (rigidbody 500 -> seletop 200 -> flexref -> emref -> clustfcc ->
seletopclusts -> caprieval).

CNS cannot handle the space in "Ryan Gustafson" and is painfully slow over /mnt/d,
so inputs are staged into a WSL-native scratch tree and only the best models are
copied back to Data/TIMP_Complexes/<arm folder>/.

Run (from WSL, inside the haddock env):
    python3 run_zn_trial.py --manifest <manifest.json> --repo-mnt /mnt/d/...
    optional: --arms Zn noZn --targets ADAM17 --jobs 8 --dry-run --force
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

HADDOCK_CMD = os.environ.get("HADDOCK3_CMD", "haddock3")
SCRATCH = Path(os.environ.get("HADDOCK_SCRATCH",
                              "/home/ryangustafson/haddock_scratch")) / "zn_trial"
ARM_DIR = {"noZn": "HADDOCK_indep_noZn", "Zn": "HADDOCK_indep_Zn",
           "ens": "HADDOCK_indep_ens"}

CFG = """
run_dir = "{run_dir}"
mode = "local"
ncores = {ncores}
postprocess = true
clean = false

molecules = [
    "{tgt}",
    "{timp}"
]

[topoaa]

[rigidbody]
sampling = {sampling}
iniseed = {seed}
ambig_fname = "{tbl}"
{unambig}
[seletop]
select = 200

[flexref]
ambig_fname = "{tbl}"
{unambig}{flexseg}
[emref]
ambig_fname = "{tbl}"
{unambig}
[clustfcc]

[seletopclusts]
top_clusters = 4
top_models = 5

[caprieval]
"""


def win_to_wsl(p) -> Path:
    """'D:\\Ryan Gustafson\\PhD-Research\\x' -> '/mnt/d/Ryan Gustafson/PhD-Research/x'."""
    s = str(p).replace("\\", "/")
    if len(s) > 1 and s[1] == ":":
        return Path(f"/mnt/{s[0].lower()}{s[2:]}")
    return Path(s)


def write_air(tbl: Path, motif: list[int], active: list[int], passive: list[int]) -> None:
    """Legacy AIRs: target zinc motif (segid A) <-> TIMP3 edge (segid B). Both arms."""
    mp = " or ".join(f"resid {i}" for i in motif)
    tp = " or ".join(f"resid {i}" for i in list(active) + list(passive))
    with open(tbl, "w") as f:
        f.write("! AIR: target zinc motif (A) <-> TIMP3 reactive edge (B)\n!\n")
        f.write("! target motif -> TIMP3 active+passive\n")
        for i in motif:
            f.write(f"assign (resid {i} and segid A) (({tp}) and segid B) 2.0 2.0 0.0\n")
        f.write("!\n! TIMP3 active -> target motif\n")
        for i in active:
            f.write(f"assign (resid {i} and segid B) (({mp}) and segid A) 2.0 2.0 0.0\n")


def write_zn_unambig(tbl: Path, zn_resid: int, his: list[int]) -> None:
    """Hold the Zn in its OWN histidine site (target-internal geometry only).

    Deliberately NOT restraining TIMP3 Cys1 to the zinc: whether the reactive edge
    finds the metal is exactly what this trial is testing.
    """
    with open(tbl, "w") as f:
        f.write("! Zn(2+) held in its own His triad — target-internal, no binder bias\n")
        for h in his:
            f.write(f"assign (segid A and resid {zn_resid} and name ZN) "
                    f"(segid A and resid {h} and (name NE2 or name ND1)) 2.1 0.3 0.3\n")


def best_model(run_dir: Path) -> Path | None:
    capri = sorted(glob.glob(str(run_dir / "*_caprieval")))
    if not capri:
        return None
    ranking = Path(capri[-1]) / "capri_ss.tsv"
    if not ranking.exists():
        return None
    lines = ranking.read_text().splitlines()
    if len(lines) < 2:
        return None
    cand = (Path(capri[-1]) / lines[1].split("\t")[0].strip()).resolve()
    return cand if cand.exists() else None


def run_pair(arm: str, t: dict, repo_mnt: Path, jobs: int,
             dry: bool, force: bool, seed: int) -> str:
    target = t["target"]
    out_dir = repo_mnt / "Data" / "TIMP_Complexes" / ARM_DIR[arm]
    out_pdb = out_dir / f"{target}_TIMP3_HADDOCK_s{seed}.pdb"
    if out_pdb.exists() and not force:
        return "skip(exists)"

    work = SCRATCH / arm / f"{target}_s{seed}"
    inp, run_dir = work / "inputs", work / "run"
    inp.mkdir(parents=True, exist_ok=True)
    if run_dir.exists():
        shutil.rmtree(run_dir, ignore_errors=True)

    # stage monomers from the D: prep into WSL-native, space-free scratch
    timp = inp / f"{target}_TIMP3.pdb"
    tgt = inp / f"{target}_target.pdb"
    if arm == "ens":
        src_timp = win_to_wsl(t["timp_pdb_ens"])
        src_tgt = win_to_wsl(t["target_pdb_ens"])
    else:
        src_timp = win_to_wsl(t["timp_pdb"])
        src_tgt = win_to_wsl(t["target_pdb_zn" if arm == "Zn" else "target_pdb"])
    for s, d in ((src_timp, timp), (src_tgt, tgt)):
        if not s.exists():
            return f"missing_input:{s.name}"
        shutil.copy(s, d)

    tbl = inp / f"{target}_air.tbl"
    write_air(tbl, t["motif_resids"], t["timp_active"], t["timp_passive"])
    unambig = ""
    if arm == "Zn":
        ztbl = inp / f"{target}_zn.tbl"
        write_zn_unambig(ztbl, t["zn_resid"], t["zn_his_resids"])
        unambig = f'unambig_fname = "{ztbl}"\n'

    # ensemble arm: more rigidbody sampling (3x3 conformer combinations) and
    # explicit semi-flexible segments so the interface can adapt during flexref
    sampling, flexseg = 500, ""
    if arm == "ens":
        sampling = 900
        ft, fi = t.get("flex_target"), t.get("flex_timp")
        if ft and fi:
            flexseg = (f"nseg1 = 1\nseg_sta_1_1 = {ft[0]}\nseg_end_1_1 = {ft[1]}\n"
                       f"nseg2 = 1\nseg_sta_2_1 = {fi[0]}\nseg_end_2_1 = {fi[1]}\n")

    cfg = inp / f"{target}_{arm}_s{seed}.cfg"
    cfg.write_text(CFG.format(run_dir=run_dir, ncores=jobs, tgt=tgt, timp=timp,
                              tbl=tbl, unambig=unambig, seed=seed,
                              sampling=sampling, flexseg=flexseg))
    if dry:
        return "dry-run(prepared)"

    log = inp / f"{target}_{arm}.log"
    with open(log, "w") as lf:
        subprocess.run(f"{HADDOCK_CMD} {cfg}", shell=True, stdout=lf,
                       stderr=subprocess.STDOUT)
    bm = best_model(run_dir)
    if bm is None:
        return "FAILED(no ranked model)"
    out_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(bm, out_pdb)
    shutil.rmtree(run_dir, ignore_errors=True)
    return f"ok -> {out_pdb.name}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--repo-mnt", required=True,
                    help="repo root as seen from WSL, e.g. '/mnt/d/Ryan Gustafson/PhD-Research'")
    ap.add_argument("--arms", nargs="*", default=["noZn", "Zn"])
    ap.add_argument("--targets", nargs="*", default=None)
    ap.add_argument("--seeds", nargs="*", type=int, default=[917, 424242, 20260716],
                    help="one HADDOCK run per seed = independent replicate")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    man = json.loads(Path(a.manifest).read_text())
    repo = Path(a.repo_mnt)
    rows = [t for t in man["targets"] if t.get("zn_status") == "ok"]
    if a.targets:
        rows = [t for t in rows if t["target"] in a.targets]
    print(f"{len(rows)} targets x {len(a.arms)} arms x {len(a.seeds)} seeds "
          f"= {len(rows) * len(a.arms) * len(a.seeds)} runs; scratch={SCRATCH}")
    for seed in a.seeds:
        for arm in a.arms:
            for t in rows:
                print(f"[seed {seed}][{arm}] {t['target']:<8}", end=" ", flush=True)
                try:
                    print(run_pair(arm, t, repo, a.jobs, a.dry_run, a.force, seed))
                except Exception as e:
                    print(f"ERROR {e}")


if __name__ == "__main__":
    main()
