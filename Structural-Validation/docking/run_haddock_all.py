"""
Dock every TIMP3 construct against every target with HADDOCK3.

Generalises Generation/run_haddock_timp_mp.py from "WT TIMP3 vs 6 proteases" to
"14 constructs x 5 targets", driving the AIR restraints off each target's
catalytic zinc motif (HExxHxxGxxH), located by sequence in the actual input
model so residue numbering never has to be hand-maintained.

Docking inputs are the monomer models from a chosen fold source (default AF3;
falls back to the repo's AF catalytic-domain CIFs / crystals for targets when
an AF3 target model isn't present). Chain A = target, chain B = construct.

Runs HADDOCK3 in the `haddock` conda env, on your side. Use --dry-run to build
every restraint file + config without executing (works anywhere), then launch
for real on the cluster.

Examples:
  python Structural-Validation/docking/run_haddock_all.py --dry-run
  python Structural-Validation/docking/run_haddock_all.py --fold-source AF3 --jobs 8
  python Structural-Validation/docking/run_haddock_all.py --pairs AB_3:ADAM17,C_16:MMP9
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO, PDBParser
from Bio.PDB.Polypeptide import is_aa

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
import config as C  # noqa: E402
import model_registry as MR  # noqa: E402
import structure_io as sio  # noqa: E402

warnings.filterwarnings("ignore")

RIGIDBODY_SAMPLING = 500
SELECT_TOP = 200


def target_model(target_id: str, target_src: str) -> Path | None:
    """Docking target from the requested source.

    Crystal   -> the extracted CD-only clean crystal (crystal/clean_targets/)
    AF3/ESM   -> that fold source's target monomer, with graceful fallback to the
                 repo AF CD CIF then the raw crystal so a track never silently
                 vanishes when a fold model isn't present yet.
    """
    if target_src == "Crystal":
        ct = MR.crystal_target(target_id)
        if ct:
            return ct
    p = MR.monomer_model(target_id, target_src) if target_src in C.FOLDERS else None
    if p:
        return p
    af = C.CRYSTAL_DIR / C.TARGETS[target_id]["af_cif"]
    if af.exists():
        return af
    cr = C.CRYSTAL_DIR / C.TARGETS[target_id]["crystal"]
    return cr if cr.exists() else None


def to_pdb_chain(src: Path, dst: Path, chain_id: str) -> None:
    """Convert any PDB/CIF (incl. minimal ESMFold2 CIFs) to a single-chain PDB
    with a forced chain id."""
    s = sio.load(src)   # tolerant loader: falls back for minimal CIFs
    model = list(s)[0]
    # keep the longest protein chain only
    best = max(model, key=lambda ch: sum(1 for r in ch if is_aa(r, standard=True)))
    for r in list(best):
        if not is_aa(r, standard=True):
            best.detach_child(r.id)
    best.id = chain_id
    for ch in list(model):
        if ch.id != chain_id:
            model.detach_child(ch.id)
    io = PDBIO()
    io.set_structure(model)
    io.save(str(dst))


def write_restraints(tbl: Path, target_motif_resids, timp_all, timp_active) -> None:
    mp_active = " or ".join(f"resid {i}" for i in target_motif_resids)
    timp_all_s = " or ".join(f"resid {i}" for i in timp_all)
    with open(tbl, "w") as f:
        f.write("! AIR: target zinc motif (A) <-> TIMP3 reactive edge (B)\n")
        for i in target_motif_resids:
            f.write(f"assign (resid {i} and segid A) (({timp_all_s}) and segid B) 2.0 2.0 0.0\n")
        for i in timp_active:
            f.write(f"assign (resid {i} and segid B) (({mp_active}) and segid A) 2.0 2.0 0.0\n")


def make_config(cfg: Path, run_dir: Path, tgt_pdb: Path, con_pdb: Path,
                tbl: Path, ncores: int) -> None:
    cfg.write_text(f"""
run_dir = "{run_dir}"
mode = "local"
ncores = {ncores}
postprocess = true
clean = false

molecules = [
    "{tgt_pdb}",
    "{con_pdb}"
]

[topoaa]

[rigidbody]
sampling = {RIGIDBODY_SAMPLING}
ambig_fname = "{tbl}"

[seletop]
select = {SELECT_TOP}

[flexref]
ambig_fname = "{tbl}"

[emref]
ambig_fname = "{tbl}"

[clustfcc]

[seletopclusts]
top_clusters = 4
top_models = 5

[caprieval]
""")


def prepare_pair(construct, target, construct_src, target_src, input_dir) -> dict | None:
    con_model = MR.monomer_model(construct["id"], construct_src)
    tgt_model = target_model(target["id"], target_src)
    status = {"construct_id": construct["id"], "target_id": target["id"],
              "construct_src": construct_src, "target_src": target_src}
    if con_model is None:
        status["prep"] = f"missing_construct_model({construct_src})"
        return status
    if tgt_model is None:
        status["prep"] = "missing_target_model"
        return status

    pair = f"{construct['id']}__{target['id']}"
    tgt_pdb = input_dir / f"{pair}_target_A.pdb"
    con_pdb = input_dir / f"{pair}_construct_B.pdb"
    to_pdb_chain(tgt_model, tgt_pdb, "A")
    to_pdb_chain(con_model, con_pdb, "B")

    # locate the zinc motif in the prepared target chain. Prefer an exact
    # substring; otherwise map the motif from the reference CD sequence onto the
    # model by alignment (handles AF-CIF / crystal fallbacks whose sequence
    # differs slightly from the registry CD).
    tchain = sio.longest_chain(sio.get_chains(tgt_pdb))
    motif_resids = sio.zinc_motif_resids(tchain)
    if not motif_resids:
        status["prep"] = "motif_not_found_in_target_model"
        return status

    cchain = sio.longest_chain(sio.get_chains(con_pdb))
    timp_active = cchain.resids[C.TIMP3_ACTIVE[0] - 1: C.TIMP3_ACTIVE[1]]
    timp_all = cchain.resids[C.TIMP3_ACTIVE[0] - 1: C.TIMP3_PASSIVE[1]]

    tbl = input_dir / f"{pair}_restraints.tbl"
    write_restraints(tbl, motif_resids, timp_all, timp_active)
    status.update(prep="ready", target_pdb=str(tgt_pdb), construct_pdb=str(con_pdb),
                  restraints=str(tbl), motif_resids=str(motif_resids))
    return status


HADDOCK_CMD = os.environ.get("HADDOCK3_CMD", "haddock3")  # run inside the env
FORCE = False   # set from --force; when False, existing best models are skipped


def run_one(pair: str, cfg: Path, run_dir: Path, best_dir: Path) -> None:
    log = cfg.with_suffix(".log")
    with open(log, "w") as lf:
        subprocess.run(f"{HADDOCK_CMD} {cfg}",
                       shell=True, stdout=lf, stderr=subprocess.STDOUT)
    import glob
    import os
    capri = sorted(glob.glob(os.path.join(str(run_dir), "*_caprieval")))
    if not capri:
        print(f"  {pair}: no caprieval output"); return
    ranking = Path(capri[-1]) / "capri_ss.tsv"
    if not ranking.exists():
        print(f"  {pair}: no capri_ss.tsv"); return
    lines = ranking.read_text().splitlines()
    if len(lines) < 2:
        print(f"  {pair}: no models ranked"); return
    best_rel = lines[1].split("\t")[0].strip()
    best = (Path(capri[-1]) / best_rel).resolve()
    if best.exists():
        best_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(best, best_dir / f"{pair}_HADDOCK.pdb")
        print(f"  {pair}: best model -> {pair}_HADDOCK.pdb")


def run_track(construct_src, target_src, constructs, targets, want,
              jobs, dry_run) -> tuple[int, int]:
    """Prepare (and optionally run) every pair for one docking track."""
    track = MR.dock_track(construct_src, target_src)
    base = C.OUT_DOCK / track
    input_dir, run_root, best_dir = base / "inputs", base / "runs", base / "best_models"
    for d in (input_dir, run_root, best_dir):
        d.mkdir(parents=True, exist_ok=True)

    print(f"\n### TRACK {track}  (construct={construct_src}, target={target_src})")
    prep_rows = []
    for c in constructs:
        for t in targets:
            if want and (c["id"], t["id"]) not in want:
                continue
            st = prepare_pair(c, t, construct_src, target_src, input_dir)
            prep_rows.append(st)
            if st.get("prep") != "ready":
                print(f"  SKIP {c['id']}__{t['id']}: {st.get('prep')}")
                continue
            pair = f"{c['id']}__{t['id']}"
            cfg = input_dir / f"{pair}_haddock3.cfg"
            run_dir = run_root / pair
            make_config(cfg, run_dir, Path(st["target_pdb"]),
                        Path(st["construct_pdb"]), Path(st["restraints"]), jobs)
            if dry_run:
                print(f"  PREP {pair}: ready (motif {st['motif_resids']})")
            elif (best_dir / f"{pair}_HADDOCK.pdb").exists() and not FORCE:
                print(f"  SKIP {pair}: already docked (resume)")
            else:
                if run_dir.exists():
                    shutil.rmtree(run_dir)
                print(f"  RUN  {pair} ...")
                run_one(pair, cfg, run_dir, best_dir)

    with open(base / "prep_manifest.csv", "w", newline="") as fh:
        keys = sorted({k for r in prep_rows for k in r})
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(prep_rows)
    ready = sum(1 for r in prep_rows if r.get("prep") == "ready")
    return ready, len(prep_rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", default=C.DOCK_MATRIX_DEFAULT,
                    choices=list(C.DOCK_MATRIX), help="preset (construct,target) tracks")
    ap.add_argument("--construct-source", choices=C.CONSTRUCT_SOURCES,
                    help="single-track override: construct fold source")
    ap.add_argument("--target-source", choices=C.TARGET_SOURCES,
                    help="single-track override: target structure source")
    ap.add_argument("--jobs", type=int, default=4, help="HADDOCK ncores per run")
    ap.add_argument("--pairs", default="", help="comma list construct:target to restrict")
    ap.add_argument("--dry-run", action="store_true",
                    help="build inputs/restraints/configs only; do not run HADDOCK")
    ap.add_argument("--force", action="store_true",
                    help="re-dock pairs even if a best model already exists")
    args = ap.parse_args()
    global FORCE
    FORCE = args.force

    reg = {r["id"]: r for r in MR.load_registry()}
    constructs = [r for r in reg.values() if r["kind"] == "construct"]
    targets = [r for r in reg.values() if r["kind"] == "target"]
    want = {tuple(p.split(":")) for p in args.pairs.split(",")} if args.pairs else None

    if args.construct_source and args.target_source:
        tracks = [(args.construct_source, args.target_source)]
    else:
        tracks = C.DOCK_MATRIX[args.matrix]

    print(f"Docking matrix: {tracks}  ({len(tracks)} track(s) x up to "
          f"{len(constructs)*len(targets)} pairs)")
    total_ready = total = 0
    for cs, ts in tracks:
        r, n = run_track(cs, ts, constructs, targets, want, args.jobs, args.dry_run)
        total_ready += r; total += n

    print(f"\n{total_ready}/{total} pairs prepared across {len(tracks)} track(s) "
          f"({'dry-run, nothing executed' if args.dry_run else 'HADDOCK executed'}).")
    print(f"Tracks under {C.OUT_DOCK.relative_to(C.REPO_ROOT)}/<construct>__<target>/")


if __name__ == "__main__":
    main()
