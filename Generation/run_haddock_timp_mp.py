import argparse
import glob
import math
import os
import shutil
import subprocess
from datetime import date

# ---------------------------------------------------------------------------
# Paths — run from the Generation/ directory; data and local dirs are ../
# ---------------------------------------------------------------------------
BASE_DIR = os.getcwd()
DATA_DIR = os.path.join(BASE_DIR, "../Data/Target_Crystal_Structures")

RUN_DATE = date.today().strftime("%Y-%m-%d")
LOCAL_DIR  = os.path.join(BASE_DIR, "../Local/haddock_timp_mp", RUN_DATE)
INPUT_DIR  = os.path.join(LOCAL_DIR, "inputs")
RUN_DIR    = os.path.join(LOCAL_DIR, "runs")
OUTPUT_DIR = os.path.join(BASE_DIR, "../Local/haddock_outputs", RUN_DATE)

# ---------------------------------------------------------------------------
# MP definitions (Name: (CIF_File, zinc-loop motif range))
# Motif ranges are the HExxHxxGxxH catalytic zinc loop in each AF catalytic
# domain CIF (sequential numbering from 1).
# ---------------------------------------------------------------------------
MPS = {
    "ADAM10": ("ADAM10cd_AF.cif", (164, 174)),
    "ADAM17": ("ADAM17cd_AF.cif", (187, 197)),
    "MMP10":  ("MMP10cd_AF.cif",  (119, 129)),
    "MMP2":   ("MMP2cd_AF.cif",   (121, 131)),
    "MMP3":   ("MMP3cd_AF.cif",   (119, 129)),
    "MMP9":   ("MMP9cd_AF.cif",   (120, 130)),
}

# TIMP3 active: CTCSP N-terminal loop (1-5) — coordinates directly with zinc
# TIMP3 passive: extended loop (6-10) — nearby context
TIMP3_ACTIVE  = (1, 5)
TIMP3_PASSIVE = (6, 10)
TIMP3_PDB = os.path.join(DATA_DIR, "TIMP3_Xray.pdb")

# Thresholds used by evaluate_complex()
SCORE_GOOD = -80.0   # HADDOCK score (more negative = better)
BSA_GOOD   = 1200.0  # buried surface area (Å²)
DIST_GOOD  = 9.0     # min CA-CA distance between motifs (Å); lower = closer

# ---------------------------------------------------------------------------
# Structure preparation
# ---------------------------------------------------------------------------

def prepare_structures():
    """Convert CIF to PDB and assign chain IDs (MP->A, TIMP3->B)."""
    for d in [INPUT_DIR, RUN_DIR, OUTPUT_DIR]:
        os.makedirs(d, exist_ok=True)

    timp3_out = os.path.join(INPUT_DIR, "TIMP3_B.pdb")
    print(f"Processing TIMP3 -> {timp3_out}")
    with open(TIMP3_PDB) as src, open(timp3_out, 'w') as dst:
        for line in src:
            if line.startswith(("ATOM", "HETATM")):
                dst.write(line[:21] + "B" + line[22:])
            else:
                dst.write(line)

    for name, (cif_file, _) in MPS.items():
        cif_path = os.path.join(DATA_DIR, cif_file)
        pdb_out  = os.path.join(INPUT_DIR, f"{name}_A.pdb")
        print(f"Processing {name} -> {pdb_out}")
        try:
            cmd = (
                f"conda run -n haddock python -c \""
                f"from Bio.PDB import MMCIFParser, PDBIO; "
                f"p = MMCIFParser(); s = p.get_structure('{name}', '{cif_path}'); "
                f"io = PDBIO(); io.set_structure(s); io.save('/tmp/tmp_{name}.pdb')\""
            )
            subprocess.run(cmd, shell=True, check=True)
            with open(f"/tmp/tmp_{name}.pdb") as src, open(pdb_out, 'w') as dst:
                for line in src:
                    if line.startswith(("ATOM", "HETATM")):
                        dst.write(line[:21] + "A" + line[22:])
                    else:
                        dst.write(line)
        except Exception as e:
            print(f"Error processing {name}: {e}")

# ---------------------------------------------------------------------------
# Restraint generation
# ---------------------------------------------------------------------------

def generate_tbl(mp_name, mp_range):
    """Write HADDOCK AIR restraint file for one MP target.

    Each active MP residue (HExxHxxGxxH loop) is ambiguously restrained to
    the combined TIMP3 active+passive region, and each TIMP3 active residue
    is ambiguously restrained back to the full MP loop.
    """
    tbl_path = os.path.join(INPUT_DIR, f"{mp_name}_restraints.tbl")

    mp_active_str = " or ".join(f"resid {i}" for i in range(mp_range[0], mp_range[1] + 1))
    timp_all_str  = " or ".join(f"resid {i}" for i in range(TIMP3_ACTIVE[0], TIMP3_PASSIVE[1] + 1))

    with open(tbl_path, 'w') as f:
        f.write("! AIR restraints: MP catalytic zinc loop (A) <-> TIMP3 N-terminal loop (B)\n!\n")
        f.write("! MP active -> TIMP3 (active + passive)\n")
        for i in range(mp_range[0], mp_range[1] + 1):
            f.write(f"assign (resid {i} and segid A) (({timp_all_str}) and segid B) 2.0 2.0 0.0\n")
        f.write("!\n! TIMP3 active -> MP active\n")
        for i in range(TIMP3_ACTIVE[0], TIMP3_ACTIVE[1] + 1):
            f.write(f"assign (resid {i} and segid B) (({mp_active_str}) and segid A) 2.0 2.0 0.0\n")

    return tbl_path

# ---------------------------------------------------------------------------
# Docking
# ---------------------------------------------------------------------------

def run_haddock(mp_name, ncores):
    """Generate HADDOCK3 config and run docking for one MP target."""
    work_dir = os.path.join(RUN_DIR, mp_name)
    if os.path.exists(work_dir):
        print(f"Cleaning up existing run directory for {mp_name}...")
        shutil.rmtree(work_dir)
    os.makedirs(work_dir, exist_ok=True)

    mp_pdb   = os.path.join(INPUT_DIR, f"{mp_name}_A.pdb")
    timp_pdb = os.path.join(INPUT_DIR, "TIMP3_B.pdb")
    tbl_path = generate_tbl(mp_name, MPS[mp_name][1])

    cfg_content = f"""
run_dir = "{work_dir}"
mode = "local"
ncores = {ncores}
postprocess = true
clean = false

molecules = [
    "{mp_pdb}",
    "{timp_pdb}"
]

[topoaa]

[rigidbody]
sampling = 500
ambig_fname = "{tbl_path}"

[seletop]
select = 200

[flexref]
ambig_fname = "{tbl_path}"

[emref]
ambig_fname = "{tbl_path}"

[clustfcc]

[seletopclusts]
top_clusters = 4
top_models = 5

[caprieval]
"""
    cfg_path = os.path.join(INPUT_DIR, f"{mp_name}_haddock3.cfg")
    with open(cfg_path, 'w') as f:
        f.write(cfg_content)

    print(f"Starting HADDOCK3 for {mp_name} ({ncores} cores)...")
    log_path = os.path.join(INPUT_DIR, f"{mp_name}_haddock3.log")
    with open(log_path, 'w') as log_file:
        subprocess.run(
            f"conda run -n haddock haddock3 {cfg_path}",
            shell=True, stdout=log_file, stderr=subprocess.STDOUT,
        )
    print(f"Finished {mp_name}. Log: {log_path}")

    # Extract best model from caprieval ranking
    try:
        capri_folders = sorted(glob.glob(os.path.join(work_dir, "*_caprieval")))
        if not capri_folders:
            print(f"No caprieval folder found for {mp_name}, skipping post-processing.")
            return

        last_capri   = capri_folders[-1]
        ranking_file = os.path.join(last_capri, "capri_ss.tsv")
        if not os.path.exists(ranking_file):
            print(f"capri_ss.tsv not found in {last_capri}")
            return

        with open(ranking_file) as f:
            lines = f.readlines()
        if len(lines) < 2:
            print(f"capri_ss.tsv has no model entries for {mp_name}")
            return

        best_model_rel  = lines[1].split('\t')[0].strip()
        best_model_path = os.path.normpath(os.path.join(last_capri, best_model_rel))
        if not os.path.exists(best_model_path):
            print(f"Best model file not found: {best_model_path}")
            return

        target_pdb = os.path.join(OUTPUT_DIR, f"{mp_name}_TIMP3_HADDOCK.pdb")
        shutil.copy(best_model_path, target_pdb)
        print(f"Best model for {mp_name} saved to {target_pdb}")

    except Exception as e:
        print(f"Post-processing failed for {mp_name}: {e}")

# ---------------------------------------------------------------------------
# Interface evaluation
# ---------------------------------------------------------------------------

def _ca_coords(pdb_path, chain, res_start, res_end):
    """Return list of (resnum, x, y, z) for CA atoms in a residue range."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            if line[21] != chain:
                continue
            resnum = int(line[22:26])
            if res_start <= resnum <= res_end:
                atoms.append((resnum, float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return atoms


def _parse_remark(pdb_path, keyword, field=0):
    """Return a float from a HADDOCK REMARK line matching keyword.

    field selects which comma-separated value to return (default 0 = first).
    Used for violations lines like '0, 0, 0, 0, 0, 0, 0'.
    """
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("REMARK") and keyword in line:
                try:
                    raw = line.split(":")[-1].strip()
                    return float(raw.split(",")[field].strip())
                except (ValueError, IndexError):
                    pass
    return None


def evaluate_complex(mp_name, pdb_path):
    """Print interface quality metrics for one docked complex.

    Metrics
    -------
    HADDOCK score   lower (more negative) is better; < -80 is good
    Buried SA (A^2) higher is better; > 1200 A^2 indicates a real interface
    Min CA-CA (A)   minimum CA-CA distance between the two motif regions;
                    < 9 A means the loops are in contact
    AIR violations  should be 0
    """
    mp_range = MPS[mp_name][1]

    score      = _parse_remark(pdb_path, "HADDOCK score:")
    bsa        = _parse_remark(pdb_path, "buried surface area:")
    violations = _parse_remark(pdb_path, "violations.:")

    mp_ca   = _ca_coords(pdb_path, "A", mp_range[0], mp_range[1])
    timp_ca = _ca_coords(pdb_path, "B", TIMP3_ACTIVE[0], TIMP3_ACTIVE[1])

    if mp_ca and timp_ca:
        min_dist = min(
            math.sqrt((a[1]-b[1])**2 + (a[2]-b[2])**2 + (a[3]-b[3])**2)
            for a in mp_ca for b in timp_ca
        )
    else:
        min_dist = None

    score_ok = score    is not None and score    < SCORE_GOOD
    bsa_ok   = bsa      is not None and bsa      > BSA_GOOD
    dist_ok  = min_dist is not None and min_dist < DIST_GOOD
    viol_ok  = violations == 0.0

    def fmt(val, fmt_str, ok):
        tag = "OK" if ok else "!!"
        return f"{val:{fmt_str}}  [{tag}]" if val is not None else "     N/A"

    print(f"\n{'='*55}")
    print(f"  {mp_name} + TIMP3")
    print(f"{'='*55}")
    print(f"  HADDOCK score    : {fmt(score,    '8.2f', score_ok)}  (threshold < {SCORE_GOOD})")
    print(f"  Buried SA (A^2)  : {fmt(bsa,      '8.1f', bsa_ok  )}  (threshold > {BSA_GOOD})")
    print(f"  Min CA-CA (A)    : {fmt(min_dist, '8.1f', dist_ok )}  (threshold < {DIST_GOOD})")
    print(f"  AIR violations   : {fmt(violations,'8.0f', viol_ok )}")

    passed = sum([score_ok, bsa_ok, dist_ok, viol_ok])
    print(f"  Overall          : {passed}/4 checks passed")
    print(f"{'='*55}")

    return {"mp": mp_name, "score": score, "bsa": bsa, "min_ca_ca": min_dist,
            "violations": violations, "checks_passed": passed}


def evaluate_all(output_dir):
    """Evaluate every *_TIMP3_HADDOCK.pdb found in output_dir."""
    pdbs = sorted(glob.glob(os.path.join(output_dir, "*_TIMP3_HADDOCK.pdb")))
    if not pdbs:
        print(f"No HADDOCK output PDBs found in {output_dir}")
        return

    results = []
    for pdb in pdbs:
        mp_name = os.path.basename(pdb).replace("_TIMP3_HADDOCK.pdb", "")
        if mp_name in MPS:
            results.append(evaluate_complex(mp_name, pdb))

    print(f"\n{'TARGET':<10} {'SCORE':>10} {'BSA':>10} {'MIN CA-CA':>10} {'CHECKS':>8}")
    print("-" * 52)
    for r in results:
        print(
            f"{r['mp']:<10} "
            f"{r['score']:>10.2f} "
            f"{r['bsa']:>10.1f} "
            f"{r['min_ca_ca']:>10.1f} "
            f"{r['checks_passed']:>6}/4"
        )

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run HADDOCK3 docking of TIMP3 against metalloprotease catalytic domains. "
                    "Run from the Generation/ directory."
    )
    parser.add_argument(
        "--ncores", type=int, default=14,
        help="Number of CPU cores for HADDOCK3 (default: 14)",
    )
    parser.add_argument(
        "--targets", nargs="+", choices=list(MPS.keys()), default=list(MPS.keys()),
        metavar="TARGET",
        help=f"Targets to dock. Choices: {list(MPS.keys())} (default: all)",
    )
    parser.add_argument(
        "--evaluate-only", metavar="OUTPUT_DIR", default=None,
        help="Skip docking and evaluate PDBs in the given directory.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if args.evaluate_only:
        evaluate_all(args.evaluate_only)
    else:
        prepare_structures()
        for mp in args.targets:
            run_haddock(mp, args.ncores)

        print("\nAll docking runs complete. Evaluating outputs...\n")
        evaluate_all(OUTPUT_DIR)
