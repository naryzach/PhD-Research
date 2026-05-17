import os
import subprocess
import glob
from datetime import date

# Run from the PhD-Research root directory; all paths are relative to it.
BASE_DIR = os.getcwd()
DATA_DIR = os.path.join(BASE_DIR, "Data/Target_Crystal_Structures")

RUN_DATE = date.today().strftime("%Y-%m-%d")
LOCAL_DIR = os.path.join(BASE_DIR, "Local/haddock_timp_mp", RUN_DATE)
INPUT_DIR = os.path.join(LOCAL_DIR, "inputs")
RUN_DIR = os.path.join(LOCAL_DIR, "runs")
OUTPUT_DIR = os.path.join(BASE_DIR, "Local/haddock_outputs", RUN_DATE)

# Ensure directories exist
for d in [INPUT_DIR, RUN_DIR, OUTPUT_DIR]:
    os.makedirs(d, exist_ok=True)

# MP definitions (Name: (CIF_File, Motif_Range))
mps = {
    "ADAM10": ("ADAM10cd_AF.cif", (164, 174)),
    "ADAM17": ("ADAM17cd_AF.cif", (187, 197)),
    "MMP10": ("MMP10cd_AF.cif", (119, 129)),
    "MMP2": ("MMP2cd_AF.cif", (121, 131)),
    "MMP3": ("MMP3cd_AF.cif", (119, 129)),
    "MMP9": ("MMP9cd_AF.cif", (120, 130)),
}

NCORES = 14

# TIMP3 active residues: CTCSP loop (1-5) — coordinates directly with catalytic zinc
# TIMP3 passive residues: extended N-terminal loop (6-10) — nearby but not directly in active site
TIMP3_ACTIVE = (1, 5)
TIMP3_PASSIVE = (6, 10)
TIMP3_PDB = os.path.join(DATA_DIR, "TIMP3_Xray.pdb")

def prepare_structures():
    """Convert CIF to PDB and re-chain."""
    # Process TIMP3
    timp3_out = os.path.join(INPUT_DIR, "TIMP3_B.pdb")
    print(f"Processing TIMP3 -> {timp3_out}")
    with open(TIMP3_PDB, 'r') as f, open(timp3_out, 'w') as out:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Change chain ID (column 22, 1-indexed) to B
                new_line = line[:21] + "B" + line[22:]
                out.write(new_line)
            else:
                out.write(line)

    # Process MPs
    for name, (cif_file, _) in mps.items():
        cif_path = os.path.join(DATA_DIR, cif_file)
        pdb_out = os.path.join(INPUT_DIR, f"{name}_A.pdb")
        print(f"Processing {name} -> {pdb_out}")
        
        # Use a simple script to extract ATOM records and set chain to A
        # For CIF, we can use a basic parser or just grep ATOM and reformat
        # Easier to use a python script with a library if available, but let's do a simple text conversion
        # since we know the format of AF cif files.
        try:
            # Actually, I'll use a temporary python script to do this properly if needed,
            # but let's try a direct conversion first.
            cmd = f"conda run -n haddock python -c \"from Bio.PDB import MMCIFParser, PDBIO; parser = MMCIFParser(); structure = parser.get_structure('{name}', '{cif_path}'); io = PDBIO(); io.set_structure(structure); io.save('/tmp/tmp.pdb')\""
            subprocess.run(cmd, shell=True, check=True)
            
            with open('/tmp/tmp.pdb', 'r') as f, open(pdb_out, 'w') as out:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        # Chain to A
                        new_line = line[:21] + "A" + line[22:]
                        out.write(new_line)
                    else:
                        out.write(line)
        except Exception as e:
            print(f"Error processing {name}: {e}")

def generate_tbl(mp_name, mp_range):
    """Generate HADDOCK TBL file for AIRs.

    Uses the standard HADDOCK ambiguous interaction restraint (AIR) format:
    each active residue of molecule 1 must be within 2+2=4 Å of any atom of
    (active + passive) residues of molecule 2, and vice versa.
    Passive residues provide one-sided restraints only (they are targets, not sources).
    """
    tbl_path = os.path.join(INPUT_DIR, f"{mp_name}_restraints.tbl")

    # The full set of MP motif residues spans the HExxHxxGxxH zinc-binding loop.
    # All are treated as active since the loop is the direct interface site.
    mp_active_str = " or ".join([f"resid {i}" for i in range(mp_range[0], mp_range[1] + 1)])

    # Combined TIMP3 target (active + passive) for MP active restraints
    timp_all_str = " or ".join(
        [f"resid {i}" for i in range(TIMP3_ACTIVE[0], TIMP3_PASSIVE[1] + 1)]
    )

    with open(tbl_path, 'w') as f:
        f.write("! AIR restraints: MP catalytic zinc loop (A) <-> TIMP3 N-terminal loop (B)\n")
        f.write("!\n")

        # Each active MP residue restrained to any atom of TIMP3 active+passive
        f.write("! MP active -> TIMP3 (active + passive)\n")
        for i in range(mp_range[0], mp_range[1] + 1):
            f.write(
                f"assign (resid {i} and segid A) "
                f"(({timp_all_str}) and segid B) 2.0 2.0 0.0\n"
            )

        f.write("!\n")
        # Each TIMP3 active residue restrained to any atom of MP active loop
        f.write("! TIMP3 active -> MP active\n")
        for i in range(TIMP3_ACTIVE[0], TIMP3_ACTIVE[1] + 1):
            f.write(
                f"assign (resid {i} and segid B) "
                f"(({mp_active_str}) and segid A) 2.0 2.0 0.0\n"
            )

    return tbl_path

def run_haddock(mp_name):
    """Generate config and run HADDOCK3."""
    work_dir = os.path.join(RUN_DIR, mp_name)
    
    # Clean up existing run directory to avoid HADDOCK3 "directory not empty" error
    if os.path.exists(work_dir):
        import shutil
        print(f"Cleaning up existing run directory for {mp_name}...")
        shutil.rmtree(work_dir)
        
    os.makedirs(work_dir, exist_ok=True)
    
    mp_pdb = os.path.join(INPUT_DIR, f"{mp_name}_A.pdb")
    timp_pdb = os.path.join(INPUT_DIR, "TIMP3_B.pdb")
    tbl_path = generate_tbl(mp_name, mps[mp_name][1])
    
    cfg_content = f"""
run_dir = "{work_dir}"
mode = "local"
ncores = {NCORES}
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
# Restraints must carry through flexible refinement so the loop interface is preserved
ambig_fname = "{tbl_path}"

[emref]
# Restraints must carry through EM refinement so the loop interface is preserved
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
        
    print(f"Starting HADDOCK3 for {mp_name}...")
    log_path = os.path.join(INPUT_DIR, f"{mp_name}_haddock3.log")
    with open(log_path, 'w') as log_file:
        subprocess.run(f"conda run -n haddock haddock3 {cfg_path}", shell=True, stdout=log_file, stderr=subprocess.STDOUT)
    print(f"Finished {mp_name}. Log: {log_path}")
    
    # Post-processing: Extract best model from the final caprieval ranking
    try:
        import shutil

        capri_folders = sorted(glob.glob(os.path.join(work_dir, "*_caprieval")))
        if not capri_folders:
            print(f"No caprieval folder found for {mp_name}, skipping post-processing.")
            return

        last_capri = capri_folders[-1]
        ranking_file = os.path.join(last_capri, "capri_ss.tsv")
        if not os.path.exists(ranking_file):
            print(f"capri_ss.tsv not found in {last_capri}")
            return

        with open(ranking_file) as f:
            lines = f.readlines()

        if len(lines) < 2:
            print(f"capri_ss.tsv has no model entries for {mp_name}")
            return

        # capri_ss.tsv column 0 is the model path relative to the caprieval folder
        best_model_rel = lines[1].split('\t')[0].strip()
        best_model_path = os.path.normpath(os.path.join(last_capri, best_model_rel))

        if not os.path.exists(best_model_path):
            print(f"Best model file not found: {best_model_path}")
            return

        target_pdb = os.path.join(OUTPUT_DIR, f"{mp_name}_TIMP3_HADDOCK.pdb")
        shutil.copy(best_model_path, target_pdb)
        print(f"Best model for {mp_name} saved to {target_pdb}")

    except Exception as e:
        print(f"Post-processing failed for {mp_name}: {e}")

if __name__ == "__main__":
    prepare_structures()
    for mp in mps:
        run_haddock(mp)
    print("All docking runs and post-processing complete.")
