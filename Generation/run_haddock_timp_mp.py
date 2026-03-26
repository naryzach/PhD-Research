import os
import subprocess
import glob

# Configuration
BASE_DIR = "/home/ryangustafson/Documents/GitHubProj/PhD-Research"
DATA_DIR = os.path.join(BASE_DIR, "Data/Target_Crystal_Structures")
LOCAL_DIR = os.path.join(BASE_DIR, "Local/haddock_timp_mp")
INPUT_DIR = os.path.join(LOCAL_DIR, "inputs")
RUN_DIR = os.path.join(LOCAL_DIR, "runs")
OUTPUT_DIR = os.path.join(BASE_DIR, "Data/HADDOCK_Outputs")

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

TIMP3_MOTIF = (1, 5)  # CTCSP loop
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
    """Generate HADDOCK TBL file for AIRs."""
    tbl_path = os.path.join(INPUT_DIR, f"{mp_name}_restraints.tbl")
    with open(tbl_path, 'w') as f:
        # Define AIRs: MP motif residues are active, TIMP3 motif residues are active.
        # Format: assign (selection1) (selection2) distance upper lower
        
        # Selection 1: MP Motif (Chain A)
        mp_resid_str = " or ".join([f"resid {i}" for i in range(mp_range[0], mp_range[1] + 1)])
        # Selection 2: TIMP3 Motif (Chain B)
        timp_resid_str = " or ".join([f"resid {i}" for i in range(TIMP3_MOTIF[0], TIMP3_MOTIF[1] + 1)])
        
        # In HADDOCK, "active" residues on protein 1 should be restrained to "all residues (active+passive) on protein 2".
        # For simplicity and strong bias, we restrain the two motifs together.
        for i in range(mp_range[0], mp_range[1] + 1):
            f.write(f"assign (resid {i} and segid A) (({timp_resid_str}) and segid B) 2.0 2.0 0.0\n")
        
        for i in range(TIMP3_MOTIF[0], TIMP3_MOTIF[1] + 1):
            f.write(f"assign (resid {i} and segid B) (({mp_resid_str}) and segid A) 2.0 2.0 0.0\n")
            
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
select = 100

[flexref]

[emref]

[clustfcc]

[seletopclusts]
top_clusters = 1
top_models = 5

[caprieval]
reference_fname = "{mp_pdb}"
"""
    cfg_path = os.path.join(INPUT_DIR, f"{mp_name}_haddock3.cfg")
    with open(cfg_path, 'w') as f:
        f.write(cfg_content)
        
    print(f"Starting HADDOCK3 for {mp_name}...")
    log_path = os.path.join(INPUT_DIR, f"{mp_name}_haddock3.log")
    with open(log_path, 'w') as log_file:
        subprocess.run(f"conda run -n haddock haddock3 {cfg_path}", shell=True, stdout=log_file, stderr=subprocess.STDOUT)
    print(f"Finished {mp_name}. Log: {log_path}")
    
    # Post-processing: Extract best model and convert to CIF
    try:
        # HADDOCK3 output structure: work_dir/07_caprieval/capri_ss.tsv (or similar)
        # However, it's easier to find the best model in the last caprieval folder.
        capri_folders = sorted(glob.glob(os.path.join(work_dir, "*_caprieval")))
        if capri_folders:
            last_capri = capri_folders[-1]
            # The top models are copied to the cluster/model folders or kept in the previous step's output.
            # HADDOCK3 usually has a 'best_models' or similar in the run folder if seletop was used.
            # Let's look for the ranking.
            ranking_file = os.path.join(last_capri, "capri_ss.tsv")
            if os.path.exists(ranking_file):
                with open(ranking_file, 'r') as f:
                    lines = f.readlines()
                    if len(lines) > 1:
                        # First line is header, second line is best model
                        best_model_rel = lines[1].split('\t')[0]
                        best_model_path = os.path.join(work_dir, best_model_rel)
                        
                        target_pdb = os.path.join(OUTPUT_DIR, f"{mp_name}_TIMP3_complex.pdb")
                        target_cif = os.path.join(OUTPUT_DIR, f"{mp_name}_TIMP3_complex.cif")
                        
                        import shutil
                        shutil.copy(best_model_path, target_pdb)
                        print(f"Best model for {mp_name} saved to {target_pdb}")
                        
                        # Convert to CIF
                        cmd = f"conda run -n haddock python -c \"from Bio.PDB import PDBParser, MMCIFIO; parser = PDBParser(); structure = parser.get_structure('{mp_name}', '{target_pdb}'); io = MMCIFIO(); io.set_structure(structure); io.save('{target_cif}')\""
                        subprocess.run(cmd, shell=True, check=True)
                        print(f"Converted {mp_name} to CIF: {target_cif}")
    except Exception as e:
        print(f"Post-processing failed for {mp_name}: {e}")

if __name__ == "__main__":
    prepare_structures()
    for mp in mps:
        run_haddock(mp)
    print("All docking runs and post-processing complete.")
