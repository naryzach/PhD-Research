import argparse
import os
import sys
import subprocess
import glob
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def get_sequence(pdb_path, chain_id=None):
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("temp", pdb_path)
    except:
        return None
        
    sequences = {}
    
    for model in structure:
        for chain in model:
            if chain_id and chain.id != chain_id:
                continue
            
            seq = []
            for residue in chain:
                if is_aa(residue, standard=True):
                    resname = residue.resname.upper()
                    seq.append(THREE_TO_ONE.get(resname, 'X'))
            if seq:
                sequences[chain.id] = "".join(seq)
    
    # Allowed ions for AlphaFold
    ALLOWED_IONS = {"MG", "ZN", "CL", "CA", "NA", "MN", "K", "FE", "CU", "CO"}
    
    # Check for Ligand Chains
    for model in structure:
        for chain in model:
            if chain.id not in sequences or len(sequences[chain.id]) == 0:
                residues = list(chain.get_residues())
                if residues:
                    residue = residues[0]
                    resname = residue.resname.strip().upper()
                    if resname in ALLOWED_IONS:
                        sequences[chain.id] = resname
                    else:
                        found_ion = None
                        for atom in residue:
                            if atom.element.upper() in ALLOWED_IONS:
                                found_ion = atom.element.upper()
                                break
                        if found_ion:
                            sequences[chain.id] = found_ion
                        else:
                            if resname != 'HOH':
                                sequences[chain.id] = resname
        break
    return sequences

def get_fixed_residues(pdb_path):
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("temp", pdb_path)
    except:
        return {}
        
    # Allowed ions list (same as above)
    ALLOWED_IONS = {"MG", "ZN", "CL", "CA", "NA", "MN", "K", "FE", "CU", "CO", "CD", "NI", "HG"}
    
    metal_atoms = []
    # Find metals
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname.strip().upper() in ALLOWED_IONS:
                    for atom in residue:
                        metal_atoms.append(atom)
                else:
                    try:
                        for atom in residue:
                            if atom.element.upper() in ALLOWED_IONS:
                                metal_atoms.append(atom)
                    except:
                        pass
                            
    if not metal_atoms:
        return {}
    
    fixed_residues = {} # chain -> list of ints
    
    # Find neighbors
    threshold = 4.0
    for model in structure:
        for chain in model:
            if chain.id not in fixed_residues:
                fixed_residues[chain.id] = []
                
            for residue in chain:
                # Skip if it is the metal itself
                if residue.resname.strip().upper() in ALLOWED_IONS:
                    continue
                if not is_aa(residue, standard=True):
                    continue
                    
                is_contact = False
                for atom in residue:
                    vec = atom.get_vector()
                    for metal in metal_atoms:
                        metal_vec = metal.get_vector()
                        if (vec - metal_vec).norm() < threshold:
                            is_contact = True
                            break
                    if is_contact: break
                
                if is_contact:
                    # MPNN uses 1-based indexing matching PDB
                    fixed_residues[chain.id].append(residue.id[1])

    return fixed_residues

def main():
    parser = argparse.ArgumentParser(description="Run LigandMPNN on Design Folder")
    parser.add_argument("--pred_dir", default="../Local/Metal_Predictions", help="Predictions Directory")
    parser.add_argument("--lmpnn_path", default="../Tools/LigandMPNN", help="Path to LigandMPNN")
    parser.add_argument("--num_seqs", type=int, default=5, help="Number of sequences to generate per design")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")
    parser.add_argument("--dry_run", action="store_true", help="Print command without executing")
    args = parser.parse_args()
    
    print(f"Scanning {args.pred_dir}...")
    pdb_files = glob.glob(os.path.join(args.pred_dir, "**", "*.pdb"), recursive=True)
    
    # Filter for RFdiffusion outputs
    design_files = [f for f in pdb_files if "rfdiffusion" in f and ("scaffold" in os.path.basename(f) or "design" in os.path.basename(f)) and "traj" not in os.path.basename(f) and "unidealized" not in f]
    
    print(f"Found {len(design_files)} designs for LigandMPNN.")

    for i, pdb in enumerate(design_files):
        # Output folder: sibling of rfdiffusion
        dirname = os.path.dirname(pdb)
        parent_dir = os.path.dirname(dirname) 
        motif_name = os.path.basename(parent_dir)
        
        # Use 'ligandmpnn' folder instead of 'proteinmpnn'
        mpnn_out = os.path.join(parent_dir, "ligandmpnn")
        os.makedirs(mpnn_out, exist_ok=True)
        
        # Determine chains to design (needed for checking expected output name if ligand present)
        seqs_dict = get_sequence(pdb)
        if seqs_dict is None: seqs_dict = {}
        all_chains = sorted(seqs_dict.keys())
        
        protein_chains = [c for c in all_chains if len(seqs_dict[c]) > 5]
        ligand_chains = [c for c in all_chains if len(seqs_dict[c]) <= 5]
        
        basename = os.path.basename(pdb).replace(".pdb", "")
        
        # Output naming
        # ...
        
        # Path parsing for logging: .../Category/Ion/Motif/rfdiffusion/file.pdb
        # parent_dir = .../Motif
        ion_dir = os.path.dirname(parent_dir)
        ion = os.path.basename(ion_dir)
        category_dir = os.path.dirname(ion_dir)
        category = os.path.basename(category_dir)

        design_chains_str = ",".join(protein_chains) # LigandMPNN comma sep
        
        # Expected output usually: {out_folder}/seqs/{basename}.fa
        expected_output = os.path.join(mpnn_out, "seqs", f"{basename}.fa")
        if os.path.exists(expected_output) and not args.overwrite:
            print(f"[{i+1}/{len(design_files)}] [{category}] [{ion}] [{motif_name}] Skipping {basename} (Already exists)")
            continue
        
        # Fixed Residues
        fixed_res_dict = get_fixed_residues(pdb)
        fixed_str = ""
        fixed_list = []
        if any(v for v in fixed_res_dict.values()):
            # Convert to "C1 C2 ..."
            for chain, res_ids in fixed_res_dict.items():
                for rid in res_ids:
                    # RID is int.
                    fixed_list.append(f"{chain}{rid}")
            fixed_str = " ".join(fixed_list)

        print(f"[{i+1}/{len(design_files)}] [{category}] [{ion}] [{motif_name}] Processing {os.path.basename(pdb)}")
        
        # Command
        # conda run -n ligandmpnn python ...
        # But we are calling from python subprocess, better to use 'conda run' or direct python path
        # User said "The pipeline will not be able to spin up ligandmpnn instead of protein mpnn"
        # Implies we should assume environment availability or use explicit conda run.
        # User created 'ligandmpnn' env.
        
        # Using sys.executable for the WRAPPER script, but the SUBPROCESS needs 'ligandmpnn' python.
        # We can try finding it or just use 'conda run -n ligandmpnn python'.
        
        # Convert to absolute paths because we will change CWD to LigandMPNN dir
        pdb_abs = os.path.abspath(pdb)
        out_abs = os.path.abspath(mpnn_out)
        
        # Command construction
        # Optimization: Check if we are already in the ligandmpnn environment
        current_env = os.environ.get('CONDA_DEFAULT_ENV')
        
        if current_env == 'ligandmpnn':
             # Use current python executable directly
             base_cmd = [sys.executable, "run.py"]
        else:
             # Use conda run
             base_cmd = ["conda", "run", "-n", "ligandmpnn", "python", "run.py"]

        cmd = base_cmd + [
            "--pdb_path", pdb_abs,
            "--out_folder", out_abs,
            "--model_type", "ligand_mpnn",
            "--chains_to_design", design_chains_str,
            "--number_of_batches", str(args.num_seqs),
            "--batch_size", "1",
            "--seed", "37",
            "--pack_side_chains", "1",
            "--number_of_packs_per_design", "1"
        ]
        
        if fixed_str:
            cmd.extend(["--fixed_residues", fixed_str])
            
        if args.dry_run:
            print(f"[Dry Run] (CWD={args.lmpnn_path}) {' '.join(cmd)}")
            continue

        try:
            # Run with cwd set to LigandMPNN path so it can find model_params/
            subprocess.run(cmd, check=True, cwd=args.lmpnn_path, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
            pass
        except subprocess.CalledProcessError as e:
            print(f"  LigandMPNN failed for {pdb}")
            print(f"  Error: {e.stderr}")

if __name__ == "__main__":
    main()
