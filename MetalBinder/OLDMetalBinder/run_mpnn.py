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
    parser = argparse.ArgumentParser(description="Run ProteinMPNN on Design Folder")
    parser.add_argument("--pred_dir", default="../Local/Metal_Predictions", help="Predictions Directory")
    parser.add_argument("--pmpnn_path", default="../Tools/ProteinMPNN", help="Path to ProteinMPNN")
    parser.add_argument("--num_seqs", type=int, default=5, help="Number of sequences to generate per design")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")
    args = parser.parse_args()
    
    print(f"Scanning {args.pred_dir}...")
    pdb_files = glob.glob(os.path.join(args.pred_dir, "**", "*.pdb"), recursive=True)
    
    # Filter for RFdiffusion outputs
    design_files = [f for f in pdb_files if "rfdiffusion" in f and ("scaffold" in os.path.basename(f) or "design" in os.path.basename(f)) and "traj" not in os.path.basename(f) and "unidealized" not in f]
    
    print(f"Found {len(design_files)} designs for ProteinMPNN.")

    for i, pdb in enumerate(design_files):
        # Output folder: sibling of rfdiffusion
        dirname = os.path.dirname(pdb)
        parent_dir = os.path.dirname(dirname) 
        motif_name = os.path.basename(parent_dir)
        
        mpnn_out = os.path.join(parent_dir, "proteinmpnn")
        os.makedirs(mpnn_out, exist_ok=True)
        
        # Determine chains to design (needed for checking expected output name if ligand present)
        seqs_dict = get_sequence(pdb)
        if seqs_dict is None: seqs_dict = {}
        all_chains = sorted(seqs_dict.keys())
        
        protein_chains = [c for c in all_chains if len(seqs_dict[c]) > 5]
        ligand_chains = [c for c in all_chains if len(seqs_dict[c]) <= 5]
        
        design_chains = protein_chains
        design_chains_str = " ".join(design_chains)
        
        basename = os.path.basename(pdb).replace(".pdb", "")
        
        # If ligand chains exist, we masquerade, so output filename changes
        if ligand_chains:
             expected_basename = f"{basename}_mpnn_input"
        else:
             expected_basename = basename
             
        expected_output = os.path.join(mpnn_out, "seqs", f"{expected_basename}.fa")
        
        if os.path.exists(expected_output) and not args.overwrite:
            # Path parsing for logging
            ion_dir = os.path.dirname(parent_dir)
            ion = os.path.basename(ion_dir)
            category_dir = os.path.dirname(ion_dir)
            category = os.path.basename(category_dir)
            
            print(f"[{i+1}/{len(design_files)}] [{category}] [{ion}] [{motif_name}] Skipping {basename} (Already exists)")
            continue
        
        run_pdb = pdb
        temp_pdb = None
        
        if ligand_chains:
            # Create temporary PDB with GLY masquerading for all ligand chains
            temp_pdb = pdb.replace(".pdb", "_mpnn_input.pdb")
            try:
                with open(pdb, 'r') as f:
                    lines = f.readlines()
                
                new_lines = []
                for line in lines:
                    if line.startswith("HETATM"):
                        chain_id = line[21]
                        if chain_id in ligand_chains:
                             try:
                                 x = float(line[30:38])
                                 y = float(line[38:46])
                                 z = float(line[46:54])
                                 gly_line = f"ATOM  {9999:5d}  CA  GLY {chain_id}   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
                                 new_lines.append(gly_line)
                             except:
                                 new_lines.append(line)
                        else:
                             new_lines.append(line)
                    else:
                        new_lines.append(line)
                
                with open(temp_pdb, 'w') as f:
                    f.writelines(new_lines)
                
                run_pdb = temp_pdb
            except Exception as e:
                print(f"Error preparing ligand masquerade for {pdb}: {e}")
                continue
        
        # Calculate fixed positions (from ORIGINAL pdb, not masqueraded one, though coords should match)
        # Using run_pdb is safer if we masqueraded ligands, but we want to know what residue numbers to fix.
        # The masquerade only changes ligands to GLY, protein residues are untouched.
        # However, get_fixed_residues relies on finding METALs.
        # If run_pdb is temp_pdb, ligands might be GLY now?
        # Wait, if we converted metals to GLY, get_fixed_residues won't find them in temp_pdb!
        # FAST FIX: Compute fixed residues from ORIGINAL pdb.
        
        fixed_res_dict = get_fixed_residues(pdb)
        fixed_jsonl = None
        
        # Only create JSONL if we have fixed residues
        import json
        if any(v for v in fixed_res_dict.values()):
            # Format: {"pdb_name": "structure_name", "fixed_chain_id": [res1, res2]}?
            # No, standard MPNN format is: {"structure_name": {"chain_id": [res_ids]}} ?
            # Let's check help output or docs.
            # Help said: "dictionary with tied positions".
            # Actually --fixed_positions_jsonl "Path to dictionary with fixed positions"
            # Format usually: {"pdb_name": {"chain_A": [1, 2, 3], "chain_B": []}}
            
            # Since we pass --pdb_path, MPNN keys by the basename (name of the parsed structure)
            # When running single PDB, often the key is the pdb filename without extension?
            # Or we can use global fixed positions?
            
            # To be safe, let's look at protein_mpnn_run.py source or standard usage.
            # Standard usage usually expects the key to match the name of the structure.
            # If we run on single pdb, name is usually derived from filename.
            # run_pdb basename if temp: 'design_0..._mpnn_input.pdb' -> 'design_0..._mpnn_input'
            
            struct_name = os.path.basename(run_pdb).replace(".pdb", "")
            
            # Construct the dict
            # fixed_res_dict is { chain: [res_ids] }
            json_content = {struct_name: fixed_res_dict}
            
            fixed_jsonl = os.path.join(mpnn_out, f"{basename}_fixed.jsonl")
            with open(fixed_jsonl, 'w') as f:
                f.write(json.dumps(json_content))
        
        cmd = [
            sys.executable, os.path.join(args.pmpnn_path, 'protein_mpnn_run.py'),
            "--pdb_path", run_pdb,
            "--out_folder", mpnn_out,
            "--num_seq_per_target", str(args.num_seqs),
            "--batch_size", "1",
            "--suppress_print", "1",
            "--pdb_path_chains", design_chains_str
        ]
        
        if fixed_jsonl:
            cmd.extend(["--fixed_positions_jsonl", fixed_jsonl])
        
        # Path parsing for logging
        ion_dir = os.path.dirname(parent_dir)
        ion = os.path.basename(ion_dir)
        category_dir = os.path.dirname(ion_dir)
        category = os.path.basename(category_dir)
        
        print(f"[{i+1}/{len(design_files)}] [{category}] [{ion}] [{motif_name}] Processing {os.path.basename(pdb)}")
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as e:
            print(f"  MPNN failed for {pdb}")
            print(f"  Error: {e.stderr}")
        finally:
            if temp_pdb and os.path.exists(temp_pdb):
                os.remove(temp_pdb)
            if fixed_jsonl and os.path.exists(fixed_jsonl):
                os.remove(fixed_jsonl)

if __name__ == "__main__":
    main()
