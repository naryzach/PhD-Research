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
        mpnn_out = os.path.join(parent_dir, "proteinmpnn")
        os.makedirs(mpnn_out, exist_ok=True)
        
        # Check if already run
        basename = os.path.basename(pdb).replace(".pdb", "")
        expected_output = os.path.join(mpnn_out, "seqs", f"{basename}.fa")
        if os.path.exists(expected_output) and not args.overwrite:
            print(f"[{i+1}/{len(design_files)}] Skipping {basename} (Already exists)")
            continue
        
        # Determine chains to design
        seqs_dict = get_sequence(pdb)
        if seqs_dict is None: seqs_dict = {}
        all_chains = sorted(seqs_dict.keys())
        
        protein_chains = [c for c in all_chains if len(seqs_dict[c]) > 5]
        ligand_chains = [c for c in all_chains if len(seqs_dict[c]) <= 5]
        
        design_chains = protein_chains
        design_chains_str = " ".join(design_chains)
        
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
        
        cmd = [
            sys.executable, os.path.join(args.pmpnn_path, 'protein_mpnn_run.py'),
            "--pdb_path", run_pdb,
            "--out_folder", mpnn_out,
            "--num_seq_per_target", str(args.num_seqs),
            "--batch_size", "1",
            "--suppress_print", "1",
            "--pdb_path_chains", design_chains_str
        ]
        
        print(f"[{i+1}/{len(design_files)}] Processing {os.path.basename(pdb)}")
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print(f"  MPNN failed for {pdb}")
        finally:
            if temp_pdb and os.path.exists(temp_pdb):
                os.remove(temp_pdb)

if __name__ == "__main__":
    main()
