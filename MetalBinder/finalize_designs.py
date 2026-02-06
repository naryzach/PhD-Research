import argparse, os, sys, subprocess, glob
import pandas as pd
import json
import re
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
    
    # Check for HETATM metals if Z not found (or always check)
    # If we restored the metal, it will be a HETATM
    for model in structure:
        for chain in model:
            if chain.id == 'Z' and 'Z' not in sequences:
                # Check if it has atoms
                has_atoms = any(True for _ in chain.get_atoms())
                if has_atoms:
                    # It's likely the metal. We use a placeholder for sequence
                    sequences['Z'] = 'X' # X marks the spot
        break

    return sequences
        
    return sequences

def main():
    parser = argparse.ArgumentParser(description="Finalize Designs and Create AlphaFold CSV")
    parser.add_argument("--pred_dir", default="../Local/Metal_Predictions", help="Predictions Directory")
    parser.add_argument("--out_csv", default="../Local/Metal_Predictions/alphafold_inputs.csv", help="Output CSV path")
    parser.add_argument("--out_json", default="../Local/Metal_Predictions/alphafold_inputs.json", help="Output JSON path for AlphaFold")
    parser.add_argument("--run_mpnn", action="store_true", help="Run ProteinMPNN on designs")
    parser.add_argument("--pmpnn_path", default="../Tools/ProteinMPNN", help="Path to ProteinMPNN")
    args = parser.parse_args()
    
    # recursive search
    print(f"Scanning {args.pred_dir}...")
    pdb_files = glob.glob(os.path.join(args.pred_dir, "**", "*.pdb"), recursive=True)
    
    # Filter for RFdiffusion outputs
    design_files = [f for f in pdb_files if "rfdiffusion" in f and ("scaffold" in os.path.basename(f) or "design" in os.path.basename(f)) and "traj" not in os.path.basename(f)]
    
    if args.run_mpnn:
        print(f"Running ProteinMPNN on {len(design_files)} designs...")
        project_root = os.path.dirname(os.path.abspath(__file__))
        
        for i, pdb in enumerate(design_files):
            # Output folder: sibling of rfdiffusion
            # .../DeNovo/ZN/Motif/rfdiffusion/scaffold.pdb
            # -> .../DeNovo/ZN/Motif/proteinmpnn/
            
            dirname = os.path.dirname(pdb)
            parent_dir = os.path.dirname(dirname) 
            mpnn_out = os.path.join(parent_dir, "proteinmpnn")
            os.makedirs(mpnn_out, exist_ok=True)
            
            # Check if already run? (Skipping check for now to ensure run)
            
            # Determine chains to design (exclude Z)
            seqs_dict = get_sequence(pdb)
            if seqs_dict is None: seqs_dict = {}
            all_chains = sorted(seqs_dict.keys())
            design_chains = [c for c in all_chains if c != 'Z']
            design_chains_str = " ".join(design_chains)
            
            # Check if we need to masquerade HETATM Z as GLY Z for MPNN?
            # If Z is 'X' (from our get_sequence update), it means it's a HETATM.
            # ProteinMPNN will ignore HETATM. To make MPNN respect steric constraints, 
            # we must convert it to a GLY anchor temporarily.
            
            run_pdb = pdb
            temp_pdb = None
            
            if 'Z' in seqs_dict and seqs_dict['Z'] == 'X':
                # Create temporary PDB with GLY Z
                temp_pdb = pdb.replace(".pdb", "_mpnn_input.pdb")
                with open(pdb, 'r') as f:
                    lines = f.readlines()
                
                new_lines = []
                for line in lines:
                    if line.startswith("HETATM") and " Z " in line:
                        # HETATM 9999 ZN   ZN Z   1 ... -> ATOM ... CA  GLY Z   1 ...
                        # Parse coords
                        try:
                            # Fixed column format? PDB is fixed width.
                            # x: 30-38, y: 38-46, z: 46-54
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            
                            gly_line = f"ATOM  {9999:5d}  CA  GLY Z   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
                            new_lines.append(gly_line)
                        except:
                            new_lines.append(line) # Fallback
                    else:
                        new_lines.append(line)
                
                with open(temp_pdb, 'w') as f:
                    f.writelines(new_lines)
                
                run_pdb = temp_pdb
            
            # ProteinMPNN run command
            cmd = [
                sys.executable, os.path.join(args.pmpnn_path, 'protein_mpnn_run.py'),
                "--pdb_path", run_pdb,
                "--out_folder", mpnn_out,
                "--num_seq_per_target", "5",
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


    # Now collect sequences
    # Prefer MPNN outputs
    data = []
    
    # Re-scan for MPNN outputs if run_mpnn was used or just check existence
    # Logic: Iterate RFdiffusion designs, look for corresponding MPNN designs.
    
    # Logic: Iterate RFdiffusion designs, look for corresponding MPNN designs.
    
    json_data = []

    for pdb in design_files:
        dirname = os.path.dirname(pdb)
        parent_dir = os.path.dirname(dirname)
        mpnn_out = os.path.join(parent_dir, "proteinmpnn")
        
        # Original Name: scaffold_0.pdb
        basename = os.path.basename(pdb).replace(".pdb", "")
        
        # ProteinMPNN outputs: seqs/{basename}.fa
        
        found_seqs = False
        
        seqs_dir = os.path.join(mpnn_out, "seqs")
        fasta_file = os.path.join(seqs_dir, f"{basename}.fa")

        if os.path.exists(fasta_file):
            try:
                with open(fasta_file, 'r') as f:
                    lines = f.readlines()
                
                # Parse FASTA and find best score
                # >T=0.1, sample=1, score=1.2345, global_score=..., seq_recovery=...
                # SEQUENCE
                
                best_score = float('inf')
                best_seq = None
                best_sample_id = None
                
                current_header = ""
                for line in lines:
                    line = line.strip()
                    if line.startswith(">"):
                        current_header = line
                    else:
                        if not line: continue
                        seq = line
                        
                        # Parse score
                        # Example header: >T=0.1, sample=1, score=1.2345, global_score=...
                        score_match = re.search(r"score=([-\d\.]+)", current_header)
                        if score_match:
                            score = float(score_match.group(1))
                            if score < best_score:
                                best_score = score
                                best_seq = seq
                                
                                # Extract sample ID for uniqueness
                                sample_match = re.search(r"sample=(\d+)", current_header)
                                if sample_match:
                                     best_sample_id = sample_match.group(1)
                                else:
                                     best_sample_id = "0"

                if best_seq:
                    # Create IDs
                    rel_path = os.path.relpath(pdb, args.pred_dir)
                    base_id = rel_path.replace("/", "_").replace(".pdb", "")
                    unique_id = f"{base_id}_pmpnn_{best_sample_id}"

                    
                    # Moved CSV append to after JSON logic to use full sequence
                    
                    
                    # Prepare JSON entry
                    chains_seqs = best_seq.split('/')
                    
                    # Re-derive chain info
                    seqs_dict = get_sequence(pdb)
                    if seqs_dict:
                        all_chains = sorted(seqs_dict.keys())
                        # Must match the logic used for --pdb_path_chains
                        design_chains = [c for c in all_chains if c != 'Z']
                    else:
                        print(f"Error parse chains for {pdb}")
                        continue

                    # Mapping: design_chains[i] -> chains_seqs[i]
                    # We iterate all_chains to preserve order in JSON (likely A, Z)
                    # And populate appropriate data
                    
                    if len(chains_seqs) != len(design_chains):
                         print(f"Warning: Design seq count mismatch for {unique_id}: {len(chains_seqs)} vs {len(design_chains)}")
                    
                    # Create a map for designed sequences
                    try:
                        design_map = {chain: seq for chain, seq in zip(design_chains, chains_seqs)}
                    except:
                        design_map = {} # Should not happen if check passed

                    try:
                        metal_type = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(pdb))))
                        if not (len(metal_type) > 0 and "site" not in metal_type):
                             metal_type = "ZN" 
                    except:
                        metal_type = "ZN"

                    json_sequences = []
                    
                    for chain_id in all_chains:
                        if chain_id == 'Z':
                             # Ligand Chain
                             json_sequences.append({
                                "ligand": {
                                    "ccdCodes": [metal_type],
                                    "count": 1
                                }
                            })
                        elif chain_id in design_map:
                             # Designed Protein Chain
                             json_sequences.append({
                                "proteinChain": {
                                    "sequence": design_map[chain_id],
                                    "count": 1
                                }
                            })
                        else:
                             # Fixed Protein Chain (original sequence)
                             json_sequences.append({
                                "proteinChain": {
                                    "sequence": seqs_dict.get(chain_id, ""),
                                    "count": 1
                                }
                            })
                    
                    json_data.append({
                        "name": unique_id,
                        "sequences": json_sequences
                    })

                    # Construct full sequence string for CSV (Chain A / Chain Z)
                    # For Ligand, we use 'G' (as in PDB) or we could use the metal code? 
                    # Users usually expect amino acids in sequence column. 'G' is safer.
                    # We will add a separate 'ligand' column.
                    
                    csv_seq_parts = []
                    for item in json_sequences:
                        if "proteinChain" in item:
                            csv_seq_parts.append(item["proteinChain"]["sequence"])
                        elif "ligand" in item:
                            csv_seq_parts.append("G") # Ligand is G in PDB
                    
                    full_seq_str = "/".join(csv_seq_parts)
                    
                    data.append({
                        "id": unique_id,
                        "sequence": full_seq_str,
                        "ligand": metal_type if metal_type else "",
                        "path": os.path.abspath(pdb), 
                        "source": "ProteinMPNN"
                    })


                    found_seqs = True
            except Exception as e:
                print(f"Error parsing {fasta_file}: {e}")
                    
        if not found_seqs:
            # Fallback to RFdiffusion backbone (Poly-Gly?)
            seqs = get_sequence(pdb)
            if seqs:
                full_seq = ":".join(seqs.values())
                rel_path = os.path.relpath(pdb, args.pred_dir)
                unique_id = rel_path.replace("/", "_").replace(".pdb", "")
                
                data.append({
                    "id": unique_id,
                    "sequence": full_seq, 
                    "path": os.path.abspath(pdb),
                    "source": "RFdiffusion_Backbone"
                })
                
                # Also add backbone to JSON? Maybe not if it's Poly-Gly/Backbone only. 
                # Converting backbone dict to list
                json_sequences = []
                for chain_id, chain_seq in seqs.items():
                     json_sequences.append({
                        "proteinChain": {
                            "sequence": chain_seq,
                            "count": 1
                        }
                    })
                
                json_data.append({
                    "name": unique_id,
                    "sequences": json_sequences
                })

    if data:
        df = pd.DataFrame(data)
        df.to_csv(args.out_csv, index=False)
        print(f"Saved {len(df)} sequences to {args.out_csv}")
        
    if json_data:
        with open(args.out_json, 'w') as f:
            json.dump(json_data, f, indent=2)
        print(f"Saved {len(json_data)} JSON entries to {args.out_json}")
    else:
        print("No designs found.")

if __name__ == "__main__":
    main()
