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
    
    # Allowed ions for AlphaFold
    ALLOWED_IONS = {"MG", "ZN", "CL", "CA", "NA", "MN", "K", "FE", "CU", "CO"}
    
    # Check for Ligand Chains (chains not in sequences or empty sequences)
    for model in structure:
        for chain in model:
            if chain.id not in sequences or len(sequences[chain.id]) == 0:
                # Check for residues to identify ligand
                residues = list(chain.get_residues())
                if residues:
                    # Use the first residue
                    residue = residues[0]
                    resname = residue.resname.strip().upper()
                    
                    # Check if resname is allowed
                    if resname in ALLOWED_IONS:
                        sequences[chain.id] = resname
                    else:
                        # Try to find an allowed ion in atoms
                        found_ion = None
                        for atom in residue:
                            if atom.element.upper() in ALLOWED_IONS:
                                found_ion = atom.element.upper()
                                break
                        
                        if found_ion:
                            sequences[chain.id] = found_ion
                        else:
                            # Fallback to resname if not HOH (Water)
                            if resname != 'HOH':
                                sequences[chain.id] = resname
        break
    
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
    design_files = [f for f in pdb_files if "rfdiffusion" in f and ("scaffold" in os.path.basename(f) or "design" in os.path.basename(f)) and "traj" not in os.path.basename(f) and "unidealized" not in f]
    
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
            
            # Determine chains to design
            seqs_dict = get_sequence(pdb)
            if seqs_dict is None: seqs_dict = {}
            all_chains = sorted(seqs_dict.keys())
            
            protein_chains = [c for c in all_chains if len(seqs_dict[c]) > 5]
            ligand_chains = [c for c in all_chains if len(seqs_dict[c]) <= 5]
            
            design_chains = protein_chains
            design_chains_str = " ".join(design_chains)
            
            # Check if we need to masquerade HETATM Ligands as GLY for MPNN?
            # ProteinMPNN will ignore HETATM. To make MPNN respect steric constraints, 
            # we must convert it to a GLY anchor temporarily.
            
            run_pdb = pdb
            temp_pdb = None
            
            if ligand_chains:
                # Create temporary PDB with GLY masquerading for all ligand chains
                temp_pdb = pdb.replace(".pdb", "_mpnn_input.pdb")
                with open(pdb, 'r') as f:
                    lines = f.readlines()
                
                new_lines = []
                for line in lines:
                    if line.startswith("HETATM"):
                        # Check if chain is in ligand_chains
                        # HETATM 1234 FE   FE  A 100 ...
                        # Chain ID is typically col 21 (0-indexed -> 21) -> line[21]
                        chain_id = line[21]
                        if chain_id in ligand_chains:
                             # Masquerade
                             try:
                                 x = float(line[30:38])
                                 y = float(line[38:46])
                                 z = float(line[46:54])
                                 # Use correct chain ID
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
    # Separate lists for categories
    data_all = []
    data_scaffold = []
    data_swap = []
    data_graft = []
    
    json_all = []
    json_scaffold = []
    json_swap = []
    json_graft = []

    for pdb in design_files:
        dirname = os.path.dirname(pdb)
        parent_dir = os.path.dirname(dirname)
        mpnn_out = os.path.join(parent_dir, "proteinmpnn")
        
        # Determine Category
        category = "Other"
        if "DeNovo" in pdb:
            category = "Scaffold"
        elif "Swaps" in pdb:
            category = "Swap"
        elif "Grafted" in pdb:
            category = "Graft"
        
        # Original Name: scaffold_0.pdb
        basename = os.path.basename(pdb).replace(".pdb", "")
        
        # ProteinMPNN outputs: seqs/{basename}.fa
        
        found_seqs = False
        
        seqs_dir = os.path.join(mpnn_out, "seqs")
        fasta_file = os.path.join(seqs_dir, f"{basename}.fa")

        if not os.path.exists(fasta_file):
             # Try _mpnn_input suffix (used during ligand masquerading)
             fasta_file = os.path.join(seqs_dir, f"{basename}_mpnn_input.fa")

        entry_data = None
        entry_json = None

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
                        # Heuristic: Protein chains have length > 5 (residues)
                        # Ligand chains have length <= 5 (e.g. "ZN", "SF4")
                        protein_chains = [c for c in all_chains if len(seqs_dict[c]) > 5]
                        ligand_chains = [c for c in all_chains if len(seqs_dict[c]) <= 5]
                        
                        # We assume we design all protein chains?
                        # Or match --pdb_path_chains logic.
                        design_chains = protein_chains

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
                        if chain_id in ligand_chains:
                             # Ligand Chain
                             ligand_code = seqs_dict[chain_id]
                             json_sequences.append({
                                "ion": {
                                    "ion": ligand_code,
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
                    
                    entry_json = {
                        "name": unique_id,
                        "sequences": json_sequences
                    }

                    # Construct full sequence string for CSV (Chain A / Chain Z)
                    # For Ligand, we use 'G' (as in PDB) or we could use the metal code? 
                    # Users usually expect amino acids in sequence column. 'G' is safer.
                    # We will add a separate 'ligand' column.
                    
                    csv_seq_parts = []
                    for item in json_sequences:
                        if "proteinChain" in item:
                            csv_seq_parts.append(item["proteinChain"]["sequence"])
                        # elif "ion" in item:
                        #     csv_seq_parts.append("G") # Ligand is G in PDB - Removed per user request
                    
                    full_seq_str = "/".join(csv_seq_parts)
                    
                    # Extract ligand from PDB content instead of directory
                    detected_ligands = set()
                    for item in json_sequences:
                        if "ion" in item:
                            detected_ligands.add(item["ion"]["ion"])
                    
                    if detected_ligands:
                        final_ligand = "+".join(sorted(detected_ligands))
                    else:
                        # Fallback to directory inference if no ligand chain found (e.g. scaffold only?)
                        final_ligand = metal_type if metal_type else ""
                        
                    entry_data = {
                        "id": unique_id,
                        "sequence": full_seq_str,
                        "ligand": final_ligand,
                        "path": os.path.abspath(pdb), 
                        "source": "ProteinMPNN"
                    }

                    found_seqs = True
            except Exception as e:
                print(f"Error parsing {fasta_file}: {e}")
                    
        if not found_seqs:
            # Fallback to RFdiffusion backbone
            seqs = get_sequence(pdb)
            if seqs:
                all_backbone_chains = sorted(seqs.keys())
                prot_chains_bb = [c for c in all_backbone_chains if len(seqs[c]) > 5]
                lig_chains_bb = [c for c in all_backbone_chains if len(seqs[c]) <= 5]
                
                # Sequence: Only protein chains
                prot_seqs = [seqs[c] for c in prot_chains_bb]
                full_seq = "/".join(prot_seqs) # Use / separator for consistency
                
                # Ligand: From detected ions
                ligands = [seqs[c] for c in lig_chains_bb]
                final_ligand = "+".join(sorted(set(ligands)))
                
                rel_path = os.path.relpath(pdb, args.pred_dir)
                unique_id = rel_path.replace("/", "_").replace(".pdb", "")
                
                entry_data = {
                    "id": unique_id,
                    "sequence": full_seq, 
                    "ligand": final_ligand,
                    "path": os.path.abspath(pdb),
                    "source": "RFdiffusion_Backbone"
                }
                
                json_sequences = []
                for chain_id in all_backbone_chains:
                    if chain_id in lig_chains_bb:
                         json_sequences.append({
                            "ion": {
                                "ion": seqs[chain_id],
                                "count": 1
                            }
                        })
                    else:
                         json_sequences.append({
                            "proteinChain": {
                                "sequence": seqs[chain_id],
                                "count": 1
                            }
                        })
                
                entry_json = {
                    "name": unique_id,
                    "sequences": json_sequences
                }
        
        # Append to master lists and category lists
        if entry_data:
            data_all.append(entry_data)
            if category == "Scaffold":
                data_scaffold.append(entry_data)
            elif category == "Swap":
                data_swap.append(entry_data)
            elif category == "Graft":
                data_graft.append(entry_data)

        if entry_json:
             json_all.append(entry_json)
             if category == "Scaffold":
                json_scaffold.append(entry_json)
             elif category == "Swap":
                json_swap.append(entry_json)
             elif category == "Graft":
                json_graft.append(entry_json)

    # Function to save CSV
    def save_csv(data_list, filename):
        if data_list:
            df = pd.DataFrame(data_list)
            df.to_csv(filename, index=False)
            print(f"Saved {len(df)} sequences to {filename}")
        else:
            print(f"No designs found for {filename}")

    # Function to save JSON
    def save_json(data_list, filename):
        if data_list:
            with open(filename, 'w') as f:
                json.dump(data_list, f, indent=2)
            print(f"Saved {len(data_list)} JSON entries to {filename}")
        else:
            print(f"No designs found for {filename}")

    # Save All
    save_csv(data_all, args.out_csv)
    save_json(json_all, args.out_json)

    # Save Categories
    base_csv = args.out_csv.replace(".csv", "")
    base_json = args.out_json.replace(".json", "")

    save_csv(data_scaffold, f"{base_csv}_scaffold.csv")
    save_json(json_scaffold, f"{base_json}_scaffold.json")

    save_csv(data_swap, f"{base_csv}_swap.csv")
    save_json(json_swap, f"{base_json}_swap.json")

    save_csv(data_graft, f"{base_csv}_graft.csv")
    save_json(json_graft, f"{base_json}_graft.json")

if __name__ == "__main__":
    main()
