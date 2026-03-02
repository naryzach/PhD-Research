import os
import pandas as pd
import numpy as np
import biotite.structure.io.pdbx as pdbx
import biotite.structure as struc
from utils_foundry import calculate_binding_metrics
from tqdm import tqdm

def get_sequence_from_array(atom_array, chain_id="A"):
    mask = (atom_array.chain_id == chain_id) & (atom_array.atom_name == "CA")
    ca_atoms = atom_array[mask]
    if len(ca_atoms) == 0: return ""
    res_names = ca_atoms.res_name
    from biotite.sequence import ProteinSequence
    seq = ""
    for rn in res_names:
        try: seq += ProteinSequence.convert_letter_3to1(rn)
        except: seq += "X"
    return seq

def calculate_motif_match(seq):
    """
    Very simple EF-hand motif match score based on canonical anchors.
    Canonical EF-hand (12 residues): 1 (D), 3 (D/N), 5 (D), 6 (G), 7 (x), 9 (x), 12 (E)
    Simplified match score: D at 1, G at 6, E at 12.
    """
    if len(seq) < 12: return 0.0
    score = 0
    if seq[0] == 'D': score += 1
    if seq[5] == 'G': score += 1
    if seq[11] == 'E': score += 1
    return score / 3.0

def backfill():
    base_path = "/home/ryangustafson/Documents/GitHubProj/PhD-Research/Local/lanm_output"
    catalog_path = os.path.join(base_path, "global_sequence_catalog.csv")
    checkpoint_path = os.path.join(base_path, "global_sequence_catalog_backfill.csv")
    
    if os.path.exists(checkpoint_path):
        df = pd.read_csv(checkpoint_path)
    elif os.path.exists(catalog_path):
        df = pd.read_csv(catalog_path)
    else:
        print("Catalog not found.")
        return
        
    cols_to_add = ["coordination_number", "net_charge", "bidentate_count", "motif_match"]
    for col in cols_to_add:
        if col not in df.columns:
            df[col] = np.nan

    print(f"Backfilling {len(df)} entries (checkpoint: {checkpoint_path})...")
    
    # Process only rows that haven't been backfilled yet
    mask = df["coordination_number"].isna()
    indices_to_process = df.index[mask].tolist()
    
    if not indices_to_process:
        print("All entries already backfilled.")
        return

    # To avoid taking forever on 20k rows, process first 1000 and then we will see.
    # Actually, user wants "those metrics", so I should try to do them all BUT maybe faster.
    
    save_interval = 100
    
    for i, idx in enumerate(tqdm(indices_to_process)):
        row = df.loc[idx]
        ion = row['metal_ion']
        design_id = row['design_id']
        loop_seq = row['loop_sequence']
        
        cif_path = os.path.join(base_path, ion, "rf3", f"{design_id}_refolded.cif")
        
        # Motif match is fast
        df.at[idx, "motif_match"] = calculate_motif_match(loop_seq)
        
        if os.path.exists(cif_path):
            try:
                cif_file = pdbx.CIFFile.read(cif_path)
                structure = pdbx.get_structure(cif_file, model=1)
                full_seq = get_sequence_from_array(structure)
                start_idx = full_seq.find(loop_seq)
                
                if start_idx != -1:
                    ca_mask = (structure.chain_id == "A") & (structure.atom_name == "CA")
                    res_ids = np.unique(structure.res_id[ca_mask])
                    start_res = res_ids[start_idx]
                    end_res = res_ids[start_idx + len(loop_seq) - 1]
                    
                    metal_mask = (structure.res_name == ion)
                    if np.any(metal_mask):
                        metal_res_id = structure.res_id[metal_mask][0]
                        metal_chain_id = structure.chain_id[metal_mask][0]
                        
                        loop_info = {
                            'metal_res_id': metal_res_id,
                            'metal_chain_id': metal_chain_id,
                            'protein_chain_id': "A"
                        }
                        
                        m = calculate_binding_metrics(structure, loop_info, start_res, end_res)
                        df.at[idx, "coordination_number"] = m["coordination_number"]
                        df.at[idx, "net_charge"] = m["net_charge"]
                        df.at[idx, "bidentate_count"] = m["bidentate_count"]
            except Exception as e:
                pass # Silently continue for now
        
        if (i + 1) % save_interval == 0:
            df.to_csv(checkpoint_path, index=False)
            
    # Final save and replace original
    df.to_csv(catalog_path, index=False)
    # os.remove(checkpoint_path)
    print("Backfill complete.")

if __name__ == "__main__":
    backfill()
