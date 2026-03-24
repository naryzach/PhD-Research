import os
import json
import numpy as np
import pandas as pd
import re
from scipy.spatial.distance import cdist
from biotite.structure.io.pdbx import CIFFile, get_structure
from biotite.structure import sasa, hbond
from biotite.sequence import ProteinSequence

OUT_BASE_DIR = "../Local/TIMP-Dashboard_output"
CATALOG_PATH = os.path.join(OUT_BASE_DIR, "global_sequence_catalog.csv")

def get_sequence_from_array(atom_array, chain_id="A"):
    mask = (atom_array.chain_id == chain_id) & (atom_array.atom_name == "CA")
    ca_atoms = atom_array[mask]
    if len(ca_atoms) == 0: return ""
    sort_idx = np.argsort(ca_atoms.res_id)
    return "".join([ProteinSequence.convert_letter_3to1(rn) if rn != "UNK" else "X" for rn in ca_atoms.res_name[sort_idx]])

def compute_biotite_interactions(atom_array, chain_A="A", chain_B="B"):
    res = {"h_bonds": 0, "salt_bridges": 0, "hydrophobic_contacts": 0, "bsa": 0.0}
    try:
        # BSA
        mask_A = atom_array.chain_id == chain_A
        mask_B = atom_array.chain_id == chain_B
        
        sasa_all = sasa(atom_array)
        sasa_A_iso = sasa(atom_array[mask_A])
        sasa_B_iso = sasa(atom_array[mask_B])
        
        total_iso = np.sum(sasa_A_iso) + np.sum(sasa_B_iso)
        total_complex = np.sum(sasa_all[mask_A]) + np.sum(sasa_all[mask_B])
        res["bsa"] = float((total_iso - total_complex) / 2.0)
        
        # Polar Contacts (H-bond approximations)
        polar_res = ["N", "O", "S"]
        polar_A = atom_array[(atom_array.chain_id == chain_A) & np.isin(atom_array.element, polar_res)]
        polar_B = atom_array[(atom_array.chain_id == chain_B) & np.isin(atom_array.element, polar_res)]
        hbond_count = 0
        if len(polar_A) > 0 and len(polar_B) > 0:
            dists = cdist(polar_A.coord, polar_B.coord)
            hbond_count = int(np.sum(dists < 3.5))
        res["h_bonds"] = hbond_count

        # Salt bridges (basic N to acidic O)
        basic_res = ["ARG", "LYS", "HIS"]
        acidic_res = ["ASP", "GLU"]
        n_atoms = atom_array[np.isin(atom_array.res_name, basic_res) & (atom_array.element == "N")]
        o_atoms = atom_array[np.isin(atom_array.res_name, acidic_res) & (atom_array.element == "O")]
        if len(n_atoms) > 0 and len(o_atoms) > 0:
            for (c1, c2) in [(chain_A, chain_B), (chain_B, chain_A)]:
                n_chain = n_atoms[n_atoms.chain_id == c1]
                o_chain = o_atoms[o_atoms.chain_id == c2]
                if len(n_chain) > 0 and len(o_chain) > 0:
                    dists = cdist(n_chain.coord, o_chain.coord)
                    res["salt_bridges"] += int(np.sum(dists < 4.0))

        # Hydrophobic Contacts
        hydro_res = ["VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"]
        c_atoms_A = atom_array[(atom_array.chain_id == chain_A) & np.isin(atom_array.res_name, hydro_res) & (atom_array.element == "C")]
        c_atoms_B = atom_array[(atom_array.chain_id == chain_B) & np.isin(atom_array.res_name, hydro_res) & (atom_array.element == "C")]
        if len(c_atoms_A) > 0 and len(c_atoms_B) > 0:
            dists = cdist(c_atoms_A.coord, c_atoms_B.coord)
            res["hydrophobic_contacts"] += int(np.sum(dists < 5.0))

    except Exception as e:
        print(f"Error computing detailed interactions: {e}")
    
    return res

def analyze_all_results():
    if not os.path.exists(CATALOG_PATH):
        print(f"Catalog not found at {CATALOG_PATH}")
        return

    df = pd.read_csv(CATALOG_PATH)
    advanced_records = []

    for idx, row in df.iterrows():
        design_id = row['design_id']
        target_name = row['target']
        combo_name = row['loop_combo']
        
        base_dir = os.path.join(OUT_BASE_DIR, target_name, combo_name, "rf3")
        cif_path = os.path.join(base_dir, f"{design_id}_refolded.cif")
        confidences_path = os.path.join(base_dir, f"{design_id}_confidences.json")
        summary_path = os.path.join(base_dir, f"{design_id}_summary_confidences.json")
        
        # Base record with existing info
        rec = row.to_dict()
        # Initialize all advanced metrics for consistent CSV columns
        rec.update({
            "h_bonds": 0,
            "salt_bridges": 0,
            "hydrophobic_contacts": 0,
            "bsa": 0.0,
            "mean_loop_plddt": row.get("plddt", 0.0),
            "mean_loop_pae": 30.0,
            "TIMP_pLDDT": row.get("plddt", 0.0), 
            "Target_pLDDT": 0.0,
            "TIMP_pTM": 0.0, 
            "Target_pTM": 0.0,
            "probability_of_binding_score": 0.0
        })
        
        if not os.path.exists(cif_path):
            advanced_records.append(rec)
            continue
            
        try:
            # Biochemical features
            import biotite.structure.io.pdbx as pdbx
            # Use get_structure which handles the file reading internally or via a generic file
            atom_array = pdbx.get_structure(pdbx.CIFFile.read(cif_path), model=1)
            inter_metrics = compute_biotite_interactions(atom_array, "A", "B")
            rec.update(inter_metrics)
            
            # Confidence features
            if os.path.exists(confidences_path):
                with open(confidences_path, 'r') as f:
                    conf = json.load(f)
                
                chain_A_seq = get_sequence_from_array(atom_array, "A")
                chain_A_len = len(chain_A_seq)
                loop_indices = set()
                
                # Identify loop residues
                for loop in ["AB", "C", "EF", "GH", "Multi"]:
                    seq_key = f"loop_{loop}_seq"
                    if seq_key in row and pd.notna(row[seq_key]) and row[seq_key] != "MISSING":
                        # Ensure sequence key exists in array, or fallback
                        match = re.search(str(row[seq_key]), chain_A_seq)
                        if match:
                            loop_indices.update(range(match.start(), match.end()))
                
                # Fallback to entire TIMP3 chain if regex targeting failed
                if not loop_indices:
                    loop_indices = set(range(chain_A_len))
                
                loop_indices = list(loop_indices)
                # Handle RF3 specific keys: atom_plddts instead of plddt
                if "plddt" not in conf and "atom_plddts" in conf:
                    conf["plddt"] = conf["atom_plddts"]
                
                if loop_indices and "pae" in conf and "plddt" in conf:
                    pae = np.array(conf["pae"])
                    plddt = np.array(conf["plddt"])
                    
                    # If plddt is per-atom, we might need to map it back to residues
                    # or just take the first N (if it's already residue-level but named differently)
                    if len(plddt.shape) == 2: plddt = plddt[0] # Handle batch dim
                    if len(pae.shape) == 3: pae = pae[0] # Handle batch dim
                    
                    # If it's still weirdly shaped (e.g. per-atom), we'll do our best
                    # Usually token_res_ids implies plddt per residue.
                    
                    loop_plddts = [plddt[i] for i in loop_indices if i < len(plddt)]
                    rec["mean_loop_plddt"] = np.mean(loop_plddts) if loop_plddts else rec["mean_loop_plddt"]
                    
                    # Chain-specific pLDDT and pTM
                    if len(plddt) >= chain_A_len:
                        rec["TIMP_pLDDT"] = float(np.mean(plddt[:chain_A_len]))
                        if len(plddt) > chain_A_len:
                            rec["Target_pLDDT"] = float(np.mean(plddt[chain_A_len:]))
                        else:
                            # DEBUG: If they match, why?
                            print(f"[DEBUG] {design_id}: plddt_len={len(plddt)}, chain_A_len={chain_A_len}")
                    
                    # Approximate per-chain pTM from PAE if available
                    if "ptm" in conf:
                        rec["ptm"] = conf["ptm"]
                        if isinstance(conf["ptm"], list) and len(conf["ptm"]) >= 2:
                            rec["TIMP_pTM"] = conf["ptm"][0]
                            rec["Target_pTM"] = conf["ptm"][1]
                    
                    # Fallback chain pTM from PAE diagonal blocks
                    if rec["TIMP_pTM"] == 0.0 and "pae" in conf:
                        # Convert mean PAE to a pTM-like score (0-1, higher is better)
                        # Approximation: exp(-PAE/10)
                        timp_pae = float(np.mean(pae[:chain_A_len, :chain_A_len]))
                        rec["TIMP_pTM"] = round(float(np.exp(-timp_pae / 10.0)), 3)
                        
                        if pae.shape[0] > chain_A_len:
                            tgt_pae = float(np.mean(pae[chain_A_len:, chain_A_len:]))
                            rec["Target_pTM"] = round(float(np.exp(-tgt_pae / 10.0)), 3)
                    
                    # Loop PAE vs Target (Chain B)
                    # Loop indices x Target Indices (chain A len to end)
                    if pae.shape[0] > chain_A_len and pae.shape[1] > chain_A_len:
                        loop_target_paes = []
                        for i in loop_indices:
                            if i < pae.shape[0]:
                                loop_target_paes.extend(pae[i, chain_A_len:])
                        rec["mean_loop_pae"] = float(np.mean(pae[loop_indices, chain_A_len:])) if len(loop_indices) > 0 else 30.0
        except Exception as e:
            print(f"Error processing {design_id}: {e}")
        
        # Calculate Probability of Binding Score (AlphaFoldServer_SA formula)
        iptm = rec.get("iptm", 0.0)
        rec["iptm"] = iptm
        
        hb_score = min(rec["h_bonds"], 15) / 15.0
        pae_score = max(0, (15.0 - rec["mean_loop_pae"]) / 15.0)
        bsa_score = min(rec["bsa"], 1000) / 1000.0
        
        comp_prob = (iptm * 0.35) + (bsa_score * 0.25) + (pae_score * 0.20) + (hb_score * 0.20)
        
        sb_bonus = min(rec["salt_bridges"] * 0.02, 0.10)
        hydro_bonus = min(rec["hydrophobic_contacts"] * 0.01, 0.05)
        
        comp_prob = (comp_prob + sb_bonus + hydro_bonus) * 100.0
        comp_prob = min(comp_prob, 100.0)
        
        rec["probability_of_binding_score"] = round(comp_prob, 2)
        
        advanced_records.append(rec)

    adv_df = pd.DataFrame(advanced_records)
    adv_df.to_csv(os.path.join(OUT_BASE_DIR, "advanced_metrics.csv"), index=False)
    print("Saved advanced_metrics.csv")
    
    # Generate Consensus Leaderboard
    # Group by loop combination and sequence to find the true winners
    leaderboard = []
    for combo, group in adv_df.groupby("loop_combo"):
        # We find top 3 by probability_of_binding_score
        top = group.sort_values("probability_of_binding_score", ascending=False).head(3)
        for _, t in top.iterrows():
            leaderboard.append(t)
            
    if leaderboard:
        lead_df = pd.DataFrame(leaderboard).sort_values("probability_of_binding_score", ascending=False)
        lead_df.to_csv(os.path.join(OUT_BASE_DIR, "consensus_ranking.csv"), index=False)
        print("Saved consensus_ranking.csv")

if __name__ == "__main__":
    analyze_all_results()
