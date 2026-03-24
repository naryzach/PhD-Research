import os
import json
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm
from biotite.structure.io.pdbx import CIFFile, get_structure as get_cif_structure
from atomworks.io.utils.io_utils import to_cif_file
from rf3.inference_engines.rf3 import RF3InferenceEngine
from rf3.utils.inference import InferenceInput

# --- CONFIGURATION ---
OUT_BASE_DIR = "../Local/TIMP-Dashboard_output"
ADVANCED_METRICS_PATH = os.path.join(OUT_BASE_DIR, "advanced_metrics.csv")
CROSS_DOCK_OUT_DIR = os.path.join(OUT_BASE_DIR, "cross_docking")
DATA_DIR = "../Data/TIMP_Complexes/AlphaFold_CIF"

TARGET_PDBS = {
    "ADAM10": "TIMP3_vs_ADAM10_AF.cif",
    "ADAM17": "TIMP3_vs_ADAM17_AF.cif",
    # "MMP2":   "TIMP3_vs_MMP2_AF.cif",
    # "MMP3":   "TIMP3_vs_MMP3_AF.cif",
    # "MMP9":   "TIMP3_vs_MMP9_AF.cif",
    # "MMP10":  "TIMP3_vs_MMP10_AF.cif"
}

TOP_K = 1 # Number of top designs PER Target-Combo to cross-dock

import biotite.structure.io.pdbx as pdbx

def load_target_structure(target_name):
    pdb_path = os.path.join(DATA_DIR, TARGET_PDBS[target_name])
    if not os.path.exists(pdb_path):
        return None
    # Use verified parsing method
    return pdbx.get_structure(pdbx.CIFFile.read(pdb_path), model=1)

def run_cross_docking():
    os.makedirs(CROSS_DOCK_OUT_DIR, exist_ok=True)
    
    # Setup Foundry Checkpoints
    checkpoint_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Tools", "foundry_checkpoints"))
    os.environ["FOUNDRY_CHECKPOINT_DIRS"] = checkpoint_dir
    
    if not os.path.exists(ADVANCED_METRICS_PATH):
        print(f"Metrics not found at {ADVANCED_METRICS_PATH}")
        return

    df = pd.read_csv(ADVANCED_METRICS_PATH)
    
    # Selection: find top designs per (target, loop_combo)
    top_designs = []
    for (target, combo), group in df.groupby(['target', 'loop_combo']):
        # Standardize target names for comparison
        clean_target = target.replace("TIMP3_vs_", "").replace("_AF", "").upper()
        top = group.sort_values("probability_of_binding_score", ascending=False).head(TOP_K)
        for _, row in top.iterrows():
            top_designs.append({
                "orig_row": row.to_dict(),
                "clean_target": clean_target,
                "full_seq": row['full_seq'],
                "design_id": row['design_id']
            })

    print(f"Selected {len(top_designs)} designs for cross-docking against {len(TARGET_PDBS)} targets.")
    
    # Prepare inference tasks
    # We need to fold [TIMP_seq_from_design] + [Target_B_structure]
    # To simplify, we use the original fold as a template but swap the target.
    # Actually, the most robust way is to just provide the sequence to RF3.
    
    rf3_engine = RF3InferenceEngine(ckpt_path='rf3', verbose=False)
    cross_dock_results = []

    for design in tqdm(top_designs, desc="Designs"):
        orig_id = design['design_id']
        orig_target = design['clean_target']
        timp_seq = design['full_seq']
        
        # We also need to know the length of TIMP3 to separate it from the target in the complex
        # TIMP3 ends at the last residue before Chain B starts.
        # In our pipeline, Chain A is TIMP, Chain B is Target.
        
        for target_name in TARGET_PDBS.keys():
            # If target_name == orig_target, we already have the result in advanced_metrics.csv
            # But for consistency, we could re-run or just copy.
            
            new_id = f"{orig_id}_vs_{target_name}_XD"
            out_folder = os.path.join(CROSS_DOCK_OUT_DIR, target_name)
            os.makedirs(out_folder, exist_ok=True)
            
            refolded_path = os.path.join(out_folder, f"{new_id}_refolded.cif")
            if os.path.exists(refolded_path):
                # Optionally skip if already done
                continue

            # Load target structure as template for the refolding
            target_struct = load_target_structure(target_name)
            if target_struct is None: continue
            
            # Identify Chain B residues
            target_b = target_struct[target_struct.chain_id == "B"]
            
            # We need an InferenceInput. 
            # We'll construct it using the sequence of TIMP and the structure of Target B.
            # However, RF3 often takes a single atom array and refolds it.
            # We'll create a dummy structure with the new TIMP sequence (as CA-only) and the target B.
            
            # Simpler: just use sequence-only input if RF3 supports it, 
            # or use the original multi-chain structure and swap the TIMP sequence.
            # Re-using the design structure is easy.
            
            # Let's find the original CIF
            orig_row = design['orig_row']
            orig_cif_path = os.path.join(OUT_BASE_DIR, orig_row['target'], orig_row['loop_combo'], "rf3", f"{orig_id}_refolded.cif")
            
            if not os.path.exists(orig_cif_path):
                print(f"Warning: Original CIF not found at {orig_cif_path}")
                continue
                
            orig_complex = pdbx.get_structure(pdbx.CIFFile.read(orig_cif_path), model=1)
            
            # Create a cross-docked complex: 
            # 1. Take TIMP from original complex (already has the design loops)
            # 2. Take Target from the NEW target_struct
            timp_part = orig_complex[orig_complex.chain_id == "A"]
            target_part = target_struct[target_struct.chain_id == "B"]
            
            # Combine them. Note: residue numbering might overlap, but RF3 handles chains.
            # We might need to ensure they don't clash too horribly initially.
            xd_complex = timp_part + target_part
            
            try:
                input_spec = InferenceInput.from_atom_array(xd_complex, example_id=new_id)
                rf3_outputs = rf3_engine.run(inputs=input_spec)
                
                if new_id in rf3_outputs:
                    res = rf3_outputs[new_id][0]
                    to_cif_file(res.atom_array, refolded_path, file_type="cif")
                    
                    # Save summary
                    sum_path = os.path.join(out_folder, f"{new_id}_summary.json")
                    with open(sum_path, 'w') as f:
                        json.dump(res.summary_confidences, f)
                    
                    # Store result record
                    rec = {
                        "design_id": orig_id,
                        "intended_target": orig_target,
                        "folded_target": target_name,
                        "plddt": res.summary_confidences.get('overall_plddt', 0),
                        "ptm": res.summary_confidences.get('ptm', 0),
                        "iptm": res.summary_confidences.get('iptm', 0)
                    }
                    cross_dock_results.append(rec)
            except Exception as e:
                print(f"Error folding {new_id}: {e}")
            finally:
                torch.cuda.empty_cache()

    # Save results
    if cross_dock_results:
        xd_df = pd.DataFrame(cross_dock_results)
        xd_df.to_csv(os.path.join(OUT_BASE_DIR, "cross_docking_metrics.csv"), index=False)
        print(f"Cross-docking complete. Saved to {os.path.join(OUT_BASE_DIR, 'cross_docking_metrics.csv')}")

if __name__ == "__main__":
    run_cross_docking()
