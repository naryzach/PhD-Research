import os
import sys
import numpy as np
import pandas as pd
import argparse
import dataclasses
import gc
import torch
import time
import re
from itertools import combinations
from tqdm import tqdm

from biotite.structure.io.pdbx import CIFFile, get_structure as get_cif_structure
from biotite.structure.io.pdb import PDBFile
from biotite.structure import AtomArray, get_residues, array
from atomworks.io.utils.io_utils import to_cif_file
from biotite.sequence import ProteinSequence
from biotite.structure import rmsd, superimpose, sasa
from atomworks.constants import PROTEIN_BACKBONE_ATOM_NAMES

# Foundry Inference imports
from lightning.fabric import seed_everything
from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine
from mpnn.inference_engines.mpnn import MPNNInferenceEngine
from rf3.inference_engines.rf3 import RF3InferenceEngine
from rf3.utils.inference import InferenceInput
from rfd3.inference.input_parsing import DesignInputSpecification

torch.set_float32_matmul_precision('medium')

import logging
logging.getLogger("transforms").setLevel(logging.ERROR)
logging.getLogger("atomworks.io").setLevel(logging.ERROR)
logging.getLogger("atomworks.ml").setLevel(logging.ERROR)
logging.getLogger("foundry").setLevel(logging.ERROR)

# --- USER CONFIGURABLE VARIABLES ---
NUM_RFD3_DESIGNS = 3
NUM_LMPNN_STRUCTURES = 10
LMPNN_SAMPLING_TEMP = 0.1

TARGET_PDBS = [
    "PDB_fold_timp3_v_adam10cd_wt_model_0.pdb",
    "PDB_fold_timp3_v_mmp10cd_wt_model_0.pdb",
    "PDB_fold_timp3_v_mmp3cd_wt_model_0.pdb",
    "PDB_fold_timp3_variant_adam17cd_wt_model_0.pdb",
    "PDB_fold_timp3_variant_mmp2cd_wt_model_0.pdb",
    "PDB_fold_timp3_variant_mmp9cd_wt_model_0.pdb"
]

LOOP_DEFINITIONS = {
    "AB": {"normal": 6, "max": 15, "pos": 30, "left": "LVK", "right": "LVY"},
    "C": {"normal": 6, "max": 15, "pos": 62, "left": "HTE", "right": "GLK"},
    "EF": {"normal": 4, "max": 10, "pos": 92, "left": "MYT", "right": "FVE"},
    "GH": {"normal": 10, "max": 20, "pos": 127, "left": "KSC", "right": "NEC"},
    "Multi": {"normal": 10, "max": 20, "pos": 143, "left": "LWT", "right": "YQS"}
}

DATA_DIR = "../Data"
OUT_BASE_DIR = "../Local/TIMP-Dashboard_output"

CHAIN_TO_DESIGN = "A" # TIMP3 chain
FIXED_CHAINS = ["B"]  # Target chain

# --- HELPER FUNCTIONS ---
def get_sequence_from_array(atom_array, chain_id="A"):
    mask = (atom_array.chain_id == chain_id) & (atom_array.atom_name == "CA")
    ca_atoms = atom_array[mask]
    if len(ca_atoms) == 0: return ""
    
    sort_idx = np.argsort(ca_atoms.res_id)
    res_names = ca_atoms.res_name[sort_idx]
    
    seq_letters = []
    for rn in res_names:
        try:
            seq_letters.append(ProteinSequence.convert_letter_3to1(rn))
        except Exception:
            seq_letters.append("X")
    return "".join(seq_letters)

def get_structure_length(file_path, chain_id):
    try:
        if file_path.endswith('.cif') or file_path.endswith('.mmcif'):
            cif_file = CIFFile.read(file_path)
            atom_array = get_cif_structure(cif_file, model=1)
        else:
            pdb_file = PDBFile.read(file_path)
            atom_array = pdb_file.get_structure()[0]
        return len(get_sequence_from_array(atom_array, chain_id))
    except Exception:
        count = 0
        if not (file_path.endswith('.cif') or file_path.endswith('.mmcif')):
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') and line.split()[4] == chain_id and line.split()[2] == "CA":
                        count += 1
        return count

def renumber_atom_array_residues(atom_array):
    new_res_ids = np.zeros(len(atom_array), dtype=int)
    for chain_id in np.unique(atom_array.chain_id):
        chain_mask = atom_array.chain_id == chain_id
        chain_res_ids = atom_array.res_id[chain_mask]
        unique_old_ids = []
        last_id = None
        for r_id in chain_res_ids:
            if r_id != last_id:
                unique_old_ids.append(r_id)
                last_id = r_id
        id_map = {old_id: new_id for new_id, old_id in enumerate(unique_old_ids, start=1)}
        new_res_ids[chain_mask] = [id_map[old_id] for old_id in chain_res_ids]
    atom_array.res_id = new_res_ids
    return atom_array

def calc_protein_protein_metrics(atom_array, chain_A, chain_B):
    mask_A = (atom_array.chain_id == chain_A) & (atom_array.element != "H")
    mask_B = (atom_array.chain_id == chain_B) & (atom_array.element != "H")
    coords_A = atom_array.coord[mask_A]
    coords_B = atom_array.coord[mask_B]
    
    metrics = {
        "contacts": 0,
        "clashes": 0,
        "centroid_distance": float('nan'),
        "interface_area": float('nan')
    }
    
    if len(coords_A) == 0 or len(coords_B) == 0:
        return metrics
        
    try:
        from scipy.spatial.distance import cdist
        dist_matrix = cdist(coords_A, coords_B)
        metrics["contacts"] = int(np.sum(dist_matrix < 5.0))
        metrics["clashes"] = int(np.sum(dist_matrix < 2.2))
    except Exception:
        pass

    centroid_A = np.mean(coords_A, axis=0)
    centroid_B = np.mean(coords_B, axis=0)
    metrics["centroid_distance"] = float(np.linalg.norm(centroid_A - centroid_B))
    
    try:
        sasa_all = sasa(atom_array)
        sasa_A_iso = sasa(atom_array[mask_A])
        sasa_B_iso = sasa(atom_array[mask_B])
        total_iso = np.sum(sasa_A_iso) + np.sum(sasa_B_iso)
        total_complex = np.sum(sasa_all[mask_A]) + np.sum(sasa_all[mask_B])
        metrics["interface_area"] = float((total_iso - total_complex) / 2.0)
    except Exception:
        pass
        
    return metrics

def calculate_heuristic_score(contacts, interface_area, clashes, centroid_distance, plddt_mean):
    interface_area = 0.0 if np.isnan(interface_area) else interface_area
    score = (contacts * 0.4 
             + interface_area * 0.002 
             - clashes * 3.0 
             - centroid_distance * 0.2)
    if not np.isnan(plddt_mean):
        score += plddt_mean * 0.5
    return score

# --- GENERATION PIPELINE ---
def main():
    seed_everything(42)
    checkpoint_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Tools", "foundry_checkpoints"))
    os.makedirs(checkpoint_dir, exist_ok=True)
    os.environ["FOUNDRY_CHECKPOINT_DIRS"] = checkpoint_dir
    os.environ["DGLBACKEND"] = "pytorch"

    # Generate all mathematical combinations of loops
    loop_keys = list(LOOP_DEFINITIONS.keys())
    all_loop_combos = []
    for r in range(1, len(loop_keys) + 1):
        all_loop_combos.extend(list(combinations(loop_keys, r)))

    print(f"Total loop combinations to test: {len(all_loop_combos)}")
    print(f"Total targets: {len(TARGET_PDBS)}")

    final_records = []

    for target_pdb in TARGET_PDBS:
        pdb_path = os.path.join(DATA_DIR, target_pdb)
        if not os.path.exists(pdb_path):
            print(f"WARNING: Skipping {target_pdb}, file does not exist at {pdb_path}.")
            continue
            
        target_name = target_pdb.replace(".pdb", "").replace(".cif", "")
        fix_chain_len = get_structure_length(pdb_path, FIXED_CHAINS[0])
        timp3_length = get_structure_length(pdb_path, CHAIN_TO_DESIGN)
        
        valid_combos = []
        for combo in all_loop_combos:
            valid = True
            for loop_name in combo:
                loop_def = LOOP_DEFINITIONS[loop_name]
                if loop_def["pos"] + loop_def["normal"] > timp3_length:
                    valid = False
                    break
            if valid:
                valid_combos.append(combo)

        for combo in valid_combos:
            combo_name = "_".join(combo)
            print(f"\n======================================")
            print(f"Processing Target: {target_name} | Loops: {combo_name}")
            print(f"======================================")

            combo_rfd3_out_dir = os.path.join(OUT_BASE_DIR, target_name, combo_name, "rfd3")
            combo_lmpnn_out_dir = os.path.join(OUT_BASE_DIR, target_name, combo_name, "lmpnn")
            combo_rf3_out_dir = os.path.join(OUT_BASE_DIR, target_name, combo_name, "rf3")
            os.makedirs(combo_rfd3_out_dir, exist_ok=True)
            os.makedirs(combo_lmpnn_out_dir, exist_ok=True)
            os.makedirs(combo_rf3_out_dir, exist_ok=True)

            selected_loops = [LOOP_DEFINITIONS[name] for name in combo]
            selected_loops.sort(key=lambda x: x["pos"])

            contig_parts = []
            current_pos = 1
            for loop in selected_loops:
                if current_pos <= loop["pos"]:
                    contig_parts.append(f"{CHAIN_TO_DESIGN}{current_pos}-{loop['pos']}")
                contig_parts.append(f"{loop['normal']}-{loop['max']}")
                current_pos = loop["pos"] + loop["normal"] + 1

            if current_pos <= timp3_length:
                contig_parts.append(f"{CHAIN_TO_DESIGN}{current_pos}-{timp3_length}")

            contig_string = ",".join(contig_parts)
            full_contig_string = f"{contig_string},/0,{FIXED_CHAINS[0]}1-{fix_chain_len}"

            min_inpaint_len = sum(loop["normal"] for loop in selected_loops)
            max_inpaint_len = sum(loop["max"] for loop in selected_loops)
            native_removed_len = sum(loop["normal"] for loop in selected_loops) 
                
            base_length = timp3_length + fix_chain_len - native_removed_len
            total_output_length_min = base_length + min_inpaint_len
            total_output_length_max = base_length + max_inpaint_len

            # --- STAGE 1: RUN RFd3 ---
            print("--- Running RFd3 ---")
            generated_arrays = []
            try:
                rfd3_config = getattr(__import__('rfd3.engine', fromlist=['RFD3InferenceConfig']), 'RFD3InferenceConfig')(
                    diffusion_batch_size=min(10, NUM_RFD3_DESIGNS), 
                    low_memory_mode=False,
                    specification={'length': f"{total_output_length_min}-{total_output_length_max}", 'contig': full_contig_string, 'extra': {}}
                )
                rfd3_engine = RFD3InferenceEngine(**dataclasses.asdict(rfd3_config))
                
                spec_input = DesignInputSpecification(
                    input=pdb_path,
                    contig=full_contig_string, 
                    length=f"{total_output_length_min}-{total_output_length_max}",
                    extra={}
                )
                
                batches_needed = NUM_RFD3_DESIGNS // rfd3_config.diffusion_batch_size + (1 if NUM_RFD3_DESIGNS % rfd3_config.diffusion_batch_size != 0 else 0)
                
                rfd3_outputs_dict = rfd3_engine.run(inputs=spec_input, n_batches=batches_needed, out_dir=None)
                if rfd3_outputs_dict:
                    global_idx = 0
                    for key, rfd3_out_list in rfd3_outputs_dict.items():
                        if not key.startswith("backbone"): continue
                        for rfd3_out in rfd3_out_list:
                            if global_idx >= NUM_RFD3_DESIGNS: break
                            design_id = f"{target_name}_{combo_name}_rfd3_{global_idx}"
                            clean_array = renumber_atom_array_residues(rfd3_out.atom_array)
                            to_cif_file(clean_array, f"{combo_rfd3_out_dir}/{design_id}.cif", file_type="cif")
                            generated_arrays.append((design_id, clean_array))
                            global_idx += 1
            except Exception as e:
                print(f"Error during RFd3 running: {e}")
            finally:
                if 'rfd3_engine' in locals(): del rfd3_engine
                torch.cuda.empty_cache()

            # --- STAGE 2: RUN LigandMPNN ---
            print(f"--- Running LigandMPNN on {len(generated_arrays)} structures ---")
            lmpnn_jobs = []
            try:
                lmpnn_engine = MPNNInferenceEngine(model_type="ligand_mpnn", is_legacy_weights=True, write_structures=False, write_fasta=False)
                
                for design_id, rfd3_array in generated_arrays:
                    aa_sequence = get_sequence_from_array(rfd3_array, CHAIN_TO_DESIGN)
                    fixed_positions_A = []
                    current_fixed_start = 1
                    current_seq_idx = 0
                    
                    for loop in selected_loops:
                        flank_left, flank_right = loop["left"], loop["right"]
                        regex_pattern = re.compile(f"{flank_left}([A-Z]+?){flank_right}")
                        match = regex_pattern.search(aa_sequence[current_seq_idx:])
                        if match:
                            match_start = current_seq_idx + match.start()
                            match_end = current_seq_idx + match.end()
                            inserted_seq = match.group(1)
                            loop_start_1idx = match_start + len(flank_left) + 1
                            fixed_positions_A.extend(range(current_fixed_start, loop_start_1idx))
                            current_fixed_start = loop_start_1idx + len(inserted_seq)
                            current_seq_idx = match_end - len(flank_right)
                    
                    new_total_length = len(aa_sequence)
                    fixed_positions_A.extend(range(current_fixed_start, new_total_length + 1))
                    
                    fixed_residues_str = [f"{CHAIN_TO_DESIGN}{pos}" for pos in fixed_positions_A]
                    b_chain_mask = (rfd3_array.chain_id == FIXED_CHAINS[0]) & (rfd3_array.atom_name == "CA")
                    b_res_ids = np.unique(rfd3_array.res_id[b_chain_mask])
                    fixed_residues_str.extend([f"{FIXED_CHAINS[0]}{pos}" for pos in b_res_ids])

                    mpnn_input_dict = {
                        "name": design_id,
                        "batch_size": NUM_LMPNN_STRUCTURES,
                        "remove_waters": True,
                        "seed": 42,
                        "fixed_residues": fixed_residues_str,
                        "sampling_temp": LMPNN_SAMPLING_TEMP
                    }
                    
                    mpnn_outputs = lmpnn_engine.run(input_dicts=[mpnn_input_dict], atom_arrays=[rfd3_array])
                    
                    for seq_idx, mpnn_out in enumerate(mpnn_outputs):
                        valid_mask = ~np.isnan(mpnn_out.atom_array.coord[:, 0])
                        lmpnn_array = mpnn_out.atom_array[valid_mask]
                        lmpnn_array = renumber_atom_array_residues(lmpnn_array)
                        
                        full_seq_designed = get_sequence_from_array(lmpnn_array, CHAIN_TO_DESIGN)
                        to_cif_file(lmpnn_array, f"{combo_lmpnn_out_dir}/{design_id}_mpnn{seq_idx}.cif", file_type="cif")
                        
                        # Compute loop sequences
                        loop_data = {}
                        curr_idx = 0
                        for name_idx, loop in enumerate(selected_loops):
                            loop_name = combo[name_idx]
                            f_left, f_right = loop["left"], loop["right"]
                            m = re.search(f"{f_left}(.*?){f_right}", full_seq_designed[curr_idx:])
                            if m:
                                seq = m.group(1)
                                loop_data[f"loop_{loop_name}_seq"] = seq
                                loop_data[f"loop_{loop_name}_length"] = len(seq)
                                curr_idx += m.end() - len(f_right)
                            else:
                                loop_data[f"loop_{loop_name}_seq"] = "MISSING"
                                loop_data[f"loop_{loop_name}_length"] = 0
                        
                        lmpnn_jobs.append({
                            "target": target_name,
                            "combo": combo_name,
                            "design_id": design_id,
                            "seq_idx": seq_idx,
                            **loop_data,
                            "full_seq": full_seq_designed,
                            "lmpnn_array": lmpnn_array,
                            "rfd3_array": rfd3_array
                        })
            except Exception as e:
                print(f"Error during LigandMPNN sequence generation: {e}")
            finally:
                if 'lmpnn_engine' in locals(): del lmpnn_engine
                torch.cuda.empty_cache()

            # --- STAGE 3: RUN RF3 VALIDATION ---
            print(f"--- Running RF3 Validations on {len(lmpnn_jobs)} sequences ---")
            try:
                rf3_engine = RF3InferenceEngine(ckpt_path='rf3', verbose=False)
                for job in lmpnn_jobs:
                    design_id = f"{job['design_id']}_mpnn{job['seq_idx']}"
                    lmpnn_array = job['lmpnn_array']
                    rfd3_array = job['rfd3_array']
                    
                    valid_atoms = ['N', 'CA', 'C', 'O', 'CB']
                    array_for_rf3 = lmpnn_array[np.isin(lmpnn_array.atom_name, valid_atoms)].copy()
                    input_structure = InferenceInput.from_atom_array(array_for_rf3, example_id=design_id)
                    rf3_outputs_dict = rf3_engine.run(inputs=input_structure)
                    
                    rf3_target_key = next((k for k in rf3_outputs_dict.keys() if design_id in k), list(rf3_outputs_dict.keys())[0] if rf3_outputs_dict else None)
                    if rf3_target_key:
                        rf3_output = rf3_outputs_dict[rf3_target_key][0]
                        rf3_atom_array = renumber_atom_array_residues(rf3_output.atom_array)
                        to_cif_file(rf3_atom_array, f"{combo_rf3_out_dir}/{design_id}_refolded.cif", file_type="cif")

                        bb_mask_rfd3 = np.isin(rfd3_array.atom_name, PROTEIN_BACKBONE_ATOM_NAMES)
                        bb_mask_rf3 = np.isin(rf3_atom_array.atom_name, PROTEIN_BACKBONE_ATOM_NAMES)
                        bb_generated = rfd3_array[bb_mask_rfd3]
                        bb_refolded = rf3_atom_array[bb_mask_rf3]
                        bb_generated = bb_generated[bb_generated.atom_name != "OXT"]
                        bb_refolded = bb_refolded[bb_refolded.atom_name != "OXT"]
                        
                        if len(bb_generated) != len(bb_refolded):
                            min_len = min(len(bb_generated), len(bb_refolded))
                            bb_generated = bb_generated[:min_len]
                            bb_refolded = bb_refolded[:min_len]

                        bb_refolded_fitted, _ = superimpose(bb_generated, bb_refolded)
                        overall_rmsd = rmsd(bb_generated, bb_refolded_fitted)
                        
                        summary = rf3_output.summary_confidences
                        plddt = summary.get('overall_plddt', 0.0)
                        ptm = summary.get('ptm', 0.0)
                        
                        metrics = calc_protein_protein_metrics(rf3_atom_array, CHAIN_TO_DESIGN, FIXED_CHAINS[0])
                        heur_score = calculate_heuristic_score(
                            contacts=metrics["contacts"], 
                            interface_area=metrics["interface_area"], 
                            clashes=metrics["clashes"], 
                            centroid_distance=metrics["centroid_distance"], 
                            plddt_mean=plddt
                        )

                        loop_kwargs = {k: v for k, v in job.items() if k.startswith("loop_")}
                        
                        final_records.append({
                            "target": job["target"],
                            "loop_combo": job["combo"],
                            "design_id": design_id,
                            **loop_kwargs,
                            "overall_rmsd": overall_rmsd,
                            "plddt": plddt,
                            "ptm": ptm,
                            "contacts": metrics["contacts"],
                            "clashes": metrics["clashes"],
                            "interface_area": metrics["interface_area"],
                            "centroid_distance": metrics["centroid_distance"],
                            "heuristic_score": heur_score,
                            "full_seq": job["full_seq"]
                        })
            except Exception as e:
                print(f"Error during RF3 validation: {e}")
            finally:
                if 'rf3_engine' in locals(): del rf3_engine
                torch.cuda.empty_cache()
                
            # Intermediately save results after each combination to be safe
            df_interim = pd.DataFrame(final_records)
            df_interim.to_csv(os.path.join(OUT_BASE_DIR, "global_sequence_catalog.csv"), index=False)

    print("\n===============================")
    print("PIPELINE COMPLETED SUCCESSFULLY")
    print("===============================")
    
if __name__ == "__main__":
    main()
