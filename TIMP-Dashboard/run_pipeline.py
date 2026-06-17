import os
import subprocess

# --- PORTABILITY: GPU-Aware Environment Setup ---
# We must detect the GPU and set DISABLE_CUEQUIVARIANCE before importing heavy ML libraries.
# cuEquivariance checks this at import time, so inside main() is too late.
try:
    # Use nvidia-smi to check for V100
    smi_out = subprocess.check_output(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"]).decode()
    if "V100" in smi_out:
        os.environ["DISABLE_CUEQUIVARIANCE"] = "1"
        print(f"Detected V100 GPU. Automatically setting DISABLE_CUEQUIVARIANCE=1 for compatibility.")
except Exception:
    pass

import sys
import json
import numpy as np
import pandas as pd
import argparse
import dataclasses
import gc
import torch
import time
import re
import subprocess
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
from rfd3.inference.input_parsing import DesignInputSpecification

torch.set_float32_matmul_precision('medium')

SAVE_ESM_STRUCTS = True   # write predicted CIFs for downstream inspection

# ── ESMFold2 helpers ──────────────────────────────────────────────────────────
# ESMFold2 (Chan Zuckerberg Biohub / Rives lab) folds binder + target as a
# two-chain complex and returns esm_iptm / esm_ptm / esm_plddt (binder mean).
# Install alongside the foundry env or in a compatible one:
#   pip install "esm @ git+https://github.com/Biohub/esm.git@main"
#   pip install transformers
# The three functions below are the complete interface — no external script needed.

def _esm_load_model(device: str = "cuda"):
    """Load ESMFold2 once and return (model, input_builder)."""
    from esm.models.esmfold2 import ESMFold2InputBuilder
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
    model = ESMFold2Model.from_pretrained("biohub/ESMFold2").to(device).eval()
    return model, ESMFold2InputBuilder()


def _esm_save_structure(res, stub: str):
    """Write the predicted complex to {stub}.cif (or .pdb). Returns path or None."""
    comp = getattr(res, "complex", None) or getattr(res, "structure", None)
    if comp is None:
        return None
    for meth, ext in (("to_mmcif", "cif"), ("to_mmcif_string", "cif"),
                      ("to_pdb", "pdb"), ("to_pdb_string", "pdb")):
        fn = getattr(comp, meth, None)
        if callable(fn):
            try:
                text = fn()
                if isinstance(text, str) and text.strip():
                    path = f"{stub}.{ext}"
                    with open(path, "w") as fh:
                        fh.write(text)
                    return path
            except Exception:
                continue
    for meth, ext in (("save_mmcif", "cif"), ("save_pdb", "pdb"), ("save", "cif")):
        fn = getattr(comp, meth, None)
        if callable(fn):
            try:
                path = f"{stub}.{ext}"
                fn(path)
                return path
            except Exception:
                continue
    return None


def _write_esm_confidences(path, plddt_full, pae_mat, ptm=None, iptm=None):
    """
    Persist ESMFold2 per-residue pLDDT (0–1) and PAE (Å) to a confidences JSON in
    the shape analyze_results.py / the dashboard read (keys: plddt, pae, ptm, iptm).
    Written only when pLDDT or PAE is available; PAE is omitted when the model does
    not expose it (downstream degrades gracefully).
    """
    if plddt_full is None and pae_mat is None:
        return
    conf = {}
    if ptm is not None and not (isinstance(ptm, float) and np.isnan(ptm)):
        conf["ptm"] = float(ptm)
    if iptm is not None and not (isinstance(iptm, float) and np.isnan(iptm)):
        conf["iptm"] = float(iptm)
    if plddt_full is not None:
        conf["plddt"] = [round(float(x), 4) for x in np.asarray(plddt_full).flatten()]
    if pae_mat is not None:
        conf["pae"] = np.round(np.asarray(pae_mat, dtype=float), 2).tolist()
    try:
        with open(path, "w") as fh:
            json.dump(conf, fh)
    except Exception as e:
        print(f"  Could not write confidences JSON {path}: {e}", flush=True)


def _esm_predict_complex(model, builder, binder_seq: str, target_seq: str,
                          binder_len: int, seed: int = 0,
                          cif_stub: str = None) -> dict:
    """
    Fold a two-chain complex and return scalar scores plus the raw per-residue
    arrays for the confidences dump:
      {esm_iptm, esm_ptm, esm_plddt, _plddt_full, _pae[, esm_cif]}
    esm_plddt is the binder-chain mean on the 0–1 scale this dashboard pipeline
    uses (matching the legacy RF3 convention).  _plddt_full is the full-complex
    per-residue pLDDT (0–1); _pae is the NxN PAE matrix (Å) or None if the model
    does not expose it.  If cif_stub is given the structure is written to
    {cif_stub}.cif/.pdb at no extra GPU cost.
    """
    from esm.models.esmfold2 import ProteinInput, StructurePredictionInput

    spi = StructurePredictionInput(sequences=[
        ProteinInput(id="A", sequence=binder_seq),
        ProteinInput(id="B", sequence=target_seq),
    ])
    res = builder.fold(model, spi, num_loops=3, num_sampling_steps=50,
                       num_diffusion_samples=1, seed=seed)

    def _get(obj, *names):
        for n in names:
            v = getattr(obj, n, None)
            if v is not None:
                return v
        return None

    iptm = _get(res, "iptm", "interface_ptm")
    ptm  = _get(res, "ptm")

    # Full per-residue pLDDT (both chains), normalized to the 0–1 scale the TIMP
    # dashboard / analyze_results pipeline uses.
    plddt_full = None
    plddt_arr = _get(res, "plddt", "plddts")
    if plddt_arr is not None:
        a = np.asarray(plddt_arr, dtype=float).flatten()
        if a.size:
            if np.nanmax(a) > 1.0:           # model emitted a 0–100 scale
                a = a / 100.0
            plddt_full = a

    # Scalar pLDDT = binder-chain (first binder_len residues) mean, 0–1.
    plddt = np.nan
    if plddt_full is not None and plddt_full.size:
        binder_part = plddt_full[:binder_len] if plddt_full.size >= binder_len else plddt_full
        plddt = float(np.mean(binder_part))

    # PAE matrix (Å), if this build of ESMFold2 exposes one.
    pae_mat = None
    pae_raw = _get(res, "pae", "predicted_aligned_error", "pae_matrix", "aligned_error")
    if pae_raw is not None:
        try:
            pm = np.asarray(pae_raw, dtype=float)
            if pm.ndim == 3:
                pm = pm[0]
            if pm.ndim == 2:
                pae_mat = pm
        except Exception:
            pae_mat = None

    def _f(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return np.nan

    out = {"esm_iptm": _f(iptm), "esm_ptm": _f(ptm), "esm_plddt": plddt,
           "_plddt_full": plddt_full, "_pae": pae_mat}
    if cif_stub:
        saved = _esm_save_structure(res, cif_stub)
        if saved:
            out["esm_cif"] = saved
    return out

import logging
import warnings
logging.getLogger("transforms").setLevel(logging.ERROR)
logging.getLogger("atomworks.io").setLevel(logging.ERROR)
logging.getLogger("atomworks.ml").setLevel(logging.ERROR)
logging.getLogger("foundry").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", module="biotite")

# --- USER CONFIGURABLE VARIABLES ---
NUM_RFD3_DESIGNS = 5
NUM_LMPNN_STRUCTURES = 5
LMPNN_SAMPLING_TEMP = 0.1

TARGET_PDBS = [
    "TIMP3_vs_ADAM10_AF.cif",
    "TIMP3_vs_MMP10_AF.cif",
    "TIMP3_vs_MMP3_AF.cif",
    "TIMP3_vs_ADAM17_AF.cif",
    "TIMP3_vs_MMP2_AF.cif",
    "TIMP3_vs_MMP9_AF.cif"
]

LOOP_DEFINITIONS = {
    "AB": {"normal": 6, "max": 15, "pos": 30, "left": "LVK", "right": "LVY"},
    "C": {"normal": 6, "max": 15, "pos": 62, "left": "HTE", "right": "GLK"},
    "EF": {"normal": 4, "max": 10, "pos": 92, "left": "MYT", "right": "FVE"},
    "GH": {"normal": 10, "max": 20, "pos": 127, "left": "KSC", "right": "NEC"},
    "Multi": {"normal": 10, "max": 20, "pos": 143, "left": "LWT", "right": "YQS"}
}

DATA_DIR = "../Data/TIMP_Complexes/AlphaFold_CIF"
OUT_BASE_DIR = "../Local/TIMP-Dashboard_output"

CHAIN_TO_DESIGN = "A" # TIMP3 chain
FIXED_CHAINS = ["B"]  # Target chain

VARY_LOOP_LENGTHS = False
if not VARY_LOOP_LENGTHS:
    for loop_combo in LOOP_DEFINITIONS.keys():
        LOOP_DEFINITIONS[loop_combo]["max"] = LOOP_DEFINITIONS[loop_combo]["normal"]

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

def calculate_heuristic_score(contacts, interface_area, clashes, centroid_distance, plddt_mean, iptm=0.0):
    interface_area = 0.0 if np.isnan(interface_area) else interface_area
    score = (contacts * 0.4 
             + interface_area * 0.002 
             - clashes * 3.0 
             - centroid_distance * 0.2)
    if not np.isnan(plddt_mean):
        score += plddt_mean * 0.5
    if iptm:
        score += iptm * 100.0
    return score

# --- GENERATION PIPELINE ---
def main():
    #seed_everything(42)
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

    for target_pdb in tqdm(TARGET_PDBS, desc="Targets", unit="target"):
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

        for combo in tqdm(valid_combos, desc=f" {target_name} loops", leave=False):
            combo_name = "_".join(combo)
            print(f"\n======================================")
            print(f"Processing Target: {target_name} | Loops: {combo_name}")
            print(f"======================================")

            combo_rfd3_out_dir = os.path.join(OUT_BASE_DIR, target_name, combo_name, "rfd3")
            combo_lmpnn_out_dir = os.path.join(OUT_BASE_DIR, target_name, combo_name, "lmpnn")
            combo_esm_out_dir = os.path.join(OUT_BASE_DIR, target_name, combo_name, "esmfold2")
            os.makedirs(combo_rfd3_out_dir, exist_ok=True)
            os.makedirs(combo_lmpnn_out_dir, exist_ok=True)
            os.makedirs(combo_esm_out_dir, exist_ok=True)

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
                # Detect GPU and set precision overrides
                precision = "bf16-mixed"
                if torch.cuda.is_available():
                    device_name = torch.cuda.get_device_name(0)
                    if "V100" in device_name:
                        precision = "16-mixed"
                    print(f"Detected GPU: {device_name}. Using precision: {precision}, DISABLE_CUEQUIVARIANCE: {os.environ.get('DISABLE_CUEQUIVARIANCE', '0')}")

                rfd3_config = getattr(__import__('rfd3.engine', fromlist=['RFD3InferenceConfig']), 'RFD3InferenceConfig')(
                    diffusion_batch_size=min(5, NUM_RFD3_DESIGNS), 
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
                
                for design_id, rfd3_array in tqdm(generated_arrays, desc="Designing Seqs", leave=False):
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
                        
                        # Initialize all possible loop columns for consistent schema
                        loop_data = {
                            "loop_AB_seq": "MISSING", "loop_AB_length": 0,
                            "loop_C_seq": "MISSING", "loop_C_length": 0,
                            "loop_EF_seq": "MISSING", "loop_EF_length": 0,
                            "loop_GH_seq": "MISSING", "loop_GH_length": 0,
                            "loop_Multi_seq": "MISSING", "loop_Multi_length": 0
                        }
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

            # --- STAGE 3: ESMFold2 VALIDATION ---
            # ESMFold2 (MSA-free, language-model-based) replaced RF3 as the fold validator.
            # Binder (chain A) + target (chain B) are folded as a two-chain complex, returning
            # esm_iptm, esm_ptm, and esm_plddt (binder-chain mean).
            # Backbone RMSD is computed vs the RFd3 chain A design.
            # Interface metrics (contacts, area, clashes) are also computed here from the
            # lmpnn_array, which retains both chains with the RFd3-placed coordinates.
            print(f"--- Running ESMFold2 Validations on {len(lmpnn_jobs)} sequences ---")

            # ── 3a. Score sequences directly via ESMFold2 (in-process) ──────────
            #   esm_scores[design_id] keeps the scalar metrics (for esm_scores.csv
            #   and Stage 3b).  The per-residue pLDDT + PAE arrays are detached and
            #   written to <design_id>_confidences.json so analyze_results.py and
            #   the dashboard can recover loop pLDDT / PAE / chain-split confidence
            #   the way they did for RF3.
            esm_scores = {}
            _esm_model = _esm_builder = None
            try:
                device = "cuda" if torch.cuda.is_available() else "cpu"
                print(f"Loading ESMFold2 on {device} ...", flush=True)
                _esm_model, _esm_builder = _esm_load_model(device)

                for job in tqdm(lmpnn_jobs, desc="ESMFold2 Scoring", leave=False):
                    design_id  = f"{job['design_id']}_mpnn{job['seq_idx']}"
                    binder_seq = get_sequence_from_array(job['lmpnn_array'], CHAIN_TO_DESIGN)
                    target_seq = get_sequence_from_array(job['lmpnn_array'], FIXED_CHAINS[0])
                    if not binder_seq or not target_seq:
                        continue
                    try:
                        stub = os.path.join(combo_esm_out_dir, design_id) if SAVE_ESM_STRUCTS else None
                        metrics = _esm_predict_complex(
                            _esm_model, _esm_builder, binder_seq, target_seq,
                            len(binder_seq), cif_stub=stub,
                        )
                        # Detach the heavy per-residue arrays and persist them as a
                        # confidences JSON; keep only scalars in esm_scores.
                        plddt_full = metrics.pop("_plddt_full", None)
                        pae_mat    = metrics.pop("_pae", None)
                        _write_esm_confidences(
                            os.path.join(combo_esm_out_dir, f"{design_id}_confidences.json"),
                            plddt_full, pae_mat,
                            ptm=metrics.get("esm_ptm"), iptm=metrics.get("esm_iptm"),
                        )
                        esm_scores[design_id] = metrics
                        print(f"  {design_id}: ipTM={metrics['esm_iptm']:.3f}  "
                              f"pLDDT={metrics['esm_plddt']:.3f}", flush=True)
                    except Exception as e:
                        print(f"  ESMFold2 failed on {design_id}: {e}", flush=True)

                if esm_scores:
                    pd.DataFrame([{"design_id": did, **m} for did, m in esm_scores.items()]).to_csv(
                        os.path.join(combo_esm_out_dir, "esm_scores.csv"), index=False)
                    print(f"ESMFold2 scored {len(esm_scores)}/{len(lmpnn_jobs)} designs.", flush=True)

            except ImportError as e:
                print(
                    f"ESMFold2 not available ({e}).\n"
                    "  pip install 'esm @ git+https://github.com/Biohub/esm.git@main' transformers\n"
                    "  pLDDT/pTM/iPTM will be NaN; interface metrics will still be computed.",
                    flush=True,
                )
            finally:
                del _esm_model, _esm_builder
                torch.cuda.empty_cache()

            # ── 3b. Per-design: RMSD from ESMFold2 structure + interface metrics ──
            bb_atom_names = set(PROTEIN_BACKBONE_ATOM_NAMES) - {"OXT"}
            try:
                for job in tqdm(lmpnn_jobs, desc="ESMFold2 Post-processing", leave=False):
                    design_id = f"{job['design_id']}_mpnn{job['seq_idx']}"
                    lmpnn_array = job['lmpnn_array']
                    rfd3_array = job['rfd3_array']

                    esm_m  = esm_scores.get(design_id, {})
                    plddt  = float(esm_m.get("esm_plddt", float("nan")))
                    ptm    = float(esm_m.get("esm_ptm",   float("nan")))
                    iptm   = float(esm_m.get("esm_iptm",  float("nan")))
                    esm_cif = esm_m.get("esm_cif")

                    # ── RMSD: RFd3 chain A backbone vs ESMFold2 monomer refold ──
                    overall_rmsd = float("nan")
                    if esm_cif and os.path.exists(str(esm_cif)):
                        try:
                            if str(esm_cif).endswith(".cif"):
                                esm_array = get_cif_structure(CIFFile.read(str(esm_cif)), model=1)
                            else:
                                esm_array = PDBFile.read(str(esm_cif)).get_structure()[0]

                            rfd3_bb = rfd3_array[
                                (rfd3_array.chain_id == CHAIN_TO_DESIGN) &
                                np.isin(rfd3_array.atom_name, list(bb_atom_names))]
                            esm_bb  = esm_array[np.isin(esm_array.atom_name, list(bb_atom_names))]

                            n = min(len(rfd3_bb), len(esm_bb))
                            if n > 0:
                                esm_bb_fitted, _ = superimpose(rfd3_bb[:n], esm_bb[:n])
                                overall_rmsd = float(rmsd(rfd3_bb[:n], esm_bb_fitted))
                        except Exception as e:
                            print(f"  RMSD computation failed for {design_id}: {e}")

                    # ── Interface metrics from lmpnn_array (both chains retained) ──
                    metrics = calc_protein_protein_metrics(lmpnn_array, CHAIN_TO_DESIGN, FIXED_CHAINS[0])
                    heur_score = calculate_heuristic_score(
                        contacts=metrics["contacts"],
                        interface_area=metrics["interface_area"],
                        clashes=metrics["clashes"],
                        centroid_distance=metrics["centroid_distance"],
                        plddt_mean=plddt,
                        iptm=iptm if not np.isnan(iptm) else 0.0
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
                        "iptm": iptm,
                        "contacts": metrics["contacts"],
                        "clashes": metrics["clashes"],
                        "interface_area": metrics["interface_area"],
                        "centroid_distance": metrics["centroid_distance"],
                        "heuristic_score": heur_score,
                        "full_seq": job["full_seq"]
                    })
            except Exception as e:
                print(f"Error during ESMFold2 post-processing: {e}")
            finally:
                torch.cuda.empty_cache()
                
            # Intermediately save results after each combination to be safe
            df_interim = pd.DataFrame(final_records)
            df_interim.to_csv(os.path.join(OUT_BASE_DIR, "global_sequence_catalog.csv"), index=False)

    print("\n===============================")
    print("PIPELINE COMPLETED SUCCESSFULLY")
    print("===============================")
    
def run_reflexive_updates():
    """Run the analysis and cross-docking scripts to update the dashboard data."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    print("\n--- Running Reflexive Post-Processing Updates ---")
    
    # 1. Structural Analysis & Best Binder Ranking
    print(f"Updating structural analysis via analyze_results.py...")
    try:
        subprocess.run([sys.executable, os.path.join(script_dir, "analyze_results.py")], check=True)
    except Exception as e:
        print(f"ERROR running analyze_results.py: {e}")
        
    # 2. Cross-Docking & Specificity Analysis
    # We only run this if cross_docking_analysis.py exists (it might be in development)
    xd_path = os.path.join(script_dir, "cross_docking_analysis.py")
    if os.path.exists(xd_path):
        print(f"Updating specificity analysis via cross_docking_analysis.py...")
        try:
            subprocess.run([sys.executable, xd_path], check=True)
        except Exception as e:
            print(f"ERROR running cross_docking_analysis.py: {e}")

    print("--- Post-Processing Updates Complete ---")
    
if __name__ == "__main__":
    main()
    run_reflexive_updates()
