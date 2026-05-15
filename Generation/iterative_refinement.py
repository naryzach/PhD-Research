import os
import time
import json
import dataclasses
import re
import pandas as pd
import numpy as np
import torch
import logging
from datetime import datetime, timedelta

# Foundry Inference imports
from lightning.fabric import seed_everything
from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine
from mpnn.inference_engines.mpnn import MPNNInferenceEngine
from rf3.inference_engines.rf3 import RF3InferenceEngine
from rf3.utils.inference import InferenceInput
from rfd3.inference.input_parsing import DesignInputSpecification
from biotite.structure.io.pdbx import CIFFile, get_structure as get_cif_structure
from biotite.structure.io.pdb import PDBFile
from biotite.structure import superimpose, rmsd
from atomworks.io.utils.io_utils import to_cif_file
from atomworks.constants import PROTEIN_BACKBONE_ATOM_NAMES
from biotite.sequence import ProteinSequence

# Configure Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Suppress noisy logs
logging.getLogger("transforms").setLevel(logging.ERROR)
logging.getLogger("atomworks.io").setLevel(logging.ERROR)
logging.getLogger("atomworks.ml").setLevel(logging.ERROR)
logging.getLogger("foundry").setLevel(logging.ERROR)

torch.set_float32_matmul_precision('medium')

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

class IterativeGenerator:
    def __init__(self, target_pdb_list, loop_names, state_path="refinement_state.json"):
        self.target_pdb_list = target_pdb_list
        self.loop_names = loop_names
        self.state_path = state_path
        self.out_base = "../Local/iterative_refinement"
        os.makedirs(self.out_base, exist_ok=True)
        
        self.loop_configs = {
            "AB": {"normal": 6, "max": 15, "pos": 30, "left": "LVK", "right": "LVY"},
            "C": {"normal": 6, "max": 15, "pos": 62, "left": "HTE", "right": "GLK"},
            "EF": {"normal": 4, "max": 10, "pos": 92, "left": "MYT", "right": "FVE"},
        }
        
        self.checkpoint_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Tools", "foundry_checkpoints"))
        os.environ["FOUNDRY_CHECKPOINT_DIRS"] = self.checkpoint_dir
        os.environ["DGLBACKEND"] = "pytorch"
        
        self.load_state()

    def load_state(self):
        if os.path.exists(self.state_path):
            try:
                with open(self.state_path, 'r') as f:
                    self.state = json.load(f)
                logger.info(f"Loaded state from {self.state_path}. Iteration: {self.state['iteration']}")
            except Exception as e:
                logger.error(f"Failed to load state: {e}. Starting fresh.")
                self.init_fresh_state()
        else:
            self.init_fresh_state()

    def init_fresh_state(self):
        self.state = {
            "iteration": 0,
            "hall_of_fame": [], # Top candidates
            "last_af3_export": datetime.now().isoformat(),
            "temperature": 0.5,
            "exported_sequences": []
        }
        self.save_state()

    def save_state(self):
        with open(self.state_path, 'w') as f:
            json.dump(self.state, f, indent=4)

    def export_for_af3(self):
        # Sort Hall of Fame by pLDDT
        sorted_hof = sorted(self.state["hall_of_fame"], key=lambda x: x.get("plddt", 0), reverse=True)
        top_15 = sorted_hof[:15]
        
        if not top_15:
            logger.warning("Hall of Fame is empty. Skipping AF3 export.")
            return

        af3_jobs = []
        for i, entry in enumerate(top_15):
            job_name = f"refine_it{self.state['iteration']}_{i:02d}"
            job = {
                "name": job_name,
                "model_seeds": [42],
                "sequences": [
                    {"proteinChain": {"sequence": entry.get("target_seq", ""), "count": 1}},
                    {"proteinChain": {"sequence": entry.get("full_seq", ""), "count": 1}}
                ]
            }
            af3_jobs.append(job)
            
        export_path = os.path.join(self.out_base, f"af3_submission_it{self.state['iteration']}.json")
        with open(export_path, 'w') as f:
            json.dump(af3_jobs, f, indent=2)
            
        logger.info(f"EXPORTED 15 sequences to {export_path}")
        self.state["last_af3_export"] = datetime.now().isoformat()
        self.save_state()
        
        print("\n" + "="*50)
        print(f"PAUSED: AlphaFold3 export generated at {export_path}")
        print("Please run these sequences on AF3 Server and provide 'af3_results.json' in this directory.")
        print("="*50 + "\n")
        
        while not os.path.exists("af3_results.json"):
            time.sleep(60)
            
        self.import_af3_results("af3_results.json")
        os.remove("af3_results.json")

    def import_af3_results(self, results_path):
        with open(results_path, 'r') as f:
            results = json.load(f)
        
        logger.info(f"Imported {len(results)} AF3 results. Updating Hall of Fame...")
        for res in results:
            res["source"] = "AF3"
            # Add or update metrics
            self.state["hall_of_fame"].append(res)
        
        # Keep only top 100 in HOF to avoid bloat
        self.state["hall_of_fame"] = sorted(self.state["hall_of_fame"], key=lambda x: x.get("plddt", 0), reverse=True)[:100]
        self.save_state()

    def run_rfd3(self, pdb_path, num_to_gen, out_dir):
        logger.info(f"Running RFd3 on {os.path.basename(pdb_path)} for {num_to_gen} backbones...")
        os.makedirs(out_dir, exist_ok=True)
        
        # Hardcoded loop params for now (TIMP3 specific)
        total_length = 121
        fixed_chains = ["B"]
        chain_to_design = "A"
        
        selected_loops = [self.loop_configs[name] for name in self.loop_names]
        selected_loops.sort(key=lambda x: x["pos"])
        
        contig_parts = []
        current_pos = 1
        for loop in selected_loops:
            if current_pos <= loop["pos"]:
                contig_parts.append(f"{chain_to_design}{current_pos}-{loop['pos']}")
            contig_parts.append(f"{loop['normal']}-{loop['max']}")
            current_pos = loop["pos"] + loop["normal"] + 1
        if current_pos <= total_length:
            contig_parts.append(f"{chain_to_design}{current_pos}-{total_length}")
        
        contig_string = ",".join(contig_parts)
        
        # Fix chain length
        try:
            from biotite.structure.io.pdb import PDBFile
            pdb_f = PDBFile.read(pdb_path)
            atom_array = pdb_f.get_structure()[0]
            fix_chain_len = len(atom_array[atom_array.chain_id == fixed_chains[0]])
        except:
            fix_chain_len = 160 # Fallback
            
        full_contig_string = f"{contig_string},/0,{fixed_chains[0]}1-{fix_chain_len}"
        
        rfd3_config = RFD3InferenceConfig(
            diffusion_batch_size=min(10, num_to_gen),
            specification={'length': f"{total_length+fix_chain_len-10}-{total_length+fix_chain_len+20}", 'contig': full_contig_string, 'extra': {}}
        )
        rfd3_engine = RFD3InferenceEngine(**dataclasses.asdict(rfd3_config))
        
        spec_input = DesignInputSpecification(
            input=pdb_path,
            contig=full_contig_string,
            length=f"{total_length+fix_chain_len-10}-{total_length+fix_chain_len+20}",
            extra={}
        )
        
        batches_needed = (num_to_gen + rfd3_config.diffusion_batch_size - 1) // rfd3_config.diffusion_batch_size
        rfd3_outputs_dict = rfd3_engine.run(inputs=spec_input, n_batches=batches_needed, out_dir=None)
        
        generated_arrays = []
        if rfd3_outputs_dict:
            global_idx = 0
            for key, rfd3_out_list in rfd3_outputs_dict.items():
                if not key.startswith("backbone"): continue
                for rfd3_out in rfd3_out_list:
                    if global_idx >= num_to_gen: break
                    clean_array = renumber_atom_array_residues(rfd3_out.atom_array)
                    generated_arrays.append(clean_array)
                    global_idx += 1
        
        del rfd3_engine
        torch.cuda.empty_cache()
        return generated_arrays

    def run_lmpnn(self, backbones, out_dir, temp=0.1):
        logger.info(f"Running LigandMPNN on {len(backbones)} backbones (T={temp})...")
        os.makedirs(out_dir, exist_ok=True)
        
        lmpnn_engine = MPNNInferenceEngine(model_type="ligand_mpnn", is_legacy_weights=True, write_structures=False, out_directory=out_dir)
        results = []
        
        for idx, array in enumerate(backbones):
            # Simplified fixed residues logic (fix everything except loops)
            aa_seq = get_sequence_from_array(array, "A")
            fixed_residues_str = [] # Placeholder: redesign all of A for now
            
            mpnn_input = {
                "name": f"design_{idx}",
                "batch_size": 1,
                "fixed_residues": fixed_residues_str,
                "sampling_temp": temp
            }
            
            try:
                mpnn_outputs = lmpnn_engine.run(input_dicts=[mpnn_input], atom_arrays=[array])
                for mp_out in mpnn_outputs:
                    results.append({
                        "array": renumber_atom_array_residues(mp_out.atom_array),
                        "parent_array": array,
                        "full_seq": get_sequence_from_array(mp_out.atom_array, "A"),
                        "target_seq": get_sequence_from_array(mp_out.atom_array, "B")
                    })
            except Exception as e:
                logger.error(f"LMPNN error: {e}")
                
        del lmpnn_engine
        torch.cuda.empty_cache()
        return results

    def run_rf3(self, candidates, out_dir):
        logger.info(f"Running RF3 validation on {len(candidates)} candidates...")
        os.makedirs(out_dir, exist_ok=True)
        
        rf3_engine = RF3InferenceEngine(ckpt_path='rf3', verbose=False)
        scored_candidates = []
        
        for idx, cand in enumerate(candidates):
            design_id = f"cand_{idx}"
            try:
                # Prepare input
                valid_atoms = ['N', 'CA', 'C', 'O', 'CB']
                array_for_rf3 = cand["array"][np.isin(cand["array"].atom_name, valid_atoms)].copy()
                input_struct = InferenceInput.from_atom_array(array_for_rf3, example_id=design_id)
                
                rf3_out_dict = rf3_engine.run(inputs=input_struct)
                rf3_key = next(iter(rf3_out_dict.keys()))
                rf3_output = rf3_out_dict[rf3_key][0]
                
                metrics = rf3_output.summary_confidences
                plddt = metrics.get('overall_plddt', 0.0)
                ptm = metrics.get('ptm', 0.0)
                
                # Calculate RMSD to backbone
                bb_gen = cand["parent_array"][np.isin(cand["parent_array"].atom_name, PROTEIN_BACKBONE_ATOM_NAMES)]
                bb_ref = rf3_output.atom_array[np.isin(rf3_output.atom_array.atom_name, PROTEIN_BACKBONE_ATOM_NAMES)]
                
                # Align and RMSD
                min_len = min(len(bb_gen), len(bb_ref))
                fitted, _ = superimpose(bb_gen[:min_len], bb_ref[:min_len])
                rmsd_val = rmsd(bb_gen[:min_len], fitted)
                
                cand_result = {
                    **cand,
                    "plddt": plddt,
                    "ptm": ptm,
                    "rmsd": rmsd_val,
                    "iteration": self.state["iteration"]
                }
                # Don't store the full array in state (too big), just metadata and sequence
                state_entry = {k: v for k, v in cand_result.items() if k not in ["array", "parent_array"]}
                scored_candidates.append(state_entry)
                
                # Save CIF
                to_cif_file(renumber_atom_array_residues(rf3_output.atom_array), os.path.join(out_dir, f"{design_id}_rf3.cif"), file_type="cif")
                
            except Exception as e:
                logger.error(f"RF3 error on {design_id}: {e}")
                
        del rf3_engine
        torch.cuda.empty_cache()
        return scored_candidates

    def run_iteration(self):
        it = self.state["iteration"]
        logger.info(f"--- Iteration {it} ---")
        it_out = os.path.join(self.out_base, f"it_{it}")
        os.makedirs(it_out, exist_ok=True)
        
        all_new_scored = []
        
        for pdb_file in self.target_pdb_list:
            pdb_path = os.path.join("../Data/TIMP_Complexes/HADDOCK_PDB", pdb_file)
            if not os.path.exists(pdb_path): continue
            
            # 1. Generate Backbones
            backbones = self.run_rfd3(pdb_path, 10, os.path.join(it_out, "rfd3"))
            
            # 2. Design Sequences
            candidates = self.run_lmpnn(backbones, os.path.join(it_out, "lmpnn"), temp=self.state["temperature"])
            
            # 3. Validate
            scored = self.run_rf3(candidates, os.path.join(it_out, "rf3"))
            all_new_scored.extend(scored)
            
        # Update HOF
        self.state["hall_of_fame"].extend(all_new_scored)
        self.state["hall_of_fame"] = sorted(self.state["hall_of_fame"], key=lambda x: x.get("plddt", 0), reverse=True)[:100]
        self.save_state()

    def main_loop(self):
        logger.info("Starting main optimization loop...")
        while True:
            # Time-based AF3 export
            last_export = datetime.fromisoformat(self.state["last_af3_export"])
            if datetime.now() - last_export > timedelta(hours=18):
                self.export_for_af3()
            
            self.run_iteration()
            
            self.state["iteration"] += 1
            self.state["temperature"] = max(0.1, self.state["temperature"] * 0.9)
            self.save_state()
            logger.info(f"Iteration {self.state['iteration']} complete. Temp: {self.state['temperature']:.2f}")

if __name__ == "__main__":
    pdbs = ["TIMP3_vs_MMP9_HADDOCK_Xray.pdb", "TIMP3_vs_MMP2_HADDOCK_Xray.pdb"]
    generator = IterativeGenerator(pdbs, ["C", "EF"])
    generator.main_loop()
