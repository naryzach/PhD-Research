import os, time, sys
import dataclasses
import re
import pandas as pd
from biotite.structure.io.pdbx import CIFFile, get_structure as get_cif_structure
from biotite.structure.io.pdb import PDBFile
from biotite.structure import superimpose, rmsd
from atomworks.io.utils.io_utils import to_cif_file
from atomworks.constants import PROTEIN_BACKBONE_ATOM_NAMES
from biotite.sequence import ProteinSequence
import numpy as np
import torch

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

def get_sequence_from_array(atom_array, chain_id="A"):
    """Helper function to extract sequence cleanly from an AtomArray."""
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

def get_pdb_length(pdb_path, chain_id):
    """Fallback sequence logic for PDB parsing."""
    try:
        pdb_file = PDBFile.read(pdb_path)
        atom_array = pdb_file.get_structure()[0]
        return len(get_sequence_from_array(atom_array, chain_id))
    except Exception:
        # Fallback to text parsing
        count = 0
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM') and line.split()[4] == chain_id and line.split()[2] == "CA":
                    count += 1
        return count

def renumber_atom_array_residues(atom_array):
    """Renumber residue IDs contiguously starting from 1 for each chain."""
    new_res_ids = np.zeros(len(atom_array), dtype=int)
    for chain_id in np.unique(atom_array.chain_id):
        chain_mask = atom_array.chain_id == chain_id
        chain_res_ids = atom_array.res_id[chain_mask]
        
        # Keep track of old residue IDs and transition them linearly
        unique_old_ids = []
        last_id = None
        for r_id in chain_res_ids:
            if r_id != last_id:
                unique_old_ids.append(r_id)
                last_id = r_id
                
        # Map old IDs to new 1-indexed contiguous IDs
        id_map = {old_id: new_id for new_id, old_id in enumerate(unique_old_ids, start=1)}
        
        # Apply the mapping
        new_res_ids[chain_mask] = [id_map[old_id] for old_id in chain_res_ids]
        
    atom_array.res_id = new_res_ids
    return atom_array

def main():
    #seed_everything(42)
    checkpoint_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Tools", "foundry_checkpoints"))
    os.makedirs(checkpoint_dir, exist_ok=True)
    os.environ["FOUNDRY_CHECKPOINT_DIRS"] = checkpoint_dir
    os.environ["DGLBACKEND"] = "pytorch"

    original_directory = os.getcwd()
    data_path = "../Data"

    output_prefix = "design"

    loop_names = ["C", "EF"]
    if isinstance(loop_names, str):
        loop_names = [loop_names]

    loop_configs = {
        "AB": {"normal": 6, "max": 15, "pos": 30, "left": "LVK", "right": "LVY"},
        "C": {"normal": 6, "max": 15, "pos": 62, "left": "HTE", "right": "GLK"},
        "EF": {"normal": 4, "max": 10, "pos": 92, "left": "MYT", "right": "FVE"},
        "GH": {"normal": 10, "max": 20, "pos": 127, "left": "KSC", "right": "NEC"},
        "Multi": {"normal": 10, "max": 20, "pos": 143, "left": "LWT", "right": "YQS"}
    }

    for name in loop_names:
        if name not in loop_configs:
            raise Exception(f"Specified Loop '{name}' is unknown")

    selected_loops = [loop_configs[name] for name in loop_names]
    selected_loops.sort(key=lambda x: x["pos"])

    chain_to_design = "A"
    fixed_chains = ["B"]
    total_length = 121

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
    # NOTE: In foundry rfd3, chain breaks are designated by adding the next chain ID without explicit /0.
    
    num_sequences_to_generate = 25

    print(f"Targeting loops: {loop_names}")
    
    pdb_file_list = [
        "TIMP3_vs_ADAM10_HADDOCK_Xray.pdb",
        "TIMP3_vs_ADAM17_HADDOCK_Xray.pdb",
        "TIMP3_vs_MMP2_HADDOCK_Xray.pdb",
        "TIMP3_vs_MMP7_HADDOCK_Xray.pdb",
        "TIMP3_vs_MMP9_HADDOCK_Xray.pdb",
        "TIMP3_vs_MMP10_HADDOCK_Xray.pdb"
    ]

    rfd3_out_base = "../Local/rfd3_output"
    lmpnn_out_base = "../Local/ligandmpnn_output"

    for pdb_complex_file_name in pdb_file_list:
        pdb_path = os.path.join(data_path, pdb_complex_file_name)
        out_name = pdb_complex_file_name.replace(".pdb", "").replace(".cif", "")
        
        rfd3_out_dir = os.path.join(rfd3_out_base, out_name)
        lmpnn_out_dir = os.path.join(lmpnn_out_base, out_name)
        
        os.makedirs(rfd3_out_dir, exist_ok=True)
        os.makedirs(lmpnn_out_dir, exist_ok=True)

        if not os.path.exists(pdb_path):
            print(f"Skipping {pdb_complex_file_name}, does not exist.")
            continue

        fix_chain_len = get_pdb_length(pdb_path, fixed_chains[0])
        print(f"B chain length: {fix_chain_len}")
        
        full_contig_string = f"{contig_string},/0,{fixed_chains[0]}1-{fix_chain_len}"
        print(f"Final Contig String: {full_contig_string}")
        
        # Calculate dynamic length constraints
        min_inpaint_len = sum(loop["normal"] for loop in selected_loops)
        max_inpaint_len = sum(loop["max"] for loop in selected_loops)
        
        native_removed_len = 0
        for loop in selected_loops:
            # Each loop replaces the context between flank_left and flank_right
            # In single-chain RFdiffusion, "pos" is the start of the loop, and "normal" was its native length
            native_removed_len += loop["normal"] 
            
        base_length = total_length + fix_chain_len - native_removed_len
        total_output_length_min = base_length + min_inpaint_len
        total_output_length_max = base_length + max_inpaint_len

        # --- RUN RFd3 ---
        print("\n--- Running RFd3 ---")
        try:
            rfd3_config = RFD3InferenceConfig(
                diffusion_batch_size=min(10, num_sequences_to_generate), 
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
            
            batches_needed = num_sequences_to_generate // rfd3_config.diffusion_batch_size + (1 if num_sequences_to_generate % rfd3_config.diffusion_batch_size != 0 else 0)
            
            st = time.time()
            rfd3_outputs_dict = rfd3_engine.run(inputs=spec_input, n_batches=batches_needed, out_dir=None)
            end = time.time()
            print(f"RFd3 finished in {(end-st)/60:.2f} minutes.")
            
            generated_arrays = []
            
            if rfd3_outputs_dict:
                global_idx = 0
                for key, rfd3_out_list in rfd3_outputs_dict.items():
                    if not key.startswith("backbone"): continue
                    
                    for rfd3_out in rfd3_out_list:
                        if global_idx >= num_sequences_to_generate: break
                        design_id = f"{output_prefix}_{global_idx}"
                        
                        clean_array = renumber_atom_array_residues(rfd3_out.atom_array)
                        
                        to_cif_file(clean_array, f"{rfd3_out_dir}/{design_id}.cif", file_type="cif")
                        generated_arrays.append((design_id, clean_array))
                        global_idx += 1
            
            del rfd3_engine
            torch.cuda.empty_cache()
            
        except Exception as e:
            print(f"Error during RFd3 running: {e}")
            continue

        # --- RUN LigandMPNN ---
        print(f"\n--- Running LigandMPNN on {len(generated_arrays)} structures ---")
        
        lmpnn_engine = MPNNInferenceEngine(model_type="ligand_mpnn", is_legacy_weights=True, write_structures=False, write_fasta=True, out_directory=lmpnn_out_dir)
        metric_records = []
        final_records = []
        
        st = time.time()
        for idx, (design_id, rfd3_array) in enumerate(generated_arrays):
            print(f"Designing {design_id}...")
            try:
                # Find fixed residues dynamically for each rfd3 array instance based on flanks
                aa_sequence = get_sequence_from_array(rfd3_array, chain_to_design)
                
                fixed_positions_A = []
                current_fixed_start = 1
                current_seq_idx = 0
                
                for loop in selected_loops:
                    flank_left = loop["left"]
                    flank_right = loop["right"]
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
                
                fixed_residues_str = [f"{chain_to_design}{pos}" for pos in fixed_positions_A]
                
                b_chain_mask = (rfd3_array.chain_id == fixed_chains[0]) & (rfd3_array.atom_name == "CA")
                b_res_ids = np.unique(rfd3_array.res_id[b_chain_mask])
                fixed_residues_str.extend([f"{fixed_chains[0]}{pos}" for pos in b_res_ids])

                mpnn_input_dict = {
                    "name": design_id,
                    "batch_size": 1, # Number of sequences to sample per structure
                    "remove_waters": True,
                    "seed": 42,
                    "fixed_residues": fixed_residues_str,
                    "sampling_temp": 0.1
                }
                
                mpnn_outputs = lmpnn_engine.run(input_dicts=[mpnn_input_dict], atom_arrays=[rfd3_array])
                
                for seq_idx, mpnn_out in enumerate(mpnn_outputs):
                    valid_mask = ~np.isnan(mpnn_out.atom_array.coord[:, 0])
                    lmpnn_array = mpnn_out.atom_array[valid_mask]
                    lmpnn_array = renumber_atom_array_residues(lmpnn_array)
                    
                    full_seq_designed = get_sequence_from_array(lmpnn_array, chain_to_design)
                    
                    to_cif_file(lmpnn_array, f"{lmpnn_out_dir}/{design_id}_mpnn{seq_idx}.cif", file_type="cif")
                    
                    # Compute loops for metadata extraction
                    loop_data = {}
                    curr_idx = 0
                    
                    for name_idx, loop in enumerate(selected_loops):
                        loop_name = loop_names[name_idx]
                        f_left = loop["left"]
                        f_right = loop["right"]
                        m = re.search(f"{f_left}(.*?){f_right}", full_seq_designed[curr_idx:])
                        
                        if m:
                            seq = m.group(1)
                            loop_data[f"loop_{loop_name}_seq"] = seq
                            loop_data[f"loop_{loop_name}_length"] = len(seq)
                            curr_idx += m.end() - len(f_right)
                        else:
                            loop_data[f"loop_{loop_name}_seq"] = "MISSING"
                            loop_data[f"loop_{loop_name}_length"] = 0
                    
                    final_records.append({
                        "design_id": design_id,
                        "seq_idx": seq_idx,
                        **loop_data,
                        "seq_recovery": float(getattr(mpnn_out, "output_dict", {}).get("sequence_recovery", 0.0)),
                        "full_seq": full_seq_designed,
                        "lmpnn_array": lmpnn_array,
                        "rfd3_array": rfd3_array
                    })

            except Exception as e:
                print(f"Error during LigandMPNN sequence generation for {design_id}: {e}")
                
        end = time.time()
        print(f"LigandMPNN finished in {(end-st)/60:.2f} minutes.")
        
        del lmpnn_engine
        torch.cuda.empty_cache()

        # --- STAGE 3: RUN RF3 VALIDATION ---
        print(f"\n--- Running RF3 Validations on {len(final_records)} sequences ---")
        rf3_engine = RF3InferenceEngine(ckpt_path='rf3', verbose=False)
        metric_records = []
        
        st = time.time()
        
        for job in final_records:
            design_id = f"{job['design_id']}_mpnn{job['seq_idx']}"
            print(f"Running RF3 validation for {design_id}...")
            try:
                # Prepare input array for structural prediction
                lmpnn_array = job['lmpnn_array']
                rfd3_array = job['rfd3_array']
                
                # Standard atoms for folding
                valid_atoms = ['N', 'CA', 'C', 'O', 'CB']
                array_for_rf3 = lmpnn_array[np.isin(lmpnn_array.atom_name, valid_atoms)].copy()
                
                input_structure = InferenceInput.from_atom_array(array_for_rf3, example_id=design_id)
                rf3_outputs_dict = rf3_engine.run(inputs=input_structure)
                
                # Grab the first available dictionary output key if specific target fails
                rf3_target_key = next((k for k in rf3_outputs_dict.keys() if design_id in k), list(rf3_outputs_dict.keys())[0] if rf3_outputs_dict else None)

                if rf3_target_key and rf3_target_key in rf3_outputs_dict:
                    rf3_output = rf3_outputs_dict[rf3_target_key][0]
                    rf3_atom_array = rf3_output.atom_array
                    
                    rf3_out_dir = os.path.join(original_directory, "../Local/rf3_output", "TIMP3_vs_ADAM17_HADDOCK_Xray")
                    os.makedirs(rf3_out_dir, exist_ok=True)
                    
                    rf3_atom_array = renumber_atom_array_residues(rf3_atom_array)
                    to_cif_file(rf3_atom_array, f"{rf3_out_dir}/{design_id}_refolded.cif", file_type="cif")

                    # CA atom superimposition for evaluation
                    bb_mask_rfd3 = np.isin(rfd3_array.atom_name, PROTEIN_BACKBONE_ATOM_NAMES)
                    bb_mask_rf3 = np.isin(rf3_atom_array.atom_name, PROTEIN_BACKBONE_ATOM_NAMES)
                    
                    bb_generated = rfd3_array[bb_mask_rfd3]
                    bb_refolded = rf3_atom_array[bb_mask_rf3]
                    
                    bb_generated = bb_generated[bb_generated.atom_name != "OXT"]
                    bb_refolded = bb_refolded[bb_refolded.atom_name != "OXT"]
                    
                    # Truncate ends if length mismatch
                    if len(bb_generated) != len(bb_refolded):
                        min_len = min(len(bb_generated), len(bb_refolded))
                        bb_generated = bb_generated[:min_len]
                        bb_refolded = bb_refolded[:min_len]

                    # Structural metrics
                    bb_refolded_fitted, _ = superimpose(bb_generated, bb_refolded)
                    overall_rmsd = rmsd(bb_generated, bb_refolded_fitted)
                    
                    summary = rf3_output.summary_confidences
                    plddt = summary.get('overall_plddt', 0.0)
                    ptm = summary.get('ptm', 0.0)
                    
                    # Dynamically collect passed loop configurations
                    loop_kwargs = {k: v for k, v in job.items() if k.startswith("loop_")}
                    
                    metric_records.append({
                        "file": f"{design_id}.cif",
                        **loop_kwargs,
                        "overall_rmsd": overall_rmsd,
                        "plddt": plddt,
                        "ptm": ptm,
                        "seq_recovery": job["seq_recovery"],
                        "full_seq": job["full_seq"]
                    })
                    
            except Exception as e:
                print(f"Error during RF3 validation for {design_id}: {e}")
                
        end = time.time()
        print(f"RF3 finished in {(end-st)/60:.2f} minutes.")
        
        del rf3_engine
        torch.cuda.empty_cache()

        # --- SUMMARIZE RESULTS ---
        print("\n--- Saving Summaries ---")
        if metric_records:
            df = pd.DataFrame(metric_records)
            df.to_csv(os.path.join(lmpnn_out_dir, "design_summary.csv"), index=False)
            
            # Find all dynamic loop sequence columns
            seq_cols = [c for c in df.columns if c.startswith("loop_") and c.endswith("_seq")]
            len_cols = [c for c in df.columns if c.startswith("loop_") and c.endswith("_length")]
            
            clean_df = df.dropna(subset=seq_cols + ["plddt"]).copy()
            if not clean_df.empty:
                # Rank primarily by RF3 Folding confidence (pLDDT is 0-100 scale, higher is better)
                best_per_loop = clean_df.sort_values("plddt", ascending=False).groupby(seq_cols, as_index=False).first()
                
                avg_stats = clean_df.groupby(len_cols, as_index=False).agg(
                    avg_plddt=("plddt", "mean"),
                    avg_rmsd=("overall_rmsd", "mean"),
                    count=("plddt", "count")
                ).sort_values(len_cols)
                
                best_per_length = best_per_loop.sort_values(len_cols + ["plddt"], ascending=[True]*len(len_cols) + [False]).groupby(len_cols, group_keys=False).head(5)
                best_per_length = best_per_length.merge(avg_stats, on=len_cols, how="left")
                best_per_length = best_per_length.sort_values(len_cols + ["plddt"], ascending=[True]*len(len_cols) + [False])
                
                best_per_length.to_csv(os.path.join(lmpnn_out_dir, "best_loops_per_length.csv"), index=False)
                avg_stats.to_csv(os.path.join(lmpnn_out_dir, "loop_length_averages.csv"), index=False)
                
                print(f"Saved summaries to {lmpnn_out_dir}")
        else:
            print("No valid metric records generated.")

if __name__ == "__main__":
    main()
