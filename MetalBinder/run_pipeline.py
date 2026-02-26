import os
import sys
import numpy as np
import pandas as pd
import argparse
import dataclasses
import gc
import torch
from tqdm import tqdm
from biotite.structure.io.pdbx import CIFFile, get_structure as get_cif_structure
from biotite.structure import AtomArray, get_residues, array
from atomworks.io.utils.io_utils import to_cif_file
from biotite.sequence import ProteinSequence

# Foundry Inference imports
from lightning.fabric import seed_everything
from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine
from mpnn.inference_engines.mpnn import MPNNInferenceEngine
from rf3.inference_engines.rf3 import RF3InferenceEngine
from rf3.utils.inference import InferenceInput
from biotite.structure import rmsd, superimpose
from atomworks.constants import PROTEIN_BACKBONE_ATOM_NAMES
from rfd3.inference.input_parsing import DesignInputSpecification

# Custom Utilities
from utils_foundry import get_ef_hand_loops, create_masked_input, mutate_metals, calculate_binding_radius

torch.set_float32_matmul_precision('medium')

def parse_args():
    parser = argparse.ArgumentParser(description="Batch RFd3 -> LMPNN -> RF3 binding loop generation.")
    parser.add_argument("--ions", nargs='+', required=True, help="List of target metal ions (e.g. EU ZN TB DY CA)")
    parser.add_argument("--num-designs", type=int, default=1, help="Number of sequences to generate per loop")
    parser.add_argument("--input", type=str, default="../Data/8FNS.cif", help="Path to input CIF structure")
    parser.add_argument("--test", action="store_true", help="Run a short test by processing only the first loop")
    parser.add_argument("--test-steps", type=int, default=10, help="Number of diffusion steps for RFd3 when in test mode (default: 10).")
    return parser.parse_args()

def get_sequence_from_array(atom_array, chain_id="A", res_ids=None):
    """Helper function to extract sequence cleanly from an AtomArray."""
    mask = (atom_array.chain_id == chain_id) & (atom_array.atom_name == "CA")
    if res_ids is not None:
        mask = mask & np.isin(atom_array.res_id, res_ids)
        
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

def main():
    args = parse_args()
    seed_everything(42)
    
    checkpoint_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Tools", "foundry_checkpoints"))
    os.makedirs(checkpoint_dir, exist_ok=True)
    os.environ["FOUNDRY_CHECKPOINT_DIRS"] = checkpoint_dir
    
    full_sequence_records = []
    
    # ---------------------------------------------------------
    # STAGE 0: Preparation & Queuing
    # ---------------------------------------------------------
    print("--- STAGE 0: Preparing Inputs ---", flush=True)
    cif_file = CIFFile.read(args.input)
    base_atom_array = get_cif_structure(cif_file, model=1)
    
    rfd3_queue = []
    
    for ion in args.ions:
        ion_dir = f"outputs/{ion}"
        os.makedirs(f"{ion_dir}/rfd3", exist_ok=True)
        os.makedirs(f"{ion_dir}/lmpnn", exist_ok=True)
        os.makedirs(f"{ion_dir}/rf3", exist_ok=True)
        
        atom_array = mutate_metals(base_atom_array.copy(), original_res_name="ND", new_res_name=ion)
        loops = get_ef_hand_loops(atom_array, metal_res_name=ion)
        
        print(f"\nFound {len(loops)} metal binding sites for {ion}.")
        if args.test: 
            print("Test mode enabled. Only processing the first metal binding site.")
            loops = loops[:1]
            
        water_mask = atom_array.res_name == "HOH"
        context_array = atom_array[~water_mask]
        
        chain_id = "A"
        ca_atoms = context_array[(context_array.atom_name == "CA") & (context_array.element == "C") & (context_array.chain_id == chain_id)]
        min_res = int(np.min(ca_atoms.res_id))
        max_res = int(np.max(ca_atoms.res_id))
        total_len = max_res - min_res + 1 
        
        contig_parts = []
        curr_res = min_res
        redesign_res_ids = []
        
        sorted_loops = sorted(loops, key=lambda x: x["start_res"])
        new_res_counter = 1
        new_loops = []
        
        for loop_info in sorted_loops:
            start_res = loop_info["start_res"]
            end_res = loop_info["end_res"]
            
            if start_res > curr_res:
                segment_len = start_res - curr_res
                contig_parts.append(f"{chain_id}{curr_res}-{start_res-1}")
                new_res_counter += segment_len
            
            loop_length = end_res - start_res + 1
            contig_parts.append(str(loop_length))
            
            new_loop = loop_info.copy()
            new_loop['start_res'] = new_res_counter
            new_loop['end_res'] = new_res_counter + loop_length - 1
            new_loops.append(new_loop)
            
            for r in range(new_loop['start_res'], new_loop['end_res'] + 1):
                redesign_res_ids.append(r)
                
            new_res_counter += loop_length
            curr_res = end_res + 1
                
        if curr_res <= max_res:
            contig_parts.append(f"{chain_id}{curr_res}-{max_res}")
            new_res_counter += (max_res - curr_res + 1)
            
        contig_str = ",".join(contig_parts)
        new_total_len = new_res_counter - 1
        
        for new_loop in new_loops:
            new_loop['metal_res_id'] = new_total_len + new_loop['metal_idx'] + 1
            new_loop['metal_chain_id'] = chain_id
        
        loop_cif_path = f"outputs/{ion}/rfd3/context_all_loops.cif"
        to_cif_file(context_array, loop_cif_path)
        
        # LOG THE STARTING CONTEXT SEQUENCE
        context_seq = get_sequence_from_array(context_array, chain_id=chain_id)
        full_sequence_records.append({
            'ion': ion,
            'phase': 'Context',
            'design_id': f"{ion}_wildtype_context",
            'full_sequence': context_seq
        })
        
        spec_input = DesignInputSpecification(
            input=loop_cif_path,
            contig=contig_str, 
            length=f"{total_len}-{total_len}",
            ligand=ion,
            select_buried={ion: "ALL"},
            extra={}
        )
        
        rfd3_queue.append({
            'ion': ion,
            'loops': new_loops,
            'spec': spec_input,
            'total_len': total_len,
            'contig_str': contig_str,
            'redesign_res_ids': redesign_res_ids,
            'min_res': 1,
            'max_res': new_total_len,
            'chain_id': chain_id
        })

    # ---------------------------------------------------------
    # STAGE 1: Backbone Generation (RFd3)
    # ---------------------------------------------------------
    print(f"\n--- STAGE 1: Running {len(rfd3_queue)} RFd3 Jobs ---", flush=True)
    mpnn_queue = []
    
    for job in tqdm(rfd3_queue, desc="RFd3 Generation", unit="design"):
        print(f"Generating simultaneous backbones for {job['ion']}...", flush=True)
        print(f"Generated contig string: {job['contig_str']}", flush=True)
        try:
            inference_sampler_kwargs = {}
            if args.test:
                test_steps = max(2, args.test_steps)
                inference_sampler_kwargs['num_timesteps'] = test_steps
                
            rfd3_config = getattr(__import__('rfd3.engine', fromlist=['RFD3InferenceConfig']), 'RFD3InferenceConfig')(
                diffusion_batch_size=args.num_designs,
                low_memory_mode=False,
                specification={'length': job['total_len'], 'contig': job['contig_str'], 'extra': {}},
                inference_sampler=inference_sampler_kwargs
            )
            rfd3_engine = RFD3InferenceEngine(**dataclasses.asdict(rfd3_config))
            rfd3_outputs_dict = rfd3_engine.run(inputs=job['spec'], n_batches=1, out_dir=None)
            del rfd3_engine
            
            if rfd3_outputs_dict:
                design_count = 0
                for key, rfd3_out_list in rfd3_outputs_dict.items():
                    if not key.startswith("backbone"): continue
                    
                    for rfd3_out in rfd3_out_list:
                        design_id = f"{job['ion']}_all_loops_design{design_count}"
                        to_cif_file(rfd3_out.atom_array, f"outputs/{job['ion']}/rfd3/{design_id}.cif", file_type="cif")
                        
                        full_sequence_records.append({
                            'ion': job['ion'],
                            'phase': 'RFd3',
                            'design_id': design_id,
                            'full_sequence': get_sequence_from_array(rfd3_out.atom_array, chain_id=job['chain_id'])
                        })
                        
                        mpnn_queue.append({
                            **job, 
                            'design_id': design_id,
                            'design_number': design_count,
                            'rfd3_atom_array': rfd3_out.atom_array
                        })
                        design_count += 1
        except Exception as e:
            print(f"Error during RFd3 generation for {job['ion']}: {e}")

    gc.collect()
    torch.cuda.empty_cache()

    # ---------------------------------------------------------
    # STAGE 2: Sequence Design (LigandMPNN)
    # ---------------------------------------------------------
    print(f"\n--- STAGE 2: Running {len(mpnn_queue)} LigandMPNN Jobs ---", flush=True)
    rf3_queue = []
    
    lmpnn_engine = MPNNInferenceEngine(model_type="ligand_mpnn", is_legacy_weights=True, write_structures=False, write_fasta=False)
    
    for job in tqdm(mpnn_queue, desc="LigandMPNN", unit="seq"):
        print(f"Running MPNN for {job['design_id']}...", flush=True)
        try:
            fixed_res_set = set(range(job['min_res'], job['max_res'] + 1)) - set(job['redesign_res_ids'])
            fixed_positions = {job['chain_id']: sorted(list(fixed_res_set))}
            mpnn_input_configs = [{"batch_size": 1, "remove_waters": True, "fixed_positions_dict": fixed_positions}]
            
            mpnn_outputs = lmpnn_engine.run(input_dicts=mpnn_input_configs, atom_arrays=[job['rfd3_atom_array']])
            
            if mpnn_outputs:
                mpnn_out = mpnn_outputs[0]
                valid_mask = ~np.isnan(mpnn_out.atom_array.coord[:, 0])
                clean_mpnn_array = mpnn_out.atom_array[valid_mask]
                
                lmpnn_design_id = f"{job['design_id']}_LMPNN"
                to_cif_file(clean_mpnn_array, f"outputs/{job['ion']}/lmpnn/{lmpnn_design_id}.cif", file_type="cif")
                
                full_sequence_records.append({
                    'ion': job['ion'],
                    'phase': 'LMPNN',
                    'design_id': lmpnn_design_id,
                    'full_sequence': get_sequence_from_array(clean_mpnn_array, chain_id=job['chain_id'])
                })
                
                rf3_queue.append({
                    **job,
                    'mpnn_atom_array': clean_mpnn_array
                })
        except Exception as e:
            print(f"Error during LMPNN for {job['design_id']}: {e}")

    del lmpnn_engine
    gc.collect()
    torch.cuda.empty_cache()

    # ---------------------------------------------------------
    # STAGE 3: Validation & Scoring (RF3)
    # ---------------------------------------------------------
    print(f"\n--- STAGE 3: Running {len(rf3_queue)} RF3 Validations ---", flush=True)
    rf3_engine = RF3InferenceEngine(ckpt_path='rf3', verbose=False)
    final_records = []
    
    for job in tqdm(rf3_queue, desc="RF3 Validations", unit="fold"):
        print(f"Running RF3 validation for {job['design_id']}...", flush=True)
        try:
            input_structure = InferenceInput.from_atom_array(job['mpnn_atom_array'], example_id=job['design_id'])
            rf3_outputs_dict = rf3_engine.run(inputs=input_structure)
            
            rf3_target_key = None
            for key in rf3_outputs_dict.keys():
                if job['design_id'] in key or key in job['design_id'] or key == job['design_id']:
                    rf3_target_key = key
                    break
            
            if not rf3_target_key and len(rf3_outputs_dict) > 0:
                rf3_target_key = list(rf3_outputs_dict.keys())[0]
                
            if rf3_target_key in rf3_outputs_dict:
                rf3_output = rf3_outputs_dict[rf3_target_key][0]
                rf3_atom_array = rf3_output.atom_array
                
                rf3_out_dir = f"outputs/{job['ion']}/rf3"
                os.makedirs(rf3_out_dir, exist_ok=True)
                to_cif_file(rf3_atom_array, f"{rf3_out_dir}/{job['design_id']}_refolded.cif", file_type="cif")
                
                # Extract RF3 refolded sequences
                full_rf3_seq = get_sequence_from_array(rf3_atom_array, chain_id=job['chain_id'])
                
                full_sequence_records.append({
                    'ion': job['ion'],
                    'phase': 'RF3',
                    'design_id': job['design_id'],
                    'full_sequence': full_rf3_seq
                })
                
                # --- Structure Verification & RMSD ---
                bb_mask_rfd3 = np.isin(job['rfd3_atom_array'].atom_name, PROTEIN_BACKBONE_ATOM_NAMES)
                bb_mask_rf3 = np.isin(rf3_atom_array.atom_name, PROTEIN_BACKBONE_ATOM_NAMES)
                
                bb_generated = job['rfd3_atom_array'][bb_mask_rfd3]
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
                plddt = summary['overall_plddt']
                
                # Iterate through each specific loop to record individual rows
                for loop_idx, loop_info in enumerate(job['loops']):
                    loop_start = loop_info['start_res']
                    loop_end = loop_info['end_res']
                    
                    # Individual Loop RMSD
                    loop_mask = (bb_generated.res_id >= loop_start) & (bb_generated.res_id <= loop_end)
                    loop_bb_generated = bb_generated[loop_mask]
                    loop_bb_refolded_fitted = bb_refolded_fitted[loop_mask]
                    
                    loop_rmsd = rmsd(loop_bb_generated, loop_bb_refolded_fitted) if len(loop_bb_generated) > 0 else float('nan')
                    
                    # Individual Binding Radius
                    rad = calculate_binding_radius(rf3_atom_array, loop_info, loop_start, loop_end)
                    
                    # Individual Loop Sequence
                    loop_res_ids = list(range(loop_start, loop_end + 1))
                    loop_seq = get_sequence_from_array(rf3_atom_array, chain_id=job['chain_id'], res_ids=loop_res_ids)
                    
                    final_records.append({
                        "metal_ion": job['ion'],
                        "design_id": job['design_id'],
                        "design_number": job['design_number'],
                        "loop_index": loop_idx + 1,
                        "loop_sequence": loop_seq,
                        "overall_rmsd": overall_rmsd,
                        "loop_rmsd": loop_rmsd,
                        "plddt": plddt,
                        "ptm": summary.get('ptm', 0),
                        "binding_radius_A": rad
                    })

        except Exception as e:
            print(f"Error during RF3 validation for {job['design_id']}: {e}")

    del rf3_engine
    gc.collect()
    torch.cuda.empty_cache()

    # ---------------------------------------------------------
    # STAGE 4: Save Catalogs
    # ---------------------------------------------------------
    if full_sequence_records:
        seq_df = pd.DataFrame(full_sequence_records)
        seq_df.to_csv("outputs/full_sequences_log.csv", index=False)
        print("\nFull sequences logged to outputs/full_sequences_log.csv")

    if final_records:
        df = pd.DataFrame(final_records)
        df.to_csv("outputs/global_sequence_catalog.csv", index=False)
        print("Batch processing complete. Metrics saved to outputs/global_sequence_catalog.csv")
    else:
        print("\nPipeline finished, but no valid metric records were generated.")

if __name__ == "__main__":
    main()