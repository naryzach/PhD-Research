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
    parser.add_argument("--ions", nargs='+', required=True, help="List of target metal ions (e.g. EU ZN TB DY)")
    parser.add_argument("--num-designs", type=int, default=1, help="Number of sequences to generate per loop")
    parser.add_argument("--input", type=str, default="../Data/8FNS.cif", help="Path to input CIF structure")
    parser.add_argument("--test", action="store_true", help="Run a short test by processing only the first loop")
    parser.add_argument("--test-steps", type=int, default=10, help="Number of diffusion steps for RFd3 when in test mode (default: 10).")
    return parser.parse_args()

def main():
    args = parse_args()
    seed_everything(42)
    
    # Point Foundry to the correct checkpoint dir relative to this script
    checkpoint_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Tools", "foundry_checkpoints"))
    os.makedirs(checkpoint_dir, exist_ok=True)
    os.environ["FOUNDRY_CHECKPOINT_DIRS"] = checkpoint_dir
    
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
        os.makedirs(f"{ion_dir}/rf3", exist_ok=True)
        
        # Mutate to current ion in the loop
        atom_array = mutate_metals(base_atom_array.copy(), original_res_name="ND", new_res_name=ion)
        loops = get_ef_hand_loops(atom_array, metal_res_name=ion)
        
        print(f"\nFound {len(loops)} metal binding sites for {ion}.")
        if args.test: 
            print("Test mode enabled. Only processing the first metal binding site.")
            loops = loops[:1]
            
        # We don't need to delete the loop residues! RFd3 ignores them if they fall into the generated gaps in the contig.
        # But we must remove waters or other unwanted heterogens just in case.
        water_mask = atom_array.res_name == "HOH"
        context_array = atom_array[~water_mask]
        
        chain_id = "A"
        ca_atoms = context_array[(context_array.atom_name == "CA") & (context_array.element == "C") & (context_array.chain_id == chain_id)]
        min_res = int(np.min(ca_atoms.res_id))
        max_res = int(np.max(ca_atoms.res_id))
        total_len = max_res - min_res + 1 
        
        # Build single multi-loop contig string for all loops
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
            # RFd3 clusters input ligands chronologically to the end of the built chain
            new_loop['metal_res_id'] = new_total_len + new_loop['metal_idx'] + 1
            new_loop['metal_chain_id'] = chain_id
        
        loop_cif_path = f"outputs/{ion}/rfd3/context_all_loops.cif"
        to_cif_file(context_array, loop_cif_path)
        
        spec_input = DesignInputSpecification(
            input=loop_cif_path,
            contig=contig_str, 
            length=f"{total_len}-{total_len}",
            ligand=ion,
            select_buried={ion: "ALL"}, # Bury all instances of the ion
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
    
    # Wrap rfd3_queue in tqdm
    for job in tqdm(rfd3_queue, desc="RFd3 Generation", unit="design"):
        print(f"Generating simultaneous backbones for {job['ion']}...", flush=True)
        print(f"Generated contig string: {job['contig_str']}", flush=True)
        
        try:
            inference_sampler_kwargs = {}
            if args.test:
                # RFd3 requires at least 2 timesteps to run the diffusion loop properly
                test_steps = max(2, args.test_steps)
                inference_sampler_kwargs['num_timesteps'] = test_steps
                print(f"Test mode active: setting RFd3 diffusion steps to {test_steps}")
                
            rfd3_config = getattr(__import__('rfd3.engine', fromlist=['RFD3InferenceConfig']), 'RFD3InferenceConfig')(
                diffusion_batch_size=args.num_designs,
                low_memory_mode=False,
                specification={'length': job['total_len'], 'contig': job['contig_str'], 'extra': {}},
                inference_sampler=inference_sampler_kwargs
            )
            rfd3_engine = RFD3InferenceEngine(
                #verbose=True,
                **dataclasses.asdict(rfd3_config)
            )
            rfd3_outputs_dict = rfd3_engine.run(inputs=job['spec'], n_batches=1, out_dir=None)
            del rfd3_engine
            
            if rfd3_outputs_dict:
                design_count = 0
                for key, rfd3_out_list in rfd3_outputs_dict.items():
                    if not key.startswith("backbone"): continue
                    
                    for rfd3_out in rfd3_out_list:
                        design_id = f"{job['ion']}_all_loops_design{design_count}"
                        to_cif_file(rfd3_out.atom_array, f"outputs/{job['ion']}/rfd3/{design_id}.cif", file_type="cif")
                        
                        mpnn_queue.append({
                            **job, 
                            'design_id': design_id,
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
    
    mpnn_engine = MPNNInferenceEngine(model_type="ligand_mpnn", is_legacy_weights=True, write_structures=False, write_fasta=False)
    
    # Wrap mpnn_queue in tqdm
    for job in tqdm(mpnn_queue, desc="MPNN Sequences", unit="seq"):
        print(f"Running MPNN for {job['design_id']}...", flush=True)
        
        try:
            # Fix all positions except the combined redesign loop residues
            fixed_res_set = set(range(job['min_res'], job['max_res'] + 1)) - set(job['redesign_res_ids'])
            fixed_res_list = sorted(list(fixed_res_set))
            
            fixed_positions = {job['chain_id']: fixed_res_list}
            mpnn_input_configs = [{"batch_size": 1, "remove_waters": True, "fixed_positions_dict": fixed_positions}]
            
            mpnn_outputs = mpnn_engine.run(input_dicts=mpnn_input_configs, atom_arrays=[job['rfd3_atom_array']])
            
            if mpnn_outputs:
                mpnn_out = mpnn_outputs[0]
                mpnn_res_starts = get_residues(mpnn_out.atom_array)[0]
                mpnn_res_names = []
                
                for rs in mpnn_res_starts:
                    residue_id = mpnn_out.atom_array.res_id[rs]
                    # Only map redesigned residues, and ensure we only count one atom (CA) per residue
                    if residue_id in job['redesign_res_ids'] and mpnn_out.atom_array.atom_name[rs] == 'CA':
                        mpnn_res_names.append(mpnn_out.atom_array.res_name[rs])
                
                seq_letters = []
                for rn in mpnn_res_names:
                    try:
                        seq_letters.append(ProteinSequence.convert_letter_3to1(rn))
                    except:
                        seq_letters.append("X")
                gen_seq = "".join(seq_letters) # This represents the concatenation of all designed loop sequences
                
                # Save LigandMPNN output independently
                lmpnn_out_dir = f"outputs/{job['ion']}/lmpnn"
                os.makedirs(lmpnn_out_dir, exist_ok=True)
                
                # Filter out NaN coordinates (like terminal OXT) before saving to prevent PyMOL parser errors
                valid_mask = ~np.isnan(mpnn_out.atom_array.coord[:, 0])
                clean_mpnn_array = mpnn_out.atom_array[valid_mask]
                
                to_cif_file(clean_mpnn_array, f"{lmpnn_out_dir}/{job['design_id']}_lmpnn.cif", file_type="cif")
                
                rf3_queue.append({
                    **job,
                    'gen_seq': gen_seq,
                    'mpnn_atom_array': mpnn_out.atom_array
                })
        except Exception as e:
            print(f"Error during MPNN for {job['design_id']}: {e}")

    del mpnn_engine
    gc.collect()
    torch.cuda.empty_cache()

    # ---------------------------------------------------------
    # STAGE 3: Validation & Scoring (RF3)
    # ---------------------------------------------------------
    print(f"\n--- STAGE 3: Running {len(rf3_queue)} RF3 Validations ---", flush=True)
    rf3_engine = RF3InferenceEngine(ckpt_path='rf3', verbose=False)
    final_records = []
    
    for job in rf3_queue:
        print(f"Running RF3 validation for {job['design_id']}...", flush=True)
        
        try:
            input_structure = InferenceInput.from_atom_array(job['mpnn_atom_array'], example_id=job['design_id'])
            rf3_outputs_dict = rf3_engine.run(inputs=input_structure)
            
            # RF3 output keys might have suffixes. Find the matching key.
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
                
                # Save RF3 specific output
                rf3_out_dir = f"outputs/{job['ion']}/rf3"
                os.makedirs(rf3_out_dir, exist_ok=True)
                to_cif_file(rf3_atom_array, f"{rf3_out_dir}/{job['design_id']}_refolded.cif", file_type="cif")
                
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
                
                loop_mask = np.isin(bb_generated.res_id, job['redesign_res_ids'])
                loop_bb_generated = bb_generated[loop_mask]
                loop_bb_refolded_fitted = bb_refolded_fitted[loop_mask]
                
                loop_rmsd = rmsd(loop_bb_generated, loop_bb_refolded_fitted) if len(loop_bb_generated) > 0 else float('nan')
                
                summary = rf3_output.summary_confidences
                plddt = summary['overall_plddt']
                
                # Average binding radius for all loops
                radii = []
                for loop_info in job['loops']:
                    rad = calculate_binding_radius(rf3_atom_array, loop_info, loop_info['start_res'], loop_info['end_res'])
                    if not np.isnan(rad): radii.append(rad)
                
                avg_binding_radius = np.mean(radii) if radii else float('nan')
                
                print(f">>> {job['design_id']} | Seq: {job['gen_seq']} | Overall RMSD: {overall_rmsd:.2f} | Avg Loop RMSD: {loop_rmsd:.2f} | Avg Radius: {avg_binding_radius:.2f} A", flush=True)
                
                final_records.append({
                    "metal_ion": job['ion'],
                    "design_id": job['design_id'],
                    "generated_sequence_concat": job['gen_seq'],
                    "overall_rmsd": overall_rmsd,
                    "avg_loop_rmsd": loop_rmsd,
                    "plddt": plddt,
                    "ptm": summary.get('ptm', 0),
                    "avg_binding_radius_A": avg_binding_radius
                })
        except Exception as e:
            print(f"Error during RF3 validation for {job['design_id']}: {e}")

    del rf3_engine
    gc.collect()
    torch.cuda.empty_cache()

    # ---------------------------------------------------------
    # STAGE 4: Save Global Catalog
    # ---------------------------------------------------------
    if final_records:
        df = pd.DataFrame(final_records)
        df.to_csv("outputs/global_sequence_catalog.csv", index=False)
        print("\nBatch processing complete. Catalog saved to outputs/global_sequence_catalog.csv")
    else:
        print("\nPipeline finished, but no valid records were generated.")

if __name__ == "__main__":
    main()