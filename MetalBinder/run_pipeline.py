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
            
        for loop_info in loops:
            context_array, start_res, end_res = create_masked_input(atom_array, loop_info)
            chain_id = loop_info.get("protein_chain_id", "A") 
            
            # Extract true protein boundaries using Alpha Carbons (filtering out Calcium ions)
            context_chain = context_array[context_array.chain_id == chain_id]
            ca_atoms = context_chain[
                (context_chain.atom_name == "CA") & 
                (context_chain.element == "C")
            ]
            
            min_res = int(np.min(ca_atoms.res_id))
            max_res = int(np.max(ca_atoms.res_id))
            total_len = max_res - min_res + 1 
            
            # Build the contig string
            contig_parts = []
            if start_res > min_res:
                contig_parts.append(f"{chain_id}{min_res}-{start_res-1}")
            
            loop_length = end_res - start_res + 1
            contig_parts.append(str(loop_length))
            
            if end_res < max_res:
                contig_parts.append(f"{chain_id}{end_res+1}-{max_res}")
                
            contig_str = ",".join(contig_parts)
            
            loop_cif_path = f"outputs/{ion}/rfd3/context_{loop_info['metal_idx']}.cif"
            to_cif_file(context_array, loop_cif_path)
            
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
                'loop_info': loop_info,
                'spec': spec_input,
                'total_len': total_len,
                'contig_str': contig_str,
                'start_res': start_res,
                'end_res': end_res,
                'min_res': min_res,
                'max_res': max_res,
                'chain_id': chain_id
            })

    # ---------------------------------------------------------
    # STAGE 1: Backbone Generation (RFd3)
    # ---------------------------------------------------------
    print(f"\n--- STAGE 1: Running {len(rfd3_queue)} RFd3 Jobs ---", flush=True)
    mpnn_queue = []
    
    # Wrap rfd3_queue in tqdm
    for job in tqdm(rfd3_queue, desc="RFd3 Generation", unit="loop"):
        print(f"Generating backbones for {job['ion']} Loop {job['loop_info']['metal_idx']}...", flush=True)
        
        try:
            rfd3_config = getattr(__import__('rfd3.engine', fromlist=['RFD3InferenceConfig']), 'RFD3InferenceConfig')(
                diffusion_batch_size=args.num_designs,
                low_memory_mode=True,
                specification={'length': job['total_len'], 'contig': job['contig_str'], 'extra': {}}
            )
            rfd3_engine = RFD3InferenceEngine(
                #verbose=True,
                **dataclasses.asdict(rfd3_config)
            )
            rfd3_outputs_dict = rfd3_engine.run(inputs=job['spec'], n_batches=1, out_dir=None)
            del rfd3_engine
            
            if rfd3_outputs_dict:
                design_count = 0
                for key, rfd3_out in rfd3_outputs_dict.items():
                    if not key.startswith("backbone"): continue
                    
                    design_id = f"{job['ion']}_loop_{job['loop_info']['metal_idx']}_design{design_count}"
                    to_cif_file(rfd3_out.atom_array, f"outputs/{job['ion']}/rfd3/{design_id}.cif", file_type="cif")
                    
                    mpnn_queue.append({
                        **job, 
                        'design_id': design_id,
                        'rfd3_atom_array': rfd3_out.atom_array
                    })
                    design_count += 1
        except Exception as e:
            print(f"Error during RFd3 generation for {job['ion']} Loop {job['loop_info']['metal_idx']}: {e}")

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
            fixed_res_list = []
            for r in range(job['min_res'], job['max_res'] + 1):
                if r < job['start_res'] or r > job['end_res']:
                    fixed_res_list.append(r)
            
            fixed_positions = {job['chain_id']: fixed_res_list}
            mpnn_input_configs = [{"batch_size": 1, "remove_waters": True, "fixed_positions_dict": fixed_positions}]
            
            mpnn_outputs = mpnn_engine.run(input_dicts=mpnn_input_configs, atom_arrays=[job['rfd3_atom_array']])
            
            if mpnn_outputs:
                mpnn_out = mpnn_outputs[0]
                mpnn_res_starts = get_residues(mpnn_out.atom_array)[0]
                mpnn_res_names = []
                
                for rs in mpnn_res_starts:
                    residue_id = mpnn_out.atom_array.res_id[rs]
                    if job['start_res'] <= residue_id <= job['end_res']:
                        mpnn_res_names.append(mpnn_out.atom_array.res_name[rs])
                
                seq_letters = [ProteinSequence.convert_letter_3to1(rn) if rn in ProteinSequence._letters_3 else "X" for rn in mpnn_res_names]
                gen_seq = "".join(seq_letters)
                
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
            
            if job['design_id'] in rf3_outputs_dict:
                rf3_output = rf3_outputs_dict[job['design_id']][0]
                rf3_atom_array = rf3_output.atom_array
                to_cif_file(rf3_atom_array, f"outputs/{job['ion']}/rf3/{job['design_id']}_refolded.cif", file_type="cif")
                
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
                
                loop_mask = (bb_generated.res_id >= job['start_res']) & (bb_generated.res_id <= job['end_res'])
                loop_bb_generated = bb_generated[loop_mask]
                loop_bb_refolded_fitted = bb_refolded_fitted[loop_mask]
                
                loop_rmsd = rmsd(loop_bb_generated, loop_bb_refolded_fitted) if len(loop_bb_generated) > 0 else float('nan')
                
                summary = rf3_output.summary_confidences
                plddt = summary['overall_plddt']
                
                binding_radius = calculate_binding_radius(rf3_atom_array, job['loop_info'], job['start_res'], job['end_res'])
                
                print(f">>> {job['design_id']} | Seq: {job['gen_seq']} | Overall RMSD: {overall_rmsd:.2f} | Loop RMSD: {loop_rmsd:.2f} | Radius: {binding_radius:.2f} A", flush=True)
                
                final_records.append({
                    "metal_ion": job['ion'],
                    "target_loop_idx": job['loop_info']['metal_idx'],
                    "design_id": job['design_id'],
                    "generated_sequence": job['gen_seq'],
                    "overall_rmsd": overall_rmsd,
                    "loop_rmsd": loop_rmsd,
                    "plddt": plddt,
                    "ptm": summary.get('ptm', 0),
                    "binding_radius_A": binding_radius
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