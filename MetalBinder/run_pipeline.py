import os
import sys
import numpy as np
import pandas as pd
from biotite.structure.io.pdbx import PDBxFile, get_structure as get_cif_structure
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

from utils_foundry import get_ef_hand_loops, create_masked_input

def main():
    seed_everything(42)
    os.makedirs("outputs", exist_ok=True)
    os.makedirs("outputs/rfd3", exist_ok=True)
    os.makedirs("outputs/rf3", exist_ok=True)
    
    cif_file = PDBxFile.read("../Data/8FNS.cif")
    atom_array = get_cif_structure(cif_file, model=1)
    
    loops = get_ef_hand_loops(atom_array, "ND")
    print(f"Found {len(loops)} metal binding sites.")
    
    records = []
    
    # Preload engines to avoid reloading weights 4 times
    print("\nLoading RFd3...")
    rfd3_engine = RFD3InferenceEngine(
        diffusion_batch_size=4,
        specification={'length': 105, 'extra': {}} # 8FNS monomer is ~105 residues
    )
    
    print("\nLoading LigandMPNN...")
    mpnn_engine = MPNNInferenceEngine(
        model_type="ligand_mpnn",
        is_legacy_weights=True,
        out_directory=None,
        write_structures=False,
        write_fasta=False,
    )
    
    print("\nLoading RF3...")
    rf3_engine = RF3InferenceEngine(ckpt_path='rf3', verbose=False)
    
    for loop_info in loops:
        print(f"\n--- Processing Metal {loop_info['metal_idx']} at {loop_info['metal_chain_id']}:{loop_info['metal_res_id']} ---")
        
        # 1. Prepare Input
        context_array, start_res, end_res = create_masked_input(atom_array, loop_info)
        print(f"Redesigning residues {start_res} to {end_res}. Context sizes: {len(context_array)} atoms.")
        
        # Determine fixed atoms string for RFd3
        # Ligand is ND
        # We want to redesign the gap.
        # So we fix everything in the context array (which has the metal and the rest of the protein)
        # RFd3 input conditioning (in code) is tricky, we might need to rely on the JSON config
        # OR we can just pass the AtomArray and hope the default `run` method handles it.
        # Actually, let's look at how RFd3 expects it in the Python API. 
        # For simplicity, we'll write a JSON config for each loop and execute via bash OR pass kwargs.
        
        # Define RFd3 Job Spec Dictionary. By passing this dict, `RFD3InferenceEngine` will 
        # ingest it the same way as `inputs=file.json` 
        # We need to redesign `start_res` to `end_res` taking care of keeping the ligand 'ND' fixed. 
        # The easiest approach is to dump the `context_array` as a temporary `.cif` file
        # and pass that as input to RFd3, and just tell it to do 105 length. 
        # If we dump the context_array containing the metal but missing the target loop, 
        # we can tell RFd3 to fix everything in the input, and generate the missing part.
        
        loop_cif_path = f"outputs/rfd3/context_{loop_info['metal_idx']}.cif"
        to_cif_file(context_array, loop_cif_path)
        
        # In RFd3, we want to design around the `ND` ligand.
        # "select_fixed_atoms": {"*": ""} fixes everything BUT the missing residues. Wait, 
        # actually, if we pass an input PDB with missing residues, we can just fix ALL provided residues.
        # Let's specify select_fixed_res or simply fix the whole input.
        # According to standard RFd3 behavior: `length` defines output length.
        # `model_type` defaults to ligand. We want to bury the metal.
        spec_dict = {
            f"loop_{loop_info['metal_idx']}": {
                "input": loop_cif_path,
                "length": f"{105}-{105}",
                "ligand": "ND",
                # "select_fixed_atoms": {"ND": ""},  
                "select_buried": {"ND": "*"},
                "extra": {}
            }
        }
        
        print("Running RFd3...")
        rfd3_outputs_dict = rfd3_engine.run(
            inputs=spec_dict,
            n_batches=1, # 1 batch of 4 (diffusion_batch_size)
            out_dir=None
        )
        # Outputs dict has { "loop_0_0": [RFd3Output, RFd3Output, ...] }
        # the key will be `loop_0_0` based on our prefix and batch logic.
        rfd3_key = f"loop_{loop_info['metal_idx']}_0"
        
        if rfd3_outputs_dict is None or rfd3_key not in rfd3_outputs_dict:
            print("Failed RFd3 generation.")
            continue
            
        rfd3_outputs = rfd3_outputs_dict[rfd3_key]
        
        for design_idx, rfd3_out in enumerate(rfd3_outputs):
            design_id = f"loop_{loop_info['metal_idx']}_design_{design_idx}"
            rfd3_atom_array = rfd3_out.atom_array
            
            # Save RFd3 backbone
            to_cif_file(rfd3_atom_array, f"outputs/rfd3/{design_id}.cif", file_type="cif")
            
            # 2. Sequence Design
            print(f"Running MPNN for {design_id}...")
            mpnn_input_configs = [{
                "batch_size": 1, # Generate 1 sequence for this specific backbone
                "remove_waters": True,
            }]
            
            mpnn_outputs = mpnn_engine.run(input_dicts=mpnn_input_configs, atom_arrays=[rfd3_atom_array])
            if not mpnn_outputs:
                print(f"Failed MPNN for {design_id}")
                continue
                
            mpnn_out = mpnn_outputs[0]
            # Convert sequences to 1-letter
            res_starts = get_residues(mpnn_out.atom_array)[0] # Biotite >= 0.40 uses differently, handle carefully
            # Let's just use string parsing if needed:
            res_names = []
            for rs in res_starts:
                res_names.append(mpnn_out.atom_array.res_name[rs])
            seq_letters = []
            for rn in res_names:
                try:
                    seq_letters.append(ProteinSequence.convert_letter_3to1(rn))
                except:
                    seq_letters.append("X") # For ligand or unnatural
            seq_1letter = "".join(seq_letters)
            
            # 3. Structure Verification
            print(f"Running RF3 validation for {design_id}...")
            # Re-fold the sequence using RF3
            input_structure = InferenceInput.from_atom_array(mpnn_out.atom_array, example_id=design_id)
            rf3_outputs_dict = rf3_engine.run(inputs=input_structure)
            
            if design_id not in rf3_outputs_dict:
                print(f"Failed RF3 validation for {design_id}")
                continue
                
            rf3_output = rf3_outputs_dict[design_id][0]
            rf3_atom_array = rf3_output.atom_array
            
            # Save RF3 validation
            to_cif_file(rf3_atom_array, f"outputs/rf3/{design_id}_refolded.cif", file_type="cif")
            
            # Calculate RMSD
            bb_generated = rfd3_atom_array[np.isin(rfd3_atom_array.atom_name, PROTEIN_BACKBONE_ATOM_NAMES)]
            bb_refolded = rf3_atom_array[np.isin(rf3_atom_array.atom_name, PROTEIN_BACKBONE_ATOM_NAMES)]
            
            # Must be same length for superimpose. They should be, but let's be safe.
            min_len = min(len(bb_generated), len(bb_refolded))
            bb_generated = bb_generated[:min_len]
            bb_refolded = bb_refolded[:min_len]
            
            bb_refolded_fitted, _ = superimpose(bb_generated, bb_refolded)
            rmsd_value = rmsd(bb_generated, bb_refolded_fitted)
            
            summary = rf3_output.summary_confidences
            plddt = summary['overall_plddt']
            
            print(f">>> RMSD: {rmsd_value:.2f} | pLDDT: {plddt:.3f}")
            
            # 4. Record Results
            records.append({
                "domain_idx": loop_info['metal_idx'],
                "design_idx": design_idx,
                "sequence": seq_1letter,
                "rmsd_vs_rfd3": rmsd_value,
                "plddt": plddt,
                "ptm": summary.get('ptm', 0)
            })
            
    # Save Catalog
    df = pd.DataFrame(records)
    print("\nFinal Results:")
    print(df)
    df.to_csv("outputs/sequence_catalog.csv", index=False)


if __name__ == "__main__":
    main()
