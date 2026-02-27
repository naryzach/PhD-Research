import os
import sys
import numpy as np
import argparse
import torch
from biotite.structure.io.pdbx import CIFFile, get_structure as get_cif_structure
from biotite.structure.io.pdb import PDBFile
from biotite.sequence import ProteinSequence

# Foundry Inference imports
from mpnn.inference_engines.mpnn import MPNNInferenceEngine

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
    parser = argparse.ArgumentParser(description="Test LigandMPNN fixed region design")
    parser.add_argument("--input", type=str, default="Data/8FNS.cif", help="Path to input CIF/PDB structure")
    parser.add_argument("--chain", type=str, default="A", help="Chain to process")
    parser.add_argument("--design-regions", type=str, nargs="+", default=["51-60"], help="List of regions to design, e.g. 51-60 80-90. All other residues will be fixed.")
    args = parser.parse_args()

    input_file = args.input
    output_dir = "outputs/test_lmpnn"
    os.makedirs(output_dir, exist_ok=True)
    chain_id = args.chain
    
    print(f"--- Loading input structure: {input_file} ---")
    if not os.path.exists(input_file):
        # Fallback to local test if run from wrong directory
        if os.path.exists(os.path.join("..", input_file)):
            input_file = os.path.join("..", input_file)
        else:
            print(f"Error: input file '{input_file}' not found.")
            return

    if input_file.endswith(".cif"):
        cif_file = CIFFile.read(input_file)
        try:
            atom_array = get_cif_structure(cif_file, model=1)
        except:
            atom_array = get_cif_structure(cif_file)
            if isinstance(atom_array, list): atom_array = atom_array[0]
    elif input_file.endswith(".pdb"):
        pdb_file = PDBFile.read(input_file)
        atom_array = pdb_file.get_structure(model=1)
    else:
        print("Unsupported file format. Please use .cif or .pdb")
        return

    # Extract original full sequence
    original_full_seq = get_sequence_from_array(atom_array, chain_id=chain_id)
    print(f"\nOriginal sequence (Chain {chain_id}):\n{original_full_seq}")
    
    # Identify all residue IDs in the chain
    mask = (atom_array.chain_id == chain_id) & (atom_array.atom_name == "CA")
    ca_atoms = atom_array[mask]
    all_res_ids = np.unique(ca_atoms.res_id)
    
    designed_res_ids = []
    for region in args.design_regions:
        try:
            start, end = map(int, region.split('-'))
            designed_res_ids.extend(list(range(start, end + 1)))
        except ValueError:
            print(f"Error: Invalid region format '{region}'. Use start-end (e.g., 51-60).")
            return
            
    designed_res_ids = sorted(list(set(designed_res_ids))) # Ensure unique and sorted
    
    # Validation
    valid_designed = [r for r in designed_res_ids if r in all_res_ids]
    if len(valid_designed) == 0:
        print(f"Error: None of the requested designed residues are in the structure.")
        print(f"Available residues in chain {chain_id}: {all_res_ids.min()}-{all_res_ids.max()}")
        return
        
    # The rest are fixed
    fixed_res_ids = [r for r in all_res_ids if r not in designed_res_ids]
    
    fixed_residues_str = [f"{chain_id}{r}" for r in fixed_res_ids]

    print(f"\nDesigned residues: {', '.join([f'{chain_id}{r}' for r in args.design_regions])}")
    print(f"Total designed residues: {len(designed_res_ids)}")
    print(f"Fixed residues: {len(fixed_res_ids)} residues will be fixed.")
    
    # Set Foundry checkpoint dir from project root
    checkpoint_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Tools", "foundry_checkpoints"))
    os.environ["FOUNDRY_CHECKPOINT_DIRS"] = checkpoint_dir

    # Initialize engine
    print("\n--- Initializing LigandMPNN Engine ---")
    engine = MPNNInferenceEngine(
        model_type="ligand_mpnn", 
        is_legacy_weights=True, 
        out_directory=output_dir,
        write_structures=True,
        write_fasta=True
    )
    
    # Setup input dictionary
    input_dict = {
        "name": "test_fixed_design",
        "batch_size": 1,
        "fixed_residues": fixed_residues_str,
        "seed": 42
    }
    
    print("\n--- Running LigandMPNN Inference ---")
    outputs = engine.run(input_dicts=[input_dict], atom_arrays=[atom_array])
    
    if not outputs:
        print("Error: No outputs from LigandMPNN.")
        return
        
    result = outputs[0]
    result_array = result.atom_array
    
    # Extract designed sequence parts
    designed_full_seq = get_sequence_from_array(result_array, chain_id=chain_id)
    
    print(f"\nGenerated sequence (Chain {chain_id}):\n{designed_full_seq}")
    
    # Verify using Biotite method
    original_fixed_seq = get_sequence_from_array(atom_array, chain_id=chain_id, res_ids=fixed_res_ids)
    designed_fixed_seq = get_sequence_from_array(result_array, chain_id=chain_id, res_ids=fixed_res_ids)
    
    print("\n--- Verification ---")
    if original_fixed_seq == designed_fixed_seq:
        print("✅ FIXED region remains unchanged.")
    else:
        print("❌ FIXED region has CHANGED!")
        
    print("\n✅ DESIGNED regions:")
    for region in args.design_regions:
        start, end = map(int, region.split('-'))
        region_res_ids = list(range(start, end + 1))
        
        orig_seq = get_sequence_from_array(atom_array, chain_id=chain_id, res_ids=region_res_ids)
        gen_seq  = get_sequence_from_array(result_array, chain_id=chain_id, res_ids=region_res_ids)
        
        if not orig_seq:
             print(f"   Region {region} (Missing from structure)")
             continue

        if orig_seq != gen_seq:
            print(f"   Region {region} (Modified)")
            print(f"     Original : {orig_seq}")
            print(f"     Generated: {gen_seq}")
        else:
            print(f"   Region {region} (Unchanged)")
            print(f"     Original : {orig_seq}")
            print(f"     Generated: {gen_seq}")

    print("\n--- Full Sequence Alignment ---")
    print(f"Original : {original_full_seq}")
    
    match_str = ""
    for o, g in zip(original_full_seq, designed_full_seq):
        if o == g:
            match_str += "|"
        else:
            match_str += " "
            
    print(f"Alignment: {match_str}")
    print(f"Generated: {designed_full_seq}")

if __name__ == "__main__":
    main()
