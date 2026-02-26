import os
import sys
import numpy as np
import torch
from biotite.structure.io.pdbx import CIFFile, get_structure as get_cif_structure
from biotite.structure import get_residues
from biotite.sequence import ProteinSequence
from Bio import PDB
from Bio.SeqUtils import seq1

# Foundry Inference imports
from mpnn.inference_engines.mpnn import MPNNInferenceEngine

def get_sequence_biotite(atom_array, chain_id="A", res_ids=None):
    """Helper function to extract sequence cleanly from an AtomArray using Biotite."""
    mask = (atom_array.chain_id == chain_id) & (atom_array.atom_name == "CA")
    if res_ids is not None:
        mask = mask & np.isin(atom_array.res_id, res_ids)
        
    ca_atoms = atom_array[mask]
    if len(ca_atoms) == 0: return ""
    
    # Ensure residues are sorted by ID
    sort_idx = np.argsort(ca_atoms.res_id)
    res_names = ca_atoms.res_name[sort_idx]
    
    seq_letters = []
    for rn in res_names:
        try:
            seq_letters.append(ProteinSequence.convert_letter_3to1(rn))
        except Exception:
            seq_letters.append("X")
            
    return "".join(seq_letters)

def get_sequence_biopython(pdb_file, chain_id="A", res_range=None):
    """Helper function to extract sequence using BioPython."""
    parser = PDB.PDBParser(QUIET=True)
    if pdb_file.endswith(".cif"):
        parser = PDB.MMCIFParser(QUIET=True)
    
    structure = parser.get_structure("test", pdb_file)
    model = structure[0]
    chain = model[chain_id]
    
    seq = []
    for residue in chain:
        if PDB.is_aa(residue):
            res_id = residue.get_id()[1]
            if res_range is not None and res_id not in res_range:
                continue
            seq.append(seq1(residue.get_resname()))
            
    return "".join(seq)

def main():
    input_file = "Data/8FNS.cif"
    output_dir = "outputs/test_lmpnn"
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"--- Loading input structure: {input_file} ---")
    cif_file = CIFFile.read(input_file)
    # Get the structure from its first (and likely only) model
    try:
        atom_array = get_cif_structure(cif_file, model=1)
    except:
        # Fallback if model index 1 fails (sometimes model is None or 0 depending on version)
        atom_array = get_cif_structure(cif_file)
        if isinstance(atom_array, list): atom_array = atom_array[0]

    # Define fixed and designed residues
    # Let's fix residues 30-50 and design 51-60 of Chain A
    fixed_res_ids = list(range(30, 51))
    designed_res_ids = list(range(51, 61))
    
    # Convert to strings as expected by the engine (e.g., "A1", "A2", ...)
    fixed_residues_str = [f"A{r}" for r in fixed_res_ids]
    designed_residues_str = [f"A{r}" for r in designed_res_ids]
    
    print(f"Fixed residues: {fixed_residues_str[0]} to {fixed_residues_str[-1]}")
    print(f"Designed residues: {designed_residues_str[0]} to {designed_residues_str[-1]}")
    
    # Extract original sequence for the fixed part
    original_fixed_seq = get_sequence_biotite(atom_array, chain_id="A", res_ids=fixed_res_ids)
    print(f"Original fixed sequence: {original_fixed_seq}")
    
    # Initialize engine
    print("\n--- Initializing LigandMPNN Engine ---")
    # Note: run_pipeline.py uses is_legacy_weights=True
    engine = MPNNInferenceEngine(
        model_type="ligand_mpnn", 
        is_legacy_weights=True, 
        out_directory=output_dir,
        write_structures=True,
        write_fasta=True
    )
    
    # Setup input dictionary
    # Based on MPNN_PER_INPUT_INFERENCE_DEFAULTS in inference.py
    input_dict = {
        "name": "test_fixed_design",
        "batch_size": 1,
        "fixed_residues": fixed_residues_str,
        # If we specify fixed_residues, others are designed by default if not specified otherwise.
        # However, to be explicit we could also set designed_residues but they are mutually exclusive in CLI.
        # In the engine, it seems we can provide either.
        "seed": 42
    }
    
    print("\n--- Running LigandMPNN Inference ---")
    outputs = engine.run(input_dicts=[input_dict], atom_arrays=[atom_array])
    
    if not outputs:
        print("Error: No outputs from LigandMPNN.")
        return
        
    result = outputs[0]
    result_array = result.atom_array
    
    # Extract result sequence for fixed part
    designed_fixed_seq = get_sequence_biotite(result_array, chain_id="A", res_ids=fixed_res_ids)
    print(f"Designed structure fixed sequence: {designed_fixed_seq}")
    
    # Verify using Biotite method
    if original_fixed_seq == designed_fixed_seq:
        print("✅ Sequence verification (Biotite): FIXED region remains unchanged.")
    else:
        print("❌ Sequence verification (Biotite): FIXED region has CHANGED!")
        print(f"Expected: {original_fixed_seq}")
        print(f"Actual:   {designed_fixed_seq}")

    # Verify using Biopython method (reading the output file)
    print("\n--- Verifying using Biopython ---")
    output_cif = os.path.join(output_dir, "test_fixed_design_b0_d0.cif")
    if os.path.exists(output_cif):
        biopython_fixed_seq = get_sequence_biopython(output_cif, chain_id="A", res_range=fixed_res_ids)
        print(f"Biopython fixed sequence from file: {biopython_fixed_seq}")
        if biopython_fixed_seq == original_fixed_seq:
            print("✅ Sequence verification (Biopython): FIXED region verified in output file.")
        else:
            print("❌ Sequence verification (Biopython): FIXED region mismatch in output file.")
    else:
        print(f"Warning: Output file {output_cif} not found for Biopython verification.")

    # Show the designed sequence
    designed_seq_part = get_sequence_biotite(result_array, chain_id="A", res_ids=designed_res_ids)
    print(f"\nDesigned region sequence: {designed_seq_part}")

if __name__ == "__main__":
    main()
