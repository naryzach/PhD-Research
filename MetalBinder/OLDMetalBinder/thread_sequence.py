import argparse
import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from openmm import app
import os

def thread_sequence(input_pdb, output_pdb, target_sequence_str):
    """
    Threads a sequence onto a PDB structure using PDBFixer.
    Replaces residues and rebuilds sidechains.
    """
    print(f"Reading {input_pdb}...")
    fixer = PDBFixer(filename=input_pdb)
    
    # target_sequence_str is likely "/"-separated chains e.g. "AAAA/GGGG"
    # or just "AAAA" if single chain.
    # PDBFixer topology: list of chains.
    
    chains = list(fixer.topology.chains())
    target_seqs = target_sequence_str.split("/")
    
    if len(target_seqs) != len(chains):
        print(f"Warning: Chain mismatch! PDB has {len(chains)} chains, Target has {len(target_seqs)} parts.")
        # Attempt to match? Or just apply to first N chains?
        # Assuming Rdfiffusion output order matches MPNN output order.
    
    # Apply mutations
    # PDBFixer.applyMutations(mutations, chain_id)
    # mutations: list of Strings "ALA-123-CYS" (Original-ResID-New)
    # Actually PDBFixer doesn't have a direct "mutate sequence" high level API easily.
    # But OpenMM Modeller does allow mutating if we delete and add?
    # PDBFixer has 'applyMutations'.
    
    # Better approach with PDBFixer:
    # working with the topology directly.
    # Actually, simpler approach:
    # 1. Identify residues in topology.
    # 2. Construct valid mutation strings for applyMutations.
    
    # Mapping amino acid codes
    one_to_three = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
        'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
        'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
        'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR',
        'X': 'ALA' # Default to ALA for unknown
    }
    
    for i, chain in enumerate(chains):
        if i >= len(target_seqs): break
        
        tgt_seq = target_seqs[i]
        
        # Iterating residues in this chain
        residues = list(chain.residues())
        
        # If lengths differ (e.g. ligand chain vs protein seq), careful.
        # Ligand chains from MPNN might be labeled "X" or short.
        # If target seq is short, maybe skip?
        
        if len(tgt_seq) != len(residues):
             print(f"Warning: Chain {chain.id} length mismatch ({len(residues)} res vs {len(tgt_seq)} seq). Skipping mutations for this chain.")
             continue
             
        mutations = []
        for j, res in enumerate(residues):
            target_aa = tgt_seq[j]
            if target_aa not in one_to_three:
                target_aa = 'X' # Unknown
            
            target_resname = one_to_three[target_aa]
            
            # Check if mutation needed (or always mutate to ensure sidechain rebuild?)
            if res.name != target_resname:
                # Format: "ALA-123-CYS"
                # ResID is res.id (which is a string like "123")
                mutations.append(f"{res.name}-{res.id}-{target_resname}")
        
        if mutations:
            print(f"Applying {len(mutations)} mutations to chain {chain.id}...")
            fixer.applyMutations(mutations, chain.id)

    # Now rebuild missing atoms (sidechains)
    print("Finding missing atoms...")
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    
    print("Adding missing atoms (sidechains)...")
    fixer.addMissingAtoms()
    
    print("Adding missing hydrogens...")
    fixer.addMissingHydrogens(7.0)
    
    print(f"Writing to {output_pdb}...")
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sequence", required=True, help="Sequence string (e.g. AAAA/GGGG)")
    args = parser.parse_args()
    
    try:
        thread_sequence(args.input, args.output, args.sequence)
    except Exception as e:
        print(f"Error threading: {e}")
        sys.exit(1)
