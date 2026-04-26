import sys
import numpy as np
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner

def get_sequence_and_residues(structure, chain_id='A'):
    residues = []
    seq = ""
    try:
        chain = structure[0][chain_id]
    except KeyError:
        return "", []
        
    for res in chain:
        if is_aa(res):
            try:
                # Use 1-letter code
                aa = seq1(res.get_resname())
                seq += aa
                residues.append(res)
            except KeyError:
                continue
    return seq, residues

def compare_pdbs_aligned(pdb1, pdb2, chain_id='A'):
    p = PDBParser(QUIET=True)
    s1 = p.get_structure("s1", pdb1)
    s2 = p.get_structure("s2", pdb2)
    
    seq1, res1 = get_sequence_and_residues(s1, chain_id)
    seq2, res2 = get_sequence_and_residues(s2, chain_id)
    
    if not seq1 or not seq2:
        print("Error reading sequences.")
        return

    # Align sequences to find correspondence
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    alignments = aligner.align(seq1, seq2)
    alignment = alignments[0] # Best alignment
    
    print(f"Sequence Alignment Score: {alignment.score}")
    # print(alignment)
    
    # Map residues based on alignment
    # alignment.indices gives indices of matches in seq1 and seq2
    # But filtering for 'identity' is better manually.
    
    atoms1 = []
    atoms2 = []
    
    # Iterate through the alignment to find matched columns
    # We can use alignment.aligned which returns list of (start, end) tuples for matches
    # but that might include mismatches if penalties allow? No, checking alignment directly is safer.
    
    # Simple approach: walk through the aligned strings
    
    # Wait, simple approach:
    pattern1, pattern2 = alignment.__str__().split('\n')[0], alignment.__str__().split('\n')[2]
    # No, alignment object printing is complex.
    
    # Use aligned indices
    # shape of alignment.coordinates is (N+1, 2) path.
    # We want matched pairs.
    
    idx1 = 0
    idx2 = 0
    
    # Get the aligned sequences (with gaps)
    # The new Bio.Align 1.8+ format:
    # We can iterate over the alignment coordinates?
    
    # Let's simple check:
    # Just map identical residues that are aligned.
    
    aligned_pairs = []
    
    # Extract aligned indices
    # alignment.indices are indices of the second sequence aligned to the first.
    # alignment.aligned is easier.
    
    match_count = 0
    
    # Use alignment.aligned to find matching segments
    # aligned is a tuple (ranges1, ranges2)
    # ranges1 is list of (start, end) in seq1
    # ranges2 is list of (start, end) in seq2
    
    ranges1, ranges2 = alignment.aligned
    
    match_count = 0
    
    for (start1, end1), (start2, end2) in zip(ranges1, ranges2):
        # The lengths must be equal for pairwise alignment segments
        length = end1 - start1
        assert length == end2 - start2
        
        for k in range(length):
            idx1 = start1 + k
            idx2 = start2 + k
            
            # Check for identity
            if seq1[idx1] == seq2[idx2]:
                 # Add CA atoms
                 # Handle missing atoms or index out of range (though seq comes from structure)
                 if idx1 < len(res1) and idx2 < len(res2):
                    r1 = res1[idx1]
                    r2 = res2[idx2]
                    if 'CA' in r1 and 'CA' in r2:
                        atoms1.append(r1['CA'])
                        atoms2.append(r2['CA'])
                        match_count += 1
                
    print(f"Aligned {match_count} identical residues (Scaffold context).")
    
    if len(atoms1) < 10:
        print("Not enough aligned residues.")
        return

    # Calculate Superimposed RMSD
    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    print(f"Whole Structure Superimposed RMSD: {sup.rms:.4f} A")
    
    if sup.rms < 1.0:
        print("PASS: Structure is preserved (RMSD < 1.0 A)")
    else:
        print("FAIL: Structure is distorted!")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python compare_pdbs.py <pdb1> <pdb2>")
        sys.exit(1)
        
    compare_pdbs_aligned(sys.argv[1], sys.argv[2])
