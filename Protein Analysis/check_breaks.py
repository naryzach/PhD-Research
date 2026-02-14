import sys
import numpy as np
from Bio.PDB import PDBParser

def check_breaks(pdb_file, threshold=2.0):
    p = PDBParser(QUIET=True)
    s = p.get_structure("check", pdb_file)
    chain = s[0]['A']
    
    residues = list(chain.get_residues())
    # Sort by residue ID just in case
    residues.sort(key=lambda r: r.id[1])
    
    breaks = []
    
    for i in range(len(residues) - 1):
        r1 = residues[i]
        r2 = residues[i+1]
        
        # Check residue numbering continuity
        if r2.id[1] != r1.id[1] + 1:
            print(f"Numbering Gap: {r1.id[1]} -> {r2.id[1]}")
            
        # Check C-N distance
        if 'C' in r1 and 'N' in r2:
            c_atom = r1['C']
            n_atom = r2['N']
            dist = c_atom - n_atom
            if dist > threshold:
                print(f"Break detected: Res {r1.id[1]}-{r2.id[1]} Dist: {dist:.2f} A")
                breaks.append((r1.id[1], r2.id[1], dist))
        else:
            print(f"Missing backbone atoms at {r1.id[1]}-{r2.id[1]}")
            
    if not breaks:
        print("No backbone breaks found.")
    else:
        print(f"Found {len(breaks)} breaks > {threshold} A.")
        
if __name__ == "__main__":
    check_breaks(sys.argv[1])
