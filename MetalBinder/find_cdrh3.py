
from Bio.PDB import PDBParser
import sys

parser = PDBParser(QUIET=True)
structure = parser.get_structure("temp", "../Local/Templates/human_VH3_IgG.pdb")

for model in structure:
    for chain in model:
        if chain.id == 'H':
            print(f"Chain {chain.id}:")
            residues = list(chain.get_residues())
            for res in residues:
                res_id = res.id[1]
                res_name = res.resname
                if 80 <= res_id <= 120:
                    print(f"  {res_id}: {res_name}")
            break
