import os
import sys
import numpy as np
import subprocess
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

class ChainSelect(Select):
    """
    Select specific chain and residues.
    """
    def __init__(self, chain_id, residues=None):
        self.chain_id = chain_id
        self.residues = residues

    def accept_chain(self, chain):
        return chain.id == self.chain_id
    
    def accept_residue(self, residue):
        if self.residues:
             return residue.id[1] in self.residues
        return True

class CleanSelect(Select):
    """
    Selects Protein and specific Metal atoms.
    """
    def __init__(self, ion_name):
        self.ion_name = ion_name.upper()
        
    def accept_residue(self, residue):
        # Keep Amino Acids
        if is_aa(residue, standard=True):
            return 1
        # Keep Target Metal (by resname or element)
        if residue.resname.strip().upper() == self.ion_name:
            return 1
        # Check atoms inside if needed
        for atom in residue:
            if atom.element.upper() == self.ion_name:
                return 1
        return 0

def find_model_weights(rf_path):
    """
    Finds the latest model weights file in the RFdiffusion2 directory.
    """
    weights_dir = os.path.join(os.path.abspath(rf_path), "rf_diffusion", "model_weights")
    if os.path.exists(weights_dir):
        pt_files = [f for f in os.listdir(weights_dir) if f.endswith(".pt") and f.startswith("RFD")]
        if pt_files:
            pt_files.sort()
            return os.path.join(weights_dir, pt_files[-1])
    return None

def get_contig_str(residues):
    """
    Generates an RFdiffusion2 compatible contig string from a list of residues (chain, res_id).
    Format: '10-60,ChainStart-End,1-20,ChainStart-End,10-60'
    """
    if not residues:
        return ""
    
    residues.sort(key=lambda x: (x[0], x[1]))
    
    segments = []
    if residues:
        current_chain = residues[0][0]
        start_res = residues[0][1]
        prev_res = residues[0][1]
        
        for i in range(1, len(residues)):
            chain, res_id = residues[i]
            if chain == current_chain and res_id == prev_res + 1:
                prev_res = res_id
            else:
                segments.append(f"{current_chain}{start_res}-{prev_res}")
                current_chain = chain
                start_res = res_id
                prev_res = res_id
        segments.append(f"{current_chain}{start_res}-{prev_res}")
    
    # Construct RFdiffusion2 contig string
    contig_parts = ["10-60"] 
    for i, seg in enumerate(segments):
        contig_parts.append(seg)
        if i < len(segments) - 1:
            contig_parts.append("1-20") 
    contig_parts.append("10-60")
    
    return f"['{','.join(contig_parts)}']"

def add_ori_token(pdb_path, ion_names):
    """
    Calculates the center of mass of the specified ions and appends an ORI HETATM token to the PDB file.
    ion_names: list of strings (e.g. ['ZN', 'CU'])
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("clean", pdb_path)
    except Exception as e:
        print(f"Error parsing for ORI addition: {e}")
        return False
        
    target_atoms = []
    ion_set = set([i.upper() for i in ion_names])
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname.strip().upper() in ion_set:
                    for atom in residue:
                        target_atoms.append(atom)
                else:
                    for atom in residue:
                        if atom.element.upper() in ion_set:
                            target_atoms.append(atom)
                            
    if target_atoms:
        coords = np.array([a.get_coord() for a in target_atoms])
        center = np.mean(coords, axis=0)
        
        with open(pdb_path, 'a') as f:
            x, y, z = center
            ori_line = f"HETATM 9999  CA  ORI Z 999    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n"
            f.write(ori_line)
            return True
    return False
