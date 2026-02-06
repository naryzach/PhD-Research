import argparse
import os
import sys
import glob
import subprocess
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
import numpy as np

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

def get_contig_and_prepare_pdb(motif_path, ion, temp_pdb_path):
    """
    Parse motif PDB, standardize Metal to Chain Z, and generate RFdiffusion contig.
    Writes a temporary standardized PDB.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("motif", motif_path)
    
    residues = []
    metal_atom = None
    
    # 1. Identify Metal and Residues
    for model in structure:
        for chain in model:
            for residue in chain:
                # Check for Metal
                # Heuristic: Check element or resname against Ion
                # Or check if HETATM
                
                is_metal = False
                for atom in residue:
                    if atom.element.upper() == ion or residue.resname.upper() == ion:
                        # Found metal
                        metal_atom = atom
                        is_metal = True
                        break
                
                if is_metal:
                    continue # specific handling later
                
                if is_aa(residue, standard=True):
                    # Store (Chain ID, ResidueID object)
                    residues.append((chain.id, residue.id[1]))

    if not residues:
        return None
        
    if not metal_atom:
         print(f"  Warning: No metal '{ion}' found in {os.path.basename(motif_path)}")
         # Continue without metal? Or Fail?
         # If no metal, we can't design a binder FOR the metal properly.
         return None

    # 2. Write Standardized PDB (Protein Chains preserved, Metal -> Chain Z Res 1)
    io = PDBIO()
    
    # We need to construct a new structure or modify existing?
    # Modifying existing in place is risky for BioPython objects if complex.
    # Easiest: Write protein parts, then append metal line manually.
    
    # Write protein parts
    class ProteinSelect(Select):
        def accept_residue(self, residue):
             return is_aa(residue, standard=True)
             
    io.set_structure(structure)
    io.save(temp_pdb_path, ProteinSelect())
    
    # Append Metal as Chain Z, Res 1
    # Note: PDBIO.save() typically writes an END record. We must insert the metal BEFORE END, or remove END.
    
    # 1. Read the saved protein-only PDB
    with open(temp_pdb_path, 'r') as f:
        lines = f.readlines()
        
    # 2. Filter out END/ENDMDL and TER
    # Use strip() to handle trailing spaces or leading spaces if any
    lines = [l for l in lines if not l.strip().startswith("END") and not l.strip().startswith("TER")]
    
    # 3. Create Metal Line

    x, y, z = metal_atom.get_coord()
    
    # Native Metal Handling
    # Use HETATM and correct residue name
    metal_lines = []
    metal_lines.append(f"HETATM{9999:5d}  {ion:<4}{ion:<3} Z   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {ion.rjust(2)}\n")

    
    print(f"    Writing Metal Lines (N, CA, C) for {ion}...")
    
    # 4. Write back
    with open(temp_pdb_path, 'w') as f:
        f.writelines(lines)
        for line in metal_lines:
            f.write(line)
        f.write("END\n")
    

        
    # 3. Build Contig
    # Group residues
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
        
    # Construct RFstring
    # Strategy: [50-50/Seg1/10-20/Seg2/50-50 Z1-1]
    contig_parts = ["5-50"] 
    for seg in segments:
        contig_parts.append(seg)
        contig_parts.append("1-20") 
    contig_parts[-1] = "5-50" # End with flexible
    
    contig_str = "/".join(contig_parts)
    # Add Metal (Chain Z, Res 1)
    # CRITICAL: Add '/0' chain break so metal is not peptide-bonded to protein
    contig = f"['{contig_str}/0 Z1-1']"
    
    return contig

def run_scaffolding_for_ion(ion, catalog_dir, out_dir, args):
    ion_cat_dir = os.path.join(catalog_dir, ion)
    ion_out_dir = os.path.join(out_dir, ion)
    
    if not os.path.exists(ion_cat_dir):
        print(f"No catalog found for {ion} in {ion_cat_dir}")
        return
    
    motifs = glob.glob(os.path.join(ion_cat_dir, "*.pdb"))
    print(f"Found {len(motifs)} motifs for {ion}.")
    
    project_root = os.path.dirname(os.path.abspath(__file__))
    
    for motif_path in motifs:
        motif_name = os.path.basename(motif_path).replace(".pdb", "")
        print(f"  Scaffolding {motif_name}...")
        
        scaffold_out_dir = os.path.join(ion_out_dir, motif_name)
        rf_out = os.path.join(scaffold_out_dir, "rfdiffusion")
        
        if os.path.exists(rf_out) and not args.overwrite:
            # Check if non-empty
            if glob.glob(os.path.join(rf_out, "*.pdb")):
                print(f"    Skipping {motif_name} (already exists)")
                continue
            
        os.makedirs(rf_out, exist_ok=True)
        
        # Prepare Input PDB with explicit metal chain & generate contig
        temp_input = os.path.join(rf_out, "input_with_metal.pdb")
        contig = get_contig_and_prepare_pdb(motif_path, ion, temp_input)
        
        if not contig:
            print(f"    Skipping {motif_name} (failed to generate contig/metal)")
            continue
            
        # Run RFdiffusion
        # Use main script (now patched with MKL fix)
        rf_script = os.path.join(args.rf_path, "scripts/run_inference.py")
        
        cmd = [
            sys.executable, rf_script,
            f"inference.input_pdb={temp_input}",
            f"inference.output_prefix={os.path.join(rf_out, 'scaffold')}",
            f"inference.num_designs={args.num_designs}",
            f"contigmap.contigs={contig}"
        ]
        
        # Inject RFDIFFUSION_ION env var for generalized metal support
        env = os.environ.copy()
        env["RFDIFFUSION_ION"] = ion # e.g. "ZN", "CU"

        if args.dry_run:
            print(f"    Dry Run RFdiffusion (Ion {ion}): {' '.join(cmd)}")
        else:
            try:
                subprocess.run(cmd, check=True, cwd=project_root, env=env)
                
                # Check output?
                # RFdiffusion usually preserves fixed input (including metal) in output.
                # But just in case, verify?
                # Usually fine.
                
                # Cleanup temp file?
                # os.remove(temp_input) 
                # Keep it for debug.
                
            except Exception as e:
                print(f"    RFdiffusion failed: {e}")

def main():
    parser = argparse.ArgumentParser(description="De Novo Scaffolding Workflow")
    parser.add_argument("--catalog_dir", default="../Local/Metal_Catalog", help="Input Catalog Directory")
    parser.add_argument("--out_dir", default="../Local/Metal_Predictions/DeNovo", help="Output Directory")
    parser.add_argument("--rf_path", default="../Tools/RFdiffusion", help="RFdiffusion Path")
    parser.add_argument("--metals", default="ZN,CU,NI,CO,MN,FE,MG,CA", help="Metals to process")
    parser.add_argument("--num_designs", type=int, default=2)
    parser.add_argument("--dry_run", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    
    args = parser.parse_args()
    
    metals = [m.strip().upper() for m in args.metals.split(',')]
    
    print(f"Starting Scaffolding Workflow for: {metals}")
    
    for ion in metals:
        run_scaffolding_for_ion(ion, args.catalog_dir, args.out_dir, args)

if __name__ == "__main__":
    main()
