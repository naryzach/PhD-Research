import argparse
import os
import sys
import glob
import subprocess
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
import argparse
import os
import sys
import glob
import subprocess
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Polypeptide import is_aa

# Import shared utilities
try:
    from utils_rfd2 import CleanSelect, find_model_weights, get_contig_str, add_ori_token
except ImportError:
    # Handle running from same directory
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from utils_rfd2 import CleanSelect, find_model_weights, get_contig_str, add_ori_token


def get_contig_and_prepare_pdb(motif_path, ion, temp_pdb_path):
    """
    Parse motif PDB, standardize input, and generate RFdiffusion2 contig.
    Writes a temporary standardized PDB with ORI token.
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("motif", motif_path)
    except Exception as e:
        print(f"Error parsing {motif_path}: {e}")
        return None
    
    residues = []
    metal_found = False
    
    # 1. Identify Metal and Residues
    detected_ligands = set()
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if is_aa(residue, standard=True):
                    # Store (Chain ID, ResidueID object)
                    residues.append((chain.id, residue.id[1]))
                
                # Check for Metal
                is_metal = False
                if residue.resname.strip().upper() == ion.upper():
                    is_metal = True
                else:
                    for atom in residue:
                        if atom.element.upper() == ion.upper():
                            is_metal = True
                            break
                if is_metal:
                    metal_found = True
                    detected_ligands.add(residue.resname.strip().upper())

    if not residues:
        return None, None
        
    if not metal_found:
         print(f"  Warning: No metal '{ion}' found in {os.path.basename(motif_path)}")

    # 2. Write Standardized PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(temp_pdb_path, CleanSelect(ion))
    
    # 3. Add ORI Token
    if not add_ori_token(temp_pdb_path, [ion]):
        print("    Warning: Could not calculate metal center for ORI token.")

    # 4. Build Contig String
    contig = get_contig_str(residues)
    
    # Return contig and comma-separated ligand names
    ligand_str = ",".join(list(detected_ligands)) if detected_ligands else ion
    return contig, ligand_str


def run_scaffolding_for_ion(ion, catalog_dir, out_dir, args):
    ion_cat_dir = os.path.join(catalog_dir, ion)
    ion_out_dir = os.path.join(out_dir, ion)
    
    if not os.path.exists(ion_cat_dir):
        print(f"No catalog found for {ion} in {ion_cat_dir}")
        return
    
    motifs = glob.glob(os.path.join(ion_cat_dir, "*.pdb"))
    print(f"Found {len(motifs)} motifs for {ion}.")
    
    # Find checkpoint file
    ckpt_path = find_model_weights(args.rf_path)
    if ckpt_path:
        print(f"  Using model weights: {os.path.basename(ckpt_path)}")
    else:
        print("  Warning: No model weights found in RFdiffusion2 directory.")
    
    project_root = os.path.dirname(os.path.abspath(__file__)) 
    
    # RFdiffusion2 script
    rf_script = os.path.join(os.path.abspath(args.rf_path), "rf_diffusion", "run_inference.py")

    for motif_path in motifs:
        motif_name = os.path.basename(motif_path).replace(".pdb", "")
        print(f"  Scaffolding {motif_name}...")
        
        scaffold_out_dir = os.path.join(ion_out_dir, motif_name)
        rf_out = os.path.join(scaffold_out_dir, "rfdiffusion")
        
        if os.path.exists(rf_out) and not args.overwrite:
            if glob.glob(os.path.join(rf_out, "scaffold*.pdb")):
                print(f"    Skipping {motif_name} (already exists)")
                continue
            
        os.makedirs(rf_out, exist_ok=True)
        
        # Prepare Input PDB
        temp_input = os.path.join(rf_out, "input_with_metal.pdb")
        contig, ligand_name = get_contig_and_prepare_pdb(motif_path, ion, temp_input)
        
        if not contig:
            print(f"    Skipping {motif_name} (failed to generate contig)")
            continue
            
        # Run RFdiffusion2
        output_prefix = os.path.join(rf_out, "scaffold")
        
        cmd = [
            sys.executable,
            rf_script,
            f"inference.input_pdb={os.path.abspath(temp_input)}",
            f"inference.output_prefix={os.path.abspath(output_prefix)}",
            f"inference.num_designs={args.num_designs}",
            f"contigmap.contigs={contig}",
            f"inference.ligand='{ligand_name}'"
        ]

        
        if ckpt_path:
            cmd.append(f"inference.ckpt_path={ckpt_path}")

        if args.dry_run:
            print(f"    Dry Run (Ion {ion}): {' '.join(cmd)}")
        else:
            try:
                # Ensure PYTHONPATH includes RFdiffusion2
                env = os.environ.copy()
                rf_abs = os.path.abspath(args.rf_path)
                if "PYTHONPATH" in env:
                    env["PYTHONPATH"] = f"{rf_abs}:{env['PYTHONPATH']}"
                else:
                    env["PYTHONPATH"] = rf_abs
                
                subprocess.run(cmd, check=True, env=env, cwd=rf_abs)
                
            except Exception as e:
                print(f"    RFdiffusion failed: {e}")

def main():
    parser = argparse.ArgumentParser(description="De Novo Scaffolding Workflow (RFdiffusion2)")
    parser.add_argument("--catalog_dir", default="../Local/Metal_Catalog", help="Input Catalog Directory")
    parser.add_argument("--out_dir", default="../Local/Metal_Predictions/DeNovo", help="Output Directory")
    parser.add_argument("--rf_path", default="../Tools/RFdiffusion2", help="RFdiffusion2 Path")
    parser.add_argument("--metals", default="ZN,CU,NI,CO,MN,FE,MG,CA", help="Metals to process")
    parser.add_argument("--num_designs", type=int, default=2)
    parser.add_argument("--dry_run", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    
    args = parser.parse_args()
    
    # Metal Selection Logic
    if args.metals.upper() == "ALL":
        # Scan catalog for all available ions
        if os.path.exists(args.catalog_dir):
            metals = [d for d in os.listdir(args.catalog_dir) 
                      if os.path.isdir(os.path.join(args.catalog_dir, d))]
            metals = [m.upper() for m in metals]
        else:
            print(f"Catalog directory not found at {args.catalog_dir}")
            return
    else:
        metals = [m.strip().upper() for m in args.metals.split(',')]
    
    print(f"Starting Scaffolding Workflow for: {metals}")
    
    for ion in metals:
        run_scaffolding_for_ion(ion, args.catalog_dir, args.out_dir, args)

if __name__ == "__main__":
    main()

