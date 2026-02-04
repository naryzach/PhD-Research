import argparse
import os
import sys
import glob
import subprocess
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

def get_contig_from_motif(motif_path):
    """
    Parse motif PDB and generate RFdiffusion contig string.
    Motif PDB contains only the residues to be fixed + metal.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("motif", motif_path)
    
    residues = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if is_aa(residue, standard=True):
                    # Store (Chain ID, ResidueID)
                    residues.append((chain.id, residue.id[1]))
                    
    if not residues:
        return None

    # Sort residues
    residues.sort(key=lambda x: (x[0], x[1]))
    
    # Group into segments
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
        
    return segments

def run_scaffolding_for_ion(ion, catalog_dir, out_dir, args):
    ion_cat_dir = os.path.join(catalog_dir, ion)
    ion_out_dir = os.path.join(out_dir, ion)
    
    if not os.path.exists(ion_cat_dir):
        print(f"No catalog found for {ion} in {ion_cat_dir}")
        return

    motifs = glob.glob(os.path.join(ion_cat_dir, "*.pdb"))
    print(f"Found {len(motifs)} motifs for {ion}.")
    
    # Calculate project root for correct relative paths if needed
    project_root = os.path.dirname(os.path.abspath(__file__))
    
    for motif_path in motifs:
        motif_name = os.path.basename(motif_path).replace(".pdb", "")
        print(f"  Scaffolding {motif_name}...")
        
        scaffold_out_dir = os.path.join(ion_out_dir, motif_name)
        rf_out = os.path.join(scaffold_out_dir, "rfdiffusion")
        
        if os.path.exists(rf_out) and not args.overwrite:
            print(f"    Skipping {motif_name} (already exists)")
            continue
            
        os.makedirs(rf_out, exist_ok=True)
        
        # 1. Generate Contig
        segments = get_contig_from_motif(motif_path)
        if not segments:
            print(f"    Skipping {motif_name} (no AA residues found)")
            continue
            
        # Strategy: [50-50/Seg1/15-25/Seg2/50-50]
        # Dynamic linker length?
        contig_parts = ["50-50"] 
        for seg in segments:
            contig_parts.append(seg)
            contig_parts.append("10-20") # Random linker
        contig_parts[-1] = "50-50"
        
        contig_str = "/".join(contig_parts)
        contig = f"['{contig_str}']"
        
        # 2. Run RFdiffusion
        rf_script = os.path.join(args.rf_path, "scripts/run_inference.py")
        
        cmd = [
            sys.executable, rf_script,
            f"inference.input_pdb={motif_path}",
            f"inference.output_prefix={os.path.join(rf_out, 'scaffold')}",
            f"inference.num_designs={args.num_designs}",
            f"contigmap.contigs={contig}"
        ]
        
        if args.dry_run:
            print(f"    Dry Run RFdiffusion: {' '.join(cmd)}")
        else:
            try:
                subprocess.run(cmd, check=True, cwd=project_root) # Run from ./ to allow relative paths if any
                
                # 3. Post-process: Add Metal
                scaffolds = glob.glob(os.path.join(rf_out, "*.pdb"))
                
                # Read metal lines
                metal_lines = []
                with open(motif_path, 'r') as f:
                    for line in f:
                        if "HETATM" in line or "ZN" in line or "CU" in line: # Simple filter, logic can be improved
                            # Check if it is the metal atom
                            if ion in line: 
                                metal_lines.append(line)
                
                for sc in scaffolds:
                    with open(sc, 'a') as f:
                        f.write("\n")
                        for line in metal_lines:
                            f.write(line)
                            
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
