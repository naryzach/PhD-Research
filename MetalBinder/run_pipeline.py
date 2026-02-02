import argparse
import os
import subprocess
import json
import glob
import sys

def run_pipeline(args):
    # Absolute paths
    base_dir = os.path.dirname(os.path.abspath(__file__))
    input_path = os.path.abspath(args.input)
    out_dir = os.path.abspath(args.out)
    
    os.makedirs(out_dir, exist_ok=True)
    
    # Check if input is file or directory
    if os.path.isdir(input_path):
        pdb_files = glob.glob(os.path.join(input_path, "*.pdb"))
    else:
        pdb_files = [input_path]
        
    print(f"Starting MetalBinder Pipeline on {len(pdb_files)} files using {sys.executable}.")
    
    for pdb_file in pdb_files:
        basename = os.path.basename(pdb_file).replace(".pdb", "")
        file_out_dir = os.path.join(out_dir, basename)
        os.makedirs(file_out_dir, exist_ok=True)
        
        print(f"\nProcessing {basename}...")
        
        # 1. Find Metal Sites
        sites_json = os.path.join(file_out_dir, "sites.json")
        pymol_scr = os.path.join(file_out_dir, "view.pml")
        
        cmd_find = [
            sys.executable, os.path.join(base_dir, "find_metal_sites.py"),
            pdb_file,
            "--out_json", sites_json,
            "--out_pymol", pymol_scr
        ]
        
        try:
            subprocess.run(cmd_find, check=True)
        except Exception as e:
            print(f"Failed to find sites for {basename}: {e}")
            continue
            
        # Check if sites found
        if not os.path.exists(sites_json):
            print("No sites JSON created.")
            continue
            
        with open(sites_json, 'r') as f:
            sites = json.load(f)
            
        if not sites:
            print(f"No metal sites found in {basename}.")
            continue
            
        # 2. Swap Metals (If requested)
        if args.swap_metals:
            print(f"Running Metal Swap for {args.swap_metals}...")
            cmd_swap = [
                sys.executable, os.path.join(base_dir, "metal_swap.py"),
                "--pdb", pdb_file,
                "--sites_json", sites_json,
                "--metals", args.swap_metals,
                "--out", os.path.join(file_out_dir, "swap"),
                "--num_designs", str(args.num_designs)
            ]
            
            if args.dry_run:
                print("Dry Run Swap:", " ".join(cmd_swap))
            else:
                try:
                    subprocess.run(cmd_swap, check=True)
                except Exception as e:
                    print(f"Swap failed: {e}")

    # 3. Catalog (If input was directory)
    motifs_dir = None
    if os.path.isdir(input_path) and args.catalog:
        print("\nBuilding Catalog...")
        catalog_csv = os.path.join(out_dir, "catalog.csv")
        motifs_dir = os.path.join(out_dir, "motifs")
        
        cmd_cat = [
            sys.executable, os.path.join(base_dir, "build_catalog.py"),
            input_path,
            catalog_csv,
            "--scaffold_dir", motifs_dir
        ]
        subprocess.run(cmd_cat)
        
        # 4. Scaffold Motifs (If requested)
        if args.scaffold_motifs and motifs_dir and os.path.exists(motifs_dir):
            print(f"\nScaffolding Motifs in {motifs_dir}...")
            motif_files = glob.glob(os.path.join(motifs_dir, "*.pdb"))
            
            for motif_file in motif_files:
                print(f"Scaffolding {os.path.basename(motif_file)}...")
                scaffold_out = os.path.join(file_out_dir, "scaffolds", os.path.basename(motif_file).replace(".pdb", ""))
                
                cmd_scaffold = [
                    sys.executable, os.path.join(base_dir, "motif_scaffold.py"),
                    "--motif_pdb", motif_file,
                    "--out", scaffold_out,
                    "--num_designs", str(args.num_designs)
                ]
                
                if args.dry_run:
                     cmd_scaffold.append("--dry_run")
                     print("Dry Run Scaffold:", " ".join(cmd_scaffold))
                     # We execute it even in dry run so the script prints its own dry output
                     subprocess.run(cmd_scaffold)
                else:
                    try:
                        subprocess.run(cmd_scaffold, check=True)
                    except Exception as e:
                        print(f"Scaffolding failed for {motif_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description="MetalBinder Master Pipeline")
    parser.add_argument("--input", required=True, help="Input PDB file or Directory")
    parser.add_argument("--out", default="../Local/MetalBinder_Results", help="Output directory")
    
    parser.add_argument("--swap_metals", help="Comma-separated list of metals to swap/design for (e.g. ZN,CU)")
    parser.add_argument("--catalog", action="store_true", help="Build a catalog (if input is directory)")
    parser.add_argument("--scaffold_motifs", action="store_true", help="Scaffold all motifs found in the cataloging step")
    
    parser.add_argument("--num_designs", type=int, default=1)
    parser.add_argument("--dry_run", action="store_true")
    
    args = parser.parse_args()
    run_pipeline(args)

if __name__ == "__main__":
    main()
