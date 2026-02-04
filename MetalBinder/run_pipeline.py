import argparse
import os
import sys
import subprocess

def run_command(cmd, dry_run=False):
    if dry_run:
        print(f"[Dry Run] {' '.join(cmd)}")
    else:
        print(f"[Exec] {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error executing command: {e}")
            sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="MetalBinder v2.0 Orchestrator")
    
    # Workflow Steps
    parser.add_argument("--download", action="store_true", help="Run Dataset Downloader")
    parser.add_argument("--catalog", action="store_true", help="Run Site Cataloging")
    parser.add_argument("--scaffold", action="store_true", help="Run De Novo Scaffolding")
    parser.add_argument("--swap", action="store_true", help="Run Metal Swapping")
    parser.add_argument("--graft", action="store_true", help="Run Antibody Grafting")
    parser.add_argument("--finalize", action="store_true", help="Run Finalization (MPNN + AF2 CSV)")
    
    parser.add_argument("--run_all", action="store_true", help="Run all steps (except Graft which needs specific args)")
    
    # Common Args
    parser.add_argument("--metals", default="ZN,CU,NI,CO,MN,FE,MG,CA", help="Metals to process")
    parser.add_argument("--dry_run", action="store_true")
    
    parser.add_argument("--num_designs", type=int, default=2, help="Number of RFdiffusion designs per motif")
    
    # Graft Specific
    parser.add_argument("--template", help="Template PDB for grafting")
    parser.add_argument("--insert_at", type=int, help="Insertion residue for grafting")
    parser.add_argument("--graft_chain", default="H", help="Template chain")
    
    # Finalize Specific
    parser.add_argument("--run_mpnn", action="store_true", help="Enable ProteinMPNN in finalization")
    
    args = parser.parse_args()
    
    base_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 1. Download
    if args.download or args.run_all:
        print("\n=== Stage 1: Download Datasets ===")
        cmd = [sys.executable, os.path.join(base_dir, "download_datasets.py"), "--metals", args.metals]
        run_command(cmd, args.dry_run)
        
    # 2. Catalog
    if args.catalog or args.run_all:
        print("\n=== Stage 2: Catalog Sites ===")
        cmd = [sys.executable, os.path.join(base_dir, "catalog_sites.py"), 
               "--pdb_dir", "../Local/Metal_PDBs", 
               "--out_dir", "../Local/Metal_Catalog",
               "--metals", args.metals]
        run_command(cmd, args.dry_run)
        
    # 3. Scaffold (De Novo)
    if args.scaffold or args.run_all:
        print("\n=== Stage 3: De Novo Scaffolding ===")
        cmd = [sys.executable, os.path.join(base_dir, "run_scaffold.py"),
               "--metals", args.metals,
               "--num_designs", str(args.num_designs)]
        if args.dry_run: cmd.append("--dry_run")
        run_command(cmd, args.dry_run)
        
    # 4. Swap
    if args.swap or args.run_all:
        print("\n=== Stage 4: Metal Swapping ===")
        cmd = [sys.executable, os.path.join(base_dir, "run_swap.py"),
               "--target_metals", args.metals,
               "--num_designs", str(args.num_designs)]
        if args.dry_run: cmd.append("--dry_run")
        run_command(cmd, args.dry_run)
        
    # 5. Graft
    if args.graft:
        print("\n=== Stage 5: Antibody Grafting ===")
        if not args.template or not args.insert_at:
            print("Error: --template and --insert_at are required for grafting.")
        else:
            # We need to iterate over catalog motifs... run_graft processes ONE motif.
            # We need a wrapper loop here? run_graft.py currently takes single --motif argument.
            # I should update run_graft.py to handle mass processing or loop here.
            # Looping here is easier for now given the design.
            import glob
            import pandas as pd
            
            # Read catalog to get motifs
            try:
                df = pd.read_csv("../Local/Metal_Catalog/master_catalog.csv")
                motifs = df['path'].tolist()
            except:
                print("Catalog not found. Skipping graft loop.")
                motifs = []
            
            # Filter motifs by metal?
            # Creating output dir: ../Local/Metal_Predictions/Grafted/{Ion}/{Motif}
            
            for motif in motifs:
                # Check metal
                # TODO: Filter by args.metals
                cmd = [sys.executable, os.path.join(base_dir, "run_graft.py"),
                       "--motif", motif,
                       "--template", args.template,
                       "--insert_at", str(args.insert_at),
                       "--chain", args.graft_chain,
                       "--out_dir", f"../Local/Metal_Predictions/Grafted/{os.path.basename(os.path.dirname(motif))}/{os.path.basename(motif).replace('.pdb','')}",
                       "--num_designs", str(args.num_designs)
                       ]
                if args.dry_run: cmd.append("--dry_run")
                run_command(cmd, args.dry_run)

    # 6. Finalize
    if args.finalize or args.run_all:
        print("\n=== Stage 6: Finalization ===")
        cmd = [sys.executable, os.path.join(base_dir, "finalize_designs.py")]
        if args.run_mpnn: cmd.append("--run_mpnn")
        run_command(cmd, args.dry_run)

    print("\nPipeline Complete.")

if __name__ == "__main__":
    main()
