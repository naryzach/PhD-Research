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
    parser = argparse.ArgumentParser(description="MetalBinder v2.0 Orchestrator (RFdiffusion2)")
    
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
    parser.add_argument("--rf_path", default="../Tools/RFdiffusion2", help="Path to RFdiffusion2")
    parser.add_argument("--dry_run", action="store_true")
    
    parser.add_argument("--num_designs", type=int, default=2, help="Number of RFdiffusion designs per motif")
    
    # Graft Specific
    parser.add_argument("--template", default="../Local/Templates/human_VH3_IgG.pdb", help="Template PDB for grafting")
    parser.add_argument("--insert_at", default="97-108", help="Insertion residue/range (e.g. 97-108)")
    parser.add_argument("--graft_chain", default="H", help="Template chain")
    
    # Finalize Specific
    parser.add_argument("--finalize_only", action="store_true", help="Deprecated: Use --finalize")
    parser.add_argument("--run_mpnn", action="store_true", help="Enable ProteinMPNN in finalization")
    
    args = parser.parse_args()
    
    base_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 1. Download (Unchanged)
    if args.download or args.run_all:
        print("\n=== Stage 1: Download Datasets ===")
        cmd = [sys.executable, os.path.join(base_dir, "download_datasets.py"), "--metals", args.metals, "--limit", str(args.num_designs)]
        run_command(cmd, args.dry_run)
        
    # 2. Catalog (Unchanged)
    if args.catalog or args.run_all:
        print("\n=== Stage 2: Catalog Sites ===")
        cmd = [sys.executable, os.path.join(base_dir, "catalog_sites.py"), 
               "--pdb_dir", "../Local/Metal_PDBs", 
               "--out_dir", "../Local/Metal_Catalog",
               "--metals", args.metals]
        run_command(cmd, args.dry_run)
        
    # 3. Scaffold (De Novo) - Updated for RFd2
    if args.scaffold or args.run_all:
        print("\n=== Stage 3: De Novo Scaffolding ===")
        cmd = [sys.executable, os.path.join(base_dir, "run_scaffold.py"),
               "--metals", args.metals,
               "--num_designs", str(args.num_designs),
               "--rf_path", args.rf_path]
        if args.dry_run: cmd.append("--dry_run")
        run_command(cmd, args.dry_run)
        
    # 4. Swap - Updated for RFd2
    if args.swap or args.run_all:
        print("\n=== Stage 4: Metal Swapping ===")
        cmd = [sys.executable, os.path.join(base_dir, "run_swap.py"),
               "--target_metals", args.metals,
               "--num_designs", str(args.num_designs),
               "--rf_path", args.rf_path]
        if args.dry_run: cmd.append("--dry_run")
        run_command(cmd, args.dry_run)
        
    # 5. Graft - Updated for RFd2
    if args.graft:
        print("\n=== Stage 5: Antibody Grafting ===")
        import glob
        import pandas as pd
        
        # Read catalog to get motifs
        try:
            # Assuming catalog location
            catalog_path = "../Local/Metal_Catalog/master_catalog.csv"
            if os.path.exists(catalog_path):
                df = pd.read_csv(catalog_path)
                motifs = df['path'].tolist()
            else:
                motifs = []
                print(f"Warning: Catalog not found at {catalog_path}")
        except Exception as e:
            print(f"Error reading catalog: {e}")
            motifs = []
        
        for motif in motifs:
            # Filter by metal? args.metals?
            # Ideally yes, but catalog usually has all.
            # Simplification: checks if motif path exists
            if not os.path.exists(motif): continue
            
            # Output dir structure
            # ../Local/Metal_Predictions/Grafted/{Ion}/{Motif}
            # infer ion from path? motif is usually .../ZN/filename.pdb
            ion = os.path.basename(os.path.dirname(motif))
            
            # Check if ion is in requested metals
            if ion not in args.metals.split(','):
                 continue
            
            out_dir = f"../Local/Metal_Predictions/Grafted/{ion}/{os.path.basename(motif).replace('.pdb','')}"
            
            cmd = [sys.executable, os.path.join(base_dir, "run_graft.py"),
                   "--motif", motif,
                   "--template", args.template,
                   "--insert_at", str(args.insert_at),
                   "--chain", args.graft_chain,
                   "--out_dir", out_dir,
                   "--num_designs", str(args.num_designs),
                   "--rf_path", args.rf_path
                   ]
            if args.dry_run: cmd.append("--dry_run")
            run_command(cmd, args.dry_run)

    # 6. Finalize (Unchanged)
    if args.finalize or args.run_all:
        print("\n=== Stage 6: Finalization ===")
        cmd = [sys.executable, os.path.join(base_dir, "finalize_designs.py")]
        if args.run_mpnn: cmd.append("--run_mpnn")
        run_command(cmd, args.dry_run)

    print("\nPipeline Complete.")

if __name__ == "__main__":
    main()
