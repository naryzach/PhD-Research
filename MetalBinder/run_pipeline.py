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
    # Changed default to None to detect user input
    parser.add_argument("--metals", default=None, help="Metals to process (Comma-separated). If unset, uses step-specific defaults.")
    parser.add_argument("--rf_path", default="../Tools/RFdiffusion2", help="Path to RFdiffusion2")
    parser.add_argument("--dry_run", action="store_true")
    
    parser.add_argument("--num_designs", type=int, default=2, help="Number of RFdiffusion designs per motif")
    
    # Download Specific
    parser.add_argument("--dataset", choices=["general", "curated", "mmp"], default="general", help="Dataset mode to download")

    
    # Graft Specific
    parser.add_argument("--template", default="../Local/Templates/human_VH3_IgG.pdb", help="Template PDB for grafting")
    parser.add_argument("--insert_at", default="93-102", help="Insertion residue/range (recommended: 93-102 for CDRH3)")
    parser.add_argument("--graft_chain", default="H", help="Template chain")
    
    # Finalize Specific
    parser.add_argument("--finalize_only", action="store_true", help="Deprecated: Use --finalize")
    parser.add_argument("--run_mpnn", action="store_true", help="Run ProteinMPNN on designs")
    parser.add_argument("--run_ligandmpnn", action="store_true", help="Run LigandMPNN on designs")
    # parser.add_argument("--select_best", action="store_true", help="Deprecated: Merged into run_allmetal3d")
    parser.add_argument("--run_allmetal3d", action="store_true", help="Run AllMetal3D (requires allmetal3d env)")
    parser.add_argument("--run_allmetal3d_ligand", action="store_true", help="Run AllMetal3D on LigandMPNN outputs (skips threading)")
    parser.add_argument("--water", action="store_true", help="Run Water3D prediction in AllMetal3D (default: False)")
    
    # Common/New Args
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")
    parser.add_argument("--linker_len", type=str, default="5-15", help="Graft linker range (default: 5-15)")
    parser.add_argument("--rigid_gaps", action="store_true", help="Disable flexibility in auto-filled gaps")
    parser.add_argument("--clash_threshold", type=float, default=3.5, help="Clash cleaning threshold in Angstroms (default 3.5)")
    parser.add_argument("--no_clash_cleanup", action="store_true", help="Disable automatic removal of clashing scaffold residues")
    
    args = parser.parse_args()
    
    base_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Logic for Metals:
    # If user provides --metals, use that for EVERYTHING.
    # If not, use defaults:
    #   Download -> Common (ZN, CU...)
    #   Catalog -> All (implicitly scans input)
    #   Scaffold -> ALL (scans catalog)
    #   Swap -> Default Target Set (ZN..CA)
    #   Graft -> ALL (scans catalog)
    
    user_metals = args.metals
    
    # 1. Download (Unchanged)
    if args.download or args.run_all:
        print("\n=== Stage 1: Download Datasets ===")
        # Default for download is the Common Set if not specified
        download_metals = user_metals if user_metals else "ZN,CU,NI,CO,MN,FE,MG,CA"
        
        cmd = [sys.executable, os.path.join(base_dir, "download_datasets.py"), 
               "--metals", download_metals, 
               "--limit", str(args.num_designs),
               "--mode", args.dataset]
        run_command(cmd, args.dry_run)
        
    # 2. Catalog (Unchanged)
    if args.catalog or args.run_all:
        print("\n=== Stage 2: Catalog Sites ===")
        # Catalog scans input dir, metals arg creates filter. 
        # If user_metals is None, we want to catalog EVERYTHING found.
        # But catalog_sites.py expects a metals list to *filter* for sites in the PDBs.
        # So we need a default list for cataloging.
        # Let's use expanded default or the user list.
        catalog_metals = user_metals if user_metals else "ZN,CU,NI,CO,MN,FE,MG,CA,LA,CE,ND,EU,GD,TB,DY,YB,LU,Y,U,URI,PU,TH,HG,CD,PB,AS"
        
        cmd = [sys.executable, os.path.join(base_dir, "catalog_sites.py"), 
               "--pdb_dir", "../Local/Metal_PDBs", 
               "--out_dir", "../Local/Metal_Catalog",
               "--metals", catalog_metals]
        run_command(cmd, args.dry_run)
        
    # 3. Scaffold (De Novo) - Updated for RFd2
    if args.scaffold or args.run_all:
        print("\n=== Stage 3: De Novo Scaffolding ===")
        # If user_metals is None -> "ALL" (scan catalog)
        scaffold_metals = user_metals if user_metals else "ALL"
        
        cmd = [sys.executable, os.path.join(base_dir, "run_scaffold.py"),
               "--metals", scaffold_metals,
               "--num_designs", str(args.num_designs),
               "--rf_path", args.rf_path]
        if args.dry_run: cmd.append("--dry_run")
        run_command(cmd, args.dry_run)
        
    # 4. Swap - Updated for RFd2
    if args.swap or args.run_all:
        print("\n=== Stage 4: Metal Swapping ===")
        # Swap targets. If user_metals is None -> Use run_swap defaults.
        cmd = [sys.executable, os.path.join(base_dir, "run_swap.py"),
               "--num_designs", str(args.num_designs),
               "--rf_path", args.rf_path]
        
        if user_metals:
             cmd.extend(["--target_metals", user_metals])
        
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
            if not os.path.exists(motif): continue
            
            # Output dir structure
            # ../Local/Metal_Predictions/Grafted/{Ion}/{Motif}
            # infer ion from path? motif is usually .../ZN/filename.pdb
            ion = os.path.basename(os.path.dirname(motif))
            
            # Check if ion is in requested metals (if user specified)
            if user_metals and ion not in user_metals.split(','):
                 continue
            
            out_dir = f"../Local/Metal_Predictions/Grafted/{ion}/{os.path.basename(motif).replace('.pdb','')}"
            
            cmd = [sys.executable, os.path.join(base_dir, "run_graft.py"),
                   "--motif", motif,
                   "--template", args.template,
                   "--insert_at", str(args.insert_at),
                   "--chain", args.graft_chain,
                   "--out_dir", out_dir,
                   "--num_designs", str(args.num_designs),
                   "--rf_path", args.rf_path,
                   "--linker_len", args.linker_len,
                   "--clash_threshold", str(args.clash_threshold)
                   ]
            if args.overwrite: cmd.append("--overwrite")
            if args.rigid_gaps: cmd.append("--rigid_gaps")
            if args.no_clash_cleanup: cmd.append("--no_clash_cleanup")
            if args.dry_run: cmd.append("--dry_run")
            run_command(cmd, args.dry_run)

    if args.run_mpnn: #or args.run_all:
        print("\n=== Stage 6: ProteinMPNN ===")
        cmd = [sys.executable, os.path.join(base_dir, "run_mpnn.py"),
               "--pred_dir", "../Local/Metal_Predictions",
               "--pmpnn_path", "../Tools/ProteinMPNN"]
        if args.overwrite: cmd.append("--overwrite")
        run_command(cmd, args.dry_run)

    # 7. LigandMPNN (New Step)
    if args.run_ligandmpnn:
        print("\n=== Stage 7: LigandMPNN ===")
        cmd = [sys.executable, os.path.join(base_dir, "run_ligandmpnn.py"),
               "--pred_dir", "../Local/Metal_Predictions",
               "--lmpnn_path", "../Tools/LigandMPNN"]
        if args.overwrite: cmd.append("--overwrite")
        run_command(cmd, args.dry_run)
    
    # 8. AllMetal3D (New Step)
    if args.run_allmetal3d:
        print("\n=== Stage 8: AllMetal3D ===")
        
        # Check if we can import allmetal3d in the current environment
        import importlib.util
        if importlib.util.find_spec("allmetal3d") is None:
            print("Error: 'allmetal3d' package not found in current environment.")
            print("Please activate the environment before running this step:")
            print("  conda activate allmetal3d")
            print("  python MetalBinder/run_pipeline.py --run_allmetal3d")
            sys.exit(1)
            
        cmd = [sys.executable, os.path.join(base_dir, "run_allmetal3d.py"),
               "--pred_dir", "../Local/Metal_Predictions"]
        if args.overwrite: cmd.append("--overwrite")
        if args.water: cmd.append("--water")
        run_command(cmd, args.dry_run)

    # 9. AllMetal3D (LigandMPNN Source)
    if args.run_allmetal3d_ligand: #or args.run_all:
        print("\n=== Stage 8b: AllMetal3D (LigandMPNN Source) ===")
        # Check if we can import allmetal3d in the current environment
        import importlib.util
        if importlib.util.find_spec("allmetal3d") is None:
            print("Error: 'allmetal3d' package not found in current environment.")
            print("Please activate the environment before running this step:")
            print("  conda activate allmetal3d")
            print("  python MetalBinder/run_pipeline.py --run_allmetal3d_ligand")
            sys.exit(1)
            
        cmd = [sys.executable, os.path.join(base_dir, "run_allmetal3d_ligand.py"),
               "--pred_dir", "../Local/Metal_Predictions"]
        if args.overwrite: cmd.append("--overwrite")
        if args.water: cmd.append("--water")
        run_command(cmd, args.dry_run)

    # 9. Finalize
    if args.finalize or args.run_all:
        print("\n=== Stage 9: Finalization ===")
        cmd = [sys.executable, os.path.join(base_dir, "finalize_designs.py"),
               "--pred_dir", "../Local/Metal_Predictions",
               "--out_csv", "../Local/Metal_Predictions/design_summary.csv",
               "--out_json", "../Local/Metal_Predictions/alphafold_inputs.json"]
        run_command(cmd, args.dry_run)

    print("\nPipeline Complete.")

if __name__ == "__main__":
    main()
