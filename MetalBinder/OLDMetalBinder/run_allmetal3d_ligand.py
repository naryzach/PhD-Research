import argparse
import os
import sys
import glob
import json
import csv
import pandas as pd

# Import AllMetal3D
try:
    from allmetal3d import predict
except ImportError:
    predict = None
    pass

def main():
    parser = argparse.ArgumentParser(description="Run AllMetal3D on LigandMPNN Packed Outputs")
    parser.add_argument("--pred_dir", default="../Local/Metal_Predictions", help="Predictions Directory")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")
    parser.add_argument("--water", action="store_true", help="Run Water3D prediction (default: False)")
    parser.add_argument("--dry_run", action="store_true", help="Print command without executing")
    args = parser.parse_args()

    if predict is None and not args.dry_run:
        print("Error: Could not import allmetal3d. Make sure you are running in the 'allmetal3d' environment.")
        sys.exit(1)

    print(f"Scanning {args.pred_dir} for LigandMPNN packed structures...")
    
    # Path pattern: .../ligandmpnn/packed/*.pdb
    packed_files = glob.glob(os.path.join(args.pred_dir, "**", "ligandmpnn", "packed", "*.pdb"), recursive=True)
    
    print(f"Found {len(packed_files)} packed structures.")

    for i, pdb_path in enumerate(packed_files):
        # Path parsing for logging
        # .../Category/Ion/Motif/ligandmpnn/packed/file.pdb
        packed_dir = os.path.dirname(pdb_path)
        ligandmpnn_dir = os.path.dirname(packed_dir)
        motif_dir = os.path.dirname(ligandmpnn_dir)
        
        motif_name = os.path.basename(motif_dir)
        ion_dir = os.path.dirname(motif_dir)
        ion = os.path.basename(ion_dir)
        category_dir = os.path.dirname(ion_dir)
        category = os.path.basename(category_dir)
        
        basename = os.path.basename(pdb_path).replace(".pdb", "")
        
        # Output directory: .../Motif/allmetal3d_ligand
        # Move up from ligandmpnn/packed to Motif level
        allmetal_out_dir = os.path.join(motif_dir, "allmetal3d_ligand")
        os.makedirs(allmetal_out_dir, exist_ok=True)
        
        # Check if already processed
        # AllMetal3D output naming: {basename}_results.json
        expected_json = os.path.join(allmetal_out_dir, f"{basename}_results.json")
        
        # Granular skip logic (same as run_allmetal3d.py)
        if os.path.exists(expected_json) and not args.overwrite:
            if args.water:
                # Check if water prediction is missing
                try:
                    with open(expected_json, 'r') as f:
                         data = json.load(f)
                    if "water_probabilities" not in data:
                         pass # Needs update
                    else:
                         print(f"[{i+1}/{len(packed_files)}] [{category}] [{ion}] [{motif_name}] Skipping {basename} (Already exists)")
                         continue
                except:
                     pass # Corrupt, re-run
            else:
                 print(f"[{i+1}/{len(packed_files)}] [{category}] [{ion}] [{motif_name}] Skipping {basename} (Already exists)")
                 continue

        print(f"[{i+1}/{len(packed_files)}] [{category}] [{ion}] [{motif_name}] Processing {basename}")
        
        if args.dry_run:
            continue
            
        try:
            model_mode = "all" if args.water else "allmetal3d"
            
            # Predict directly on the packed PDB
            paths, results = predict(
                pdb_path,
                output_dir=allmetal_out_dir,
                models=model_mode
            )
            
            # Save results to JSON
            with open(expected_json, 'w') as f:
                json.dump(results, f, indent=4)
                
            # Create summary CSV similar to original script
            csv_path = os.path.join(allmetal_out_dir, f"{basename}_probs.csv")
            
            # Flatten results for CSV
            rows = []
            for probefile_key, res_data in results.items():
                 for residue, probs in res_data.items():
                      row = {"residue": residue}
                      row.update(probs)
                      rows.append(row)
            
            if rows:
                pd.DataFrame(rows).to_csv(csv_path, index=False)
                
                # Print best match
                df = pd.DataFrame(rows)
                # Drop residue col
                ion_cols = [c for c in df.columns if c != "residue"]
                if ion_cols:
                    max_val = df[ion_cols].max().max()
                    print(f"  -> Max Probability: {max_val:.3f}")

        except Exception as e:
            print(f"  Failed to run AllMetal3D on {basename}: {e}")
            import traceback
            traceback.print_exc()

if __name__ == "__main__":
    main()
