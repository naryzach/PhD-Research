import os
import glob
import pandas as pd
import sys
import argparse
import warnings

# Suppress warnings from torch/gradio/etc
warnings.filterwarnings("ignore")

# Import AllMetal3D
try:
    from allmetal3d.utils.main import predict_cli
except ImportError:
    print("Error: Could not import allmetal3d. Please run this script within the allmetal3d conda environment.")
    sys.exit(1)

# Labels as defined in allmetal3d/utils/main.py
LABELS = ['Alkali', 'MG', 'CA', 'ZN', 'NonZNTM', 'NoMetal']

def main():
    parser = argparse.ArgumentParser(description="Batch run AllMetal3D and summarize results with detailed probabilities.")
    parser.add_argument("--input_dir", default="../Local/AllMetal_input", help="Directory containing input PDB files")
    parser.add_argument("--output_dir", default="../Local/Allmetal_output", help="Directory to save outputs")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing runs")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_dir):
        print(f"Error: Input directory {args.input_dir} does not exist.")
        return

    os.makedirs(args.output_dir, exist_ok=True)
    
    pdb_files = glob.glob(os.path.join(args.input_dir, "*.pdb"))
    if not pdb_files:
        print(f"No PDB files found in {args.input_dir}")
        return
        
    print(f"Found {len(pdb_files)} PDB files to process.")
    
    summary_data = []

    for pdb_file in pdb_files:
        basename = os.path.splitext(os.path.basename(pdb_file))[0]
        pdb_out_dir = os.path.join(args.output_dir, basename)
        
        # Check if already run (check for _metals.pdb)
        metals_pdb = os.path.join(pdb_out_dir, f"{basename}_metals.pdb")
        
        # If overwrite is False and output exists, we might want to skip execution but RE-READ the results?
        # Re-reading detailed probabilities from PDB is hard because they aren't stored there (only top choice).
        # So we MUST re-run prediction to get granular data if we don't have it saved elsewhere.
        # Since the user specifically asked for granularity, let's assume we need to run it unless we save a separate JSON/CSV per run.
        # For now, let's just re-run if overwrite is True OR if we want to update the summary.
        # Actually, let's force run if overwrite is True, otherwise if output exists, we skip it. 
        # BUT if we skip, we can't get the granular data for the summary unless we stored it.
        # The previous version didn't store it. So for THIS run, we should probably process everything.
        # However, to save time on subsequent runs, we should save the detailed results to a file.
        
        csv_out_path = os.path.join(pdb_out_dir, f"{basename}_detailed.csv")
        
        if os.path.exists(csv_out_path) and not args.overwrite:
            print(f"Skipping {basename} (Detailed output exists, reading from file)")
            # Read back the detailed CSV to add to summary
            try:
                df = pd.read_csv(csv_out_path)
                summary_data.extend(df.to_dict('records'))
            except Exception as e:
                print(f"Error reading {csv_out_path}: {e}")
            continue

        print(f"Processing {basename}...")
        os.makedirs(pdb_out_dir, exist_ok=True)
        
        try:
            # predict_cli arguments: 
            # pdb, models, pthreshold=0.1, threshold=7, batch_size=20, mode="fast", central_residue=None, radius=8, output_dir = ".", backend="multiprocessing"
            # It returns a tuple. results is at index 5.
            
            # Use 'allmetal3d' model string
            ret = predict_cli(
                pdb=pdb_file, 
                models="allmetal3d", 
                output_dir=pdb_out_dir
            )
            
            # ret structure: (probefile, cubefile, probefile_water, cubefile_water, results)
            # predict_cli returns 5 elements (lines 609 in main.py)
            results = ret[4]
            
            if not results:
                print(f"  No sites found for {basename}.")
                continue

            row_list = []
            for res in results:
                # res is a dict:
                # {
                #   "index": int,
                #   "location_confidence": float,
                #   "probabilities_identity": [float, ...], # Corresponds to LABELS
                #   "probabilities_geometry": [float, ...],
                #   "close_residues": [...]
                # }
                
                probs = res["probabilities_identity"]
                
                row = {
                    "Structure": basename,
                    "Site_Index": res["index"],
                    "Location_Prob": res["location_confidence"],
                    "File": metals_pdb # Path where the PDB representation is
                }
                
                # Add individual probabilities
                for label, p in zip(LABELS, probs):
                    row[label] = p
                
                # Identify the "Top Identity" which matches the max probability
                # Note: The code logic uses argmax of the LOGITS/Output to define identity, 
                # but probabilities_identity is the softmaxed probability? 
                # Let's assume the highest probability here corresponds to the identity.
                # Actually, `predict_identity` code:
                # identity = f"{l_metal[o[0][i].argmax()]} ..."
                # probabilities_identity = [round(x,2) for x in o[0][row.Site].tolist()]
                # checks out.
                
                # Let's also add the Top Identity explicitly
                max_prob_idx = max(range(len(probs)), key=probs.__getitem__)
                row["Top_Identity"] = LABELS[max_prob_idx]
                
                row_list.append(row)
                summary_data.append(row)
            
            # Save individual detailed CSV
            if row_list:
                df_detailed = pd.DataFrame(row_list)
                df_detailed.to_csv(csv_out_path, index=False)
                print(f"  Saved detailed results to {csv_out_path}")

        except Exception as e:
            print(f"  Error processing {basename}: {e}")
            import traceback
            traceback.print_exc()

    # Generate Summary
    if summary_data:
        df = pd.DataFrame(summary_data)
        
        # Sort by Structure then Location Prob
        df = df.sort_values(by=["Structure", "Location_Prob"], ascending=[True, False])
        
        # Reorder columns slightly for readability
        cols = ["Structure", "Site_Index", "Location_Prob", "Top_Identity"] + LABELS + ["File"]
        # Ensure all cols exist
        cols = [c for c in cols if c in df.columns]
        df = df[cols]
        
        output_csv = os.path.join(args.output_dir, "binding_summary.csv")
        df.to_csv(output_csv, index=False)
        print(f"\nSummary saved to {output_csv}")
        print(df.head())
    else:
        print("No data extracted for summary.")

if __name__ == "__main__":
    main()
