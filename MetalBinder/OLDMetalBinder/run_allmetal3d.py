import argparse
import os
import glob
import re
import subprocess
import sys

# Import AllMetal3D
# Must process this import carefully as it relies on PyTorch/CUDA
try:
    from allmetal3d import predict
except ImportError:
    # If run outside the environment, this might fail.
    # We will handle it, but for now assuming correct env.
    predict = None
    pass

# Import thread_sequence
# Assuming MetalBinder is in python path or we append it
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from thread_sequence import thread_sequence

def main():
    parser = argparse.ArgumentParser(description="Run AllMetal3D on Best Designs")
    parser.add_argument("--pred_dir", default="../Local/Metal_Predictions", help="Predictions Directory")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")
    parser.add_argument("--water", action="store_true", help="Run Water3D prediction (default: False)")
    args = parser.parse_args()
    
    if predict is None:
        print("Error: Could not import allmetal3d. Make sure you are running in the 'allmetal3d' environment.")
        sys.exit(1)
    
    print(f"Scanning {args.pred_dir} for ProteinMPNN outputs...")
    
    # Find all sequence files
    # pattern: .../proteinmpnn/seqs/*.fa
    seq_files = glob.glob(os.path.join(args.pred_dir, "**", "proteinmpnn", "seqs", "*.fa"), recursive=True)
    
    print(f"Found {len(seq_files)} MPNN output files.")
    
    for fasta_file in seq_files:
        # Identify original PDB
        # fasta: .../Type/Ion/Motif/proteinmpnn/seqs/design_name.fa
        # original: .../Type/Ion/Motif/rfdiffusion/design_name.pdb
        
        basename = os.path.basename(fasta_file).replace(".fa", "").replace("_mpnn_input", "")
        dirname = os.path.dirname(os.path.dirname(os.path.dirname(fasta_file)))
        rfdiffusion_pdb = os.path.join(dirname, "rfdiffusion", f"{basename}.pdb")
        
        if not os.path.exists(rfdiffusion_pdb):
            continue
            
        # Parse FASTA for best score
        best_score = float('inf')
        best_seq = None
        
        try:
            with open(fasta_file, 'r') as f:
                lines = f.readlines()
            
            current_header = ""
            for line in lines:
                line = line.strip()
                if line.startswith(">"):
                    current_header = line
                else:
                    if not line: continue
                    # Parse score
                    score_match = re.search(r"score=([-\d\.]+)", current_header)
                    if score_match:
                        score = float(score_match.group(1))
                        if score < best_score:
                            best_score = score
                            best_seq = line
        except Exception as e:
            print(f"Error parsing {fasta_file}: {e}")
            continue

        if not best_seq:
            print(f"No valid sequence found in {fasta_file}")
            continue

        # Prepare AllMetal3D paths
        allmetal_dir = os.path.join(dirname, "allmetal3d")
        os.makedirs(allmetal_dir, exist_ok=True)
        
        threaded_pdb = os.path.join(allmetal_dir, f"{basename}_threaded.pdb")
        allmetal_out_dir = os.path.join(allmetal_dir, f"{basename}_am3d")
        os.makedirs(allmetal_out_dir, exist_ok=True)
        
        print(f"Processing {basename}...")
        
        # 1. Thread Sequence
        try:
            if not os.path.exists(threaded_pdb) or args.overwrite:
                print(f"  Threading sequence (Score: {best_score:.2f})...")
                thread_sequence(rfdiffusion_pdb, threaded_pdb, best_seq)
            else:
                print(f"  Using existing threaded structure.")
        except Exception as e:
            print(f"  Error threading {basename}: {e}")
            continue

        # 2. Run AllMetal3D
        try:
            # Check if prediction is already done
            expected_results = os.path.join(allmetal_out_dir, f"{basename}_results.json")
            skip_prediction = False
            
            if os.path.exists(expected_results) and not args.overwrite:
                # Load existing results to check consistency (e.g. water vs no-water)
                try:
                    import json
                    with open(expected_results, 'r') as f:
                        existing_data = json.load(f)
                    
                    # If user wants water, but existing data says water_predicted=False, must RERUN.
                    existing_water = existing_data.get("water_predicted", False)
                    
                    if args.water and not existing_water:
                        print(f"  Existing results lack water prediction. Rerunning...")
                        skip_prediction = False
                    else:
                        print(f"  Skipping AllMetal3D (Results exist).")
                        skip_prediction = True
                        
                except Exception as e:
                    print(f"  Error reading existing results, will rerun: {e}")
                    skip_prediction = False
            
            if not skip_prediction:
                model_mode = "all" if args.water else "allmetal3d"
                print(f"  Running AllMetal3D API (models={model_mode})...")
                
                # Predict
                # Returns: probefile, cubefile, probefile_water, cubefile_water, results
                probefile, cubefile, probefile_water, cubefile_water, results = predict(
                    threaded_pdb,
                    output_dir=allmetal_out_dir,
                    models=model_mode
                )
                
                results_wrapper = {
                    "results": results,
                    "dataset": basename,
                    "water_predicted": args.water
                }
    
                with open(expected_results, 'w') as f:
                    import json
                    json.dump(results_wrapper, f, indent=2, default=str)
                    
                print(f"  AllMetal3D completed. Results saved to {expected_results}")

        except Exception as e:
             print(f"  Error running AllMetal3D API for {basename}: {e}")
             import traceback
             traceback.print_exc()

if __name__ == "__main__":
    main()
