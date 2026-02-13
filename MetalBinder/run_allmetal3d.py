import argparse
import os
import glob
import re
import subprocess
import sys

# Import thread_sequence
# Assuming MetalBinder is in python path or we append it
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from thread_sequence import thread_sequence

def main():
    parser = argparse.ArgumentParser(description="Run AllMetal3D on Best Designs")
    parser.add_argument("--pred_dir", default="../Local/Metal_Predictions", help="Predictions Directory")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")
    args = parser.parse_args()
    
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
        # Target: .../Type/Ion/Motif/allmetal3d/name_threaded.pdb
        
        allmetal_dir = os.path.join(dirname, "allmetal3d")
        os.makedirs(allmetal_dir, exist_ok=True)
        
        threaded_pdb = os.path.join(allmetal_dir, f"{basename}_threaded.pdb")
        allmetal_out_dir = os.path.join(allmetal_dir, f"{basename}_am3d")
        os.makedirs(allmetal_out_dir, exist_ok=True)
        
        # Check if done
        if os.path.exists(threaded_pdb) and os.path.exists(allmetal_out_dir) and not args.overwrite:
            # Check if directory is empty? 
            if os.listdir(allmetal_out_dir):
                print(f"Skipping {basename} (Already done)")
                continue

        print(f"Processing {basename}...")
        
        # 1. Thread Sequence (Build Sidechains)
        try:
            print(f"  Threading sequence (Score: {best_score:.2f})...")
            thread_sequence(rfdiffusion_pdb, threaded_pdb, best_seq)
        except Exception as e:
            print(f"  Error threading {basename}: {e}")
            continue
            
        # 2. Run AllMetal3D
        try:
            print(f"  Running AllMetal3D...")
            cmd = [
                "allmetal3d",
                "-i", threaded_pdb,
                "-o", allmetal_out_dir
            ]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print(f"  AllMetal3D completed.")
        except subprocess.CalledProcessError:
             print(f"  AllMetal3D failed for {basename}")
        except Exception as e:
             print(f"  Error running AllMetal3D for {basename}: {e}")

if __name__ == "__main__":
    main()
