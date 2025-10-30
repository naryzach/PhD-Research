#!/usr/bin/env python

"""
Runs ProteinMPNN to design a specific motif on a ligand chain
while in complex with a target.

This script automates the creation of the '--fixed_positions_jsonl'
file based on a user-provided motif range.
"""

import argparse
import json
import subprocess
import tempfile
import os
import sys
from typing import Set, Dict, List

def get_ligand_residues(pdb_path: str, chain_id: str) -> Set[int]:
    """
    Parses a PDB file to get all residue numbers for a specific chain.
    
    Args:
        pdb_path: Path to the PDB file.
        chain_id: The chain ID to parse.
    
    Returns:
        A set of integer residue numbers found for that chain.
    """
    residues: Set[int] = set()
    
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                # Check for ATOM records belonging to the specified chain
                if line.startswith("ATOM") and line[21] == chain_id:
                    try:
                        # PDB format: residue number is in columns 23-26
                        res_num = int(line[22:26].strip())
                        residues.add(res_num)
                    except ValueError:
                        # Skip lines with improperly formatted residue numbers
                        continue
    except FileNotFoundError:
        print(f"Error: PDB file not found at {pdb_path}", file=sys.stderr)
        sys.exit(1)
    
    if not residues:
        print(f"Error: No residues found for chain '{chain_id}' in {pdb_path}.", file=sys.stderr)
        print("Please check your --ligand_chain ID and PDB file.", file=sys.stderr)
        sys.exit(1)
        
    return residues

def run_motif_design(args: argparse.Namespace):
    """
    Generates the fixed positions JSONL file and runs the ProteinMPNN command.
    """
    
    # 1. Get all residues for the ligand chain from the PDB
    print(f"Parsing PDB file: {args.pdb_file} for chain {args.ligand_chain}...")
    all_residues = get_ligand_residues(args.pdb_file, args.ligand_chain)
    print(f"Found {len(all_residues)} total residues for chain {args.ligand_chain}.")

    # 2. Determine fixed vs. designable residues
    design_start, design_end = args.design_motif
    
    # Create a set of all residues in the designable motif range
    designable_residues = set(range(design_start, design_end + 1))
    
    # Validate that the design motif is actually part of the PDB chain
    if not designable_residues.issubset(all_residues):
        missing = designable_residues - all_residues
        print(f"Warning: The design motif {design_start}-{design_end} contains residues "
              f"not found in chain {args.ligand_chain}: {sorted(list(missing))}", file=sys.stderr)
        
        # Filter designable residues to only those that actually exist
        designable_residues = designable_residues.intersection(all_residues)
        
        if not designable_residues:
            print("Error: No valid designable residues found in the specified motif.", file=sys.stderr)
            sys.exit(1)

    # Calculate fixed residues by subtracting the designable set from the total
    fixed_residues = all_residues - designable_residues
    
    print(f"Designing {len(designable_residues)} residues: {sorted(list(designable_residues))}")
    print(f"Fixing {len(fixed_residues)} residues.")

    # 3. Create a temporary JSONL file for fixed positions
    fixed_positions_data: Dict[str, List[int]] = {args.ligand_chain: sorted(list(fixed_residues))}
    
    temp_jsonl_path = None
    try:
        # Create a named temporary file that subprocess can access
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.jsonl', encoding='utf-8') as f:
            json.dump(fixed_positions_data, f)
            f.write('\n')  # .jsonl requires a newline
            temp_jsonl_path = f.name
        
        print(f"Created temporary fixed positions file: {temp_jsonl_path}")

        # Get the absolute path to the ProteinMPNN script's *directory*
        # We will run the subprocess from this directory.
        protein_mpnn_dir = os.path.dirname(os.path.abspath(args.protein_mpnn_path))
        
        # Convert all user-provided paths to absolute paths
        # This ensures they are found correctly when the working directory is changed
        abs_pdb_path = os.path.abspath(args.pdb_file)
        abs_output_folder = os.path.abspath(args.output_folder)
        # temp_jsonl_path is already absolute from tempfile.NamedTemporaryFile

        # 4. Build and run the ProteinMPNN command
        # Combine ligand and target chains for the --pdb_path_chains argument
        all_chains = f"{args.ligand_chain} {args.target_chains}"
        
        command = [
            "python",
            args.protein_mpnn_path,
            "--pdb_path", abs_pdb_path,           # Use absolute path
            "--out_folder", abs_output_folder,    # Use absolute path
            "--fixed_positions_jsonl", temp_jsonl_path, # Already absolute
            "--pdb_path_chains", all_chains,
            "--num_seq_per_target", str(args.num_sequences),
            "--sampling_temp", str(args.temperature)
        ]

        print("\n" + "="*30)
        print("Running ProteinMPNN command:")
        # Join for a readable command
        print(" ".join(command))
        print(f"Running in directory: {protein_mpnn_dir}") # Add print for clarity
        print("="*30 + "\n")

        # Execute the command
        subprocess.run(command, check=True, cwd=protein_mpnn_dir)
        
        print("\nProteinMPNN run completed successfully.")
        print(f"Output sequences are in: {args.output_folder}")

    except subprocess.CalledProcessError as e:
        print(f"\nError running ProteinMPNN. Return code: {e.returncode}", file=sys.stderr)
        print(f"Command was: {' '.join(command)}", file=sys.stderr)
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}", file=sys.stderr)
    finally:
        # 5. Clean up the temporary file
        if temp_jsonl_path and os.path.exists(temp_jsonl_path):
            os.remove(temp_jsonl_path)
            print(f"Cleaned up temporary file: {temp_jsonl_path}")

def main():
    """Parses command-line arguments and starts the design process."""
    parser = argparse.ArgumentParser(
        description="Run ProteinMPNN to design a specific motif on a ligand in complex with a target.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    pdb_location = "Data/TIMP3_vs_ADAM17_X_ray.pdb"

    # --- Required Arguments ---
    req_group = parser.add_argument_group('Required Arguments')
    req_group.add_argument("--pdb_file", type=str, default=pdb_location,
                           help="Path to the PDB file containing the ligand-target complex.")
    req_group.add_argument("--ligand_chain", type=str, default="B",
                           help="Chain ID of the ligand you want to modify (e.g., 'A').")
    req_group.add_argument("--target_chains",  type=str, default="A",
                           help="Chain ID(s) of the fixed target protein, space-separated (e.g., 'B' or 'B C').")
    req_group.add_argument("--design_motif", type=int, default=(30, 36), nargs=2,
                           metavar=('START_RES', 'END_RES'),
                           help="Start and end residue numbers of the 6-amino-acid motif to design (inclusive).")
    req_group.add_argument("--protein_mpnn_path", type=str, default="../Tools/ProteinMPNN/protein_mpnn_run.py",
                           help="Full path to the 'protein_mpnn_run.py' script from the ProteinMPNN repo.")
    req_group.add_argument("--output_folder", type=str, default="Local/",
                           help="Path to the folder where output sequences will be saved.")
    
    # --- Optional Arguments ---
    opt_group = parser.add_argument_group('Optional Arguments')
    opt_group.add_argument("--num_sequences", type=int, default=10,
                           help="Number of sequences to generate (default: 10).")
    opt_group.add_argument("--temperature", type=float, default=0.1,
                           help="Sampling temperature for sequence generation (default: 0.1).")

    args = parser.parse_args()

    # --- Validation ---
    start, end = args.design_motif
    if start > end:
        print(f"Error: --design_motif START ({start}) must be less than or equal to END ({end}).", file=sys.stderr)
        sys.exit(1)
        
    motif_length = end - start + 1
    if motif_length != 6:
        print(f"Warning: The specified motif range {start}-{end} contains {motif_length} residues, not 6.")

    if not os.path.exists(args.output_folder):
        print(f"Creating output folder: {args.output_folder}")
        os.makedirs(args.output_folder, exist_ok=True)

    run_motif_design(args)

if __name__ == "__main__":
    main()
