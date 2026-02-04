import argparse
import os
import glob
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_one, is_aa

def get_sequence(pdb_path, chain_id=None):
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("temp", pdb_path)
    except:
        return None
        
    sequences = {}
    
    for model in structure:
        for chain in model:
            if chain_id and chain.id != chain_id:
                continue
            
            seq = []
            for residue in chain:
                if is_aa(residue, standard=True):
                    try:
                        seq.append(three_to_one(residue.resname))
                    except:
                        seq.append("X")
            if seq:
                sequences[chain.id] = "".join(seq)
        break # First model
        
    return sequences

def main():
    parser = argparse.ArgumentParser(description="Finalize Designs and Create AlphaFold CSV")
    parser.add_argument("--pred_dir", default="../Local/Metal_Predictions", help="Predictions Directory")
    parser.add_argument("--run_mpnn", action="store_true", help="Run ProteinMPNN on designs")
    parser.add_argument("--lmpnn_path", default="../Tools/LigandMPNN", help="Path to LigandMPNN")
    args = parser.parse_args()
    
    # recursive search
    print(f"Scanning {args.pred_dir}...")
    pdb_files = glob.glob(os.path.join(args.pred_dir, "**", "*.pdb"), recursive=True)
    
    # Filter for RFdiffusion outputs
    design_files = [f for f in pdb_files if "rfdiffusion" in f and ("scaffold" in os.path.basename(f) or "design" in os.path.basename(f))]
    
    if args.run_mpnn:
        print(f"Running ProteinMPNN on {len(design_files)} designs...")
        project_root = os.path.dirname(os.path.abspath(__file__))
        
        for i, pdb in enumerate(design_files):
            # Output folder: sibling of rfdiffusion
            # .../DeNovo/ZN/Motif/rfdiffusion/scaffold.pdb
            # -> .../DeNovo/ZN/Motif/ligandmpnn/
            
            dirname = os.path.dirname(pdb)
            parent_dir = os.path.dirname(dirname) 
            mpnn_out = os.path.join(parent_dir, "ligandmpnn")
            os.makedirs(mpnn_out, exist_ok=True)
            
            # Check if already run
            # LigandMPNN output naming: {input_name}_mpnn_...
            # We can just run it.
            
            cmd = [
                sys.executable, os.path.join(args.lmpnn_path, 'run.py'),
                "--pdb_path", dirname, # RFdiffusion output folder
                "--out_folder", mpnn_out,
                "--model_type", "ligand_mpnn",
                "--num_seq_per_target", "2",
                "--batch_size", "1"
            ]
            
            print(f"[{i+1}/{len(design_files)}] Processing {os.path.basename(pdb)}")
            try:
                subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except subprocess.CalledProcessError:
                print(f"  MPNN failed for {pdb} (Check dependencies: pip install dm-tree prody ml_collections)")
                continue

    # Now collect sequences
    # Prefer MPNN outputs
    data = []
    
    # Re-scan for MPNN outputs if run_mpnn was used or just check existence
    # Logic: Iterate RFdiffusion designs, look for corresponding MPNN designs.
    
    for pdb in design_files:
        dirname = os.path.dirname(pdb)
        parent_dir = os.path.dirname(dirname)
        mpnn_out = os.path.join(parent_dir, "ligandmpnn")
        
        # Original Name: scaffold_0.pdb
        basename = os.path.basename(pdb).replace(".pdb", "")
        
        # MPNN outputs: scaffold_0.pdb (it copies it?) or scaffold_0_seqs...
        # LigandMPNN outputs .pdb files with sequences.
        # Usually {original_name}_1.pdb, {original_name}_2.pdb etc.
        
        found_seqs = False
        
        if os.path.exists(mpnn_out):
            mpnn_pdbs = glob.glob(os.path.join(mpnn_out, f"{basename}*.pdb"))
            for mp_pdb in mpnn_pdbs:
                seqs = get_sequence(mp_pdb)
                if seqs:
                    # Filter for Poly-Glycine (if MPNN failed or mapped back incorrectly?)
                    # Assuming MPNN always produces valid seqs.
                    full_seq = ":".join(seqs.values())
                    
                    # Create ID
                    rel_path = os.path.relpath(mp_pdb, args.pred_dir)
                    unique_id = rel_path.replace("/", "_").replace(".pdb", "")
                    
                    data.append({
                        "id": unique_id,
                        "sequence": full_seq,
                        "path": os.path.abspath(mp_pdb),
                        "source": "ProteinMPNN"
                    })
                    found_seqs = True
                    
        if not found_seqs:
            # Fallback to RFdiffusion backbone (Poly-Gly?)
            seqs = get_sequence(pdb)
            if seqs:
                full_seq = ":".join(seqs.values())
                rel_path = os.path.relpath(pdb, args.pred_dir)
                unique_id = rel_path.replace("/", "_").replace(".pdb", "")
                
                data.append({
                    "id": unique_id,
                    "sequence": full_seq, 
                    "path": os.path.abspath(pdb),
                    "source": "RFdiffusion_Backbone"
                })

    if data:
        df = pd.DataFrame(data)
        df.to_csv(args.out_csv, index=False)
        print(f"Saved {len(df)} sequences to {args.out_csv}")
    else:
        print("No designs found.")

if __name__ == "__main__":
    main()
