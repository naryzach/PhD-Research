import sys
import subprocess
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

def run_scaffold(args):
    motif_path = os.path.abspath(args.motif_pdb)
    out_dir = os.path.abspath(args.out)
    rf_path = os.path.abspath(args.rf_path)
    lmpnn_path = os.path.abspath(args.lmpnn_path)
    
    os.makedirs(out_dir, exist_ok=True)
    
    # --- RFdiffusion ---
    rf_out = os.path.join(out_dir, "rfdiffusion")
    os.makedirs(rf_out, exist_ok=True)
    
    # Scaffolding strategy
    # Parse the input PDB to find the actual residues to fix
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("motif", motif_path)
    
    residues = []
    hetatms_lines = [] # To store metal/ligand lines to append later
    
    # Read raw lines for HETATMs to append later (simple text processing)
    with open(motif_path, 'r') as f:
        for line in f:
            if line.startswith("HETATM") or line.startswith("CONECT"):
                 hetatms_lines.append(line)
    
    for model in structure:
        for chain in model:
            for residue in chain:
                # Check if standard amino acid
                if is_aa(residue, standard=True):
                    residues.append((chain, residue.id[1]))
                # Else: it's likely the metal or ligand, we ignore it for the CONTIG but keep it in the PDB context

    
    if not residues:
        print(f"Error: No residues found in {motif_path}")
        return

    # Group into contiguous segments
    # Sort by Chain, then Residue ID
    residues.sort(key=lambda x: (x[0], x[1]))
    
    segments = []
    if residues:
        current_chain = residues[0][0]
        start_res = residues[0][1]
        prev_res = residues[0][1]
        
        for i in range(1, len(residues)):
            chain, res_id = residues[i]
            if chain == current_chain and res_id == prev_res + 1:
                # Contiguous
                prev_res = res_id
            else:
                # Break
                segments.append(f"{current_chain}{start_res}-{prev_res}")
                current_chain = chain
                start_res = res_id
                prev_res = res_id
        # Add last segment
        segments.append(f"{current_chain}{start_res}-{prev_res}")
    
    # Construct Contig
    # Strategy: [50 / Seg1 / 10-20 / Seg2 / 50]
    # We will put flexible loopers between segments.
    
    if args.contig:
        contig = args.contig
    else:
        # Build dynamic contig
        contig_parts = ["50"] # Start with N-term scaffold
        
        for seg in segments:
            contig_parts.append(seg)
            contig_parts.append("10-20") # Linker between segments
            
        contig_parts[-1] = "50" # Replace last linker with C-term scaffold
        
        # Format: ['Part/Part/Part']
        contig_str = "/".join(contig_parts)
        contig = f"['{contig_str}']"
        print(f"Auto-generated contig: {contig}")

    print(f"Running RFdiffusion with contig: {contig}")
    
    cmd_rf = [
        sys.executable, os.path.join(rf_path, 'scripts/run_inference.py'),
        f"inference.input_pdb={motif_path}",
        f"inference.output_prefix={os.path.join(rf_out, 'scaffold')}",
        f"inference.num_designs={args.num_designs}",
        f"contigmap.contigs={contig}"
    ]
    
    if args.dry_run:
        print("Dry Run RFdiffusion:", " ".join(cmd_rf))
    else:
        try:
            subprocess.run(cmd_rf, check=True)
        except subprocess.CalledProcessError as e:
            print(f"RFdiffusion failed: {e}")
            return # Don't proceed to LigandMPNN if scaffold failed

    # --- Post-Processing: Merge Metals back ---
    # RFdiffusion output PDBs likely don't have the metal if it wasn't in the contig as 'fixed'.
    # We need to append the HETATM lines from the original motif file to the generated scaffolds.
    scaffold_pdbs = [os.path.join(rf_out, f) for f in os.listdir(rf_out) if f.endswith(".pdb")]
    print(f"Merging metals into {len(scaffold_pdbs)} scaffolds...")
    
    for pdb in scaffold_pdbs:
        with open(pdb, 'a') as f:
            f.write("\n")
            for line in hetatms_lines:
                f.write(line)

    # --- LigandMPNN ---
    lmpnn_out = os.path.join(out_dir, "ligandmpnn")
    os.makedirs(lmpnn_out, exist_ok=True)
    
    print("Running LigandMPNN...")
    
    cmd_lmpnn = [
        sys.executable, os.path.join(lmpnn_path, 'run.py'),
        "--pdb_path", rf_out,
        "--out_folder", lmpnn_out,
        "--model_type", "ligand_mpnn",
        "--num_seq_per_target", "2",
        "--batch_size", "1"
    ]
    
    if args.dry_run:
        print("Dry Run LigandMPNN:", " ".join(cmd_lmpnn))
    else:
        try:
            subprocess.run(cmd_lmpnn, check=True)
        except subprocess.CalledProcessError as e:
            print(f"LigandMPNN failed: {e}")

def main():
    parser = argparse.ArgumentParser(description="Automated Motif Scaffolding")
    parser.add_argument("--motif_pdb", required=True, help="Input Motif PDB (Metal + Residues)")
    parser.add_argument("--motif_chain", default="A", help="Chain ID of the motif in the input PDB")
    parser.add_argument("--contig", help="RFdiffusion contig string (e.g. '[50/A1-10/50]')")
    parser.add_argument("--out", default="./motif_scaffold_output", help="Output directory")
    parser.add_argument("--num_designs", type=int, default=10)
    parser.add_argument("--dry_run", action="store_true")
    
    parser.add_argument("--rf_path", default="../Tools/RFdiffusion", help="Path to RFdiffusion")
    parser.add_argument("--lmpnn_path", default="../Tools/LigandMPNN", help="Path to LigandMPNN")
    
    args = parser.parse_args()
    run_scaffold(args)

if __name__ == "__main__":
    main()
