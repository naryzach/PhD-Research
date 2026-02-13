import argparse, os, sys, glob, json, re
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

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
                    resname = residue.resname.upper()
                    seq.append(THREE_TO_ONE.get(resname, 'X'))
            if seq:
                sequences[chain.id] = "".join(seq)
    
    ALLOWED_IONS = {"MG", "ZN", "CL", "CA", "NA", "MN", "K", "FE", "CU", "CO"}
    
    for model in structure:
        for chain in model:
            if chain.id not in sequences or len(sequences[chain.id]) == 0:
                residues = list(chain.get_residues())
                if residues:
                    residue = residues[0]
                    resname = residue.resname.strip().upper()
                    if resname in ALLOWED_IONS:
                        sequences[chain.id] = resname
                    else:
                        found_ion = None
                        for atom in residue:
                            if atom.element.upper() in ALLOWED_IONS:
                                found_ion = atom.element.upper()
                                break
                        if found_ion:
                            sequences[chain.id] = found_ion
                        else:
                            if resname != 'HOH':
                                sequences[chain.id] = resname
        break
    return sequences

def parse_allmetal3d_metals(pdb_path):
    sites = []
    if not os.path.exists(pdb_path):
        return sites
        
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("HETATM"):
                try:
                    prob = float(line[60:66])
                    atom_name = line[12:16].strip()
                    element = line[76:78].strip()
                    if not element: element = atom_name
                    
                    sites.append({'ion': element, 'prob': prob})
                except:
                    pass
    return sites

def main():
    parser = argparse.ArgumentParser(description="Finalize Designs and Create Summary CSV")
    parser.add_argument("--pred_dir", default="Local/Metal_Predictions", help="Predictions Directory")
    parser.add_argument("--out_csv", default="Local/Metal_Predictions/design_summary.csv", help="Output CSV path")
    parser.add_argument("--out_json", default="Local/Metal_Predictions/alphafold_inputs.json", help="Output JSON path")
    args = parser.parse_args()
    
    data_all = []
    json_all = []
    
    print(f"Scanning {args.pred_dir} for ProteinMPNN outputs...")
    # Scan for ANY .fa files in a proteinmpnn/seqs directory
    seq_files = glob.glob(os.path.join(args.pred_dir, "**", "proteinmpnn", "seqs", "*.fa"), recursive=True)
    
    # Also try just searching for ANY .fa file and filtering by path
    if len(seq_files) == 0:
        print("Recursive glob failed, trying broader scan...")
        all_fa = glob.glob(os.path.join(args.pred_dir, "**", "*.fa"), recursive=True)
        seq_files = [f for f in all_fa if "proteinmpnn" in f and "seqs" in f]

    # Group by design
    design_map = {}
    for f in seq_files:
        filename = os.path.basename(f)
        if filename.endswith("_mpnn_input.fa"):
            basename = filename.replace("_mpnn_input.fa", "")
            type_ = "input"
        else:
            basename = filename.replace(".fa", "")
            type_ = "output"
            
        # Basename is NOT unique across directories (e.g. scaffold_0 exists in many folders).
        # We must key by the unique design identifier which includes the path structure.
        # But we want to group input and output for the SAME design.
        # They should be in the same directory.
        
        dir_path = os.path.dirname(f)
        unique_key = os.path.join(dir_path, basename)
        
        if unique_key not in design_map:
            design_map[unique_key] = {}
        design_map[unique_key][type_] = f
        
    final_files = []
    print(f"Found {len(design_map)} unique designs via MPNN files (raw count {len(seq_files)}).")
    
    for basename, files in design_map.items():
        if "output" in files:
            final_files.append(files["output"])
        elif "input" in files:
            final_files.append(files["input"])
            
    mpnn_files = final_files
    print(f"Processing {len(mpnn_files)} designs...")
    
    for fasta_file in mpnn_files:
         basename = os.path.basename(fasta_file).replace(".fa", "")
         if fasta_file.endswith("_mpnn_input.fa"):
             basename = os.path.basename(fasta_file).replace("_mpnn_input.fa", "")
             
         mpnn_dir = os.path.dirname(fasta_file) 
         proteinmpnn_conn_dir = os.path.dirname(mpnn_dir) 
         parent_dir = os.path.dirname(proteinmpnn_conn_dir) # .../Motif
         
         # Identify original PDB
         rfdiffusion_base = basename
         
         rfdiffusion_pdb = os.path.join(parent_dir, "rfdiffusion", f"{rfdiffusion_base}.pdb")
         if not os.path.exists(rfdiffusion_pdb):
             rfdiffusion_pdb = os.path.join(parent_dir, "rfdiffusion", "unidealized", f"{rfdiffusion_base}.pdb")
         
         if not os.path.exists(rfdiffusion_pdb):
             # Try stricter unidealized check or differing basename?
             # Sometimes MPNN modifies name? No, usually not.
             pass
         
         if not os.path.exists(rfdiffusion_pdb):
             # Final check: recursive search in motif dir for this pdb?
             # Might be slow but accurate
             found = glob.glob(os.path.join(parent_dir, "**", f"{rfdiffusion_base}.pdb"), recursive=True)
             if found:
                 rfdiffusion_pdb = found[0]
         
         if not os.path.exists(rfdiffusion_pdb):
             # print(f"Skipping {basename}: PDB not found")
             continue
             
         original_pdb = os.path.abspath(rfdiffusion_pdb)
         design_id = basename
         
         # Determine Metadata from Parent Dir structure
         # Path: .../Type/Ion/Motif/
         
         rel_path = os.path.relpath(original_pdb, args.pred_dir)
         parts = rel_path.split(os.sep)
         
         category = "Unknown"
         ion_type = "Unknown"
         motif_name = "Unknown"
         
         # Robust metadata extraction relative to pred_dir or known structure
         # We know parent_dir is the Motif dir
         try:
             motif_name = os.path.basename(parent_dir)
             ion_dir_path = os.path.dirname(parent_dir)
             ion_type = os.path.basename(ion_dir_path)
             cat_dir_path = os.path.dirname(ion_dir_path)
             category_dir = os.path.basename(cat_dir_path)
             
             if "DeNovo" in category_dir: category = "Scaffold"
             elif "Grafted" in category_dir: category = "Graft"
             elif "Swaps" in category_dir: category = "Swap"
         except:
             pass

         # 1. Parse ProteinMPNN Sequence and Score
         mpnn_seq = None
         mpnn_score = None
         mpnn_design_idx = 0 # Default if input file
         
         try:
             with open(fasta_file, 'r') as f:
                 lines = f.readlines()
             best_s = float('inf')
             
             # Scan headers to find best score
             for i, line in enumerate(lines):
                 if line.startswith(">"):
                     header = line
                     # Example header: >design_0_sample=1 score=1.234
                     # Global seq index? 
                     # Usually MPNN output has multiple sequences.
                     # We want the index of the BEST one.
                     
                     # Check if it's an MPNN output header
                     score_match = re.search(r"score=([-\d\.]+)", header)
                     if score_match:
                         score = float(score_match.group(1))
                         if score < best_s:
                             best_s = score
                             # Get sequence from next line
                             if i + 1 < len(lines):
                                 mpnn_seq = lines[i+1].strip()
                                 mpnn_score = score
                                 
                                 # Try to extract sample index
                                 # >...sample=5... or just use loop index? 
                                 # Standard MPNN: sample=5
                                 sample_match = re.search(r"sample=(\d+)", header)
                                 if sample_match:
                                     mpnn_design_idx = int(sample_match.group(1))
                                 else:
                                     # Fallback: create an index based on order
                                     # Not ideal but works.
                                     pass
         except:
             pass
         
         full_sequence = mpnn_seq
         if not full_sequence:
             # Extract from PDB
             seqs = get_sequence(original_pdb)
             if seqs:
                 full_sequence = "/".join(seqs.values())
        
         # 2. Check AllMetal3D Results
         allmetal_dir = os.path.join(parent_dir, "allmetal3d", f"{basename}_am3d")   
         am3d_metals_pdb = os.path.join(allmetal_dir, f"{basename}_threaded_metals.pdb")
         
         am3d_sites = parse_allmetal3d_metals(am3d_metals_pdb)
         
         am3d_best_prob = 0.0
         am3d_best_ion = "None"
         am3d_count = len(am3d_sites)
         
         if am3d_sites:
             am3d_sites.sort(key=lambda x: x['prob'], reverse=True)
             am3d_best_prob = am3d_sites[0]['prob']
             am3d_best_ion = am3d_sites[0]['ion']
        
         entry = {
             "id": design_id,
             "category": category,
             "target_ion": ion_type,
             "motif": motif_name,
             "mpnn_score": mpnn_score,
             "am3d_count": am3d_count,
             "am3d_best_prob": am3d_best_prob,
             "am3d_best_ion": am3d_best_ion,
             "sequence": full_sequence,
             "mpnn_design_idx": mpnn_design_idx
         }
         data_all.append(entry)

         # JSON Logic (Simplified for brevity but retaining structure)
         seqs_dict = get_sequence(original_pdb)
         if not seqs_dict: seqs_dict = {}
         all_chains = sorted(seqs_dict.keys())
         protein_chains = [c for c in all_chains if len(seqs_dict[c]) > 5]
         ligand_chains = [c for c in all_chains if len(seqs_dict[c]) <= 5]
         
         json_sequences = []
         if mpnn_seq:
            mpnn_parts = mpnn_seq.split('/')
            if len(mpnn_parts) == len(protein_chains):
                 for i, chain_id in enumerate(protein_chains):
                     json_sequences.append({"proteinChain": {"sequence": mpnn_parts[i], "count": 1}})
            else:
                 for chain_id in protein_chains:
                     json_sequences.append({"proteinChain": {"sequence": seqs_dict[chain_id], "count": 1}})
         else:
             for chain_id in protein_chains:
                 json_sequences.append({"proteinChain": {"sequence": seqs_dict[chain_id], "count": 1}})

         for chain_id in ligand_chains:
             json_sequences.append({"ion": {"ion": seqs_dict[chain_id], "count": 1}})
            
         entry_json = {
            "name": design_id,
            "sequences": json_sequences,
            "category": category # Temporary for splitting
         }
         json_all.append(entry_json)

    # Save CSV
    if data_all:
        df = pd.DataFrame(data_all)
        # Select columns as requested
        cols = ["id","category","target_ion","motif","mpnn_score","mpnn_design_idx","am3d_count","am3d_best_prob","am3d_best_ion","sequence"]
        df = df[cols]
        
        df.to_csv(args.out_csv, index=False)
        print(f"Summary saved to {args.out_csv}")
        
        # Split Categories
        for cat in df['category'].unique():
             sub_df = df[df['category'] == cat]
             base_fn = os.path.basename(args.out_csv).replace(".csv", "")
             out_fn = f"{base_fn}_{cat}.csv" if cat != "Unknown" else f"{base_fn}_Unknown.csv"
             # Map requested names
             if cat == "Scaffold": out_fn = f"{base_fn}_DeNovo.csv"
             elif cat == "Graft": out_fn = f"{base_fn}_Grafted.csv"
             elif cat == "Swap": out_fn = f"{base_fn}_Swaps.csv"
             
             out_path = os.path.join(os.path.dirname(args.out_csv), out_fn)
             sub_df.to_csv(out_path, index=False)
             print(f"  Saved subset: {out_path}")
             
    else:
        print("No designs found.")

    # Save JSON
    if json_all:
        # Clean up category from main list before saving, or keep it?
        # Standard: remove internal helper key
        clean_json = [{k:v for k,v in x.items() if k!='category'} for x in json_all]
        
        with open(args.out_json, 'w') as f:
            json.dump(clean_json, f, indent=2)
        print(f"JSON saved to {args.out_json}")
        
        # Split JSONs
        # Group by category
        from collections import defaultdict
        cats = defaultdict(list)
        for entry in json_all:
            cat = entry.get('category', 'Unknown')
            # Clean entry
            clean_entry = {k:v for k,v in entry.items() if k!='category'}
            cats[cat].append(clean_entry)
            
        for cat, entries in cats.items():
             base_fn = os.path.basename(args.out_json).replace(".json", "")
             out_fn = f"{base_fn}_{cat}.json"
             # Map requested names
             if cat == "Scaffold": out_fn = f"{base_fn}_DeNovo.json"
             elif cat == "Graft": out_fn = f"{base_fn}_Grafted.json"
             elif cat == "Swap": out_fn = f"{base_fn}_Swaps.json"
             
             out_path = os.path.join(os.path.dirname(args.out_json), out_fn)
             with open(out_path, 'w') as f:
                 json.dump(entries, f, indent=2)
             print(f"  Saved JSON subset: {out_path}")

if __name__ == "__main__":
    main()
