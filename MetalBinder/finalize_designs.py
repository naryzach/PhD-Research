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

def get_allmetal3d_metrics(parent_dir, basename, source):
    # Defaults
    metrics = {
        "am3d_count": 0,
        "am3d_best_prob": 0.0,
        "am3d_best_ion": "None"
    }
    
    # Paths
    json_path = None
    pdb_path = None
    
    if source == 'ligandmpnn':
        # .../Motif/allmetal3d_ligand/basename_results.json
        json_path = os.path.join(parent_dir, "allmetal3d_ligand", f"{basename}_results.json")
    else:
        # .../Motif/allmetal3d/basename_am3d/basename_results.json
        # OR .../basename_threaded_metals.pdb
        am3d_sub = os.path.join(parent_dir, "allmetal3d", f"{basename}_am3d")
        json_path = os.path.join(am3d_sub, f"{basename}_results.json")
        pdb_path = os.path.join(am3d_sub, f"{basename}_threaded_metals.pdb")

    # Try JSON first (Cleanest)
    if json_path and os.path.exists(json_path):
        try:
            with open(json_path, 'r') as f:
                data = json.load(f)
            
            # wrapper detection (run_allmetal3d usually wraps in "results")
            results = data.get("results", data) 
            
            # Flatten to find max
            # results = {path: {res: {ion: prob}}}
            all_sites = []
            
            # Handle potential different structures
            if isinstance(results, dict):
                for path_k, res_v in results.items():
                    if isinstance(res_v, dict):
                        for res_k, ion_v in res_v.items():
                            if isinstance(ion_v, dict):
                                for ion, prob in ion_v.items():
                                    all_sites.append({'ion': ion, 'prob': prob})
                            
            if all_sites:
                all_sites.sort(key=lambda x: x['prob'], reverse=True)
                metrics["am3d_count"] = len(all_sites)
                metrics["am3d_best_prob"] = all_sites[0]['prob']
                metrics["am3d_best_ion"] = all_sites[0]['ion']
                
            return metrics
        except:
            pass # Fallback

    # Try PDB (Legacy / ProteinMPNN threaded)
    if pdb_path and os.path.exists(pdb_path):
         sites = parse_allmetal3d_metals(pdb_path)
         if sites:
             sites.sort(key=lambda x: x['prob'], reverse=True)
             metrics["am3d_count"] = len(sites)
             metrics["am3d_best_prob"] = sites[0]['prob']
             metrics["am3d_best_ion"] = sites[0]['ion']
             
    return metrics

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
    
    # Scan for LigandMPNN FASTAs (in seqs folder)
    # Pattern: .../ligandmpnn/seqs/*.fa
    lmpnn_seq_files = glob.glob(os.path.join(args.pred_dir, "**", "ligandmpnn", "seqs", "*.fa"), recursive=True)
    
    print(f"Found {len(seq_files)} ProteinMPNN FASTAs.")
    print(f"Found {len(lmpnn_seq_files)} LigandMPNN FASTAs.")

    # Group by design
    design_map = {}
    
    # Process ProteinMPNN Files
    for f in seq_files:
        filename = os.path.basename(f)
        if filename.endswith("_mpnn_input.fa"):
            basename = filename.replace("_mpnn_input.fa", "")
            type_ = "input"
        else:
            basename = filename.replace(".fa", "")
            type_ = "output"
            
        dir_path = os.path.dirname(f)
        unique_key = os.path.join(dir_path, basename)
        
        if unique_key not in design_map:
            design_map[unique_key] = {'source': 'proteinmpnn'}
        design_map[unique_key][type_] = f

    # Process LigandMPNN Files
    for f in lmpnn_seq_files:
        # For LigandMPNN, one FASTA contains multiple designs (usually 1-N)
        # We treat the FASTA as the "output" file source
        filename = os.path.basename(f)
        basename = filename.replace(".fa", "")
        dir_path = os.path.dirname(f)
        unique_key = os.path.join(dir_path, basename)
        
        if unique_key not in design_map:
            design_map[unique_key] = {'source': 'ligandmpnn'}
        design_map[unique_key]['output'] = f
        
    print(f"Processing unique design sets...")
    
    for unique_key, data in design_map.items():
         source = data.get('source', 'proteinmpnn')
         
         if source == 'proteinmpnn':
             if 'output' in data: design_file = data['output']
             elif 'input' in data: design_file = data['input']
             else: continue
             
             basename = os.path.basename(design_file).replace(".fa", "").replace("_mpnn_input", "")
             mpnn_dir = os.path.dirname(design_file)
             conn_dir = os.path.dirname(mpnn_dir) 
             parent_dir = os.path.dirname(conn_dir) # .../Motif
             
             # Identify original PDB (Common)
             rfdiffusion_base = basename
             rfdiffusion_pdb = os.path.join(parent_dir, "rfdiffusion", f"{rfdiffusion_base}.pdb")
             if not os.path.exists(rfdiffusion_pdb):
                 rfdiffusion_pdb = os.path.join(parent_dir, "rfdiffusion", "unidealized", f"{rfdiffusion_base}.pdb")
             if not os.path.exists(rfdiffusion_pdb):
                 found = glob.glob(os.path.join(parent_dir, "**", f"{rfdiffusion_base}.pdb"), recursive=True)
                 if found: rfdiffusion_pdb = found[0]
                 
             if not os.path.exists(rfdiffusion_pdb): continue
             original_pdb = os.path.abspath(rfdiffusion_pdb)
             
             # Metadata extraction
             motif_name, ion_type, category = "Unknown", "Unknown", "Unknown"
             try:
                 motif_name = os.path.basename(parent_dir)
                 ion_dir_path = os.path.dirname(parent_dir)
                 ion_type = os.path.basename(ion_dir_path)
                 cat_dir_path = os.path.dirname(ion_dir_path)
                 category_dir = os.path.basename(cat_dir_path)
                 if "DeNovo" in category_dir: category = "Scaffold"
                 elif "Grafted" in category_dir: category = "Graft"
                 elif "Swaps" in category_dir: category = "Swap"
             except: pass
             
             # Single Best Parse
             mpnn_seq = None
             mpnn_score = None
             mpnn_design_idx = 0
             try:
                 with open(design_file, 'r') as f:
                     lines = f.readlines()
                 best_s = float('inf')
                 for i, line in enumerate(lines):
                     if line.startswith(">"):
                         header = line
                         score_match = re.search(r"score=([-\d\.]+)", header)
                         if score_match:
                             score = float(score_match.group(1))
                             if score < best_s:
                                 best_s = score
                                 if i + 1 < len(lines):
                                     mpnn_seq = lines[i+1].strip()
                                     mpnn_score = score
                                     sample_match = re.search(r"sample=(\d+)", header)
                                     if sample_match: mpnn_design_idx = int(sample_match.group(1))
             except: pass
             
             full_sequence = mpnn_seq
             if not full_sequence:
                 seqs = get_sequence(original_pdb)
                 if seqs: full_sequence = "/".join(seqs.values())
                 
             am3d_metrics = get_allmetal3d_metrics(parent_dir, basename, source)
             
             entry = {
                 "id": basename,
                 "source": source,
                 "category": category,
                 "target_ion": ion_type,
                 "motif": motif_name,
                 "mpnn_score": mpnn_score,
                 "am3d_count": am3d_metrics["am3d_count"],
                 "am3d_best_prob": am3d_metrics["am3d_best_prob"],
                 "am3d_best_ion": am3d_metrics["am3d_best_ion"],
                 "sequence": full_sequence,
                 "mpnn_design_idx": mpnn_design_idx
             }
             data_all.append(entry)
             
             # JSON Logic for ProteinMPNN
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
                "name": basename,
                "sequences": json_sequences,
                "category": category
             }
             json_all.append(entry_json)
             
         elif source == 'ligandmpnn':
             design_file = data['output'] # FASTA
             basename = os.path.basename(design_file).replace(".fa", "")
             mpnn_dir = os.path.dirname(design_file) # .../seqs
             conn_dir = os.path.dirname(mpnn_dir) # .../ligandmpnn
             parent_dir = os.path.dirname(conn_dir) # .../Motif
             
             # Identify original PDB (same logic)
             rfdiffusion_base = basename
             rfdiffusion_pdb = os.path.join(parent_dir, "rfdiffusion", f"{rfdiffusion_base}.pdb")
             if not os.path.exists(rfdiffusion_pdb):
                 rfdiffusion_pdb = os.path.join(parent_dir, "rfdiffusion", "unidealized", f"{rfdiffusion_base}.pdb")
             if not os.path.exists(rfdiffusion_pdb):
                 found = glob.glob(os.path.join(parent_dir, "**", f"{rfdiffusion_base}.pdb"), recursive=True)
                 if found: rfdiffusion_pdb = found[0]
             
             # Metadata
             motif_name, ion_type, category = "Unknown", "Unknown", "Unknown"
             try:
                 motif_name = os.path.basename(parent_dir)
                 ion_dir_path = os.path.dirname(parent_dir)
                 ion_type = os.path.basename(ion_dir_path)
                 cat_dir_path = os.path.dirname(ion_dir_path)
                 category_dir = os.path.basename(cat_dir_path)
                 if "DeNovo" in category_dir: category = "Scaffold"
                 elif "Grafted" in category_dir: category = "Graft"
                 elif "Swaps" in category_dir: category = "Swap"
             except: pass

             # Parse ALL sequences in FASTA
             try:
                 with open(design_file, 'r') as f:
                     content = f.read()
                 
                 # Split by '>'
                 entries = content.split('>')
                 for chunk in entries:
                     if not chunk.strip(): continue
                     header_line = chunk.strip().split('\n')[0]
                     
                     id_match = re.search(r"id=(\d+)", header_line)
                     conf_match = re.search(r"overall_confidence=([-\d\.]+)", header_line)
                     
                     if id_match:
                         idx = int(id_match.group(1))
                         score = float(conf_match.group(1)) if conf_match else 0.0
                         
                         # Check for corresponding PDB/AllMetal3D result
                         sub_basename = f"{basename}_packed_{idx}_1"
                         packed_pdb_path = os.path.join(parent_dir, "ligandmpnn", "packed", f"{sub_basename}.pdb")
                         
                         am3d_metrics = get_allmetal3d_metrics(parent_dir, sub_basename, source)
                         
                         entry_id = f"{basename}_{idx}"
                         
                         # Get sequence from Packed PDB (Robusted for JSON)
                         full_sequence = ""
                         seqs_dict = {}
                         if os.path.exists(packed_pdb_path):
                             seqs_dict = get_sequence(packed_pdb_path)
                             if seqs_dict: full_sequence = "/".join(seqs_dict.values())
                         else:
                             # Fallback to FASTA seq (might be incomplete chain info)
                             full_sequence = "".join(chunk.strip().split('\n')[1:])

                         entry = {
                             "id": entry_id,
                             "source": source,
                             "category": category,
                             "target_ion": ion_type,
                             "motif": motif_name,
                             "mpnn_score": score, 
                             "am3d_count": am3d_metrics["am3d_count"],
                             "am3d_best_prob": am3d_metrics["am3d_best_prob"],
                             "am3d_best_ion": am3d_metrics["am3d_best_ion"],
                             "sequence": full_sequence,
                             "mpnn_design_idx": idx
                         }
                         data_all.append(entry)
                         
                         # JSON Logic for LigandMPNN
                         if not seqs_dict: seqs_dict = {}
                         all_chains = sorted(seqs_dict.keys())
                         # Heuristic: shorter count usually ligand?
                         # Or check standard AA.
                         protein_chains = []
                         ligand_chains = []
                         
                         # Use helper logic from main
                         from Bio.PDB.Polypeptide import is_aa 
                         # Actually get_sequence filters non-standard for seqs.
                         # If length > 5 likely protein.
                         protein_chains = [c for c in all_chains if len(seqs_dict[c]) > 5]
                         ligand_chains = [c for c in all_chains if len(seqs_dict[c]) <= 5]

                         json_sequences = []
                         for chain_id in protein_chains:
                             json_sequences.append({"proteinChain": {"sequence": seqs_dict[chain_id], "count": 1}})
                         for chain_id in ligand_chains:
                             json_sequences.append({"ion": {"ion": seqs_dict[chain_id], "count": 1}})
                         
                         entry_json = {
                            "name": entry_id,
                            "sequences": json_sequences,
                            "category": category
                         }
                         json_all.append(entry_json)
                     
             except Exception as e:
                 print(f"Error parsing LigandMPNN FASTA {design_file}: {e}")

    # Save CSV
    if data_all:
        df = pd.DataFrame(data_all)
        # Select columns as requested
        cols = ["id","source","category","target_ion","motif","mpnn_score","mpnn_design_idx","am3d_count","am3d_best_prob","am3d_best_ion","sequence"]
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
