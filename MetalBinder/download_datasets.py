import argparse
import os
import requests
import json
import shutil
import time
from Bio.PDB import MMCIFParser, PDBParser

# --- CONFIGURATION FROM METAL_MINER ---

# 1. GENERAL MINING: Broad search for these ions
# Note: Expanded from original download_datasets.py
DEFAULT_METALS = ["ZN", "CU", "NI", "CO", "MN", "FE", "MG", "CA"]

TARGET_METALS_CATEGORIES = {
    "Common": ["FE", "CU", "ZN", "MN", "NI", "CO", "MG", "CA"],
    "Rare_Earths": ["LA", "CE", "ND", "EU", "GD", "TB", "DY", "YB", "LU", "Y"], 
    "Actinides": ["U", "URI", "PU", "TH"], 
    "Toxic": ["HG", "CD", "PB", "AS"]
}

# 2. HAND-PICKED: The "Greatest Hits" of metal binding
HAND_PICKED_PDBS = {
    "Siderophores_Iron": ["1EFD", "6Z2N", "1FHA", "2IAH"],
    "Copper_Machinery": ["4AZU", "1PLC", "7ZC3", "1CC8"],
    "Metallothioneins_Sponges": ["4MT2", "1MHU", "2FJ5", "1DME"],
    "Lanthanide_Specialists": ["6MI5", "8DQ2", "9C8Z", "7CCO"],
    "Actinide_Engineered": ["4GLW", "4FZP"]
}

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"

def search_rcsb(query_payload):
    """Helper to send queries to RCSB."""
    try:
        response = requests.post(SEARCH_URL, json=query_payload)
        if response.status_code == 200:
            result = response.json()
            if 'result_set' in result and len(result['result_set']) > 0:
                return [entry['identifier'] for entry in result['result_set']]
        return []
    except Exception as e:
        print(f"Search Error: {e}")
        return []

def get_pdb_ids_for_ligand(ligand_id, limit=10):
    """Standard metal search (General Mining)."""
    # Uses the improved logic from metal_miner.py (resolution restriction)
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                        "operator": "exact_match",
                        "value": ligand_id
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal",
                        "value": 2.5
                    }
                }
            ]
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": limit}}
    }
    return search_rcsb(query)

def get_metalloproteinases(limit=50):
    """Targets Matrix Metalloproteinases (MMPs) and ADAMs."""
    print(f"\n--- HUNTING METALLOPROTEINASES (MMPs & ADAMs) ---")
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "group",
                    "logical_operator": "or",
                    "nodes": [
                        {"type": "terminal", "service": "full_text", "parameters": {"value": "Matrix Metalloproteinase"}},
                        {"type": "terminal", "service": "full_text", "parameters": {"value": "ADAM17"}},
                        {"type": "terminal", "service": "full_text", "parameters": {"value": "ADAM10"}},
                        {"type": "terminal", "service": "full_text", "parameters": {"value": "MMP-9"}}
                    ]
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                        "operator": "exact_match",
                        "value": "ZN" # Ensure they are the active metallo-forms
                    }
                }
            ]
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": limit}}
    }
    return search_rcsb(query)

def download_file(pdb_id, out_dir, file_format="cif"):
    """Downloads a single PDB/CIF file."""
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    pdb_id = pdb_id.upper()
    ext = "cif" if file_format == "cif" else "pdb"
    filename = f"{pdb_id}.{ext}"
    final_path = os.path.join(out_dir, filename)
    
    if os.path.exists(final_path) and os.path.getsize(final_path) > 0:
        return final_path # Exists
        
    url = f"https://files.rcsb.org/download/{filename}"
    try:
        r = requests.get(url)
        if r.status_code == 200:
            with open(final_path, 'wb') as f:
                f.write(r.content)
            # print(f"  Downloaded {filename}")
            return final_path
        else:
            print(f"  Failed to download {filename} (Status: {r.status_code})")
            return None
    except Exception as e:
        print(f"  Error downloading {filename}: {e}")
        return None

def sort_into_ions(file_path, base_out_dir, valid_metals=None):
    """
    Parses the valid PDB/CIF file to find which metals it contains,
    then symlinks/copies it to the appropriate Ion folders.
    """
    if not file_path or not os.path.exists(file_path):
        return

    # Use defaults if not provided
    if valid_metals is None:
        valid_metals = DEFAULT_METALS

    basename = os.path.basename(file_path)
    
    # Parse file to find metals
    found_metals = set()
    
    # Fast parsing attempt (text based) to avoid loading full structure if possible?
    # Actually, MMCIF/PDB parsing is safer.
    
    try:
        if file_path.endswith(".cif"):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
            
        structure = parser.get_structure("temp", file_path)
        
        for atom in structure.get_atoms():
            element = atom.element.upper()
            if element in valid_metals:
                found_metals.add(element)
                
    except Exception as e:
        print(f"  Warning: Could not parse {basename} for sorting: {e}")
        return

    if not found_metals:
        print(f"  Warning: No valid metals found in {basename} (checked: {valid_metals})")
        return

    # Sort
    for metal in found_metals:
        metal_dir = os.path.join(base_out_dir, metal)
        if not os.path.exists(metal_dir):
            os.makedirs(metal_dir)
            
        dest_path = os.path.join(metal_dir, basename)
        
        # Copy if not exists
        if not os.path.exists(dest_path):
            shutil.copy2(file_path, dest_path)
            # print(f"  -> Sorted into {metal}")

def main():
    parser = argparse.ArgumentParser(description="Download Metal Binding PDBs/CIFs (v2.0 Integrated)")
    parser.add_argument("--mode", choices=["general", "curated", "mmp"], default="general", 
                        help="Mining mode: 'general' (search by ion), 'curated' (hand-picked), 'mmp' (metalloproteinases)")
    parser.add_argument("--metals", default=",".join(DEFAULT_METALS), help="Comma-separated list of metals (for General/Sorting)")
    parser.add_argument("--out", default="../Local/Metal_PDBs", help="Base output directory")
    parser.add_argument("--limit", type=int, default=10, help="Max number of files per metal (General mode only)")
    parser.add_argument("--format", choices=["cif", "pdb"], default="cif", help="File format to download")
    
    args = parser.parse_args()
    
    # Flatten all categories into a single list for sorting safety
    MINING_METALS = [m.strip().upper() for m in args.metals.split(',')]
    
    # Create a superset of all known metals for sorting
    ALL_KNOWN_METALS = set(DEFAULT_METALS)
    for cat_metals in TARGET_METALS_CATEGORIES.values():
        ALL_KNOWN_METALS.update(cat_metals)
    
    # Ensure base output directory exists
    if not os.path.exists(args.out):
        os.makedirs(args.out)

    downloaded_files = []

    # --- EXECUTE MODES ---
    
    if args.mode == 'general':
        print(f"--- GENERAL MINING (Limit: {args.limit}) ---")
        for metal in MINING_METALS:
            print(f"\nScanning for {metal}...")
            # For general mining, we trust the search result, so we download DIRECTLY to the metal folder
            # This saves parsing time, although parsing is safer. 
            # Let's stick to the direct download for efficiency as per original design, 
            # BUT we should ensure the files are valid.
            
            ids = get_pdb_ids_for_ligand(metal, args.limit)
            if not ids:
                print(f"No structures found for {metal}")
                continue
                
            print(f"Found {len(ids)} structures. Downloading...")
            metal_dir = os.path.join(args.out, metal)
            
            for pdb_id in ids:
                f = download_file(pdb_id, metal_dir, args.format)
                if f: downloaded_files.append(f)
                
    elif args.mode == 'curated':
        print(f"--- CURATED DATASETS ---")
        # Download to a temporary "Received" folder first, then sort
        recv_dir = os.path.join(args.out, "_Received_Curated")
        
        for category, pdb_list in HAND_PICKED_PDBS.items():
            print(f"Processing Category: {category}")
            for pdb_id in pdb_list:
                f = download_file(pdb_id, recv_dir, args.format)
                if f: downloaded_files.append(f)
        
        print("\nSorting Curated files by Ion...")
        for f in downloaded_files:
            # Use ALL known metals for sorting curated/mmp sets to catch exotic ones
            sort_into_ions(f, args.out, ALL_KNOWN_METALS)
            
        # Cleanup temporary folder
        if os.path.exists(recv_dir):
            shutil.rmtree(recv_dir)
            print(f"Cleaned up temporary directory: {recv_dir}")
        
    elif args.mode == 'mmp':
        print(f"--- METALLOPROTEINASES (MMP) ---")
        ids = get_metalloproteinases()
        print(f"Found {len(ids)} MMP structures.")
        
        recv_dir = os.path.join(args.out, "_Received_MMP")
        
        for pdb_id in ids:
            f = download_file(pdb_id, recv_dir, args.format)
            if f: downloaded_files.append(f)
            
        print("\nSorting MMP files by Ion...")
        for f in downloaded_files:
            # Use ALL known metals for sorting curated/mmp sets
            sort_into_ions(f, args.out, ALL_KNOWN_METALS)

        # Cleanup temporary folder
        if os.path.exists(recv_dir):
            shutil.rmtree(recv_dir)
            print(f"Cleaned up temporary directory: {recv_dir}")

    print("\nDataset Download & Sorting Complete.")

if __name__ == "__main__":
    main()

