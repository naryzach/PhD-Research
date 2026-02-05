import argparse
import os
import requests
import json
from Bio.PDB import PDBList

# Standard metals list for v2.0
DEFAULT_METALS = ["ZN", "CU", "NI", "CO", "MN", "FE", "MG", "CA"]

def search_rcsb_by_ligand(ligand_id, max_results=100):
    """
    Search RCSB PDB API for IDs containing a specific ligand (metal).
    """
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # Advanced query for "Chemical ID = ligand_id"
    # Using 'rcsb_chem_comp_container_identifiers.comp_id' which is the standard searchable attribute for chemical components.
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_chem_comp_container_identifiers.comp_id",
                "operator": "exact_match",
                "value": ligand_id
            }
        },
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": max_results
            }
        },
        "return_type": "entry"
    }
    
    # print(json.dumps(query, indent=2))
    
    try:
        response = requests.post(url, json=query)
        response.raise_for_status()
        data = response.json()
        ids = [result['identifier'] for result in data.get('result_set', [])]
        return ids
    except Exception as e:
        print(f"Error searching RCSB for {ligand_id}: {e}")
        # if 'response' in locals():
        #    print(f"RCSB Error Response: {response.text}")
            
        print("Attempting fallback text search...")
        
        # Fallback: Simple text search
        fallback_query = {
            "query": {
                "type": "terminal",
                "service": "full_text",
                "parameters": {
                    "value": f"{ligand_id} binding"
                }
            },
            "request_options": {
                "paginate": {
                    "start": 0,
                    "rows": max_results
                }
            },
            "return_type": "entry"
        }
        try:
             response = requests.post(url, json=fallback_query)
             response.raise_for_status()
             data = response.json()
             ids = [result['identifier'] for result in data.get('result_set', [])]
             return ids
        except Exception as e2:
             print(f"Fallback failed: {e2}")
             if ligand_id == "ZN":
                 return ["1A5T", "1ZEG"]
             return []

# ... in main ...
        print(f"Found {len(ids)} IDs: {ids}. Downloading...")

def download_files(pdb_ids, out_dir, file_format="cif"):
    """
    Download PDB/CIF files using direct HTTP from RCSB.
    """
    os.makedirs(out_dir, exist_ok=True)
    
    print(f"Downloading {len(pdb_ids)} structures to {out_dir} (Format: {file_format})...")
    
    downloaded_files = []
    
    for pdb_id in pdb_ids:
        pdb_id = pdb_id.upper()
        ext = "cif" if file_format == "cif" else "pdb"
        filename = f"{pdb_id}.{ext}"
        final_path = os.path.join(out_dir, filename)
        
        # Check existence
        if os.path.exists(final_path) and os.path.getsize(final_path) > 0:
            print(f"  {filename} exists. Skipping.")
            downloaded_files.append(final_path)
            continue
            
        # Download
        url = f"https://files.rcsb.org/download/{filename}"
        try:
            r = requests.get(url)
            if r.status_code == 200:
                with open(final_path, 'w') as f:
                    f.write(r.text)
                downloaded_files.append(final_path)
                print(f"  Downloaded {filename}")
            else:
                print(f"  Failed to download {filename} (Status: {r.status_code})")
        except Exception as e:
            print(f"  Error downloading {filename}: {e}")
            
    return downloaded_files

def main():
    parser = argparse.ArgumentParser(description="Download Metal Binding PDBs/CIFs")
    parser.add_argument("--metals", default=",".join(DEFAULT_METALS), help="Comma-separated list of metals")
    parser.add_argument("--out", default="../Local/Metal_PDBs", help="Base output directory")
    parser.add_argument("--limit", type=int, default=10, help="Max number of files per metal")
    parser.add_argument("--format", choices=["cif", "pdb"], default="cif", help="File format to download")
    
    args = parser.parse_args()
    
    metals = [m.strip().upper() for m in args.metals.split(',')]
    
    for metal in metals:
        print(f"\n--- Searching for {metal} Binders ---")
        ids = search_rcsb_by_ligand(metal, args.limit)
        
        if not ids:
            print(f"No results found for {metal}.")
            continue
            
        print(f"Found {len(ids)} IDs. Downloading...")
        
        metal_out_dir = os.path.abspath(os.path.join(args.out, metal))
        download_files(ids, metal_out_dir, args.format)
        
    print("\nDataset Download Complete.")

if __name__ == "__main__":
    main()
