import argparse
import os
import requests
import json
from Bio.PDB import PDBList

def search_rcsb(query_str, max_results=20):
    """
    Search RCSB PDB API for IDs matching the query.
    Falls back to a known list if API fails.
    """
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # Simpler query
    query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": query_str
            }
        },
        "request_options": {
            "return_all_hits": False,
            "results_max": max_results
        },
        "return_type": "entry"
    }
    
    try:
        response = requests.post(url, json=query)
        response.raise_for_status()
        data = response.json()
        ids = [result['identifier'] for result in data.get('result_set', [])]
        return ids
    except Exception as e:
        print(f"Error searching RCSB: {e}")
        # Fallback for testing "Zinc"
        if "Zinc" in query_str or "zinc" in query_str:
            print("Using fallback list for Zinc...")
            return ["1A5T", "1ZEG", "2CB3", "3H4R", "4V2S"][:max_results]
        return []

def download_pdbs(pdb_ids, out_dir):
    """
    Download PDB files using BioPython.
    """
    os.makedirs(out_dir, exist_ok=True)
    pdbl = PDBList(verbose=False)
    
    print(f"Downloading {len(pdb_ids)} structures to {out_dir}...")
    
    # BioPython downloads to weird directory structures usually, let's control it.
    # actually PDBList.retrieve_pdb_file puts it in `out_dir` but might rename it `pdbXXXX.ent`.
    
    downloaded_files = []
    
    for pdb_id in pdb_ids:
        try:
            # Download as .ent file
            filename = pdbl.retrieve_pdb_file(pdb_id, pdir=out_dir, file_format="pdb")
            
            # PDBList usually renames file to 'pdbXXXX.ent', return value is the path.
            # We want 'XXXX.pdb'
            final_path = os.path.join(out_dir, f"{pdb_id}.pdb")
            
            if os.path.exists(filename):
                # Rename if needed (BioPython might handle this slightly differently across versions)
                # Usually returns path to file.
                if filename.endswith(".ent"):
                    os.rename(filename, final_path)
                    downloaded_files.append(final_path)
                elif filename.endswith(".pdb"):
                     downloaded_files.append(filename)
                
        except Exception as e:
            print(f"Failed to download {pdb_id}: {e}")

    return downloaded_files

def main():
    parser = argparse.ArgumentParser(description="Download Metal Binding PDBs form RCSB")
    parser.add_argument("--query", default="Zinc binding", help="Text query for RCSB")
    parser.add_argument("--out", default="../Local/Metal_PDBs", help="Output directory")
    parser.add_argument("--limit", type=int, default=10, help="Max number of PDBs")
    
    args = parser.parse_args()
    
    print(f"Searching for: {args.query}")
    ids = search_rcsb(args.query, args.limit)
    
    if ids:
        print(f"Found {len(ids)} IDs: {', '.join(ids[:5])}...")
        download_pdbs(ids, args.out)
    else:
        print("No results found.")

if __name__ == "__main__":
    main()
