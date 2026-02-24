import os
import requests
import json
import time

# --- CONFIGURATION ---

# 1. GENERAL MINING: Broad search for these ions
TARGET_METALS = {
    "Common": ["FE", "CU", "ZN", "MN", "NI", "CO", "MG", "CA"],
    "Rare_Earths": ["LA", "CE", "ND", "EU", "GD", "TB", "DY", "YB", "LU", "Y"], 
    "Actinides": ["U", "URI", "PU", "TH"], 
    "Toxic": ["HG", "CD", "PB", "AS"]
}

# 2. HAND-PICKED: The "Greatest Hits" of metal binding
# These will be downloaded into a special "Curated" folder.
HAND_PICKED_PDBS = {
    "Siderophores_Iron": [
        "1EFD", # FhuD (Hydroxamate siderophore binding)
        "6Z2N", # PfeA (Enterobactin receptor - strongest iron binder)
        "1FHA", # Ferritin (Iron storage cage)
        "2IAH"  # FhuA (Ferrichrome transporter)
    ],
    "Copper_Machinery": [
        "4AZU", # Azurin (Electron transfer)
        "1PLC", # Plastocyanin
        "7ZC3", # Atox1 (Copper chaperone)
        "1CC8"  # Cytochrome c Oxidase (Cu center)
    ],
    "Metallothioneins_Sponges": [
        "4MT2", # Rat Metallothionein-2
        "1MHU", # Human Metallothionein-2
        "2FJ5", # Human Metallothionein-3
        "1DME"  # Crab Metallothionein-1
    ],
    "Lanthanide_Specialists": [
        "6MI5", # Lanmodulin (Yttrium bound)
        "8DQ2", # Lanmodulin (Lanthanum bound)
        "9C8Z", # LanD (New lanthanide chaperone)
        "7CCO"  # Engineered Lanthanide Binding Tag (LBT)
    ],
    "Actinide_Engineered": [
        "4GLW", # Super Uranyl-binding Protein (SUP)
        "4FZP"  # SUP variant
    ]
}

# Output settings
BASE_DIR = "protein_metal_mine_v2"
LIMIT_PER_METAL = 5 
LIMIT_MMP = 15 # Higher limit for the specialized MMP search

# URLs
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DOWNLOAD_URL = "https://files.rcsb.org/download/{}.cif"

def download_file(pdb_id, folder_path):
    """Downloads the .cif file for a given PDB ID."""
    url = DOWNLOAD_URL.format(pdb_id)
    file_path = os.path.join(folder_path, f"{pdb_id}.cif")
    
    if os.path.exists(file_path):
        print(f"    [Skipped] {pdb_id} exists.")
        return

    try:
        r = requests.get(url)
        if r.status_code == 200:
            with open(file_path, 'wb') as f:
                f.write(r.content)
            print(f"    [Downloaded] {pdb_id}")
        else:
            print(f"    [Failed] {pdb_id} (Status: {r.status_code})")
    except Exception as e:
        print(f"    [Error] {pdb_id}: {e}")

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

def get_pdb_ids_for_ligand(ligand_id):
    """Standard metal search."""
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
        "request_options": {"paginate": {"start": 0, "rows": LIMIT_PER_METAL}}
    }
    return search_rcsb(query)

def get_metalloproteinases():
    """
    Specifically targets Matrix Metalloproteinases (MMPs) and ADAMs 
    that contain Zinc (their active cofactor).
    """
    print(f"\n--- HUNTING METALLOPROTEINASES (MMPs & ADAMs) ---")
    
    # Complex query: (Title contains MMP or ADAM) AND (Has Zinc)
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
        "request_options": {"paginate": {"start": 0, "rows": LIMIT_MMP}}
    }
    
    ids = search_rcsb(query)
    
    # Create folder
    mmp_dir = os.path.join(BASE_DIR, "Specialist_MMPs")
    if not os.path.exists(mmp_dir):
        os.makedirs(mmp_dir)
        
    print(f"Found {len(ids)} Metalloproteinases with Zinc.")
    for pdb_id in ids:
        download_file(pdb_id, mmp_dir)

def main():
    if not os.path.exists(BASE_DIR):
        os.makedirs(BASE_DIR)

    # 1. Download Hand-Picked Curated Set
    print("--- DOWNLOADING CURATED SET ---")
    for category, pdb_list in HAND_PICKED_PDBS.items():
        print(f"Processing Curated: {category}")
        cat_dir = os.path.join(BASE_DIR, "Curated_Best_Of", category)
        if not os.path.exists(cat_dir):
            os.makedirs(cat_dir)
        
        for pdb_id in pdb_list:
            download_file(pdb_id, cat_dir)

    # 2. Download Metalloproteinases
    get_metalloproteinases()

    # 3. General Mining (as before)
    print("\n--- GENERAL METAL MINING ---")
    for category, metals in TARGET_METALS.items():
        print(f"Category: {category}")
        for metal in metals:
            metal_dir = os.path.join(BASE_DIR, "General_Mining", category, metal)
            if not os.path.exists(metal_dir):
                os.makedirs(metal_dir)
            
            pdb_ids = get_pdb_ids_for_ligand(metal)
            if pdb_ids:
                print(f"  {metal}: Found {len(pdb_ids)} structures.")
                for pdb_id in pdb_ids:
                    download_file(pdb_id, metal_dir)
                    time.sleep(0.2)
            else:
                print(f"  {metal}: No structures found.")

    print("\n--- COMPLETE ---")
    print(f"Check '{BASE_DIR}' for Curated, Specialist, and General folders.")

if __name__ == "__main__":
    main()