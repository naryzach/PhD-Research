import argparse
import os
import glob
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

def get_parser(file_path):
    if file_path.endswith(".cif"):
        return MMCIFParser(QUIET=True)
    else:
        return PDBParser(QUIET=True)

class SiteSelect(Select):
    def __init__(self, residues, metal_atom):
        self.residues = residues # List of (chain, id)
        self.metal_atom = metal_atom

    def accept_residue(self, residue):
        # Keep if in list
        if (residue.get_parent().id, residue.id[1]) in self.residues:
            return 1
        # Keep if it is the metal (residue level check is tricky for single atom, usually metal is its own residue)
        if residue == self.metal_atom.get_parent():
            return 1
        return 0

def find_sites_in_structure(structure, metal_types, cutoff=6.0):
    """
    Find metal sites.
    cutoff: 6.0A to capture sufficient context (first and second shell).
    """
    sites = []
    
    atoms = list(structure.get_atoms())
    metals = [a for a in atoms if a.element.upper() in metal_types]
    
    for metal in metals:
        # Find neighbors
        nearby = []
        metal_res = metal.get_parent()
        
        for atom in atoms:
            if atom == metal: continue
            if atom.get_parent() == metal_res: continue # Skip same residue
            
            # Distance check
            dist = atom - metal
            if dist < cutoff:
                # Check if AA
                res = atom.get_parent()
                if is_aa(res, standard=True):
                    nearby.append(res)
        
        # Filter unique residues
        nearby = list(set(nearby))
        
        if len(nearby) >= 2: # At least 2 residues (bidentate+)
            site_data = {
                'metal_id': f"{metal_res.get_parent().id}_{metal_res.id[1]}",
                'metal_element': metal.element.upper(),
                'residues': [(r.get_parent().id, r.id[1], r.resname) for r in nearby],
                'metal_atom_obj': metal 
            }
            sites.append(site_data)
            
    return sites

def save_motif(structure, site, out_path):
    # Select residues and metal
    # Note: Bio.PDB structure is shared, so we need a selector class
    io = PDBIO()
    io.set_structure(structure)
    
    # Residue identifier list for Selector
    target_res = [(r[0], r[1]) for r in site['residues']]
    selector = SiteSelect(target_res, site['metal_atom_obj'])
    
    io.save(out_path, selector)

def catalog_directory(input_dir, output_dir, metals):
    catalog_data = []
    
    # Support both .pdb and .cif
    files = glob.glob(os.path.join(input_dir, "*"))
    files = [f for f in files if f.endswith(".pdb") or f.endswith(".cif")]
    
    print(f"Scanning {len(files)} files in {input_dir}...")
    
    for file_path in files:
        basename = os.path.basename(file_path).replace(".pdb", "").replace(".cif", "")
        parser = get_parser(file_path)
        
        try:
            structure = parser.get_structure(basename, file_path)
        except Exception as e:
            print(f"Failed to parse {basename}: {e}")
            continue
            
        sites = find_sites_in_structure(structure, metals)
        
        if sites:
            # Create motif dir for this Ion if processed individually? 
            # Or assume input_dir is Ion specific.
            # Let's save to output_dir directly.
            
            for i, site in enumerate(sites):
                motif_id = f"{basename}_{site['metal_element']}_site{i}"
                motif_path = os.path.join(output_dir, f"{motif_id}.pdb")
                
                # Save motif PDB
                try:
                    save_motif(structure, site, motif_path)
                except Exception as e:
                    print(f"Failed to save motif {motif_id}: {e}")
                    continue
                
                # Add to Catalog
                res_str = ";".join([f"{r[0]}{r[1]}({r[2]})" for r in site['residues']])
                catalog_data.append({
                    'motif_id': motif_id,
                    'source_pdb': basename,
                    'ion': site['metal_element'],
                    'num_residues': len(site['residues']),
                    'residues': res_str,
                    'path': os.path.abspath(motif_path)
                })
                
    return catalog_data

def main():
    parser = argparse.ArgumentParser(description="Catalog Metal Binding Sites")
    parser.add_argument("--pdb_dir", required=True, help="Input directory containing PDB/CIFs (can be root of Ion folders)")
    parser.add_argument("--out_dir", required=True, help="Output directory for Catalog and Motifs")
    parser.add_argument("--metals", default="ZN,CU,NI,CO,MN,FE,MG,CA", help="Metals to filter for")
    
    args = parser.parse_args()
    metals_list = [m.strip().upper() for m in args.metals.split(',')]
    
    # Check if pdb_dir has subdirectories (Ion folders)
    # We want to maintain structure: Local/Metal_Catalog/{Ion}/...
    
    # If pdb_dir has 'ZN', 'CU' etc inside, we recurse.
    # Simpler: Just walk.
    
    full_catalog = []
    
    # If users point to 'Local/Metal_PDBs', we iterate subfolders matching metals
    subdirs = [d for d in os.listdir(args.pdb_dir) if os.path.isdir(os.path.join(args.pdb_dir, d))]
    
    # Filter subdirs if they match our metal list? Or just process all.
    # Let's process all subdirs that look like valid folders.
    
    # If no subdirs, process root.
    dirs_to_process = []
    if not subdirs:
        dirs_to_process = [(args.pdb_dir, "Unknown")]
    else:
        for d in subdirs:
            dirs_to_process.append((os.path.join(args.pdb_dir, d), d)) # path, ion_name
            
    for inp, ion_name in dirs_to_process:
        print(f"Processing folder: {ion_name}")
        
        # Output folder for motifs
        motif_out = os.path.join(args.out_dir, ion_name)
        os.makedirs(motif_out, exist_ok=True)
        
        # Run catalog
        # Filter: Only look for the specific ion if the folder is named after it? 
        # Or look for *all* metals in *all* folders?
        # Safe bet: Look for all requested metals.
        
        data = catalog_directory(inp, motif_out, metals_list)
        full_catalog.extend(data)
        
    # Save CSV
    if full_catalog:
        df = pd.DataFrame(full_catalog)
        csv_path = os.path.join(args.out_dir, "master_catalog.csv")
        df.to_csv(csv_path, index=False)
        print(f"Catalog saved to {csv_path} with {len(df)} entries.")
    else:
        print("No sites found.")

if __name__ == "__main__":
    main()
