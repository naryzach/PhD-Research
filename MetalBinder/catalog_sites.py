import argparse
import os
import glob
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select, NeighborSearch
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
    Find metal sites using KD-tree for speed.
    cutoff: 6.0A.
    """
    sites = []
    
    atoms = list(structure.get_atoms())
    metals = [a for a in atoms if a.element.upper() in metal_types]
    
    if not metals:
        return []
        
    # Build KD-tree for all atoms
    ns = NeighborSearch(atoms)
    
    for metal in metals:
        metal_res = metal.get_parent()
        
        # Search for nearby residues (level='R' returns residues)
        nearby_residues = ns.search(metal.coord, cutoff, level='R')
        
        # Filter:
        # 1. Remove the metal's own residue
        # 2. Keep only Amino Acids (is_aa)
        
        valid_neighbors = []
        for res in nearby_residues:
            if res == metal_res:
                continue
            if is_aa(res, standard=True):
                valid_neighbors.append(res)
        
        if len(valid_neighbors) >= 2:
            site_data = {
                'metal_id': f"{metal_res.get_parent().id}_{metal_res.id[1]}",
                'metal_element': metal.element.upper(),
                'residues': [(r.get_parent().id, r.id[1], r.resname) for r in valid_neighbors],
                'metal_atom_obj': metal 
            }
            sites.append(site_data)
            
    return sites

from Bio.PDB import Structure, Model, Chain, Residue, Atom

def save_motif(structure, site, out_path):
    # Create a new structure for the motif to handle remapping of chain IDs
    # (since some source PDBs have chain IDs > 1 char which breaks PDB format)
    
    new_struct = Structure.Structure("motif")
    model = Model.Model(0)
    new_struct.add(model)
    
    # Track used chain IDs to assign new single-char IDs
    # We will simply assign A, B, C... to the chains involved
    chain_map = {} # old_chain_id -> new_chain_id
    used_new_ids = set()
    possible_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
    next_id_idx = 0
    
    # Collect all atoms to save: residues + metal
    atoms_to_save = []
    
    # Add metal
    atoms_to_save.append(site['metal_atom_obj'])
    
    # Add residues
    # We need to find the specific residues in the structure again or just use the ones we found
    # The 'site' dict stores residue info, but we need the actual objects 
    # Luckily 'site' was constructed from 'structure', so we can look them up, 
    # BUT we didn't store residue objects, only IDs in 'residues' list.
    # Wait, 'find_sites_in_structure' logic:
    # 'residues': [(r.get_parent().id, r.id[1], r.resname) for r in nearby]
    # We need to re-fetch these residues from the structure.
    
    target_res_ids = set([(r[0], r[1]) for r in site['residues']]) # (chain_id, res_seq)
    
    # Helper to get residue object
    def get_residue(chain_id, res_seq):
        chain = structure[0][chain_id] # Assuming model 0
        return chain[res_seq]

    residue_objs = []
    try:
        for chain_id, res_seq, _ in site['residues']:
            residue_objs.append(get_residue(chain_id, res_seq))
            
        # Also need the metal's residue (parent)
        metal_res = site['metal_atom_obj'].get_parent()
        if metal_res not in residue_objs:
             residue_objs.append(metal_res)
             
    except KeyError:
        # Fallback if lookup fails (shouldn't happen)
        print("Warning: Could not lookup residues for motif.")
        return

    # Now reconstruct structure
    # We group by chain to keep residue connectivity if possible, 
    # but for a motif (disjoint residues usually), it matters less.
    # However, to be nice, let's map source chains to new chains.
    
    for res in residue_objs:
        old_chain_id = res.get_parent().id
        
        if old_chain_id not in chain_map:
            if next_id_idx < len(possible_ids):
                new_id = possible_ids[next_id_idx]
                next_id_idx += 1
            else:
                new_id = "X" # Fallback if we somehow have >36 chains in a tiny motif
            chain_map[old_chain_id] = new_id
            
            # Create new chain
            new_chain = Chain.Chain(new_id)
            model.add(new_chain)
            
        new_chain = model[chain_map[old_chain_id]]
        
        # Copy residue
        new_res = res.copy() # This copies atoms too
        new_chain.add(new_res)
        
    # Save
    io = PDBIO()
    io.set_structure(new_struct)
    io.save(out_path)

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
