import argparse
import os
import glob
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select
import warnings
from find_metal_sites import find_metal_sites  # Import the function from our sibling script

# Suppress BioPython warnings
warnings.filterwarnings("ignore")

class SiteSelect(Select):
    def __init__(self, residues, metal_chain, metal_id):
        self.residues = residues # List of dicts with chain/resid
        self.metal_chain = metal_chain
        self.metal_id = metal_id

    def accept_residue(self, residue):
        # Accept the metal
        if residue.get_parent().id == self.metal_chain and residue.id[1] == self.metal_id:
            return 1
        
        # Accept contact residues
        for res in self.residues:
            if residue.get_parent().id == res['chain'] and residue.id[1] == res['resid']:
                return 1
        return 0

def build_catalog(pdb_dir, out_file, scaffold_dir=None):
    results = []
    
    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))
    print(f"Scanning {len(pdb_files)} PDB files...")
    
    if scaffold_dir:
        os.makedirs(scaffold_dir, exist_ok=True)
    
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    
    for pdb_path in pdb_files:
        try:
            sites = find_metal_sites(pdb_path)
            if not sites:
                continue
            
            basename = os.path.basename(pdb_path).replace(".pdb", "")
            
            for i, site in enumerate(sites):
                entry = {
                    "pdb": basename,
                    "metal": site['metal'],
                    "metal_chain": site['chain'],
                    "metal_resid": site['residue_id'],
                    "num_contacts": len(site['nearby_residues']),
                    "residues": ";".join([f"{r['chain']}{r['resid']}({r['resname']})" for r in site['nearby_residues']])
                }
                results.append(entry)
                
                # Save Motif PDB if requested
                if scaffold_dir:
                    structure = parser.get_structure("temp", pdb_path)
                    io.set_structure(structure)
                    
                    motif_name = f"{basename}_{site['metal']}_site{i}.pdb"
                    motif_path = os.path.join(scaffold_dir, motif_name)
                    
                    selector = SiteSelect(
                        site['nearby_residues'], 
                        site['chain'], 
                        site['residue_id']
                    )
                    io.save(motif_path, select=selector)
                    
        except Exception as e:
            print(f"Error processing {pdb_path}: {e}")
            
    df = pd.DataFrame(results)
    df.to_csv(out_file, index=False)
    print(f"Catalog saved to {out_file}. Found {len(df)} sites.")

def main():
    parser = argparse.ArgumentParser(description="Build Metal Binder Catalog")
    parser.add_argument("pdb_dir", help="Directory containing PDBs to scan")
    parser.add_argument("out_csv", help="Output CSV path")
    parser.add_argument("--scaffold_dir", help="Directory to save extracted motifs for simple scaffolding")
    
    args = parser.parse_args()
    build_catalog(args.pdb_dir, args.out_csv, args.scaffold_dir)

if __name__ == "__main__":
    main()
