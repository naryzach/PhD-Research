import argparse
import os
import json
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch, Selection, PDBIO, Select

def find_metal_sites(pdb_file, cutoff=4.0, metals=None, output_dir=None):
    if metals is None:
        metals = ["ZN", "CU", "FE", "MG", "CA", "NI", "CO", "MN"]
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    
    metal_atoms = []
    protein_atoms = []
    
    for atom in structure.get_atoms():
        # Check if it is a heteroatom residue (often metals are HETATMs)
        res_id_full = atom.get_parent().id
        is_hetatm = res_id_full[0] != " "
        
        # Determine if it is likely a metal
        is_metal = False
        
        # 1. Trust Element field if populated
        if atom.element.strip():
            if atom.element.upper() in metals:
                is_metal = True
        else:
            # 2. Fallback to name checking only if element is missing
            name_clean = ''.join(filter(str.isalpha, atom.name.upper()))
            if name_clean in metals:
                # Disambiguate CA if needed (though is_hetatm usually handles it)
                if name_clean == "CA" and not is_hetatm:
                    is_metal = False
                else:
                    is_metal = True
        
        # Debugging metals
        # if is_hetatm:
        #      print(f"DEBUG: {atom.name} Elem:{atom.element} Metal:{is_metal} Flag:{res_id_full}")

        if is_metal and is_hetatm:
             metal_atoms.append(atom)
        elif not is_hetatm: # Standard residues have empty heteroflag
            protein_atoms.append(atom)
            
    if not metal_atoms:
        print(f"No metal atoms found in {pdb_file}")
        return []

    ns = NeighborSearch(protein_atoms)
    
    sites = []
    
    for metal in metal_atoms:
        site_info = {
            "metal": metal.element,
            "chain": metal.get_parent().get_parent().id,
            "residue_id": metal.get_parent().id[1],
            "atom_name": metal.name,
            "nearby_residues": []
        }
        
        nearby = ns.search(metal.coord, cutoff, level="R")
        for res in nearby:
            # Calculate distance specific to atomic contacts if needed, but residue level is good for selection
            min_dist = min([atom - metal for atom in res])
            
            res_info = {
                "chain": res.get_parent().id,
                "resname": res.get_resname(),
                "resid": res.id[1],
                "distance": float(min_dist)
            }
            site_info["nearby_residues"].append(res_info)
            
        sites.append(site_info)
        
    return sites

def generate_pymol_script(sites, pdb_file, output_path):
    with open(output_path, 'w') as f:
        f.write(f"load {os.path.abspath(pdb_file)}, protein\n")
        f.write("hide everything, protein\n")
        f.write("show cartoon, protein\n")
        f.write("color gray80, protein\n")
        
        for i, site in enumerate(sites):
            metal_sel = f"chain {site['chain']} and resi {site['residue_id']} and name {site['atom_name']}"
            f.write(f"select metal_{i}, {metal_sel}\n")
            f.write(f"show spheres, metal_{i}\n")
            f.write(f"color marine, metal_{i}\n")
            
            res_sels = []
            for res in site['nearby_residues']:
                res_sels.append(f"(chain {res['chain']} and resi {res['resid']}")
            
            if res_sels:
                combined_sel = " or ".join(res_sels)
                f.write(f"select site_{i}_residues, {combined_sel}\n")
                f.write(f"show sticks, site_{i}_residues\n")
                f.write(f"util.cbag site_{i}_residues\n")
                f.write(f"dist metal_{i}_contacts, metal_{i}, site_{i}_residues, mode=2\n")

def main():
    parser = argparse.ArgumentParser(description="Find metal binding sites in a PDB file.")
    parser.add_argument("pdb_file", help="Path to input PDB file")
    parser.add_argument("--cutoff", type=float, default=4.0, help="Distance cutoff for finding residues")
    parser.add_argument("--metals", nargs="+", help="Method names to look for (default: ZN CU FE MG CA NI CO MN)")
    parser.add_argument("--out_json", help="Path to save JSON report")
    parser.add_argument("--out_pymol", help="Path to save PyMOL visualization script")
    
    args = parser.parse_args()
    
    sites = find_metal_sites(args.pdb_file, args.cutoff, args.metals)
    
    print(json.dumps(sites, indent=2))
    
    if args.out_json:
        with open(args.out_json, 'w') as f:
            json.dump(sites, f, indent=2)
            
    if args.out_pymol:
        generate_pymol_script(sites, args.pdb_file, args.out_pymol)

if __name__ == "__main__":
    main()
