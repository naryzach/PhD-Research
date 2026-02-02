import argparse
import os
import sys
import json
import subprocess
import shutil
from Bio.PDB import PDBParser, PDBIO, Select

def swap_metal_in_pdb(input_pdb, output_pdb, target_metal, metal_site_info):
    """
    Reads the input PDB, finds the metal atom at the specified site, 
    changes its element/name to target_metal, and saves the new PDB.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("temp", input_pdb)
    
    metal_chain = metal_site_info['metal_chain']
    metal_resid = int(metal_site_info['metal_resid'])
    
    swapped = False
    for atom in structure.get_atoms():
        # Match the metal
        if (atom.get_parent().id[1] == metal_resid and 
            atom.get_parent().get_parent().id == metal_chain):
            
            # Simple swap: change element and name
            # Note: For proper forcefields, might need more, but for RFdiffusion inference, element is key.
            atom.element = target_metal.upper()
            atom.name = target_metal.upper() + "+2" # Generic ion name
            swapped = True
            
    if not swapped:
        print(f"Warning: Could not find metal at {metal_chain}{metal_resid} to swap.")
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    return swapped

def run_metal_swap(args):
    # Absolute paths
    pdb_path = os.path.abspath(args.pdb)
    base_out_dir = os.path.abspath(args.out)
    rf_path = os.path.abspath(args.rf_path)
    lmpnn_path = os.path.abspath(args.lmpnn_path)
    
    # Load sites
    with open(args.sites_json, 'r') as f:
        sites = json.load(f)
    
    # Filter sites if user requested specific one? 
    # For now, let's assume the user wants to process *all* valid sites in the file unless specified.
    # Or just the first one for simplicity.
    target_site = sites[0] # Default to first found site
    print(f"Targeting metal site: {target_site['metal']} at {target_site['chain']}{target_site['residue_id']}")
    
    # Define contiguous mask based on binding residues
    # We want to redesign the residues *interacting* with the metal.
    # The JSON has 'nearby_residues'.
    # We need to turn this into a contig string for RFdiffusion.
    # RFdiffusion contig logic: [FIXED/DESIGN/FIXED]
    # If we have scattered residues (e.g. H166, H170, H176), we can't easily redesign *just* those dots efficiently with contigmap alone 
    # without making many segments. 
    # Better approach: Find the bounding box or 'loops' covering these residues.
    
    # Simple logic: Group residues into ranges.
    # e.g. 166, 170, 176 -> Range 166-176.
    residues = sorted([int(r['resid']) for r in target_site['nearby_residues']])
    
    if not residues:
        print("No residues found for this site.")
        return

    # Add padding?
    start_res = max(1, min(residues) - args.padding)
    end_res = max(residues) + args.padding
    chain = target_site['chain'] # Assuming single chain binding for now
    
    print(f"Defined design region: Chain {chain}, Residues {start_res}-{end_res}")
    
    metals_to_test = args.metals.split(',')
    
    for metal in metals_to_test:
        metal = metal.strip().upper()
        print(f"\n--- Processing Target Metal: {metal} ---")
        
        metal_out_dir = os.path.join(base_out_dir, metal)
        os.makedirs(metal_out_dir, exist_ok=True)
        
        # 1. Swap Metal in PDB
        swapped_pdb = os.path.join(metal_out_dir, f"input_{metal}.pdb")
        site_info_for_swap = {
           'metal_chain': target_site['chain'], # Note: found_sites puts metal chain here. Verify if it's 'metal chain' or 'protein chain' in JSON.
           'metal_resid': target_site['residue_id'] 
        }
        # Actually in find_metal_sites, 'chain' is the metal's chain.
        
        swap_metal_in_pdb(pdb_path, swapped_pdb, metal, site_info_for_swap)
        
        # 2. Run RFdiffusion (Partial Diffusion / Inpainting)
        rf_out = os.path.join(metal_out_dir, "rfdiffusion")
        os.makedirs(rf_out, exist_ok=True)
        
        design_len = end_res - start_res + 1
        # Contig: FIXED_BEFORE / DESIGN_LEN / FIXED_AFTER
        # [A1-(start-1)/len/A(end+1)-]
        # We need to handle the metal chain too.
        
        metal_chain_str = f" {site_info_for_swap['metal_chain']}{site_info_for_swap['metal_resid']}-{site_info_for_swap['metal_resid']}"
        # Note: If metal is part of protein chain, this logic changes. Assuming HetAtom separate chain or same chain?
        # If 'B' is metal chain.
        # Contig: [A1-start/len/A-end B1-1]
        
        protein_chain = target_site['nearby_residues'][0]['chain'] # Get protein chain from residues
        
        contig = f"[{protein_chain}1-{start_res-1}/{design_len}/{protein_chain}{end_res+1}-{metal_chain_str}]"
        
        print(f"Running RFdiffusion with contig: {contig}")
        
        cmd_rf = [
            sys.executable, os.path.join(rf_path, 'scripts/run_inference.py'),
            f"inference.input_pdb={swapped_pdb}",
            f"inference.output_prefix={os.path.join(rf_out, 'design')}",
            f"inference.num_designs={args.num_designs}",
            f"contigmap.contigs={contig}"
        ]
        
        # Execute RFdiffusion
        if not args.dry_run:
            try:
                subprocess.run(cmd_rf, check=True)
            except subprocess.CalledProcessError as e:
                print(f"RFdiffusion failed for {metal}: {e}")
                continue
        else:
            print(f"Dry run: {' '.join(cmd_rf)}")

        # 3. Run LigandMPNN
        lmpnn_out = os.path.join(metal_out_dir, "ligandmpnn")
        os.makedirs(lmpnn_out, exist_ok=True)
        
        print(f"Running LigandMPNN for {metal}...")
        
        cmd_lmpnn = [
            sys.executable, os.path.join(lmpnn_path, 'run.py'),
            "--pdb_path", rf_out,
            "--out_folder", lmpnn_out,
            "--model_type", "ligand_mpnn",
            "--num_seq_per_target", "2",
            "--batch_size", "1"
        ]
        
        if not args.dry_run:
            try:
                subprocess.run(cmd_lmpnn, check=True)
            except subprocess.CalledProcessError as e:
                print(f"LigandMPNN failed for {metal}: {e}")
                continue
        else:
            print(f"Dry run: {' '.join(cmd_lmpnn)}")
            
    print("\nAll Metal Swaps Completed.")

def main():
    parser = argparse.ArgumentParser(description="Automated Metal Swap & Redesign")
    parser.add_argument("--pdb", required=True, help="Input PDB file")
    parser.add_argument("--sites_json", required=True, help="JSON file from find_metal_sites.py")
    parser.add_argument("--metals", default="ZN,CU,NI,CO,MN,FE,MG,CA", help="Comma-separated list of metals to generate binders for")
    parser.add_argument("--metal_chain", help="Chain ID of the metal in the input PDB (e.g. B, Z)")
    parser.add_argument("--out", default="../Local/Metal_Outputs", help="Output directory")
    parser.add_argument("--padding", type=int, default=2, help="Padding residues around the binding site to redesign")
    parser.add_argument("--num_designs", type=int, default=1, help="Number of RFdiffusion designs per metal")
    parser.add_argument("--dry_run", action="store_true", help="Print commands without running")
    
    # Defaults
    parser.add_argument("--rf_path", default="../Tools/RFdiffusion", help="Path to RFdiffusion")
    parser.add_argument("--lmpnn_path", default="../Tools/LigandMPNN", help="Path to LigandMPNN")
    
    args = parser.parse_args()
    run_metal_swap(args)

if __name__ == "__main__":
    main()
