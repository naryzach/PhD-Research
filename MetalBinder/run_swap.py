import argparse
import os
import sys
import pandas as pd
import subprocess
import glob
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

# Import shared utilities
try:
    from utils_rfd2 import find_model_weights, add_ori_token, get_contig_str
except ImportError:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from utils_rfd2 import find_model_weights, add_ori_token, get_contig_str

def get_parser(file_path):
    if file_path.endswith(".cif"):
        return MMCIFParser(QUIET=True)
    return PDBParser(QUIET=True)

def get_chain_residues(structure, chain_id):
    """Get sorted list of standard AA residues in the chain."""
    residues = []
    for model in structure:
        if chain_id in model:
            for r in model[chain_id]:
                if is_aa(r, standard=True):
                    residues.append(r.id[1])
            break 
    return sorted(list(set(residues)))

def run_swap(args):
    if not os.path.exists(args.catalog_csv):
        print(f"Catalog CSV not found: {args.catalog_csv}")
        return

    catalog = pd.read_csv(args.catalog_csv)
    target_metals = [m.strip().upper() for m in args.target_metals.split(',')]
    
    project_root = os.path.dirname(os.path.abspath(__file__))
    
    print(f"Loaded {len(catalog)} sites from catalog.")
    
    # Pre-find checkpoints
    ckpt_path = find_model_weights(args.rf_path)
    if ckpt_path:
        print(f"Using model weights: {os.path.basename(ckpt_path)}")
    
    rf_script = os.path.join(os.path.abspath(args.rf_path), "rf_diffusion", "run_inference.py")
    
    for idx, row in catalog.iterrows():
        source_pdb = row['source_pdb'] # E.g. "1ZEG"
        source_ion = row['ion']        # E.g. "ZN"
        motif_id = row['motif_id']     # E.g. "1ZEG_ZN_site0"
        
        # Locate Source PDB
        src_path_pdb = os.path.join(args.pdb_dir, source_ion, f"{source_pdb}.pdb")
        src_path_cif = os.path.join(args.pdb_dir, source_ion, f"{source_pdb}.cif")
        
        if os.path.exists(src_path_pdb):
            src_path = src_path_pdb
        elif os.path.exists(src_path_cif):
            src_path = src_path_cif
        else:
            print(f"Source file not found for {source_pdb} ({source_ion})")
            continue
            
        print(f"\nProcessing {motif_id} (Source: {source_pdb})")
        
        parser = get_parser(src_path)
        try:
            structure = parser.get_structure(source_pdb, src_path)
        except:
             print("Failed to parse structure")
             continue

        # Remap Chains (Detect AAA -> A)
        chain_map = {}
        used_chains = set()
        
        # First pass: identify existing single-char chains
        for model in structure:
            for chain in model:
                if len(chain.id) == 1:
                    used_chains.add(chain.id)
                    
        # Second pass: remap long chains
        import string
        available_chars = [c for c in string.ascii_uppercase if c not in used_chains and c != 'Z']
        
        for model in structure:
            for chain in model:
                original_id = chain.id
                if len(original_id) > 1:
                    if original_id in chain_map:
                        new_id = chain_map[original_id]
                    else:
                        if not available_chars:
                            print(f"Error: Run out of single-character chain IDs for {source_pdb}")
                            break # Critical failure
                        new_id = available_chars.pop(0)
                        chain_map[original_id] = new_id
                        used_chains.add(new_id)
                    
                    print(f"  Remapping chain {original_id} -> {new_id}")
                    chain.id = new_id
        
        # Update internal maps if remapping occurred
        if chain_map:
             # structure is updated in place
             pass

        # Parse residues to find Target Chain
        res_str = row['residues']
        chains = set()
        residue_ids = []
        
        # Parse logic (A52;B23)
        # Parse logic (A52;B23)
        import re
        for r in res_str.split(';'):
            # Use non-greedy match for chain to allow digits in residue to be captured by group 2
            match = re.match(r"([A-Za-z0-9]+?)(\d+)", r.split('(')[0])
            if match:
                raw_chain = match.group(1)
                res_num = int(match.group(2))
                
                # Apply remapping
                final_chain = chain_map.get(raw_chain, raw_chain)
                
                chains.add(final_chain)
                residue_ids.append((final_chain, res_num))
            
        if not chains:
             print(f"  Could not parse chain from {res_str}")
             continue
             
        target_chains = sorted(list(chains))
        
        # Identify residues to design (mask)
        # Store as (chain, resnum)
        design_residues = set()
        padding = 1
        for r_c, r_n in residue_ids:
            for k in range(r_n - padding, r_n + padding + 1):
                design_residues.add((r_c, k))
        
        rf_contig_list = []
        
        # Build contig for each chain
        for chain_id in target_chains:
            chain_res_list = get_chain_residues(structure, chain_id)
            if not chain_res_list: continue
                
            contig_parts = []
            start_segment = None
            prev_r = None
            design_block_len = 0
            
            for r in chain_res_list:
                is_designable = (chain_id, r) in design_residues
                is_gap = (prev_r is not None) and (r != prev_r + 1)
                
                if is_gap:
                    if start_segment is not None:
                         contig_parts.append(f"{chain_id}{start_segment}-{prev_r}")
                         start_segment = None
                    if design_block_len > 0:
                        contig_parts.append(f"{design_block_len}-{design_block_len}")
                        design_block_len = 0
                
                if is_designable:
                    if start_segment is not None:
                        contig_parts.append(f"{chain_id}{start_segment}-{prev_r}")
                        start_segment = None
                    design_block_len += 1
                else:
                    if design_block_len > 0:
                        contig_parts.append(f"{design_block_len}-{design_block_len}")
                        design_block_len = 0
                    if start_segment is None:
                        start_segment = r
                    prev_r = r
                    
            if start_segment is not None:
                 contig_parts.append(f"{chain_id}{start_segment}-{prev_r}")
            if design_block_len > 0:
                 contig_parts.append(f"{design_block_len}-{design_block_len}")
            
            # Add padding around each chain? Or just once?
            # Original code added 10-10 padding.
            # "['10-10,{contig_raw},10-10']"
            contig_str = ",".join(contig_parts)
            rf_contig_list.append(f"10-10,{contig_str},10-10")

        # Identify Metal Atom
        atom_source_metal = None
        for model in structure:
            for chain in model:
                for residue in chain:
                     for atom in residue:
                          if atom.element.upper() == source_ion:
                               atom_source_metal = atom
                               break
                     if atom_source_metal: break
                if atom_source_metal: break
        
        if not atom_source_metal:
             print("  No metal found in structure to swap.")
             continue

        # For each target metal
        for metal in target_metals:
            out_dir = os.path.join(args.out_dir, metal, motif_id)
            os.makedirs(out_dir, exist_ok=True)
            rf_out = os.path.join(out_dir, "rfdiffusion")
            os.makedirs(rf_out, exist_ok=True)

            if not args.overwrite and glob.glob(os.path.join(rf_out, "design*.pdb")):
                print(f"    Skipping {motif_id} -> {metal} (already exists)")
                continue
            
            swapped_pdb_path = os.path.join(rf_out, "input_swapped.pdb")
            
            # Write PDB (Protein + Swapped Metal)
            with open(swapped_pdb_path, 'w') as f:
                io = PDBIO()
                class ChainSelect(Select):
                    def accept_chain(self, chain):
                        return chain.id in target_chains
                io.set_structure(structure)
                io.save(os.path.join(rf_out, "temp_prot.pdb"), ChainSelect())
                
                with open(os.path.join(rf_out, "temp_prot.pdb"), 'r') as temp:
                    for line in temp:
                        if not line.startswith("END") and not line.startswith("TER"):
                            f.write(line)
                os.remove(os.path.join(rf_out, "temp_prot.pdb"))
                
                # Write Swapped Metal
                x, y, z = atom_source_metal.get_coord()
                # Use Chain Z
                pdb_line = f"HETATM{9999:5d} {metal:>4s} {metal:>3s} Z   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {metal:>2s}\n"
                f.write(pdb_line)
                
                # Write ORI Token (Centered on Metal)
                ori_line = f"HETATM 9998  CA  ORI Z 999    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n"
                f.write(ori_line)
                
                f.flush()
                os.fsync(f.fileno())
                
            # Verify ORI presence
            with open(swapped_pdb_path, 'r') as f:
                 if "ORI" not in f.read():
                     print(f"Error: ORI token failed to write to {swapped_pdb_path}")
                     continue
            
            # Add ORI (Manually done above)
            # add_ori_token(swapped_pdb_path, [metal])
            
            # Contig: append reference to Metal Chain Z
            # RFd2 Contig string: List of strings
            current_contigs = list(rf_contig_list)
            current_contigs.append("Z1-1")
            
            # Format as string representation of list
            # "['...', '...', 'Z1-1']"
            contig_final_str = str(current_contigs).replace("'", '"')
            
            cmd = [
                sys.executable, rf_script,
                f"inference.input_pdb={os.path.abspath(swapped_pdb_path)}",
                f"inference.output_prefix={os.path.abspath(os.path.join(rf_out, 'design'))}",
                f"inference.num_designs={args.num_designs}",
                f"contigmap.contigs={contig_final_str}",
                f"inference.ligand='{metal}'"
            ]
            
            if ckpt_path:
                cmd.append(f"inference.ckpt_path={ckpt_path}")
            
            if args.dry_run:
                print(f"    Dry Run RFdiffusion ({metal}): {' '.join(cmd)}")
            else:
                try:
                    env = os.environ.copy()
                    rf_abs = os.path.abspath(args.rf_path)
                    if "PYTHONPATH" in env:
                        env["PYTHONPATH"] = f"{rf_abs}:{env['PYTHONPATH']}"
                    else:
                        env["PYTHONPATH"] = rf_abs
                    
                    subprocess.run(cmd, check=True, env=env, cwd=rf_abs)
                    # print(f"DEBUG: Skipping RFdiffusion. Input generated at {swapped_pdb_path}")
                except Exception as e:
                    print(f"    RFdiffusion failed: {e}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog_csv", default="../Local/Metal_Catalog/master_catalog.csv")
    parser.add_argument("--pdb_dir", default="../Local/Metal_PDBs")
    parser.add_argument("--out_dir", default="../Local/Metal_Predictions/Swaps")
    parser.add_argument("--rf_path", default="../Tools/RFdiffusion2")
    parser.add_argument("--target_metals", default="ZN,CU,NI,CO,MN,FE,MG,CA", help="Metals to swap TO")
    parser.add_argument("--num_designs", type=int, default=2)
    parser.add_argument("--dry_run", action="store_true")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output.")
    
    args = parser.parse_args()
    run_swap(args)

if __name__ == "__main__":
    main()
