import argparse
import os
import sys
import pandas as pd
import subprocess
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

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
    catalog = pd.read_csv(args.catalog_csv)
    target_metals = [m.strip().upper() for m in args.target_metals.split(',')]
    
    project_root = os.path.dirname(os.path.abspath(__file__))
    
    print(f"Loaded {len(catalog)} sites from catalog.")
    
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

        # Parse residues to find Target Chain and Designable Spots
        res_str = row['residues']
        chains = set()
        residue_ids = []
        
        for r in res_str.split(';'):
            # Parse A52(GLN) or 52(GLN) if no chain char?
            # PDB ID often has Chain.
            p1 = r.split('(')[0] # A52
            
            # Simple heuristic: First char is chain if alphabetic, else assume empty chain ' '?
            # Actually catalogue_sites.py outputs A52, so Chain is A.
            # If chain was ' ', it implies it might be just numbers?
            # Let's assume catalog is consistent: {Chain}{ResNum}.
            
            # Robust parsing:
            import re
            match = re.match(r"([A-Za-z0-9]+)(\d+)", p1)
            if match:
                chain_char = match.group(1)
                res_num = int(match.group(2))
                chains.add(chain_char)
                residue_ids.append((chain_char, res_num))
            else:
                # Fallback
                 pass
            
        if len(chains) > 1:
            print(f"  Skipping multi-chain site {motif_id}")
            continue
        if not chains:
             print(f"  Could not parse chain from {res_str}")
             continue
             
        target_chain = list(chains)[0]
        
        # Define Designable Residues (Binding Site +/- padding)
        design_residues = set()
        padding = 1
        for r_c, r_n in residue_ids:
            for k in range(r_n - padding, r_n + padding + 1):
                design_residues.add(k)
        
        chain_res_list = get_chain_residues(structure, target_chain)
        if not chain_res_list:
            print(f"  No residues found in chain {target_chain}")
            continue
            
        # Build Contig
        contig_parts = []
        start_segment = None
        prev_r = None
        
        design_block_len = 0 # To track length of contiguous design blocks
        
        for r in chain_res_list:
            is_designable = r in design_residues
            is_gap = (prev_r is not None) and (r != prev_r + 1)
            
            if is_gap:
                if start_segment is not None:
                     contig_parts.append(f"{target_chain}{start_segment}-{prev_r}")
                     start_segment = None
                if design_block_len > 0:
                    contig_parts.append(f"{design_block_len}-{design_block_len}")
                    design_block_len = 0
            
            if is_designable:
                if start_segment is not None:
                    contig_parts.append(f"{target_chain}{start_segment}-{prev_r}")
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
             contig_parts.append(f"{target_chain}{start_segment}-{prev_r}")
        if design_block_len > 0:
             contig_parts.append(f"{design_block_len}-{design_block_len}")

        contig_str = "/".join(contig_parts)
        
        # Metal Handling
        # We need to find the metal atom and ensure it is preserved as a separate chain (Z).
        # We will extract the metal atom, give it Chain Z, Res ID 1, and save it in input_swapped.pdb.
        
        atom_source_metal = None
        # Find the metal associated with this site.
        # Catalog has 'ion' (e.g. ZN).
        # We search atoms in structure.
        
        # Heuristic: Find ZN atom closest to the binding residues?
        # Or just pick the first ZN?
        # Let's pick the first one matching the source_ion for now, to support mass processing.
        # But we must ensure we only output ONE metal if we are swapping one site.
        
        # Actually, if we just want to swap metals, we can just replace the element in standard PDB?
        # But RFdiffusion needs to "see" the metal to design around it.
        # Best practice: Provide metal as a fixed input.
        
        # Extract protein chain + 1 Metal Atom -> input_swapped.pdb
        
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
             
        # Create Custom Structure for RFdiffusion
        # contain target_chain + Metal (remapped to Chain Z)
        
        # Metal Block
        metal_res_seq = 1
        metal_chain_id = 'Z'
        
        # Contig: ProteinParts/0 Z1-1
        contig = f"['{contig_str}/0 {metal_chain_id}{metal_res_seq}-{metal_res_seq}']"
        
        for metal in target_metals:
            out_dir = os.path.join(args.out_dir, metal, motif_id)
            os.makedirs(out_dir, exist_ok=True)
            rf_out = os.path.join(out_dir, "rfdiffusion")
            os.makedirs(rf_out, exist_ok=True)
            
            swapped_pdb_path = os.path.join(rf_out, "input_swapped.pdb")
            
            # Write PDB manually or via PDBIO
            with open(swapped_pdb_path, 'w') as f:
                # Write Protein Chain
                io = PDBIO()
                class ChainSelect(Select):
                    def accept_chain(self, chain):
                        return chain.id == target_chain
                
                io.set_structure(structure)
                # We can't easily pipe IO to file handle with other writes in BioPython < 1.78 sometimes?
                # Let's save protein to temp first.
                io.save(os.path.join(rf_out, "temp_prot.pdb"), ChainSelect())
                
                with open(os.path.join(rf_out, "temp_prot.pdb"), 'r') as temp:
                    f.write(temp.read())
                
                os.remove(os.path.join(rf_out, "temp_prot.pdb"))
                
                # Write Metal Atom as HETATM in Chain Z
                # Formatting PDB line manually
                # HETATM 1234 ZN   ZN Z   1      12.345  67.890  -5.432  1.00 20.00          ZN
                x, y, z = atom_source_metal.get_coord()
                pdb_line = f"HETATM{9999:5d} {metal:>4s} {metal:>3s} {metal_chain_id}{metal_res_seq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {metal:>2s}\n"
                f.write(pdb_line)
                
            # Run RFdiffusion
            rf_script = os.path.join(args.rf_path, "scripts/run_inference.py")
            cmd = [
                sys.executable, rf_script,
                f"inference.input_pdb={swapped_pdb_path}",
                f"inference.output_prefix={os.path.join(rf_out, 'design')}",
                f"inference.num_designs={args.num_designs}",
                f"contigmap.contigs={contig}"
            ]
            
            if args.dry_run:
                print(f"    Dry Run RFdiffusion ({metal}): {' '.join(cmd)}")
            else:
                try:
                    subprocess.run(cmd, check=True, cwd=project_root)
                except Exception as e:
                    print(f"    RFdiffusion failed: {e}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog_csv", default="../Local/Metal_Catalog/master_catalog.csv")
    parser.add_argument("--pdb_dir", default="../Local/Metal_PDBs")
    parser.add_argument("--out_dir", default="../Local/Metal_Predictions/Swaps")
    parser.add_argument("--rf_path", default="../Tools/RFdiffusion")
    parser.add_argument("--target_metals", default="CU,NI,CO,FE")
    parser.add_argument("--num_designs", type=int, default=2)
    parser.add_argument("--dry_run", action="store_true")
    
    args = parser.parse_args()
    run_swap(args)

if __name__ == "__main__":
    main()
