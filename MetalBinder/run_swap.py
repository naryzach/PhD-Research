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

def residues_to_ranges(residue_list, chain_id):
    """Convert list of residue numbers to RFdiffusion ranges."""
    if not residue_list:
        return []
    ranges = []
    start = residue_list[0]; prev = residue_list[0]
    for r in residue_list[1:]:
        if r != prev + 1:
            ranges.append(f"{chain_id}{start}-{prev}")
            start = r
        prev = r
    ranges.append(f"{chain_id}{start}-{prev}")
    return ranges

def run_swap(args):
    catalog = pd.read_csv(args.catalog_csv)
    target_metals = [m.strip().upper() for m in args.target_metals.split(',')]
    
    project_root = os.path.dirname(os.path.abspath(__file__))
    
    print(f"Loaded {len(catalog)} sites from catalog.")
    
    for idx, row in catalog.iterrows():
        source_pdb = row['source_pdb']
        source_ion = row['ion']
        motif_id = row['motif_id']
        path_in_catalog = row['path']
        
        # Locate Source PDB
        # Assumes PDBs are in default location: Local/Metal_PDBs/{Ion}/{ID}.cif or .pdb
        # The catalog doesn't store full source path, only source PDB ID + Ion.
        # We need to construct it.
        
        # Try both extensions
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

        # Parse residues from catalog string: "A52(GLN);A60(GLY)"
        # We need to identify the protein chain involved to build the contig.
        # Catalog only has 'residues' string.
        # Extract Chain IDs.
        res_str = row['residues']
        chains = set()
        residue_ids = []
        
        for r in res_str.split(';'):
            # Parse A52(GLN) -> Chain A, Res 52
            # Simple regex-like parsing
            # A52(GLN)
            p1 = r.split('(')[0] # A52
            chain_char = p1[0]
            res_num = int(p1[1:])
            chains.add(chain_char)
            residue_ids.append((chain_char, res_num))
            
        if len(chains) > 1:
            print(f"  Skipping multi-chain site {motif_id} (complex to inpaint)")
            continue
            
        target_chain = list(chains)[0]
        
        # Build Gap-Aware "Keep" regions
        # We want to redesign LOOPs around the site.
        # Actually RFdiffusion "Inpainting" usually works by masking the region to design.
        # Here we want to keep the whole protein but redesign the neighborhood of the metal?
        # Or just swap the metal?
        # Metal swap usually implies geometry optimization (RFdiffusion) around the metal.
        # So we "mask" the binding site residues to be redesigned?
        # The prompt says: "Inpaint loops".
        # Let's design 5 residues around the binding residues.
        
        # Identify residues to redesign (exclude from 'Fixed' contig)
        # We define a "Binding Pocket" range.
        
        # Actually, let's keep it simple: Fixed Scaffold. Just redesign the immediate residues?
        # RFdiffusion is better at redesigning loops.
        # Let's simply fix the whole protein EXCEPT the binding residues +/- padding?
        
        # Logic: 
        # 1. Identify all residues in chain.
        # 2. Identify "Designable" residues (The binding residues + padding).
        # 3. Everything else is "Fixed".
        
        design_residues = set()
        padding = 1
        for r_c, r_n in residue_ids:
            for k in range(r_n - padding, r_n + padding + 1):
                design_residues.add(k)
                
        chain_res_list = get_chain_residues(structure, target_chain)
        if not chain_res_list:
            continue
            
        # Build Fixed Ranges
        fixed_ranges = []
        start = chain_res_list[0]
        prev = chain_res_list[0]
        current_segment = []
        
        # Iterate and build ranges that are NOT in design_residues
        segments = []
        
        # Using boolean mask approach is easier?
        # Iterate sorted chain_res_list
        # If r is NOT in design_residues -> add to current fixed block
        # If r IS in design_residues -> close current fixed block, start new one after.
        
        # But we also need to respect gaps in the PDB itself (missing density).
        
        # RFdiffusion input: [Fixed/10-10/Fixed...] or [Fixed/10/Fixed] (length)
        # If we want to redesign existing residues, we just OMIT them from the fixed contig and specifying length.
        # E.g. A1-49/3/A53-100 (Redesigns 50-52 as 3 residues).
        
        contig_parts = []
        
        start_segment = None
        prev_r = None
        
        for r in chain_res_list:
            is_designable = r in design_residues
            
            # Check for gap in PDB
            is_gap = (prev_r is not None) and (r != prev_r + 1)
            
            if is_gap:
                # Close any open segment
                if start_segment is not None:
                     contig_parts.append(f"{target_chain}{start_segment}-{prev_r}")
                     start_segment = None
            
            if is_designable:
                # If we were tracking a fixed segment, close it
                if start_segment is not None:
                    contig_parts.append(f"{target_chain}{start_segment}-{prev_r}")
                    start_segment = None
                
                # Add design step (1 residue)
                # We add '1-1' to indicate 1 residue to be built? Or multiple?
                # If we want to preserve length, we add '1-1'.
                contig_parts.append("1-1") 
            else:
                # Fixed residue
                if start_segment is None:
                    start_segment = r
                prev_r = r
                
        # Close last segment
        if start_segment is not None:
             contig_parts.append(f"{target_chain}{start_segment}-{prev_r}")
             
        # Compact "1-1/1-1/1-1" to "3-3" logic? 
        # RFdiffusion is fine with lists.
        # But "1-1/1-1" is creating multiple segments.
        # Ideally we group designable blocks.
        
        # Simplify: Join contig_parts with '/'
        # Clean up repeated design blocks? 
        # For now, it works.
        
        contig_str = "/".join(contig_parts)
        contig = f"['{contig_str} Z0-0']" # Append Metal Chain Z (we will create it)
        
        # Run RFdiffusion for each target metal
        for metal in target_metals:
            out_dir = os.path.join(args.out_dir, metal, motif_id)
            os.makedirs(out_dir, exist_ok=True)
            rf_out = os.path.join(out_dir, "rfdiffusion")
            os.makedirs(rf_out, exist_ok=True)
            
            # Prepare Input PDB with Changed Metal
            # We need to create a temporary input PDB where the metal atom is swapped.
            # Using Site Info to find the metal atom.
            # Actually, extract the site metal? The source PDB has a metal.
            # We assume the metal found in catalog is the one to swap.
            # We need to find the specific metal atom corresponding to 'site'.
            # Catalog doesn't store atom serial number.
            # But we can find the metal near the residues.
            
            # For simplicity in mass processing:
            # Just take the first metal of type 'source_ion' in the structure?
            # Or use the catalog logic again to identify the metal?
            # Re-finding is safer.
            
            # Skip actual swap file generation for dry run speed if needed, but logic is fast.
            # Let's assume we output to `swapped_input.pdb`.
            
            swapped_pdb_path = os.path.join(rf_out, "input_swapped.pdb")
            
            io = PDBIO()
            # Swap metal element in structure object (in memory)
            # Find closest metal to residues
            # ... (Simplification: We assume single metal or we just swap ALL metals of that type?)
            
            # Let's swap ALL metals of source_ion to target_metal for the input PDB.
            # This is easiest.
            
            for atom in structure.get_atoms():
                if atom.element.upper() == source_ion:
                    atom.element = metal
                    atom.name = metal
            
            io.set_structure(structure)
            io.save(swapped_pdb_path)
            
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
