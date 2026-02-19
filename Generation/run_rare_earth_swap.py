
import argparse
import os
import sys
import numpy as np
import subprocess
import subprocess
import glob
import pandas as pd
from Bio.PDB import MMCIFParser, PDBIO, Select, PDBParser
from Bio.PDB.Polypeptide import is_aa

# Add current directory to path to import utils
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
try:
    from utils_rfd2 import find_model_weights
except ImportError:
    pass

def get_structure(file_path):
    if file_path.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
        return parser.get_structure(os.path.basename(file_path), file_path)
    else:
        parser = PDBParser(QUIET=True)
        return parser.get_structure(os.path.basename(file_path), file_path)

def get_neighbors(target_atom, atoms, radius=5.0):
    """Find residues with atoms within radius of target_atom."""
    neighbors = set()
    target_coord = target_atom.get_coord()
    for atom in atoms:
        if np.linalg.norm(atom.get_coord() - target_coord) < radius:
            parent = atom.get_parent()
            if is_aa(parent, standard=True):
                neighbors.add((parent.get_parent().id, parent.id[1])) # (Chain, ResNum)
    return sorted(list(neighbors))

def parse_fa_header(header):
    # >design_0-atomized-bb-False, id=1, T=0.1, seed=37, overall_confidence=0.4527, ligand_confidence=0.7263, seq_rec=0.4190
    data = {}
    parts = header.strip().lstrip('>').split(',')
    
    # First part is name
    data['name'] = parts[0]
    
    for part in parts[1:]:
        if '=' in part:
            k, v = part.strip().split('=')
            try:
                data[k] = float(v)
            except ValueError:
                data[k] = v
    return data

def generate_summary(root_dir, out_csv="rare_earth_summary.csv"):
    print(f"\n--- Generating Summary ---")
    results = []
    
    # Iterate over metal_site directories
    # Only look for folders matching METAL_SiteX pattern
    site_dirs = glob.glob(os.path.join(root_dir, "*_*"))
    
    print(f"  Found {len(site_dirs)} potential site directories in {root_dir}")

    for site_dir in site_dirs:
        if not os.path.isdir(site_dir):
            continue
            
        dirname = os.path.basename(site_dir)
        if '_' not in dirname:
            continue
            
        parts = dirname.split('_')
        if len(parts) < 2: 
             continue
        
        metal = parts[0]
        site = "_".join(parts[1:]) 
        
        # Find LigandMPNN outputs
        fa_files = glob.glob(os.path.join(site_dir, "ligandmpnn", "seqs", "*.fa"))
        
        if not fa_files:
            continue
            
        for fa_path in fa_files:
            basename = os.path.basename(fa_path).replace('.fa', '')
            
            # Read FASTA
            with open(fa_path, 'r') as f:
                lines = f.readlines()
                
            current_header = None
            for line in lines:
                line = line.strip()
                if line.startswith('>'):
                    current_header = parse_fa_header(line)
                else:
                    if current_header and 'id' in current_header:
                        row = {
                            'Metal': metal,
                            'Site': site,
                            'Design_Base': basename,
                            'Seq_ID': current_header.get('id'),
                            'Overall_Confidence': current_header.get('overall_confidence', np.nan),
                            'Ligand_Confidence': current_header.get('ligand_confidence', np.nan),
                            'Seq_Recovery': current_header.get('seq_rec', np.nan),
                            'Sequence': line
                        }
                        results.append(row)

    if results:
        df = pd.DataFrame(results)
        cols = ['Metal', 'Site', 'Design_Base', 'Seq_ID', 'Overall_Confidence', 'Ligand_Confidence', 'Seq_Recovery', 'Sequence']
        remaining = [c for c in df.columns if c not in cols]
        df = df[cols + remaining]
        
        out_path = os.path.join(root_dir, out_csv)
        df.to_csv(out_path, index=False)
        print(f"  Summary written to {out_path}")
    else:
        print("  No results found to summarize.")

def run_rare_earth_swap(args):
    # Setup Paths
    input_cif = os.path.abspath(args.input_cif)
    if not os.path.exists(input_cif):
        print(f"Error: Input file not found: {input_cif}")
        return

    # Load Structure
    structure = get_structure(input_cif)
    
    # Identify ND atoms
    nd_atoms = []
    all_atoms = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    all_atoms.append(atom)
                    if atom.element.upper() == args.binding_metal.upper():
                        nd_atoms.append(atom)
        break # Only first model

    if not nd_atoms:
        print(f"No {args.binding_metal} atoms found in structure.")
        return

    print(f"Found {len(nd_atoms)} {args.binding_metal} atoms.")
    
    # Target Metals
    target_metals = [m.strip() for m in args.metals.split(',')]
    
    # RFdiffusion Path
    rf_script = os.path.join(os.path.abspath(args.rf_path), "rf_diffusion", "run_inference.py")
    ckpt_path = find_model_weights(args.rf_path)
    
    # LigandMPNN Path
    # Assuming run_ligandmpnn.py logic or direct call. 
    # Use the existing wrapper or call the customized script?
    # I'll call the user's existing run_ligandmpnn.py if possible, or replicate the call.
    # Replicating the call is safer to avoid dependency on the wrapper's specific args.
    
    mpnn_script = os.path.join(os.path.abspath(args.mpnn_path), "run.py") 
    # Wait, looking at file list, there is a `run_ligandmpnn.py` in current dir.
    # I should use that or see how it calls the tool.
    # For now, I will use a direct subprocess call to the installed LigandMPNN if I can find it,
    # OR use the `run_ligandmpnn.py` wrapper. 
    # Better: Use the `run_ligandmpnn.py` wrapper to ensure consistency.
    local_mpnn_wrapper = os.path.abspath("run_ligandmpnn.py")

    for i, nd_atom in enumerate(nd_atoms):
        site_id = f"Site{i+1}"
        print(f"\n--- Processing {site_id} ---")
        
        # 1. Identify Binding Residues
        # These will be FIXED.
        binding_residues = get_neighbors(nd_atom, all_atoms, radius=args.radius)
        if not binding_residues:
            print(f"  No binding residues found for {site_id}")
            continue
            
        print(f"  Binding Residues: {binding_residues}")
        
        # 2. Define Contigs
        # Logic: Fixed residues are fixed. Gaps between them are loops to be redesigned.
        # But we need to include the whole protein structure as input?
        # If we provide input_pdb, we define contigs as the WHOLE chain, and use `fixed_residues` or `contig` logic?
        # RFdiffusion Partial Diffusion: 
        #   Provide input_pdb.
        #   contigs="[10-100]" (keep all? No, that just defines length).
        #   If we want to KEEP coordinates of binding sites and diffuse others, we use `contigmap.contigs` + `inference.input_pdb`.
        #   Actually, to FIX residues, we need to specify them in the contig string as fixed segments 
        #   OR rely on the fact that if we provide the PDB, we can specify which parts to keep.
        #   
        #   Correct approach for "Partial Diffusion" / "Inpainting":
        #   contigmap.contigs=["A1-19", "10-10", "A30-40"] 
        #      -> Keeps A1-19, Generates 10 residues, Keeps A30-40.
        #   
        #   The user wants to "fill in the loop domains".
        #   I will assume this means: KEEP the binding residues. REGENERATE everything else (or at least the loops between them).
        #   
        #   To keep it simple and robust:
        #   I'll map the chain. 
        #   Identify chunks of "Fixed" (Binding) and "Variable" (Non-binding).
        #   Construct contig string: "A1-19/A20-20/A21-29..."
        #   But purely fixed residues might detach if loops are generated from scratch.
        #   
        #   Better Strategy:
        #   Keep the binding residues fixed.
        #   Keep the REST of the protein fixed?
        #   If I simply "fill in loop domains around the new ions", maybe I should only regenerate the immediate loops connecting the binding sites?
        #   
        #   Let's check the residue numbers. If they are close (e.g. 20, 22, 24), the loops are tiny (1 residue).
        #   If they are far, the loop is long.
        #   
        #   I will treat binding residues as FIXED (source coordinates).
        #   I will treat everything else as FLEXIBLE (redesign).
        #   However, completely redesigning the Whole scaffold might be too much.
        #   
        #   Alternative: Keep everything fixed EXCEPT +/- N residues around binding sites?
        #   No, usually you keep the motif and redesign the scaffold.
        #   
        #   Let's go with: Keep Binding Residues FIXED. Keep the rest of the chain as Template (Noise? or Fixed?).
        #   "Fill in the loop domains" -> Redesign loops.
        #   
        #   I will parse the sequence from N to C.
        #   If residue is Binding: Output "A<res>-<res>".
        #   If residue is Non-Binding:
        #       Check if user wants it fixed or generated.
        #       "Fill in" implies generation.
        #       I will generate a length match for the gap.
        #       So "N-N" (length) instead of "A<start>-<end>" (sequence/structure).
        #       
        #   So:
        #   binding_res: {20, 22, 24}
        #   Chain A (1-100).
        #   Contig: "19-19" (Gen 1-19) + "A20-20" (Fix 20) + "1-1" (Gen 21) + "A22-22" (Fix 22) ...
        #   
        #   This effectively destroys the original scaffold and rebuilds it around the binding sites.
        #   This seems essentially like "Motif Scaffolding".
        #   
        #   Is that what the user wants? "fill in the loop domains around the new ions".
        #   Yes, implies the loops adapt to the new ion.
        
        # Get Chain IDs
        chains = sorted(list(set([r[0] for r in binding_residues])))
        
        contig_specs = []
        
        # Build contig string for the FIRST chain involved (Assuming single chain for now for simplicity, or handle multiples)
        # 8FNS is likely a single chain or homomer.
        # If multiple chains, we join them with "/0 " (break).
        
        chain_sequences = {} # ChainID -> List of (ResNum, is_binding)
        
        # reconstruct chain sequence from PDB
        for model in structure:
            for chain in model:
                if chain.id not in chains: continue
                
                res_list = []
                for r in chain:
                    if is_aa(r, standard=True):
                        full_id = (chain.id, r.id[1])
                        is_bind = full_id in binding_residues
                        res_list.append( (r.id[1], is_bind) )
                
                chain_sequences[chain.id] = sorted(res_list, key=lambda x: x[0])
        
        contig_parts = []
        for chain_id in chains:
            seq = chain_sequences.get(chain_id, [])
            if not seq: continue
            
            # Group into segments
            current_type = None # 'fixed' or 'gen'
            current_start = None
            current_end = None
            current_len = 0
            
            for r_num, is_bind in seq:
                type_ = 'fixed' if is_bind else 'gen'
                
                if type_ != current_type:
                    # Flush previous
                    if current_type == 'fixed':
                         # A<start>-<end>
                         contig_parts.append(f"{chain_id}{current_start}-{current_end}")
                    elif current_type == 'gen':
                         # <len>-<len>
                         contig_parts.append(f"{current_len}-{current_len}")
                    
                    current_type = type_
                    current_start = r_num
                    current_len = 0
                
                current_end = r_num
                current_len += 1
            
            # Flush last
            if current_type == 'fixed':
                 contig_parts.append(f"{chain_id}{current_start}-{current_end}")
            elif current_type == 'gen':
                 contig_parts.append(f"{current_len}-{current_len}")
            
            contig_parts.append("0") # Chain break
            
        if contig_parts and contig_parts[-1] == "0":
            contig_parts.pop()
            
        contig_str = "[" + "/".join(contig_parts) + "]" # RFd2 uses list-like string? No, just keys.
        # Actually RFdiffusion expects: '["contig1", "contig2"]' or 'contig1/contig2'.
        # Actually standard CLI: 'contigmap.contigs=[...]'
        # My utils_rfd2 might have helpers? No, I'll just format it manually.
        # The list format is safest: '["10-10", "A10-15"]'
        
        # Reformulate contig_parts to be a list of strings
        # Currently they are raw strings.
        final_contig_list_str = str(contig_parts).replace("'", '"').replace("0", "") # 0 is separate handling?
        # Actually if I have multiple chains, I need multiple entries in the list?
        # Or one entry with breaks?
        # RFd2: contigmap.contigs=['10-10/A5-10/10-10']  (Single quotes inside list)
        # So I will join with "/" and put in a single list item?
        # Or multiple items for multiple chains?
        # Usually one string with / is enough for multi-chain output if doing binder design?
        # I'll stick to a single string with steps.
        
        final_contig_str = ",".join(contig_parts).replace(",,", ",") # Simplify
        
        # Final Format for CLI: "['spec']"
        cli_contig = f"['{final_contig_str}']"
        
        
        # 3. Process Per Metal
        for metal in target_metals:
            print(f"  > Swapping to {metal}")
            
            out_subdir = os.path.join(args.out_dir, "RareEarthSwaps", f"{metal}_{site_id}")
            rf_out_dir = os.path.join(out_subdir, "rfdiffusion")
            os.makedirs(rf_out_dir, exist_ok=True)
            
            # Check if designs already exist
            if not args.overwrite and glob.glob(os.path.join(rf_out_dir, "design*.pdb")):
                print(f"    Skipping {site_id} -> {metal} (already exists)")
                continue
            
            # Create Input PDB (Swapped Metal)
            nd_pos = nd_atom.get_coord()
            
            input_pdb = os.path.join(rf_out_dir, "input_swapped.pdb")
            
            with open(input_pdb, 'w') as f:
                io = PDBIO()
                class ResSelect(Select):
                    def accept_residue(self, residue):
                         return is_aa(residue, standard=True)
                
                io.set_structure(structure)
                io.save(os.path.join(rf_out_dir, "temp.pdb"), ResSelect())
                
                # Copy Protein
                with open(os.path.join(rf_out_dir, "temp.pdb"), 'r') as temp:
                    for line in temp:
                         if line.startswith("ATOM"):
                             f.write(line)
                os.remove(os.path.join(rf_out_dir, "temp.pdb"))
                
                # Write Metal
                # Chain Z, Res 1
                f.write(f"HETATM 9999 {metal:>4} {metal:>3} Z   1    {nd_pos[0]:8.3f}{nd_pos[1]:8.3f}{nd_pos[2]:8.3f}  1.00 20.00          {metal:>2}\n")
                
                # Write ORI Token (Centered on Metal)
                f.write(f"HETATM 9998  CA  ORI Z 999    {nd_pos[0]:8.3f}{nd_pos[1]:8.3f}{nd_pos[2]:8.3f}  1.00 20.00           C  \n")
            
            
            # Run RFdiffusion
            # We need to specify the Ligand.
            
            cmd = [
                sys.executable, rf_script,
                f"inference.input_pdb={os.path.abspath(input_pdb)}",
                f"inference.output_prefix={os.path.abspath(os.path.join(rf_out_dir, 'design'))}",
                f"inference.num_designs={args.num_designs}",
                f"contigmap.contigs={cli_contig}",
                f"inference.ligand={metal}"
            ]
            
            if ckpt_path:
                 cmd.append(f"inference.ckpt_path={ckpt_path}")

            print(f"    Running RFdiffusion...")
            print(" ".join(cmd))
            
            try:
                env = os.environ.copy()
                rf_root = os.path.abspath(args.rf_path)
                env["PYTHONPATH"] = f"{rf_root}:{env.get('PYTHONPATH', '')}"
                
                # Run with env and cwd so RFdiffusion finds modules. Output to console.
                subprocess.run(cmd, check=True, env=env, cwd=rf_root)
                print("    RFdiffusion Complete.")
            except subprocess.CalledProcessError as e:
                print(f"    RFdiffusion Failed.")
                if e.stderr:
                    print(e.stderr.decode())
                continue
            
            # Run LigandMPNN
            print("    Running LigandMPNN...")
            
            # Check for designs
            designs = glob.glob(os.path.join(rf_out_dir, "design*.pdb"))
            if not designs:
                print("    No designs found.")
                continue
                
            mpnn_cmd = [
                sys.executable, local_mpnn_wrapper,
                "--pred_dir", out_subdir, # Pass the parent site dir, wrapper will find rfdiffusion/*.pdb
                "--lmpnn_path", args.mpnn_path
            ]
            
            try:
                subprocess.run(mpnn_cmd, check=True, stdout=subprocess.DEVNULL)
                print("    LigandMPNN Complete.")
            except subprocess.CalledProcessError as e:
                print(f"    LigandMPNN Failed.")
    
    # Generate Summary after all processing
    generate_summary(os.path.join(args.out_dir, "RareEarthSwaps"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_cif", default="../Local/Metal_PDBs/LA/8FNS.cif")
    parser.add_argument("--out_dir", default="../Local/Metal_Predictions")
    parser.add_argument("--rf_path", default="../Tools/RFdiffusion2")
    parser.add_argument("--mpnn_path", default="../Tools/LigandMPNN") # Adjust if needed
    parser.add_argument("--metals", default="LA,ND,SM,GD,TB,Y")
    parser.add_argument("--binding_metal", default="ND", help="Metal element in input PDB to define binding sites")
    parser.add_argument("--radius", type=float, default=6.0, help="Radius to define binding site")
    parser.add_argument("--num_designs", type=int, default=1)
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output.")
    
    args = parser.parse_args()
    run_rare_earth_swap(args)

if __name__ == "__main__":
    main()
