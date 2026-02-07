import argparse
import os
import sys
import numpy as np
import subprocess
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

# Import shared utilities
try:
    from utils_rfd2 import find_model_weights, add_ori_token
except ImportError:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from utils_rfd2 import find_model_weights, add_ori_token

def get_parser(file_path):
    if file_path.endswith(".cif"):
        return MMCIFParser(QUIET=True)
    return PDBParser(QUIET=True)

def get_center_of_mass(entity):
    atoms = list(entity.get_atoms())
    coords = np.array([a.get_coord() for a in atoms])
    return coords.mean(axis=0)

def translate_entity(entity, vector):
    for atom in entity.get_atoms():
        atom.set_coord(atom.get_coord() + vector)

def run_graft(args):
    project_root = os.path.dirname(os.path.abspath(__file__))
    
    # 1. Load Template
    parser = get_parser(args.template)
    template_struct = parser.get_structure("template", args.template)
    template_chain = template_struct[0][args.chain]
    
    # Find insertion point / replacement range
    try:
        if "-" in str(args.insert_at):
            start_str, end_str = args.insert_at.split("-")
            start_res = int(start_str)
            end_res = int(end_str)
        else:
            start_res = int(args.insert_at)
            end_res = int(args.insert_at) # Single residue replacement
            
        print(f"Grafting at {args.chain}:{start_res}-{end_res}")
        
        # Calculate target center (COM of replaced region)
        coords = []
        for i in range(start_res, end_res + 1):
            if i in template_chain:
                res = template_chain[i]
                if 'CA' in res:
                    coords.append(res['CA'].get_coord())
                else:
                    coords.append(get_center_of_mass(res))
        
        if coords:
            target_center = np.mean(coords, axis=0)
        else:
            print(f"Error: No residues found in range {start_res}-{end_res} in chain {args.chain}")
            return

    except Exception as e:
        print(f"Error parsing insert_at '{args.insert_at}': {e}")
        return

    # 2. Load Motif
    parser_m = get_parser(args.motif)
    motif_struct = parser_m.get_structure("motif", args.motif)
    motif_model = motif_struct[0] # Assume one model
    
    # Calculate Motif Center
    motif_center = get_center_of_mass(motif_model)
    
    # 3. Align Motif to Template (Translation only)
    translation = target_center - motif_center
    translate_entity(motif_model, translation)
    
    # 4. Merge Structures
    motif_chain_id = list(motif_model.child_dict.keys())[0]
    if motif_chain_id == args.chain:
        print(f"Warning: Motif chain {motif_chain_id} conflicts with Template chain {args.chain}. Renaming Motif to 'M'.")
        motif_model[motif_chain_id].id = 'M'
        motif_chain_id = 'M'
        
    merged_path = os.path.join(args.out_dir, "grafted_input.pdb")
    os.makedirs(args.out_dir, exist_ok=True)
    
    io = PDBIO()
    
    with open(merged_path, 'w') as out_f:
        # Write Template (Filtered for target chain)
        with open(args.template, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    chain_id = line[21]
                    if chain_id == args.chain:
                        out_f.write(line)
        out_f.write("TER\n")
        
        # Write Motif
        temp_motif = os.path.join(args.out_dir, "temp_motif_transformed.pdb")
        io.set_structure(motif_struct)
        io.save(temp_motif)
        
        with open(temp_motif, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    out_f.write(line)
        
        os.remove(temp_motif)
        
    print(f"Merged structure saved to {merged_path}")
    
    # Add ORI Token (Centered on Motif Metal if present, or Motif Center)
    # We use add_ori_token but need to know the ion.
    # We can infer ion from the Motif PDB or passed args?
    # run_graft didn't take ion arg.
    # Let's try to add ORI centered on ANY HETATM in Motif?
    # Or just use the 'translation' target point (which is where we put the motif).
    
    # Actually, add_ori_token takes ion_names. If we don't know, we can't search.
    # But we know where we put the motif: `target_center`.
    # Let's manually append ORI at `target_center`.
    
    with open(merged_path, 'a') as f:
        x, y, z = target_center
        ori_line = f"HETATM 9998  CA  ORI Z 999    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n"
        f.write(ori_line)
        f.flush()
        os.fsync(f.fileno())
    
    # Verify ORI presence
    with open(merged_path, 'r') as f:
         if "ORI" not in f.read():
             print("Error: ORI token failed to write to grafted_input.pdb")
             return

        
    # 5. Build Contig
    # [Template_Before, 10-20, Motif, 10-20, Template_After]
    
    # Get Template max res
    res_list = sorted([r.id[1] for r in template_chain if is_aa(r)])
    max_res = res_list[-1]
    
    contig_parts = []
    
    # Template Before
    if start_res > 1:
        contig_parts.append(f"{args.chain}1-{start_res - 1}")
    
    contig_parts.append("10-20") # Loop 1
    
    # Motif (Assume Chain M is fully preserved)
    # Motif (Assume Chain M is fully preserved but maybe discontinuous)
    m_res = sorted([r.id[1] for r in motif_struct[0][motif_chain_id] if is_aa(r)])
    
    if m_res:
         # Identify contiguous segments
         segments = []
         if not m_res:
             pass
         else:
             current_segment = [m_res[0]]
             for r in m_res[1:]:
                 if r == current_segment[-1] + 1:
                     current_segment.append(r)
                 else:
                     segments.append(current_segment)
                     current_segment = [r]
             segments.append(current_segment)
         
         for i, seg in enumerate(segments):
             start, end = seg[0], seg[-1]
             contig_parts.append(f"{motif_chain_id}{start}-{end}")
             
             # Add loop between segments if not last
             if i < len(segments) - 1:
                  contig_parts.append("6-15") # Flexible linker between fixed binding residues

    
    contig_parts.append("10-20") # Loop 2
    
    # Template After
    if end_res < max_res:
        contig_parts.append(f"{args.chain}{end_res + 1}-{max_res}")
        
    # RFdiffusion2 Format: ['part,part,part']
    contig = f"['{','.join(contig_parts)}']"
    
    # 6. Run RFdiffusion2
    rf_out = os.path.join(args.out_dir, "rfdiffusion")
    os.makedirs(rf_out, exist_ok=True)
    
    rf_script = os.path.join(os.path.abspath(args.rf_path), "rf_diffusion", "run_inference.py")
    
    ckpt_path = find_model_weights(args.rf_path)
    if ckpt_path:
        print(f"Using model weights: {os.path.basename(ckpt_path)}")
        
    output_prefix = os.path.join(rf_out, "design")
    
    # Detect Ligand/Metal in Motif
    ligand_name = None
    common_metals = {'ZN', 'CU', 'FE', 'NI', 'CO', 'MN', 'MG', 'CA', 'SF4', 'FES'}
    
    for model in motif_struct:
        for chain in model:
            for residue in chain:
                if residue.resname.strip().upper() in common_metals:
                    ligand_name = residue.resname.strip().upper()
                    break
            if ligand_name: break
        if ligand_name: break
    
    cmd = [
        sys.executable, rf_script,
        f"inference.input_pdb={os.path.abspath(merged_path)}",
        f"inference.output_prefix={os.path.abspath(output_prefix)}",
        f"inference.num_designs={args.num_designs}",
        f"contigmap.contigs={contig}"
    ]
    
    if ligand_name:
        print(f"Detected ligand in motif: {ligand_name}")
        cmd.append(f"inference.ligand='{ligand_name}'")
    
    if ckpt_path:
        cmd.append(f"inference.ckpt_path={ckpt_path}")
    
    if args.dry_run:
        print(f"Dry Run RFdiffusion: {' '.join(cmd)}")
    else:
        try:
            # PYTHONPATH
            env = os.environ.copy()
            rf_abs = os.path.abspath(args.rf_path)
            if "PYTHONPATH" in env:
                env["PYTHONPATH"] = f"{rf_abs}:{env['PYTHONPATH']}"
            else:
                env["PYTHONPATH"] = rf_abs
            
            subprocess.run(cmd, check=True, env=env, cwd=rf_abs)
        except Exception as e:
            print(f"RFdiffusion failed: {e}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--motif", required=True, help="Input Motif PDB")
    parser.add_argument("--template", default="../Local/Templates/human_VH3_IgG.pdb", help="Input Template PDB (Antibody). Default: VH3 IgG")
    parser.add_argument("--out_dir", required=True, help="Output directory")
    parser.add_argument("--chain", default="H", help="Template Chain to graft into")
    parser.add_argument("--insert_at", default="97-108", help="Residue number or range to replace (e.g. 97-108). Default: 97-108 (CDRH3)")
    
    parser.add_argument("--rf_path", default="../Tools/RFdiffusion2")
    parser.add_argument("--num_designs", type=int, default=2)
    parser.add_argument("--dry_run", action="store_true")
    
    args = parser.parse_args()
    run_graft(args)

if __name__ == "__main__":
    main()
