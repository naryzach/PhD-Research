import argparse
import os
import sys
import copy
import numpy as np
import subprocess
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select, Superimposer, Vector, rotaxis
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB import NeighborSearch

# Import shared utilities
try:
    from utils_rfd2 import find_model_weights, add_ori_token
except ImportError:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from utils_rfd2 import find_model_weights, add_ori_token

from Bio.PDB import PDBList

def rotate_around_axis(entity, axis_start, axis_end, angle_deg):
    """
    Rotates entity around an axis defined by two points.
    angle_deg: Angle in degrees
    axis_start, axis_end: Bio.PDB.Vector or np.array coordinates
    """
    angle_rad = np.radians(angle_deg)
    
    # Ensure vectors
    try:
        if not isinstance(axis_start, Vector): axis_start = Vector(axis_start)
        if not isinstance(axis_end, Vector): axis_end = Vector(axis_end)
    except:
        axis_start = Vector(list(axis_start))
        axis_end = Vector(list(axis_end))
    
    axis = axis_end - axis_start
    axis.normalize()
    
    # Rotation matrix (returns 3x3)
    rot = rotaxis(angle_rad, axis)
    
    # Apply to atoms: translate -> rotate -> translate_back
    start_arr = axis_start.get_array()
    
    for atom in entity.get_atoms():
        c = atom.get_coord()
        c_centered = c - start_arr
        c_rotated = np.dot(rot, c_centered) 
        c_final = c_rotated + start_arr
        atom.set_coord(c_final)

def count_clashes(entity_mobile, entity_static, threshold=2.5, exclude_chain_id=None):
    """
    Counts clashes between mobile entity and static entity.
    Default threshold 2.5A for sidechains roughly.
    """
    # Build NS on static
    atoms_static = []
    for a in entity_static.get_atoms():
         # Exclude if chain ID matches
         chain_id = a.get_parent().get_parent().id
         if exclude_chain_id and chain_id == exclude_chain_id:
             continue
         atoms_static.append(a)
    
    if not atoms_static: return 0
    
    ns = NeighborSearch(atoms_static)
    clashes = 0
    for atom in entity_mobile.get_atoms():
         if ns.search(atom.get_coord(), threshold, level='A'):
             clashes += 1
    return clashes

def ensure_template_exists(template_path, default_pdb_id="1HZH"):
    """
    Checks if template exists. If not, downloads default PDB (1HZH),
    extracts Chain H, and saves it to template_path.
    """
    if os.path.exists(template_path):
        return template_path, False
        
    print(f"Template {template_path} not found. Downloading default {default_pdb_id}...")
    
    # Ensure directory exists
    pdir = os.path.dirname(os.path.abspath(template_path))
    os.makedirs(pdir, exist_ok=True)
    
    # Download using urllib (HTTP) because PDBList (FTP) might be blocked
    import urllib.request
    url = f"https://files.rcsb.org/download/{default_pdb_id.upper()}.pdb"
    ent_file = os.path.join(pdir, f"{default_pdb_id.lower()}.pdb")
    
    try:
        urllib.request.urlretrieve(url, ent_file)
        print(f"Downloaded {url} to {ent_file}")
    except Exception as e:
        if os.path.exists(ent_file): os.remove(ent_file)
        raise FileNotFoundError(f"Failed to download {default_pdb_id} from {url}: {e}")

    if not os.path.exists(ent_file):
        raise FileNotFoundError(f"Failed to download {default_pdb_id} to {ent_file}")

    # Extract Chain H and save to template_path
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(default_pdb_id, ent_file)
    
    io = PDBIO()
    class ChainHSelect(Select):
        def accept_chain(self, chain):
            return chain.id == 'H'
            
    io.set_structure(structure)
    io.save(template_path, ChainHSelect())
    
    # Clean up raw download if different
    if os.path.abspath(ent_file) != os.path.abspath(template_path):
        os.remove(ent_file)
        
    print(f"Saved default template to {template_path}")
    return template_path, True

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

def detect_cdrh3(structure, chain_id):
    """
    Heuristic detection of CDRH3 using conserved Cys (~92) and Trp (~103).
    Returns (start, end) residue numbers (exclusive of anchors).
    """
    try:
        chain = structure[0][chain_id]
        sorted_res = sorted([r for r in chain if is_aa(r)], key=lambda x: x.id[1])
        
        # Find last Cys before 105 (allow for numbering shifts)
        cys_res = None
        for r in sorted_res:
            if r.resname == 'CYS' and r.id[1] < 105:
                 cys_res = r
                 
        # Find first Trp after Cys
        trp_res = None
        if cys_res:
            cys_index = sorted_res.index(cys_res)
            for r in sorted_res[cys_index+1:]:
                 if r.resname == 'TRP':
                      trp_res = r
                      break
                      
        if cys_res and trp_res:
             print(f"Detected CDRH3 anchors: CYS {cys_res.id[1]} - TRP {trp_res.id[1]}")
             return cys_res.id[1] + 1, trp_res.id[1] - 1
             
    except Exception as e:
        print(f"Warning: CDRH3 detection failed: {e}")
        
    return None

def run_graft(args):
    project_root = os.path.dirname(os.path.abspath(__file__))
    
    # 1. Load Template
    args.template, downloaded = ensure_template_exists(args.template)
    parser = get_parser(args.template)
    template_struct = parser.get_structure("template", args.template)
    
    # Check/Detect range
    if downloaded:
        detected_range = detect_cdrh3(template_struct, args.chain)
        if detected_range:
            start, end = detected_range
            print(f"Automatically detected CDRH3 range for default template: {start}-{end}")
            args.insert_at = f"{start}-{end}"
    
    template_chain = next(iter(template_struct))[args.chain]
    
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
    
    # 3. Align Motif to Template
    # Try Superimposer (Stem alignment) first
    aligned = False
    try:
        # Anchors: start_res - 1, end_res + 1
        anchor_ids = [start_res - 1, end_res + 1]
        template_anchors = []
        for rid in anchor_ids:
            if rid in template_chain and 'CA' in template_chain[rid]:
                template_anchors.append(template_chain[rid]['CA'])
        
        # Stems: First and Last CA of motif
        motif_residues = sorted([r for r in motif_model.get_residues() if is_aa(r)], key=lambda x: x.id[1])
        motif_stems = []
        if len(motif_residues) >= 2:
            if 'CA' in motif_residues[0]: motif_stems.append(motif_residues[0]['CA'])
            if 'CA' in motif_residues[-1]: motif_stems.append(motif_residues[-1]['CA'])
            
        if len(template_anchors) == 2 and len(motif_stems) == 2:
            print(f"Aligning Motif stems to Template anchors {anchor_ids}...")
            sup = Superimposer()
            sup.set_atoms(template_anchors, motif_stems)
            sup.apply(motif_model.get_atoms())
            print(f"Alignment RMSD: {sup.rms:.3f}")
            aligned = True
            
            # Recalculate target center based on aligned motif
            motif_center = get_center_of_mass(motif_model)
            target_center = motif_center 

            # --- Rotation Optimization ---
            print("Optimizing rotation around stem axis to minimize steric clashes...")
            
            # Prepare Template for clash check (excluding replaced region)
            # We clone the template structure to modify it for clash checking
            temp_struct_clash = template_struct.copy()
            temp_model_clash = temp_struct_clash[0]
            temp_chain_clash = temp_model_clash[args.chain]
            
            # Remove residues in the insertion range from the clash check template
            # (they will be removed anyway, so they shouldn't count as clashes)
            ids_to_remove_clash = [r.id for r in temp_chain_clash if is_aa(r) and start_res <= r.id[1] <= end_res]
            for rid in ids_to_remove_clash:
                temp_chain_clash.detach_child(rid)
            
            best_angle = 0
            min_clashes = float('inf')
            max_dist = -float('inf')
            
            # Calculate Scaffold Center (Approximate, using clash model)
            scaffold_coords = [a.get_coord() for a in temp_model_clash.get_atoms()]
            if scaffold_coords:
                 scaffold_center = np.mean(scaffold_coords, axis=0)
            else:
                 scaffold_center = np.array([0,0,0])

            # Calculate Graft Base (Midpoint of anchors)
            # Use original anchor coordinates (not modified by rotation)
            # template_anchors contains the CA atoms of start-1 and end+1
            anchor_coords = [a.get_coord() for a in template_anchors]
            if len(anchor_coords) == 2:
                graft_base = np.mean(anchor_coords, axis=0)
            else:
                graft_base = scaffold_center # Fallback
            
            # Outward Vector: Core -> Graft Base
            outward_vec = graft_base - scaffold_center
            norm = np.linalg.norm(outward_vec)
            if norm > 0.001: outward_vec /= norm
            
            # Initial state score
            clashes = count_clashes(motif_model, temp_model_clash, threshold=3.5)
            motif_now_center = get_center_of_mass(motif_model)
            motif_vec = motif_now_center - graft_base
            proj = np.dot(motif_vec, outward_vec)
            
            print(f"  Angle 0 deg: {clashes} clashes, Projection: {proj:.2f}")
            
            # Record scores: (angle, clashes, projection)
            scores = [(0, clashes, proj)]

            # Optimization Loop
            step_deg = 10
            last_angle = 0
            
            # Define axis atoms 
            stem_start_vec = motif_stems[0].get_vector()
            stem_end_vec = motif_stems[-1].get_vector()

            for angle in range(step_deg, 360, step_deg):
                # Rotate by step
                rotate_around_axis(motif_model, stem_start_vec, stem_end_vec, step_deg)
                
                clashes = count_clashes(motif_model, temp_model_clash, threshold=3.5)
                
                # Calculate Projection
                motif_now_center = get_center_of_mass(motif_model)
                motif_vec = motif_now_center - graft_base
                proj = np.dot(motif_vec, outward_vec)
                
                scores.append((angle, clashes, proj))
                last_angle = angle
            
            # Selection Logic:
            # 1. Filter for minimal clashes
            min_c = min(s[1] for s in scores)
            
            # 2. Filter candidates with min_c
            candidates = [s for s in scores if s[1] == min_c]
            
            # 3. Select max projection among candidates
            best_score = max(candidates, key=lambda x: x[2])
            
            best_angle = best_score[0]
            best_clashes = best_score[1]
            best_proj = best_score[2]
            
            print(f"Best rotation found: {best_angle} degrees with {best_clashes} clashes and outward projection {best_proj:.2f}.")
            
            # Apply Best Rotation
            # Current state is at 'last_angle' (e.g. 350)
            correction_angle = best_angle - last_angle
            
            if correction_angle != 0:
                print(f"Applying correction rotation of {correction_angle} degrees...")
                rotate_around_axis(motif_model, stem_start_vec, stem_end_vec, correction_angle)
            
            # Re-calculate center
            motif_center = get_center_of_mass(motif_model)
            target_center = motif_center
            # --- End Optimization ---
        else:
            print("Insufficient atoms for stem alignment (need 2 anchors + 2 stems).")
            
    except Exception as e:
        print(f"Alignment failed: {e}")
        
    if not aligned:
        print("Falling back to Center-of-Mass translation...")
        motif_center = get_center_of_mass(motif_model)
        translation = target_center - motif_center
        translate_entity(motif_model, translation)
    
    # 4. Merge Structures
    
    # Create merged structure first so we have the base
    merged_structure = template_struct.copy()
    merged_model = merged_structure[0]
    merged_chain = merged_model[args.chain]
    
    # Remove replaced residues from template chain in merged model
    ids_to_remove = []
    for res in merged_chain:
        if is_aa(res) and start_res <= res.id[1] <= end_res:
             ids_to_remove.append(res.id)
    
    for rid in ids_to_remove:
        merged_chain.detach_child(rid)
    print(f"Removed {len(ids_to_remove)} residues ({start_res}-{end_res}) from template chain {args.chain}")

    # Identify Metal in Motif (for clash cleaning later)
    metal_atoms_in_motif = []
    common_metals = {'ZN', 'CU', 'FE', 'NI', 'CO', 'MN', 'MG', 'CA', 'SF4', 'FES', 'CD'}
    
    # Merge ALL chains from Motif into Template
    # Logic: Rename chains if conflict
    existing_chains = set([c.id for c in merged_model])
    
    # Re-do merge loop robustly
    
    # Re-do merge loop robustly
    motif_chains = list(motif_model.get_chains()) # Get list
    for m_chain in motif_chains:
        # Save metal atoms object refs (they move with chain)
        for res in m_chain:
            if res.resname.strip().upper() in common_metals:
                 for atom in res:
                     metal_atoms_in_motif.append(atom)
        
        old_id = m_chain.id
        motif_model.detach_child(old_id)
        
        new_id = old_id
        if new_id in existing_chains:
             # Find new ID
             import string
             for char in "MNOPQRSTUVWXYZABCDEFGHIJKL":
                 if char not in existing_chains:
                     new_id = char
                     break
             print(f"Renaming Motif chain {old_id} to {new_id}")
             m_chain.id = new_id
             
        merged_model.add(m_chain)
        # Update motif_chain_id for contig generation if it was the main chain
        # Heuristic: if this chain has the most residues, it's likely the motif protein chain
        # But `motif_chain_id` used later is for the protein part.
    
    # Re-identify the main motif protein chain for contig generation
    # Warning: this assumes one main protein chain in motif.
    # We need to find which chain ID corresponds to the protein part we aligned.
    # The `motif_stems` came from `motif_model`. The chain ID might have changed.
    # Let's find the chain that contains `motif_stems[0]`? 
    # Atoms conform to hierarchy. `motif_stems[0].get_parent().get_parent().id`
    
    motif_chain_id = 'M' # Default fallback
    if motif_stems:
         motif_chain_id = motif_stems[0].get_parent().get_parent().id
    else:
         # Fallback search
         for c in merged_model:
             if c.id not in [args.chain]: # not template
                 motif_chain_id = c.id
                 break

    merged_path = os.path.join(args.out_dir, "grafted_input.pdb")
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Clash Cleaning
    if metal_atoms_in_motif:
        print(f"Checking for clashes with {len(metal_atoms_in_motif)} metal atoms...")
        ns = NeighborSearch(list(merged_model.get_atoms()))
        clash_radius = 3.5
        
        atoms_to_remove = []
        residues_to_clean = set()
        
        for metal_atom in metal_atoms_in_motif:
            close_atoms = ns.search(metal_atom.get_coord(), clash_radius, level='A')
            for atom in close_atoms:
                # Check if atom belongs to Template Chain
                parent_chain = atom.get_parent().get_parent()
                if parent_chain.id == args.chain:
                    residues_to_clean.add(atom.get_parent())
        
        if residues_to_clean:
            print(f"Cleaning {len(residues_to_clean)} clashing residues from template...")
            for res in residues_to_clean:
                 # Check if we should remove (don't remove anchors?)
                 # For now, remove everything close to metal to clear pocket
                 try:
                     merged_model[args.chain].detach_child(res.id)
                     print(f"  Removed clashing residue {args.chain} {res.id[1]}")
                 except:
                     pass
    
    io = PDBIO()
    io.set_structure(merged_structure)
    io.save(merged_path)
        
    print(f"Merged structure saved to {merged_path}")
    
    # Add ORI Token (Centered on Metal if possible)
    # If we have metal atoms, use their center
    if metal_atoms_in_motif:
         coords = np.array([a.get_coord() for a in metal_atoms_in_motif])
         target_center = np.mean(coords, axis=0) # Update center for ORI
         
    with open(merged_path, 'a') as f:
        x, y, z = target_center
        ori_line = f"HETATM 9998  CA  ORI Z 999    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n"
        f.write(ori_line)
        f.flush()
        # os.fsync(f.fileno()) # Optional
    
    # Verify ORI presence
    with open(merged_path, 'r') as f:
         if "ORI" not in f.read():
             print("Error: ORI token failed to write to grafted_input.pdb")
             return

        
    # 5. Build Contig
    # [Template_Before, 10-20, Motif, 10-20, Template_After]
    
    # Helper for contiguous segments
    def get_contig_segs(residues):
        if not residues: return []
        segments = []
        current = [residues[0]]
        for r in residues[1:]:
            if r == current[-1] + 1:
                current.append(r)
            else:
                segments.append(current)
                current = [r]
        segments.append(current)
        return [(s[0], s[-1]) for s in segments]

    contig_parts = []
    
    # helper to get merged chain
    merged_model = next(iter(merged_structure))
    merged_chain = merged_model[args.chain]
    
    # Template Before (Handle Gaps)
    # Get all AA residues < start_res FROM MERGED CHAIN
    before_res = sorted([r.id[1] for r in merged_chain if is_aa(r) and r.id[1] < start_res])
    if before_res:
         segs = get_contig_segs(before_res)
         for start, end in segs:
             contig_parts.append(f"{args.chain}{start}-{end}")
    
    contig_parts.append("10-20") # Loop 1
    
    # Motif (Assume Chain M is fully preserved)
    # Need to find the motif chain in the merged model
    motif_chain = None
    motif_chain_id = None
    for c in merged_model:
        if c.id != args.chain:
            motif_chain = c
            motif_chain_id = c.id
            break
            
    if motif_chain:
        m_res = sorted([r.id[1] for r in motif_chain if is_aa(r)])
        if m_res:
             segs = get_contig_segs(m_res)
             for i, (start, end) in enumerate(segs):
                 contig_parts.append(f"{motif_chain_id}{start}-{end}")
                 if i < len(segs) - 1:
                      contig_parts.append("6-15") # Linker
    else:
        print("Error: Could not find motif chain in merged model for contig generation.")
        return

    contig_parts.append("10-20") # Loop 2
    
    # Template After (Handle Gaps)
    # Get all AA residues > end_res FROM MERGED CHAIN
    after_res = sorted([r.id[1] for r in merged_chain if is_aa(r) and r.id[1] > end_res])
    if after_res:
         segs = get_contig_segs(after_res)
         for start, end in segs:
             contig_parts.append(f"{args.chain}{start}-{end}")
        
    # RFdiffusion2 Format: ['part,part,part']
    contig = f"['{','.join(contig_parts)}']"
    
    # 6. Run RFdiffusion2
    rf_out = os.path.join(args.out_dir, "rfdiffusion")
    os.makedirs(rf_out, exist_ok=True)
    
    output_prefix = os.path.join(rf_out, "design")
    
    # Check for existing outputs
    import glob
    existing_designs = glob.glob(f"{output_prefix}*.pdb")
    if len(existing_designs) >= args.num_designs and not args.overwrite:
        print(f"Skipping RFdiffusion: Found {len(existing_designs)} existing designs in {rf_out}. Use --overwrite to re-run.")
        return

    rf_script = os.path.join(os.path.abspath(args.rf_path), "rf_diffusion", "run_inference.py")
    
    ckpt_path = find_model_weights(args.rf_path)
    if ckpt_path:
        print(f"Using model weights: {os.path.basename(ckpt_path)}")
        
    # Detect Ligand/Metal in Merged Structure
    ligand_name = None
    ligand_block = ""
    common_metals = {'ZN', 'CU', 'FE', 'NI', 'CO', 'MN', 'MG', 'CA', 'SF4', 'FES', 'CD'}
    
    # Scan merged_structure because motif atoms are moved there
    for model in merged_structure:
        for chain in model:
            for residue in chain:
                if residue.resname.strip().upper() in common_metals:
                    ligand_name = residue.resname.strip().upper()
                    break
            if ligand_name: break
        if ligand_name: break
            
    if ligand_name:
        # Extract ligand lines from merged_structure
        class LigandSelect(Select):
             def accept_residue(self, residue):
                 return residue.resname.strip().upper() == ligand_name
        
        io_lig = PDBIO()
        io_lig.set_structure(merged_structure)
        temp_lig_file = os.path.join(args.out_dir, "temp_ligand.pdb")
        io_lig.save(temp_lig_file, LigandSelect())
        
        with open(temp_lig_file, 'r') as f:
            lines = f.readlines()
            # Filter for HETATM/ATOM to be safe
            ligand_block = "".join([l for l in lines if l.startswith("HETATM") or l.startswith("ATOM")])
        
        if os.path.exists(temp_lig_file):
             os.remove(temp_lig_file)

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
        # In dry run, we can't append to non-existent output, but we can print the block
        if ligand_block:
            print("Ligand block to be appended:")
            print(ligand_block)
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
            
            # Post-processing: Append ligand to outputs
            if ligand_block:
                designs = glob.glob(f"{output_prefix}*.pdb")
                for d in designs:
                    # Check if ligand already appended to avoid duplication on re-runs (though we skip re-runs now)
                    with open(d, 'r') as f:
                        content = f.read()
                    if ligand_name not in content: # heuristic check
                         with open(d, 'a') as f:
                            f.write("\n" + ligand_block)
                         print(f"Appended ligand {ligand_name} to {os.path.basename(d)}")

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
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing designs")
    
    args = parser.parse_args()
    run_graft(args)

if __name__ == "__main__":
    main()
