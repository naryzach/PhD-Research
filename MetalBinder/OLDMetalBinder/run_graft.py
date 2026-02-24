import argparse
import os
import sys
import copy
import numpy as np
import subprocess
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select, Superimposer, Vector, rotaxis, Structure, Model, Chain
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
    if metal_atoms_in_motif and not args.no_clash_cleanup:
        print(f"Checking for clashes with {len(metal_atoms_in_motif)} metal atoms (Threshold: {args.clash_threshold} A)...")
        ns = NeighborSearch(list(merged_model.get_atoms()))
        clash_radius = args.clash_threshold
        
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
    elif args.no_clash_cleanup:
        print("Skipping clash cleaning (disabled by user).")
    
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

            
    # 5. Renumber Merged Structure (Crucial for RFdiffusion)
    # RFdiffusion works best with 1-indexed, contiguous residues.
    # We will renumber the MERGED model (template + graft) to 1..N.
    # This avoids insertion code issues and duplicate IDs.
    
    print("Renumbering merged structure to contiguous 1..N...")
    new_res_num = 1
    # We need to track the mapping to generate the contig string correctly
    # Old ID -> New ID
    # But "Old ID" is (chain, resseq, icode).
    # The merged model has:
    # 1. Template Chain (args.chain) with a GAP where graft goes.
    # 2. Motif Chain (motif_chain_id).
    # We want to merge them into ONE chain if possible?
    # If we keep them as separate chains in PDB structures, we can still renumber them.
    # But for a SINGLE CHAIN OUTPUT from RFdiffusion, we should ideally have a single input chain.
    # Let's Move motif residues into args.chain!
    
    main_chain = merged_model[args.chain]
    if motif_chain_id and motif_chain_id != args.chain:
        motif_chain = merged_model[motif_chain_id]
        print(f"Merging motif chain {motif_chain_id} into {args.chain}...")
        for atom in list(motif_chain.get_atoms()):
            # We need to detach from old and add to new?
            # Easier to just move residues.
            pass # PDB structure modification is tricky inplace.
            
    # Simpler approach: Renumber everything in place, keeping chains distinct but numbers unique?
    # No, for `/` syntax, they usually imply same chain or connected.
    # Let's just renumber the residues and update the `merged_model` object.
    # And we need to know that 'Before' residues map to 1..X
    # 'Motif' residues map to Y..Z
    # 'After' residues map to ...
    
    # Let's collect all residues in desired order:
    # Before -> Motif -> After
    # But they are currently in different chains/gap.
    
    # Actually, we can just rely on the existing IDs for now but use `/` to connect them.
    # The `H1-52/2-6/M10-20` syntax works even if M is a different chain ID.
    # It tells RFdiffusion: "Take seq from H:1-52, then build 2-6 residues, then take seq from M:10-20".
    # Result will be a single chain.
    # So we DON'T strictly need to merge chains in the PDB.
    # WE DO need to handle duplicate 52/52A if they exist.
    
    # Scan for duplicates in merged_chain
    # If 52 and 52A exist, they will appear as `r.id[1]=52` twice.
    # BioPython parses 52A as `id=(' ', 52, 'A')`.
    # My `get_contig_segs` uses `r.id[1]`.
    # It ignores icode!
    # Fix `get_contig_segs` to handle full IDs or renumber.
    
    # Let's just fix the CONTIG GENERATION to handle insertion codes?
    # RFdiffusion understands `H52` but maybe not `H52A`.
    # If PDB has 52A, RFdiffusion might fail or ignore it.
    # Renumbering is the safest bet.
    
    # Implementation:
    # 1. Create a clean new chain.
    # 2. Add 'Before' residues from Template.
    # 3. Add 'Motif' residues.
    # 4. Add 'After' residues from Template.
    # 5. Save this new single-chain PDB.
    # 6. Generate contig `1-N/link/N+M...`.
    
    # This is complex to implement correctly in one go.
    # Let's stick to the `/` fix first.
    # And fix `get_contig_segs` to NOT duplicate 52.
    
    pass 

    # 5. Build Contig
    # [Template_Before, 10-20, Motif, 10-20, Template_After]
    
    
    
    
    # 5. Renumber Merged Structure & Build Contig
    # RFdiffusion works best with 1-indexed, contiguous residues.
    # We will renumber the MERGED model (template + graft) to 1..N.
    # This avoids insertion code issues and duplicate IDs.
    
    print("Renumbering merged structure to contiguous 1..N...")
    
    template_chain_obj = merged_model[args.chain]
    
    # Collect protein residues to add in order: Before -> Motif -> After
    residues_to_add = [] 
    
    # 1. Template Before
    all_temp_res = sorted(list(template_chain_obj.get_residues()), key=lambda r: r.id[1])
    temp_before = [r for r in all_temp_res if is_aa(r) and r.id[1] < start_res]
    for r in temp_before: residues_to_add.append((r, r.id))
    
    # 2. Motif
    motif_len = 0
    if motif_chain_id:
        motif_chain_obj = merged_model[motif_chain_id]
        motif_res_obj = sorted(list(motif_chain_obj.get_residues()), key=lambda r: r.id[1])
        motif_aa = [r for r in motif_res_obj if is_aa(r)]
        motif_len = len(motif_aa)
        
        # Create set for boundary checking later
        motif_residues_set = set(motif_aa)
        
        # Create set for boundary checking later
        motif_residues_set = set(motif_aa)
        
        for r in motif_aa: residues_to_add.append((r, r.id))
        
    # 3. Template After
    temp_after = [r for r in all_temp_res if is_aa(r) and r.id[1] > end_res]
    for r in temp_after: residues_to_add.append((r, r.id))
    
    # Prep for renumbering
    renumbered_structure = Structure.Structure("renumbered")
    renumbered_model = Model.Model(0)
    renumbered_structure.add(renumbered_model)
    renumbered_chain = Chain.Chain("A") # Always Chain A for protein
    renumbered_model.add(renumbered_chain)
    
    contig_parts = []
    current_segment_start_new = 1
    current_segment_end_new = 0 
    
    idx_switch_to_motif = len(temp_before)
    idx_switch_to_after = len(temp_before) + motif_len
    
    new_res_counter = 1
    
    for i, (r, old_id) in enumerate(residues_to_add):
        # 1. Renumber
        # Note: detach_child might fail if we already moved it? 
        # But we are building a NEW structure object usually.
        # But `r` is the same object.
        # We must clone it to be safe, or just move it.
        # Moving it is fine as long as we don't need old structure anymore.
        
        # FIX: Must detach from parent BEFORE changing ID to avoid collision with siblings!
        if r.get_parent():
            r.get_parent().detach_child(old_id)
            
        r.id = (' ', new_res_counter, ' ')
        renumbered_chain.add(r)
        
        # 2. Contig Logic
        if i == idx_switch_to_motif:
            if current_segment_end_new >= current_segment_start_new:
                 contig_parts.append(f"A{current_segment_start_new}-{current_segment_end_new}")
            contig_parts.append(args.linker_len) # Linker to Motif
            current_segment_start_new = new_res_counter
            
        elif i == idx_switch_to_after:
             if current_segment_end_new >= current_segment_start_new:
                 contig_parts.append(f"A{current_segment_start_new}-{current_segment_end_new}")
             contig_parts.append(args.linker_len) # Linker to Scaffold
             current_segment_start_new = new_res_counter
             
        if i > 0:
            # Check for Gap using Geometry (Robust to numbering weirdness)
            prev_res, prev_id = residues_to_add[i-1]
            curr_res, curr_id = residues_to_add[i]
            
            # Calculate distance between C of prev and N of curr
            # Note: At the Motif<->Template boundaries, we inserted linkers via loop logic (lines 642/649)
            # But here we are iterating residues_to_add which is the MERGED list.
            # Wait, residues_to_add is constructed sequentially.
            # Does it include the linkers? NO. Linkers are added to `contig_parts` string, not to the PDB structure yet.
            # RFdiffusion builds the linkers.
            # So `residues_to_add` contains the Fixed residues.
            # If `prev` and `curr` are far apart, we MUST insert a gap in `contig_parts`.
            
            # Check if we are crossing a linker boundary
            # If so, `current_segment_end_new` was updated in the linker block?
            # actually the linker block (lines 642/650) appends to `contig_parts`.
            # But `residues_to_add` is the list of RESIDUES.
            # We are iterating `residues_to_add`.
            # If `prev` was Template and `curr` is Motif -> Linker logic handles it (is_motif check).
            # If `prev` and `curr` are both Template -> Gap check needed.
            
            gap_needed = False
            gap_size_est = 0
            
            is_prev_motif = (prev_res in motif_residues_set)
            is_curr_motif = (curr_res in motif_residues_set)
            
            # Only check geometry if we are NOT crossing the graft boundaries (which are handled by linker_len)
            if is_prev_motif == is_curr_motif:
                 try:
                     c_atom = prev_res['C']
                     n_atom = curr_res['N']
                     dist = c_atom - n_atom
                     if dist > 2.0: # Threshold for peptide bond approx 1.33A
                         gap_needed = True
                         gap_size_est = int(dist / 3.0) # Approx 3A per residue stretch
                         if gap_size_est < 1: gap_size_est = 1
                 except KeyError:
                     # Missing atoms? Assume gap if numbering jumps?
                     # Fallback to numbering
                     if curr_id[1] - prev_id[1] > 1:
                         gap_needed = True
                         gap_size_est = curr_id[1] - prev_id[1] - 1

            if gap_needed:
                 # Close current segment
                 if current_segment_end_new >= current_segment_start_new:
                     contig_parts.append(f"A{current_segment_start_new}-{current_segment_end_new}")
                 
                 # Add gap
                 gap_start = gap_size_est
                 gap_end = gap_size_est + 2 if not args.rigid_gaps else gap_size_est
                 print(f"Detected geometric gap {dist:.2f}A. Auto-filling {gap_start}-{gap_end}.")
                 contig_parts.append(f"{gap_start}-{gap_end}")
                 
                 current_segment_start_new = new_res_counter
        
        current_segment_end_new = new_res_counter
        new_res_counter += 1
        
    if current_segment_end_new >= current_segment_start_new:
        contig_parts.append(f"A{current_segment_start_new}-{current_segment_end_new}")

    # Handle Metals (Move to Chain B)
    metal_contigs = []
    if metal_atoms_in_motif:
         print("Moving metals to Chain B...")
         metal_chain = Chain.Chain("B")
         renumbered_model.add(metal_chain)
         
         # These atoms might be in residues we already moved? 
         # No, we filtered `is_aa(r)`. Metals are HETATMs.
         # So they are still in `merged_model` old chains.
         
         # We need to find them. `metal_atoms_in_motif` has the atoms.
         # Get their parents (Residues).
         added_metals = set()
         metal_res_counter = 1
         
         for atom in metal_atoms_in_motif:
             res = atom.get_parent()
             # ID might be (H_Z, 100, )
             # Check if we already added this residue
             if res in added_metals: continue
             
             # Clone/Move
             res.detach_parent()
             
             # FIX 2: Set hetflag to 'H_XYZ' so PDBIO writes HETATM.
             # RFdiffusion parser (aa2num) crashes if it sees 'ATOM ... ZN'.
             # It ignores HETATM lines when parse_hetatom=False (default for target).
             hetflag = f"H_{res.resname.strip()}"
             res.id = (hetflag, metal_res_counter, ' ')
             metal_chain.add(res)
             
             # Note: Contig string B1-1 etc is NOT needed (and handled above)
             pdb_res_num = metal_res_counter
             
             # We track them just in case, but we don't put them in valid_contigs
             metal_contigs.append(f"B{metal_res_counter}-{metal_res_counter}")
             
             added_metals.add(res)
             metal_res_counter += 1

    # Replace merged structure
    merged_structure = renumbered_structure
    merged_model = renumbered_model
    merged_chain = renumbered_chain # Chain A
    args.chain = "A" 
    
    # Save
    merged_path = os.path.join(args.out_dir, "grafted_input.pdb")
    io = PDBIO()
    io.set_structure(merged_structure)
    io.save(merged_path)
    print(f"Renumbered structure saved to {merged_path}")
    
    # Add ORI (re-center if needed, but relative coords conserved)
    if metal_atoms_in_motif:
         coords = np.array([a.get_coord() for a in metal_atoms_in_motif])
         target_center = np.mean(coords, axis=0)
         
    with open(merged_path, 'a') as f:
        x, y, z = target_center
        ori_line = f"HETATM 9998  CA  ORI Z 999    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n"
        f.write(ori_line)

    
    
    # Build Final Contig String
    # Protein parts joined by ',' (single chain segments)
    # RFdiffusion contigs.py splits by ',' to get subcons. 
    # It does NOT split by '/'. Using '/' caused ValueError in length calculation.
    protein_contig = ",".join(contig_parts)
    
    # NOTE: We do NOT include the metal chain (B) in the contig string.
    # RFdiffusion treats items in `contigs` as protein chains to be diffused or fixed.
    # If we include 'B1-1' where B1 is ZN, it tries to parse ZN as an amino acid and fails (KeyError 'ZN').
    # Instead, we rely on `inference.ligand` to tell RFdiffusion to respect the ligand in the input PDB.
    # The metal is preserved in the PDB (Chain B), so RFdiffusion will see it.
    
    contig_str = protein_contig

    contig = f"['{contig_str}']"
    
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
    
    if args.overwrite and existing_designs:
        print(f"Overwriting: Removing {len(existing_designs)} existing designs in {rf_out}...")
        for f in existing_designs:
            try:
                os.remove(f)
                # Also remove trb file
                trb = f.replace(".pdb", ".trb")
                if os.path.exists(trb): os.remove(trb)
            except Exception as e:
                print(f"Error removing {f}: {e}")

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
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output directory")
    parser.add_argument("--linker_len", type=str, default="5-15", help="Length range for linkers connecting motif to scaffold (default: 5-15)")
    parser.add_argument("--rigid_gaps", action="store_true", help="Disable flexibility in auto-filled gaps (force fixed length)")
    parser.add_argument("--clash_threshold", type=float, default=3.5, help="Distance threshold (Angstroms) for cleaning clashing scaffold residues (default: 3.5)")
    parser.add_argument("--no_clash_cleanup", action="store_true", help="Disable automatic removal of clashing scaffold residues")
    
    args = parser.parse_args()
    run_graft(args)

if __name__ == "__main__":
    main()
