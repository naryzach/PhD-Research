import numpy as np
import biotite.structure as struc


def get_ef_hand_loops(atom_array, window_size=15, metal_res_name="ND"):
    """
    Finds all specified metal ions in the structure and extracts their associated loops.
    EF hand loops are typically 12 residues long. We will mask out an arbitrary
    window of `window_size` residues around the nearest neighbor atom to the ion for redesign.
    """
    metal_atoms = atom_array[atom_array.res_name == metal_res_name]
    loops = []
    
    # Pre-defined explicit EF hand ranges for 8FNS (Chain A)
    ef_ranges = [
        (35, 46),
        (59, 70),
        (84, 95),
        (108, 119)
    ]
    
    for i, metal in enumerate(metal_atoms):
        if i < len(ef_ranges):
            start_res, end_res = ef_ranges[i]
        else:
            # Fallback for extra metals, if any
            start_res, end_res = 1, 10
            
        print(f"Metal {i+1} (Chain {metal.chain_id}, ResID {metal.res_id}) -> Assigned to EF Hand loop {start_res}-{end_res}")
        
        loops.append({
            "metal_idx": i,
            "metal_res_name": metal.res_name,
            "metal_res_id": metal.res_id,
            "metal_chain_id": metal.chain_id,
            "start_res": start_res,
            "end_res": end_res,
            "protein_chain_id": "A", # 8FNS is chain A
            "coord": metal.coord
        })
    return loops


def mutate_metals(atom_array, original_res_name, new_res_name):
    """
    Finds all metals with original_res_name and mutates their atom_name, res_name, and element to new_res_name.
    Returns the modified AtomArray.
    """
    metal_mask = atom_array.res_name == original_res_name
    
    if np.any(metal_mask):
        print(f"Mutating {np.sum(metal_mask)} '{original_res_name}' ions to '{new_res_name}'...")
        atom_array.res_name[metal_mask] = new_res_name
        # Typically the atom name matches the element/residue name, keeping it simple
        atom_array.atom_name[metal_mask] = new_res_name
        # Keep element clean (remove any numbers if user passes something like EU3)
        element = ''.join([c for c in new_res_name if c.isalpha()]).upper()
        atom_array.element[metal_mask] = element
        
    return atom_array


def create_masked_input(atom_array, loop_info):
    """
    Creates an input AtomArray for RFd3 where the loop residues are removed.
    """
    protein_chain = loop_info["protein_chain_id"]
    start_res = loop_info["start_res"]
    end_res = loop_info["end_res"]
    
    # Keep atoms that are NOT in the loop we want to redesign
    mask = ~(
        (atom_array.chain_id == protein_chain) & 
        (atom_array.res_id >= start_res) & 
        (atom_array.res_id <= end_res)
    )
    
    # Ensure our specific metal is kept
    metal_mask = (atom_array.chain_id == loop_info["metal_chain_id"]) & (atom_array.res_id == loop_info["metal_res_id"])
    mask = mask | metal_mask
    
    # Remove waters
    water_mask = atom_array.res_name == "HOH"
    mask = mask & ~water_mask
    
    # Create the fixed context
    context_array = atom_array[mask]
    return context_array, start_res, end_res

def calculate_binding_radius(atom_array, loop_info, start_res, end_res):
    """
    Calculates the effective binding radius between the target metal ion 
    and the nearest coordinating oxygen atoms in the designed loop.
    """
    # 1. Locate the specific metal ion in the structure
    metal_mask = (atom_array.res_id == loop_info['metal_res_id']) & \
                 (atom_array.chain_id == loop_info['metal_chain_id'])
    
    metal_atoms = atom_array[metal_mask]
    
    if len(metal_atoms) == 0:
        print(f"Warning: Metal ion not found at {loop_info['metal_chain_id']}:{loop_info['metal_res_id']}")
        return float('nan')
        
    # Get the 3D coordinates of the metal ion
    metal_coord = metal_atoms.coord[0]
    
    # 2. Isolate the newly designed loop residues
    loop_mask = (atom_array.chain_id == loop_info["protein_chain_id"]) & \
                (atom_array.res_id >= start_res) & \
                (atom_array.res_id <= end_res)
    
    loop_atoms = atom_array[loop_mask]
    
    if len(loop_atoms) == 0:
        return float('nan')
        
    # 3. Filter for likely coordinating atoms
    # Lanthanides and EF-hands primarily coordinate via Oxygen (O, OD1, OD2, OE1, OE2)
    oxygen_mask = np.char.startswith(loop_atoms.atom_name.astype(str), 'O')
    coordinating_atoms = loop_atoms[oxygen_mask]
    
    # Fallback to all heavy atoms if no oxygens are present (e.g., if sidechains failed to pack)
    if len(coordinating_atoms) == 0:
        coordinating_atoms = loop_atoms
        
    # 4. Calculate Euclidean distances from the metal to all coordinating atoms
    distances = np.linalg.norm(coordinating_atoms.coord - metal_coord, axis=1)
    
    # 5. Calculate the radius of the coordination sphere
    # Lanthanides typically have a coordination number of 7 to 9.
    # We take the mean of the 6 closest contacts to establish the primary binding radius.
    closest_distances = np.sort(distances)[:6]
    
    if len(closest_distances) > 0:
        return np.mean(closest_distances)
    else:
        return float('nan')

