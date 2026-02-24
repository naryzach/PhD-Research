import numpy as np
import biotite.structure as struc


def get_ef_hand_loops(atom_array, window_size=15, metal_res_name="ND"):
    """
    Finds all 'ND' ions in the structure and extracts their associated loops.
    EF hand loops are typically 12 residues long. We will mask out an arbitrary
    window of `window_size` residues around the nearest neighbor atom to the ion for redesign.
    """
    metal_atoms = atom_array[atom_array.res_name == metal_res_name]
    loops = []
    
    for i, metal in enumerate(metal_atoms):
        # Find nearest protein atom to this metal to find the center of the binding loop
        protein_atoms = atom_array[
            np.isin(atom_array.res_name, ["HOH", metal_res_name], invert=True)
        ]
        
        if len(protein_atoms) == 0:
            continue
            
        distances = np.linalg.norm(protein_atoms.coord - metal.coord, axis=1)
        closest_idx = np.argmin(distances)
        closest_atom = protein_atoms[closest_idx]
        
        # Center the loop around this residue
        center_res_id = closest_atom.res_id
        chain_id = closest_atom.chain_id
        
        print(f"Metal {i+1} (Chain {metal.chain_id}, ResID {metal.res_id}) -> Closest protein res: Chain {chain_id}, ResID {center_res_id}")
        
        loops.append({
            "metal_idx": i,
            "metal_res_id": metal.res_id,
            "metal_chain_id": metal.chain_id,
            "center_res_id": center_res_id,
            "protein_chain_id": chain_id,
            "window_size": window_size,
            "coord": metal.coord
        })
    return loops


def create_masked_input(atom_array, loop_info):
    """
    Creates an input AtomArray for RFd3 where the loop residues are removed.
    """
    protein_chain = loop_info["protein_chain_id"]
    center_res_id = loop_info["center_res_id"]
    window_size = loop_info["window_size"]
    
    start_res = max(1, center_res_id - window_size // 2)
    end_res = center_res_id + window_size // 2
    
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
