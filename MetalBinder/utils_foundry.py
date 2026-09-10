import numpy as np
import biotite.structure as struc


def get_ef_hand_loops(atom_array, window_size=15, metal_res_name="ND",
                       ef_ranges=None, protein_chain_id="A"):
    """
    Finds all specified metal ions in the structure and extracts their associated loops.
    EF hand loops are typically 12 residues long. We will mask out an arbitrary
    window of `window_size` residues around the nearest neighbor atom to the ion for redesign.

    `ef_ranges` / `protein_chain_id` let this be reused for a different template
    structure than 8FNS (e.g. Calmodulin) without changing the 8FNS default
    behavior for any existing caller -- omit both to get the original
    hardcoded 8FNS (chain A) ranges exactly as before.
    """
    metal_atoms = atom_array[atom_array.res_name == metal_res_name]
    loops = []

    if ef_ranges is None:
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
            "protein_chain_id": protein_chain_id,
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

# Shannon effective ionic radii (six-coordinate, Angstrom; Shannon, Acta Cryst. 1976).
# Fe/Mn/Co use high-spin values (dominant in biological contexts). All lanthanides
# taken at their common +3 oxidation state, matching METAL_SMILES in score_with_chai1.py.
IONIC_RADII_CN6 = {
    "CA": 1.00, "MG": 0.72, "ZN": 0.74, "CU": 0.73, "FE": 0.78, "FE3": 0.645,
    "MN": 0.83, "CO": 0.745, "NI": 0.69,
    "LA": 1.032, "CE": 1.01, "PR": 0.99, "ND": 0.983, "PM": 0.97, "SM": 0.958,
    "EU": 0.947, "GD": 0.938, "TB": 0.923, "DY": 0.912, "HO": 0.901, "ER": 0.89,
    "TM": 0.88, "YB": 0.868, "LU": 0.861, "Y": 0.90, "SC": 0.745,
}
OXYGEN_RADIUS_CN6 = 1.40  # O2- effective ionic radius, six-coordinate


def ideal_mo_distance(ion: str) -> float:
    """Sum-of-ionic-radii estimate of the ideal metal-oxygen bond length for `ion`.
    Returns NaN for an ion not in IONIC_RADII_CN6 rather than guessing."""
    r = IONIC_RADII_CN6.get(ion.upper())
    return r + OXYGEN_RADIUS_CN6 if r is not None else float("nan")


def calculate_binding_metrics(atom_array, loop_info, start_res, end_res, cutoff=3.1):
    """
    Calculates advanced binding metrics for the target metal ion.
    Returns a dictionary of metrics.

    `cutoff` (coordination-sphere cutoff, Angstrom) defaults to a single fixed
    3.1 A for all ions, which is coarse: it can't distinguish coordination
    number differences that fall inside the gap between e.g. a lanthanide's
    ideal M-O distance (~2.3-2.4 A) and Ca's (~2.4 A) vs. the cutoff itself.
    For ion-aware analysis, pass `cutoff=ideal_mo_distance(ion) + tolerance`
    (e.g. +0.4 A) instead — see ideal_mo_distance() above.
    """
    metrics = {
        "binding_radius_A": float('nan'),
        "coordination_number": 0,
        "net_charge": 0,
        "bidentate_count": 0
    }

    # 1. Locate the specific metal ion
    metal_mask = (atom_array.res_id == loop_info['metal_res_id']) & \
                 (atom_array.chain_id == loop_info['metal_chain_id'])
    metal_atoms = atom_array[metal_mask]
    if len(metal_atoms) == 0:
        return metrics
    metal_coord = metal_atoms.coord[0]

    # 2. Isolate the loop residues
    loop_mask = (atom_array.chain_id == loop_info["protein_chain_id"]) & \
                (atom_array.res_id >= start_res) & \
                (atom_array.res_id <= end_res)
    loop_atoms = atom_array[loop_mask]
    if len(loop_atoms) == 0:
        return metrics

    # 3. Calculate Net Charge of the loop
    # Asp (D) and Glu (E) are -1, Lys (K) and Arg (R) are +1, His (H) is ~0.1 (simplified)
    # We focus on the formal charges.
    # struc.get_residues() returns (res_ids, res_names) in that order — this was
    # previously unpacked backwards, so `rn` held integer residue ids and never
    # matched the ASP/GLU/LYS/ARG string checks below, silently giving net_charge=0
    # for every single loop regardless of actual sequence.
    res_ids, res_names = struc.get_residues(loop_atoms)
    charge = 0
    for rn in res_names:
        if rn in ['ASP', 'GLU']: charge -= 1
        elif rn in ['LYS', 'ARG']: charge += 1
    metrics["net_charge"] = charge

    # 4. Filter for coordinating Oxygens
    # EF-hands coordinate via Oxygen (O, OD1, OD2, OE1, OE2)
    potential_coord_mask = np.isin(loop_atoms.atom_name, ['O', 'OD1', 'OD2', 'OE1', 'OE2'])
    coordinating_atoms = loop_atoms[potential_coord_mask]
    
    if len(coordinating_atoms) == 0:
        return metrics

    # 5. Calculate Distances
    distances = np.linalg.norm(coordinating_atoms.coord - metal_coord, axis=1)

    # Coordination Sphere Cutoff (see `cutoff` param docstring above)
    cn_mask = distances <= cutoff
    metrics["coordination_number"] = int(np.sum(cn_mask))
    
    # 6. Binding Radius (Mean of closest 6 contacts)
    closest_distances = np.sort(distances)[:6]
    if len(closest_distances) > 0:
        metrics["binding_radius_A"] = np.mean(closest_distances)

    # 7. Bidentate Detection
    # If both oxygens of an Asp/Glu are within the cutoff, it's bidentate.
    bidentate_count = 0
    for res_id in np.unique(coordinating_atoms.res_id):
        res_atoms = coordinating_atoms[coordinating_atoms.res_id == res_id]
        # Check for OD1/OD2 (Asp) or OE1/OE2 (Glu)
        sidechain_ox_mask = np.isin(res_atoms.atom_name, ['OD1', 'OD2', 'OE1', 'OE2'])
        sc_oxygens = res_atoms[sidechain_ox_mask]
        if len(sc_oxygens) >= 2:
            sc_distances = np.linalg.norm(sc_oxygens.coord - metal_coord, axis=1)
            if np.all(sc_distances <= cutoff):
                bidentate_count += 1
    metrics["bidentate_count"] = bidentate_count

    return metrics

