import os
import argparse
from string import ascii_uppercase
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select, Chain

# Metals and their formal charges
METAL_IONS = {
    "ZN": 2, "MG": 2, "MN": 2, "FE": 2, "CU": 2, "CA": 2, "CO": 2,
    "NI": 2, "CD": 2, "K": 1, "NA": 1, "HG": 2, "SR": 2, "BA": 2
}

# Standard amino-acid 3-letter codes
STANDARD_AA3 = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
}

WATER_NAMES = {"HOH", "WAT", "H2O"}


class CleanSelect(Select):
    """Return True for atoms we want to write.
    Rules:
      - Always remove waters
      - Keep protein ATOM records
      - Keep metal HETATM records only if keep_hetatms==True
      - **HADDOCK FIX:** Only keep ONE atom per metal residue
      - Remove any other HETATM (non-standard cofactor)
      - Respect altloc: only keep altloc 'A' or single-position atoms
    """
    def __init__(self, keep_hetatms=False):
        super().__init__()
        self.keep_hetatms = keep_hetatms
        # --- ADD THIS ---
        # This set will store the renumbered IDs of metal residues
        # we have already written an atom for (e.g., ('H_', 1, ' '))
        self.seen_metal_residues = set()

    def accept_atom(self, atom):
        residue = atom.get_parent()
        resname = residue.get_resname().strip().upper()
        base = "".join(ch for ch in resname if ch.isalpha())

        # remove water
        if base in WATER_NAMES:
            return False

        # --- REVISED METAL LOGIC ---
        if base in METAL_IONS:
            if not self.keep_hetatms:
                return False

            # HADDOCK Error Fix: Check if we've already kept an atom for this
            # renumbered metal residue.
            if residue.id in self.seen_metal_residues:
                # If yes, reject this atom (e.g., it's an altloc 'B'
                # or a duplicate non-disordered atom)
                return False

            # This is the first atom we've seen for this residue.
            # We must still check for altlocs on it.
            if atom.is_disordered() and atom.get_altloc() != "A":
                # It's disordered, but not 'A', so we reject it.
                # We don't add it to 'seen' because we haven't kept it.
                return False

            # This is an atom we are keeping.
            # Add its residue ID to the 'seen' set to block all others.
            self.seen_metal_residues.add(residue.id)
            return True
        # --- END OF REVISED LOGIC ---

        # skip any other hetero residue (non-protein/non-metal)
        if residue.id[0] != " ":
            return False

        # altloc check for standard protein residues
        if atom.is_disordered() and atom.get_altloc() != "A":
            return False

        return True


def _post_process_for_haddock(filepath, metal_charges_dict):
    """
    Reads the output PDB file, finds HETATM lines for metals,
    and reformats the atom name to include the charge (e.g., ZN+2).
    This is a post-processing step to bypass Biopython's PDBIO writer.
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        processed_lines = []
        for line in lines:
            if line.startswith("HETATM"):
                # Atom name written by Biopython (cols 13-16)
                atom_name_from_pdb = line[12:16].strip()
                # Res name written by Biopython (cols 18-20)
                res_name_from_pdb = line[17:20].strip()

                # Use the atom name (or res name as fallback) as the key
                base_element = None
                if atom_name_from_pdb in metal_charges_dict:
                    base_element = atom_name_from_pdb
                elif res_name_from_pdb in metal_charges_dict:
                    base_element = res_name_from_pdb

                if base_element:
                    charge = metal_charges_dict.get(base_element)
                    if charge is not None:
                        # HADDOCK name, e.g., "ZN+2"
                        haddock_name = f"{base_element}+{charge}"
                        # Format to 4 chars, left-justified, e.g., "ZN+2"
                        haddock_name_formatted = "{:<4}".format(haddock_name)
                        
                        # Splice the new name into the line
                        line = line[:12] + haddock_name_formatted + line[16:]
            
            processed_lines.append(line)

        # Write the modified lines back to the same file
        with open(filepath, 'w') as f:
            f.writelines(processed_lines)
    
    except Exception as e:
        print(f"Error during HADDOCK post-processing: {e}")


def clean_and_renumber(
    input_file, output_file,
    start_res=None, end_res=None,
    chains_to_keep=None, keep_hetatms=False,
    haddock=False
):
    """Clean and renumber PDB/mmCIF for HADDOCK with strict metal handling."""
    ext = os.path.splitext(input_file)[1].lower()
    if ext == ".cif":
        parser = MMCIFParser(QUIET=True)
    elif ext == ".pdb":
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError("Unsupported input format: must be .pdb or .cif")

    structure = parser.get_structure("structure", input_file)
    model = structure[0]

    # ... (Code for filtering chains, trimming, etc. - no changes needed here) ...
    # Filter chains if requested (use original chain IDs)
    if chains_to_keep:
        chains_to_keep = [c.strip().upper() for c in chains_to_keep]
        for chain in list(model):
            if chain.id.upper() not in chains_to_keep:
                model.detach_child(chain.id)
        print(f"Keeping chains: {', '.join(chains_to_keep)}")

    # Rename protein chains sequentially from 'A'
    chain_labels = iter(ascii_uppercase)
    protein_chains = []
    for chain in list(model):
        protein_chains.append(chain)
    for chain in protein_chains:
        chain.id = next(chain_labels, "Z") # Assign A, B, C...

    # Optionally trim residues
    if start_res is not None and end_res is not None:
        for chain in list(model):
            residues_to_keep = [res for res in chain if start_res <= res.id[1] <= end_res]
            for res in list(chain):
                if res not in residues_to_keep:
                    try:
                        chain.detach_child(res.id)
                    except Exception:
                        pass

    # --- Handle residues: collect metals, remove nonstandard ---
    metals_to_move = []
    metals_chain_id = "Z" # placeholder
    for chain in list(model):
        for residue in list(chain):
            orig_resname = residue.get_resname().strip().upper()
            base = "".join(ch for ch in orig_resname if ch.isalpha())

            if base in WATER_NAMES: # remove waters
                try: chain.detach_child(residue.id)
                except Exception: pass
                continue
            
            if base in METAL_IONS: # extract metals
                metals_to_move.append((residue, chain.id))
                try: chain.detach_child(residue.id)
                except Exception: pass
                continue
            
            if base not in STANDARD_AA3: # remove nonstandard AAs
                try: chain.detach_child(residue.id)
                except Exception: pass
                continue

    # --- **FIX**: Renumber protein residues CONTINUOUSLY starting at 1 ---
    current_residue_number = 1
    for chain in model:
        # Get all protein residues in this chain, in order
        protein_residues = [r for r in chain.get_unpacked_list() if r.id[0] == " "]
        if not protein_residues:
            continue # Skip chains that are now empty (e.g., were all-metal)

        for residue in protein_residues:
            residue.id = (" ", current_residue_number, " ")
            current_residue_number += 1
    
    # current_residue_number now holds the *next* available residue ID
    # (e.g., if last res was 159, this is 160)
    # --- End of FIX ---

    # --- Create a new chain for metals if any and if keep_hetatms True ---
    if metals_to_move and keep_hetatms:
        metals_chain_id = next(chain_labels, "Z")
        metals_chain = Chain.Chain(metals_chain_id)
        
        # --- **FIX**: Use the continuous residue number from above ---
        metal_residue_counter = current_residue_number
        
        # --- HADDOCK FIX: Dictionary to count metal types ---
        metal_type_counts = {}
        unique_suffix_chars = "123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"

        for residue, orig_chain in metals_to_move:
            orig_resname = residue.get_resname().strip().upper()
            base = "".join(ch for ch in orig_resname if ch.isalpha()) 

            # Create unique residue name
            count = metal_type_counts.get(base, 0) + 1
            metal_type_counts[base] = count

            if count > len(unique_suffix_chars):
                new_resname = f"{base[0]}{count:02d}" 
            else:
                suffix = unique_suffix_chars[count-1]
                new_resname = f"{base[:2]}{suffix}" 
            
            new_resname = "{:<3}".format(new_resname[:3])

            # --- **FIX**: Set residue as hetero AND use continuous renumbering ---
            residue.id = ('H_', metal_residue_counter, ' ')
            metal_residue_counter += 1 # Increment the global counter
            # --- End of FIX ---

            residue.resname = new_resname 

            for atom in residue:
                atom.name = base
                atom.element = base 

            metals_chain.add(residue)
        
        model.add(metals_chain)

    # --- Save using the STANDARD PDBIO ---
    io = PDBIO() 
    io.set_structure(structure)
    io.save(output_file, select=CleanSelect(keep_hetatms=keep_hetatms))

    # --- HADDOCK Post-Processing Step ---
    if haddock and keep_hetatms and metals_to_move:
        print("   - Applying HADDOCK post-processing for metal charges & serial numbers...")
        _post_process_for_haddock(output_file, METAL_IONS)

    print(f"\nCleaned structure saved to: {output_file}")
    print("   - Output format: PDB")
    print("   - Protein chains renamed sequentially starting from 'A'")
    if metals_to_move:
        print(f"   - Metals extracted into chain {metals_chain_id if keep_hetatms else '(removed)'}")
    print("   - All residues renumbered continuously across chains") # <-- Updated message
    if haddock and metals_to_move:
        print("   - Metal residue names and charges formatted for HADDOCK")

def batch_process(directory, start_res=None, end_res=None, chains_to_keep=None, keep_hetatms=False, haddock=False): # <-- NEW FLAG
    files = [f for f in os.listdir(directory) if f.lower().endswith((".pdb", ".cif"))]
    if not files:
        print("No PDB or CIF files found in the specified directory.")
        return

    output_dir = os.path.join(directory, "cleaned_outputs")
    os.makedirs(output_dir, exist_ok=True)

    for file in files:
        input_path = os.path.join(directory, file)
        base_name = os.path.splitext(file)[0]
        chain_suffix = ""
        if chains_to_keep:
            chain_suffix = "_" + "_".join([c.strip().upper() for c in chains_to_keep])
        
        # Add haddock suffix if requested
        hetatm_suffix = "_withMetals" if keep_hetatms else ""
        if haddock and keep_hetatms:
             hetatm_suffix = "_HADDOCK"

        output_name = f"{base_name}{chain_suffix}{hetatm_suffix}_cleaned.pdb"
        output_path = os.path.join(output_dir, output_name)

        print(f"\nProcessing {file}...")
        try:
            clean_and_renumber(
                input_file=input_path,
                output_file=output_path,
                start_res=start_res,
                end_res=end_res,
                chains_to_keep=chains_to_keep,
                keep_hetatms=keep_hetatms,
                haddock=haddock # <-- PASS FLAG
            )
        except Exception as e:
            print(f"Error processing {file}: {e}")

    print(f"\nBatch processing complete. Cleaned files saved to: {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean and renumber PDB/mmCIF structures for HADDOCK.")
    parser.add_argument("input", help="Input PDB/CIF file or directory")
    parser.add_argument("output", nargs="?", help="Output cleaned PDB file (ignored in batch mode)")
    parser.add_argument("--start", type=int, help="Start residue number to keep (optional)")
    parser.add_argument("--end", type=int, help="End residue number to keep (optional)")
    parser.add_argument("--chains", nargs="+", help="List of chain IDs to keep (e.g. --chains A B)")
    parser.add_argument("--keep-hetatms", action="store_true", help="Keep and format metal ions")
    
    # --- NEW ARGUMENT ---
    parser.add_argument("--haddock", action="store_true", help="Format metal HETATM records for HADDOCK (e.g., ZN+2). Implies --keep-hetatms.")
    
    parser.add_argument("--batch", action="store_true", help="Process all .pdb and .cif files in a directory")

    args = parser.parse_args()

    # If --haddock is used, it must also --keep-hetatms
    if args.haddock:
        args.keep_hetatms = True

    if args.batch:
        batch_process(
            directory=args.input,
            start_res=args.start,
            end_res=args.end,
            chains_to_keep=args.chains,
            keep_hetatms=args.keep_hetatms,
            haddock=args.haddock
        )
    else:
        if not args.output:
            raise ValueError("Output file must be specified in single-file mode.")
        clean_and_renumber(
            input_file=args.input,
            output_file=args.output,
            start_res=args.start,
            end_res=args.end,
            chains_to_keep=args.chains,
            keep_hetatms=args.keep_hetatms,
            haddock=args.haddock 
        )