import os
import argparse
from string import ascii_uppercase
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select

# Common biologically relevant metal ions and their formal charges
METAL_IONS = {
    "ZN": 2, "MG": 2, "MN": 2, "FE": 2, "CU": 2, "CA": 2, "CO": 2,
    "NI": 2, "CD": 2, "K": 1, "NA": 1, "HG": 2, "SR": 2, "BA": 2
}


class CleanSelect(Select):
    """Selects which atoms to keep based on cleanup and --keep-hetatms flag."""
    def __init__(self, keep_hetatms=False):
        super().__init__()
        self.keep_hetatms = keep_hetatms

    def accept_atom(self, atom):
        parent = atom.get_parent()
        resname = parent.get_resname().strip().upper()

        # Always skip water molecules
        if resname in {"HOH", "WAT", "H2O"}:
            return False

        # Keep metal ions if requested
        if self.keep_hetatms and resname in METAL_IONS:
            return True

        # Skip all other HETATMs unless they’re allowed metals
        if parent.id[0] != " ":
            return False

        # Skip alternate locations except 'A'
        if atom.is_disordered() and atom.get_altloc() != "A":
            return False

        return True


def clean_and_renumber(
    input_file, output_file,
    start_res=None, end_res=None,
    chains_to_keep=None, keep_hetatms=False
):
    """Clean, trim, and renumber a PDB/mmCIF structure for HADDOCK."""
    ext = os.path.splitext(input_file)[1].lower()
    if ext == ".cif":
        parser = MMCIFParser(QUIET=True)
    elif ext == ".pdb":
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported input format: {ext}")

    structure = parser.get_structure("structure", input_file)
    model = structure[0]

    # Filter chains if requested
    if chains_to_keep:
        chains_to_keep = [c.strip().upper() for c in chains_to_keep]
        for chain in list(model):
            if chain.id.upper() not in chains_to_keep:
                model.detach_child(chain.id)
        print(f"Keeping chains: {', '.join(chains_to_keep)}")

    # Rename chains sequentially from A
    chain_labels = iter(ascii_uppercase)
    for chain in model:
        chain.id = next(chain_labels, "Z")

    # Trim residue span if requested
    if start_res is not None and end_res is not None:
        for chain in list(model):
            residues_to_keep = [
                res for res in chain
                if start_res <= res.id[1] <= end_res
            ]
            for res in list(chain):
                if res not in residues_to_keep:
                    chain.detach_child(res.id)

    # Renumber residues starting at 1 for each chain
    for chain in model:
        new_id = 1
        for residue in chain.get_unpacked_list():
            residue.id = (" ", new_id, " ")
            new_id += 1

    # --- Normalize metal ions to HETATM for HADDOCK ---
    if keep_hetatms:
        for chain in model:
            for residue in list(chain):
                resname = residue.get_resname().strip().upper()
                if resname in METAL_IONS:
                    charge = METAL_IONS[resname]
                    new_resname = f"{resname}{charge}+"
                    residue.resname = new_resname
                    for atom in residue:
                        # Force them to be treated as heteroatoms
                        atom.parent.id = ("H_", atom.parent.id[1], atom.parent.id[2])

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file, select=CleanSelect(keep_hetatms=keep_hetatms))

    print(f"\nCleaned structure saved to: {output_file}")
    print("   - Output format: PDB")
    print("   - Chains renamed sequentially starting from 'A'")
    if start_res and end_res:
        print(f"   - Residues kept: {start_res}–{end_res} → renumbered 1–{end_res - start_res + 1}")
    print(f"   - Metal ions kept: {'Yes' if keep_hetatms else 'No'}")
    print("   - Waters and other heteroatoms removed")
    print("   - Nonstandard atoms, altlocs, and naming fixed")


def batch_process(directory, start_res=None, end_res=None, chains_to_keep=None, keep_hetatms=False):
    """Process all .pdb and .cif files in a directory."""
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
        hetatm_suffix = "_withMetals" if keep_hetatms else ""

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
            )
        except Exception as e:
            print(f"Error processing {file}: {e}")

    print(f"\nBatch processing complete. Cleaned files saved to: {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Clean and renumber PDB/mmCIF structures for HADDOCK."
    )
    parser.add_argument("input", help="Input PDB/CIF file or directory")
    parser.add_argument("output", nargs="?", help="Output cleaned PDB file (ignored in batch mode)")
    parser.add_argument("--start", type=int, help="Start residue number to keep (optional)")
    parser.add_argument("--end", type=int, help="End residue number to keep (optional)")
    parser.add_argument("--chains", nargs="+", help="List of chain IDs to keep (e.g. --chains A B)")
    parser.add_argument("--keep-hetatms", action="store_true",
                        help="Keep and correctly format metal ions for HADDOCK")
    parser.add_argument("--batch", action="store_true",
                        help="Process all .pdb and .cif files in a directory")

    args = parser.parse_args()

    if args.batch:
        batch_process(
            directory=args.input,
            start_res=args.start,
            end_res=args.end,
            chains_to_keep=args.chains,
            keep_hetatms=args.keep_hetatms,
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
        )
