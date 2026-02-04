import argparse
import os
import sys
import numpy as np
import subprocess
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Structure, Model, Chain
from Bio.PDB.Polypeptide import is_aa

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
    
    # Find insertion point
    target_res = template_chain[args.insert_at]
    target_center = target_res['CA'].get_coord() if 'CA' in target_res else get_center_of_mass(target_res)
    
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
    # We move Motif atoms into a new chain (or keep same if distinct)
    # RFdiffusion treats separate chains in input PDB as inputs?
    # Yes.
    # Let's verify Motif Chain ID.
    motif_chain_id = list(motif_model.child_dict.keys())[0]
    if motif_chain_id == args.chain:
        print(f"Warning: Motif chain {motif_chain_id} conflicts with Template chain {args.chain}. Renaming Motif to 'M'.")
        motif_model[motif_chain_id].id = 'M'
        motif_chain_id = 'M'
        
    # Create Merged Structure
    merged_path = os.path.join(args.out_dir, "grafted_input.pdb")
    os.makedirs(args.out_dir, exist_ok=True)
    
    io = PDBIO()
    
    # Merge logic: Add chains from Motif to Template Model?
    # BioPython structures are hierarchical.
    # Copy template model.
    merged_struct = template_struct
    # Add motif chain(s)
    # Need to detach from motif struct first?
    pass # BioPython entities can be reparented?
    
    # Easier: Just write both to file.
    with open(merged_path, 'w') as out_f:
        # Write Template
        with open(args.template, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    chain_id = line[21]
                    if chain_id == args.chain:
                        out_f.write(line)
        out_f.write("TER\n")
        
        # Write Motif (using PDBIO for transformed coords)
        # We need to save the transformed motif to a temp file first
        temp_motif = os.path.join(args.out_dir, "temp_motif_transformed.pdb")
        io.set_structure(motif_struct)
        io.save(temp_motif)
        
        with open(temp_motif, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # Ensure chain matching?
                    # The IO saved it with 'M' if we renamed it.
                    out_f.write(line)
        
        os.remove(temp_motif)
        
    print(f"Merged structure saved to {merged_path}")
    
    # 5. Build Contig
    # [Template_Before / Loop / Motif / Loop / Template_After]
    # Lengths:
    # Template Before: 1 to (insert_at - 1)
    # Loop: 10-20
    # Motif: Fixed (Length of Motif)
    # Loop: 10-20
    # Template After: (insert_at + 1) to End
    
    # Get Template max res
    res_list = sorted([r.id[1] for r in template_chain if is_aa(r)])
    max_res = res_list[-1]
    
    contig_parts = []
    
    # Template Before
    if args.insert_at > 1:
        contig_parts.append(f"{args.chain}1-{args.insert_at - 1}")
    
    contig_parts.append("10-20") # Loop 1
    
    # Motif (Assume Chain M or whatever it is)
    # Motif might be disjoint!
    # We need to use logic from run_scaffold to get segments?
    # Simplification: Assume Motif is Chain 'M' and we fix ALL of it.
    # Get Motif max res
    m_res = sorted([r.id[1] for r in motif_struct[0][motif_chain_id] if is_aa(r)])
    if m_res:
         start, end = m_res[0], m_res[-1]
         contig_parts.append(f"{motif_chain_id}{start}-{end}")
    
    contig_parts.append("10-20") # Loop 2
    
    # Template After
    if args.insert_at < max_res:
        contig_parts.append(f"{args.chain}{args.insert_at + 1}-{max_res}")
        
    contig_parts.append("Z0-0") # Metal? Assuming it is in Motif HETATM and RFdiffusion keeps it? 
    # RFdiffusion needs "Inpaint Seq" or just contig?
    # Handling metal in contig:
    # If the metal is Chain Z in the Merged PDB, we append it.
    
    contig = f"['{'/'.join(contig_parts)}']"
    
    # 6. Run RFdiffusion
    rf_out = os.path.join(args.out_dir, "rfdiffusion")
    os.makedirs(rf_out, exist_ok=True)
    
    rf_script = os.path.join(args.rf_path, "scripts/run_inference.py")
    cmd = [
        sys.executable, rf_script,
        f"inference.input_pdb={merged_path}",
        f"inference.output_prefix={os.path.join(rf_out, 'design')}",
        f"inference.num_designs={args.num_designs}",
        f"contigmap.contigs={contig}"
    ]
    
    if args.dry_run:
        print(f"Dry Run RFdiffusion: {' '.join(cmd)}")
    else:
        try:
            subprocess.run(cmd, check=True, cwd=project_root)
        except Exception as e:
            print(f"RFdiffusion failed: {e}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--motif", required=True, help="Input Motif PDB")
    parser.add_argument("--template", required=True, help="Input Template PDB (Antibody)")
    parser.add_argument("--out_dir", required=True, help="Output directory")
    parser.add_argument("--chain", default="H", help="Template Chain to graft into")
    parser.add_argument("--insert_at", type=int, required=True, help="Residue number to insert at")
    
    parser.add_argument("--rf_path", default="../Tools/RFdiffusion")
    parser.add_argument("--num_designs", type=int, default=2)
    parser.add_argument("--dry_run", action="store_true")
    
    args = parser.parse_args()
    run_graft(args)

if __name__ == "__main__":
    main()
