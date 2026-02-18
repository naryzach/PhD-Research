import os
import glob
import pandas as pd
import sys
import argparse
import warnings
import numpy as np

# Suppress warnings from torch/gradio/etc
warnings.filterwarnings("ignore")

# Import AllMetal3D
try:
    from allmetal3d.utils.main import predict_cli
except ImportError:
    print("Error: Could not import allmetal3d. Please run this script within the allmetal3d conda environment.")
    sys.exit(1)

# Labels as defined in allmetal3d/utils/main.py
LABELS = ['Alkali', 'MG', 'CA', 'ZN', 'NonZNTM', 'NoMetal']

def main():
    parser = argparse.ArgumentParser(description="Batch run AllMetal3D and summarize results with detailed probabilities.")
    parser.add_argument("--input_dir", default="../Local/AllMetal_input", help="Directory containing input PDB files")
    parser.add_argument("--output_dir", default="../Local/Allmetal_output", help="Directory to save outputs")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing runs")
    
    parser.add_argument("--insert_metal", default=None, help="Metal to insert (e.g., 'ZN') or 'auto'. Defaults to None (no insertion).")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_dir):
        print(f"Error: Input directory {args.input_dir} does not exist.")
        return

    os.makedirs(args.output_dir, exist_ok=True)
    if args.insert_metal:
         holo_dir_likely = os.path.join(args.output_dir, "holo_structures")
         holo_dir_unlikely = os.path.join(args.output_dir, "holo_structures_unlikely")
         os.makedirs(holo_dir_likely, exist_ok=True)
         os.makedirs(holo_dir_unlikely, exist_ok=True)
    
    pdb_files = glob.glob(os.path.join(args.input_dir, "*.pdb")) + glob.glob(os.path.join(args.input_dir, "*.cif"))
    if not pdb_files:
        print(f"No PDB or CIF files found in {args.input_dir}")
        return
        
    print(f"Found {len(pdb_files)} structure files to process.")
    
    summary_data = []

    for pdb_file in pdb_files:
        basename = os.path.splitext(os.path.basename(pdb_file))[0]
        pdb_out_dir = os.path.join(args.output_dir, basename)
        
        # Check if already run (check for _metals.pdb)
        metals_pdb = os.path.join(pdb_out_dir, f"{basename}_metals.pdb")
        
        csv_out_path = os.path.join(pdb_out_dir, f"{basename}_detailed.csv")
        
        if os.path.exists(csv_out_path) and not args.overwrite:
            print(f"Skipping prediction for {basename} (Detailed output exists)")
            try:
                # Load CSV for summary
                df = pd.read_csv(csv_out_path)
                records = df.to_dict('records')
                summary_data.extend(records)
                
                # If we need to insert metal, we have to find the best site from these records
                if args.insert_metal:
                    # Find best site (highest Location_Prob)
                    best_site = max(records, key=lambda x: x['Location_Prob'])
                    
                    # Determine Likely/Unlikely
                    # Likely: Location_Prob >= 50.0 AND Top_Identity != 'NoMetal'
                    is_likely = (best_site['Location_Prob'] >= 50.0) and (best_site['Top_Identity'] != 'NoMetal')
                    target_dir = holo_dir_likely if is_likely else holo_dir_unlikely

                    # We need coordinates. Read them from metals_pdb
                    # metals_pdb line: HETATM    1  X     X X  1 ...
                    # Site_Index in CSV should match atom serial number or order?
                    # In my previous code: Site_Index comes from `res["index"]`.
                    # In `write_probefile`, `i+1` is used for serial. So Site_Index 1 = Atom 1.
                    
                    coord_line = None
                    if os.path.exists(metals_pdb):
                        with open(metals_pdb, 'r') as f:
                            for line in f:
                                if line.startswith("HETATM") or line.startswith("ATOM"):
                                    try:
                                        serial = int(line[6:11].strip())
                                        if serial == best_site['Site_Index']:
                                            coord_line = line
                                            break
                                    except:
                                        pass
                    
                    if coord_line:
                        # Extract coordinates
                        x = float(coord_line[30:38])
                        y = float(coord_line[38:46])
                        z = float(coord_line[46:54])
                        
                        # Determine metal identity
                        if args.insert_metal.lower() == 'auto':
                            metal_id = best_site.get('Top_Identity', 'ZN')
                            if metal_id == 'NoMetal': metal_id = 'ZN' # Fallback
                            if metal_id == 'Alkali': metal_id = 'NA' # Default for alkali
                            if metal_id == 'NonZNTM': metal_id = 'CU' # Default
                        else:
                            metal_id = args.insert_metal.upper()
                            
                        # Insert Metal
                        insert_metal_into_structure(pdb_file, target_dir, basename, x, y, z, metal_id)

            except Exception as e:
                print(f"Error handling existing files for {basename}: {e}")
            continue

        print(f"Processing {basename}...")
        os.makedirs(pdb_out_dir, exist_ok=True)
        
        try:
            # Predict
            ret = predict_cli(
                pdb=pdb_file, 
                models="allmetal3d", 
                output_dir=pdb_out_dir
            )
            
            results = ret[4]
            
            if not results:
                print(f"  No sites found for {basename}.")
                continue

            row_list = []
            best_res = None
            max_p = -1.0
            best_identity = "NoMetal"
            
            for res in results:
                probs = res["probabilities_identity"]
                p_loc = res["location_confidence"]
                
                # Identify the "Top Identity" which matches the max probability
                max_prob_idx = max(range(len(probs)), key=probs.__getitem__)
                top_identity = LABELS[max_prob_idx]

                # Tracking best for insertion
                if p_loc > max_p:
                    max_p = p_loc
                    best_res = res
                    best_identity = top_identity

                row = {
                    "Structure": basename,
                    "Site_Index": res["index"],
                    "Location_Prob": p_loc,
                    "File": metals_pdb 
                }
                
                for label, p in zip(LABELS, probs):
                    row[label] = p
                
                row["Top_Identity"] = top_identity
                
                row_list.append(row)
                summary_data.append(row)
            
            if row_list:
                df_detailed = pd.DataFrame(row_list)
                df_detailed.to_csv(csv_out_path, index=False)
                print(f"  Saved detailed results to {csv_out_path}")
                
            # Metal Insertion (Fresh Run)
            if args.insert_metal and best_res:
                # Determine Likely/Unlikely
                # Likely: Location_Prob >= 50.0 AND Top_Identity != 'NoMetal'
                # Note: best_res has 'location_confidence' and we derived 'best_identity'
                is_likely = (max_p >= 50.0) and (best_identity != 'NoMetal')
                target_dir = holo_dir_likely if is_likely else holo_dir_unlikely

                # We need coordinates. `results` object DOES NOT contain coordinates directly?
                # Wait, `res` has "close_residues" but coordinate data might be lost?
                # Actually, predict_cli writes the probefile. We can read coordinates from there as before.
                
                coord_line = None
                if os.path.exists(metals_pdb):
                     with open(metals_pdb, 'r') as f:
                        for line in f:
                            if line.startswith("HETATM") or line.startswith("ATOM"):
                                try:
                                    serial = int(line[6:11].strip())
                                    if serial == best_res['index']:
                                        coord_line = line
                                        break
                                except:
                                    pass
                
                if coord_line:
                    x = float(coord_line[30:38])
                    y = float(coord_line[38:46])
                    z = float(coord_line[46:54])
                    
                    if args.insert_metal.lower() == 'auto':
                         # Use best_res probabilities
                         probs = best_res["probabilities_identity"]
                         max_prob_idx = max(range(len(probs)), key=probs.__getitem__)
                         metal_id = LABELS[max_prob_idx]
                         
                         if metal_id == 'NoMetal': metal_id = 'ZN'
                         if metal_id == 'Alkali': metal_id = 'NA'
                         if metal_id == 'NonZNTM': metal_id = 'CU'
                    else:
                        metal_id = args.insert_metal.upper()
                        
                    insert_metal_into_structure(pdb_file, target_dir, basename, x, y, z, metal_id)

        except Exception as e:
            print(f"  Error processing {basename}: {e}")
            import traceback
            traceback.print_exc()

    # Generate Summary
    if summary_data:
        df = pd.DataFrame(summary_data)
        df = df.sort_values(by=["Structure", "Location_Prob"], ascending=[True, False])
        cols = ["Structure", "Site_Index", "Location_Prob", "Top_Identity"] + LABELS + ["File"]
        cols = [c for c in cols if c in df.columns]
        df = df[cols]
        output_csv = os.path.join(args.output_dir, "binding_summary.csv")
        df.to_csv(output_csv, index=False)
        print(f"\nSummary saved to {output_csv}")
        print(df.head())
    else:
        print("No data extracted for summary.")

def insert_metal_into_structure(original_path, output_dir, basename, x, y, z, metal_id):
    """
    Reads original PDB/CIF and writes a new PDB with the inserted metal.
    Converts CIF to PDB format if input is CIF (simple atom extraction).
    """
    out_name = f"{basename}_holo_{metal_id}.pdb"
    out_path = os.path.join(output_dir, out_name)
    
    try:
        atom_lines = []
        is_cif = original_path.endswith(".cif")
        
        with open(original_path, 'r') as f:
            lines = f.readlines()
            
        if is_cif:
            # Very basic CIF to PDB converter for ATOM/HETATM sites
            # Assumes loop_ structure for _atom_site.
            # This is fragile but works for RFdiffusion/standard CIFs usually.
            # Ideally use Bio.PDB but keeping it dependency-light/fast if possible.
            # Or just filter for "ATOM" lines if they start with it (PDB-style CIF?)
            # No, CIF uses `_atom_site.group_PDB` etc. 
            # Let's try to use Bio.PDB for robustness if we have it.
            # allmetal3d env has Bio.PDB.
            
            from Bio.PDB.MMCIFParser import MMCIFParser
            from Bio.PDB import PDBIO
            from Bio.PDB.Atom import Atom
            
            parser = MMCIFParser()
            structure = parser.get_structure("struct", original_path)
            
        else:
            from Bio.PDB import PDBParser, PDBIO
            from Bio.PDB.Atom import Atom
            parser = PDBParser()
            structure = parser.get_structure("struct", original_path)

        # Create new Metal Atom
        # Chain X, Residue 1
        model = structure[0]
        
        # Create a new chain for the metal if possible, or append to existing?
        # RFdiffusion usually expects single chain or specific format?
        # Let's append to a new Chain 'M' for Metal
        from Bio.PDB.Chain import Chain
        from Bio.PDB.Residue import Residue
        
        metal_chain_id = 'M'
        # Check if M exists
        if metal_chain_id in model:
            metal_chain_id = 'Z' # Fallback
            
        chain = Chain(metal_chain_id)
        residue = Residue(('H_M', 1, ' '), metal_id, ' ')
        
        # Atom name, coord, B-factor, Occupancy, Altloc, Fullname, Serial
        # Atom(name, coord, bfactor, occupancy, altloc, fullname, serial_number, element=None)
        
        new_atom = Atom(
            name=metal_id,
            coord=np.array([x, y, z], dtype='f'),
            bfactor=0.0,
            occupancy=1.0,
            altloc=' ',
            fullname=f"{metal_id:>4}",
            serial_number=9999,
            element=metal_id
        )
        
        residue.add(new_atom)
        chain.add(residue)
        model.add(chain)
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(out_path)
        print(f"  Generated holo structure: {out_path}")
        
    except Exception as e:
        print(f"  Error generating holo structure for {basename}: {e}")
        # Fallback to simple text append if Bio.PDB fails or for PDBs
        if original_path.endswith(".pdb"):
             with open(original_path, 'r') as f:
                content = f.read()
             # Remove END
             lines = [l for l in content.splitlines() if not l.startswith("END")]
             
             # HETATM 9999 ZN    ZN M   1      12.345  67.890  12.345  1.00  0.00          ZN
             het_line = f"HETATM 9999 {metal_id:>2}    {metal_id:>3} M   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {metal_id:>2}  "
             lines.append(het_line)
             lines.append("END")
             
             with open(out_path, 'w') as f:
                 f.write('\n'.join(lines))
             print(f"  Generated holo structure (text fallback): {out_path}")

if __name__ == "__main__":
    main()

