# %% [markdown]
# # 🧬 ProteinMPNN Input Generator
# This notebook creates properly formatted JSONL files for partial redesign.
# %%
import os
import json
from Bio.PDB import PDBParser
import subprocess
import shutil
import re
import glob
from itertools import combinations
# %% [markdown]
# ## Configuration

'''
pdb_files = ["PDB_fold_timp3_v_adam10cd_wt_model_1.pdb", 
             "PDB_fold_timp3_v_mmp10cd_wt_model_1.pdb", 
             "PDB_fold_timp3_variant_adam17cd_wt_model_1.pdb",
             "PDB_fold_timp3_variant_mmp2cd_wt_model_1.pdb",
             "PDB_fold_timp3_variant_mmp9cd_wt_model_1.pdb"]'''
pdb_files = ["TIMP3_vs_ADAM10_HADDOCK_Xray.pdb",
             "TIMP3_vs_ADAM17_HADDOCK_Xray.pdb",
             "TIMP3_vs_MMP2_HADDOCK_Xray.pdb",
             "TIMP3_vs_MMP7_HADDOCK_Xray.pdb",
             "TIMP3_vs_MMP9_HADDOCK_Xray.pdb",
             "TIMP3_vs_MMP10_HADDOCK_Xray.pdb"]
loops_options = ["AB", "C", "EF", "GH", "Multi"]  # Options: "AB", "C", "EF", "GH", "Multi"

for pdb_name in pdb_files:
    print(f"Processing PDB file: {pdb_name}")
    for loop_name in loops_options:
        print(f"  Designing loop: {loop_name}")
        # %%
        # --- Path to ProteinMPNN directory ---
        proteinmpnn_repo = "../Tools/ProteinMPNN"

        # --- Input PDB file ---
        data_path = "../Data"
        pdb_path_orig = os.path.join(data_path, pdb_name)
        pdb_basename = os.path.basename(pdb_path_orig).replace(".pdb", "")

        # --- Output directory (will be created if not exists) ---
        output_dir = "../Local/ProteinMPNN_Designs"
        out_parse = os.path.join(output_dir, "parsed_chains.jsonl")
        out_assign = os.path.join(output_dir, "assigned_chains.jsonl")
        out_fixed_dict = os.path.join(output_dir, "fixed_positions.jsonl")

        # --- Move PDB to a solo directory ---
        solo_pdb_path = os.path.join(output_dir, "solo_pdb")
        os.makedirs(solo_pdb_path, exist_ok=True)
        shutil.copy(pdb_path_orig, solo_pdb_path)
        pdb_path = os.path.join(solo_pdb_path, os.path.basename(pdb_path_orig))
        print(f"PDB path: {pdb_path}")

        # --- Chains and residues to redesign ---
        chains_to_design = ["A"]
        loops_to_design = [loop_name]  # Options: "AB", "C", "EF", "GH", "Multi"
        residues_to_redesign = []
        if "AB" in loops_to_design:
            residues_to_redesign += [i for i in range(31, 37)]
        if "C" in loops_to_design:
            residues_to_redesign += [i for i in range(63, 69)]
        if "EF" in loops_to_design:
            residues_to_redesign += [i for i in range(93, 97)]
        if "GH" in loops_to_design:
            residues_to_redesign += [i for i in range(128, 138)]
        if "Multi" in loops_to_design:
            residues_to_redesign += [i for i in range(144, 154)]
        if loops_to_design == None:
            residues_to_redesign += [i for i in range(1, 188)]

        # --- Number of designed sequences to generate ---
        num_seq_per_target = 1000

        # --- Sampling temperature (0.1–0.2 is typical) ---
        sampling_temp = 0.2

        # --- Helper Scripts ---
        help_scritps = "../Tools/ProteinMPNN/helper_scripts"
        help_parse = os.path.join(help_scritps, "parse_multiple_chains.py")
        help_assign_chain = os.path.join(help_scritps, "assign_fixed_chains.py")
        help_fixed_dict = os.path.join(help_scritps, "make_fixed_positions_dict.py")


        original_directory = os.getcwd() 

        name_parse = pdb_name.split("_")
        target = name_parse[2]
        ligand = name_parse[0]

        # %% [markdown]
        # ## Setup

        # %%
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Derive output prefix from PDB filename (no extension)
        output_prefix = os.path.splitext(os.path.basename(pdb_path))[0]

        # Load structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("input", pdb_path)

        # Gather residue numbers for each chain
        chain_residues = {}
        for model in structure:
            for chain in model:
                chain_id = chain.id
                resnums = [res.id[1] for res in chain if res.id[0] == " "]
                chain_residues[chain_id] = resnums

        chains_to_fix = [c for c in chain_residues.keys() if c not in chains_to_design]
        print("Chains found:", list(chain_residues.keys()))

        # %%
        parse_command = [
            "python",
            help_parse.replace("../", ""),
            "--input_path", solo_pdb_path.replace("../", ""),
            "--output_path", out_parse.replace("../", "")
        ]

        print(" ".join(parse_command))

        os.chdir("..")
        result = subprocess.run(parse_command, capture_output=True, text=True)
        os.chdir(original_directory)

        print(result.stdout)
        print(result.stderr)

        # %%
        parsed_chains = {}

        AA_lookup = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
            'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
            'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
            'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }
        for model in structure:
            for chain in model:
                tmp_chain = []
                for res in chain:
                    tmp_chain.append(AA_lookup.get(res.resname, 'X'))  # 'X' for unknown residues
                parsed_chains[chain.id] = ''.join(tmp_chain)
                print(f"Chain ID: {chain.id}")
                print(f"Sequence: {parsed_chains[chain.id]}")

        # %% [markdown]
        # 
        # ### Build fixed_positions dictionary
        # For each chain, fix everything except residues_to_redesign (for chosen chains)
        # For other chains, fix all residues (i.e. no redesign)

        # %%
        fixed_positions = {output_prefix: {}}

        for chain_id, resnums in chain_residues.items():
            if chain_id in chains_to_design:
                fixed_positions[output_prefix][chain_id] = [
                    r for r in resnums if r not in residues_to_redesign
                ]
            else:
                fixed_positions[output_prefix][chain_id] = []

        # Write fixed_positions.jsonl
        fixed_positions_path = os.path.join(output_dir, f"{output_prefix}_fixed_positions.jsonl")
        with open(fixed_positions_path, "w") as f:
            f.write(json.dumps(fixed_positions))
        print(f"Saved fixed positions file: {fixed_positions_path}")

        # %% [markdown]
        # ### Create chain_id JSONL

        # %%
        assign_command = [
            "python",
            help_assign_chain.replace("../", ""),
            "--input_path", out_parse.replace("../", ""),
            "--output_path", out_assign.replace("../", ""),
            "--chain_list", str(chains_to_design).replace("[", "").replace("]", "").replace("'", "")
        ]

        print(" ".join(assign_command))

        os.chdir("..")
        result = subprocess.run(assign_command, capture_output=True, text=True)
        os.chdir(original_directory)

        print(result.stdout)
        print(result.stderr)

        chain_id_path = out_assign

        # %% [markdown]
        # ### Setup Summary

        # %%
        print("\n=== CONFIGURATION SUMMARY ===")
        print(f"PDB path: {pdb_path}")
        print(f"ProteinMPNN repo: {proteinmpnn_repo}")
        print(f"Output directory: {output_dir}")
        print(f"Chains found: {list(chain_residues.keys())}")
        print(f"Design chains: {chains_to_design}")
        print(f"Design residues: {residues_to_redesign}")
        print(f"Fixed positions file: {fixed_positions_path}")
        print(f"Chain IDs file: {chain_id_path}")

        # %% [markdown]
        # ## Run ProteinMPNN

        # %%
        """ This section is not implemented currently.
        fix_command = [
            "python",
            help_fixed_dict.replace("../", ""),
            "--input_path", out_parse.replace("../", ""),
            "--output_path", out_fixed_dict.replace("../", ""),
            "--chain_list", str(chains_to_design).replace("[", "").replace("]", "").replace("'", ""),
            "--position_list", str(residues_to_redesign).replace("[", "").replace("]", "").replace("'", "").replace(",", "")
        ]

        print(" ".join(fix_command))

        os.chdir("..")
        result = subprocess.run(fix_command, capture_output=True, text=True)
        os.chdir(original_directory)

        print(result.stdout)
        print(result.stderr)

        fixed_positions_path = out_fixed_dict
        """

        # %%
        # Construct the full command
        run_command = [
            "python",
            os.path.join(proteinmpnn_repo.replace("../", ""), "protein_mpnn_run.py"),
            "--pdb_path", pdb_path.replace("../", ""),
            #"--jsonl_path", out_parse.replace("../", ""),
            "--out_folder", output_dir.replace("../", ""),
            "--chain_id_jsonl", chain_id_path.replace("../", ""),
            "--fixed_positions_jsonl", fixed_positions_path.replace("../", ""),
            "--num_seq_per_target", str(num_seq_per_target),
            "--sampling_temp", str(sampling_temp),
        ]

        print("\n=== RUNNING PROTEINMPNN ===")
        print(" ".join(run_command))

        # Run the command and stream output live
        os.chdir("..") 
        result = subprocess.run(run_command, capture_output=True, text=True)
        os.chdir(original_directory)

        print("\n=== STDOUT ===")
        print(result.stdout)

        print("\n=== STDERR ===")
        print(result.stderr)

        print("\nProteinMPNN run complete. Check your output directory for results.")

        # %% [markdown]
        # ## Analyze the Data

        # %%

        # NOTE: This only allows for one fixed chain. Modify as needed for multiple fixed chains.

        fixed_chain_id = chains_to_fix[0]  # Assuming only one fixed chain for simplicity
        output_json_file = output_dir + f"/alphafold_batch_{ligand}_versus_{target}_{"".join(loops_to_design)}_ProteinMPNN.json"

        print(f"Starting post-processing...")

        fixed_seq = parsed_chains[fixed_chain_id]
        print(f"Found fixed chain sequence for Chain '{fixed_chain_id}' (Length: {len(fixed_seq)}).")


        seqs_dir = os.path.join(output_dir, 'seqs')
        fasta_files = glob.glob(os.path.join(seqs_dir, '*.fa*')) # Finds .fa, .fasta, etc.

        if not fasta_files:
            print(f"!! ERROR: No FASTA files found in directory: {seqs_dir}")
            print(f"!!        Please check your 'output_dir' variable.")
            raise FileNotFoundError(f"No FASTA file in {seqs_dir}")

        # Use the first FASTA file found in the directory
        mpnn_output_file = seqs_dir + f"/{pdb_basename}.fa"
        print(f"Reading designs from: {mpnn_output_file}")

        with open(mpnn_output_file, 'r') as f:
            mpnn_output_data = f.read()

        parsed_sequences = []

        all_entries = mpnn_output_data.strip().split('\n>')
        designed_entries = all_entries[1:] 

        if not designed_entries:
            print(f"!! ERROR: No designed sequences found in the file. Only the native sequence was present.")
            raise ValueError("No designed sequences to parse.")

        for entry in designed_entries:
            lines = entry.split('\n')
            if len(lines) < 2:
                continue # Skip any malformed entries
                
            header_line = lines[0].strip()
            sequence = lines[1].strip()
            
            # Extract the score
            score_match = re.search(r"score=([\d\.]+)", header_line)
            
            if score_match:
                score = float(score_match.group(1))
                parsed_sequences.append((score, sequence, header_line))

        # Sort by score (lower is better for NLL/perplexity)
        sorted_designs = sorted(parsed_sequences, key=lambda x: x[0], reverse=False)

        print(f"Parsed {len(sorted_designs)} total designs. Now filtering for unique sequences...")
        unique_designs = []
        seen_sequences = set()

        for (score, seq, header) in sorted_designs:
            if seq not in seen_sequences:
                unique_designs.append((score, seq, header))
                seen_sequences.add(seq)

        # Get the top 10 designs
        print(f"Found {len(unique_designs)} unique designs after filtering.")
        top_10_designs = unique_designs[:10]
        alphafold_jobs = []

        print("\n" + "=" * 70)
        print(f"Top {len(top_10_designs)} Designs (Sorted by score, lower is better):")
        print("-" * 70)
        print(f"{'Rank':<5} | {'Score':<8} | {'Designed Loop':<15} | {'Header Info':<30}")
        print("-" * 70)

        for i, (score, seq, header) in enumerate(top_10_designs):
            rank = i + 1
            
            # Extract the redesigned loop to print for verification
            extracted_loop = "".join([seq[i-1] for i in residues_to_redesign])

            # Create the job name
            job_name = f"{ligand}_versus_{target}_{extracted_loop}_ProteinMPNN_{rank:02d}"

            # Create the job object for AlphaFold
            job_object = {
                "name": job_name,
                "sequences": [
                    { "proteinChain": { "sequence": parsed_chains[fixed_chain_id].replace('X', ''), "count": 1 } }, # Fixed Chain
                    { "proteinChain": { "sequence": seq, "count": 1 } }                      # Designed Chain
                ]
            }
            
            # Add this job object to our main list
            alphafold_jobs.append(job_object)

            print(f"{rank:<5} | {score:<8.4f} | {extracted_loop:<15} | {header[:30]}...")

        # --- Write the final JSON file ---
        with open(output_json_file, 'w') as f:
            json.dump(alphafold_jobs, f, indent=2)

        print("\n" + "=" * 70)
        print(f"Success! Created '{output_json_file}' with {len(top_10_designs)} jobs.")
        print("   This file is ready for batch processing with AlphaFold/ColabFold.")

        final_file_dest = os.path.join(output_dir, "FA_files")
        os.makedirs(final_file_dest, exist_ok=True)
        shutil.move(f"{os.path.join(output_dir, "seqs", os.path.basename(pdb_path.replace(".pdb", "")))}.fa", 
                    f"{os.path.join(final_file_dest, "")}ProteinMPNN_{ligand}_versus_{target}_{'-'.join(loops_to_design)}_Sequences.fa")
        
    os.remove(pdb_path)