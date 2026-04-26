import os, time, subprocess, json, re, shutil

# 3-letter to 1-letter amino acid code map
aa_map = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
        'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def get_aa_sequence(pdb_file, chain_letter):
    current_sequence = []
    with open(pdb_file, 'r') as f:
        for line in f:
            line_items = line.split()
            if line.startswith('ATOM') and line_items[4] == chain_letter and line_items[2] == "CA":
                res_name = line_items[3]
                if res_name in aa_map:
                    current_sequence.append(aa_map[res_name])
    return "".join(current_sequence)

original_directory = os.getcwd()

data_path = "../Data"
rf_diff_path = "../Tools/RFdiffusion/scripts"
protein_mpnn_path = "../Tools/ProteinMPNN/"

pmpnn_out_dir = "../Local/proteinmpnn_output"
output_prefix = "design"

loop_names = ["C", "EF"]
if isinstance(loop_names, str):
    loop_names = [loop_names]

loop_configs = {
    "AB": {"normal": 6, "max": 15, "pos": 30, "left": "LVK", "right": "LVY"},
    "C": {"normal": 6, "max": 15, "pos": 62, "left": "HTE", "right": "GLK"},
    "EF": {"normal": 4, "max": 10, "pos": 92, "left": "MYT", "right": "FVE"},
    "GH": {"normal": 10, "max": 20, "pos": 127, "left": "KSC", "right": "NEC"},
    "Multi": {"normal": 10, "max": 20, "pos": 143, "left": "LWT", "right": "YQS"}
}

for name in loop_names:
    if name not in loop_configs:
        raise Exception(f"Specified Loop '{name}' is unknown")

selected_loops = [loop_configs[name] for name in loop_names]
selected_loops.sort(key=lambda x: x["pos"])

chain_to_design = "A"
fixed_chains = ["B"]
total_length = 121

contig_parts = []
current_pos = 1

for loop in selected_loops:
    if current_pos <= loop["pos"]:
        contig_parts.append(f"{chain_to_design}{current_pos}-{loop['pos']}")
    contig_parts.append(f"{loop['normal']}-{loop['max']}")
    current_pos = loop["pos"] + loop["normal"] + 1

if current_pos <= total_length:
    contig_parts.append(f"{chain_to_design}{current_pos}-{total_length}")

contig_string = "/".join(contig_parts)
num_sequences_to_generate = 25
mp_seqs = 25

print(original_directory)
print(contig_string)


# %% [markdown]
# ### Run RFdiffusion

'''
pdb_file_list = ["PDB_fold_timp3_v_adam10cd_wt_model_0.pdb", 
             "PDB_fold_timp3_v_mmp10cd_wt_model_0.pdb", 
             "PDB_fold_timp3_variant_adam17cd_wt_model_0.pdb",
             "PDB_fold_timp3_variant_mmp2cd_wt_model_0.pdb",
             "PDB_fold_timp3_variant_mmp9cd_wt_model_0.pdb"
             ]
'''
pdb_file_list = [#"TIMP3_vs_ADAM10_HADDOCK_Xray.pdb",
             "TIMP3_vs_ADAM17_HADDOCK_Xray.pdb",
             #"TIMP3_vs_MMP2_HADDOCK_Xray.pdb",
             #"TIMP3_vs_MMP7_HADDOCK_Xray.pdb",
             #"TIMP3_vs_MMP9_HADDOCK_Xray.pdb",
             #"TIMP3_vs_MMP10_HADDOCK_Xray.pdb"
             ]


for pdb_complex_file_name in pdb_file_list:
    pdb_path = os.path.join(data_path, pdb_complex_file_name)
    output_dir = os.path.join("../Local/rfdiffusion_output", pdb_complex_file_name[:-4])

    fix_chain_len = len(get_aa_sequence(pdb_path, fixed_chains[0]))
    print(fix_chain_len)

    if not pdb_path:
        print("Cannot run RFdiffusion without a scaffold PDB file.")
        raise Exception("No PDB File")

    print("Preparing to run RFdiffusion...")

    # Construct the command for RFdiffusion
    run_command = [
        "python",
        os.path.join(rf_diff_path.replace('../', ''), "run_inference.py"),
        f"inference.output_prefix={os.path.join(output_dir.replace('../', ''), output_prefix)}",
        f"inference.input_pdb={pdb_path.replace('../', '')}",
        f'contigmap.contigs=[{contig_string}/0 {fixed_chains[0]}1-{fix_chain_len}]', #-- /0 is chain break
        "diffuser.T=20", # default 50, short sequence
        f"inference.num_designs={num_sequences_to_generate}",
    ]

    print("Running RFdiffusion to generate novel loops and structures...")
    print(" ".join(run_command))

    # Run the command and stream output live
    os.chdir("..") 
    st = time.time()
    result = subprocess.run(run_command, capture_output=True, text=True)
    end = time.time()
    os.chdir(original_directory)

    print("--- STDOUT ---")
    print(result.stdout)

    print("--- STDERR ---")
    print(result.stderr)

    print(f"RFdiffusion finished in {(end-st)/60:.2f} minutes.")



for pdb_complex_file_name in pdb_file_list:
    pdb_path = os.path.join(data_path, pdb_complex_file_name)
    output_dir = os.path.join("../Local/rfdiffusion_output", pdb_complex_file_name[:-4])

    fix_chain_len = len(get_aa_sequence(pdb_path, fixed_chains[0]))
    print(fix_chain_len)

    if not pdb_path:
        print("Cannot run RFdiffusion without a scaffold PDB file.")
        raise Exception("No PDB File")
    # %%
    # --- Parse PDB Results to Extract Sequences ---
    print("parsing generated PDBs to extract sequences...")

    generated_sequences = set()
    generated_loops = set()
    chain_id_to_extract = contig_string.split('/')[0][0] # Get chain from contig

    for i in range(num_sequences_to_generate):
        pdb_file_name = f"{output_prefix}_{i}.pdb"
        pdb_file = os.path.join(output_dir, pdb_file_name)
        if os.path.exists(pdb_file):
            current_sequence = get_aa_sequence(pdb_file, chain_id_to_extract)
            generated_sequences.add(current_sequence)
            
            loop_seqs = []
            curr_idx = 0
            for loop in selected_loops:
                flank_left = loop["left"]
                flank_right = loop["right"]
                regex_pattern = re.compile(f"{flank_left}([A-Z]+?){flank_right}")
                match = regex_pattern.search(current_sequence[curr_idx:])
                if match:
                    loop_seqs.append(match.group(1))
                    curr_idx += match.end() - len(flank_right)
                else:
                    loop_seqs.append("MISSING")
            
            if "MISSING" not in loop_seqs:
                generated_loops.add("_".join(loop_seqs))

    print(f"\nExtracted {len(generated_sequences)} unique, novel loop sequences.")
    print("Here are a few examples (full):")
    for seq in list(generated_sequences)[:5]:
        print(f"   - {seq}")


    # %% [markdown]
    # ### ProteinMPNN

    # %%
    pmpnn_out_dir_full = os.path.join(pmpnn_out_dir, pdb_complex_file_name[:-4])
    os.makedirs(pmpnn_out_dir_full, exist_ok=True)

    print(pmpnn_out_dir_full)

    for i in range(num_sequences_to_generate):
        # Create JSONL helper files for ProteinMPNN
        chain_json = {f"{output_prefix}_{i}": [[c for c in fixed_chains], [str(chain_to_design)]]}
        with open(f"{pmpnn_out_dir_full}/{output_prefix}_chain_B.jsonl", "w") as f:
            json.dump(chain_json, f)

        pdb_file_name = f"{output_prefix}_{i}.pdb"
        pdb_file = os.path.join(output_dir, pdb_file_name)
        if os.path.exists(pdb_file):
            aa_sequence = get_aa_sequence(pdb_file, chain_to_design)
        else:
            raise FileNotFoundError
        print(aa_sequence)

        fixed_positions = []
        current_fixed_start = 1
        current_seq_idx = 0
        
        for loop in selected_loops:
            flank_left = loop["left"]
            flank_right = loop["right"]
            regex_pattern = re.compile(f"{flank_left}([A-Z]+?){flank_right}")
            match = regex_pattern.search(aa_sequence[current_seq_idx:])
            if match:
                match_start = current_seq_idx + match.start()
                match_end = current_seq_idx + match.end()
                inserted_seq = match.group(1)
                loop_start_1idx = match_start + len(flank_left) + 1
                
                fixed_positions.extend(range(current_fixed_start, loop_start_1idx))
                current_fixed_start = loop_start_1idx + len(inserted_seq)
                current_seq_idx = match_end - len(flank_right)
        
        new_total_length = len(aa_sequence)
        fixed_positions.extend(range(current_fixed_start, new_total_length + 1))

        # Freeze all residues except the inpainted region
        fixed_json = {
            f"{output_prefix}_{i}": {
                str(chain_to_design): fixed_positions
            }
        }
        for c in fixed_chains:
            fixed_json[f"{output_prefix}_{i}"][c] = [i for i in range(1,fix_chain_len+1)]
        with open(f"{pmpnn_out_dir_full}/fixed.jsonl", "w") as f:
            json.dump(fixed_json, f)

        if not output_dir:
            print("Cannot run ProteinMPNN without a scaffold files.")
            raise Exception("No RFdiffusion Files")

        print("Preparing to run ProteinMPNN...")

        # Construct the command for RFdiffusion
        run_command = [
            "python",
            os.path.join(protein_mpnn_path.replace('../', ''), "protein_mpnn_run.py"),
            "--pdb_path", pdb_file.replace('../', ''),
            "--out_folder", f"{pmpnn_out_dir_full.replace('../', '')}",
            "--chain_id_jsonl", f"{pmpnn_out_dir_full.replace('../', '')}/{output_prefix}_chain_B.jsonl",
            "--fixed_positions_jsonl", f"{pmpnn_out_dir_full.replace('../', '')}/fixed.jsonl",
            "--num_seq_per_target", f"{mp_seqs}",
            "--sampling_temp", "0.1"
        ]

        print("Running ProteinMPNN to generate sequences...")
        print(" ".join(run_command))

        # Run the command and stream output live
        os.chdir("..") 
        st = time.time()
        result = subprocess.run(run_command, capture_output=True, text=True)
        end = time.time()
        os.chdir(original_directory)

        print("--- STDOUT ---")
        print(result.stdout)

        print("--- STDERR ---")
        print(result.stderr)

        print(f"ProteinMPNN finished in {(end-st)/60:.2f} minutes.")

        os.makedirs(os.path.join(pmpnn_out_dir_full, "FA_files"), exist_ok=True)
        shutil.move(os.path.join(pmpnn_out_dir_full, "seqs", f"{output_prefix}_{i}.fa"), 
                    os.path.join(pmpnn_out_dir_full, "FA_files", f"ProMPNN_from_RFd_{i}.fa"))



    # %% [markdown]
    # ### Analyze ProteinMPNN

    # %%
    import re
    import csv
    from pathlib import Path
    import pandas as pd

    # ---- USER SETTINGS ----
    input_dir = Path(pmpnn_out_dir_full, "FA_files")  # directory containing your design files
    output_csv = "design_summary.csv"
    top_n = 5 # number of top unique loops to keep per loop length

    # ---- PARSE FILES ----
    records = []

    for file in input_dir.glob("*.fa"):
        with open(file) as f:
            lines = [line.strip() for line in f if line.strip()]
        
        for i in range(0, len(lines), 2):
            header = lines[i]
            seq = lines[i+1]
            if "/" in seq:
                seq = seq.split("/")[0]
            
            # extract metadata
            entry = {"file": file.name, "full_seq": seq}
            parts = header.strip(">").split(", ")
            for p in parts:
                if "=" in p:
                    key, val = p.split("=", 1)
                    entry[key.strip()] = val.strip()
                else:
                    entry["design_type"] = p.strip()
            
            # extract loop region between flanks
            # extract loop region between flanks
            loop_data = {}
            curr_idx = 0
            
            for name_idx, loop in enumerate(selected_loops):
                loop_name = loop_names[name_idx]
                f_left = loop["left"]
                f_right = loop["right"]
                m = re.search(f"{f_left}(.*?){f_right}", seq[curr_idx:])
                
                if m:
                    extracted_seq = m.group(1)
                    loop_data[f"loop_{loop_name}_seq"] = extracted_seq
                    loop_data[f"loop_{loop_name}_length"] = len(extracted_seq)
                    curr_idx += m.end() - len(f_right)
                else:
                    loop_data[f"loop_{loop_name}_seq"] = "MISSING"
                    loop_data[f"loop_{loop_name}_length"] = 0
                    
            entry.update(loop_data)
            records.append(entry)

    # ---- MAKE DATAFRAME ----
    df = pd.DataFrame(records)

    # Convert numeric columns
    for col in ["score", "global_score", "seq_recovery", "T"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # ---- REORDER COLUMNS ----
    base_columns = ["file", "sample"]
    loop_columns = []
    for ln in loop_names:
        loop_columns.extend([f"loop_{ln}_seq", f"loop_{ln}_length"])
    
    other_columns = [
        "score", "global_score", "seq_recovery",
        "design_type", "full_seq", "fixed_chains", "designed_chains",
        "model_name", "git_hash", "seed", "T"
    ]
    
    column_order = base_columns + loop_columns + other_columns
    df = df[[c for c in column_order if c in df.columns]]

    # ---- SAVE FULL SUMMARY ----
    df.to_csv(os.path.join(pmpnn_out_dir_full, output_csv), index=False)
    print(f"Saved {len(df)} entries to {output_csv}")

    # ---- SUMMARIZE BEST UNIQUE LOOPS ----
    seq_cols = [c for c in df.columns if c.startswith("loop_") and c.endswith("_seq")]
    len_cols = [c for c in df.columns if c.startswith("loop_") and c.endswith("_length")]
    
    clean_df = df.dropna(subset=seq_cols + ["score"]).copy()

    if not clean_df.empty:
        # Keep best (lowest score) per unique loop sequence array
        best_per_loop = (
            clean_df.sort_values("score")
            .groupby(seq_cols, as_index=False)
            .first()
        )

        # Compute averages per loop length array
        avg_stats = (
            clean_df.groupby(len_cols, as_index=False)
            .agg(
                avg_score=("score", "mean"),
                avg_seq_recovery=("seq_recovery", "mean"),
                count=("score", "count")
            )
            .sort_values(len_cols)
        )

        # Select top N best loops per length configuration
        best_per_length = (
            best_per_loop.sort_values(len_cols + ["score"], ascending=[True]*len(len_cols) + [True])
            .groupby(len_cols, group_keys=False)
            .head(top_n)
        )

        # Merge in average stats
        best_per_length = best_per_length.merge(avg_stats, on=len_cols, how="left")

        # Sort final summary
        best_per_length = best_per_length.sort_values(len_cols + ["score"], ascending=[True]*len(len_cols) + [True])

        # ---- SAVE SUMMARIES ----
        best_per_length.to_csv(os.path.join(pmpnn_out_dir_full, "best_loops_per_length.csv"), index=False)
        avg_stats.to_csv(os.path.join(pmpnn_out_dir_full, "loop_length_averages.csv"), index=False)

        print(f"Saved best {top_n} unique loops per length to best_loops_per_length.csv")
        print("Saved loop length averages to loop_length_averages.csv")
    else:
        print("No valid scores generated, summaries skipped.")
