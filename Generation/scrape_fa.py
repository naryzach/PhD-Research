import os
import re
import pandas as pd
from collections import Counter

# === CONFIG ===
input_dir = "../Local/ProteinMPNN_Designs/FA_files" 
output_dir = os.path.join(input_dir, "summaries")
os.makedirs(output_dir, exist_ok=True)

# Regex patterns
filename_pattern = re.compile(r"ProteinMPNN_(.+?)_versus_(.+?)_(.+?)_Sequences\.fa")
temp_pattern = re.compile(r"T\s*=\s*([0-9.]+)")
score_pattern = re.compile(r"score\s*=\s*([0-9.]+)")


# === HELPER FUNCTIONS ===
def parse_fasta(filepath):
    """Parse FASTA file -> list of (header, sequence)."""
    entries = []
    with open(filepath, "r") as f:
        header, seq = None, []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header and seq:
                    entries.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header and seq:
            entries.append((header, "".join(seq)))
    return entries


def extract_loop_sequence(seq, loop):
    """Extract loop region depending on loop label."""
    loop_seq = seq[:10]  # default
    if loop.upper() == "AB":
        loop_seq = seq[30:36]
    elif loop.upper() == "C":
        loop_seq = seq[62:68]
    elif loop.upper() == "EF":
        loop_seq = seq[92:96]
    elif loop.upper() == "GH":
        loop_seq = seq[127:137]
    elif loop.lower() == "multi":
        loop_seq = seq[143:153]
    return loop_seq


def summarize_fasta(filepath):
    """Summarize a single FASTA file into a DataFrame."""
    filename = os.path.basename(filepath)
    match = filename_pattern.match(filename)
    if not match:
        print(f"Skipping {filename}: filename pattern not recognized.")
        return None

    ligand, target, loop = match.groups()
    entries = parse_fasta(filepath)

    # Extract metadata from headers
    seq_info = []
    for header, seq in entries:
        temp_match = temp_pattern.search(header)
        score_match = score_pattern.search(header)
        temperature = float(temp_match.group(1)) if temp_match else None
        score = float(score_match.group(1)) if score_match else None
        seq_info.append((seq, temperature, score))

    # Count unique sequences
    seq_counts = Counter([seq for seq, _, _ in seq_info])

    data = []
    for seq, count in seq_counts.items():
        temps = [t for s, t, _ in seq_info if s == seq and t is not None]
        scores = [sc for s, _, sc in seq_info if s == seq and sc is not None]
        mean_temp = round(sum(temps) / len(temps), 3) if temps else None
        mean_score = round(sum(scores) / len(scores), 3) if scores else None
        loop_seq = extract_loop_sequence(seq, loop)

        data.append({
            "File": filename,
            "Ligand": ligand,
            "Target": target,
            "Loop": loop,
            "Loop_Sequence": loop_seq,
            "Occurrences": count,
            "Temperature": mean_temp,
            "Score": mean_score,
            "Sequence": seq
        })

    df = pd.DataFrame(data)
    df = df.sort_values(by="Occurrences", ascending=False)

    # Save per-file CSV
    csv_path = os.path.join(output_dir, filename.replace(".fa", "_summary.csv"))
    df.to_csv(csv_path, index=False)

    # Save compact LaTeX table (without full sequence)
    latex_path = os.path.join(output_dir, filename.replace(".fa", "_summary.tex"))
    latex_cols = ["Ligand", "Target", "Loop", "Loop_Sequence", "Occurrences", "Score"]
    with open(latex_path, "w") as f:
        f.write(df[latex_cols].to_latex(index=False, escape=False))

    print(f"Processed {filename}: {len(df)} unique sequences.")
    return df


# === MAIN ===
all_dfs = []
for file in os.listdir(input_dir):
    if file.endswith(".fa"):
        df = summarize_fasta(os.path.join(input_dir, file))
        if df is not None:
            all_dfs.append(df)

# === COMBINE ALL FILES ===
if all_dfs:
    combined_df = pd.concat(all_dfs, ignore_index=True)
    combined_df = combined_df.sort_values(by=["Ligand", "Target", "Loop", "Occurrences"], ascending=[True, True, True, False])

    combined_csv = os.path.join(output_dir, "Combined_ProteinMPNN_Summary.csv")
    combined_tex = os.path.join(output_dir, "Combined_ProteinMPNN_Summary.tex")

    combined_df.to_csv(combined_csv, index=False)
    
    latex_cols = ["Ligand", "Target", "Loop", "Loop_Sequence", "Occurrences", "Score"]
    with open(combined_tex, "w") as f:
        f.write(combined_df[latex_cols].to_latex(index=False, escape=False))

    print(f"\nCombined summary saved:")
    print(f"  --> {combined_csv}")
    print(f"  --> {combined_tex}")
else:
    print("No valid .fa files found.")
