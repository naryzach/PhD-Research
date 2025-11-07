import os
import re
import pandas as pd
from collections import Counter

# === CONFIG ===
input_dir = "../Local/ProteinMPNN_Designs/FA_files"  # <-- Change this to your directory path
output_dir = os.path.join(input_dir, "summaries")
os.makedirs(output_dir, exist_ok=True)

# Regex pattern to extract ligand, target, and loop from filename
filename_pattern = re.compile(r"ProteinMPNN_(.+?)_versus_(.+?)_(.+?)_Sequences\.fa")

# Function to read fasta files and return headers + sequences
def parse_fasta(filepath):
    entries = []
    with open(filepath, "r") as f:
        header = None
        seq = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header and seq:
                    entries.append((header, "".join(seq)))
                header = line[1:]  # remove ">"
                seq = []
            else:
                seq.append(line)
        if header and seq:
            entries.append((header, "".join(seq)))
    return entries

# Function to summarize sequences
def summarize_fasta(filepath):
    filename = os.path.basename(filepath)
    match = filename_pattern.match(filename)
    if not match:
        print(f"Skipping {filename}: filename pattern not recognized.")
        return

    ligand, target, loop = match.groups()
    entries = parse_fasta(filepath)

    # Count unique sequences
    sequences = [seq for _, seq in entries]
    seq_counts = Counter(sequences)

    # Prepare DataFrame
    rows = []
    for seq, count in seq_counts.items():
        loop_seq = seq[:10]
        if loop.upper() == "AB":
            loop_seq = seq[30:36]
        if loop.upper() == "C":
            loop_seq = seq[62:68]
        if loop.upper() == "EF":
            loop_seq = seq[92:96]
        if loop.upper() == "GH":
            loop_seq = seq[127:137]
        if loop == "Multi":
            loop_seq = seq[143:153]
        rows.append({
            "Ligand": ligand,
            "Target": target,
            "Loop": loop,
            "Sequence": seq,
            "Occurrences": count,
            "Loop_Sequence": loop_seq
        })

    df = pd.DataFrame(rows)
    df = df.sort_values(by="Occurrences", ascending=False)

    # Save CSV
    csv_path = os.path.join(output_dir, f"{filename.replace('.fa', '_summary.csv')}")
    df.to_csv(csv_path, index=False)

    # Save LaTeX table
    latex_path = os.path.join(output_dir, f"{filename.replace('.fa', '_summary.tex')}")
    with open(latex_path, "w") as f:
        f.write(df.to_latex(index=False, escape=False))

    print(f"Processed {filename}: {len(df)} unique sequences found.")
    print(f"Saved -> {csv_path}")
    print(f"Saved -> {latex_path}")

# === MAIN ===
for file in os.listdir(input_dir):
    if file.endswith(".fa"):
        summarize_fasta(os.path.join(input_dir, file))
