import os
import csv
import numpy as np
import pandas as pd
import glob

# === CODON TRANSLATION HELPERS ===
def translate(codon):
    codon_dict = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCA": "S", "TCG": "S", "TCC": "S", "TCT": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "-", "TAG": "-",
        "TGC": "C", "TGT": "C", "TGA": "-", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACA": "T", "ACT": "T", "ACC": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }
    return codon_dict.get(codon, "")

def complement(base):
    complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return complement_dict.get(base, base)

def reverse_complement(sequence):
    return ''.join(complement(base) for base in sequence[::-1])

def translate_sequence(sequence, frame=0):
    return ''.join(translate(sequence[i:i+3]) for i in range(frame, len(sequence) - 2, 3))

# === FASTQ READING & PROCESSING ===
def read_fastq(fastq_path):
    """
    Reads a FASTQ file and returns a list of (ID, Sequence, Quality).
    """
    records = []
    with open(fastq_path, 'r') as fq:
        while True:
            id_line = fq.readline().strip()
            if not id_line:
                break
            seq = fq.readline().strip()
            fq.readline()  # skip '+'
            qual = fq.readline().strip()
            seq_id = id_line[1:] if id_line.startswith('@') else id_line
            records.append((seq_id, seq, qual))
    return records

def phred_score(qual_str):
    """
    Convert ASCII-encoded Phred33 string to average quality score.
    """
    scores = [ord(c) - 33 for c in qual_str]
    return np.mean(scores) if scores else np.nan

def process_fastq(fastq_path, translate_flag=True):
    """
    Process FASTQ → DataFrame with translation & quality metrics.
    """
    records = read_fastq(fastq_path)
    if not records:
        print(f"No reads found in {fastq_path}")
        return None

    df = pd.DataFrame(records, columns=["ID", "Sequence", "Quality"])
    df["Seq_Length"] = df["Sequence"].str.len()
    df["Avg_Quality"] = df["Quality"].apply(phred_score)

    # Remove ambiguous sequences
    df = df[~df["Sequence"].str.contains('N')]

    if translate_flag:
        print("Translating sequences (6 frames)...")
        for i in range(1, 4):
            df[f'Frame_fwd_{i}'] = df["Sequence"].apply(lambda s: translate_sequence(s, frame=i-1))
            df[f'Frame_rev_{i}'] = df["Sequence"].apply(lambda s: translate_sequence(reverse_complement(s), frame=i-1))
    else:
        print("Skipping translation (--no-translate enabled)")

    return df

def save_to_csv(df, output_path):
    """
    Save DataFrame to plain CSV.
    """
    df.to_csv(output_path + ".csv", index=False)
    print(f"Saved to {output_path}.csv")

def batch_convert_fastq(input_folder, translate_flag=True):
    """
    Convert all FASTQ files in a folder.
    """
    fastq_files = glob.glob(os.path.join(input_folder, "*.fastq")) + \
                  glob.glob(os.path.join(input_folder, "*.fq"))

    if not fastq_files:
        print("No FASTQ files found in the folder.")
        return

    for fastq_file in fastq_files:
        print(f"Processing {fastq_file}...")
        df = process_fastq(fastq_file, translate_flag)
        if df is not None:
            output_path = os.path.splitext(fastq_file)[0]
            save_to_csv(df, output_path)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Convert FASTQ to CSV with optional translation & quality info.")
    parser.add_argument("path", help="Path to a FASTQ file or folder.")
    parser.add_argument("--no-translate", action="store_true", help="Skip amino acid translation for faster CSV generation.")
    args = parser.parse_args()

    translate_flag = not args.no_translate

    if os.path.isdir(args.path):
        batch_convert_fastq(args.path, translate_flag)
    else:
        df = process_fastq(args.path, translate_flag)
        if df is not None:
            output_path = os.path.splitext(args.path)[0]
            save_to_csv(df, output_path)
