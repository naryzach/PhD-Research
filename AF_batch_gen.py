import os
import json
from pathlib import Path

# --- CONFIGURATION ---
INPUT_DIR = "output_fastas"           # Folder containing your FASTA files
OUTPUT_JSON = "alphafold_jobs.json"  # Output file for batch submission

# --- FUNCTION DEFINITIONS ---

def parse_fasta(fasta_path):
    """Parse a FASTA file and return a list of (header, sequence) tuples."""
    sequences = []
    header = None
    seq_lines = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    sequences.append((header, ''.join(seq_lines)))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            sequences.append((header, ''.join(seq_lines)))

    return sequences


def fasta_to_job(fasta_path):
    """Convert one FASTA file into a JSON AlphaFold job entry."""
    sequences = parse_fasta(fasta_path)
    chains = []
    for header, seq in sequences:
        chains.append({
            "proteinChain": {
                "sequence": seq,
                "count": 1
            }
        })
    
    job = {
        "name": Path(fasta_path).stem,
        "sequences": chains,
        "modelSeeds": []
    }
    return job


def main():
    fasta_files = sorted(Path(INPUT_DIR).glob("*.fasta"))
    jobs = []

    if not fasta_files:
        print(f"No FASTA files found in {INPUT_DIR}")
        return

    for fasta in fasta_files:
        job = fasta_to_job(fasta)
        jobs.append(job)
        print(f"Added job for {fasta.name} ({len(job['sequences'])} sequences)")

    with open(OUTPUT_JSON, "w") as out:
        json.dump(jobs, out, indent=2)

    print(f"\nDone. Wrote {len(jobs)} jobs to {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
