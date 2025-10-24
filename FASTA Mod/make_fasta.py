import requests
import sys
import os

# --- Configuration ---
# UniProt IDs for the human proteins
UNIPROT_IDS = {
    "TIMP3": "P35625",
    "MMP9": "P14780",
    "MMP2": "P08253",
    "MMP10": "P09238",
    "ADAM10": "O14672",
    "ADAM17": "P78536"
}

# The target residues for replacement in TIMP3 (P35625).
# TIMP3 mature protein starts at residue 24 (cleavage of signal peptide 1-23)
# Start index in mature sequence (0-based)
REPLACEMENT_START_INDEX = 62
REPLACEMENT_LENGTH = 6

# Define variant sequences
VARIANT_SEQUENCES = [
    "AGESNA", "AGESNC", "AGESTA", "AGESTC", "ASESNA", 
    "ASESNC", "ASESTA", "ASESTC", "YSEDIC", "YSEDID", 
    "YSEDMC", "YSEDMD", "YSEDPC", "YSEDPD", "YKEDIC", 
    "YKEDID", "YKEDMC", "YKEDMD", "YKEDPC", "YKEDPD",
    "ASESLC" # Wild type
]

# --- Helper Functions ---

def fetch_sequence(uniprot_id: str) -> str:
    """Fetches the protein sequence for a given UniProt ID."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status() # Raise an exception for bad status codes (4xx or 5xx)

        # Sequences are returned in FASTA format, which includes a header line starting with '>'
        fasta_content = response.text.split('\n')
        # Filter out header and join the sequence lines
        sequence = "".join(line.strip() for line in fasta_content if not line.startswith('>'))

        if not sequence:
            raise ValueError(f"No sequence found for UniProt ID {uniprot_id}.")

        return sequence
    except requests.exceptions.RequestException as e:
        print(f"Error fetching sequence for {uniprot_id} from UniProt: {e}", file=sys.stderr)
        return ""
    except ValueError as e:
        print(f"Error processing FASTA content for {uniprot_id}: {e}", file=sys.stderr)
        return ""

def write_fasta_file(filename: str, sequences: dict):
    """Writes a dictionary of sequences to a FASTA file."""
    try:
        dir = list(sequences.keys())[1].split("_")[0]
        os.makedirs(dir, exist_ok=True)
        with open(dir + "/" + filename, 'w') as f:
            for header, seq in sequences.items():
                # Write header line
                f.write(f">{header}\n")
                # Write sequence (can be wrapped, but simple AlphaFold usage often takes unwrapped)
                f.write(f"{seq}\n")
        with open(dir + "/" + f"complex_{filename}", 'w') as f:
            f.write(f">{'+'.join(sequences.keys())}\n")
            # Write sequence (can be wrapped, but simple AlphaFold usage often takes unwrapped)
            f.write(f"{':'.join(sequences.values())}\n")
        print(f"Successfully created file: {filename}")
    except IOError as e:
        print(f"Error writing file {filename}: {e}", file=sys.stderr)

# --- Main Logic ---

def generate_fasta_files():
    """Fetches sequences, creates the variants, and generates all required FASTA files."""
    
    # Fetch all required sequences (only once)
    print("Fetching protein sequences from UniProt...")
    full_sequences = {}
    for name, uid in UNIPROT_IDS.items():
        full_sequences[name] = fetch_sequence(uid)
        if not full_sequences[name]:
            print("Aborting file generation due to sequence fetch failure.", file=sys.stderr)
            return

    # Check for TIMP3 full sequence length for validation (P35625 is 211 residues)
    timp3_full_seq = full_sequences["TIMP3"]
    if len(timp3_full_seq) != 211:
         print(f"Warning: TIMP3 sequence length is {len(timp3_full_seq)}, expected 211. Proceeding anyway.")
         
    # Extract the Mature TIMP3 Sequence (remove signal peptide 1-23)
    TIMP3_MATURE_WT = timp3_full_seq[23:] 
    
    # Validate the segment we are about to replace
    segment_to_replace = TIMP3_MATURE_WT[REPLACEMENT_START_INDEX : REPLACEMENT_START_INDEX + REPLACEMENT_LENGTH]
    print(f"Wild-type 6-residue segment being replaced: {segment_to_replace}")
    print(f"Total number of variant sequences to process: {len(VARIANT_SEQUENCES)}\n")
    
    # Iterate through all variants and generate files
    print("Generating FASTA files for AlphaFold-Multimer:")
    for i, variant_sequence in enumerate(VARIANT_SEQUENCES):
        print(f"--- Processing Variant {i+1}/{len(VARIANT_SEQUENCES)}: {variant_sequence} ---")
        
        # Construct the Variant TIMP3 Sequence
        TIMP3_VARIANT_SEQ = (
            TIMP3_MATURE_WT[:REPLACEMENT_START_INDEX] + 
            variant_sequence + 
            TIMP3_MATURE_WT[REPLACEMENT_START_INDEX + REPLACEMENT_LENGTH:]
        )
        
        if variant_sequence == "ASESLC": # WT sequence
            timp3_variant_header = f"TIMP3_VARIANT_WT_HUMAN|P35625"
        else: 
            timp3_variant_header = f"TIMP3_VARIANT_{variant_sequence}_HUMAN|P35625"

        # Define sequences for the two AlphaFold runs
        for name, uid in UNIPROT_IDS.items():
            if name == "TIMP3":
                continue
            target_header = f"{name}_HUMAN|{uid}"
            # Run 1: TIMP3 Variant and target
            fasta_target = {
                timp3_variant_header: TIMP3_VARIANT_SEQ,
                target_header: full_sequences[name]
            }

            # Write the files
            # Using the variant sequence in the filename for easy identification
            if variant_sequence == "ASESLC":
                variant_sequence = "WT"
            write_fasta_file(f"TIMP3_v_{name}_C_{variant_sequence}.fasta", fasta_target)
        
    print(f"\nBatch generation complete. {len(VARIANT_SEQUENCES) * (len(UNIPROT_IDS) - 1)} FASTA files\
           ({len(VARIANT_SEQUENCES)} variants * {len(UNIPROT_IDS) - 1} targets) are ready for AlphaFold-Multimer.")


if __name__ == "__main__":
    # Ensure requests library is available (though typically it is in these environments)
    try:
        import requests
        generate_fasta_files()
    except ImportError:
        print("Error: The 'requests' library is required. Please ensure it is installed.", file=sys.stderr)
