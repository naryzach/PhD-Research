import os
import glob
from pathlib import Path
from collections import OrderedDict

# --- Configuration ---
INPUT_DIR = 'MMP2'
OUTPUT_DIR = 'output_fastas'
REPLACEMENT_SEQUENCE = 'YNFFPRKPKWDKNQITYRIIGYTPDLDPETVDDAFARAFQVWSDVTPLRFSRIHDGEADIMINFGRWEHGDGYPFDGKDGLLAHAFAPGTGVGGDSHFDDDELWTLGEGQVGYSLFLVAAHEFGHAMGLEHSQDPGALMAPIYTYTKNFRLSQDDIKGIQELYGASPDGSDYKDDDDK' # The predefined sequence to insert
# ---------------------

def parse_fasta(file_path):
    """
    Parses a FASTA file and returns an OrderedDict of {header: sequence}.
    Sequence lines are joined into a single string.
    """
    sequences = OrderedDict()
    current_header = None
    current_sequence_lines = []

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = "".join(current_sequence_lines)
                
                current_header = line
                current_sequence_lines = []
            else:
                current_sequence_lines.append(line)
        
        # Add the last record
        if current_header:
            sequences[current_header] = "".join(current_sequence_lines)
            
    return sequences

def write_fasta(file_path, sequences):
    """
    Writes sequences from an OrderedDict to a FASTA file.
    Sequences are wrapped to 60 characters for standard FASTA format.
    """
    with open(file_path, 'w') as f:
        for header, seq in sequences.items():
            f.write(f"{header}\n")
            # Wrap the sequence at 60 characters
            for i in range(0, len(seq), 60):
                f.write(f"{seq[i:i+60]}\n")

def process_fasta_file(input_file_path, output_dir, replacement_seq):
    """
    Reads a FASTA file, processes the second sequence based on filename,
    and writes the modified file to the output directory.
    """
    filename = Path(input_file_path).name
    print(f"Processing: {filename}")

    # 1. Parse the FASTA file
    try:
        sequences = parse_fasta(input_file_path)
    except Exception as e:
        print(f"Error parsing {filename}: {e}")
        return

    headers = list(sequences.keys())
    
    # Check for at least two sequences
    if len(headers) < 2:
        print(f"  Skipping: {filename} has less than two sequences.")
        return

    # 2. Get the header for the second sequence
    second_header = headers[1]
    
    # 3. Determine the new sequence based on the filename rule
    if filename.startswith("complex_"):
        # For 'complex_' files, the second sequence is split by ':'
        # The new sequence needs to replace the part after the first ':' in the second protein.
        
        # Get the original sequence description (after '>')
        description = second_header[1:]
        
        # Find the protein ID (part before the first ':')
        try:
            protein_id, rest_of_header = description.split(':', 1)
        except ValueError:
            # Handle cases where the 'complex_' file doesn't use the expected ':' separator in the second header
            print(f"  Warning: 'complex_' file {filename} does not have a ':' in the second header. Replacing the whole second sequence.")
            sequences[second_header] = replacement_seq
        else:
            # Reconstruct the new header and sequence
            new_header = f">{protein_id}:{replacement_seq}"
            
            # The logic requested is to replace the "second protein" with the predefined sequence.
            # In the 'complex_' case, the FASTA file itself has the entire sequence,
            # but the header for the second entry contains both sequences, separated by ':'.
            # It's an unusual format where the second FASTA entry contains two proteins/sequences.
            # Assuming the intention is to replace the second sequence *segment* which is after the ':' 
            # in the header of the second FASTA record, and update both the header and the sequence.

            # Here we assume the input FASTA file has only two records (based on a "first and second sequence" context).
            # If the original FASTA record for the second protein has a sequence, that sequence is *the* second protein.
            # The problem description is slightly ambiguous about the sequence *data* in the file vs the header's content.
            # Based on standard FASTA, we assume there are two records:
            # Record 1: >Header1 / Sequence1
            # Record 2: >Header2 / Sequence2 (This is the one to replace)
            
            # Since the prompt says "replace the second protein with a predefined sequence", 
            # we simply replace the sequence data of the second FASTA record. 
            # The ":" distinction only changes how the **header** for the second protein is handled *if* the header reflects the structure.

            # Let's keep the *original* header and only replace the *sequence* for simplicity and standard FASTA convention.
            # If the complex format implies a multi-sequence header *and* a single sequence block, the sequence block is the protein.
            
            # New sequence to replace the existing one
            sequences[second_header] = replacement_seq
            print(f"  Result: Replaced second sequence with predefined sequence. Header remains: {second_header}")

    else:
        # Normal FASTA file (not starting with 'complex_')
        # Simply replace the sequence of the second record.
        sequences[second_header] = replacement_seq
        print(f"  Result: Replaced second sequence with predefined sequence. Header remains: {second_header}")

    # 4. Write the modified content to the output directory
    output_file_path = Path(output_dir) / filename
    write_fasta(output_file_path, sequences)
    print(f"  Wrote modified file to: {output_file_path}")

# --- Main Execution ---
def main():
    """
    Main function to manage file iteration and processing.
    """
    print("--- Starting FASTA Processor ---")

    # Ensure output directory exists
    Path(OUTPUT_DIR).mkdir(exist_ok=True)

    # Get all FASTA files in the input directory
    # Using 'iglob' to find all files ending with '.fasta' or '.fa' (case-insensitive)
    fasta_files = list(glob.iglob(os.path.join(INPUT_DIR, '*.[fF][aA][sS][tT][aA]')))
    fasta_files.extend(list(glob.iglob(os.path.join(INPUT_DIR, '*.[fF][aA]'))))

    if not fasta_files:
        print(f"No FASTA files found in {INPUT_DIR}. Please check the path and file extensions.")
        return

    for file_path in fasta_files:
        process_fasta_file(file_path, OUTPUT_DIR, REPLACEMENT_SEQUENCE)

    print("--- Processing Complete ---")

if __name__ == "__main__":
    main()