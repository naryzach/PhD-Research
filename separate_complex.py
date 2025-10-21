import os
import argparse

def process_fasta_files(input_folder, output_folder):
    """
    Reads FASTA files from an input folder, splits sequences separated by a colon,
    and writes them to new FASTA files in an output folder.
    """
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Created output directory: {output_folder}")

    # Loop through all files in the input folder
    for filename in os.listdir(input_folder):
        # Construct the full file path
        input_file_path = os.path.join(input_folder, filename)

        # Process only files, not subdirectories
        if os.path.isfile(input_file_path):
            try:
                with open(input_file_path, 'r') as f_in:
                    lines = f_in.readlines()

                # Basic validation for FASTA format
                if not lines or not lines[0].startswith('>'):
                    print(f"Warning: Skipping '{filename}' as it doesn't appear to be a FASTA file.")
                    continue
                
                # Get the header and the combined sequence line
                header = lines[0].strip().lstrip('>')
                combined_sequence = lines[1].strip()

                # Split the sequence by the colon
                sequences = combined_sequence.split(':')

                if len(sequences) != 2:
                    print(f"Warning: Could not find two colon-separated sequences in '{filename}'. Skipping.")
                    continue

                # Prepare the new content for the output file
                seq1 = sequences[0]
                seq2 = sequences[1]
                
                # Create a new filename for the output
                output_filename = f"{os.path.splitext(filename)[0]}_separated.fasta"
                output_file_path = os.path.join(output_folder, output_filename)

                # Write the two new sequences to the output file
                with open(output_file_path, 'w') as f_out:
                    f_out.write(f">{header}_seq1\n")
                    f_out.write(f"{seq1}\n")
                    f_out.write(f">{header}_seq2\n")
                    f_out.write(f"{seq2}\n")
                
                print(f"Successfully processed '{filename}' -> '{output_filename}'")

            except Exception as e:
                print(f"Error processing file {filename}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split colon-separated sequences in a folder of FASTA files."
    )
    parser.add_argument("input_folder", help="The path to the folder containing your FASTA files.")
    parser.add_argument("output_folder", help="The path to the folder where separated files will be saved.")
    
    args = parser.parse_args()
    
    process_fasta_files(args.input_folder, args.output_folder)