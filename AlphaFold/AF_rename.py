import os
import re

# Define the directory where your zip files are located
target_directory = './AlphaFold'

# Define the regex pattern to match and capture the parts of the filename
# Group 1: ligand (e.g., timp3)
# Group 2: variant (e.g., agestc)
# Group 3: target (e.g., mmp2)
pattern = re.compile(r"^fold_([a-zA-Z0-9]+)_variant_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)\.zip$", re.IGNORECASE)

print(f"Starting file renaming in: {os.path.abspath(target_directory)}\n")

renamed_count = 0
for filename in os.listdir(target_directory):
    if filename.endswith(".zip"):
        match = pattern.match(filename)
        
        if match:
            # Extract the captured groups
            ligand = match.group(1)
            variant = match.group(2)
            target = match.group(3)
            
            # Construct the new filename: {ligand}_variant_{target}_{variant}.zip
            new_filename = f"{ligand}_variant_{target}_{variant}.zip"
            
            # Full paths for renaming
            old_path = os.path.join(target_directory, filename)
            new_path = os.path.join(target_directory, new_filename)
            
            # Check if the new filename is different from the old one
            if old_path != new_path:
                try:
                    os.rename(old_path, new_path)
                    print(f"Renamed: {filename} -> {new_filename}")
                    renamed_count += 1
                except OSError as e:
                    print(f"Error renaming {filename}: {e}")
            else:
                # This should rarely happen but handles cases where old_path == new_path
                # based on the pattern match logic.
                pass 

print(f"\nCompleted. Total files renamed: {renamed_count}")