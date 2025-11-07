import os
import json
from pathlib import Path

def swap_sequences_in_json(file_path, make_backup=True):
    """Swaps the two protein sequences in each job entry of a JSON file."""
    with open(file_path, "r") as f:
        data = json.load(f)

    modified = False
    for job in data:
        sequences = job.get("sequences", [])
        if len(sequences) == 2:
            sequences[0], sequences[1] = sequences[1], sequences[0]
            modified = True

    if modified:
        if make_backup:
            backup_path = file_path.with_suffix(file_path.suffix + ".bak")
            os.replace(file_path, backup_path)
        with open(file_path, "w") as f:
            json.dump(data, f, indent=2)
        print(f"✅ Swapped sequences in {file_path.name}")
    else:
        print(f"⚠️ No swap performed in {file_path.name} (not exactly 2 sequences)")

def process_folder(folder_path):
    """Processes all JSON files in a folder."""
    folder = Path(folder_path)
    json_files = list(folder.glob("*.json"))

    if not json_files:
        print("No JSON files found in the folder.")
        return

    for json_file in json_files:
        swap_sequences_in_json(json_file)

if __name__ == "__main__":
    folder = input("Enter the path to your JSON folder: ").strip()
    process_folder(folder)
