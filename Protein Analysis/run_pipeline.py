import os
import subprocess
import sys

# Determine paths relative to this script
script_dir = os.path.dirname(os.path.abspath(__file__))
local_dir = os.path.abspath(os.path.join(script_dir, "..", "Local"))
analyze_script = os.path.join(script_dir, "analyze_fcs.py")
aggregate_script = os.path.join(script_dir, "aggregate_analysis.py")

print(f"Targeting data directory: {local_dir}")

# Find all _Renamed folders in the Local directory
if not os.path.exists(local_dir):
    print(f"Error: Could not find directory {local_dir}")
    sys.exit(1)

renamed_folders = [
    f.path for f in os.scandir(local_dir) 
    if f.is_dir() and f.name.endswith("_Renamed")
]

if not renamed_folders:
    print("No '_Renamed' folders found.")
else:
    print(f"Found {len(renamed_folders)} target folder(s):")
    for f in renamed_folders: 
        print(f"  - {os.path.basename(f)}")
    
    # 1. Run analyze_fcs.py on each isolated folder
    for folder in renamed_folders:
        print(f"\n{'='*60}\n[PROCESSING] {os.path.basename(folder)}\n{'='*60}")
        try:
            subprocess.run(["python", analyze_script, "-i", folder], check=True, cwd=script_dir)
        except subprocess.CalledProcessError:
            print(f"Warning: Script failed on {folder}. Proceeding to next...")
            
    # 2. Run aggregate_analysis.py
    print(f"\n{'='*60}\n[STARTING AGGREGATE ANALYSIS]\n{'='*60}")
    if os.path.exists(aggregate_script):
        try:
            subprocess.run(["python", aggregate_script], check=True, cwd=script_dir)
            print("\n✔️ PIPELINE COMPLETE ✔️")
        except subprocess.CalledProcessError:
            print("Error: Aggregate analysis failed.")
    else:
        print(f"Error: {aggregate_script} not found. Skipping aggregation.")
