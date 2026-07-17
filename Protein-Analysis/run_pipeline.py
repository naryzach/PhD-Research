import os
import sys
import argparse
import subprocess

# Determine paths relative to this script
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)
analyze_script = os.path.join(script_dir, "analyze_fcs.py")
aggregate_script = os.path.join(script_dir, "aggregate_analysis.py")

# Per-folder expression/binding channels come from channel_manifest.csv (default: His).
import aggregate_analysis

parser = argparse.ArgumentParser(description="Run per-folder FCS analysis then aggregate.")
parser.add_argument("--data_dir", help="Directory containing the *_Renamed trial folders "
                                        "(default: ../Local). Point this at the raw flow-data root, "
                                        "e.g. 'D:/Ryan Gustafson/Flow Cytometry Data'.")
parser.add_argument("--only", help="Comma-separated folder names to (re)process instead of all "
                                   "(e.g. for fixing just the folders with corrected channels).")
args = parser.parse_args()

data_dir = os.path.abspath(args.data_dir) if args.data_dir else os.path.abspath(os.path.join(script_dir, "..", "Local"))
only = set(s.strip() for s in args.only.split(",")) if args.only else None

print(f"Targeting data directory: {data_dir}")

if not os.path.exists(data_dir):
    print(f"Error: Could not find directory {data_dir}")
    sys.exit(1)

renamed_folders = [
    f.path for f in os.scandir(data_dir)
    if f.is_dir() and f.name.endswith("_Renamed")
]
if only:
    renamed_folders = [f for f in renamed_folders if os.path.basename(f) in only]

# Persistent exclusions from excluded_folders.csv (skip without deleting).
_excluded = aggregate_analysis.load_excluded_folders()
_skipped = [f for f in renamed_folders if aggregate_analysis.is_excluded(os.path.basename(f))]
renamed_folders = [f for f in renamed_folders if not aggregate_analysis.is_excluded(os.path.basename(f))]
if _skipped:
    print(f"Skipping {len(_skipped)} folder(s) listed in excluded_folders.csv:")
    for f in _skipped:
        name = aggregate_analysis._norm_folder_name(os.path.basename(f))
        print(f"  - {os.path.basename(f)} ({_excluded.get(name, '')})")

if not renamed_folders:
    print("No '_Renamed' folders found (matching filter)." if only else "No '_Renamed' folders found.")
    sys.exit(0)

print(f"Found {len(renamed_folders)} target folder(s):")
for f in renamed_folders:
    print(f"  - {os.path.basename(f)}")

# 1. Run analyze_fcs.py on each isolated folder with its manifest-defined channels
for folder in renamed_folders:
    folder_name = os.path.basename(folder)
    print(f"\n{'='*60}\n[PROCESSING] {folder_name}\n{'='*60}")

    expr_ch, bind_ch, tag = aggregate_analysis.channels_for_folder(folder_name)
    print(f"Channels (tag={tag}): Expr = {expr_ch}, Bind = {bind_ch}"
          + ("  [channel_manifest.csv]" if tag != aggregate_analysis.DEFAULT_TAG else "  [default His]"))

    cmd = ["python", analyze_script, "-i", folder,
           "--expr_channel", expr_ch, "--bind_channel", bind_ch]
    try:
        subprocess.run(cmd, check=True, cwd=script_dir)
    except subprocess.CalledProcessError:
        print(f"Warning: Script failed on {folder}. Proceeding to next...")

# 2. Run aggregate_analysis.py over the same data directory
print(f"\n{'='*60}\n[STARTING AGGREGATE ANALYSIS]\n{'='*60}")
if os.path.exists(aggregate_script):
    agg_cmd = ["python", aggregate_script, "--data_dir", data_dir]
    try:
        subprocess.run(agg_cmd, check=True, cwd=script_dir)
        print("\n=== PIPELINE COMPLETE ===")
    except subprocess.CalledProcessError:
        print("Error: Aggregate analysis failed.")
else:
    print(f"Error: {aggregate_script} not found. Skipping aggregation.")
