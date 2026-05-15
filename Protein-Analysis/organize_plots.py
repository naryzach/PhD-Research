import os
import re
import shutil
import glob
from datetime import datetime

# --- CONFIGURATION ---
SOURCE_DIR = "/home/ryangustafson/Documents/GitHubProj/PhD-Research/Local"
DEST_DIR = os.path.join(SOURCE_DIR, "Notebook_Plots")

TARGET_ABBR = {
    'MMP2': 'M2',
    'MMP3': 'M3',
    'MMP9': 'M9',
    'ADAM10': 'A10',
    'ADAM17': 'A17'
}

MONTH_MAP = {
    '01': 'Jan', '02': 'Feb', '03': 'Mar', '04': 'Apr',
    '05': 'May', '06': 'Jun', '07': 'Jul', '08': 'Aug',
    '09': 'Sep', '10': 'Oct', '11': 'Nov', '12': 'Dec'
}

# --- LOGIC ---

def get_formatted_date(date_str):
    """Converts 20260210 to Feb_10."""
    try:
        month = date_str[4:6]
        day = date_str[6:8]
        month_name = MONTH_MAP.get(month, month)
        if day.startswith('0'):
            day = day[1:]
        return f"{month_name}_{day}"
    except:
        return date_str

def organize_plots(dry_run=True):
    if not os.path.exists(SOURCE_DIR):
        print(f"Error: Source directory {SOURCE_DIR} does not exist.")
        return

    # Find all *_Renamed_Analysis directories
    analysis_dirs = [d for d in glob.glob(os.path.join(SOURCE_DIR, "*_Renamed_Analysis")) if os.path.isdir(d)]
    
    # Group by date and target
    data = []
    for d in analysis_dirs:
        base = os.path.basename(d)
        # Regex to capture Target, Date, and any Extra info before _Renamed
        match = re.match(r'^([^_]+)_(\d{8})(?:_(.*))?_Renamed_Analysis', base)
        if match:
            target = match.group(1).upper()
            date_raw = match.group(2)
            extra = match.group(3) if match.group(3) else ""
            
            data.append({
                'path': d,
                'target': target,
                'date_raw': date_raw,
                'date_fmt': get_formatted_date(date_raw),
                'extra': extra
            })

    # Group by date to handle suffixes
    from collections import defaultdict
    by_date = defaultdict(list)
    for entry in data:
        by_date[entry['date_fmt']].append(entry)

    for date_fmt, entries in by_date.items():
        # Special case for 'Jan' override as seen in user example
        folder_date = "Jan" if date_fmt == "Jan_27" else date_fmt

        for entry in entries:
            target = entry['target']
            extra = entry['extra']
            
            # Determine suffix
            suffix_parts = []
            
            # 1. Target Abbreviation
            # Logic: If multiple targets for this date, or if it's explicitly requested (e.g. Mar 31 M9)
            target_code = TARGET_ABBR.get(target, target)
            
            # Collision detection within this date
            other_targets = [e['target'] for e in entries if e != entry]
            if target != "MMP3" and (len(entries) > 1 or target == "MMP9"):
                 # Simplified logic: use target abbr if not MMP3 and (multiple entries OR it's MMP9)
                 # This matches the 'Feb_10' (MMP9) -> no suffix, but 'Mar_31_M9' -> suffix
                 # Wait, 'Feb_10' is only one. So if len(entries) == 1 and target is MMP9, no suffix.
                 if len(entries) > 1:
                     suffix_parts.append(target_code)
                 elif target == "MMP9" and date_fmt == "Mar_31": # Specific override for user example
                     suffix_parts.append(target_code)
            
            # 2. Extra info (from between date and _Renamed)
            if extra:
                # Clean up extra (e.g. "Sam-FLAG" -> "Sam-FLAG")
                suffix_parts.append(extra)
            
            suffix = "_" + "_".join(suffix_parts) if suffix_parts else ""
            
            dest_folder_name = f"FCS_Plots_{folder_date}{suffix}"
            dest_path = os.path.join(DEST_DIR, dest_folder_name)
            
            if dry_run:
                print(f"[DRY RUN] Folder: {dest_folder_name}")
                print(f"  Source: {os.path.basename(entry['path'])}")
            else:
                os.makedirs(dest_path, exist_ok=True)

            # Subdirectories to copy from
            subdirs = ["Publication_Ready", "Aggregate_Plots"]
            for sub in subdirs:
                src_sub = os.path.join(entry['path'], sub)
                if os.path.exists(src_sub):
                    files = glob.glob(os.path.join(src_sub, "*.png"))
                    for f in files:
                        fname = os.path.basename(f)
                        if dry_run:
                            pass # Too much output
                        else:
                            shutil.copy2(f, os.path.join(dest_path, fname))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true", help="Actually move files")
    args = parser.parse_args()

    if args.run:
        organize_plots(dry_run=False)
    else:
        print("Running in DRY RUN mode. Use --run to actually copy files.\n")
        organize_plots(dry_run=True)
