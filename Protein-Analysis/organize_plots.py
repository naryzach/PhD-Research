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
            subdirs = {
                "Publication_Ready": "*.png",
                "Aggregate_Plots": "*.png",
                "Positive_Distributions": "*_PosDist.png",
                "Filtered_Histograms": "*_Bind_Filtered.png",
                "Filtered_Histograms/Binding": "*_Bind_Filtered.png",
                "Positive_Ratios": "*.png"
            }
            
            copied_files = set() # Use set to avoid duplicates if found in multiple places
            
            # 1. Root files (like Gating_Strategy)
            root_patterns = ["Gating_Strategy_*.png"]
            for pat in root_patterns:
                files = glob.glob(os.path.join(entry['path'], pat))
                for f in files:
                    fname = os.path.basename(f)
                    if not dry_run:
                        shutil.copy2(f, os.path.join(dest_path, fname))
                    copied_files.add(fname)

            # 2. Subdirectories
            for sub, pattern in subdirs.items():
                src_sub = os.path.join(entry['path'], sub)
                if os.path.exists(src_sub):
                    files = glob.glob(os.path.join(src_sub, pattern))
                    for f in files:
                        fname = os.path.basename(f)
                        if not dry_run:
                            shutil.copy2(f, os.path.join(dest_path, fname))
                        copied_files.add(fname)

            if not dry_run:
                generate_tex_template(dest_path, entry['date_raw'], target, extra, list(copied_files))

def generate_tex_template(dest_path, date_raw, target, extra, files):
    """Generates a .tex file with figure blocks for easy notebook inclusion."""
    tex_path = os.path.join(dest_path, "figures.tex")
    
    # Target abbreviation for filename consistency
    target_abbr_map = {'MMP2': 'MMP2', 'MMP3': 'MMP3', 'MMP9': 'MMP9', 'ADAM10': 'ADAM10', 'ADAM17': 'ADAM17'}
    target_name = target_abbr_map.get(target, target)
    target_full = f"{target_name}({extra})" if extra else target_name
    
    # Format date as YYYY-MM-DD
    try:
        caption_date = f"{date_raw[0:4]}-{date_raw[4:6]}-{date_raw[6:8]}"
    except:
        caption_date = date_raw

    # Define categories and their specific files/rules
    # Using an OrderedDict to maintain sequence from user example
    from collections import OrderedDict
    categories = OrderedDict([
        ("Controls", []),
        ("ABs", []),
        ("Cs", []),
        ("ABCs", []),
        ("Aggregate Expression and Binding", []),
        ("Negative Control Comparison Metrics", []),
        ("Positive Control Comparison Metrics", []),
        ("Positive Distributions", []),
        ("Filtered APC Plots", []),
        ("Gating Strategy", [])
    ])
    
    # Sort files to ensure deterministic order
    for f in sorted(files):
        # 1. Gating
        if "Gating_Strategy" in f:
            categories["Gating Strategy"].append(f)
        # 2. Aggregates
        elif "Aggregate_Ridgeline" in f:
            # Only include _All ridgelines as per user request
            if "_All.png" in f:
                categories["Aggregate Expression and Binding"].append(f)
        elif "Aggregate_FoldChange" in f:
            categories["Negative Control Comparison Metrics"].append(f)
        elif "Aggregate_Norm" in f:
            categories["Positive Control Comparison Metrics"].append(f)
        # 3. Distributions and Filtered
        elif "_PosDist" in f:
            # Only include PosDist for positive controls (TIMP or Positive Control)
            if "TIMP" in f.upper() or "POSITIVE CONTROL" in f.upper():
                categories["Positive Distributions"].append(f)
        elif "_Bind_Filtered" in f:
            categories["Filtered APC Plots"].append(f)
        # 4. Individuals (Controls, ABs, Cs, ABCs)
        elif "_iso.png" in f:
            if "Negative Control" in f or "TIMP" in f:
                categories["Controls"].append(f)
            elif f.startswith("AB "):
                categories["ABs"].append(f)
            elif f.startswith("C "):
                categories["Cs"].append(f)
            elif f.startswith("ABC "):
                categories["ABCs"].append(f)

    with open(tex_path, 'w') as t:
        t.write(f"\\subsection*{{Data \\& Results ({target_full})}}\n")
        
        rel_folder = os.path.join("Resources", os.path.basename(dest_path))
        
        for cat, cat_files in categories.items():
            if not cat_files:
                continue
            
            t.write(f"\t\\begin{{figure}}[H]\n")
            t.write(f"\t\t\\centering\n")
            
            # Determine width based on user example
            width = "0.45"
            if "Aggregate Expression and Binding" in cat:
                width = "0.7"
            elif "Comparison Metrics" in cat:
                width = "0.6"
            elif "Gating Strategy" in cat:
                width = "0.6"
            
            for f in cat_files:
                t.write(f"\t\t\\includegraphics[width={width}\\textwidth]{{{rel_folder}/{f}}}\n")
            
            t.write(f"\t\t\\caption{{TwistBio expression flow {caption_date} ({target_full}) -- {cat}}}\n")
            t.write(f"\t\\end{{figure}}\n\n")

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
