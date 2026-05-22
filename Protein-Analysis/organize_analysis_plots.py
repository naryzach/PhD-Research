import os
import re
import shutil
import glob
from collections import defaultdict, OrderedDict

# --- CONFIGURATION ---
SOURCE_DIR = "../Local"
DEST_DIR = os.path.join(SOURCE_DIR, "Notebook_Analysis_Plots")

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

# Aggregate files to include alongside the _analysis.png files.
# Searched first in Aggregate_Plots/ subdir, then in the root of the analysis dir.
AGGREGATE_FILES = [
    "Aggregate_Norm_Bind_Med_(Expr+).png",
    "Aggregate_Norm_Pos_Med_Ratio.png",
]

# --- LOGIC ---

def get_formatted_date(date_str):
    """Converts 20260403 to Apr_3."""
    try:
        month = date_str[4:6]
        day = date_str[6:8]
        month_name = MONTH_MAP.get(month, month)
        if day.startswith('0'):
            day = day[1:]
        return f"{month_name}_{day}"
    except Exception:
        return date_str


def find_source_dirs():
    """Return all TARGET_DATE[_Extra]_Renamed_Analysis dirs that have Individual_Plots/_analysis.png files."""
    # Same pattern as organize_plots.py: TARGET_DATE[_Extra]_Renamed_Analysis
    pattern = re.compile(r'^([^_]+)_(\d{8})(?:_(.*))?_Renamed_Analysis$')
    results = []
    for entry in sorted(os.scandir(SOURCE_DIR), key=lambda e: e.name):
        if not entry.is_dir():
            continue
        m = pattern.match(entry.name)
        if not m:
            continue
        individual_plots_dir = os.path.join(entry.path, "Individual_Plots")
        has_analysis = bool(glob.glob(os.path.join(individual_plots_dir, "*_analysis.png")))
        if has_analysis:
            results.append({
                'path': entry.path,
                'target': m.group(1).upper(),
                'date_raw': m.group(2),
                'date_fmt': get_formatted_date(m.group(2)),
                'extra': m.group(3) or "",
            })
    return results


def find_aggregate_files(dir_path):
    """Find the two aggregate PNGs from Positive_Ratios/ (where analyze_fcs.py saves them)."""
    pos_ratio_dir = os.path.join(dir_path, "Positive_Ratios")
    found = []
    for fname in AGGREGATE_FILES:
        candidate = os.path.join(pos_ratio_dir, fname)
        if os.path.exists(candidate):
            found.append(candidate)
        else:
            print(f"  Warning: {fname} not found in {pos_ratio_dir}")
    return found


def organize_analysis_plots(dry_run=True):
    if not os.path.exists(SOURCE_DIR):
        print(f"Error: Source directory {SOURCE_DIR} does not exist.")
        return

    data = find_source_dirs()
    if not data:
        print("No *_Renamed_Analysis directories with Individual_Plots/*_analysis.png files found.")
        return

    # Group entries by date to decide whether to add a target suffix
    by_date = defaultdict(list)
    for entry in data:
        by_date[entry['date_fmt']].append(entry)

    for entry in data:
        date_fmt = entry['date_fmt']
        target = entry['target']
        extra = entry['extra']
        peers = by_date[date_fmt]

        suffix_parts = []
        if len(peers) > 1:
            suffix_parts.append(TARGET_ABBR.get(target, target))
        if extra:
            suffix_parts.append(extra)
        suffix = ("_" + "_".join(suffix_parts)) if suffix_parts else ""

        dest_folder_name = f"Analysis_Plots_{date_fmt}{suffix}"
        dest_path = os.path.join(DEST_DIR, dest_folder_name)

        # Collect files to copy
        individual_dir = os.path.join(entry['path'], "Individual_Plots")
        analysis_files = sorted(glob.glob(os.path.join(individual_dir, "*_analysis.png")))
        aggregate_files = find_aggregate_files(entry['path'])
        all_files = analysis_files + aggregate_files

        if dry_run:
            print(f"\n[DRY RUN] {dest_folder_name}/")
            print(f"  Source: {os.path.basename(entry['path'])}")
            for f in all_files:
                print(f"    {os.path.basename(f)}")
        else:
            os.makedirs(dest_path, exist_ok=True)
            for f in all_files:
                shutil.copy2(f, os.path.join(dest_path, os.path.basename(f)))
            print(f"Copied {len(all_files)} files -> {dest_folder_name}/")
            generate_tex_template(dest_path, entry['date_raw'], target, extra, all_files)


def generate_tex_template(dest_path, date_raw, target, extra, files):
    """Generates a .tex file with figure blocks for notebook inclusion."""
    tex_path = os.path.join(dest_path, "figures.tex")

    try:
        caption_date = f"{date_raw[0:4]}-{date_raw[4:6]}-{date_raw[6:8]}"
    except Exception:
        caption_date = date_raw

    target_full = f"{target}({extra})" if extra else target

    categories = OrderedDict([
        ("Controls", []),
        ("TIMP (Positive Control)", []),
        ("ABs", []),
        ("Cs", []),
        ("ABCs", []),
        ("Aggregate Metrics", []),
    ])

    for f in sorted(files, key=os.path.basename):
        fname = os.path.basename(f)

        if fname in AGGREGATE_FILES:
            categories["Aggregate Metrics"].append(fname)
        elif fname.lower().endswith("_analysis.png"):
            stem = re.sub(r'_analysis\.png$', '', fname, flags=re.IGNORECASE)
            if re.match(r'^Negative Control', stem, re.IGNORECASE):
                categories["Controls"].append(fname)
            elif re.match(r'^TIMP', stem, re.IGNORECASE):
                categories["TIMP (Positive Control)"].append(fname)
            elif re.match(r'^AB ', stem, re.IGNORECASE):
                categories["ABs"].append(fname)
            elif re.match(r'^C ', stem, re.IGNORECASE):
                categories["Cs"].append(fname)
            elif re.match(r'^ABC ', stem, re.IGNORECASE):
                categories["ABCs"].append(fname)
            else:
                categories["Controls"].append(fname)

    rel_folder = os.path.join("Resources", os.path.basename(dest_path))

    with open(tex_path, 'w') as t:
        t.write(f"\\subsection*{{Analysis Plots ({target_full}, {caption_date})}}\n\n")

        for cat, cat_files in categories.items():
            if not cat_files:
                continue

            width = "0.45" if cat != "Aggregate Metrics" else "0.6"

            t.write(f"\\begin{{figure}}[H]\n")
            t.write(f"\t\\centering\n")
            for fname in cat_files:
                t.write(f"\t\\includegraphics[width={width}\\textwidth]{{{rel_folder}/{fname}}}\n")
            t.write(f"\t\\caption{{{target_full} analysis {caption_date} -- {cat}}}\n")
            t.write(f"\\end{{figure}}\n\n")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Organize Individual_Plots/*_analysis.png files into dated Notebook_Analysis_Plots folders."
    )
    parser.add_argument("--run", action="store_true", help="Actually copy files (omit for dry run)")
    args = parser.parse_args()

    if args.run:
        organize_analysis_plots(dry_run=False)
    else:
        print("Running in DRY RUN mode. Use --run to actually copy files.\n")
        organize_analysis_plots(dry_run=True)
