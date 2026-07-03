import os
import re
import shutil
import glob
import argparse
from collections import defaultdict

"""
split_fcs.py

Splits a single raw acquisition folder (containing FCS files for MULTIPLE
targets) into one '_Renamed' folder per target, ready for analyze_fcs.py.

Raw file naming convention (per the flow cytometer export):
    "<Well> <Sample>-<Target>(<Source>).fcs"   e.g. "B01 AB 6-ADAM10(Abcam).fcs"
    "<Well> NC.fcs"                             e.g. "G12 NC.fcs"

For each file we:
  1. Strip the leading well id (everything before the first space).
  2. Group sample files by their target token (text after the last '-').
  3. Emit one folder per target named "<Base>_<DATE>_<Source>_Renamed"
     (e.g. "ADAM10_20260701_Abcam_Renamed"). Targets with no "(Source)"
     become "<Target>_<DATE>_Renamed".
  4. Copy the negative control(s) (files starting with "NC" / "Negative
     Control") into EVERY target folder, renaming "NC" -> "Negative Control".
  5. Resolve name collisions (e.g. duplicate replicate wells, or 2+ NCs) by
     appending " 1", " 2", ... exactly like rename_fcs.py.

Dash-less files that are not negative controls (e.g. "EtOH", "DECON") are
skipped, matching how prior _Renamed folders were produced.
"""

NEG_CONTROL_PATTERNS = ["NC", "Negative Control"]


def strip_well(filename):
    """Drop the leading well id. Returns the name (with extension) or None."""
    parts = filename.split(' ', 1)
    if len(parts) < 2:
        return None
    return parts[1]


def normalize_nc(base_name):
    """Map a leading 'NC' token to 'Negative Control'."""
    if base_name.startswith("Negative Control"):
        return base_name
    if base_name.startswith("NC"):
        return base_name.replace("NC", "Negative Control", 1)
    return base_name


def is_neg_control(base_name):
    return any(base_name.startswith(p) for p in NEG_CONTROL_PATTERNS)


def parse_target(base_name):
    """Return the target token (text after the last '-'), or None if no dash."""
    if '-' in base_name:
        return base_name.rsplit('-', 1)[1].strip()
    return None


def target_base_name(target, date):
    """'ADAM10(Abcam)' + '20260701' -> 'ADAM10_20260701_Abcam' (no _Renamed)."""
    m = re.match(r'^(.*?)\s*\(([^)]*)\)\s*$', target)
    if m:
        return f"{m.group(1).strip()}_{date}_{m.group(2).strip()}"
    return f"{target}_{date}"


def target_to_folder(target, date):
    """'ADAM10(Abcam)' + '20260701' -> 'ADAM10_20260701_Abcam_Renamed'."""
    return f"{target_base_name(target, date)}_Renamed"


def parse_well_column(filename):
    """'B01 AB 6-...' -> '01'. Returns the column string, or None."""
    well = filename.split(' ', 1)[0]
    m = re.match(r'^[A-Za-z]+(\d+)$', well)
    return m.group(1) if m else None


def extract_date(input_path):
    """Pull the leading 8-digit YYYYMMDD from the raw folder name."""
    name = os.path.basename(os.path.normpath(input_path))
    m = re.match(r'(\d{8})', name)
    return m.group(1) if m else name


def copy_with_collision(items, dest_dir):
    """items: list of (src_path, base_name, ext). Numbers duplicates."""
    groups = defaultdict(list)
    for it in items:
        groups[it[1].lower().strip()].append(it)

    count = 0
    for _, group in groups.items():
        collision = len(group) > 1
        for i, (src, base, ext) in enumerate(group):
            final = f"{base} {i+1}{ext}" if collision else f"{base}{ext}"
            try:
                shutil.copy2(src, os.path.join(dest_dir, final))
                count += 1
            except Exception as e:
                print(f"  Error copying {final}: {e}")
    return count


def main():
    parser = argparse.ArgumentParser(
        description="Split a raw multi-target FCS folder into per-target _Renamed folders.")
    parser.add_argument("-i", "--input", required=True,
                        help="Raw acquisition folder containing mixed-target .fcs files.")
    parser.add_argument("-o", "--output", default=None,
                        help="Parent directory for the _Renamed folders (default: input's parent).")
    parser.add_argument("--date", default=None,
                        help="Override the YYYYMMDD tag (default: parsed from the input folder name).")
    parser.add_argument("--nc", nargs='+', default=None,
                        help="Explicit negative-control .fcs file(s) to use for EVERY target folder, "
                             "overriding any NC files auto-detected in the input. Each is copied in as "
                             "'Negative Control.fcs' (numbered if more than one).")
    parser.add_argument("--by-column", action="store_true",
                        help="Split each target's plate columns into separate trial folders "
                             "(lowest column = _T1, next = _T2, ...) instead of collapsing "
                             "replicate wells into one folder.")
    args = parser.parse_args()

    src_path = os.path.abspath(args.input)
    parent = os.path.abspath(args.output) if args.output else os.path.dirname(src_path)
    date = args.date or extract_date(src_path)

    fcs_files = glob.glob(os.path.join(src_path, "*.fcs"))
    if not fcs_files:
        print(f"No FCS files found in {src_path}")
        return

    print(f"Found {len(fcs_files)} FCS files. Date tag: {date}")

    # Partition into per-target sample lists + a shared negative-control list.
    # In --by-column mode we key on (target_base, column) so we can assign
    # trial tags after seeing every column present for each target.
    target_items = defaultdict(list)   # folder_name -> [(src, base, ext), ...]
    tc_items = defaultdict(list)       # (target_base, column) -> [items]  (by-column)
    base_columns = defaultdict(set)    # target_base -> {columns}          (by-column)
    control_items = []                 # [(src, base, ext), ...]  (NC -> Negative Control)
    skipped = []

    for file_path in fcs_files:
        filename = os.path.basename(file_path)
        stripped = strip_well(filename)
        if stripped is None:
            print(f"Skipping {filename}: no well id (no space found).")
            skipped.append(filename)
            continue

        base_name, ext = os.path.splitext(stripped)

        if is_neg_control(base_name):
            control_items.append((file_path, normalize_nc(base_name), ext))
            continue

        target = parse_target(base_name)
        if target is None:
            print(f"Skipping {filename}: not a target sample or control.")
            skipped.append(filename)
            continue

        if args.by_column:
            col = parse_well_column(filename)
            if col is None:
                print(f"Skipping {filename}: could not parse plate column from well id.")
                skipped.append(filename)
                continue
            tbase = target_base_name(target, date)
            tc_items[(tbase, col)].append((file_path, base_name, ext))
            base_columns[tbase].add(col)
        else:
            folder = target_to_folder(target, date)
            target_items[folder].append((file_path, base_name, ext))

    # Resolve by-column groups into trial-tagged folders.
    if args.by_column:
        for tbase, cols in base_columns.items():
            for idx, col in enumerate(sorted(cols, key=int), start=1):
                folder = f"{tbase}_T{idx}_Renamed"
                target_items[folder] = tc_items[(tbase, col)]

    if not target_items:
        print("No target sample files identified. Nothing to do.")
        return

    # Override the negative control(s) if explicit files were supplied.
    if args.nc:
        control_items = []
        for p in args.nc:
            p = os.path.abspath(p)
            if not os.path.exists(p):
                print(f"Warning: --nc file not found, skipping: {p}")
                continue
            control_items.append((p, "Negative Control", os.path.splitext(p)[1]))
        print(f"Overriding negative control(s) with {len(control_items)} supplied file(s): "
              f"{[os.path.basename(p) for p in args.nc]}")

    print(f"\nIdentified {len(target_items)} target(s) and {len(control_items)} negative control file(s):")
    for folder in sorted(target_items):
        print(f"  - {folder}: {len(target_items[folder])} samples")

    # Create each folder, copy its samples, then distribute the controls.
    for folder, items in sorted(target_items.items()):
        dest_dir = os.path.join(parent, folder)
        if os.path.exists(dest_dir):
            print(f"\nCleaning existing directory: {dest_dir}")
            shutil.rmtree(dest_dir)
        os.makedirs(dest_dir)

        n_samples = copy_with_collision(items, dest_dir)
        n_ctrl = copy_with_collision(control_items, dest_dir)
        print(f"[{folder}] {n_samples} samples + {n_ctrl} controls copied.")

    if skipped:
        print(f"\nSkipped {len(skipped)} non-target/non-control file(s): {skipped}")

    print("\nSplit complete.")


if __name__ == "__main__":
    main()
