import os
import sys
import csv
import glob
import re
import argparse

# ---- INSTALL DEPENDENCIES IF MISSING ----
def install(package):
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

required_packages = ['matplotlib', 'pandas', 'numpy', 'scipy', 'statsmodels']


for pkg in required_packages:
    try:
        __import__(pkg)
    except ImportError:
        print(f"Installing {pkg}...")
        install(pkg)

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Use non-interactive backend for headless execution
plt.switch_backend('Agg')

from scipy import stats

def run_anova_p(df, group_col, value_col):
    """Returns the ANOVA p-value for a given value_col grouped by group_col."""
    if df.empty: return np.nan
    groups = [group[value_col].values for name, group in df.groupby(group_col)]
    # Filter out groups with < 1 data point
    groups = [g for g in groups if len(g) > 0]
    if len(groups) < 2:
        return np.nan
    try:
        f_stat, p_val = stats.f_oneway(*groups)
        return p_val
    except:
        return np.nan

def run_tukey_posthoc(df, group_col, value_col):
    """Runs Tukey HSD post-hoc test and returns the summary as a string."""
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    if df.empty: return "No data"
    
    # Filter groups with < 1 data point
    counts = df.groupby(group_col)[value_col].count()
    valid_groups = counts[counts > 0].index.tolist()
    if len(valid_groups) < 2:
        return "Insufficient groups for Post-hoc"
    
    df_f = df[df[group_col].isin(valid_groups)]
    
    try:
        tukey = pairwise_tukeyhsd(endog=df_f[value_col], groups=df_f[group_col], alpha=0.05)
        return tukey
    except Exception as e:
        return None

def format_tukey_sig_only(tukey, context_label):
    """Returns a string listing only significant pairs from a Tukey result."""
    if tukey is None: return ""
    import pandas as pd
    df = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
    sig = df[df['reject'] == True]
    if sig.empty: return ""
    
    out = [f"--- {context_label} ---"]
    for _, row in sig.iterrows():
        out.append(f"  {row['group1']} vs {row['group2']}: diff={row['meandiff']:.4f}, p-adj={row['p-adj']:.4f}")
    return "\n".join(out) + "\n"



def standardize_construct(name):
    """
    Standardize construct names to: AB, C, ABC, or TIMP followed by a number.
    Ignore everything after the number.
    """
    match = re.search(r'(ABC|AB|C|TIMP)\s*(\d+)', name, re.IGNORECASE)
    if match:
        prefix = match.group(1).upper()
        number = match.group(2)
        return f"{prefix} {number}"
    return None

# ---- CONTROL IDENTIFICATION PATTERNS ----
# Shared by the CLI (main) and the importable core used by the FCS viewer.
POS_CTRL_PATTERNS = ["POSITIVE CONTROL", "TIMP 3"]
NEG_CTRL_PATTERNS = ["NC", "NEGATIVE CONTROL"]


# ---- VENDOR / SOURCE MANIFEST ----
# Maps the short source tokens that show up in sample names and trial-folder names
# (e.g. 'Sino', 'AbC', 'Mas', 'F6') onto a canonical vendor + provenance. Edit
# Protein-Analysis/vendor_manifest.csv to add or correct mappings.
_VENDOR_MAP_CACHE = None


def _vendor_manifest_path():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "vendor_manifest.csv")


def _norm_token(token):
    """Normalize a raw source token to a lookup key: leading chunk, uppercased.
    'Mas_G' -> 'MAS', 'Arly-FLAG' -> 'ARLY', 'Sino' -> 'SINO'."""
    parts = re.split(r'[\s_\-()]+', str(token or '').strip())
    return (parts[0] if parts and parts[0] else '').upper()


def load_vendor_map(path=None, force=False):
    """Load the token -> {Vendor, Provenance} manifest (cached)."""
    global _VENDOR_MAP_CACHE
    if _VENDOR_MAP_CACHE is not None and path is None and not force:
        return _VENDOR_MAP_CACHE
    mapping = {}
    p = path or _vendor_manifest_path()
    try:
        with open(p, newline='') as f:
            for row in csv.DictReader(f):
                key = _norm_token(row.get('Token', ''))
                vendor = (row.get('Vendor', '') or '').strip()
                if key and vendor:
                    mapping[key] = {
                        'Vendor': vendor,
                        'Provenance': (row.get('Provenance', '') or '').strip() or 'Unknown',
                    }
    except FileNotFoundError:
        pass
    except Exception:
        pass
    if path is None:
        _VENDOR_MAP_CACHE = mapping
    return mapping


def canonical_vendor(token):
    """Canonical vendor name for a raw token, or None if it is not a usable source
    (empty or purely numeric). Unmapped tokens fall back to their leading chunk so
    new vendors still surface (add them to vendor_manifest.csv to canonicalize)."""
    key = _norm_token(token)
    if not key or key.isdigit():
        return None
    vm = load_vendor_map()
    if key in vm:
        return vm[key]['Vendor']
    parts = re.split(r'[\s_\-()]+', str(token).strip())
    return parts[0] if parts and parts[0] else None


def source_provenance(vendor):
    """Return 'Purchased' / 'In-house' / 'Unknown' for a canonical vendor or token."""
    vm = load_vendor_map()
    key = _norm_token(vendor)
    if key in vm:
        return vm[key]['Provenance']
    for entry in vm.values():
        if entry['Vendor'].upper() == str(vendor or '').upper():
            return entry['Provenance']
    return 'Unknown'


def derive_source(raw_name, folder_source=None):
    """Best-effort canonical vendor for a sample.

    Prefers the parenthetical in the sample name (e.g. 'AB 6-MMP9(Sino)' -> 'Sino'),
    falling back to the folder-derived token, canonicalized via vendor_manifest.csv.
    Purely numeric parentheticals (replicate labels like '(1)') are ignored; returns
    'Unknown' when no usable source is found.
    """
    tok = None
    m = re.search(r'\(([^)]+)\)', str(raw_name or ''))
    if m:
        s = m.group(1).strip()
        if s and not s.isdigit():
            tok = s
    if tok is None and folder_source:
        tok = str(folder_source).strip()
    if not tok:
        return "Unknown"
    return canonical_vendor(tok) or "Unknown"


# ---- FOLDER -> VENDOR MANIFEST ----
# Authoritative per-folder vendor/provenance (each trial folder is one target prep
# from one vendor), sourced from folder_vendor_manifest.csv. This wins over the
# per-sample token guess because it captures RG-confirmed overrides & notebook rules.
_FOLDER_MAP_CACHE = None


def _folder_manifest_path():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "folder_vendor_manifest.csv")


def _norm_folder_name(name):
    """Normalize a trial-folder name to match the manifest's Folder column: drop a
    trailing '_Analysis' (analysis dirs) so both raw and analysis folders map."""
    n = os.path.basename(str(name or '').rstrip("/\\"))
    if n.endswith("_Analysis"):
        n = n[:-len("_Analysis")]
    return n


def load_folder_map(path=None, force=False):
    """Load the folder -> {Vendor, Provenance} manifest (cached)."""
    global _FOLDER_MAP_CACHE
    if _FOLDER_MAP_CACHE is not None and path is None and not force:
        return _FOLDER_MAP_CACHE
    mapping = {}
    p = path or _folder_manifest_path()
    try:
        with open(p, newline='') as f:
            for row in csv.DictReader(f):
                folder = _norm_folder_name(row.get('Folder', ''))
                vendor = (row.get('Vendor', '') or '').strip()
                if folder and vendor:
                    mapping[folder] = {
                        'Vendor': vendor,
                        'Provenance': (row.get('Provenance', '') or '').strip() or 'Unknown',
                    }
    except FileNotFoundError:
        pass
    except Exception:
        pass
    if path is None:
        _FOLDER_MAP_CACHE = mapping
    return mapping


def folder_vendor(folder_name):
    """Authoritative canonical vendor for a trial folder, or None if not listed."""
    entry = load_folder_map().get(_norm_folder_name(folder_name))
    return entry['Vendor'] if entry else None


def resolve_vendor(folder_name=None, raw_name=None):
    """Resolve a sample's vendor: the folder manifest is authoritative; otherwise
    fall back to the per-sample token / folder token via vendor_manifest.csv."""
    if folder_name:
        v = folder_vendor(folder_name)
        if v:
            return v
    return derive_source(raw_name, folder_source(folder_name) if folder_name else None)


# ---- CHANNEL / TAG MANIFEST ----
# Which expression/binding channels a trial folder should be analyzed with.
# Some targets use a FLAG tag (different detection channels) instead of His.
# Defaults to His (FITC-A expression / APC-A binding) when a folder is not listed
# in channel_manifest.csv. Edit that file to add or correct per-folder channels.
DEFAULT_EXPR_CHANNEL = "FITC-A"
DEFAULT_BIND_CHANNEL = "APC-A"
DEFAULT_TAG = "His"
# Fallback channels per tag when a manifest row gives a Tag but leaves channels blank.
TAG_CHANNEL_DEFAULTS = {
    "HIS": ("FITC-A", "APC-A"),
    "FLAG": ("APC-A", "PE-A"),
}
_CHANNEL_MAP_CACHE = None


def _channel_manifest_path():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "channel_manifest.csv")


def load_channel_map(path=None, force=False):
    """Load the folder -> {Tag, Expr, Bind} channel manifest (cached)."""
    global _CHANNEL_MAP_CACHE
    if _CHANNEL_MAP_CACHE is not None and path is None and not force:
        return _CHANNEL_MAP_CACHE
    mapping = {}
    p = path or _channel_manifest_path()
    try:
        with open(p, newline='') as f:
            for row in csv.DictReader(f):
                folder = _norm_folder_name(row.get('Folder', ''))
                if not folder:
                    continue
                tag = (row.get('Tag', '') or '').strip() or DEFAULT_TAG
                expr = (row.get('Expr_Channel', '') or '').strip()
                bind = (row.get('Bind_Channel', '') or '').strip()
                # Fill missing channels from the tag default.
                d_expr, d_bind = TAG_CHANNEL_DEFAULTS.get(tag.upper(), (DEFAULT_EXPR_CHANNEL, DEFAULT_BIND_CHANNEL))
                mapping[folder] = {
                    'Tag': tag,
                    'Expr': expr or d_expr,
                    'Bind': bind or d_bind,
                }
    except FileNotFoundError:
        pass
    except Exception:
        pass
    if path is None:
        _CHANNEL_MAP_CACHE = mapping
    return mapping


def channels_for_folder(folder_name):
    """Return (expr_channel, bind_channel, tag) for a trial folder, from
    channel_manifest.csv. Falls back to His (FITC-A / APC-A) when not listed."""
    entry = load_channel_map().get(_norm_folder_name(folder_name)) if folder_name else None
    if entry:
        return entry['Expr'], entry['Bind'], entry['Tag']
    return DEFAULT_EXPR_CHANNEL, DEFAULT_BIND_CHANNEL, DEFAULT_TAG


def folder_source(dir_name):
    """Extract the manufacturer/source token from a trial directory name, e.g.
    'MMP9_20260701_Sino_T1_Renamed' -> 'Sino'. Returns None when absent."""
    parts = str(dir_name).split("_")
    if len(parts) >= 3:
        cand = parts[2]
        if cand and cand not in ("Renamed", "Analysis") and not re.match(r'^T\d+$', cand):
            return cand
    return None


def get_local_dir(local_dir=None):
    """Resolve the ../Local data directory relative to this script unless given."""
    if local_dir is not None:
        return local_dir
    base_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(base_dir, "../Local"))


def get_analysis_dirs(local_dir=None):
    """Return the sorted list of Local/*_Analysis trial directories (full paths)."""
    local_dir = get_local_dir(local_dir)
    return sorted([d for d in glob.glob(os.path.join(local_dir, "*_Analysis")) if os.path.isdir(d)])


def aggregate_records(records, target, date, pos_patterns=None, neg_patterns=None, folder_name=None):
    """Aggregate a trial's per-file summary_stats records into aggregate rows.

    `records` is an iterable of dict rows matching the summary_stats.csv schema
    (values may be strings, as from csv.DictReader, or numbers, as computed live
    by the FCS viewer). Mirrors the schema written to aggregate_summary.csv.
    Negative controls are dropped; positive controls are kept as their own rows.

    pos_patterns / neg_patterns override the control-identification keywords
    (case-insensitive); they default to the module constants so CLI behavior is
    unchanged. The viewer passes its "PC Keywords" setting here.
    """
    records = list(records)
    rows = []

    pos_patterns = [str(p).upper() for p in (pos_patterns or POS_CTRL_PATTERNS)]
    neg_patterns = [str(p).upper() for p in (neg_patterns or NEG_CTRL_PATTERNS)]

    # Trial-level QC from the positive control(s)
    trial_failed = False
    trial_failed_reason = ""
    pos_ctrl_dps = []
    for row in records:
        raw_name = row.get('Filename', '')
        if any(ctrl in raw_name.upper() for ctrl in pos_patterns):
            try:
                pos_ctrl_dps.append(float(row.get('Double+ %', 0)))
            except (ValueError, TypeError):
                continue
    if not pos_ctrl_dps:
        trial_failed = True
        trial_failed_reason = "No Positive Control found"
    elif max(pos_ctrl_dps) < 2.0:
        trial_failed = True
        trial_failed_reason = f"Pos Ctrl Double+ % too low ({max(pos_ctrl_dps):.1f}%)"

    # Per-construct metrics
    for row in records:
        raw_name = row.get('Filename', '')
        is_nc = any(ctrl in raw_name.upper() for ctrl in neg_patterns)
        is_pc = any(ctrl in raw_name.upper() for ctrl in pos_patterns)

        if is_nc:
            continue

        standard_name = standardize_construct(raw_name)
        if not standard_name and is_pc:
            standard_name = raw_name

        if not standard_name:
            continue

        try:
            num_match = re.search(r'(\d+)', standard_name)
            construct_num = num_match.group(1) if num_match else None

            row_target = str(row.get('Target', '') or '').strip().upper()
            if not row_target and target:
                row_target = target
            elif not row_target:
                row_target = "UNKNOWN"

            pos_med_ratio = float(row.get('Pos Med Ratio', 0))
            norm_med_ratio = float(row.get('Norm Pos Med Ratio', 0))
            pos_mean_ratio = float(row.get('Pos Mean Ratio', 0))
            norm_mean_ratio = float(row.get('Norm Pos Mean Ratio', 0))
            double_pos = float(row.get('Double+ %', 0))
            expr_pos = float(row.get('Expr+ %', 0))
            gated_events = int(float(row.get('Gated Events', 0)))

            low_expression = expr_pos < 1.0
            low_events = gated_events < 500

            bind_med_expr = float(row.get('Bind Med (Expr+)', 0))
            norm_bind_med_expr = float(row.get('Norm Bind Med (Expr+)', 0))
            bind_mean_expr = float(row.get('Bind Mean (Expr+)', 0))
            norm_bind_mean_expr = float(row.get('Norm Bind Mean (Expr+)', 0))
            expr_med_bind = float(row.get('Expr Med (Bind+)', 0))
            norm_expr_med_bind = float(row.get('Norm Expr Med (Bind+)', 0))

            bind_eff = 0
            for col in row:
                if "Binding Efficiency (DP/" in col:
                    bind_eff = float(row.get(col, 0))
                    break
            iwb_index = float(row.get('Intensity-Weighted Binding Index', 0))
            norm_iwb_index = float(row.get('Norm Intensity-Weighted Binding Index', 0))

            rows.append({
                'Target': row_target,
                'Date': date,
                'Source': resolve_vendor(folder_name, raw_name),
                'Construct': standard_name,
                'Construct Num': construct_num,
                'Raw Name': raw_name,
                'Pos Med Ratio': pos_med_ratio,
                'Norm Median Ratio': norm_med_ratio,
                'Pos Mean Ratio': pos_mean_ratio,
                'Norm Mean Ratio': norm_mean_ratio,
                'Double+ %': double_pos,
                'Bind Med (Expr+)': bind_med_expr,
                'Norm Bind Med (Expr+)': norm_bind_med_expr,
                'Bind Mean (Expr+)': bind_mean_expr,
                'Norm Bind Mean (Expr+)': norm_bind_mean_expr,
                'Expr Med (Bind+)': expr_med_bind,
                'Norm Expr Med (Bind+)': norm_expr_med_bind,
                'Binding Efficiency': bind_eff,
                'Intensity-Weighted Binding Index': iwb_index,
                'Norm Intensity-Weighted Binding Index': norm_iwb_index,
                'Expr+ %': expr_pos,
                'Gated Events': gated_events,
                'Low Expression': low_expression,
                'Low Events': low_events,
                'Trial Failed': trial_failed,
                'Trial Failed Reason': trial_failed_reason
            })
        except (ValueError, KeyError):
            continue

    return rows


def _parse_summary_stats(csv_path, target, date, verbose=True, folder_name=None):
    """Read one summary_stats.csv and aggregate it (CLI/Local path)."""
    try:
        with open(csv_path, 'r') as f:
            records = list(csv.DictReader(f))
    except Exception as e:
        if verbose:
            print(f"  Error reading {csv_path}: {e}", flush=True)
        return []
    return aggregate_records(records, target, date, folder_name=folder_name)

    return rows


def build_aggregate_data(exclude_dates=None, exclude_dirs=None, local_dir=None, verbose=True):
    """Parse every Local/*_Analysis/summary_stats.csv into a list of row dicts.

    The returned rows match the aggregate_summary.csv schema but nothing is
    written to disk. This is the importable core shared by the CLI and the FCS
    viewer's live-recalculation mode, so trials can be excluded without having to
    regenerate the static aggregate/selectivity sheets.

    exclude_dates: iterable of 'YYYYMMDD' strings to skip.
    exclude_dirs:  iterable of trial directory names or full paths to skip,
                   e.g. 'MMP2_20240510_Analysis'.
    """
    exclude_dates = set(exclude_dates or [])
    exclude_dirs = set(os.path.basename(str(d).rstrip("/\\")) for d in (exclude_dirs or []))

    local_dir = get_local_dir(local_dir)
    analysis_dirs = get_analysis_dirs(local_dir)
    if verbose:
        print(f"Found {len(analysis_dirs)} analysis directories.", flush=True)

    all_data = []
    for d in analysis_dirs:
        dir_name = os.path.basename(d)

        if dir_name in exclude_dirs:
            if verbose:
                print(f"  Skipping {dir_name}: excluded by user selection.", flush=True)
            continue

        if verbose:
            print(f"Processing {d}...", flush=True)

        # Parse target and date from directory name
        match = re.match(r'^([^_]+)_(\d{8})', dir_name)
        if not match:
            match_dt = re.match(r'^(\d{8})_(\d{6})_Analysis$', dir_name)
            if match_dt:
                date = match_dt.group(1)
                target = None
            else:
                if verbose:
                    print(f"  Skipping {dir_name}: Could not parse target and date.", flush=True)
                continue
        else:
            target = match.group(1).upper()
            date = match.group(2)

        if date in exclude_dates:
            if verbose:
                print(f"  Skipping {dir_name}: Date {date} is in exclude list.", flush=True)
            continue

        csv_path = os.path.join(d, "summary_stats.csv")
        if not os.path.exists(csv_path):
            if verbose:
                print(f"  Skipping {dir_name}: summary_stats.csv not found.", flush=True)
            continue

        all_data.extend(_parse_summary_stats(csv_path, target, date, verbose=verbose,
                                             folder_name=dir_name))

    return all_data


def build_selectivity_summary(df, metric):
    """Build the per-(Construct, Target) selectivity summary for a single metric.

    Limited to constructs tested against more than one target. Mirrors
    perform_selectivity_analysis so the viewer can compute selectivity live
    (with exclusions applied) instead of reading selectivity_summary.csv.
    Returns columns: Construct, Target, Mean, StdDev, Count, SEM, 95% CI, ANOVA_p.
    """
    import pandas as pd

    cols = ['Construct', 'Target', 'Mean', 'StdDev', 'Count', 'SEM', '95% CI', 'ANOVA_p']
    if df is None or df.empty or metric not in df.columns:
        return pd.DataFrame(columns=cols)

    stats_df = df.groupby(['Construct', 'Target'])[metric].agg(['mean', 'std', 'count']).reset_index()
    stats_df.columns = ['Construct', 'Target', 'Mean', 'StdDev', 'Count']
    stats_df['SEM'] = stats_df.apply(lambda r: r['StdDev'] / (r['Count'] ** 0.5) if r['Count'] > 1 else 0, axis=1)
    stats_df['95% CI'] = stats_df['SEM'] * 1.96

    anova_results = []
    for construct in stats_df['Construct'].unique():
        c_df = df[df['Construct'] == construct]
        anova_results.append({'Construct': construct, 'ANOVA_p': run_anova_p(c_df, 'Target', metric)})
    stats_df = stats_df.merge(pd.DataFrame(anova_results), on='Construct', how='left')

    target_counts = stats_df.groupby('Construct')['Target'].count()
    multi_target_constructs = target_counts[target_counts > 1].index.tolist()
    return stats_df[stats_df['Construct'].isin(multi_target_constructs)].copy()


def plot_stats(df, title, output_path, y_col='Binding Ratio', y_label='Binding Efficiency (Bind/Expr)', color_col=None, color_label=None):
    """Helper to generate a bar plot with optional color-coding."""
    if df.empty:
        return
        
    # Calculate statistics manually
    stats = {}
    color_stats = {}
    for _, row in df.iterrows():
        c = row['Construct']
        r = row[y_col]
        if c not in stats:
            stats[c] = []
            if color_col:
                color_stats[c] = []
        stats[c].append(r)
        if color_col:
            color_stats[c].append(row[color_col])

    sorted_constructs = sorted(stats.keys())
    means = []
    stds = []
    colors = []
    
    # Define colormap for color coding
    cmap = plt.get_cmap('viridis')
    
    for c in sorted_constructs:
        ratios = stats[c]
        n = len(ratios)
        mean = sum(ratios) / n
        means.append(mean)
        if n > 1:
            variance = sum((x - mean) ** 2 for x in ratios) / (n - 1)
            sd = variance ** 0.5
            sem = sd / (n ** 0.5)
            # 95% Confidence Interval (approx 1.96 * SEM)
            # Using 1.96 for large N, but could use t-distribution if preferred.
            # Here we use 1.96 as a standard approximation.
            stds.append(1.96 * sem)
        else:
            stds.append(0)

            
        if color_col:
            c_vals = color_stats[c]
            c_mean = sum(c_vals) / len(c_vals)
            colors.append(c_mean)

    plt.figure(figsize=(14, 8))
    
    # Global Font Sizes for Presentations/Journals
    plt.rcParams.update({
        'font.size': 12,
        'axes.titlesize': 20,
        'axes.labelsize': 18,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 14,
        'figure.titlesize': 24
    })
    
    if color_col:
        # Normalize colors for the colormap
        norm = plt.Normalize(min(colors) if colors else 0, max(colors) if colors else 100)
        bar_colors = [cmap(norm(v)) for v in colors]
        bars = plt.bar(sorted_constructs, means, yerr=stds, 
                capsize=5, color=bar_colors, edgecolor='black', alpha=0.9)
        
        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=plt.gca())
        cbar.set_label(color_label if color_label else color_col, fontsize=16)
    else:
        plt.bar(sorted_constructs, means, yerr=stds, 
                capsize=5, color='skyblue', edgecolor='black', alpha=0.8)
    
    plt.title(title, fontsize=22, pad=20)
    plt.xlabel('Sample', fontsize=18)
    plt.ylabel(y_label, fontsize=18)
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved plot to {output_path}", flush=True)

def main():
    # Paths
    parser = argparse.ArgumentParser(description="Aggregate FCS Analysis Results")
    parser.add_argument("--exclude_dates", help="Comma-separated list of dates to exclude (e.g., 20240510,20240511)")
    parser.add_argument("--data_dir", help="Directory containing the *_Analysis trial folders "
                                           "(default: ../Local). Point this at the raw flow-data root.")
    parser.add_argument("--output_dir", help="Where to write aggregate outputs "
                                            "(default: ../Local/Aggregate_FCS_Analysis).")
    args = parser.parse_args()

    exclude_dates = []
    if args.exclude_dates:
        exclude_dates = [d.strip() for d in args.exclude_dates.split(",")]

    base_dir = os.path.dirname(os.path.abspath(__file__))
    default_local = os.path.abspath(os.path.join(base_dir, "../Local"))
    local_dir = os.path.abspath(args.data_dir) if args.data_dir else default_local
    output_dir = os.path.abspath(args.output_dir) if args.output_dir else os.path.join(default_local, "Aggregate_FCS_Analysis")
    reference_csv = os.path.join(base_dir, "TIMP3_variants_names.csv")
    
    # 0. Define metrics and control patterns
    # Metrics list controls both aggregate plots and selectivity analysis
    # Aggregate metrics (Normalized) - used for within-target/trial comparisons
    aggregate_metrics = [
        {
            "y_col": "Norm Median Ratio",
            "y_label": "Normalized Median Binding Ratio (vs Pos Ctrl)",
            "title_prefix": "Binding to Expression Ratio",
            "folder": "Norm_Median_Ratio"
        },
        {
            "y_col": "Norm Bind Med (Expr+)",
            "y_label": "Normalized Median APC for FITC+ Cells",
            "title_prefix": "Binding of Expressed Cells",
            "folder": "Norm_Bind_Med_Expr_Positive"
        },
        {
            "y_col": "Binding Efficiency",
            "y_label": "Binding Efficiency (Double Positive / Expr+)",
            "title_prefix": "Binding Efficiency",
            "folder": "Binding_Efficiency"
        },
        {
            "y_col": "Norm Intensity-Weighted Binding Index",
            "y_label": "Normalized IWB Index (vs Pos Ctrl)",
            "title_prefix": "IWB Index",
            "folder": "Norm_IWB_Index"
        }
    ]

    # Selectivity metrics (Raw/Non-normalized) - used for cross-target comparisons
    selectivity_metrics = [
        {
            "y_col": "Pos Med Ratio",
            "y_label": "Median Binding Ratio (Bind/Expr)",
            "title_prefix": "Binding to Expression Ratio",
            "folder": "Median_Ratio"
        },
        {
            "y_col": "Bind Med (Expr+)",
            "y_label": "Median APC for FITC+ Cells",
            "title_prefix": "Binding of Expressed Cells",
            "folder": "Bind_Med_Expr_Positive"
        },
        {
            "y_col": "Binding Efficiency",
            "y_label": "Binding Efficiency (Double Positive / Expr+)",
            "title_prefix": "Binding Efficiency",
            "folder": "Binding_Efficiency"
        },
        {
            "y_col": "Intensity-Weighted Binding Index",
            "y_label": "IWB Index",
            "title_prefix": "IWB Index",
            "folder": "IWB_Index"
        }
    ]


    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}", flush=True)
    
    # Load reference expectations and sequences
    ref_dict = {}
    if os.path.exists(reference_csv):
        try:
            with open(reference_csv, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    construct = row.get('Construct', '').strip()
                    if construct:
                        ref_dict[construct] = {
                            'AB Sequence': row.get('AB Sequence', ''),
                            'C Sequence': row.get('C Sequence', ''),
                            'Expected': row.get('Expected', '')
                        }
        except Exception as e:
            print(f"Warning: Could not load reference CSV: {e}", flush=True)

    # Build aggregate row data via the importable core (shared with the FCS viewer)
    all_data = build_aggregate_data(exclude_dates=exclude_dates, local_dir=local_dir)

    if not all_data:
        print("No valid data points found.", flush=True)
        return

    # Create aggregate list sorted by Target and Construct
    all_data.sort(key=lambda x: (x['Target'], x['Construct']))

    # Save aggregate CSV
    output_csv = os.path.join(output_dir, "aggregate_summary.csv")
    fieldnames = ['Target', 'Date', 'Source', 'Construct', 'Raw Name', 'Pos Med Ratio', 'Norm Median Ratio', 'Pos Mean Ratio', 'Norm Mean Ratio', 'Double+ %', 'Expr+ %', 'Gated Events', 'Bind Med (Expr+)', 'Norm Bind Med (Expr+)', 'Bind Mean (Expr+)', 'Norm Bind Mean (Expr+)', 'Expr Med (Bind+)', 'Norm Expr Med (Bind+)', 'Binding Efficiency', 'Intensity-Weighted Binding Index', 'Norm Intensity-Weighted Binding Index', 'Low Expression', 'Low Events', 'Trial Failed', 'Trial Failed Reason']
    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for d in all_data:
            writer.writerow({k: d[k] for k in fieldnames})
    print(f"Saved aggregate summary to {output_csv}", flush=True)

    print("Generating plots...", flush=True)
    global_df = pd.DataFrame(all_data)
    
    # Use only valid trials and samples with sufficient expression/events for the generated bar graphs
    global_df_filtered = global_df[(global_df['Trial Failed'] == False) & (global_df['Low Expression'] == False) & (global_df['Low Events'] == False)]
    
    # Initialize master log
    master_sig_log = []

    for metric in aggregate_metrics:

        metric_dir = os.path.join(output_dir, metric["folder"])
        if not os.path.exists(metric_dir):
            os.makedirs(metric_dir)
            
        # Global Plot (using the given metric, color-coded by Double+ %)
        plot_stats(
            global_df_filtered, 
            f"{metric['title_prefix']} - All Targets", 
            os.path.join(metric_dir, f"aggregate_colorcoded.png"),
            y_col=metric['y_col'],
            y_label=metric['y_label'],
            color_col='Double+ %',
            color_label='Double Positive % (QC)'
        )
        
        # Individual Plots per Target
        targets = sorted(global_df_filtered['Target'].unique())
        for t in targets:
            target_df = global_df_filtered[global_df_filtered['Target'] == t]
            plot_stats(
                target_df,
                f"{metric['title_prefix']}: {t}",
                os.path.join(metric_dir, f"plot_{t}.png"),
                y_col=metric['y_col'],
                y_label=metric['y_label'],
                color_col='Double+ %',
                color_label='Double Positive % (QC)'
            )
            
            # RUN ANOVA across constructs for this target
            p_val = run_anova_p(target_df, 'Construct', metric['y_col'])
            if not np.isnan(p_val):
                print(f"    ANOVA (All Constructs for {t}): p = {p_val:.4e}", flush=True)
                # We can save this to a log or a specific CSV
                anova_log_path = os.path.join(metric_dir, f"anova_results_{t}.txt")
                with open(anova_log_path, 'w') as f:
                    f.write(f"Target: {t}\nMetric: {metric['y_col']}\nANOVA p-value: {p_val:.4e}\n")
                    
                    if p_val < 0.05:
                        f.write("\nPost-hoc (Tukey HSD) Analysis:\n")
                        tukey_res = run_tukey_posthoc(target_df, 'Construct', metric['y_col'])
                        if tukey_res:
                            f.write(str(tukey_res))
                            sig_str = format_tukey_sig_only(tukey_res, f"Target: {t} | Metric: {metric['y_col']}")
                            if sig_str:
                                master_sig_log.append(sig_str)




    # 2. Cross-Target Comparison Table
    print("Generating cross-target summary...", flush=True)
    cross_data = []
    constructs = global_df['Construct'].unique()
    
    for c in constructs:
        c_subset = global_df[global_df['Construct'] == c]
        
        valid_subset = c_subset[(c_subset['Trial Failed'] == False) & (c_subset['Low Expression'] == False) & (c_subset['Low Events'] == False)]
        failed_subset = c_subset[(c_subset['Trial Failed'] == True) | (c_subset['Low Expression'] == True) | (c_subset['Low Events'] == True)]
        
        # Get overall expression dates for this construct (only from valid trials)
        expression_dates_list = sorted(valid_subset['Date'].unique())
        expression_dates_str = f"{len(expression_dates_list)} dates: {', '.join(expression_dates_list)}" if expression_dates_list else ""
        
        # Determine trial dates per target (only from valid trials)
        target_counts = {}
        for tgt in ['MMP2', 'MMP3', 'MMP9', 'ADAM10', 'ADAM17']:
            tgt_dates_list = sorted(valid_subset[valid_subset['Target'] == tgt]['Date'].unique())
            tgt_dates_count = len(tgt_dates_list)
            target_counts[tgt] = f"{tgt_dates_count} dates: {', '.join(tgt_dates_list)}" if tgt_dates_count > 0 else ""
            
        # Collect comments
        comments = []
        
        if not failed_subset.empty:
            for _, row in failed_subset.drop_duplicates(['Date', 'Target']).iterrows():
                if row['Trial Failed']:
                    reason = row.get('Trial Failed Reason', 'invalid Pos Ctrl')
                    comments.append(f"Excluded {row['Date']} ({row['Target']}): {reason}")
                if row['Low Expression']:
                    comments.append(f"Excluded {row['Date']} ({row['Target']}) due to low expression ({row['Expr+ %']:.1f}%)")
                if row['Low Events']:
                    comments.append(f"Excluded {row['Date']} ({row['Target']}) due to few events ({row['Gated Events']})")
                
        # Check for unreliable double+ % in own sample (only in valid trials to avoid double reporting)
        unreliable_trials = valid_subset[valid_subset['Double+ %'] < 0.1]
        if not unreliable_trials.empty:
            for _, row in unreliable_trials.iterrows():
                comments.append(f"Low Double+ (<0.1%) on {row['Date']} ({row['Target']})")
                
        comment_str = "; ".join(comments)
        
        ref_info = ref_dict.get(c, {})
        
        cross_data.append({
            'Construct': c,
            'AB Sequence': ref_info.get('AB Sequence', ''),
            'C Sequence': ref_info.get('C Sequence', ''),
            'Expected': ref_info.get('Expected', ''),
            'Expression': expression_dates_str,
            'MMP2 Trials': target_counts['MMP2'],
            'MMP3 Trials': target_counts['MMP3'],
            'MMP9 Trials': target_counts['MMP9'],
            'ADAM10 Trials': target_counts['ADAM10'],
            'ADAM17 Trials': target_counts['ADAM17'],
            'Comments': comment_str
        })
        
    def natural_sort_key(d):
        import re
        return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', d['Construct'])]
    
    cross_data.sort(key=natural_sort_key)
    
    output_cross_csv = os.path.join(output_dir, "cross_target_summary.csv")
    cross_fieldnames = ['Construct', 'AB Sequence', 'C Sequence', 'Expected', 'Expression', 'MMP2 Trials', 'MMP3 Trials', 'MMP9 Trials', 'ADAM10 Trials', 'ADAM17 Trials', 'Comments']
    with open(output_cross_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=cross_fieldnames)
        writer.writeheader()
        writer.writerows(cross_data)
    print(f"Saved cross-target summary to {output_cross_csv}", flush=True)

    # 3. Selectivity Analysis
    # Pass the filtered global dataframe (valid trials only) and the selectivity metrics list (Raw)
    perform_selectivity_analysis(global_df_filtered, output_dir, selectivity_metrics, master_sig_log)


def perform_selectivity_analysis(df, output_dir, metrics, master_sig_log):
    """
    Identifies constructs tested against multiple targets and generates

    comparison plots and summaries for each specified metric.
    """
    print("Performing selectivity analysis...", flush=True)
    
    base_selectivity_dir = os.path.join(output_dir, "Selectivity_Analysis")
    if not os.path.exists(base_selectivity_dir):
        os.makedirs(base_selectivity_dir)
        print(f"  Created selectivity directory: {base_selectivity_dir}", flush=True)
        
    for metric_info in metrics:
        metric = metric_info["y_col"]
        folder = metric_info["folder"]
        label = metric_info.get("y_label", metric)
        title_prefix = metric_info.get("title_prefix", "Selectivity Analysis")
        
        print(f"  Analyzing selectivity for: {metric}...", flush=True)
        
        metric_dir = os.path.join(base_selectivity_dir, folder)
        if not os.path.exists(metric_dir):
            os.makedirs(metric_dir)

        # 1. Aggregate data by Construct and Target
        # Calculate Mean, StdDev, and Count for each (Construct, Target) pair
        stats_df = df.groupby(['Construct', 'Target'])[metric].agg(['mean', 'std', 'count']).reset_index()
        stats_df.columns = ['Construct', 'Target', 'Mean', 'StdDev', 'Count']
        
        # Calculate SEM and 95% CI
        stats_df['SEM'] = stats_df.apply(lambda row: row['StdDev'] / (row['Count']**0.5) if row['Count'] > 1 else 0, axis=1)
        stats_df['95% CI'] = stats_df['SEM'] * 1.96

        
        # Calculate ANOVA p-value for each construct across targets
        anova_results = []
        for construct in stats_df['Construct'].unique():
            c_df = df[df['Construct'] == construct]
            p = run_anova_p(c_df, 'Target', metric)
            anova_results.append({'Construct': construct, 'ANOVA_p': p})
        
        anova_df = pd.DataFrame(anova_results)
        stats_df = stats_df.merge(anova_df, on='Construct', how='left')

        
        # 2. Identify multi-target constructs
        target_counts = stats_df.groupby('Construct')['Target'].count()
        multi_target_constructs = target_counts[target_counts > 1].index.tolist()
        
        if not multi_target_constructs:
            print(f"    No constructs found with multiple targets for {metric} comparison.", flush=True)
            continue

        # Filter stats for these constructs
        selectivity_df = stats_df[stats_df['Construct'].isin(multi_target_constructs)].copy()
        
        # Save CSV
        csv_path = os.path.join(metric_dir, "selectivity_summary.csv")
        selectivity_df.to_csv(csv_path, index=False)
        print(f"    Saved selectivity summary to {csv_path}", flush=True)
        
        # 3. Plotting - Grouped Bar Chart
        # Use pandas plotting for grouped bars
        pivot_df = selectivity_df.pivot(index='Construct', columns='Target', values='Mean')
        pivot_std = selectivity_df.pivot(index='Construct', columns='Target', values='StdDev').fillna(0)
        
        if not pivot_df.empty:
            pivot_ci = selectivity_df.pivot(index='Construct', columns='Target', values='95% CI').fillna(0)
            plt.figure(figsize=(16, 10))
            ax = pivot_df.plot(kind='bar', yerr=pivot_ci, figsize=(16, 10), capsize=4, edgecolor='black', alpha=0.8)
            
            plt.title(f"{title_prefix}: Across Targets (Mean ± 95% CI)", fontsize=22, pad=20)

            plt.xlabel("Construct", fontsize=18)
            plt.ylabel(f"Mean {label}", fontsize=18)
            plt.xticks(rotation=45, ha='right')
            plt.legend(title="Target", fontsize=12, title_fontsize=14, loc='upper left', bbox_to_anchor=(1, 1))
            plt.grid(axis='y', linestyle='--', alpha=0.7)
            plt.tight_layout()
            
            plot_path = os.path.join(metric_dir, "selectivity_comparison_all.png")
            plt.savefig(plot_path, dpi=300)
            plt.close()
            print(f"    Saved grouped comparison plot to {plot_path}", flush=True)
        
        # 4. Individual plots for each multi-target construct
        indiv_dir = os.path.join(metric_dir, "Individual_Comparisons")
        if not os.path.exists(indiv_dir):
            os.makedirs(indiv_dir)
            
        for construct in multi_target_constructs:
            c_df = selectivity_df[selectivity_df['Construct'] == construct]
            
            plt.figure(figsize=(10, 6))
            # Use a consistent color map for targets
            colors = plt.cm.viridis(np.linspace(0, 0.8, len(c_df)))
            bars = plt.bar(c_df['Target'], c_df['Mean'], yerr=c_df['95% CI'], capsize=5, 
                           color=colors, edgecolor='black', alpha=0.9)
            
            plt.title(f"{title_prefix}: {construct}\n(Error Bars: 95% CI)", fontsize=20, pad=15)

            plt.xlabel("Target", fontsize=16)
            plt.ylabel(f"Mean {label}", fontsize=16)
            plt.grid(axis='y', linestyle='--', alpha=0.5)
            
            # Add value labels on top of bars
            for bar in bars:
                height = bar.get_height()
                plt.text(bar.get_x() + bar.get_width()/2., height + (max(c_df['Mean']) * 0.02),
                        f'{height:.2f}', ha='center', va='bottom', fontsize=11, fontweight='bold')
                
            plt.tight_layout()
            safe_name = construct.replace(" ", "_").replace("/", "-")
            plt.savefig(os.path.join(indiv_dir, f"selectivity_{safe_name}.png"), dpi=200)
            plt.close()

            
            # SAVE ANOVA + Tukey HSD text output for this construct
            anova_dir = os.path.join(metric_dir, "ANOVA_Results")
            if not os.path.exists(anova_dir):
                os.makedirs(anova_dir)
            
            p_val = c_df['ANOVA_p'].iloc[0]
            anova_txt_path = os.path.join(anova_dir, f"anova_{safe_name}.txt")
            with open(anova_txt_path, 'w') as f:
                f.write(f"Construct: {construct}\n")
                f.write(f"Metric: {metric}\n")
                f.write(f"ANOVA p-value (Across Targets): {p_val:.4e}\n")
                
                if not np.isnan(p_val) and p_val < 0.05:
                    f.write("\nPost-hoc (Tukey HSD) Analysis across Targets:\n")
                    # Get raw trials for this construct
                    c_trials = df[df['Construct'] == construct]
                    tukey_res = run_tukey_posthoc(c_trials, 'Target', metric)
                    if tukey_res:
                        f.write(str(tukey_res))
                        sig_str = format_tukey_sig_only(tukey_res, f"Construct: {construct} | Metric: {metric}")
                        if sig_str:
                            master_sig_log.append(sig_str)


        
        print(f"    Saved {len(multi_target_constructs)} individual comparison plots to {indiv_dir}", flush=True)

    # 4. Save Master Significant Results
    master_log_path = os.path.join(output_dir, "significant_tukey_summary.txt")
    with open(master_log_path, 'w') as f:
        f.write("==================================================\n")
        f.write("   MASTER SUMMARY OF SIGNIFICANT DIFFERENCES\n")
        f.write("==================================================\n\n")
        if master_sig_log:
            f.write("\n".join(master_sig_log))
        else:
            f.write("No statistically significant pairwise differences found across any analysis.\n")
    print(f"Saved master significant results to {master_log_path}", flush=True)

if __name__ == "__main__":

    main()
