import os
import sys
import csv
import glob
import re

# ---- INSTALL DEPENDENCIES IF MISSING ----
def install(package):
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

required_packages = ['matplotlib', 'pandas', 'numpy']
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
        mean = sum(ratios) / len(ratios)
        means.append(mean)
        if len(ratios) > 1:
            variance = sum((x - mean) ** 2 for x in ratios) / (len(ratios) - 1)
            stds.append(variance ** 0.5)
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
    base_dir = os.path.dirname(os.path.abspath(__file__))
    local_dir = os.path.abspath(os.path.join(base_dir, "../Local"))
    output_dir = os.path.join(local_dir, "Aggregate_Analysis")
    reference_csv = os.path.join(base_dir, "TIMP3_variants_names.csv")
    
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}", flush=True)
    
    # Load reference expectations
    expectations = {}
    if os.path.exists(reference_csv):
        try:
            with open(reference_csv, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    num = row.get('Number')
                    if num:
                        expectations[num] = row.get('Expectation', 'N/A')
        except Exception as e:
            print(f"Warning: Could not load reference CSV: {e}", flush=True)

    # Find all *_Analysis directories
    analysis_dirs = sorted([d for d in glob.glob(os.path.join(local_dir, "*_Analysis")) if os.path.isdir(d)])
    print(f"Found {len(analysis_dirs)} analysis directories.", flush=True)
    
    if not analysis_dirs:
        print(f"No analysis directories found in {local_dir}", flush=True)
        return

    all_data = []

    for d in analysis_dirs:
        print(f"Processing {d}...", flush=True)
        dir_name = os.path.basename(d)
        
        # Parse target and date from directory name
        match = re.match(r'^([^_]+)_(\d{8})', dir_name)
        if not match:
            print(f"  Skipping {dir_name}: Could not parse target and date.", flush=True)
            continue
            
        target = match.group(1).upper()
        date = match.group(2)
        
        csv_path = os.path.join(d, "summary_stats.csv")
        if not os.path.exists(csv_path):
            print(f"  Skipping {dir_name}: summary_stats.csv not found.", flush=True)
            continue
            
        try:
            with open(csv_path, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    raw_name = row['Filename']
                    
                    # Filter out controls
                    if any(ctrl in raw_name.upper() for ctrl in ["NC", "NEGATIVE CONTROL", "POSITIVE CONTROL"]):
                        continue
                    
                    standard_name = standardize_construct(raw_name)
                    if standard_name:
                        try:
                            # Extract number for expectation mapping
                            num_match = re.search(r'(\d+)', standard_name)
                            construct_num = num_match.group(1) if num_match else None
                            
                            # Extract required metrics
                            norm_ratio = float(row.get('Norm Pos Med Ratio', 0))
                            double_pos = float(row.get('Double+ %', 0))
                            
                            # New Filtered Metrics
                            norm_bind_med_expr = float(row.get('Norm Bind Med (Expr+)', 0))
                            norm_expr_med_bind = float(row.get('Norm Expr Med (Bind+)', 0))
                            
                            all_data.append({
                                'Target': target,
                                'Date': date,
                                'Construct': standard_name,
                                'Construct Num': construct_num,
                                'Raw Name': raw_name,
                                'Norm Median Ratio': norm_ratio,
                                'Double+ %': double_pos,
                                'Norm Bind Med (Expr+)': norm_bind_med_expr,
                                'Norm Expr Med (Bind+)': norm_expr_med_bind
                            })
                        except (ValueError, KeyError):
                            continue
        except Exception as e:
            print(f"  Error processing {csv_path}: {e}", flush=True)

    if not all_data:
        print("No valid data points found.", flush=True)
        return

    # Create aggregate list sorted by Target and Construct
    all_data.sort(key=lambda x: (x['Target'], x['Construct']))

    # Save aggregate CSV
    output_csv = os.path.join(output_dir, "aggregate_summary.csv")
    fieldnames = ['Target', 'Date', 'Construct', 'Raw Name', 'Norm Median Ratio', 'Double+ %', 'Norm Bind Med (Expr+)', 'Norm Expr Med (Bind+)']
    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for d in all_data:
            writer.writerow({k: d[k] for k in fieldnames})
    print(f"Saved aggregate summary to {output_csv}", flush=True)

    # 1. Plotting
    print("Generating plots...", flush=True)
    global_df = pd.DataFrame(all_data)
    
    # Global Plot (using Normalized Median Ratio, color-coded by Double+ %)
    plot_stats(
        global_df, 
        'Normalized Binding Ratio (Color-coded by Double+ %)', 
        os.path.join(output_dir, "aggregate_binding_colorcoded.png"),
        y_col='Norm Median Ratio',
        y_label='Normalized Median Binding Ratio (vs Pos Ctrl)',
        color_col='Double+ %',
        color_label='Double Positive % (QC)'
    )
    
    # Individual Plots per Target
    targets = sorted(global_df['Target'].unique())
    for t in targets:
        target_df = global_df[global_df['Target'] == t]
        plot_stats(
            target_df,
            f'Binding to Expression Ratio: {t}',
            os.path.join(output_dir, f"binding_plot_{t}.png"),
            y_col='Norm Median Ratio',
            y_label='Normalized Median Binding Ratio (vs Pos Ctrl)',
            color_col='Double+ %',
            color_label='Double Positive % (QC)'
        )

    # 2. Cross-Target Comparison Table
    print("Generating cross-target summary...", flush=True)
    cross_data = []
    constructs = global_df['Construct'].unique()
    
    for c in constructs:
        c_subset = global_df[global_df['Construct'] == c]
        c_num = c_subset['Construct Num'].iloc[0]
        
        m2_subset = c_subset[c_subset['Target'] == "MMP2"]
        m3_subset = c_subset[c_subset['Target'] == "MMP3"]
        m9_subset = c_subset[c_subset['Target'] == "MMP9"]
        
        m2_dates = m2_subset['Date'].nunique()
        m3_dates = m3_subset['Date'].nunique()
        m9_dates = m9_subset['Date'].nunique()
        
        m2_mean = m2_subset['Norm Median Ratio'].mean() if not m2_subset.empty else 0
        m3_mean = m3_subset['Norm Median Ratio'].mean() if not m3_subset.empty else 0
        m9_mean = m9_subset['Norm Median Ratio'].mean() if not m9_subset.empty else 0
        
        # Use MMP2 or MMP3 for the ratio
        m23_mean = m2_mean if not m2_subset.empty else m3_mean
        ratio_val = m23_mean / m9_mean if m9_mean > 0 else 0
        
        # Determine "Seen" based on Norm Median Ratio
        overall_max = max(m2_mean, m3_mean, m9_mean)
        if overall_max > 0.6: seen = "High"
        elif overall_max > 0.2: seen = "Medium"
        else: seen = "Low"
        
        cross_data.append({
            'Construct': c,
            'MMP2 Trials': f"{m2_dates} dates" if m2_dates > 0 else "",
            'MMP3 Trials': f"{m3_dates} dates" if m3_dates > 0 else "",
            'MMP9 Trials': f"{m9_dates} dates" if m9_dates > 0 else "",
            'MMP2/9 Ratio': f"{ratio_val:.2f}",
            'Expected': expectations.get(c_num, 'N/A'),
            'Seen': seen
        })
        
    output_cross_csv = os.path.join(output_dir, "cross_target_summary.csv")
    cross_fieldnames = ['Construct', 'MMP2 Trials', 'MMP3 Trials', 'MMP9 Trials', 'MMP2/9 Ratio', 'Expected', 'Seen']
    with open(output_cross_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=cross_fieldnames)
        writer.writeheader()
        writer.writerows(cross_data)
    print(f"Saved cross-target summary to {output_cross_csv}", flush=True)

if __name__ == "__main__":
    main()
