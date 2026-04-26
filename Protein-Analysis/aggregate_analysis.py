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
    output_dir = os.path.join(local_dir, "Aggregate_FCS_Analysis")
    reference_csv = os.path.join(base_dir, "TIMP3_variants_names.csv")
    
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
            
        trial_failed = False
        try:
            with open(csv_path, 'r') as f:
                reader = csv.DictReader(f)
                pos_ctrl_dps = []
                for row in reader:
                    raw_name = row['Filename']
                    if any(ctrl in raw_name.upper() for ctrl in ["POSITIVE CONTROL", "TIMP"]):
                        pos_ctrl_dps.append(float(row.get('Double+ %', 0)))
                
                # If no positive controls exist or max is < 2.0%, fail trial
                if not pos_ctrl_dps or max(pos_ctrl_dps) < 2.0:
                    trial_failed = True
        except Exception as e:
            print(f"  Error checking Pos Ctrls for {csv_path}: {e}", flush=True)
            
        try:
            with open(csv_path, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    raw_name = row['Filename']
                    # Filter out negative controls, keep constructs and positive controls
                    is_nc = any(ctrl in raw_name.upper() for ctrl in ["NC", "NEGATIVE CONTROL"])
                    is_pc = any(ctrl in raw_name.upper() for ctrl in ["POSITIVE CONTROL", "PC"])
                    
                    if is_nc:
                        continue
                    
                    standard_name = standardize_construct(raw_name)
                    
                    # If it's a known PC pattern or "TIMP" (exactly), and not a variant, include it
                    if not standard_name:
                        if is_pc or raw_name.upper() == "TIMP":
                             standard_name = raw_name
                    
                    if standard_name:
                        try:
                            # Extract number for expectation mapping
                            num_match = re.search(r'(\d+)', standard_name)
                            construct_num = num_match.group(1) if num_match else None
                            
                            # Extract required metrics
                            norm_med_ratio = float(row.get('Norm Pos Med Ratio', 0))
                            norm_mean_ratio = float(row.get('Norm Pos Mean Ratio', 0))
                            double_pos = float(row.get('Double+ %', 0))
                            expr_pos = float(row.get('Expr+ %', 0))
                            gated_events = int(float(row.get('Gated Events', 0)))
                            
                            # QC Flags
                            low_expression = expr_pos < 1.0
                            low_events = gated_events < 500
                            
                            # New Filtered Metrics
                            norm_bind_med_expr = float(row.get('Norm Bind Med (Expr+)', 0))
                            norm_bind_mean_expr = float(row.get('Norm Bind Mean (Expr+)', 0))
                            norm_expr_med_bind = float(row.get('Norm Expr Med (Bind+)', 0))
                            
                            all_data.append({
                                'Target': target,
                                'Date': date,
                                'Construct': standard_name,
                                'Construct Num': construct_num,
                                'Raw Name': raw_name,
                                'Norm Median Ratio': norm_med_ratio,
                                'Norm Mean Ratio': norm_mean_ratio,
                                'Double+ %': double_pos,
                                'Norm Bind Med (Expr+)': norm_bind_med_expr,
                                'Norm Bind Mean (Expr+)': norm_bind_mean_expr,
                                'Norm Expr Med (Bind+)': norm_expr_med_bind,
                                'Expr+ %': expr_pos,
                                'Gated Events': gated_events,
                                'Low Expression': low_expression,
                                'Low Events': low_events,
                                'Trial Failed': trial_failed
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
    fieldnames = ['Target', 'Date', 'Construct', 'Raw Name', 'Norm Median Ratio', 'Norm Mean Ratio', 'Double+ %', 'Expr+ %', 'Gated Events', 'Norm Bind Med (Expr+)', 'Norm Bind Mean (Expr+)', 'Norm Expr Med (Bind+)', 'Low Expression', 'Low Events', 'Trial Failed']
    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for d in all_data:
            writer.writerow({k: d[k] for k in fieldnames})
    print(f"Saved aggregate summary to {output_csv}", flush=True)

    # 1. Plotting
    print("Generating plots...", flush=True)
    global_df = pd.DataFrame(all_data)
    
    # Use only valid trials and samples with sufficient expression/events for the generated bar graphs
    global_df_filtered = global_df[(global_df['Trial Failed'] == False) & (global_df['Low Expression'] == False) & (global_df['Low Events'] == False)]
    
    metrics_to_plot = [
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
        }
    ]

    for metric in metrics_to_plot:
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
                    comments.append(f"Excluded {row['Date']} ({row['Target']}) due to invalid Pos Ctrl")
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
    # Pass the filtered global dataframe (valid trials only)
    perform_selectivity_analysis(global_df_filtered, output_dir)

def perform_selectivity_analysis(df, output_dir):
    """
    Identifies constructs tested against multiple targets and generates
    comparison plots and summaries.
    """
    print("Performing selectivity analysis...", flush=True)
    
    selectivity_dir = os.path.join(output_dir, "Selectivity_Analysis")
    if not os.path.exists(selectivity_dir):
        os.makedirs(selectivity_dir)
        print(f"  Created selectivity directory: {selectivity_dir}", flush=True)
        
    # 1. Aggregate data by Construct and Target
    # We use "Norm Median Ratio" as the primary metric
    metric = "Norm Median Ratio"
    
    # Calculate Mean, StdDev, and Count for each (Construct, Target) pair
    stats_df = df.groupby(['Construct', 'Target'])[metric].agg(['mean', 'std', 'count']).reset_index()
    stats_df.columns = ['Construct', 'Target', 'Mean', 'StdDev', 'Count']
    
    # Calculate SEM (Standard Error of the Mean) for error bars
    stats_df['SEM'] = stats_df.apply(lambda row: row['StdDev'] / (row['Count']**0.5) if row['Count'] > 1 else 0, axis=1)
    
    # 2. Identify multi-target constructs
    target_counts = stats_df.groupby('Construct')['Target'].count()
    multi_target_constructs = target_counts[target_counts > 1].index.tolist()
    
    if not multi_target_constructs:
        print("  No constructs found with multiple targets for comparison.", flush=True)
        return

    # Filter stats for these constructs
    selectivity_df = stats_df[stats_df['Construct'].isin(multi_target_constructs)].copy()
    
    # Save CSV
    csv_path = os.path.join(selectivity_dir, "selectivity_summary.csv")
    selectivity_df.to_csv(csv_path, index=False)
    print(f"  Saved selectivity summary to {csv_path}", flush=True)
    
    # 3. Plotting - Grouped Bar Chart
    # Use pandas plotting for grouped bars
    pivot_df = selectivity_df.pivot(index='Construct', columns='Target', values='Mean')
    pivot_sem = selectivity_df.pivot(index='Construct', columns='Target', values='SEM').fillna(0)
    
    if not pivot_df.empty:
        plt.figure(figsize=(16, 10))
        ax = pivot_df.plot(kind='bar', yerr=pivot_sem, figsize=(16, 10), capsize=4, edgecolor='black', alpha=0.8)
        
        plt.title("Selectivity Analysis: Binding Across Targets", fontsize=22, pad=20)
        plt.xlabel("Construct", fontsize=18)
        plt.ylabel(f"Mean {metric}", fontsize=18)
        plt.xticks(rotation=45, ha='right')
        plt.legend(title="Target", fontsize=12, title_fontsize=14, loc='upper left', bbox_to_anchor=(1, 1))
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()
        
        plot_path = os.path.join(selectivity_dir, "selectivity_comparison_all.png")
        plt.savefig(plot_path, dpi=300)
        plt.close()
        print(f"  Saved grouped comparison plot to {plot_path}", flush=True)
    
    # 4. Individual plots for each multi-target construct
    indiv_dir = os.path.join(selectivity_dir, "Individual_Comparisons")
    if not os.path.exists(indiv_dir):
        os.makedirs(indiv_dir)
        
    for construct in multi_target_constructs:
        c_df = selectivity_df[selectivity_df['Construct'] == construct]
        
        plt.figure(figsize=(10, 6))
        # Use a consistent color map for targets
        colors = plt.cm.viridis(np.linspace(0, 0.8, len(c_df)))
        bars = plt.bar(c_df['Target'], c_df['Mean'], yerr=c_df['SEM'], capsize=5, 
                       color=colors, edgecolor='black', alpha=0.9)
        
        plt.title(f"Target Selectivity: {construct}", fontsize=20, pad=15)
        plt.xlabel("Target", fontsize=16)
        plt.ylabel(f"Mean {metric}", fontsize=16)
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
    
    print(f"  Saved {len(multi_target_constructs)} individual comparison plots to {indiv_dir}", flush=True)

if __name__ == "__main__":
    main()
