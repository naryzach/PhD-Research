import os
import sys
import glob
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy import stats

# Filters warnings from pandas/matplotlib
warnings.filterwarnings("ignore")

# ---- INSTALL DEPENDENCIES IF MISSING ----
def install(package):
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

required_packages = ['fcsparser', 'seaborn', 'scipy', 'matplotlib', 'pandas', 'numpy']
for pkg in required_packages:
    try:
        __import__(pkg)
    except ImportError:
        print(f"Installing {pkg}...")
        install(pkg)

import fcsparser

# ---- CONFIGURATION ----
# Paths relative to "DNA Analysis" folder
DATA_DIR = r"../Local/BD6_Flow_Renamed"
OUTPUT_DIR = r"../Local/FCS_Analysis_Results_2"
NEG_CONTROL_PATTERN = "Negative Control" # Starts with this
POS_CONTROL_PATTERNS = ["Positive Control", "TIMP"] # Starts with any of these

# Channels
CH_FSC_A = 'FSC-A'
CH_SSC_A = 'SSC-A'
CH_FSC_H = 'FSC-H'
CH_SSC_H = 'SSC-H'
CH_FITC = 'FITC-A' # Expression
CH_APC = 'APC-A'   # Binding

# Plot settings
SNS_STYLE = "whitegrid"
sns.set_style(SNS_STYLE)

# ---- HELPER FUNCTIONS ----

def load_fcs(file_path):
    """Loads an FCS file and returns the dataframe."""
    try:
        # Suppress warnings from fcsparser if needed
        meta, data = fcsparser.parse(file_path, reformat_meta=True)
        return meta, data
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None, None

def gate_scatter(df, plot_ax=None):
    """
    Gates debris based on FSC-A vs SSC-A.
    """
    # Simple heuristic gate
    mask = (df[CH_FSC_A] > 10000) & (df[CH_SSC_A] > 2000) & \
           (df[CH_FSC_A] < df[CH_FSC_A].max() * 0.99)
    
    if plot_ax:
        # Density plot
        plot_ax.hexbin(df[CH_FSC_A], df[CH_SSC_A], gridsize=100, cmap='inferno', mincnt=1, bins='log')
        plot_ax.set_xlabel("FSC-A")
        plot_ax.set_ylabel("SSC-A")
        plot_ax.set_title("Scatter Gate (Debris Removal)")
    
    return df[mask], mask

def gate_singlets(df, plot_ax=None):
    """
    Gates singlets based on FSC-A vs FSC-H.
    """
    if len(df) == 0:
        return df, []

    ratio = df[CH_FSC_H] / df[CH_FSC_A]
    median_ratio = ratio.median()
    
    # Allow +/- 20% deviation
    mask = (ratio > median_ratio * 0.8) & (ratio < median_ratio * 1.2)
    
    if plot_ax:
        plot_ax.hexbin(df[CH_FSC_A], df[CH_FSC_H], gridsize=100, cmap='inferno', mincnt=1, bins='log')
        # Show gate lines roughly
        x_vals = np.linspace(df[CH_FSC_A].min(), df[CH_FSC_A].max(), 100)
        plot_ax.plot(x_vals, x_vals * median_ratio * 1.2, 'g--', linewidth=1)
        plot_ax.plot(x_vals, x_vals * median_ratio * 0.8, 'g--', linewidth=1)
        plot_ax.set_xlabel("FSC-A")
        plot_ax.set_ylabel("FSC-H")
        plot_ax.set_title("Singlet Gate")
        
    return df[mask], mask

def transform_logicle(data, channels):
    """Simple log10 transform for plotting, handling negatives safely (clipping)."""
    df_trans = data.copy()
    for ch in channels:
        df_trans[ch] = np.log10(np.clip(df_trans[ch], 1, None))
    return df_trans

def main():
    # 1. Setup
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    iso_dir = os.path.join(OUTPUT_DIR, "Publication_Plots")
    if not os.path.exists(iso_dir):
        os.makedirs(iso_dir)

    comp_dir = os.path.join(OUTPUT_DIR, "Comparisons_vs_PosCtrl")
    if not os.path.exists(comp_dir):
        os.makedirs(comp_dir)
        
    # Find files
    search_path = os.path.join(os.getcwd(), DATA_DIR)
    
    print(f"Searching in: {search_path}")
    all_files = glob.glob(os.path.join(search_path, "*.fcs"))
    
    if not all_files:
        print("No FCS files found!")
        return

    # Semantic Channel Mapping
    CH_EXPR = CH_FITC # Expression (Y-axis)
    CH_BIND = CH_APC  # Binding (X-axis)

    # 2. Process Controls
    # Negative Controls
    neg_files = [f for f in all_files if os.path.basename(f).startswith(NEG_CONTROL_PATTERN)]
    print(f"Found {len(neg_files)} negative control files.")
    
    neg_dfs = []
    for f in neg_files:
        _, d = load_fcs(f)
        if d is not None:
            d_scat, _ = gate_scatter(d)
            d_sing, _ = gate_singlets(d_scat)
            neg_dfs.append(d_sing)
            
    if neg_dfs:
        neg_concat = pd.concat(neg_dfs)
        thresh_expr = np.percentile(neg_concat[CH_EXPR], 99)
        thresh_bind = np.percentile(neg_concat[CH_BIND], 99)
        
        neg_mfi_expr = neg_concat[CH_EXPR].median()
        neg_mfi_bind = neg_concat[CH_BIND].median()
        neg_rsd_bind = stats.median_abs_deviation(neg_concat[CH_BIND], scale='normal')
    else:
        print("Warning: No negative controls found matching pattern. Using default.")
        thresh_expr = 1000
        thresh_bind = 1000
        neg_mfi_expr = 1
        neg_mfi_bind = 1
        neg_rsd_bind = 1
        neg_concat = None

    # Positive Controls
    pos_files = [f for f in all_files if any(os.path.basename(f).startswith(p) for p in POS_CONTROL_PATTERNS)]
    print(f"Found {len(pos_files)} positive control files (Patterns: {POS_CONTROL_PATTERNS}).")
    
    pos_dfs = []
    for f in pos_files:
        _, d = load_fcs(f)
        if d is not None:
            d_scat, _ = gate_scatter(d)
            d_sing, _ = gate_singlets(d_scat)
            pos_dfs.append(d_sing)
            
    pos_mfi_bind = 1.0
    pos_concat = None
    if pos_dfs:
        pos_concat = pd.concat(pos_dfs)
        pos_mfi_expr = pos_concat[CH_EXPR].median()
        pos_mfi_bind = pos_concat[CH_BIND].median()
        print(f"Pos Ctrl Stats: Binding (APC) Median={pos_mfi_bind:.2f}")

    # 3. Analyze All Samples
    summary_stats = []
    ridge_data = [] 
    
    # Pre-calculate log data for controls for plotting (Binding = APC)
    neg_log_bind = np.log10(np.clip(neg_concat[CH_BIND].sample(n=min(len(neg_concat), 5000)), 1, None)) if neg_concat is not None else None
    pos_log_bind = np.log10(np.clip(pos_concat[CH_BIND].sample(n=min(len(pos_concat), 5000)), 1, None)) if pos_concat is not None else None

    for f in all_files:
        fname = os.path.basename(f)
        clean_name = os.path.splitext(fname)[0]
        
        print(f"Processing {clean_name}...")
        
        meta, df = load_fcs(f)
        if df is None: continue
        
        # Gates
        df_scat, _ = gate_scatter(df)
        df_sing, _ = gate_singlets(df_scat)
        
        # Stats
        count_total = len(df)
        count_gated = len(df_sing)
        
        if count_gated == 0:
            print(f"Warning: No events gated for {clean_name}")
            continue

        # Positivity
        expr_pos = df_sing[CH_EXPR] > thresh_expr
        bind_pos = df_sing[CH_BIND] > thresh_bind
        double_pos = expr_pos & bind_pos
        
        pct_expr = expr_pos.mean() * 100
        pct_bind = bind_pos.mean() * 100
        pct_double = double_pos.mean() * 100
        
        mfi_expr = df_sing[CH_EXPR].median()
        mfi_bind = df_sing[CH_BIND].median()
        
        # Formatting Ridge Data (Binding = APC)
        sub = df_sing.sample(n=min(len(df_sing), 2000))
        t_vals = np.log10(np.clip(sub[CH_BIND], 1, None))
        for v in t_vals:
            ridge_data.append({"Sample": clean_name, "LogBinding": v})

        # Advanced Stats
        fc_expr = mfi_expr / max(neg_mfi_expr, 1)
        fc_bind = mfi_bind / max(neg_mfi_bind, 1)
        
        # % of Pos Ctrl (Binding)
        pct_of_pos_mfi = (mfi_bind / pos_mfi_bind) * 100 if pos_mfi_bind > 0 else 0
        
        # Binding Efficiency (Binding normalized to Expression = Binding / Expression)
        # Avoid div by zero
        binding_eff = mfi_bind / mfi_expr if mfi_expr > 1 else 0
        
        si_bind = (mfi_bind - neg_mfi_bind) / (2 * neg_rsd_bind) if neg_rsd_bind > 0 else 0
        
        summary_stats.append({
            "Filename": clean_name,
            "Total Events": count_total,
            "Gated Events": count_gated,
            "% Gated": (count_gated/count_total)*100,
            "Expr+ %": pct_expr,
            "Bind+ %": pct_bind,
            "Double+ %": pct_double,
            "Expr MFI": mfi_expr,
            "Bind MFI": mfi_bind,
            "Expr Fold Change": fc_expr,
            "Bind Fold Change": fc_bind,
            "Bind Stain Index": si_bind,
            "% of Pos Ctrl": pct_of_pos_mfi,
            "Binding Efficiency": binding_eff
        })
        
        # --- PLOTTING ---
        
        # 1. Main 4-panel Plot
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Scatter
        gate_scatter(df, plot_ax=axes[0,0])
        gate_singlets(df_scat, plot_ax=axes[0,1])
        
        # Density (Expression Y vs Binding X) - User requested FITC on Y (Expr), APC on X (Bind)
        t_expr = np.log10(np.clip(df_sing[CH_EXPR], 1, None))
        t_bind = np.log10(np.clip(df_sing[CH_BIND], 1, None))
        
        # Plotting: x=t_bind (APC), y=t_expr (FITC)
        axes[1,0].hexbin(t_bind, t_expr, gridsize=100, cmap='jet', mincnt=1, bins='log')
        axes[1,0].axvline(np.log10(max(1, thresh_bind)), color='k', linestyle='--')
        axes[1,0].axhline(np.log10(max(1, thresh_expr)), color='k', linestyle='--')
        axes[1,0].set_xlabel("Log10 APC (Binding)")
        axes[1,0].set_ylabel("Log10 FITC (Expression)")
        axes[1,0].set_title(f"Double+: {pct_double:.1f}%")
        
        # Histogram Overlay (vs Neg) - Binding (APC)
        sns.kdeplot(t_bind, fill=True, ax=axes[1,1], color='g', label='Sample')
        if neg_log_bind is not None:
            sns.kdeplot(neg_log_bind, fill=True, ax=axes[1,1], color='gray', alpha=0.3, label='Neg Ctrl')
        axes[1,1].axvline(np.log10(max(1, thresh_bind)), color='r', linestyle='--')
        axes[1,1].set_title(f"Binding Dist (FC: {fc_bind:.1f}x)")
        axes[1,1].legend()
        
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, f"{clean_name}_analysis.png"))
        plt.close()
        
        # 2. Publication Plot (Expr Y vs Bind X)
        plt.figure(figsize=(6, 6))
        plt.hexbin(t_bind, t_expr, gridsize=100, cmap='jet', mincnt=1, bins='log')
        plt.axvline(np.log10(max(1, thresh_bind)), color='k', linestyle='--')
        plt.axhline(np.log10(max(1, thresh_expr)), color='k', linestyle='--')
        plt.xlabel("Log10 APC-A (Binding)")
        plt.ylabel("Log10 FITC-A (Expression)")
        plt.title(f"{clean_name}")
        
        # Add a text box for statistics (cleaner than title)
        textstr = '\n'.join((
            f'Double+: {pct_double:.1f}%',
            f'Bind+: {pct_bind:.1f}%',
            f'Expr+: {pct_expr:.1f}%'
        ))
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        plt.gca().text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
                verticalalignment='top', bbox=props)
                
        plt.tight_layout()
        plt.savefig(os.path.join(iso_dir, f"{clean_name}_iso.png"), dpi=200)
        plt.close()
        
        # 3. Comparison vs Positive Control (Histogram - Binding)
        plt.figure(figsize=(8, 5))
        if neg_log_bind is not None:
            sns.kdeplot(neg_log_bind, fill=True, color='lightgray', label='Neg Ctrl', alpha=0.5)
        if pos_log_bind is not None:
            # Line only for positive control
            sns.kdeplot(pos_log_bind, color='blue', label='Pos Ctrl', linestyle='--', linewidth=2)
            
        sns.kdeplot(t_bind, fill=True, color='green', label=clean_name, alpha=0.4)
        
        plt.xlabel("Log10 APC-A (Binding)")
        plt.title(f"Binding Comparison: {clean_name}\n(% of Pos: {pct_of_pos_mfi:.1f}%)")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(comp_dir, f"Comp_vs_Pos_{clean_name}.png"))
        plt.close()

    # Save Stats
    df_stats = pd.DataFrame(summary_stats)
    df_stats.to_csv(os.path.join(OUTPUT_DIR, "summary_stats.csv"), index=False)
    print(f"\nAnalysis complete. Saved to {os.path.join(OUTPUT_DIR, 'summary_stats.csv')}")
    
    # ---- AGGREGATE PLOTS ----
    generate_aggregate_plots(df_stats, ridge_data, thresh_bind, OUTPUT_DIR)
    
    # Report
    generate_report(df_stats)

def generate_aggregate_plots(df_stats, ridge_data, thresh_bind, output_dir):
    # Re-implementing the aggregation plots to ensure they exist
    print("Generating aggregate plots...")
    
    if not ridge_data:
        return

    df_ridge = pd.DataFrame(ridge_data)

    # 0. Ridgeline Plot (Binding) - Grouped by Prefix
    # Identify Groups (First word of filename)
    df_ridge['Group'] = df_ridge['Sample'].apply(lambda x: x.split()[0] if ' ' in x else 'Misc')
    groups = df_ridge['Group'].unique()
    
    for group in groups:
        print(f"Generating Ridgeline for Group: {group}...")
        df_group = df_ridge[df_ridge['Group'] == group]
        
        # Filter Middle 95% PER SAMPLE
        df_filtered = pd.DataFrame()
        for sample in df_group['Sample'].unique():
            d_sample = df_group[df_group['Sample'] == sample]
            low = d_sample['LogBinding'].quantile(0.025)
            high = d_sample['LogBinding'].quantile(0.975)
            d_filt = d_sample[(d_sample['LogBinding'] >= low) & (d_sample['LogBinding'] <= high)]
            df_filtered = pd.concat([df_filtered, d_filt])
        
        if df_filtered.empty: continue

        # Sort within group
        group_stats = df_stats[df_stats['Filename'].isin(df_filtered['Sample'].unique())]
        order_mfi = group_stats.sort_values("Bind MFI", ascending=False)["Filename"].tolist()
        
        # Dynamic X-limits
        x_min = df_filtered["LogBinding"].min()
        x_max = df_filtered["LogBinding"].max()
        x_pad = (x_max - x_min) * 0.1
        x_lims = (x_min - x_pad, x_max + x_pad*2)

        pal = sns.cubehelix_palette(len(order_mfi), rot=-.25, light=.7)
        
        g = sns.FacetGrid(df_filtered, row="Sample", hue="Sample", aspect=10, height=0.8, palette=pal, row_order=order_mfi)
        g.map(sns.kdeplot, "LogBinding", clip_on=False, fill=True, alpha=1, lw=1.5, bw_adjust=.2)
        g.map(sns.kdeplot, "LogBinding", clip_on=False, color="w", lw=2, bw_adjust=.2)
        g.map(plt.axhline, y=0, lw=2, clip_on=False)
        
        def label(x, color, label):
            ax = plt.gca()
            ax.text(0, .2, label, fontweight="bold", color=color,
                    ha="left", va="center", transform=ax.transAxes, fontsize=10)
            ax.set_facecolor('none')
        
        g.map(label, "LogBinding")
        
        g.set(xlim=x_lims)
        g.fig.subplots_adjust(hspace=-.5)
        g.set_titles("")
        g.set(yticks=[])
        g.despine(bottom=True, left=True)
        g.set_xlabels(f"Log10 APC-A (Binding) - Group {group}")
        
        t_thresh = np.log10(max(1, thresh_bind))
        for ax in g.axes.flatten():
            ax.axvline(t_thresh, color='r', alpha=0.5, linestyle='--')
            
        out_name = f"Aggregate_Ridgeline_Binding_{group}.png"
        g.savefig(os.path.join(output_dir, out_name))
        plt.close()

    # 1. Ridgeline (All Samples)
    out_name = os.path.join(output_dir, "Aggregate_Ridgeline_All.png")
    
    # Sort by MFI
    order = df_stats.sort_values("Bind MFI", ascending=False)["Filename"].tolist()
    pal = sns.cubehelix_palette(len(order), rot=-.25, light=.7)
    
    g = sns.FacetGrid(df_ridge, row="Sample", hue="Sample", aspect=15, height=0.6, palette=pal, row_order=order)
    g.map(sns.kdeplot, "LogBinding", clip_on=False, fill=True, alpha=1, lw=1.5, bw_adjust=.2)
    g.map(sns.kdeplot, "LogBinding", clip_on=False, color="w", lw=2, bw_adjust=.2)
    g.map(plt.axhline, y=0, lw=2, clip_on=False)
    
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color,
                ha="left", va="center", transform=ax.transAxes, fontsize=9)
        ax.set_facecolor('none')
        
    g.map(label, "LogBinding")
    
    g.set_titles("")
    g.set(yticks=[])
    g.despine(bottom=True, left=True)
    g.fig.subplots_adjust(hspace=-.5)
    g.set_xlabels("Log10 APC-A (Binding)")
    
    # Ref Line
    t_thresh = np.log10(max(1, thresh_bind))
    for ax in g.axes.flatten():
        ax.axvline(t_thresh, color='r', alpha=0.5, linestyle='--')
        
    g.savefig(out_name)
    plt.close()

    # 2. Fold Change Bar Plot
    plt.figure(figsize=(12, 6))
    # Sort for bar plot
    df_sorted_fc = df_stats.sort_values("Bind Fold Change", ascending=False)
    sns.barplot(data=df_sorted_fc, x="Filename", y="Bind Fold Change", palette="viridis")
    plt.xticks(rotation=90)
    plt.axhline(1, color='k', linestyle='--', label='No Change')
    plt.title("Binding Fold Change (vs Neg Control)")
    plt.ylabel("Fold Change (Sample MFI / Neg MFI)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Aggregate_FoldChange_Bar.png"))
    plt.close()

    # 3. Stain Index Bar Plot
    plt.figure(figsize=(12, 6))
    df_sorted_si = df_stats.sort_values("Bind Stain Index", ascending=False)
    sns.barplot(data=df_sorted_si, x="Filename", y="Bind Stain Index", palette="magma")
    plt.xticks(rotation=90)
    plt.title("Binding Stain Index (Separation Quality)")
    plt.ylabel("Stain Index (Delta Median / 2*rSD)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Aggregate_StainIndex_Bar.png"))
    plt.close()
    
    # 4. Bar plot of Double Positives
    plt.figure(figsize=(12, 6))
    sns.barplot(data=df_stats, x="Filename", y="Double+ %", palette="Blues_d")
    plt.xticks(rotation=90)
    plt.title("Double Positive Population % by Sample")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Aggregate_DoublePos_Bar.png"))
    plt.close()
    
    # 5. Scatter MFI
    # User requested: APC on X (Binding), FITC on Y (Expr)
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=df_stats, x="Bind MFI", y="Expr MFI", hue="Filename", s=100)
    plt.title("Mean Fluorescence Intensity: Binding vs Expression")
    plt.xlabel("Binding MFI (APC)")
    plt.ylabel("Expression MFI (FITC)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Aggregate_MFI_Scatter.png"))
    plt.close()

def generate_report(df):
    """Generates a text report summarizing the findings."""
    
    report_path = os.path.join(OUTPUT_DIR, "Analysis_Report.md")
    
    # Identify top performers
    top_binding = df.sort_values("Bind MFI", ascending=False).iloc[0]
    top_expr = df.sort_values("Expr MFI", ascending=False).iloc[0]
    
    # Group analysis
    df['Group'] = df['Filename'].apply(lambda x: x.split()[0] if ' ' in x else 'Misc')
    groups = df['Group'].unique()
    
    with open(report_path, "w") as f:
        f.write("# FCS Analysis Report\n\n")
        f.write(f"**Total Samples**: {len(df)}\n")
        f.write(f"**Date**: {pd.Timestamp.now().strftime('%Y-%m-%d')}\n")
        f.write(f"**Average Gating Efficiency**: {df['% Gated'].mean():.1f}%\n\n")
        
        f.write("## 1. Executive Summary\n")
        f.write(f"- **Highest Binding (APC)**: {top_binding['Filename']} (MFI: {top_binding['Bind MFI']:.0f}, %+: {top_binding['Bind+ %']:.1f}%)\n")
        f.write(f"- **Highest Expression (FITC)**: {top_expr['Filename']} (MFI: {top_expr['Expr MFI']:.0f}, %+: {top_expr['Expr+ %']:.1f}%)\n")
        f.write("\n")
        
        f.write("## 2. Group Performance\n")
        for g in groups:
            sub = df[df['Group'] == g]
            if len(sub) == 0: continue
            avg_fc_bind = sub['Bind Fold Change'].mean()
            avg_fc_expr = sub['Expr Fold Change'].mean()
            avg_eff = sub['Binding Efficiency'].mean()
            avg_pct_pos = sub['% of Pos Ctrl'].mean()
            avg_si = sub['Bind Stain Index'].mean()
            
            f.write(f"### Group {g}\n")
            f.write(f"- **Samples**: {len(sub)}\n")
            f.write(f"- **Avg Binding Fold Change**: {avg_fc_bind:.2f}x\n")
            f.write(f"- **Avg Expression Fold Change**: {avg_fc_expr:.2f}x\n")
            f.write(f"- **Avg Stain Index**: {avg_si:.2f}\n")
            f.write(f"- **Avg Binding Efficiency (Binding / Expression ratio)**: {avg_eff:.2f}\n")
            f.write(f"- **Avg % of Pos Ctrl**: {avg_pct_pos:.1f}%\n")
            f.write(f"- **Observation**: {stats_observation(avg_fc_bind, avg_fc_expr, avg_eff)}\n\n")

def stats_observation(fc_bind, fc_expr, efficiency):
    obs = []
    
    # Binding Analysis
    if fc_bind < 1.1:
        obs.append("No significant binding detected.")
    elif fc_bind < 1.5:
        obs.append("Weak/Minimal binding signal.")
    elif fc_bind < 5.0:
        obs.append("Moderate binding observed.")
    else:
        obs.append("Strong binding signal detected.")
        
    # Expression Analysis
    if fc_expr < 1.1:
        obs.append("Expression is at background levels (Low display level).")
    elif fc_expr < 1.5:
        obs.append("Low expression detected.")
    elif fc_expr < 3.0:
        obs.append("Moderate expression levels.")
    else:
        obs.append("High surface expression.")

    # Efficiency/Context
    if fc_bind > 1.5 and fc_expr > 1.5:
        # Both present
        if efficiency > 1.0: # Arbitrary high efficiency threshold, depends on gain settings
            obs.append("High binding relative to expression (Potentially high affinity).")
        elif efficiency < 0.2:
            obs.append("Binding is low relative to expression (Potentially low affinity).")
    elif fc_bind > 1.5 and fc_expr < 1.2:
        obs.append("Binding detected despite low expression (Type II error or non-specific?).")
    elif fc_bind < 1.2 and fc_expr > 2.0:
        obs.append("High expression but no binding (Non-binder).")
        
    return " ".join(obs)

if __name__ == "__main__":
    main()
