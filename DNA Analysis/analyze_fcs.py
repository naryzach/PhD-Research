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
DATA_DIR = r"../Local/TwistBio_FCS"
OUTPUT_DIR = r"../Local/FCS_Analysis_Results"
NEG_CONTROL_PATTERN = "Negative Control" # Starts with this

# Channels
CH_FSC_A = 'FSC-A'
CH_SSC_A = 'SSC-A'
CH_FSC_H = 'FSC-H'
CH_SSC_H = 'SSC-H'
CH_FITC = 'FITC-A' # Binding
CH_APC = 'APC-A'   # Expression

# Plot settings
SNS_STYLE = "whitegrid"
sns.set_style(SNS_STYLE)

# ---- HELPER FUNCTIONS ----

def load_fcs(file_path):
    """Loads an FCS file and returns the dataframe."""
    try:
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
        
    # Find files
    search_path = os.path.join(os.getcwd(), DATA_DIR)
    
    print(f"Searching in: {search_path}")
    all_files = glob.glob(os.path.join(search_path, "*.fcs"))
    
    if not all_files:
        print("No FCS files found!")
        return

    # 2. Process Negative Controls
    # Strict matching to start with "Negative Control" (case sensitive usually, but Windows matches anyway)
    # This avoids "0 Negative Control" if we check startswith
    neg_files = [f for f in all_files if os.path.basename(f).startswith(NEG_CONTROL_PATTERN)]
    print(f"Found {len(neg_files)} negative control files (Pattern: '{NEG_CONTROL_PATTERN}').")
    
    neg_dfs = []
    
    for f in neg_files:
        _, d = load_fcs(f)
        if d is not None:
            d_scat, _ = gate_scatter(d)
            d_sing, _ = gate_singlets(d_scat)
            neg_dfs.append(d_sing)
            
    if neg_dfs:
        neg_concat = pd.concat(neg_dfs)
        thresh_fitc = np.percentile(neg_concat[CH_FITC], 99)
        thresh_apc = np.percentile(neg_concat[CH_APC], 99)
        print(f"Thresholds (99% Neg): FITC={thresh_fitc:.2f}, APC={thresh_apc:.2f}")
    else:
        print("Warning: No negative controls found matching pattern. Using default.")
        thresh_fitc = 1000
        thresh_apc = 1000

    # Calculate Negative Control Stats for Fold Change / Stain Index
    if neg_dfs:
        # Use simple median for MFI
        neg_mfi_fitc = neg_concat[CH_FITC].median()
        neg_mfi_apc = neg_concat[CH_APC].median()
        
        # Robust SD for Stain Index (MAD * 1.4826)
        # Handle zeros/negatives by shifting if necessary, but usually distinct population is fine
        neg_rsd_fitc = stats.median_abs_deviation(neg_concat[CH_FITC], scale='normal')
        neg_rsd_apc = stats.median_abs_deviation(neg_concat[CH_APC], scale='normal')
        
        print(f"Neg Ctrl Stats: FITC Median={neg_mfi_fitc:.2f}, rSD={neg_rsd_fitc:.2f}")
    else:
        neg_mfi_fitc = 1
        neg_mfi_apc = 1
        neg_rsd_fitc = 1
        neg_rsd_apc = 1
        print("Warning: Using default Neg Stats (1.0)")

    # 3. Analyze All Samples
    summary_stats = []
    
    # Store subsampled data for Ridgeline plot
    # DataFrame with columns: [Value, SampleName]
    ridge_data = [] 
    
    for f in all_files:
        fname = os.path.basename(f)
        # Clean name for display/saving (remove .fcs extension)
        clean_name = os.path.splitext(fname)[0]
        print(f"Processing {clean_name}...")
        
        meta, df = load_fcs(f)
        if df is None: continue
        
        # Gates
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # 1. Scatter Gate
        df_scat, _ = gate_scatter(df, plot_ax=axes[0,0])
        
        # 2. Singlet Gate
        df_sing, _ = gate_singlets(df_scat, plot_ax=axes[0,1])
        
        # Subsample for Ridgeline (Binding - FITC)
        # Log transform for layout
        if len(df_sing) > 0:
            sub = df_sing.sample(n=min(len(df_sing), 2000))
            t_vals = np.log10(np.clip(sub[CH_FITC], 1, None))
            for v in t_vals:
                ridge_data.append({"Sample": clean_name, "LogFITC": v})

        # Stats
        count_total = len(df)
        count_gated = len(df_sing)
        
        # Positivity
        fitc_pos = df_sing[CH_FITC] > thresh_fitc
        apc_pos = df_sing[CH_APC] > thresh_apc
        double_pos = fitc_pos & apc_pos
        
        pct_fitc = fitc_pos.mean() * 100
        pct_apc = apc_pos.mean() * 100
        pct_double = double_pos.mean() * 100
        
        mfi_fitc = df_sing[CH_FITC].median()
        mfi_apc = df_sing[CH_APC].median()
        
        # Advanced Stats
        # Fold Change (FC) = MFI_Sample / MFI_Neg
        fc_fitc = mfi_fitc / max(neg_mfi_fitc, 1) # Avoid div by zero
        fc_apc = mfi_apc / max(neg_mfi_apc, 1)
        
        # Stain Index (SI) = (MFI_Pos - MFI_Neg) / (2 * rSD_Neg) 
        # Here we treat the whole sample median as MFI_Pos for the metric, 
        # or we could use the MFI of the positive population. 
        # Standard "Stain Index" usually compares positive peak vs negative peak.
        # We'll use Sample Median vs Neg Median to show overall shift.
        si_fitc = (mfi_fitc - neg_mfi_fitc) / (2 * neg_rsd_fitc) if neg_rsd_fitc > 0 else 0
        
        summary_stats.append({
            "Filename": clean_name,
            "Total Events": count_total,
            "Gated Events": count_gated,
            "% Gated": (count_gated/count_total)*100 if count_total > 0 else 0,
            "FITC+ %": pct_fitc,
            "APC+ %": pct_apc,
            "Double+ %": pct_double,
            "FITC MFI": mfi_fitc,
            "APC MFI": mfi_apc,
            "FITC Fold Change": fc_fitc,
            "FITC Stain Index": si_fitc
        })
        
        # 3. Fluorescence Plot (Density)
        # Log transform
        t_fitc = np.log10(np.clip(df_sing[CH_FITC], 1, None))
        t_apc = np.log10(np.clip(df_sing[CH_APC], 1, None))
        t_thresh_fitc = np.log10(max(1, thresh_fitc))
        t_thresh_apc = np.log10(max(1, thresh_apc))
        
        # Use Hexbin for density
        # SWAPPED: X = Expression (APC), Y = Binding (FITC)
        hb = axes[1,0].hexbin(t_apc, t_fitc, gridsize=100, cmap='jet', mincnt=1, bins='log')
        # Add quadrants (Lines swapped too: Vline=APC thresh, Hline=FITC thresh)
        axes[1,0].axvline(t_thresh_apc, color='black', linestyle='--')
        axes[1,0].axhline(t_thresh_fitc, color='black', linestyle='--')
        
        axes[1,0].set_xlabel("Log10 APC-A (Expression)")
        axes[1,0].set_ylabel("Log10 FITC-A (Binding)")
        axes[1,0].set_title(f"Binding vs Expression (Density)\nDouble+: {pct_double:.1f}%")
        # Can add colorbar if needed, but might clutter 4-panel
        
        # 4. Histogram Overlay (Binding focus) - FITC is now Y in main plot, but usually histogram is 1D. Keep Binding (FITC).
        sns.kdeplot(t_fitc, fill=True, ax=axes[1,1], color='g', label='Sample')
        if neg_dfs:
            neg_sing_sample = neg_concat.sample(n=min(len(neg_concat), 5000))
            t_neg_fitc = np.log10(np.clip(neg_sing_sample[CH_FITC], 1, None))
            sns.kdeplot(t_neg_fitc, fill=True, ax=axes[1,1], color='gray', alpha=0.3, label='Neg Ctrl')
            
        axes[1,1].axvline(t_thresh_fitc, color='r', linestyle='--')
        axes[1,1].set_xlabel("Log10 FITC-A")
        axes[1,1].set_title(f"Binding Dist (FC: {fc_fitc:.1f}x)")
        axes[1,1].legend()
        
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, f"{clean_name}_analysis.png"))
        plt.close()

        # 5. Isolated Publication Plot (Binding vs Expression)
        # Clean, single figure, larger fonts
        plt.figure(figsize=(8, 7))
        # Hexbin with good contrast
        hb = plt.hexbin(t_apc, t_fitc, gridsize=100, cmap='jet', mincnt=1, bins='log')
        
        # Add Gates (Thick black dashed lines)
        plt.axvline(t_thresh_apc, color='black', linestyle='--', linewidth=1.5)
        plt.axhline(t_thresh_fitc, color='black', linestyle='--', linewidth=1.5)
        
        # Labels and Title
        plt.xlabel("Log10 APC-A (Expression)", fontsize=12, fontweight='bold')
        plt.ylabel("Log10 FITC-A (Binding)", fontsize=12, fontweight='bold')
        plt.title(f"{clean_name} Analysis", fontsize=14, fontweight='bold')
        
        # Add a text box for statistics (cleaner than title)
        textstr = '\n'.join((
            f'Double+: {pct_double:.1f}%',
            f'FITC+: {pct_fitc:.1f}%',
            f'APC+: {pct_apc:.1f}%'
        ))
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        plt.gca().text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
                verticalalignment='top', bbox=props)
        
        plt.tight_layout()
        plt.savefig(os.path.join(iso_dir, f"{clean_name}_iso_binding_vs_expression.png"), dpi=300)
        plt.close()

    # Save Stats
    df_stats = pd.DataFrame(summary_stats)
    df_stats.to_csv(os.path.join(OUTPUT_DIR, "summary_stats.csv"), index=False)
    print(f"\nAnalysis complete. Saved to {os.path.join(OUTPUT_DIR, 'summary_stats.csv')}")
    
    # ---- AGGREGATE PLOTS ----
    
    # 1. Ridgeline Plot (Binding) - Grouped by Prefix
    if ridge_data:
        df_ridge = pd.DataFrame(ridge_data)
        
        # Identify Groups (First word of filename)
        # Exclude Control if present in ridge_data (though usually we skip it in the loop or it has unique name)
        df_ridge['Group'] = df_ridge['Sample'].apply(lambda x: x.split()[0] if ' ' in x else 'Misc')
        groups = df_ridge['Group'].unique()
        
        for group in groups:
            print(f"Generating Ridgeline for Group: {group}...")
            df_group = df_ridge[df_ridge['Group'] == group]
            
            # Filter Middle 95% PER SAMPLE to cut tails
            # This makes the distribution tighter and removes outliers
            df_filtered = pd.DataFrame()
            for sample in df_group['Sample'].unique():
                d_sample = df_group[df_group['Sample'] == sample]
                low = d_sample['LogFITC'].quantile(0.025)
                high = d_sample['LogFITC'].quantile(0.975)
                d_filt = d_sample[(d_sample['LogFITC'] >= low) & (d_sample['LogFITC'] <= high)]
                df_filtered = pd.concat([df_filtered, d_filt])
            
            if df_filtered.empty: continue

            # Sort within group
            # Get subset of summary stats for this group to sort
            group_stats = df_stats[df_stats['Filename'].isin(df_filtered['Sample'].unique())]
            order_mfi = group_stats.sort_values("FITC MFI", ascending=False)["Filename"].tolist()
            
            # Dynamic X-limits based on filtered data
            x_min = df_filtered["LogFITC"].min()
            x_max = df_filtered["LogFITC"].max()
            x_pad = (x_max - x_min) * 0.1
            x_lims = (x_min - x_pad, x_max + x_pad*2) # Extra right padding for labels

            # Setup Palette
            pal = sns.cubehelix_palette(len(order_mfi), rot=-.25, light=.7)
            
            # Setup Grid
            # Height and Aspect tuned for visibility
            g = sns.FacetGrid(df_filtered, row="Sample", hue="Sample", aspect=10, height=0.8, palette=pal, row_order=order_mfi)
            
            # Densities
            # fill=True (replacement for shade), alpha=1 makes the density opaque to hide potential grid lines behind it
            # But we must make the AXES background transparent so the rectangle doesn't hide the plot above
            g.map(sns.kdeplot, "LogFITC", clip_on=False, fill=True, alpha=1, lw=1.5, bw_adjust=.5)
            g.map(sns.kdeplot, "LogFITC", clip_on=False, color="w", lw=2, bw_adjust=.5)
            g.map(plt.axhline, y=0, lw=2, clip_on=False)
            
            # Labels
            def label(x, color, label):
                ax = plt.gca()
                ax.text(0, .2, label, fontweight="bold", color=color,
                        ha="left", va="center", transform=ax.transAxes, fontsize=10)
                # CRITICAL FIX: Make the background transparent so it doesn't cover the plot above it
                ax.set_facecolor('none')
            
            g.map(label, "LogFITC")
            
            # Adjust Layout
            g.set(xlim=x_lims)
            g.fig.subplots_adjust(hspace=-.5) # Overlap amount
            g.set_titles("")
            g.set(yticks=[])
            g.despine(bottom=True, left=True)
            g.set_xlabels(f"Log10 FITC-A (Binding) - Group {group}")
            
            # Ref Line
            t_thresh = np.log10(max(1, thresh_fitc))
            for ax in g.axes.flatten():
                ax.axvline(t_thresh, color='r', alpha=0.5, linestyle='--')
                # Ensure visibility of the density on top of the grid
                # But since we remove axes background, the grid might disappear or be weird.
                # Actually, despine(left=True) removes Y axis. 
                # We want the density (alpha=1) to cover the one behind it, which works by draw order.
                
            out_name = f"Aggregate_Ridgeline_Binding_{group}.png"
            g.savefig(os.path.join(OUTPUT_DIR, out_name))
            plt.close()

        # 1b. Ridgeline Plot - ALL SAMPLES
        print("Generating Ridgeline for ALL samples...")
        # Filter Middle 95% PER SAMPLE (for all samples)
        df_filtered_all = pd.DataFrame()
        for sample in df_ridge['Sample'].unique():
            d_sample = df_ridge[df_ridge['Sample'] == sample]
            low = d_sample['LogFITC'].quantile(0.025)
            high = d_sample['LogFITC'].quantile(0.975)
            d_filt = d_sample[(d_sample['LogFITC'] >= low) & (d_sample['LogFITC'] <= high)]
            df_filtered_all = pd.concat([df_filtered_all, d_filt])
        
        if not df_filtered_all.empty:
            # Sort all
            order_mfi_all = df_stats.sort_values("FITC MFI", ascending=False)["Filename"].tolist()
            
            # Dynamic X-limits
            x_min = df_filtered_all["LogFITC"].min()
            x_max = df_filtered_all["LogFITC"].max()
            x_pad = (x_max - x_min) * 0.1
            x_lims = (x_min - x_pad, x_max + x_pad*2)

            pal = sns.cubehelix_palette(len(order_mfi_all), rot=-.25, light=.7)
            
            # Taller aspect for many samples
            g = sns.FacetGrid(df_filtered_all, row="Sample", hue="Sample", aspect=15, height=0.6, palette=pal, row_order=order_mfi_all)
            
            g.map(sns.kdeplot, "LogFITC", clip_on=False, fill=True, alpha=1, lw=1.5, bw_adjust=.5)
            g.map(sns.kdeplot, "LogFITC", clip_on=False, color="w", lw=2, bw_adjust=.5)
            g.map(plt.axhline, y=0, lw=2, clip_on=False)
            
            def label(x, color, label):
                ax = plt.gca()
                ax.text(0, .2, label, fontweight="bold", color=color,
                        ha="left", va="center", transform=ax.transAxes, fontsize=9)
                ax.set_facecolor('none') # Transparency fix
            
            g.map(label, "LogFITC")
            
            g.set(xlim=x_lims)
            g.fig.subplots_adjust(hspace=-.6) # Tighter overlap for many samples
            g.set_titles("")
            g.set(yticks=[])
            g.despine(bottom=True, left=True)
            g.set_xlabels("Log10 FITC-A (Binding) - All Samples")
            
            t_thresh = np.log10(max(1, thresh_fitc))
            for ax in g.axes.flatten():
                ax.axvline(t_thresh, color='r', alpha=0.5, linestyle='--')
                
            g.savefig(os.path.join(OUTPUT_DIR, "Aggregate_Ridgeline_Binding_All.png"))
            plt.close()

    # 2. Fold Change Bar Plot
    plt.figure(figsize=(12, 6))
    # Sort for bar plot
    df_sorted = df_stats.sort_values("FITC Fold Change", ascending=False)
    sns.barplot(data=df_sorted, x="Filename", y="FITC Fold Change", palette="viridis")
    plt.xticks(rotation=90)
    plt.axhline(1, color='k', linestyle='--', label='No Change')
    plt.title("Binding Fold Change (vs Neg Control)")
    plt.ylabel("Fold Change (Sample MFI / Neg MFI)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "Aggregate_FoldChange_Bar.png"))
    plt.close()

    # 3. Stain Index Bar Plot
    plt.figure(figsize=(12, 6))
    df_sorted = df_stats.sort_values("FITC Stain Index", ascending=False)
    sns.barplot(data=df_sorted, x="Filename", y="FITC Stain Index", palette="magma")
    plt.xticks(rotation=90)
    plt.title("Binding Stain Index (Separation Quality)")
    plt.ylabel("Stain Index (Delta Median / 2*rSD)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "Aggregate_StainIndex_Bar.png"))
    plt.close()
    
    # 1. Bar plot of Double Positives
    plt.figure(figsize=(12, 6))
    sns.barplot(data=df_stats, x="Filename", y="Double+ %")
    plt.xticks(rotation=90)
    plt.title("Double Positive Population % by Sample")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "Aggregate_DoublePos_Bar.png"))
    plt.close()
    
    # 2. Scatter MFI
    plt.figure(figsize=(10, 6))
    # Keep X=FITC, Y=APC for MFI scatter? Or swap? User said "Binding on Y" for the main plot.
    # Usually scatter plots follow the same logic. Let's make MFI plot consistent: X=Expression(APC), Y=Binding(FITC).
    sns.scatterplot(data=df_stats, x="APC MFI", y="FITC MFI", hue="Filename", s=100)
    plt.title("Mean Fluorescence Intensity: Binding vs Expression")
    plt.ylabel("Binding MFI (FITC)")
    plt.xlabel("Expression MFI (APC)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "Aggregate_MFI_Scatter.png"))
    plt.close()
    
    # 3. Generate Text Report
    generate_report(df_stats)

def generate_report(df):
    """Generates a text report summarizing the findings."""
    
    report_path = os.path.join(OUTPUT_DIR, "Analysis_Report.md")
    
    # Identify top performers
    top_binding = df.sort_values("FITC MFI", ascending=False).iloc[0]
    top_expr = df.sort_values("APC MFI", ascending=False).iloc[0]
    
    # Group analysis
    df['Group'] = df['Filename'].apply(lambda x: x.split()[0] if ' ' in x else 'Misc')
    groups = df['Group'].unique()
    
    with open(report_path, "w") as f:
        f.write("# FCS Analysis Report\n\n")
        f.write(f"**Total Samples Processed**: {len(df)}\n")
        f.write(f"**Average Gating Efficiency**: {df['% Gated'].mean():.1f}%\n\n")
        
        f.write("## 1. Executive Summary\n")
        f.write(f"- **Highest Binding (FITC)**: {top_binding['Filename']} (MFI: {top_binding['FITC MFI']:.0f}, %+: {top_binding['FITC+ %']:.1f}%)\n")
        f.write(f"- **Highest Expression (APC)**: {top_expr['Filename']} (MFI: {top_expr['APC MFI']:.0f}, %+: {top_expr['APC+ %']:.1f}%)\n")
        f.write("\n")
        
        f.write("## 2. Group Performance\n")
        for g in groups:
            sub = df[df['Group'] == g]
            if len(sub) == 0: continue
            avg_fc = sub['FITC Fold Change'].mean()
            avg_si = sub['FITC Stain Index'].mean()
            f.write(f"### Group {g}\n")
            f.write(f"- **Samples**: {len(sub)}\n")
            f.write(f"- **Avg Binding Fold Change**: {avg_fc:.2f}x\n")
            f.write(f"- **Avg Stain Index**: {avg_si:.2f}\n")
            f.write(f"- **Key Observation**: {stats_observation(avg_fc)}\n\n")
            
        f.write("## 3. Plot Guide & Interpretation\n")
        f.write("### A. Ridgeline Plots (Aggregate_Ridgeline_*.png)\n")
        f.write("- **What it shows**: The distribution of binding signal (FITC) for every sample.\n")
        f.write("- **How to read**: Peaks further to the right indicate stronger binding. 'Wider' peaks suggest heterogeneous binding.\n")
        f.write("- **Findings**: Look for samples where the main hill is shifted clearly to the right of the red dashed line (Negative Control).\n\n")
        
        f.write("### B. Fold Change Bar Plot (Aggregate_FoldChange_Bar.png)\n")
        f.write("- **What it shows**: Ratio of Sample MFI to Negative Control MFI.\n")
        f.write("- **Threshold**: Values > 1.0 indicate signal above background. Values > 2.0 vary strong binding.\n\n")
        
        f.write("### C. Isolated Plots (Publication_Plots/)\n")
        f.write("- **What it shows**: Individual Binding vs Expression density for each sample.\n")
        f.write("- **Usage**: Use these for manuscripts. They include clean gates and percentage statistics.\n")
        
    print(f"Report generated: {report_path}")

def stats_observation(fc):
    if fc < 1.0: return "Signal suppressed or below background."
    if fc < 1.2: return "Minimal binding signal observed."
    if fc < 1.5: return "Moderate binding interaction."
    return "Strong binding signal detected."

if __name__ == "__main__":
    main()
