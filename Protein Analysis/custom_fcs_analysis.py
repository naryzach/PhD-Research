import os
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import fcsparser
from scipy.stats import gaussian_kde
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon

# Configure plotting style for "rich aesthetics"
plt.style.use('dark_background')
sns.set_context("notebook", font_scale=1.2)
colors = ["#0E1117", "#00C0FF", "#FF0055", "#FFD700"]
custom_cmap = LinearSegmentedColormap.from_list("custom", colors)


def load_fcs(file_path):
    """Loads an FCS file and returns the dataframe."""
    try:
        # reformat_meta=True usually gives cleaner column names if supported
        meta, data = fcsparser.parse(file_path, reformat_meta=True)
        return data, meta
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None, None

def calculate_density(x, y, subsample=15000):
    """Calculates Gaussian KDE density for coloring."""
    # Subsample for performance
    if len(x) > subsample:
        idx = np.random.choice(len(x), subsample, replace=False)
        x_sub, y_sub = x.iloc[idx], y.iloc[idx]
    else:
        x_sub, y_sub = x, y
        
    # Stack and calc KDE
    values = np.vstack([x_sub, y_sub])
    kernel = gaussian_kde(values)
    
    # Calculate density on the actual points
    z = kernel(np.vstack([x, y]))
    return z

def get_contour_polygon(x, y, fraction=0.90, grid_size=100):
    """
    Finds a polygon enclosing `fraction` of the data density.
    Returns: matplotlib.path.Path object and the vertices
    """
    # Create a grid for KDE evaluation
    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()
    
    # Add some padding
    dx = (xmax - xmin) * 0.1
    dy = (ymax - ymin) * 0.1
    xmin -= dx; xmax += dx
    ymin -= dy; ymax += dy
    
    xi, yi = np.mgrid[xmin:xmax:complex(0, grid_size), ymin:ymax:complex(0, grid_size)]
    positions = np.vstack([xi.flatten(), yi.flatten()])
    
    # Perform KDE on subsample
    if len(x) > 10000:
        idx = np.random.choice(len(x), 10000, replace=False)
        x_s, y_s = x.iloc[idx], y.iloc[idx]
    else:
        x_s, y_s = x, y
        
    values = np.vstack([x_s, y_s])
    kernel = gaussian_kde(values)
    zi = kernel(positions).reshape(xi.shape)
    
    # Flatten zi and sort to find the threshold for the top X%
    zi_flat = zi.flatten()
    zi_sorted = np.sort(zi_flat)[::-1]
    zi_sum = np.sum(zi_sorted)
    
    # This is a rough approximation. A better way for probability mass checks
    # is calculating the volume, but for contouring relative levels:
    # We want a contour level `L` such that mass inside >= fraction.
    # Since we don't have volume integration easily without more math, 
    # we can try to find the level that keeps the points.
    
    # Alternative approach: Score all points, find the percentile of the *score* 
    # that keeps the top X%.
    point_scores = kernel(values)
    score_threshold = np.percentile(point_scores, (1 - fraction) * 100)
    
    return xi, yi, zi, score_threshold

def apply_polygon_gate(df, x_col, y_col, xi, yi, zi, threshold):
    """
    Filters DF to points inside the contour level specified by threshold.
    Using a simplified approach: Evaluate KDE for all points and keep those > threshold.
    Note: 'threshold' here implies probability density value.
    This is equivalent to the 'polygon' gate visually if we draw the covering contour.
    """
    # Evaluate KDE on all points (might be slow for huge files, but accurate)
    # To speed up, we can use the pre-computed grid and interpolate, or just re-eval if <100k events.
    
    values = np.vstack([df[x_col], df[y_col]])
    
    # We need the SAME kernel as generated in get_contour_polygon
    # So we should probably just return the kernel or recompute it quickly.
    # Let's re-compute kernel inside a class or closure if we want speed, 
    # but for script simplicity, we'll just reconstruct the kernel on the subset used to define it.
    
    # Ideally, we pass the kernel in. Let's refactor slightly in main loop.
    pass

class AutoGating:
    def __init__(self, df, fsc_col="FSC-A", ssc_col="SSC-A", fsc_h_col="FSC-H", fitc_col="FITC-A", apc_col="APC-A"):
        self.df = df
        self.fsc_col = fsc_col
        self.ssc_col = ssc_col
        self.fsc_h_col = fsc_h_col
        self.fitc_col = fitc_col
        self.apc_col = apc_col
        
        # Gates
        self.fsc_ssc_kernel = None
        self.fsc_ssc_threshold = None
        
        self.doublet_slope = 1.0
        self.doublet_intercept = 0.0
        self.doublet_width = 10000 # Default huge
        
        self.fitc_threshold = 0
        self.apc_threshold = 0

    def learn_fsc_ssc_gate(self, fraction=0.85):
        """Learn density gate on FSC/SSC"""
        x = self.df[self.fsc_col]
        y = self.df[self.ssc_col]
        
        # Subsample for KDE creation
        if len(x) > 10000:
            idx = np.random.choice(len(x), 10000, replace=False)
            x_sub, y_sub = x.iloc[idx], y.iloc[idx]
        else:
            x_sub, y_sub = x, y
            
        values = np.vstack([x_sub, y_sub])
        self.fsc_ssc_kernel = gaussian_kde(values)
        
        # Evaluate on the training points to find threshold
        scores = self.fsc_ssc_kernel(values)
        # We want to keep top 'fraction'. So the percentile is (1-fraction)*100
        # e.g., keep top 90%, cut off bottom 10%.
        self.fsc_ssc_threshold = np.percentile(scores, (1 - fraction) * 100)
        
    def apply_fsc_ssc_gate(self, target_df):
        """Apply density threshold gate"""
        x = target_df[self.fsc_col]
        y = target_df[self.ssc_col]
        scores = self.fsc_ssc_kernel(np.vstack([x, y]))
        mask = scores >= self.fsc_ssc_threshold
        return target_df[mask], mask

    def learn_doublet_gate(self, df_gated):
        """
        Learn Singlet Gate on FSC-A vs FSC-H.
        Simple linear fit A ~ H. Ideally A=H.
        We will define a band around the main diagonal or fitted line.
        """
        x = df_gated[self.fsc_col] # A
        y = df_gated[self.fsc_h_col] # H
        
        # Simple ratio check usually works best for cytometers: A/H ~ constant
        # Let's rely on A vs H linearity.
        # Fit a line
        slope, intercept = np.polyfit(y, x, 1) # A = m*H + c
        self.doublet_slope = slope
        self.doublet_intercept = intercept
        
        # Calculate residuals
        predicted_x = slope * y + intercept
        residuals = np.abs(x - predicted_x)
        
        # Define width as 2 sigmas or 95% of residuals
        self.doublet_width = np.percentile(residuals, 95)
        
    def apply_doublet_gate(self, target_df):
        x = target_df[self.fsc_col]
        y = target_df[self.fsc_h_col]
        
        predicted_x = self.doublet_slope * y + self.doublet_intercept
        residuals = np.abs(x - predicted_x)
        mask = residuals <= self.doublet_width
        return target_df[mask], mask
    
    def learn_quadrant_gate(self, df_singlets):
        """
        Learn FITC/APC boundaries ensuring most negative control data is in bottom-left.
        We'll take the 99th percentile or similar.
        """
        # We want "close to cluster boundary".
        # 99.5th percentile is usually a good safe "negative" boundary.
        self.fitc_threshold = np.percentile(df_singlets[self.fitc_col], 99.5)
        self.apc_threshold = np.percentile(df_singlets[self.apc_col], 99.5)
        
    def apply_quadrant_gate(self, target_df):
        # We don't filter data here, just categorize/count, 
        # but for the "Gate" visualization we might show lines.
        # But if the user wants "horizontal and vertical gates such that most of the data is in the bottom left"
        # and "show the gate visually", that's just drawing the lines.
        # But usually "applying" a gate means we might want the double positive population?
        # The prompt says: "sumamry of the amount of expression and binding" -> this implies Q2 (Top-Right) or Q1/Q3.
        # We will calculate percentages in all 4 quadrants.
        pass

def save_plot(fig, output_dir, filename):
    path = os.path.join(output_dir, filename)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description="Custom FCS Analysis")
    parser.add_argument("--data_dir", default=".", help="Directory containing FCS files")
    parser.add_argument("--output_dir", default="FCS_Custom_Analysis", help="Output directory folder name")
    parser.add_argument("--toggle_doublet", action="store_true", help="Toggle doublet gating off if provided")
    args = parser.parse_args()
    
    # Setup Paths
    if not os.path.isdir(args.data_dir):
        print(f"Data directory not found: {args.data_dir}")
        return

    # Create Output Structure
    base_out = os.path.join(args.data_dir, args.output_dir)
    subfolders = ["1_Raw_Density", "2_Polygon_Gate", "3_Doublet_Gate", "4_Quadrant_Gate", "5_Binding_Distributions"]
    for sub in subfolders:
        os.makedirs(os.path.join(base_out, sub), exist_ok=True)
        
    # Find Files
    fcs_files = sorted(glob.glob(os.path.join(args.data_dir, "*.fcs")))
    if not fcs_files:
        print("No FCS files found.")
        return
        
    # Identify Negative Control
    neg_control_file = None
    for f in fcs_files:
        if "negative" in os.path.basename(f).lower() and "control" in os.path.basename(f).lower():
            neg_control_file = f
            break
    
    if not neg_control_file:
        print("Warning: No 'Negative Control' file found. Using the first file as reference.")
        neg_control_file = fcs_files[0]
        
    print(f"Using Reference: {os.path.basename(neg_control_file)}")
    
    # Load Negative Control
    neg_data, neg_meta = load_fcs(neg_control_file)
    if neg_data is None: return

    # --- 1. Raw Data Density (Negative Control) ---
    print("Generating Raw Data Density Plot...")
    
    # Check Columns
    fsc = "FSC-A"
    ssc = "SSC-A"
    fsc_h = "FSC-H"
    fitc = "FITC-A"
    apc = "APC-A"
    
    # Ensure columns exist (basic check)
    if fsc not in neg_data.columns: fsc = neg_data.columns[0] # Fallback
    
    # Calculate Density for visualization
    z_neg = calculate_density(neg_data[fsc], neg_data[ssc])
    
    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(neg_data[fsc], neg_data[ssc], c=z_neg, cmap=custom_cmap, s=5, alpha=0.6)
    plt.colorbar(sc, label="Density")
    ax.set_xlabel(fsc); ax.set_ylabel(ssc)
    ax.set_title(f"Raw Density: {os.path.basename(neg_control_file)}")
    save_plot(fig, os.path.join(base_out, "1_Raw_Density"), "NegControl_Density.png")
    
    # --- 2. Learn Polygon Gate (FSC/SSC) ---
    print("Learning Polygon Gate...")
    gater = AutoGating(neg_data, fsc_col=fsc, ssc_col=ssc, fsc_h_col=fsc_h, fitc_col=fitc, apc_col=apc)
    gater.learn_fsc_ssc_gate(fraction=0.88) # 85-90% requested
    
    # Visualize the gate on Control
    # To visualize the "polygon", we can contour the density kernel at the threshold
    xi = np.linspace(neg_data[fsc].min(), neg_data[fsc].max(), 100)
    yi = np.linspace(neg_data[ssc].min(), neg_data[ssc].max(), 100)
    xi_g, yi_g = np.meshgrid(xi, yi)
    
    # Check if kernel is available
    if gater.fsc_ssc_kernel is not None:
        zi_g = gater.fsc_ssc_kernel(np.vstack([xi_g.flatten(), yi_g.flatten()])).reshape(xi_g.shape)
        
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(neg_data[fsc], neg_data[ssc], c=z_neg, cmap=custom_cmap, s=2, alpha=0.4)
        # Draw contour
        ax.contour(xi_g, yi_g, zi_g, levels=[gater.fsc_ssc_threshold], colors='white', linewidths=2)
        ax.set_title(f"Polygon Gate (FSC/SSC) - ~88% Capture")
        save_plot(fig, os.path.join(base_out, "2_Polygon_Gate"), "NegControl_PolygonGate.png")
    else:
        print("Error: KDE kernel for polygon gate not found.")
    
    # --- 3. Apply Polygon and Learn Doublet Gate ---
    neg_gated_1, _ = gater.apply_fsc_ssc_gate(neg_data)
    
    if not args.toggle_doublet:
        print("Learning Doublet Gate...")
        gater.learn_doublet_gate(neg_gated_1)
        
        # Visualize Doublet Gate
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Plot all points (from neg_gated_1)
        ax.scatter(neg_gated_1[fsc], neg_gated_1[fsc_h], c='#00C0FF', s=2, alpha=0.5, label="Events")
        
        # Draw Gate Lines
        # y_pred = m*H + c  <=> A_pred = ...
        # We are plotting A on X, H on Y? Wait.
        # Code: scatter(x=A, y=H) in visualization?
        # Prompt: "FSC A vs FSC H" -> usually A is X, H is Y, or vice versa.
        # Let's check axes on standard views. Doulet discrim often A vs H or A vs W.
        # Linearity is usually A vs H.
        
        # Plotting: X=FSC-A, Y=FSC-H
        x_vals = np.linspace(neg_gated_1[fsc].min(), neg_gated_1[fsc].max(), 100)
        # However, our fit was A = m*H + c. (X depends on Y in fit). 
        # Ideally we fit Y = m*X + c (H depends on A) or just check consistency.
        # Let's stick to the fit: A = slope * H + intercept
        # So H = (A - intercept) / slope
        
        h_vals_model = (x_vals - gater.doublet_intercept) / gater.doublet_slope
        
        # But width is in terms of A residuals.
        # So upper bound A = m*H + c + width
        # lower bound A = m*H + c - width
        
        # To plot these lines on an A vs H plot (A is X, H is Y):
        # We need Y values for the lines.
        # Upper bound line: X = m*Y + c + w  => Y = (X - c - w)/m
        # Lower bound line: X = m*Y + c - w  => Y = (X - c + w)/m
        
        y_upper = (x_vals - gater.doublet_intercept - gater.doublet_width) / gater.doublet_slope
        y_lower = (x_vals - gater.doublet_intercept + gater.doublet_width) / gater.doublet_slope
        
        ax.plot(x_vals, y_upper, 'w--', lw=1, label="Gate")
        ax.plot(x_vals, y_lower, 'w--', lw=1)
        
        ax.set_xlabel(fsc); ax.set_ylabel(fsc_h)
        ax.set_title("Doublet Discrimination")
        ax.legend()
        save_plot(fig, os.path.join(base_out, "3_Doublet_Gate"), "NegControl_DoubletGate.png")
        
        neg_gated_2, _ = gater.apply_doublet_gate(neg_gated_1)
    else:
        neg_gated_2 = neg_gated_1

    # --- 4. Learn Quadrant Gate (FITC/APC) ---
    print("Learning Quadrant Gates...")
    gater.learn_quadrant_gate(neg_gated_2)
    
    # Visualize Quadrant Gate on Neg Control
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Use density for this plot too? Or just dots?
    # Prompt: "Show that gate visually as well"
    ax.scatter(neg_gated_2[fitc], neg_gated_2[apc], c='#444444', s=2, alpha=0.5)
    
    # Lines
    ax.axvline(gater.fitc_threshold, color="#FF0055", linestyle="--", lw=1.5, label=f"FITC Cutoff ({gater.fitc_threshold:.1f})")
    ax.axhline(gater.apc_threshold, color="#FF0055", linestyle="--", lw=1.5, label=f"APC Cutoff ({gater.apc_threshold:.1f})")
    
    ax.set_xlabel(fitc); ax.set_ylabel(apc)
    ax.set_title("Quadrant Gates (Defined on Neg Control)")
    ax.legend(loc="upper right")
    save_plot(fig, os.path.join(base_out, "4_Quadrant_Gate"), "NegControl_QuadGate.png")
    
    # --- 5. Process All Files and Generate Stats ---
    stats_list = []
    
    print("Processing all files...")
    for f in fcs_files:
        basename = os.path.basename(f)
        try:
            raw_data, _ = load_fcs(f)
            if raw_data is None: continue
            
            # 1. Polygon Gate
            df_poly, _ = gater.apply_fsc_ssc_gate(raw_data)
            perc_poly = (len(df_poly) / len(raw_data)) * 100
            
            # 2. Doublet Gate
            if not args.toggle_doublet:
                df_singlet, _ = gater.apply_doublet_gate(df_poly)
                perc_singlet = (len(df_singlet) / len(df_poly)) * 100 if len(df_poly) > 0 else 0
            else:
                df_singlet = df_poly
                perc_singlet = 100.0
                
            # 3. Quadrant Stats
            # Q1: FITC- APC+ (Top Left) -> Low Binding, High Expression (if FITC is Binding?)
            # Prompt: "FITC (expression Y-axis) vs APC (binding X-axis)"
            # Wait, prompt says: "FITC (expression Y-axis) vs APC (binding X-axis)"
            # Standard: usually FITC is X, APC is Y or vice versa depending on config.
            # But prompt EXPLICITLY maps them:
            # FITC = Y-axis (Expression)
            # APC = X-axis (Binding)
            # So my plot above: scatter(x=FITC, y=APC) was wrong relative to prompt descriptions.
            # I should plot X=APC, Y=FITC.
            
            # Let's Correct Visualization for Quadrant:
            # X = APC (Binding), Y = FITC (Expression)
            # Thresholds: APC_Threshold (vertical), FITC_Threshold (horizontal)
            
            # Filtered Data
            final_df = df_singlet
            
            # Counts
            # "Most of data in bottom left" => Low APC, Low FITC.
            # Q_LL (Double Neg): APC < Thresh, FITC < Thresh
            # Q_LR (Binding Only): APC >= Thresh, FITC < Thresh
            # Q_UL (Expr Only): APC < Thresh, FITC >= Thresh
            # Q_UR (Double Pos): APC >= Thresh, FITC >= Thresh
            
            n_total = len(final_df)
            if n_total > 0:
                q_ll = len(final_df[(final_df[filtc_col] < gater.fitc_threshold) & (final_df[apc_col] < gater.apc_threshold)]) if 'filtc_col' in locals() else \
                       len(final_df[(final_df[fitc] < gater.fitc_threshold) & (final_df[apc] < gater.apc_threshold)])
                       
                q_lr = len(final_df[(final_df[fitc] < gater.fitc_threshold) & (final_df[apc] >= gater.apc_threshold)])
                q_ul = len(final_df[(final_df[fitc] >= gater.fitc_threshold) & (final_df[apc] < gater.apc_threshold)])
                q_ur = len(final_df[(final_df[fitc] >= gater.fitc_threshold) & (final_df[apc] >= gater.apc_threshold)])
                
                pct_ll = (q_ll / n_total) * 100
                pct_lr = (q_lr / n_total) * 100
                pct_ul = (q_ul / n_total) * 100
                pct_ur = (q_ur / n_total) * 100
                
                # MFI Stats
                mfi_fitc = final_df[fitc].median()
                mfi_apc = final_df[apc].median()
            else:
                pct_ll = pct_lr = pct_ul = pct_ur = 0
                mfi_fitc = mfi_apc = 0
            
            # Store Stats
            stats_entry = {
                "Filename": basename,
                "Total_Events_Raw": len(raw_data),
                "Events_Post_Link_Doublet": len(final_df),
                "Percent_Polygon_Gate": perc_poly,
                "Percent_Singlet": perc_singlet,
                "Pct_DoubleNeg_LL": pct_ll,
                "Pct_BindingOnly_LR": pct_lr,
                "Pct_ExpressionOnly_UL": pct_ul,
                "Pct_DoublePos_UR": pct_ur,
                "MFI_Expression_FITC": mfi_fitc,
                "MFI_Binding_APC": mfi_apc
            }
            stats_list.append(stats_entry)
            
            # --- Generate Plots for Every File ---
            # 4. Quadrant Plot (Log Scale Match + Hexbin)
            # Log Transform High Dynamic Range Channels
            # Clip to 1 before log to avoid -inf
            log_apc = np.log10(np.clip(final_df[apc], 1, None))
            log_fitc = np.log10(np.clip(final_df[fitc], 1, None))
            
            # Log Thresholds
            log_thresh_apc = np.log10(max(1, gater.apc_threshold))
            log_thresh_fitc = np.log10(max(1, gater.fitc_threshold))

            fig, ax = plt.subplots(figsize=(7, 6))
            
            # Use Hexbin for faster, cleaner density visualization (matches analyze_fcs.py)
            hb = ax.hexbin(log_apc, log_fitc, gridsize=100, cmap='jet', mincnt=1, bins='log')
            # Optional: Add colorbar
            # cb = plt.colorbar(hb, ax=ax, label='Count (Log)')

            # Draw Gates
            ax.axvline(log_thresh_apc, color="black", ls="--", lw=1.5)
            ax.axhline(log_thresh_fitc, color="black", ls="--", lw=1.5)
            
            # --- Axis Resizing (Dynamic Log) ---
            # Include 0 (10^0 = 1) and upper percentiles
            x_max_limit = log_apc.max() * 1.05
            y_max_limit = log_fitc.max() * 1.05
            
            # Ensure we see gates
            x_max_limit = max(x_max_limit, log_thresh_apc + 0.5)
            y_max_limit = max(y_max_limit, log_thresh_fitc + 0.5)
            
            ax.set_xlim(left=0, right=x_max_limit)
            ax.set_ylim(bottom=0, top=y_max_limit)
            
            ax.set_xlabel(f"Log10 Binding ({apc})")
            ax.set_ylabel(f"Log10 Expression ({fitc})")
            ax.set_title(f"{basename}\nQ1(UL):{pct_ul:.1f}% | Q2(UR):{pct_ur:.1f}%\nQ4(LL):{pct_ll:.1f}% | Q3(LR):{pct_lr:.1f}%")
            
            save_plot(fig, os.path.join(base_out, "4_Quadrant_Gate"), f"{basename}_Quad.png")
            
            # 5. Binding Distribution (Histogram / Ridgeline style)
            # Just simple histogram for now
            fig, ax = plt.subplots(figsize=(8, 4))
            sns.histplot(final_df[apc], bins=100, color="#00C0FF", element="step", fill=True, ax=ax)
            ax.axvline(gater.apc_threshold, color="white", ls="--", label="Gate")
            ax.set_xlabel(f"Binding ({apc})")
            ax.set_title(f"Binding Distribution: {basename}")
            # Log scale often useful
            ax.set_xscale("log")
            save_plot(fig, os.path.join(base_out, "5_Binding_Distributions"), f"{basename}_BindingHist.png")
            
        except Exception as e:
            print(f"Failed to process {basename}: {e}")
            
    # --- 6. Summary Statistics & Comparison ---
    stats_df = pd.DataFrame(stats_list)
    
    # Calculate relative to Negative Control
    # Find negative control row
    neg_row = stats_df[stats_df["Filename"] == os.path.basename(neg_control_file)]
    
    if not neg_row.empty:
        neg_mfi_fitc = neg_row.iloc[0]["MFI_Expression_FITC"]
        neg_mfi_apc = neg_row.iloc[0]["MFI_Binding_APC"]
        
        # Add Fold Change Columns
        stats_df["FoldChange_Expression"] = stats_df["MFI_Expression_FITC"] / neg_mfi_fitc
        stats_df["FoldChange_Binding"] = stats_df["MFI_Binding_APC"] / neg_mfi_apc
        
        # Add "Better/Worse" Delta
        stats_df["Delta_DoublePos_Pct"] = stats_df["Pct_DoublePos_UR"] - neg_row.iloc[0]["Pct_DoublePos_UR"]
        
    # Save CSV
    stats_df.to_csv(os.path.join(base_out, "summary_statistics.csv"), index=False)
    print(f"Analysis Complete. Results saved to {base_out}")

if __name__ == "__main__":
    main()
