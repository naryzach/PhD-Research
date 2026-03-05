import os
import sys
import glob
import argparse
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

required_packages = ['fcsparser', 'seaborn', 'scipy', 'matplotlib', 'pandas', 'numpy', 'scikit-image']
for pkg in required_packages:
    try:
        if pkg == 'scikit-image':
            __import__('skimage')
        else:
            __import__(pkg)
    except ImportError:
        print(f"Installing {pkg}...")
        install(pkg)

import fcsparser
from scipy.stats import gaussian_kde
from matplotlib.path import Path
from skimage import measure

# ---- CONFIGURATION ----
NEG_CONTROL_PATTERNS = ["NC", "Negative Control"] # Starts with any of these
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

from scipy.spatial import ConvexHull

def simplify_polygon_vw(points, num_points=5):
    """
    Simplifies a convex polygon to exactly `num_points` vertices
    using the Visvalingam-Whyatt area-based minimization algorithm.
    This guarantees the minimal possible bounded area while remaining convex 
    and strictly inside the original hull!
    """
    pts = list(points)
    # Ensure no duplicates at the end
    if np.allclose(pts[0], pts[-1]):
        pts.pop()
        
    while len(pts) > num_points:
        min_area = float('inf')
        min_idx = -1
        for i in range(len(pts)):
            p_prev = pts[i - 1]
            p_curr = pts[i]
            p_next = pts[(i + 1) % len(pts)]
            # Triangle area
            area = 0.5 * abs(p_prev[0]*(p_curr[1] - p_next[1]) + 
                             p_curr[0]*(p_next[1] - p_prev[1]) + 
                             p_next[0]*(p_prev[1] - p_curr[1]))
            if area < min_area:
                min_area = area
                min_idx = i
        pts.pop(min_idx)
    return np.array(pts)

def learn_pentagon_gate(df, fsc_col, ssc_col, fraction=0.90):
    """
    Learns a 5-point polygon (pentagon) gate on FSC vs SSC.
    Actively rejects points near the origin (debris) before contouring 
    to ensure the resulting polygon is tightly wrapped around the core events.
    Captures approximately `fraction` (e.g., 90%) of the events,
    Returns the matplotlib Path object representing the pentagon.
    """
    # 1. Reject obvious origin debris BEFORE learning density
    # This ensures the KDE isn't skewed down towards (0,0)
    origin_filter = (df[fsc_col] > 5000) & (df[ssc_col] > 1000)
    core_df = df[origin_filter]
    
    if len(core_df) < 100:
        print("Warning: Too few events after debris filtering. Using all data.")
        core_df = df
        
    x = core_df[fsc_col].values
    y = core_df[ssc_col].values
    
    # 2. Filter upper extreme outliers so KDE is focused tightly on main clusters
    # This prevents the boundary from getting stretched to 10^7
    fsc_99 = np.percentile(x, 99)
    ssc_99 = np.percentile(y, 99)
    core_idx = (x < fsc_99) & (y < ssc_99)
    
    x_core = x[core_idx]
    y_core = y[core_idx]
    
    # Subsample for KDE performance
    n_fit = min(10000, len(x_core))
    idx_fit = np.random.choice(len(x_core), n_fit, replace=False)
    fit_points = np.vstack([x_core[idx_fit], y_core[idx_fit]])
    
    kernel = gaussian_kde(fit_points)
    
    # Score a representative sample to quickly find the threshold
    # Using the real data points themselves means the gate perfectly wraps the data
    all_x = df[fsc_col].values
    all_y = df[ssc_col].values
    
    n_score = min(20000, len(all_x))
    idx_score = np.random.choice(len(all_x), n_score, replace=False)
    score_points = np.vstack([all_x[idx_score], all_y[idx_score]])
    
    scores_sub = kernel(score_points)
    
    target_captured = fraction
    low_pct = 0.0
    high_pct = 100.0
    
    best_path = None
    best_captured = 0
    
    # Binary search for the right contour threshold
    for _ in range(25): # Binary search depth 
        mid_pct = (low_pct + high_pct) / 2
        
        # We want the threshold such that we drop the lowest (mid_pct) scores.
        # e.g., if mid_pct = 10, we keep top 90%
        thresh = np.percentile(scores_sub, mid_pct)
        
        mask_sub = scores_sub >= thresh
        if np.sum(mask_sub) < 10:
            high_pct = mid_pct
            continue
            
        in_pts = score_points.T[mask_sub]
        try:
            hull = ConvexHull(in_pts)
            hull_pts = in_pts[hull.vertices]
        except Exception:
            # Collinear or too few points
            high_pct = mid_pct
            continue
            
        # Simplify to exactly 5 points using VW area minimization!
        pent_pts = simplify_polygon_vw(hull_pts, 5)
        pent_pts = np.vstack((pent_pts, pent_pts[0])) # close it
        
        path = Path(pent_pts)
        
        # Test how many of the exact data points fall in this path
        captured = path.contains_points(score_points.T).mean()
        
        if abs(captured - target_captured) < abs(best_captured - target_captured):
            best_captured = captured
            best_path = path

        if captured > target_captured:
            # We captured too much. To shrink the area, we need a TIGHTER cluster,
            # which means dropping MORE points, which means a HIGHER percentile threshold.
            low_pct = mid_pct
        else:
            high_pct = mid_pct
            
        if abs(captured - target_captured) < 0.005:
            break
            
    if best_path is None:
        print("Warning: Could not form a proper polygon. Using a wide bounding box.")
        xmin, xmax = x.min(), x.max()
        ymin, ymax = y.min(), y.max()
        poly_verts = [(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)]
        best_path = Path(poly_verts)

    exact_captured = best_path.contains_points(np.vstack([all_x, all_y]).T).mean()
    print(f"Pentagon generated: captures ~{exact_captured*100:.1f}% events.")
    
    return best_path


def apply_polygon_gate(df, path, fsc_col, ssc_col, plot_ax=None, title="Polygon Gate"):
    """
    Applies a matplotlib Path polygon gate to the dataframe.
    """
    if len(df) == 0:
        return df, []
        
    points = np.vstack([df[fsc_col], df[ssc_col]]).T
    mask = path.contains_points(points)
    
    if plot_ax:
        # Density plot 
        plot_ax.hexbin(df[fsc_col], df[ssc_col], gridsize=100, cmap='inferno', mincnt=1, bins='log')
        plot_ax.set_xlabel(fsc_col)
        plot_ax.set_ylabel(ssc_col)
        plot_ax.set_title(f"{title} ({mask.sum()}/{len(df)} events)")
        
        # Overlay Polygon
        import matplotlib.patches as patches
        patch = patches.PathPatch(path, facecolor='none', edgecolor='cyan', lw=2, linestyle='--')
        plot_ax.add_patch(patch)
        
    return df[mask]


def transform_logicle(data, channels):
    """Simple log10 transform for plotting, handling negatives safely (clipping)."""
    df_trans = data.copy()
    for ch in channels:
        df_trans[ch] = np.log10(np.clip(df_trans[ch], 1, None))
    return df_trans

def main():
    parser = argparse.ArgumentParser(description="Analyze FCS files and generate gating plots and statistics.")
    parser.add_argument("-i", "--input", required=True, help="Path to input directory containing .fcs files.")
    parser.add_argument("-o", "--output", default=None, help="Path to output directory. Defaults to '<input>_Analysis'.")
    args = parser.parse_args()

    input_dir = os.path.abspath(args.input)
    if args.output:
        output_dir = os.path.abspath(args.output)
    else:
        # Generate output directory automatically
        base = os.path.basename(os.path.normpath(input_dir))
        output_dir = os.path.join(os.path.dirname(input_dir), f"{base}_Analysis")

    # 1. Setup
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # Updated Folder Structure
    individual_dir = os.path.join(output_dir, "Individual_Plots")
    if not os.path.exists(individual_dir):
        os.makedirs(individual_dir)

    pub_dir = os.path.join(output_dir, "Publication_Ready")
    if not os.path.exists(pub_dir):
        os.makedirs(pub_dir)
        
    agg_dir = os.path.join(output_dir, "Aggregate_Plots")
    if not os.path.exists(agg_dir):
        os.makedirs(agg_dir)

    comp_bind_dir = os.path.join(output_dir, "Comparisons_vs_PosCtrl", "Binding")
    if not os.path.exists(comp_bind_dir):
        os.makedirs(comp_bind_dir)
        
    comp_expr_dir = os.path.join(output_dir, "Comparisons_vs_PosCtrl", "Expression")
    if not os.path.exists(comp_expr_dir):
        os.makedirs(comp_expr_dir)
        
    pos_dist_dir = os.path.join(output_dir, "Positive_Distributions")
    if not os.path.exists(pos_dist_dir):
        os.makedirs(pos_dist_dir)
        
    pos_ratio_dir = os.path.join(output_dir, "Positive_Ratios")
    if not os.path.exists(pos_ratio_dir):
        os.makedirs(pos_ratio_dir)
        
    # Find files
    search_path = input_dir
    
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
    neg_files = [f for f in all_files if any(os.path.basename(f).startswith(p) for p in NEG_CONTROL_PATTERNS)]
    print(f"Found {len(neg_files)} negative control files.")
    
    # Learn Pentagon Gate from the First Negative Control
    pentagon_path = None
    if neg_files:
        print("Learning pentagon gate from Negative Control...")
        _, d_neg_learn = load_fcs(neg_files[0])
        if d_neg_learn is not None:
            # Drop pure zeros or negatives for KDE stability if needed, 
            # but usually the origin points are handled by KDE
            pentagon_path = learn_pentagon_gate(d_neg_learn, CH_FSC_A, CH_SSC_A, fraction=0.90)
            
    if pentagon_path is None:
        print("Warning: Could not learn pentagon gate from NC. Cannot proceed accurately.")
        return

    # Process all Negative Controls with learned gate
    neg_dfs = []
    for f in neg_files:
        _, d = load_fcs(f)
        if d is not None:
            d_gated = apply_polygon_gate(d, pentagon_path, CH_FSC_A, CH_SSC_A)
            d_gated = d_gated[(d_gated[CH_EXPR] > 0) & (d_gated[CH_BIND] > 0)]
            neg_dfs.append(d_gated)
            
    if neg_dfs:
        neg_concat = pd.concat(neg_dfs)
        # Learn Quadrant Gate Thresholds (99.9% NC in lower left)
        # i.e., threshold is 99.9th percentile of NC
        thresh_expr = np.percentile(neg_concat[CH_EXPR], 99.9)
        thresh_bind = np.percentile(neg_concat[CH_BIND], 99.9)
        
        neg_mfi_expr = neg_concat[CH_EXPR].median()
        neg_mfi_bind = neg_concat[CH_BIND].median()
        neg_rsd_bind = stats.median_abs_deviation(neg_concat[CH_BIND], scale='normal')
        neg_rsd_expr = stats.median_abs_deviation(neg_concat[CH_EXPR], scale='normal')

        # --- GATING VISUALIZATION (Representative Neg Control) ---
        print("Generating Gating Strategy Plot overlay on NC...")
        _, d_demo = load_fcs(neg_files[0])
        if d_demo is not None:
            fig, ax = plt.subplots(1, 1, figsize=(6, 5))
            
            # Scatter/Density Gate Overlay
            d_gated_demo = apply_polygon_gate(d_demo, pentagon_path, CH_FSC_A, CH_SSC_A, plot_ax=ax, title="Pentagon Gate on NC")
            d_gated_demo = d_gated_demo[(d_gated_demo[CH_EXPR] > 0) & (d_gated_demo[CH_BIND] > 0)]
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "Gating_Strategy_NegCtrl.png"))
            plt.close()
            
            # Pure Density Gate Overlay (KDE)
            fig_kde, ax_kde = plt.subplots(1, 1, figsize=(6, 5))
            sns.kdeplot(x=d_demo[CH_FSC_A], y=d_demo[CH_SSC_A], fill=True, cmap="mako", ax=ax_kde, log_scale=(True, True))
            import matplotlib.patches as patches
            patch_kde = patches.PathPatch(pentagon_path, facecolor='none', edgecolor='cyan', lw=2, linestyle='--')
            ax_kde.add_patch(patch_kde)
            ax_kde.set_xlabel(CH_FSC_A)
            ax_kde.set_ylabel(CH_SSC_A)
            ax_kde.set_title(f"Pentagon Gate on NC (Density)")
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "Gating_Strategy_NegCtrl_Density.png"))
            plt.close()
            
            # Scatter/Density for Quad gate on NC
            fig2, ax2 = plt.subplots(1, 1, figsize=(6, 6))
            t_bind_nc = np.log10(np.clip(d_gated_demo[CH_BIND], 1, None))
            t_expr_nc = np.log10(np.clip(d_gated_demo[CH_EXPR], 1, None))
            
            ax2.hexbin(t_bind_nc, t_expr_nc, gridsize=100, cmap='jet', mincnt=1, bins='log')
            ax2.axvline(np.log10(max(thresh_bind, 1)), color='r', linestyle='--', label='99.9% Bind Thresh')
            ax2.axhline(np.log10(max(thresh_expr, 1)), color='r', linestyle='--', label='99.9% Expr Thresh')
            ax2.set_xlabel("Log10 APC-A (Binding)")
            ax2.set_ylabel("Log10 FITC-A (Expression)")
            ax2.set_title("99.9% Thresholds on Negative Control")
            ax2.legend(loc='lower left')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "Gating_Strategy_NegCtrl_Quad.png"))
            plt.close()
            
            # 1D Density for Expression on NC
            fig_ex, ax_ex = plt.subplots(1, 1, figsize=(6, 4))
            sns.kdeplot(t_expr_nc, fill=True, ax=ax_ex, color='purple')
            ax_ex.axvline(np.log10(max(thresh_expr, 1)), color='r', linestyle='--', label=f'99.9% Thresh: {thresh_expr:.0f}')
            ax_ex.set_xlabel("Log10 FITC-A (Expression)")
            ax_ex.set_title("Negative Control Expression Density")
            ax_ex.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "Gating_Strategy_NegCtrl_Expr_Density.png"))
            plt.close()

            # 1D Density for Binding on NC
            fig_bi, ax_bi = plt.subplots(1, 1, figsize=(6, 4))
            sns.kdeplot(t_bind_nc, fill=True, ax=ax_bi, color='g')
            ax_bi.axvline(np.log10(max(thresh_bind, 1)), color='r', linestyle='--', label=f'99.9% Thresh: {thresh_bind:.0f}')
            ax_bi.set_xlabel("Log10 APC-A (Binding)")
            ax_bi.set_title("Negative Control Binding Density")
            ax_bi.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "Gating_Strategy_NegCtrl_Bind_Density.png"))
            plt.close()

    else:
        print("Error: No negative controls found matching pattern. Required for dynamic gating.")
        return

    # Positive Controls
    pos_files = [f for f in all_files if any(os.path.basename(f).startswith(p) for p in POS_CONTROL_PATTERNS)]
    print(f"Found {len(pos_files)} positive control files (Patterns: {POS_CONTROL_PATTERNS}).")
    
    pos_dfs = []
    for f in pos_files:
        _, d = load_fcs(f)
        if d is not None:
            d_gated = apply_polygon_gate(d, pentagon_path, CH_FSC_A, CH_SSC_A)
            d_gated = d_gated[(d_gated[CH_EXPR] > 0) & (d_gated[CH_BIND] > 0)]
            pos_dfs.append(d_gated)
            
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
    
    # Pre-calculate log data for controls for plotting
    neg_log_bind = np.log10(np.clip(neg_concat[CH_BIND].sample(n=min(len(neg_concat), 5000)), 1, None)) if neg_concat is not None else None
    pos_log_bind = np.log10(np.clip(pos_concat[CH_BIND].sample(n=min(len(pos_concat), 5000)), 1, None)) if pos_concat is not None else None
    
    neg_log_expr = np.log10(np.clip(neg_concat[CH_EXPR].sample(n=min(len(neg_concat), 5000)), 1, None)) if neg_concat is not None else None
    pos_log_expr = np.log10(np.clip(pos_concat[CH_EXPR].sample(n=min(len(pos_concat), 5000)), 1, None)) if pos_concat is not None else None

    for f in all_files:
        fname = os.path.basename(f)
        clean_name = os.path.splitext(fname)[0]
        
        print(f"Processing {clean_name}...")
        
        meta, df = load_fcs(f)
        if df is None: continue
        
        # Gating
        df_sing = apply_polygon_gate(df, pentagon_path, CH_FSC_A, CH_SSC_A)
        df_sing = df_sing[(df_sing[CH_EXPR] > 0) & (df_sing[CH_BIND] > 0)]
        
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
        
        # MFI Stats
        mfi_expr = df_sing[CH_EXPR].median() if len(df_sing) > 0 else 0
        mfi_bind = df_sing[CH_BIND].median() if len(df_sing) > 0 else 0
        
        # Formatting Ridge Data (Binding & Expression)
        if len(df_sing) > 0:
            sub = df_sing.sample(n=min(len(df_sing), 2000))
            t_vals_bind = np.log10(np.clip(sub[CH_BIND], 1, None))
            t_vals_expr = np.log10(np.clip(sub[CH_EXPR], 1, None))
            
            for i in range(len(sub)):
                ridge_data.append({
                    "Sample": clean_name, 
                    "LogBinding": t_vals_bind.iloc[i],
                    "LogExpression": t_vals_expr.iloc[i]
                })

        # Advanced Stats
        fc_expr = mfi_expr / max(neg_mfi_expr, 1)
        fc_bind = mfi_bind / max(neg_mfi_bind, 1)
        
        # % of Pos Ctrl (Binding)
        pct_of_pos_mfi = (mfi_bind / pos_mfi_bind) * 100 if pos_mfi_bind > 0 else 0
        
        # Fold Change vs Positive Control (Ratio)
        fc_vs_pos = mfi_bind / pos_mfi_bind if pos_mfi_bind > 0 else 0
        
        # Binding Efficiency (Binding normalized to Expression = Binding / Expression)
        binding_eff = mfi_bind / mfi_expr if mfi_expr > 1 else 0
        
        si_bind = (mfi_bind - neg_mfi_bind) / (2 * neg_rsd_bind) if neg_rsd_bind > 0 else 0
        si_expr = (mfi_expr - neg_mfi_expr) / (2 * neg_rsd_expr) if neg_rsd_expr > 0 else 0
        
        # Positive Population Analysis
        df_pos_expr = df_sing[expr_pos]
        df_pos_bind = df_sing[bind_pos]
        
        pos_mean_expr = df_pos_expr[CH_EXPR].mean() if len(df_pos_expr) > 0 else 0
        pos_med_expr = df_pos_expr[CH_EXPR].median() if len(df_pos_expr) > 0 else 0
        
        pos_mean_bind = df_pos_bind[CH_BIND].mean() if len(df_pos_bind) > 0 else 0
        pos_med_bind = df_pos_bind[CH_BIND].median() if len(df_pos_bind) > 0 else 0
        
        # Ratios (Bind / Expr) for Positive events
        pos_mean_ratio = pos_mean_bind / pos_mean_expr if pos_mean_expr > 0 else 0
        pos_med_ratio = pos_med_bind / pos_med_expr if pos_med_expr > 0 else 0

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
            "Bind FC vs Pos Ctrl": fc_vs_pos,
            "Bind Stain Index": si_bind,
            "Expr Stain Index": si_expr,
            "% of Pos Ctrl": pct_of_pos_mfi,
            "Binding Efficiency": binding_eff,
            "Pos Mean Bind": pos_mean_bind,
            "Pos Med Bind": pos_med_bind,
            "Pos Mean Expr": pos_mean_expr,
            "Pos Med Expr": pos_med_expr,
            "Pos Mean Ratio": pos_mean_ratio,
            "Pos Med Ratio": pos_med_ratio
        })
        
        # --- POSITIVE DISTRIBUTIONS PLOT ---
        if len(df_pos_expr) > 0 or len(df_pos_bind) > 0:
            fig, axes = plt.subplots(1, 2, figsize=(10, 4))
            
            if len(df_pos_expr) > 0:
                t_pos_expr = np.log10(np.clip(df_pos_expr[CH_EXPR], 1, None))
                sns.histplot(t_pos_expr, ax=axes[0], color='purple', stat='count', element='step', fill=True)
                
                # We plot means and medians 
                m_mean_ex = np.log10(max(1, pos_mean_expr))
                m_med_ex = np.log10(max(1, pos_med_expr))
                
                axes[0].axvline(m_mean_ex, color='k', linestyle='--', label=f'Mean: {pos_mean_expr:.0f}')
                axes[0].axvline(m_med_ex, color='b', linestyle=':', label=f'Median: {pos_med_expr:.0f}')
                axes[0].set_title("Positive Expression Dist")
                axes[0].set_xlabel("Log10 FITC-A (Expression)")
                axes[0].legend()
            
            if len(df_pos_bind) > 0:
                t_pos_bind = np.log10(np.clip(df_pos_bind[CH_BIND], 1, None))
                sns.histplot(t_pos_bind, ax=axes[1], color='g', stat='count', element='step', fill=True)
                
                m_mean_bi = np.log10(max(1, pos_mean_bind))
                m_med_bi = np.log10(max(1, pos_med_bind))
                
                axes[1].axvline(m_mean_bi, color='k', linestyle='--', label=f'Mean: {pos_mean_bind:.0f}')
                axes[1].axvline(m_med_bi, color='b', linestyle=':', label=f'Median: {pos_med_bind:.0f}')
                axes[1].set_title("Positive Binding Dist")
                axes[1].set_xlabel("Log10 APC-A (Binding)")
                axes[1].legend()
                
            plt.suptitle(f"{clean_name} - Positive Events Distributions")
            plt.tight_layout()
            plt.savefig(os.path.join(pos_dist_dir, f"{clean_name}_PosDist.png"))
            plt.close()

        # --- PLOTTING ---
        
        # 1. Main 4-panel Plot
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        
        # Panel 1: Original Scatter Overlay with Pentagon Gate
        apply_polygon_gate(df, pentagon_path, CH_FSC_A, CH_SSC_A, plot_ax=axes[0,0], title="Pentagon Gating (from NC)")
        
        # Log Transformations
        t_expr = np.log10(np.clip(df_sing[CH_EXPR], 1, None))
        t_bind = np.log10(np.clip(df_sing[CH_BIND], 1, None))
        log_thresh_bind = np.log10(max(1, thresh_bind))
        log_thresh_expr = np.log10(max(1, thresh_expr))
        
        # Panel 2: Expression Distribution (FITC Histogram)
        sns.kdeplot(t_expr, fill=True, ax=axes[0,1], color='purple', label='Sample')
        if neg_log_expr is not None:
            sns.kdeplot(neg_log_expr, fill=True, ax=axes[0,1], color='gray', alpha=0.3, label='Neg Ctrl')
        axes[0,1].axvline(log_thresh_expr, color='r', linestyle='--')
        pct_expr_left = 100 - pct_expr
        axes[0,1].text(0.05, 0.95, f"Unexpressed: {pct_expr_left:.1f}%\nExpressed: {pct_expr:.1f}%", 
                       transform=axes[0,1].transAxes, fontsize=10, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        axes[0,1].set_xlabel("Log10 FITC (Expression)")
        axes[0,1].set_title(f"Expression Dist (FC: {fc_expr:.1f}x)")
        axes[0,1].legend()

        # Panel 3: Binding vs Expression (Hexbin)
        axes[1,0].hexbin(t_bind, t_expr, gridsize=100, cmap='jet', mincnt=1, bins='log')
        axes[1,0].axvline(log_thresh_bind, color='k', linestyle='--')
        axes[1,0].axhline(log_thresh_expr, color='k', linestyle='--')
        axes[1,0].set_xlabel("Log10 APC (Binding)")
        axes[1,0].set_ylabel("Log10 FITC (Expression)")
        pct_ll = 100 - pct_expr - pct_bind + pct_double # roughly, or explicitly calculated:
        expr_pos_only = (df_sing[CH_EXPR] >= thresh_expr) & (df_sing[CH_BIND] < thresh_bind)
        bind_pos_only = (df_sing[CH_EXPR] < thresh_expr) & (df_sing[CH_BIND] >= thresh_bind)
        double_neg = (df_sing[CH_EXPR] < thresh_expr) & (df_sing[CH_BIND] < thresh_bind)
        
        pct_ll_v = double_neg.mean() * 100
        pct_lr_v = bind_pos_only.mean() * 100
        pct_ul_v = expr_pos_only.mean() * 100
        pct_ur_v = double_pos.mean() * 100
        
        quadrant_str = f"Q1 (UL): {pct_ul_v:.1f}% | Q2 (UR): {pct_ur_v:.1f}%\nQ3 (LL): {pct_ll_v:.1f}% | Q4 (LR): {pct_lr_v:.1f}%"
        axes[1,0].set_title(f"Quadrant Plot\n{quadrant_str}")
        
        # Panel 4: Binding Distribution (APC Histogram)
        sns.kdeplot(t_bind, fill=True, ax=axes[1,1], color='g', label='Sample')
        if neg_log_bind is not None:
            sns.kdeplot(neg_log_bind, fill=True, ax=axes[1,1], color='gray', alpha=0.3, label='Neg Ctrl')
        axes[1,1].axvline(log_thresh_bind, color='r', linestyle='--')
        pct_bind_left = 100 - pct_bind
        axes[1,1].text(0.05, 0.95, f"Unbound: {pct_bind_left:.1f}%\nBound: {pct_bind:.1f}%", 
                       transform=axes[1,1].transAxes, fontsize=10, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        axes[1,1].set_xlabel("Log10 APC (Binding)")
        axes[1,1].set_title(f"Binding Dist (FC: {fc_bind:.1f}x)")
        axes[1,1].legend()
        
        plt.tight_layout()
        plt.savefig(os.path.join(individual_dir, f"{clean_name}_analysis.png"))
        plt.close()
        
        # 2. Publication Plot
        if len(df_sing) > 0:
            plt.figure(figsize=(6, 6))
            plt.hexbin(t_bind, t_expr, gridsize=100, cmap='jet', mincnt=1, bins='log')
            plt.axvline(log_thresh_bind, color='k', linestyle='--')
            plt.axhline(log_thresh_expr, color='k', linestyle='--')
            plt.xlabel("Log10 APC-A (Binding)")
            plt.ylabel("Log10 FITC-A (Expression)")
            plt.title(f"{clean_name}")
            
            textstr = '\n'.join((
                f'Double+: {pct_ur_v:.1f}%',
                f'Bind+ Only (LR): {pct_lr_v:.1f}%',
                f'Expr+ Only (UL): {pct_ul_v:.1f}%'
            ))
            props = dict(boxstyle='round', facecolor='white', alpha=0.8)
            plt.gca().text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
                    verticalalignment='top', bbox=props)
                    
            plt.tight_layout()
            plt.savefig(os.path.join(pub_dir, f"{clean_name}_iso.png"), dpi=200)
            plt.close()
            
            # 3. Comparison vs Pos (Binding)
            plt.figure(figsize=(8, 5))
            if neg_log_bind is not None:
                sns.kdeplot(neg_log_bind, fill=True, color='lightgray', label='Neg Ctrl', alpha=0.5)
            if pos_log_bind is not None:
                sns.kdeplot(pos_log_bind, color='blue', label='Pos Ctrl', linestyle='--', linewidth=2)
                
            sns.kdeplot(t_bind, fill=True, color='green', label=clean_name, alpha=0.4)
            plt.axvline(log_thresh_bind, color='r', linestyle='--', label='99.9% NC Thresh')
            
            plt.xlabel("Log10 APC-A (Binding)")
            plt.title(f"Binding Comparison: {clean_name}\n(Ratio vs Pos: {fc_vs_pos:.2f})")
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(comp_bind_dir, f"Comp_Bind_{clean_name}.png"))
            plt.close()

            # 4. Comparison vs Pos (Expression)
            plt.figure(figsize=(8, 5))
            if neg_log_expr is not None:
                sns.kdeplot(neg_log_expr, fill=True, color='lightgray', label='Neg Ctrl', alpha=0.5)
            if pos_log_expr is not None:
                sns.kdeplot(pos_log_expr, color='blue', label='Pos Ctrl', linestyle='--', linewidth=2)
                
            sns.kdeplot(t_expr, fill=True, color='purple', label=clean_name, alpha=0.4)
            plt.axvline(log_thresh_expr, color='r', linestyle='--', label='99.9% NC Thresh')
            
            plt.xlabel("Log10 FITC-A (Expression)")
            plt.title(f"Expression Comparison: {clean_name}\n(FC: {fc_expr:.1f}x)")
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(comp_expr_dir, f"Comp_Expr_{clean_name}.png"))
            plt.close()

    # Post-Calculation Normalization
    df_stats = pd.DataFrame(summary_stats)
    
    # 5. Calculate Normalized Mean and Median Ratios
    # Identify Positive Control to normalize ratios
    pos_mask = df_stats['Filename'].apply(lambda x: any(p in x for p in POS_CONTROL_PATTERNS))
    if pos_mask.any():
        pos_mean_rat_ref = df_stats.loc[pos_mask, "Pos Mean Ratio"].median() # Use median of PosCtrls if multiple
        pos_med_rat_ref = df_stats.loc[pos_mask, "Pos Med Ratio"].median()
    else:
        # Fallback if no Pos Ctrl: Normalize by maximum ratio found to keep scale 0-1
        pos_mean_rat_ref = df_stats["Pos Mean Ratio"].max()
        pos_med_rat_ref = df_stats["Pos Med Ratio"].max()
        
    pos_mean_rat_ref = max(pos_mean_rat_ref, 1e-9) # Avoid division by zero
    pos_med_rat_ref = max(pos_med_rat_ref, 1e-9)

    df_stats["Norm Pos Mean Ratio"] = df_stats["Pos Mean Ratio"] / pos_mean_rat_ref
    df_stats["Norm Pos Med Ratio"] = df_stats["Pos Med Ratio"] / pos_med_rat_ref
    
    # Save Stats
    df_stats.to_csv(os.path.join(output_dir, "summary_stats.csv"), index=False)
    print(f"\nAnalysis complete. Saved to {os.path.join(output_dir, 'summary_stats.csv')}")
    
    # ---- AGGREGATE PLOTS ----
    generate_aggregate_plots(df_stats, ridge_data, thresh_bind, thresh_expr, agg_dir)
    
    # Report
    generate_report(df_stats, output_dir)

def generate_aggregate_plots(df_stats, ridge_data, thresh_bind, thresh_expr, output_dir):
    # Re-implementing the aggregation plots to ensure they exist
    print("Generating aggregate plots...")
    
    if not ridge_data:
        return

    df_ridge = pd.DataFrame(ridge_data)

    # 0. Ridgeline Plot (Binding & Expression) - Grouped by Prefix
    # Identify Groups (First word of filename)
    df_ridge['Group'] = df_ridge['Sample'].apply(lambda x: x.split()[0] if ' ' in x else 'Misc')
    groups = df_ridge['Group'].unique()
    
    # Iterate for both Binding and Expression
    metrics = {
        "Binding": {"col": "LogBinding", "mfi": "Bind MFI", "thresh": thresh_bind, "color": "cubehelix", "x_lab": "Log10 APC-A (Binding)"},
        "Expression": {"col": "LogExpression", "mfi": "Expr MFI", "thresh": thresh_expr, "color": "mako", "x_lab": "Log10 FITC-A (Expression)"}
    }

    for m_name, m_info in metrics.items():
        for group in groups: # Grouped Plots
            # print(f"Generating {m_name} Ridgeline for Group: {group}...")
            df_group = df_ridge[df_ridge['Group'] == group]
            
            # Filter Middle 95% PER SAMPLE
            df_filtered = pd.DataFrame()
            for sample in df_group['Sample'].unique():
                d_sample = df_group[df_group['Sample'] == sample]
                low = d_sample[m_info["col"]].quantile(0.025)
                high = d_sample[m_info["col"]].quantile(0.975)
                d_filt = d_sample[(d_sample[m_info["col"]] >= low) & (d_sample[m_info["col"]] <= high)]
                df_filtered = pd.concat([df_filtered, d_filt])
            
            if df_filtered.empty: continue

            # Sort within group
            group_stats = df_stats[df_stats['Filename'].isin(df_filtered['Sample'].unique())]
            order_mfi = group_stats.sort_values(m_info["mfi"], ascending=False)["Filename"].tolist()
            
            # Dynamic X-limits
            x_min = df_filtered[m_info["col"]].min()
            x_max = df_filtered[m_info["col"]].max()
            x_pad = (x_max - x_min) * 0.1
            x_lims = (x_min - x_pad, x_max + x_pad*2)

            pal = sns.cubehelix_palette(len(order_mfi), rot=-.25, light=.7) if m_info["color"] == "cubehelix" else sns.color_palette("mako", len(order_mfi))
            
            g = sns.FacetGrid(df_filtered, row="Sample", hue="Sample", aspect=10, height=0.8, palette=pal, row_order=order_mfi)
            g.map(sns.kdeplot, m_info["col"], clip_on=False, fill=True, alpha=1, lw=1.5, bw_adjust=.2)
            g.map(sns.kdeplot, m_info["col"], clip_on=False, color="w", lw=2, bw_adjust=.2)
            g.map(plt.axhline, y=0, lw=2, clip_on=False)
            
            def label(x, color, label):
                ax = plt.gca()
                ax.text(0, .2, label, fontweight="bold", color=color,
                        ha="left", va="center", transform=ax.transAxes, fontsize=10)
                ax.set_facecolor('none')
            
            g.map(label, m_info["col"])
            
            g.set(xlim=x_lims)
            g.fig.subplots_adjust(hspace=-.5)
            g.set_titles("")
            g.set(yticks=[])
            g.despine(bottom=True, left=True)
            g.set_xlabels(f"{m_info['x_lab']} - Group {group}")
            
            t_thresh = np.log10(max(1, m_info["thresh"]))
            for ax, name in zip(g.axes.flatten(), order_mfi):
                ax.axvline(t_thresh, color='r', alpha=0.5, linestyle='--')
                # Calculate % on right
                d_sample_all = df_group[df_group['Sample'] == name]
                pct_right = (d_sample_all[m_info["col"]] > t_thresh).mean() * 100
                
                # Check metric title to change Ridgeline labels
                label_text = f"Expressed: {pct_right:.1f}%" if m_name == "Expression" else f"Bound: {pct_right:.1f}%"
                
                ax.text(0.95, 0.2, label_text, 
                       transform=ax.transAxes, fontsize=9, color='red',
                       ha='right', va='bottom', fontweight='bold')
                
            out_name = f"Aggregate_Ridgeline_{m_name}_{group}.png"
            g.savefig(os.path.join(output_dir, out_name))
            plt.close()

        # All Samples Plot
        out_name = os.path.join(output_dir, f"Aggregate_Ridgeline_{m_name}_All.png")
        order = df_stats.sort_values(m_info["mfi"], ascending=False)["Filename"].tolist()
        pal = sns.cubehelix_palette(len(order), rot=-.25, light=.7) if m_info["color"] == "cubehelix" else sns.color_palette("mako", len(order))
        
        g = sns.FacetGrid(df_ridge, row="Sample", hue="Sample", aspect=15, height=0.6, palette=pal, row_order=order)
        g.map(sns.kdeplot, m_info["col"], clip_on=False, fill=True, alpha=1, lw=1.5, bw_adjust=.2)
        g.map(sns.kdeplot, m_info["col"], clip_on=False, color="w", lw=2, bw_adjust=.2)
        g.map(plt.axhline, y=0, lw=2, clip_on=False)
        
        def label(x, color, label):
            ax = plt.gca()
            ax.text(0, .2, label, fontweight="bold", color=color,
                    ha="left", va="center", transform=ax.transAxes, fontsize=9)
            ax.set_facecolor('none')
            
        g.map(label, m_info["col"])
        g.set_titles("")
        g.set(yticks=[])
        g.despine(bottom=True, left=True)
        g.fig.subplots_adjust(hspace=-.5)
        g.set_xlabels(m_info['x_lab'])
        
        t_thresh = np.log10(max(1, m_info["thresh"]))
        for ax, name in zip(g.axes.flatten(), order):
            ax.axvline(t_thresh, color='r', alpha=0.5, linestyle='--')
            d_sample_all = df_ridge[df_ridge['Sample'] == name]
            pct_right = (d_sample_all[m_info["col"]] > t_thresh).mean() * 100
            
            label_text = f"Expressed: {pct_right:.1f}%" if m_name == "Expression" else f"Bound: {pct_right:.1f}%"
            
            ax.text(0.95, 0.2, label_text, 
                   transform=ax.transAxes, fontsize=9, color='red',
                   ha='right', va='bottom', fontweight='bold')
            
        g.savefig(out_name)
        plt.close()

    # 2. Fold Change Bar Plot (Binding & Expression)
    for m_name, col_name, title, hue_col in [("Binding", "Bind Fold Change", "Binding Fold Change (vs Neg)", "Bind+ %"), 
                                             ("Expression", "Expr Fold Change", "Expression Fold Change (vs Neg)", "Expr+ %")]:
        plt.figure(figsize=(12, 6))
        df_sorted = df_stats.sort_values(col_name, ascending=False)
        
        # Determine max value for normalization based on Positive Controls
        pos_ctrl_mask = df_stats['Filename'].apply(lambda x: any(p in x for p in POS_CONTROL_PATTERNS))
        if pos_ctrl_mask.any():
            vmax = df_stats.loc[pos_ctrl_mask, hue_col].max()
        else:
            vmax = df_stats[hue_col].max()
        vmax = max(vmax, 1.0) # Ensure no divide-by-zero
        
        norm = plt.Normalize(0, vmax)
        cmap = plt.cm.get_cmap("viridis" if m_name == "Binding" else "magma")
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        
        # Map values to colors
        colors = [cmap(norm(val)) for val in df_sorted[hue_col]]
        
        sns.barplot(data=df_sorted, x="Filename", y=col_name, palette=colors)
        
        plt.xticks(rotation=90)
        plt.axhline(1, color='k', linestyle='--', label='No Change')
        plt.title(title)
        plt.ylabel("Fold Change")
        cbar = plt.colorbar(sm, ax=plt.gca())
        cbar.set_label(f"Positive Events (%)")
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"Aggregate_FoldChange_{m_name}.png"))
        plt.close()

    # 3. Stain Index Bar Plot (Binding & Expression)
    for m_name, col_name, title, hue_col in [("Binding", "Bind Stain Index", "Binding Stain Index", "Bind+ %"), 
                                             ("Expression", "Expr Stain Index", "Expression Stain Index", "Expr+ %")]:
        plt.figure(figsize=(12, 6))
        df_sorted = df_stats.sort_values(col_name, ascending=False)
        
        pos_ctrl_mask = df_stats['Filename'].apply(lambda x: any(p in x for p in POS_CONTROL_PATTERNS))
        if pos_ctrl_mask.any():
            vmax = df_stats.loc[pos_ctrl_mask, hue_col].max()
        else:
            vmax = df_stats[hue_col].max()
        vmax = max(vmax, 1.0)
        
        norm = plt.Normalize(0, vmax)
        cmap = plt.cm.get_cmap("viridis" if m_name == "Binding" else "magma")
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        
        colors = [cmap(norm(val)) for val in df_sorted[hue_col]]
        
        sns.barplot(data=df_sorted, x="Filename", y=col_name, palette=colors)
        
        plt.xticks(rotation=90)
        plt.title(title)
        plt.ylabel("Stain Index")
        cbar = plt.colorbar(sm, ax=plt.gca())
        cbar.set_label(f"Positive Events (%)")
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"Aggregate_StainIndex_{m_name}.png"))
        plt.close()
    
    # 4. Bar plot of Double Positives
    plt.figure(figsize=(12, 6))
    df_sorted = df_stats.sort_values("Double+ %", ascending=False)
    
    pos_ctrl_mask = df_stats['Filename'].apply(lambda x: any(p in x for p in POS_CONTROL_PATTERNS))
    if pos_ctrl_mask.any():
        vmax = df_stats.loc[pos_ctrl_mask, "Double+ %"].max()
    else:
        vmax = df_stats["Double+ %"].max()
    vmax = max(vmax, 1.0)
    
    norm = plt.Normalize(0, vmax)
    cmap = plt.cm.get_cmap("Blues")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    
    colors = [cmap(norm(val)) for val in df_sorted["Double+ %"]]
    
    sns.barplot(data=df_sorted, x="Filename", y="Double+ %", palette=colors)
    
    plt.xticks(rotation=90)
    plt.title("Double Positive Population % by Sample")
    cbar = plt.colorbar(sm, ax=plt.gca())
    cbar.set_label("Double+ (%)")
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Aggregate_DoublePos_Bar.png"))
    plt.close()
    
    # 5. Scatter MFI
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=df_stats, x="Bind MFI", y="Expr MFI", hue="Filename", s=100)
    plt.axvline(thresh_bind, color='r', linestyle='--', label='99.9% NC Bind Thresh')
    plt.axhline(thresh_expr, color='r', linestyle='--', label='99.9% NC Expr Thresh')
    plt.title("Mean Fluorescence Intensity: Binding vs Expression")
    plt.xlabel("Binding MFI (APC)")
    plt.ylabel("Expression MFI (FITC)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Aggregate_MFI_Scatter.png"))
    plt.close()

    # 6. Normalized Positive Ratios (Binding Effectiveness)
    # Correcting the path to be at the top level of FCS_Analysis_Results_3, rather than inside Aggregate_Plots
    pos_ratio_dir = os.path.join(os.path.dirname(output_dir), "Positive_Ratios")
    if not os.path.exists(pos_ratio_dir):
        os.makedirs(pos_ratio_dir)
    pos_ctrl_mask = df_stats['Filename'].apply(lambda x: any(p in x for p in POS_CONTROL_PATTERNS))
    if pos_ctrl_mask.any():
        vmax_dp = df_stats.loc[pos_ctrl_mask, "Double+ %"].max()
    else:
        vmax_dp = df_stats["Double+ %"].max()
    vmax_dp = max(vmax_dp, 1.0)
        
    for m_col, title in [("Pos Mean Ratio", "Raw Pos Mean Ratio (Bind/Expr)"),
                         ("Pos Med Ratio", "Raw Pos Median Ratio (Bind/Expr)"),
                         ("Norm Pos Mean Ratio", "Normalized Pos Mean Ratio (Bind/Expr vs PosCtrl)"), 
                         ("Norm Pos Med Ratio", "Normalized Pos Median Ratio (Bind/Expr vs PosCtrl)")]:
        plt.figure(figsize=(12, 6))
        df_sorted = df_stats.sort_values(m_col, ascending=False)
        
        # Color scale based on Double+ % to give more context than just height
        norm = plt.Normalize(0, vmax_dp)
        cmap = plt.cm.get_cmap("viridis")
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        
        colors = [cmap(norm(val)) for val in df_sorted["Double+ %"]]
        
        sns.barplot(data=df_sorted, x="Filename", y=m_col, palette=colors)
        
        if "Norm" in m_col:
            plt.axhline(1, color='k', linestyle='--', label='Positive Control Baseline')
            
        plt.xticks(rotation=90)
        plt.title(title)
        plt.ylabel("Ratio (Sample / PosCtrl)" if "Norm" in m_col else "Ratio (Bind / Expr)")
        
        cbar = plt.colorbar(sm, ax=plt.gca())
        cbar.set_label("Double+ Population (%)")
        
        plt.tight_layout()
        plt.savefig(os.path.join(pos_ratio_dir, f"Aggregate_{m_col.replace(' ', '_')}.png"))
        plt.close()

def generate_report(df, output_dir):
    """Generates a text report summarizing the findings."""
    
    report_path = os.path.join(output_dir, "Analysis_Report.md")
    
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
        f.write(f"**Average Gating Efficiency**: {df['% Gated'].mean():.1f}%\n")
        f.write("\n> **Note**: All analyses are performed on **Singlet** events only. The Singlet Gate removes doublets (clumps of cells) to ensure accurate MFI quantification per cell.\n\n")
        
        f.write("## 1. Executive Summary\n")
        f.write(f"- **Highest Binding (APC)**: {top_binding['Filename']} (MFI: {top_binding['Bind MFI']:.0f}, Ratio vs Pos: {top_binding['Bind FC vs Pos Ctrl']:.2f})\n")
        f.write(f"- **Highest Expression (FITC)**: {top_expr['Filename']} (MFI: {top_expr['Expr MFI']:.0f}, %+: {top_expr['Expr+ %']:.1f}%)\n")
        
        # Best Ratio Effectiveness
        top_ratio = df.sort_values("Norm Pos Med Ratio", ascending=False).iloc[0]
        f.write(f"- **Highest Normalized Positive Ratio (Effectiveness)**: {top_ratio['Filename']} ({top_ratio['Norm Pos Med Ratio']:.2f}x Pos Ctrl)\n")
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
            avg_vs_pos = sub['Bind FC vs Pos Ctrl'].mean()
            avg_norm_ratio = sub['Norm Pos Med Ratio'].mean()
            
            f.write(f"### Group {g}\n")
            f.write(f"- **Samples**: {len(sub)}\n")
            f.write(f"- **Avg Binding Fold Change (vs Neg)**: {avg_fc_bind:.2f}x\n")
            f.write(f"- **Avg Expression Fold Change (vs Neg)**: {avg_fc_expr:.2f}x\n")
            f.write(f"- **Avg Binding Ratio (vs Pos)**: {avg_vs_pos:.2f} (1.0 = Same as Pos)\n")
            f.write(f"- **Avg Binding Efficiency (Binding / Expression ratio)**: {avg_eff:.2f}\n")
            f.write(f"- **Avg Normalized Pos Median Ratio**: {avg_norm_ratio:.2f}\n")
            f.write(f"- **Avg Stain Index**: {avg_si:.2f}\n")
            f.write(f"- **Avg % of Pos Ctrl**: {avg_pct_pos:.1f}%\n")
            f.write(f"- **Observation**: {stats_observation(avg_fc_bind, avg_fc_expr, avg_eff, avg_vs_pos)}\n\n")
            
        f.write("## 3. Explanation of Metrics\n")
        f.write("- **Fold Change (vs Neg)**: The Median MFI of the sample divided by the Median MFI of the Negative Control. A value of 1.0 means it is identical to background.\n")
        f.write("- **Binding Efficiency**: The ratio of Binding Fold Change to Expression Fold Change. Helps identify if a high binding signal is just due to abnormally high expression.\n")
        f.write("- **Stain Index**: A measure of separation between the positive population and the negative population. Calculated as: `(Sample Median - Neg Median) / (2 * robust SD of Neg)`.\n")
        f.write("- **Pos Mean / Median Ratio (Raw)**: For events strictly above the 99.9% Negative Control threshold, this is the ratio of their Binding level to their Expression level (`Pos Bind / Pos Expr`). Validates how effectively the expressed protein binds the target.\n")
        f.write("- **Normalized Pos Mean / Median Ratio**: The Raw Ratio scaled against the Positive Control's ratio. A value of `1.0` means the sample's binding-to-expression effectiveness perfectly matches the Positive Control. `>1.0` is better than Pos Ctrl, `<1.0` is worse.\n")
        f.write("- **Double+ %**: The percentage of *all* events in the sample that fell into the upper-right quadrant (i.e., they expressed *and* bound above the 99.9% Negative Control thresholds).\n")

def stats_observation(fc_bind, fc_expr, efficiency, vs_pos=0):
    obs = []
    
    # Comparison vs Positive Control
    if vs_pos > 1.2:
        obs.append("**Better binder than Positive Control.**")
    elif vs_pos > 0.8:
        obs.append("Binding comparable to Positive Control.")
    elif vs_pos > 0.1:
        obs.append("Binding weaker than Positive Control.")
    
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
        if efficiency > 1.0: 
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

