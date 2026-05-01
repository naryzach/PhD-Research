import os
import sys
import glob
import argparse
import warnings
import shutil
from collections import defaultdict

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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy import stats
from scipy.stats import gaussian_kde
from matplotlib.path import Path
from scipy.spatial import ConvexHull

warnings.filterwarnings("ignore")

# ---- CONFIGURATION ----
NEG_CONTROL_PATTERNS = ["NC", "Negative Control"]
SKIP_PATTERNS = ["EtOH", "Decon"]

# Channels
CH_FSC_A = 'FSC-A'
CH_SSC_A = 'SSC-A'
CH_FITC = 'FITC-A' # Expression (TIMP3)
CH_APC = 'APC-A'   # Binding (His-tagged)
CH_PE = 'PE-A'     # Binding (FLAG-tagged)

# Origin Debris Filter
MIN_FSC = 500000
MIN_SSC = 20000
NC_CUTOFF_PERCENTILE = 99.5
UPPER_FILTER_PERCENTILE = 95
GATE_TARGET_FRACTION = 0.90

# Plot settings
SNS_STYLE = "whitegrid"
sns.set_style(SNS_STYLE)
plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 11,
    'figure.titlesize': 20
})

# ---- HELPER FUNCTIONS ----

def load_fcs(file_path):
    try:
        meta, data = fcsparser.parse(file_path, reformat_meta=True)
        return meta, data
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None, None

def simplify_polygon_vw(points, num_points=5):
    pts = list(points)
    if np.allclose(pts[0], pts[-1]):
        pts.pop()
    while len(pts) > num_points:
        min_area = float('inf')
        min_idx = -1
        for i in range(len(pts)):
            p_prev = pts[i - 1]
            p_curr = pts[i]
            p_next = pts[(i + 1) % len(pts)]
            area = 0.5 * abs(p_prev[0]*(p_curr[1] - p_next[1]) + 
                             p_curr[0]*(p_next[1] - p_prev[1]) + 
                             p_next[0]*(p_prev[1] - p_curr[1]))
            if area < min_area:
                min_area = area
                min_idx = i
        pts.pop(min_idx)
    return np.array(pts)

def learn_pentagon_gate(df, fsc_col, ssc_col, fraction=GATE_TARGET_FRACTION):
    origin_filter = (df[fsc_col] > MIN_FSC) & (df[ssc_col] > MIN_SSC)
    core_df = df[origin_filter]
    if len(core_df) < 100:
        core_df = df
    x = core_df[fsc_col].values
    y = core_df[ssc_col].values
    fsc_upper = np.percentile(x, UPPER_FILTER_PERCENTILE)
    ssc_upper = np.percentile(y, UPPER_FILTER_PERCENTILE)
    core_idx = (x < fsc_upper) & (y < ssc_upper)
    x_core, y_core = x[core_idx], y[core_idx]
    
    n_fit = min(10000, len(x_core))
    idx_fit = np.random.choice(len(x_core), n_fit, replace=False)
    fit_points = np.vstack([x_core[idx_fit], y_core[idx_fit]])
    kernel = gaussian_kde(fit_points)
    
    all_x, all_y = core_df[fsc_col].values, core_df[ssc_col].values
    n_score = min(20000, len(all_x))
    idx_score = np.random.choice(len(all_x), n_score, replace=False)
    score_points = np.vstack([all_x[idx_score], all_y[idx_score]])
    scores_sub = kernel(score_points)
    
    low_pct, high_pct = 0.0, 100.0
    best_path, best_captured = None, 0
    
    for _ in range(25):
        mid_pct = (low_pct + high_pct) / 2
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
            high_pct = mid_pct
            continue
        pent_pts = simplify_polygon_vw(hull_pts, 5)
        pent_pts = np.vstack((pent_pts, pent_pts[0]))
        path = Path(pent_pts)
        captured = path.contains_points(score_points.T).mean()
        if abs(captured - fraction) < abs(best_captured - fraction):
            best_captured = captured
            best_path = path
        if captured > fraction:
            low_pct = mid_pct
        else:
            high_pct = mid_pct
        if abs(captured - fraction) < 0.005:
            break
            
    if best_path is None:
        xmin, xmax = x.min(), x.max()
        ymin, ymax = y.min(), y.max()
        poly_verts = [(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)]
        best_path = Path(poly_verts)
    return best_path

def apply_polygon_gate(df, path, fsc_col, ssc_col, plot_ax=None, title="Polygon Gate"):
    if len(df) == 0: return df
    points = np.vstack([df[fsc_col], df[ssc_col]]).T
    poly_mask = path.contains_points(points)
    origin_mask = (df[fsc_col] > MIN_FSC) & (df[ssc_col] > MIN_SSC)
    mask = poly_mask & origin_mask
    if plot_ax:
        plot_ax.hexbin(df[fsc_col], df[ssc_col], gridsize=100, cmap='inferno', mincnt=1, bins='log')
        plot_ax.axvline(MIN_FSC, color='gray', linestyle=':', alpha=0.5)
        plot_ax.axhline(MIN_SSC, color='gray', linestyle=':', alpha=0.5)
        import matplotlib.patches as patches
        patch = patches.PathPatch(path, facecolor='none', edgecolor='cyan', lw=2, linestyle='--')
        plot_ax.add_patch(patch)
        plot_ax.set_title(f"{title} ({mask.sum()}/{len(df)} evts)")
    return df[mask]

def classify_sample(filename):
    if any(s in filename for s in SKIP_PATTERNS):
        return "skip", None, None, None
    if any(n in filename for n in NEG_CONTROL_PATTERNS):
        return "nc", "NC", None, None
        
    tag = "unknown"
    if "ADAM10" in filename or "Arly" in filename:
        tag = "flag"
    elif "ADAM17" in filename or "Sino" in filename or "Enzo" in filename:
        tag = "his"
        
    # Extract metadata: TIMP 3-TARGET(CONC-SOURCE)
    target = "Unknown"
    conc = "Unknown"
    source = "Unknown"
    
    if "-" in filename and "(" in filename and ")" in filename:
        # e.g. "TIMP 3-ADAM17(Low-FT).fcs"
        try:
            target_part = filename.split("-", 1)[1].split("(")[0].strip()
            target = target_part
            details = filename.split("(")[1].split(")")[0]
            if "-" in details:
                conc, source = details.split("-", 1)
            else:
                source = details
                if "Enzo" in source:
                    conc = "Med"
                else:
                    conc = "None"
                    
            if conc == "Med": conc = "Medium"
            elif conc == "Hi": conc = "High"
            
            if "FXN" in source:
                source = source.replace("FXN", "Fraction ")
            if "FT" in source:
                source = source.replace("FT", "Flow Through")
        except:
            pass
            
    return tag, target, conc, source

def copy_and_rename_files(input_dir, output_dir):
    renamed_dir = os.path.join(output_dir, "Renamed_FCS")
    os.makedirs(renamed_dir, exist_ok=True)
    
    file_info = []
    
    for f in glob.glob(os.path.join(input_dir, "*.fcs")):
        fname = os.path.basename(f)
        tag, target, conc, source = classify_sample(fname)
        if tag == "skip": continue
        
        # Strip well prefix (e.g., "A04 TIMP 3-..." -> "TIMP 3-...")
        parts = fname.split(" ", 1)
        clean_name = parts[1] if len(parts) > 1 else fname
        
        if tag == "nc":
            clean_name = clean_name.replace("NC", "Negative_Control")
            dest_subdir = os.path.join(renamed_dir, "Controls")
        elif tag == "his":
            dest_subdir = os.path.join(renamed_dir, "His_Tagged")
        elif tag == "flag":
            dest_subdir = os.path.join(renamed_dir, "Flag_Tagged")
        else:
            dest_subdir = os.path.join(renamed_dir, "Unknown")
            
        os.makedirs(dest_subdir, exist_ok=True)
        dest_path = os.path.join(dest_subdir, clean_name)
        shutil.copy2(f, dest_path)
        
        file_info.append({
            "OriginalPath": f,
            "CleanName": clean_name.replace(".fcs", ""),
            "Tag": tag,
            "Target": target,
            "Conc": conc,
            "Source": source,
            "DestPath": dest_path
        })
        
    return pd.DataFrame(file_info)

def process_directory(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    print("Pre-processing and classifying files...")
    df_files = copy_and_rename_files(input_dir, output_dir)
    
    if df_files.empty:
        print("No FCS files found to process.")
        return
        
    # Directories
    dirs = {
        "His_Tagged": os.path.join(output_dir, "His_Tagged"),
        "Flag_Tagged": os.path.join(output_dir, "Flag_Tagged"),
        "Controls": os.path.join(output_dir, "Controls"),
        "Combined": os.path.join(output_dir, "Combined")
    }
    for tag_grp in ["His_Tagged", "Flag_Tagged", "Controls"]:
        os.makedirs(os.path.join(dirs[tag_grp], "Individual_Plots"), exist_ok=True)
        if tag_grp != "Controls":
            os.makedirs(os.path.join(dirs[tag_grp], "Concentration_Comparison"), exist_ok=True)
            os.makedirs(os.path.join(dirs[tag_grp], "Aggregate_Plots"), exist_ok=True)
    os.makedirs(dirs["Combined"], exist_ok=True)

    # 1. Gate Learning on NC
    nc_files = df_files[df_files["Tag"] == "nc"]["DestPath"].tolist()
    if not nc_files:
        print("Error: No Negative Controls found.")
        return
        
    print("Learning gate from NC...")
    _, d_nc_learn = load_fcs(nc_files[0])
    pentagon_path = learn_pentagon_gate(d_nc_learn, CH_FSC_A, CH_SSC_A)
    
    # 2. Gate all NCs
    nc_dfs = []
    for f in nc_files:
        _, d = load_fcs(f)
        if d is not None and CH_FITC in d.columns:
            d_gated = apply_polygon_gate(d, pentagon_path, CH_FSC_A, CH_SSC_A)
            d_gated = d_gated[(d_gated[CH_FITC] > 0) & (d_gated[CH_APC] > 0) & (d_gated[CH_PE] > 0)]
            nc_dfs.append(d_gated)
    nc_concat = pd.concat(nc_dfs) if nc_dfs else None
    
    # 3. Calculate Spillover (FITC into PE) using His samples
    # His samples have FITC (TIMP3) and APC (His). They have NO PE.
    # Therefore, any PE signal is pure spillover from FITC.
    print("Calculating Spillover Correction (FITC -> PE) from His-tagged samples...")
    his_files = df_files[df_files["Tag"] == "his"]["DestPath"].tolist()
    his_dfs = []
    for f in his_files:
        _, d = load_fcs(f)
        if d is not None:
            d_gated = apply_polygon_gate(d, pentagon_path, CH_FSC_A, CH_SSC_A)
            d_gated = d_gated[(d_gated[CH_FITC] > 0) & (d_gated[CH_APC] > 0) & (d_gated[CH_PE] > 0)]
            his_dfs.append(d_gated)
            
    alpha = 0.0
    if his_dfs:
        his_concat = pd.concat(his_dfs)
        # Find pure FITC expressing cells: high FITC, low APC (so minimal APC->PE spillover)
        thresh_expr_rough = np.percentile(nc_concat[CH_FITC].dropna(), 99.0)
        thresh_apc_rough = np.percentile(nc_concat[CH_APC].dropna(), 90.0)
        
        pure_fitc = his_concat[(his_concat[CH_FITC] > thresh_expr_rough) & (his_concat[CH_APC] < thresh_apc_rough)]
        if len(pure_fitc) > 50:
            med_pe = pure_fitc[CH_PE].median()
            med_fitc = pure_fitc[CH_FITC].median()
            alpha = med_pe / med_fitc if med_fitc > 0 else 0
            print(f"Spillover Coefficient (Alpha): {alpha:.4f}")
            
            # Plot spillover reference
            plt.figure(figsize=(6, 5))
            plt.hexbin(np.log10(np.clip(pure_fitc[CH_FITC], 1, None)), 
                       np.log10(np.clip(pure_fitc[CH_PE], 1, None)), 
                       gridsize=50, cmap='viridis', bins='log')
            plt.xlabel("Log10 FITC-A")
            plt.ylabel("Log10 PE-A (Uncorrected)")
            plt.title(f"Spillover Reference (His-tagged FITC+ APC-)\nAlpha = {alpha:.4f}")
            plt.savefig(os.path.join(output_dir, "Spillover_Correction_Reference.png"))
            plt.close()
        else:
            print("Warning: Could not find enough pure FITC+ APC- cells to calculate spillover. Alpha=0.")
            
    # Apply correction to NCs to calculate accurate thresholds
    CH_PE_CORR = 'PE-A-Corrected'
    nc_concat[CH_PE_CORR] = nc_concat[CH_PE] - alpha * nc_concat[CH_FITC]
    
    thresh_fitc = np.percentile(nc_concat[CH_FITC].dropna(), NC_CUTOFF_PERCENTILE)
    thresh_apc = np.percentile(nc_concat[CH_APC].dropna(), NC_CUTOFF_PERCENTILE)
    thresh_pe = np.percentile(nc_concat[CH_PE_CORR].dropna(), NC_CUTOFF_PERCENTILE)
    
    neg_mfi_fitc = nc_concat[CH_FITC].median()
    neg_mfi_apc = nc_concat[CH_APC].median()
    neg_mfi_pe = nc_concat[CH_PE_CORR].median()
    
    # 4. Analyze all Samples
    summary_stats = []
    ridge_data_his = []
    ridge_data_flag = []
    
    sample_files = df_files[df_files["Tag"].isin(["his", "flag", "nc"])]
    
    for idx, row in sample_files.iterrows():
        clean_name = row["CleanName"]
        tag = row["Tag"]
        print(f"Processing {clean_name} ({tag})...")
        
        _, df = load_fcs(row["DestPath"])
        if df is None: continue
        
        df_sing = apply_polygon_gate(df, pentagon_path, CH_FSC_A, CH_SSC_A)
        df_sing = df_sing[(df_sing[CH_FITC] > 0) & (df_sing[CH_APC] > 0) & (df_sing[CH_PE] > 0)]
        df_sing[CH_PE_CORR] = df_sing[CH_PE] - alpha * df_sing[CH_FITC]
        
        if len(df_sing) == 0: continue
        
        if tag == "his" or tag == "nc":
            ch_bind = CH_APC
            thresh_bind = thresh_apc
            neg_mfi_bind = neg_mfi_apc
            tag_dir = dirs["Controls"] if tag == "nc" else dirs["His_Tagged"]
            bind_lbl = "APC-A"
            tag_name = "NC" if tag == "nc" else "His (APC)"
        else:
            ch_bind = CH_PE_CORR
            thresh_bind = thresh_pe
            neg_mfi_bind = neg_mfi_pe
            tag_dir = dirs["Flag_Tagged"]
            bind_lbl = "PE-A (Corrected)"
            tag_name = "FLAG (PE)"
            
        # Drop <= 0 binding values to prevent a massive spike at 0 on log plots due to clipping
        df_sing = df_sing[df_sing[ch_bind] > 0]
        
        count_gated = len(df_sing)
        if count_gated == 0: continue
        
        expr_pos = df_sing[CH_FITC] > thresh_fitc
        bind_pos = df_sing[ch_bind] > thresh_bind
        double_pos = expr_pos & bind_pos
        
        mfi_fitc = df_sing[CH_FITC].median()
        mfi_bind = df_sing[ch_bind].median()
        
        fc_expr = mfi_fitc / max(neg_mfi_fitc, 1)
        fc_bind = mfi_bind / max(neg_mfi_bind, 1)
        
        bind_med_expr_pos = df_sing[expr_pos][ch_bind].median() if expr_pos.sum() > 0 else 0
        
        if tag != "nc":
            summary_stats.append({
                "Filename": clean_name,
                "Tag": tag_name,
                "Target": row["Target"],
                "Conc": row["Conc"],
                "Source": row["Source"],
                "Gated Events": count_gated,
                "Expr+ %": expr_pos.mean() * 100,
                "Bind+ %": bind_pos.mean() * 100,
                "Double+ %": double_pos.mean() * 100,
                "Expr MFI": mfi_fitc,
                "Bind MFI": mfi_bind,
                "Bind Med (Expr+)": bind_med_expr_pos,
                "Expr Fold Change": fc_expr,
                "Bind Fold Change": fc_bind
            })
            
            # Prepare ridgeline data
            sub = df_sing.sample(n=min(len(df_sing), 2000))
            for _, s_row in sub.iterrows():
                d = {
                    "Sample": clean_name,
                    "LogExpression": np.log10(max(s_row[CH_FITC], 1)),
                    "LogBinding": np.log10(max(s_row[ch_bind], 1))
                }
                if tag == "his": ridge_data_his.append(d)
                else: ridge_data_flag.append(d)
            
        # Plotting
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        apply_polygon_gate(df, pentagon_path, CH_FSC_A, CH_SSC_A, plot_ax=axes[0,0], title="Gating")
        
        t_expr = np.log10(np.clip(df_sing[CH_FITC], 1, None))
        t_bind = np.log10(np.clip(df_sing[ch_bind], 1, None))
        
        sns.kdeplot(t_expr, fill=True, ax=axes[0,1], color='purple')
        axes[0,1].axvline(np.log10(max(thresh_fitc, 1)), color='r', linestyle='--')
        axes[0,1].set_xlabel("Log10 FITC-A (Expression)")
        axes[0,1].set_title(f"Expression (FITC+ {expr_pos.mean()*100:.1f}%)")
        
        axes[1,0].hexbin(t_bind, t_expr, gridsize=100, cmap='jet', mincnt=1, bins='log')
        axes[1,0].axvline(np.log10(max(thresh_bind, 1)), color='k', linestyle='--')
        axes[1,0].axhline(np.log10(max(thresh_fitc, 1)), color='k', linestyle='--')
        axes[1,0].set_xlabel(f"Log10 {bind_lbl} (Binding)")
        axes[1,0].set_ylabel("Log10 FITC-A (Expression)")
        axes[1,0].set_title(f"Double+ {double_pos.mean()*100:.1f}%")
        
        sns.kdeplot(t_bind, fill=True, ax=axes[1,1], color='g')
        axes[1,1].axvline(np.log10(max(thresh_bind, 1)), color='r', linestyle='--')
        axes[1,1].set_xlabel(f"Log10 {bind_lbl} (Binding)")
        axes[1,1].set_title(f"Binding (Bind+ {bind_pos.mean()*100:.1f}%)")
        
        plt.suptitle(f"{clean_name}", fontsize=16)
        plt.tight_layout()
        plt.savefig(os.path.join(tag_dir, "Individual_Plots", f"{clean_name}_analysis.png"))
        plt.close()
        
    df_stats = pd.DataFrame(summary_stats)
    if df_stats.empty:
        print("No stats generated.")
        return
        
    df_stats.to_csv(os.path.join(output_dir, "summary_stats.csv"), index=False)
    
    # 5. Concentration Comparison Plots
    print("Generating concentration comparison plots...")
    conc_order = ["None", "Low", "Medium", "High"]
    df_stats['Conc_Cat'] = pd.Categorical(df_stats['Conc'], categories=conc_order, ordered=True)
    
    for tag_grp, d_tag in [("His (APC)", "His_Tagged"), ("FLAG (PE)", "Flag_Tagged")]:
        d_sub = df_stats[df_stats["Tag"] == tag_grp]
        if d_sub.empty: continue
        
        targets = d_sub["Target"].unique()
        for tgt in targets:
            d_tgt = d_sub[d_sub["Target"] == tgt].dropna(subset=['Conc_Cat']).sort_values('Conc_Cat')
            if d_tgt.empty or len(d_tgt['Conc'].unique()) <= 1:
                continue
                
            plt.figure(figsize=(8, 5))
            sns.lineplot(data=d_tgt, x="Conc_Cat", y="Bind Med (Expr+)", hue="Source", marker="o", markersize=10, linewidth=2)
            plt.title(f"Concentration Effect on Binding: {tgt} ({tag_grp})")
            plt.ylabel("Median Binding MFI (in FITC+ cells)")
            plt.xlabel("Target Concentration")
            plt.tight_layout()
            plt.savefig(os.path.join(dirs[d_tag], "Concentration_Comparison", f"{tgt}_ConcentrationEffect.png"), dpi=300)
            plt.close()

    # 6. Aggregate Bar Plot (Combined)
    plt.figure(figsize=(14, 7))
    # Create a descriptive label combining Filename and Tag (e.g. "Sample Name (His)" or "Sample Name (FLAG)")
    df_stats['Sample_Label'] = df_stats['Filename'] + " (" + df_stats['Tag'].apply(lambda x: x.split()[0]) + ")"
    df_sorted = df_stats.sort_values(["Target", "Conc_Cat", "Filename"], ascending=[True, True, True])
    sns.barplot(data=df_sorted, x="Sample_Label", y="Double+ %", hue="Target", dodge=False, palette="Set2")
    plt.xticks(rotation=90)
    plt.title("Double Positive Population % (Grouped by Target & Concentration)")
    plt.tight_layout()
    plt.savefig(os.path.join(dirs["Combined"], "All_Targets_DoublePos_Bar.png"), dpi=300)
    plt.close()
    
    # 7. Ridgelines
    print("Generating ridgeline plots...")
    def plot_ridge(data, metric_col, title, out_path, x_lbl, t_thresh, stats_df):
        if not data: return
        df_r = pd.DataFrame(data)
        
        df_ridge_filtered = pd.DataFrame()
        for sample in df_r['Sample'].unique():
            d_sample = df_r[df_r['Sample'] == sample]
            low = d_sample[metric_col].quantile(0.025)
            high = d_sample[metric_col].quantile(0.975)
            d_filt = d_sample[(d_sample[metric_col] >= low) & (d_sample[metric_col] <= high)]
            df_ridge_filtered = pd.concat([df_ridge_filtered, d_filt])

        if df_ridge_filtered.empty: return
        
        # Sort by Target -> Conc_Cat -> Filename using the stats_df
        ordered_samples = stats_df.sort_values(["Target", "Conc_Cat", "Filename"], ascending=[True, True, True])["Filename"].tolist()
        order = [s for s in ordered_samples if s in df_ridge_filtered['Sample'].unique()]
        
        x_min_all = df_ridge_filtered[metric_col].min()
        x_max_all = df_ridge_filtered[metric_col].max()
        x_pad_all = (x_max_all - x_min_all) * 0.1
        x_lims_all = (x_min_all - x_pad_all, x_max_all + x_pad_all*2)

        pal = sns.cubehelix_palette(len(order), rot=-.25, light=.7)
        g = sns.FacetGrid(df_ridge_filtered, row="Sample", hue="Sample", aspect=15, height=0.6, palette=pal, row_order=order, hue_order=order)
        g.map(sns.kdeplot, metric_col, clip_on=False, fill=True, alpha=1, lw=1.5, bw_adjust=.2)
        g.map(sns.kdeplot, metric_col, clip_on=False, color="w", lw=2, bw_adjust=.2)
        g.map(plt.axhline, y=0, lw=2, clip_on=False)
        
        def label(x, color, label):
            ax = plt.gca()
            ax.text(0, .2, label, fontweight="bold", color=color, ha="left", va="center", transform=ax.transAxes, fontsize=9)
            ax.set_facecolor('none')
            
        g.map(label, metric_col)
        g.set_titles("")
        g.set(yticks=[], xlim=x_lims_all)
        g.set_ylabels("")
        g.despine(bottom=True, left=True)
        g.fig.subplots_adjust(hspace=-.5, top=0.88)
        g.set_xlabels(x_lbl)
        g.fig.suptitle(title, fontsize=22)
        
        for ax, name in zip(g.axes.flatten(), order):
            ax.axvline(np.log10(max(1, t_thresh)), color='r', alpha=0.5, linestyle='--')
            d_sample_all = df_r[df_r['Sample'] == name]
            pct_right = (d_sample_all[metric_col] > np.log10(max(1, t_thresh))).mean() * 100
            
            label_text = f"Expressed: {pct_right:.1f}%" if "FITC" in x_lbl else f"Bound: {pct_right:.1f}%"
            ax.text(0.95, 0.2, label_text, transform=ax.transAxes, fontsize=9, color='red', ha='right', va='bottom', fontweight='bold')
            
        g.savefig(out_path, dpi=300)
        plt.close()

    plot_ridge(ridge_data_his, "LogBinding", "His-Tagged Binding (APC)", os.path.join(dirs["His_Tagged"], "Aggregate_Plots", "Ridgeline_Binding.png"), "Log10 APC", thresh_apc, df_stats)
    plot_ridge(ridge_data_his, "LogExpression", "His-Tagged Expression (FITC)", os.path.join(dirs["His_Tagged"], "Aggregate_Plots", "Ridgeline_Expression.png"), "Log10 FITC", thresh_fitc, df_stats)
    plot_ridge(ridge_data_flag, "LogBinding", "FLAG-Tagged Binding (PE, Corrected)", os.path.join(dirs["Flag_Tagged"], "Aggregate_Plots", "Ridgeline_Binding.png"), "Log10 PE (Corr)", thresh_pe, df_stats)
    plot_ridge(ridge_data_flag, "LogExpression", "FLAG-Tagged Expression (FITC)", os.path.join(dirs["Flag_Tagged"], "Aggregate_Plots", "Ridgeline_Expression.png"), "Log10 FITC", thresh_fitc, df_stats)

    # 8. Markdown Report
    with open(os.path.join(output_dir, "Analysis_Report.md"), "w") as f:
        f.write("# TIMP3 Target Trial FCS Analysis Report\n\n")
        f.write(f"**Spillover Coefficient (FITC -> PE)**: {alpha:.4f}\n")
        f.write("*(Applied automatically to FLAG-tagged samples to remove FITC bleed-through from the PE channel)*\n\n")
        
        top_his = df_stats[df_stats["Tag"] == "His (APC)"].sort_values("Double+ %", ascending=False)
        top_flag = df_stats[df_stats["Tag"] == "FLAG (PE)"].sort_values("Double+ %", ascending=False)
        
        f.write("## Top Binders: His-Tagged (APC)\n")
        for i, row in top_his.head(3).iterrows():
            f.write(f"- **{row['Filename']}**: {row['Double+ %']:.1f}% Double+ | Bind MFI: {row['Bind MFI']:.0f}\n")
            
        f.write("\n## Top Binders: FLAG-Tagged (PE)\n")
        for i, row in top_flag.head(3).iterrows():
            f.write(f"- **{row['Filename']}**: {row['Double+ %']:.1f}% Double+ | Bind MFI: {row['Bind MFI']:.0f}\n")
            
    print("Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    process_directory(args.input, args.output)
