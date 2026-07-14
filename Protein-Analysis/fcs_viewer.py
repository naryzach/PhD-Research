import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import fcsparser
import s3fs
import io
import os
import re
import sys
import glob
from scipy import stats
from scipy.stats import gaussian_kde
from scipy.spatial import ConvexHull
from matplotlib.path import Path
import tempfile

# Make the sibling aggregate_analysis module importable so the viewer can
# recompute aggregate/selectivity data live instead of reading the static CSVs.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import aggregate_analysis
    AGG_MODULE_OK = True
    _AGG_IMPORT_ERROR = ""
except Exception as _agg_err:  # noqa: BLE001 - surface any import failure in the UI
    aggregate_analysis = None
    AGG_MODULE_OK = False
    _AGG_IMPORT_ERROR = str(_agg_err)

# Shared explanation for the Tukey HSD post-hoc tables (shown under the dropdowns).
TUKEY_NOTE = (
    "These rows are **Tukey HSD post-hoc pairwise comparisons**, shown only when the "
    "one-way **ANOVA** across the groups is significant (p < 0.05). "
    "**Diff** (shown as *meandiff* in some tables) is the difference between the two "
    "groups' mean values of the metric — specifically *mean(group 2) − mean(group 1)* "
    "in that metric's own units, so the sign tells you which group is higher and the "
    "magnitude is the effect size. **p-adj** is the multiple-comparison–adjusted "
    "p-value for that specific pair (significant if < 0.05)."
)

# ---- PAGE CONFIGURATION ----
st.set_page_config(
    page_title="FCS Viewer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ---- CUSTOM CSS FOR RICH AESTHETICS ----
st.markdown("""
    <style>
    /* Main container and background adjustments */
    .stApp {
        background-color: #0E1117;
    }
    
    /* Headers */
    h1, h2, h3 {
        color: #00C0FF !important;
        font-family: 'Segoe UI', sans-serif;
    }
    
    /* Sidebar */
    [data-testid="stSidebar"] {
        background-color: #161B22;
        border-right: 1px solid #30363D;
    }
    
    /* Metrics */
    [data-testid="stMetricValue"] {
        color: #00C0FF;
    }
    
    /* Tabs */
    .stTabs [data-baseweb="tab-list"] {
        gap: 10px;
    }
    .stTabs [data-baseweb="tab"] {
        height: 50px;
        white-space: pre-wrap;
        background-color: #161B22;
        border-radius: 5px;
        color: #C9D1D9;
        font-weight: 600;
    }
    .stTabs [data-baseweb="tab"]:hover {
        background-color: #21262D;
        color: #58A6FF;
    }
    .stTabs [aria-selected="true"] {
        background-color: #00C0FF !important;
        color: #0E1117 !important;
    }
    </style>
""", unsafe_allow_html=True)

# Gating parameters are now managed via the sidebar UI.

# ---- DATA BACKEND INITIALIZATION ----

@st.cache_resource
def get_fs():
    """Initializes a connection to R2 if secrets are present, else returns None (local mode)."""
    if "r2" not in st.secrets:
        return None, "Missing [r2] section in secrets.toml"
        
    try:
        import s3fs
        r2_config = st.secrets["r2"]
        fs = s3fs.S3FileSystem(
            key=r2_config["access_key_id"],
            secret=r2_config["secret_access_key"],
            endpoint_url=r2_config["endpoint_url"],
            client_kwargs={'region_name': 'auto'},
            use_listings_cache=False
        )
        return fs, r2_config.get("bucket_name", "Unknown")
    except ImportError:
        return None, "❌ s3fs not found. Run pip install s3fs"
    except Exception as e:
        return None, f"❌ Connection Error: {e}"

fs, BUCKET = get_fs()
fs_status = "✅ Connected to R2" if fs else BUCKET # BUCKET contains error msg when fs is None

# --- Data Source Selection State ---
if "data_mode" not in st.session_state:
    st.session_state["data_mode"] = "Cloud (R2)" if fs else "Local"

# (Connection Debug moved to sidebar bottom)

# ---- HELPER FUNCTIONS ----

@st.cache_data
def run_anova_p(df, group_col, value_col):
    """Returns the ANOVA p-value for a given value_col grouped by group_col."""
    if df.empty: return np.nan
    groups = [group[value_col].values for name, group in df.groupby(group_col)]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) < 2:
        return np.nan
    try:
        f_stat, p_val = stats.f_oneway(*groups)
        return p_val
    except:
        return np.nan

def run_tukey_summary(df, group_col, value_col):
    """Runs Tukey HSD post-hoc test and returns a dataframe of significant pairs."""
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    if df.empty: return None
    
    df_f = df.dropna(subset=[value_col])
    # Tukey HSD requires groups with > 1 data point for variance estimation
    counts = df_f.groupby(group_col)[value_col].count()
    valid_groups = counts[counts > 1].index.tolist()
    
    if len(valid_groups) < 2:
        return None
        
    df_f = df_f[df_f[group_col].isin(valid_groups)]
    
    try:
        tukey = pairwise_tukeyhsd(endog=df_f[value_col], groups=df_f[group_col], alpha=0.05)
        summary_df = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
        significant = summary_df[summary_df['reject'] == True]
        return significant
    except:
        return None

def load_fcs(file_path):


    """Loads an FCS file and returns metadata and dataframe."""
    try:
        # Check if this is an R2 path
        if fs and BUCKET and file_path.startswith(BUCKET):
            # Create a temporary file with .fcs suffix
            with tempfile.NamedTemporaryFile(suffix=".fcs", delete=False) as tmp:
                tmp_path = tmp.name
            
            try:
                # Download the R2 object to the temporary local file
                fs.get(file_path, tmp_path)
                # Parse the local temporary file
                meta, data = fcsparser.parse(tmp_path, reformat_meta=True)
                return meta, data
            finally:
                # Clean up the temporary file
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)
        else:
            # Standard local file loading
            meta, data = fcsparser.parse(file_path, reformat_meta=True)
            return meta, data
    except Exception as e:
        st.error(f"Error loading {file_path}: {e}")
        return None, None

def get_fcs_folders(base_dir):
    """Finds all subdirectories containing .fcs files."""
    if fs and base_dir.startswith(BUCKET):
        try:
            # R2 doesn't have true directories, but we can browse via prefixes
            all_files = fs.find(base_dir)
            folders = set()
            for f in all_files:
                if f.lower().endswith('.fcs'):
                    folders.add(os.path.dirname(f))
            return sorted(list(folders))
        except:
            return []
            
    if not os.path.exists(base_dir):
        return []
    folders = []
    for root, _, files in os.walk(base_dir):
        if any(f.lower().endswith('.fcs') for f in files):
            folders.append(root)
    return sorted(folders)

def get_files(directory):
    """Finds FCS files in a specific directory."""
    if fs and directory.startswith(BUCKET):
        try:
            # We want files directly in this "directory" (prefix)
            files = fs.ls(directory, detail=False)
            fcs_files = [f for f in files if f.lower().endswith('.fcs')]
            return sorted(fcs_files)
        except:
            return []

    if not os.path.exists(directory):
        return []
    files = glob.glob(os.path.join(directory, "*.fcs"))
    return sorted(files)

def cloud_exists(path):
    """Checks if a path exists locally or in R2."""
    if fs and BUCKET and path.startswith(BUCKET):
        return fs.exists(path)
    return os.path.exists(path)

def cloud_read_csv(path):
    """Reads a CSV file from local disk or R2."""
    if fs and BUCKET and path.startswith(BUCKET):
        with fs.open(path, 'rb') as f:
            return pd.read_csv(f)
    return pd.read_csv(path)

# Control-identification patterns for on-the-fly summary-stat computation
# (mirrors analyze_fcs.py so recalculated stats match the saved sheets).
RECALC_NEG_PATTERNS = ["NC", "Negative Control"]
RECALC_POS_PATTERNS = ["Positive Control", "TIMP 3"]

@st.cache_data(show_spinner=False)
def discover_trial_folders(base_path):
    """Trial folders (local OR R2) containing .fcs files whose names look like a
    trial, e.g. the *_Renamed folders (Target_Date...). Works for both sources."""
    try:
        folders = get_fcs_folders(base_path)
    except Exception:
        folders = []
    out = [f for f in folders if re.match(r'^([^_]+)_(\d{8})', os.path.basename(f.rstrip("/")))]
    return sorted(out)

@st.cache_data(show_spinner=False)
def compute_folder_summary_records(directory, expr_ch, bind_ch, nc_pct, min_fsc, min_ssc, upper_pct, gate_frac):
    """Compute analyze_fcs-equivalent per-file summary_stats records for one raw
    trial folder (local or R2), on the fly from the FCS files.

    Returns a list of dict records in the summary_stats.csv schema (the subset
    consumed by aggregate_analysis). A negative control is required in the folder
    for gating/thresholds; returns [] otherwise, matching analyze_fcs.py.
    """
    files = get_files(directory)
    if not files:
        return []

    neg_files = [f for f in files if any(os.path.basename(f).startswith(p) for p in RECALC_NEG_PATTERNS)]
    if not neg_files:
        return []

    # Learn the pentagon gate from the first NC
    _, d_nc = load_fcs(neg_files[0])
    if d_nc is None or "FSC-A" not in d_nc.columns or "SSC-A" not in d_nc.columns:
        return []
    if expr_ch not in d_nc.columns or bind_ch not in d_nc.columns:
        return []
    pent_path, _ = learn_pentagon_gate(d_nc, "FSC-A", "SSC-A", gate_frac, min_fsc, min_ssc, upper_pct)

    # Concatenate gated negative controls to derive thresholds
    neg_dfs = []
    for f in neg_files:
        _, d = load_fcs(f)
        if d is None:
            continue
        dg = apply_polygon_gate(d, pent_path, "FSC-A", "SSC-A", min_fsc, min_ssc)
        if expr_ch in dg.columns and bind_ch in dg.columns:
            dg = dg[(dg[expr_ch] > 0) & (dg[bind_ch] > 0)]
            if len(dg) > 0:
                neg_dfs.append(dg)
    if not neg_dfs:
        return []
    neg_concat = pd.concat(neg_dfs)
    thresh_expr = np.percentile(neg_concat[expr_ch], nc_pct)
    thresh_bind = np.percentile(neg_concat[bind_ch], nc_pct)

    # Per-file metrics (subset of analyze_fcs.py needed for aggregation)
    records = []
    lbl_expr = expr_ch.replace("-A", "")
    for f in files:
        _, df = load_fcs(f)
        if df is None or expr_ch not in df.columns or bind_ch not in df.columns:
            continue
        dg = apply_polygon_gate(df, pent_path, "FSC-A", "SSC-A", min_fsc, min_ssc)
        dg = dg[(dg[expr_ch] > 0) & (dg[bind_ch] > 0)]
        count_gated = len(dg)
        if count_gated == 0:
            continue

        expr_pos = dg[expr_ch] > thresh_expr
        bind_pos = dg[bind_ch] > thresh_bind
        double_pos = expr_pos & bind_pos
        pct_expr = expr_pos.mean() * 100
        pct_double = double_pos.mean() * 100

        d_pe = dg[expr_pos]
        d_pb = dg[bind_pos]
        pos_mean_expr = d_pe[expr_ch].mean() if len(d_pe) > 0 else 0
        pos_med_expr = d_pe[expr_ch].median() if len(d_pe) > 0 else 0
        pos_mean_bind = d_pb[bind_ch].mean() if len(d_pb) > 0 else 0
        pos_med_bind = d_pb[bind_ch].median() if len(d_pb) > 0 else 0
        pos_mean_ratio = pos_mean_bind / pos_mean_expr if pos_mean_expr > 0 else 0
        pos_med_ratio = pos_med_bind / pos_med_expr if pos_med_expr > 0 else 0

        bind_mean_expr_pos = d_pe[bind_ch].mean() if len(d_pe) > 0 else 0
        bind_med_expr_pos = d_pe[bind_ch].median() if len(d_pe) > 0 else 0
        expr_mean_bind_pos = d_pb[expr_ch].mean() if len(d_pb) > 0 else 0
        expr_med_bind_pos = d_pb[expr_ch].median() if len(d_pb) > 0 else 0

        bind_eff_ratio = pct_double / pct_expr if pct_expr > 0 else 0

        d_dp = dg[double_pos]
        if len(d_dp) > 0:
            sig_bind = (d_dp[bind_ch] - thresh_bind).clip(lower=0.1)
            sig_expr = (d_dp[expr_ch] - thresh_expr).clip(lower=0.1)
            iwb_index = (sig_bind / sig_expr).median()
        else:
            iwb_index = 0

        records.append({
            "Filename": os.path.splitext(os.path.basename(f))[0],
            "Gated Events": count_gated,
            "Expr+ %": pct_expr,
            "Double+ %": pct_double,
            "Pos Mean Ratio": pos_mean_ratio,
            "Pos Med Ratio": pos_med_ratio,
            "Bind Mean (Expr+)": bind_mean_expr_pos,
            "Bind Med (Expr+)": bind_med_expr_pos,
            "Expr Mean (Bind+)": expr_mean_bind_pos,
            "Expr Med (Bind+)": expr_med_bind_pos,
            f"Binding Efficiency (DP/{lbl_expr}+)": bind_eff_ratio,
            "Intensity-Weighted Binding Index": iwb_index,
        })

    # Raw per-file records only (no Norm columns). Normalization is applied
    # separately so it can respond to the PC-keyword setting without reloading FCS.
    return records

def _normalize_folder_records(records, pos_patterns):
    """Add the Norm* columns to one folder's records, normalized to the positive
    control (median of PC rows, else max). Cheap; recomputed whenever the PC
    keyword setting changes, so it need not reload the FCS files."""
    if not records:
        return records
    pos_list = [str(p).upper() for p in pos_patterns] if pos_patterns else [p.upper() for p in RECALC_POS_PATTERNS]
    df_stats = pd.DataFrame(records)
    pos_mask = df_stats["Filename"].apply(lambda x: any(p in str(x).upper() for p in pos_list))

    def _norm(col):
        ref = df_stats.loc[pos_mask, col].median() if pos_mask.any() else df_stats[col].max()
        ref = max(float(ref) if pd.notna(ref) else 0.0, 1e-9)
        return df_stats[col] / ref

    df_stats["Norm Pos Mean Ratio"] = _norm("Pos Mean Ratio")
    df_stats["Norm Pos Med Ratio"] = _norm("Pos Med Ratio")
    df_stats["Norm Bind Mean (Expr+)"] = _norm("Bind Mean (Expr+)")
    df_stats["Norm Bind Med (Expr+)"] = _norm("Bind Med (Expr+)")
    df_stats["Norm Expr Mean (Bind+)"] = _norm("Expr Mean (Bind+)")
    df_stats["Norm Expr Med (Bind+)"] = _norm("Expr Med (Bind+)")
    df_stats["Norm Intensity-Weighted Binding Index"] = _norm("Intensity-Weighted Binding Index")
    return df_stats.to_dict("records")

def recalc_global_aggregate(base_path, excluded_trials_tuple, expr_ch, bind_ch,
                            nc_pct, min_fsc, min_ssc, upper_pct, gate_frac,
                            pos_patterns=None, progress_cb=None):
    """Recompute the global aggregate trial table on the fly from raw FCS in the
    selected data source (Local or R2). Summary stats are computed per trial
    folder, normalized to the positive control (using the PC-keyword setting),
    then aggregated via aggregate_analysis. Excluded folders are skipped.
    Returns a DataFrame equivalent to Aggregate_FCS_Analysis/aggregate_summary.csv."""
    if not AGG_MODULE_OK:
        return pd.DataFrame()

    folders = discover_trial_folders(base_path)
    excluded = set(excluded_trials_tuple)
    rows = []
    total = len(folders)
    for i, folder in enumerate(folders):
        name = os.path.basename(folder.rstrip("/"))
        if name in excluded:
            continue
        m = re.match(r'^([^_]+)_(\d{8})', name)
        if not m:
            continue
        tgt, date = m.group(1).upper(), m.group(2)
        if progress_cb:
            progress_cb(i, total, name)
        # Raw stats are cached on gating+channels; normalization + aggregation
        # are cheap and re-run each time so PC-keyword changes take effect.
        recs = compute_folder_summary_records(folder, expr_ch, bind_ch, nc_pct, min_fsc, min_ssc, upper_pct, gate_frac)
        recs = _normalize_folder_records(recs, pos_patterns)
        if recs:
            rows.extend(aggregate_analysis.aggregate_records(
                recs, tgt, date, pos_patterns=pos_patterns, folder_name=name))
    if progress_cb:
        progress_cb(total, total, "")
    return pd.DataFrame(rows) if rows else pd.DataFrame()

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

@st.cache_data(show_spinner=False)
def learn_pentagon_gate(df, fsc_col="FSC-A", ssc_col="SSC-A", fraction=0.9, min_fsc=500000, min_ssc=20000, upper_pct=95):
    origin_filter = (df[fsc_col] > min_fsc) & (df[ssc_col] > min_ssc)
    core_df = df[origin_filter]
    
    if len(core_df) < 100:
        core_df = df
        
    x = core_df[fsc_col].values
    y = core_df[ssc_col].values
    
    fsc_upper = np.percentile(x, upper_pct)
    ssc_upper = np.percentile(y, upper_pct)
    core_idx = (x < fsc_upper) & (y < ssc_upper)
    
    x_core = x[core_idx]
    y_core = y[core_idx]
    
    n_fit = min(10000, len(x_core))
    if n_fit == 0:
        xmin, xmax = df[fsc_col].min(), df[fsc_col].max()
        ymin, ymax = df[ssc_col].min(), df[ssc_col].max()
        verts = np.array([[xmin, ymin], [xmin, ymax], [xmax, ymax], [xmax, ymin], [xmin, ymin]])
        return Path(verts), verts

    idx_fit = np.random.choice(len(x_core), n_fit, replace=False)
    fit_points = np.vstack([x_core[idx_fit], y_core[idx_fit]])
    
    kernel = gaussian_kde(fit_points)
    
    all_x = core_df[fsc_col].values
    all_y = core_df[ssc_col].values
    
    n_score = min(20000, len(all_x))
    idx_score = np.random.choice(len(all_x), n_score, replace=False)
    score_points = np.vstack([all_x[idx_score], all_y[idx_score]])
    
    scores_sub = kernel(score_points)
    
    target_captured = fraction
    low_pct = 0.0
    high_pct = 100.0
    
    best_path = None
    best_verts = None
    best_captured = 0
    
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
        
        if abs(captured - target_captured) < abs(best_captured - target_captured):
            best_captured = captured
            best_path = path
            best_verts = pent_pts

        if captured > target_captured:
            low_pct = mid_pct
        else:
            high_pct = mid_pct
            
        if abs(captured - target_captured) < 0.005:
            break
            
    if best_path is None:
        xmin, xmax = x.min(), x.max()
        ymin, ymax = y.min(), y.max()
        best_verts = np.array([[xmin, ymin], [xmin, ymax], [xmax, ymax], [xmax, ymin], [xmin, ymin]])
        best_path = Path(best_verts)

    return best_path, best_verts

def apply_polygon_gate(df, path, fsc_col="FSC-A", ssc_col="SSC-A", min_fsc=500000, min_ssc=20000):
    if len(df) == 0:
        return df
    points = np.vstack([df[fsc_col], df[ssc_col]]).T
    poly_mask = path.contains_points(points)
    origin_mask = (df[fsc_col] > min_fsc) & (df[ssc_col] > min_ssc)
    return df[poly_mask & origin_mask]

@st.cache_data(show_spinner=False)
def analyze_folder_controls(directory, fsc_col="FSC-A", ssc_col="SSC-A", expr_col="FITC-A", bind_col="APC-A", nc_percentile=99.5, min_fsc=500000, min_ssc=20000, upper_pct=95, gate_fraction=0.9):
    all_files = get_files(directory)
    neg_files = [f for f in all_files if any(os.path.basename(f).startswith(p) for p in ["NC", "Negative Control"])]
    pos_files = [f for f in all_files if any(os.path.basename(f).startswith(p) for p in ["Positive Control", "TIMP 3"])]
    
    results = {
        "has_nc": False,
        "has_pc": False,
        "pentagon_path": None,
        "pentagon_verts": None,
        "thresh_expr": None,
        "thresh_bind": None,
        "neg_mfi_expr": None,
        "neg_mfi_bind": None,
        "pos_mfi_expr": None,
        "pos_mfi_bind": None,
        "neg_sample": None,
        "pos_sample": None,
    }
    
    if neg_files:
        _, d_neg_learn = load_fcs(neg_files[0])
        if d_neg_learn is not None and fsc_col in d_neg_learn.columns and ssc_col in d_neg_learn.columns:
            path, verts = learn_pentagon_gate(d_neg_learn, fsc_col, ssc_col, gate_fraction, min_fsc, min_ssc, upper_pct)
            results["pentagon_path"] = path
            results["pentagon_verts"] = verts
            
            neg_dfs = []
            for f in neg_files:
                _, d = load_fcs(f)
                if d is not None:
                    d_gated = apply_polygon_gate(d, path, fsc_col, ssc_col)
                    if expr_col in d_gated.columns and bind_col in d_gated.columns:
                        d_gated = d_gated[(d_gated[expr_col] > 0) & (d_gated[bind_col] > 0)]
                        neg_dfs.append(d_gated)
            
            if neg_dfs:
                neg_concat = pd.concat(neg_dfs)
                if len(neg_concat) > 0:
                    results["thresh_expr"] = np.percentile(neg_concat[expr_col], nc_percentile)
                    results["thresh_bind"] = np.percentile(neg_concat[bind_col], nc_percentile)
                    results["neg_mfi_expr"] = neg_concat[expr_col].median()
                    results["neg_mfi_bind"] = neg_concat[bind_col].median()
                    results["neg_sample"] = neg_concat.sample(n=min(5000, len(neg_concat)))
                    results["has_nc"] = True

    if pos_files and results["has_nc"]:
        pos_dfs = []
        for f in pos_files:
            _, d = load_fcs(f)
            if d is not None:
                d_gated = apply_polygon_gate(d, results["pentagon_path"], fsc_col, ssc_col, min_fsc, min_ssc)
                if expr_col in d_gated.columns and bind_col in d_gated.columns:
                    d_gated = d_gated[(d_gated[expr_col] > 0) & (d_gated[bind_col] > 0)]
                    pos_dfs.append(d_gated)
                
        if pos_dfs:
            pos_concat = pd.concat(pos_dfs)
            if len(pos_concat) > 0:
                results["pos_mfi_expr"] = pos_concat[expr_col].median()
                results["pos_mfi_bind"] = pos_concat[bind_col].median()
                results["pos_sample"] = pos_concat.sample(n=min(5000, len(pos_concat)))
                results["has_pc"] = True

    return results

@st.cache_data(show_spinner=False)
def get_folder_aggregate_stats(directory, _ctrls, fsc_col="FSC-A", ssc_col="SSC-A", expr_col="FITC-A", bind_col="APC-A", def_thresh_expr=1.0, def_thresh_bind=1.0, min_fsc=500000, min_ssc=20000, upper_pct=95):
    all_files = get_files(directory)
    stats_list = []
    
    pent_path = _ctrls.get("pentagon_path")
    if not pent_path:
        return pd.DataFrame(), pd.DataFrame()
        
    neg_mfi_expr = _ctrls.get("neg_mfi_expr", 1)
    neg_mfi_bind = _ctrls.get("neg_mfi_bind", 1)
    
    ridge_data = []

    for f in all_files:
        meta, df = load_fcs(f)
        if df is None: continue
        
        df_sing = apply_polygon_gate(df, pent_path, fsc_col, ssc_col, min_fsc, min_ssc)
        if len(df_sing) > 0 and expr_col in df_sing.columns and bind_col in df_sing.columns:
            df_sing = df_sing[(df_sing[expr_col] > 0) & (df_sing[bind_col] > 0)]
            
        if len(df_sing) == 0:
            continue
            
        expr_pos = df_sing[expr_col] > def_thresh_expr
        bind_pos = df_sing[bind_col] > def_thresh_bind
        double_pos = expr_pos & bind_pos
        
        pct_double = double_pos.mean() * 100
        
        mfi_expr = df_sing[expr_col].median()
        mfi_bind = df_sing[bind_col].median()
        
        fc_expr = mfi_expr / max(neg_mfi_expr, 1)
        fc_bind = mfi_bind / max(neg_mfi_bind, 1)
        
        pos_expr_mask = df_sing[expr_pos]
        pos_bind_mask = df_sing[bind_pos]
        
        pos_med_expr = pos_expr_mask[expr_col].median() if len(pos_expr_mask) > 0 else 0
        pos_med_bind = pos_bind_mask[bind_col].median() if len(pos_bind_mask) > 0 else 0
        pos_med_ratio = (pos_med_bind / pos_med_expr) if pos_med_expr > 0 else 0
        
        pos_mean_expr = pos_expr_mask[expr_col].mean() if len(pos_expr_mask) > 0 else 0
        pos_mean_bind = pos_bind_mask[bind_col].mean() if len(pos_bind_mask) > 0 else 0
        pos_mean_ratio = (pos_mean_bind / pos_mean_expr) if pos_mean_expr > 0 else 0
        
        clean_name = os.path.basename(f)
        stats_list.append({
            "Filename": clean_name,
            "Double+ %": pct_double,
            "Expr MFI": mfi_expr,
            "Bind MFI": mfi_bind,
            "Bind Fold Change": fc_bind,
            "Expr Fold Change": fc_expr,
            "Pos Med Ratio": pos_med_ratio,
            "Pos Mean Ratio": pos_mean_ratio,
            "Bind Med (Expr+)": pos_med_bind,
            "Expr Med (Bind+)": pos_med_expr,
            "Bind Mean (Expr+)": pos_mean_bind,
            "Expr Mean (Bind+)": pos_mean_expr
        })
        
        sub = df_sing.sample(n=min(len(df_sing), 1000))
        for _, row in sub.iterrows():
            ridge_data.append({
                "Sample": clean_name,
                "LogBinding": np.log10(max(row[bind_col], 1)),
                "LogExpression": np.log10(max(row[expr_col], 1))
            })
            
    df_stats = pd.DataFrame(stats_list)
    df_ridge = pd.DataFrame(ridge_data)
    
    if _ctrls["has_pc"] and not df_stats.empty:
        pos_mask = df_stats['Filename'].apply(lambda x: any(p in x for p in ["Positive Control", "TIMP 3"]))
        if pos_mask.any():
            p_med_rat_ref = df_stats.loc[pos_mask, "Pos Med Ratio"].median()
            p_mean_rat_ref = df_stats.loc[pos_mask, "Pos Mean Ratio"].median()
            p_bind_med_ref = df_stats.loc[pos_mask, "Bind Med (Expr+)"].median()
            p_expr_med_ref = df_stats.loc[pos_mask, "Expr Med (Bind+)"].median()
            p_bind_mean_ref = df_stats.loc[pos_mask, "Bind Mean (Expr+)"].median()
            p_expr_mean_ref = df_stats.loc[pos_mask, "Expr Mean (Bind+)"].median()
        else:
            p_med_rat_ref = df_stats["Pos Med Ratio"].max()
            p_mean_rat_ref = df_stats["Pos Mean Ratio"].max()
            p_bind_med_ref = df_stats["Bind Med (Expr+)"].max()
            p_expr_med_ref = df_stats["Expr Med (Bind+)"].max()
            p_bind_mean_ref = df_stats["Bind Mean (Expr+)"].max()
            p_expr_mean_ref = df_stats["Expr Mean (Bind+)"].max()

        df_stats["Norm Median Ratio"] = df_stats["Pos Med Ratio"] / max(p_med_rat_ref, 1e-9)
        df_stats["Norm Mean Ratio"] = df_stats["Pos Mean Ratio"] / max(p_mean_rat_ref, 1e-9)
        df_stats["Norm Bind Med (Expr+)"] = df_stats["Bind Med (Expr+)"] / max(p_bind_med_ref, 1e-9)
        df_stats["Norm Expr Med (Bind+)"] = df_stats["Expr Med (Bind+)"] / max(p_expr_med_ref, 1e-9)
        df_stats["Norm Bind Mean (Expr+)"] = df_stats["Bind Mean (Expr+)"] / max(p_bind_mean_ref, 1e-9)
        df_stats["Norm Expr Mean (Bind+)"] = df_stats["Expr Mean (Bind+)"] / max(p_expr_mean_ref, 1e-9)
        
    elif not df_stats.empty:
        for m in ["Pos Med Ratio", "Pos Mean Ratio", "Bind Med (Expr+)", "Expr Med (Bind+)", "Bind Mean (Expr+)", "Expr Mean (Bind+)"]:
            col_name = f"Norm {m}"
            if m == "Pos Med Ratio": col_name = "Norm Median Ratio"
            if m == "Pos Mean Ratio": col_name = "Norm Mean Ratio"
            df_stats[col_name] = df_stats[m] / max(df_stats[m].max(), 1e-9)
        
    return df_stats, df_ridge

@st.cache_data
def calculate_density(x, y, bins=100):
    """
    Computes binned density for scatter plots, more consistent with flow software
    and 'analyze_fcs.py' hexbin logic.
    """
    if len(x) == 0:
        return np.array([])
        
    counts, x_edges, y_edges = np.histogram2d(x, y, bins=bins)
    
    # Map each point to its corresponding bin index
    ix = np.searchsorted(x_edges, x) - 1
    iy = np.searchsorted(y_edges, y) - 1
    
    # Clip indices to stay within bin bounds
    ix = np.clip(ix, 0, bins - 1)
    iy = np.clip(iy, 0, bins - 1)
    
    # Retrieve the counts and apply log-scaling for the 'flow' look
    z = counts[ix, iy]
    return np.log10(z + 1)

def render_scatter(df, key_prefix):
    col1, col2 = st.columns([1, 3])
    
    with col1:
        st.markdown(f"**Settings ({key_prefix})**")
        x_axis = st.selectbox("X Axis", df.columns, index=df.columns.get_loc("FSC-A") if "FSC-A" in df.columns else 0, key=f"{key_prefix}_x")
        y_axis = st.selectbox("Y Axis", df.columns, index=df.columns.get_loc("SSC-A") if "SSC-A" in df.columns else 1, key=f"{key_prefix}_y")
        
        log_x = st.checkbox("Log X", value=False, key=f"{key_prefix}_logx")
        log_y = st.checkbox("Log Y", value=False, key=f"{key_prefix}_logy")
        
        subsample_n = st.slider("Subsample", 1000, max(1000, len(df)), min(5000, len(df)), step=1000, key=f"{key_prefix}_sub")
        
        color_mode = st.radio("Color Mode", ["Density", "Single Color", "Channel"], horizontal=True, key=f"{key_prefix}_color")
        if color_mode == "Channel":
            color_ch = st.selectbox("Color Channel", df.columns, key=f"{key_prefix}_c_ch")
        
        show_gate = st.checkbox("Show Gates", value=False, key=f"{key_prefix}_gate")
        if show_gate:
            gate_x = st.number_input(f"Gate X ({x_axis})", value=0.0, key=f"{key_prefix}_gx")
            gate_y = st.number_input(f"Gate Y ({y_axis})", value=0.0, key=f"{key_prefix}_gy")

        with st.expander("Manual Axis Limits"):
            use_manual_limits = st.checkbox("Enable Manual Limits", key=f"{key_prefix}_man_lim")
            l_col1, l_col2 = st.columns(2)
            
            x_min_d, x_max_d = float(df[x_axis].min()), float(df[x_axis].max())
            y_min_d, y_max_d = float(df[y_axis].min()), float(df[y_axis].max())
            
            with l_col1:
                x_min = st.number_input(f"X Min", value=x_min_d, key=f"{key_prefix}_xmin")
                x_max = st.number_input(f"X Max", value=x_max_d, key=f"{key_prefix}_xmax")
            with l_col2:
                y_min = st.number_input(f"Y Min", value=y_min_d, key=f"{key_prefix}_ymin")
                y_max = st.number_input(f"Y Max", value=y_max_d, key=f"{key_prefix}_ymax")
                
    with col2:
        if use_manual_limits:
            stats_df = df[
                (df[x_axis] >= x_min) & (df[x_axis] <= x_max) &
                (df[y_axis] >= y_min) & (df[y_axis] <= y_max)
            ]
        else:
            stats_df = df
            
        if len(stats_df) > subsample_n:
            plot_df = stats_df.sample(n=subsample_n, random_state=42).copy()
        else:
            plot_df = stats_df.copy()
            
        if len(plot_df) == 0:
            st.warning("No events in current range.")
            return

        if log_x:
            plot_df[x_axis] = plot_df[x_axis].apply(lambda v: np.log10(max(v, 1e-1)))
            plot_xlabel = f"Log10 {x_axis}"
            gx_val = np.log10(max(gate_x, 1e-1)) if show_gate else 0
        else:
            plot_xlabel = x_axis
            gx_val = gate_x if show_gate else 0
            
        if log_y:
            plot_df[y_axis] = plot_df[y_axis].apply(lambda v: np.log10(max(v, 1e-1)))
            plot_ylabel = f"Log10 {y_axis}"
            gy_val = np.log10(max(gate_y, 1e-1)) if show_gate else 0
        else:
            plot_ylabel = y_axis
            gy_val = gate_y if show_gate else 0

        if color_mode == "Density":
            with st.spinner("Recalculating relative density..."):
                if len(plot_df) > 5:
                    z = calculate_density(plot_df[x_axis], plot_df[y_axis])
                    plot_df["Density"] = z
                    fig = px.scatter(plot_df, x=x_axis, y=y_axis, color="Density",
                                     color_continuous_scale="Jet", template="plotly_dark",
                                     opacity=0.5, render_mode='webgl')
                else:
                     fig = px.scatter(plot_df, x=x_axis, y=y_axis,
                                      template="plotly_dark", opacity=0.6, render_mode='webgl')
            
            # Draw pentagon gate if axes match FSC/SSC
            if show_gate and 'selected_folder' in globals() and 'ctrls' in globals():
                 # We need to check if the axes currently selected match the ones used for gating
                 # Usually FSC-A and SSC-A. Since we are in render_scatter, we check session state or ctrls
                 p_verts = ctrls.get("pentagon_verts")
                 if p_verts is not None and x_axis == "FSC-A" and y_axis == "SSC-A":
                     path_parts = [f"{p[0]},{p[1]}" for p in p_verts]
                     path_str = "M " + " L ".join(path_parts) + " Z"
                     fig.add_shape(type="path", path=path_str, line=dict(color='cyan', dash='dash', width=3), layer='above')
        elif color_mode == "Channel":
             fig = px.scatter(plot_df, x=x_axis, y=y_axis, color=color_ch,
                              color_continuous_scale="Viridis", template="plotly_dark",
                              opacity=0.6, render_mode='webgl')
        else:
            fig = px.scatter(plot_df, x=x_axis, y=y_axis, 
                             template="plotly_dark", opacity=0.6, render_mode='webgl')
            fig.update_traces(marker=dict(color='#00C0FF'))

        if show_gate:
            total = len(plot_df)
            q1 = len(plot_df[(plot_df[x_axis] < gx_val) & (plot_df[y_axis] >= gy_val)])
            q2 = len(plot_df[(plot_df[x_axis] >= gx_val) & (plot_df[y_axis] >= gy_val)])
            q3 = len(plot_df[(plot_df[x_axis] >= gx_val) & (plot_df[y_axis] < gy_val)])
            q4 = len(plot_df[(plot_df[x_axis] < gx_val) & (plot_df[y_axis] < gy_val)])
            
            fig.add_vline(x=gx_val, line_dash="dash", line_color="white")
            fig.add_hline(y=gy_val, line_dash="dash", line_color="white")
            
            st.code(f"""
            Q1 (UL): {q1/total*100:.1f}% | Q2 (UR): {q2/total*100:.1f}%
            Q4 (LL): {q4/total*100:.1f}% | Q3 (LR): {q3/total*100:.1f}%
            """)
            
        fig.update_layout(
            autosize=True, 
            height=500,
            xaxis_title=plot_xlabel,
            yaxis_title=plot_ylabel,
            margin=dict(l=0, r=0, t=0, b=0)
        )
        st.plotly_chart(fig, width='stretch', key=f"scatter_plot_{key_prefix}")

# ---- SIDEBAR ----
st.sidebar.title("🧬 FCS Viewer")
st.sidebar.markdown("Explore raw flow cytometry data interactively.")

if st.session_state["data_mode"] == "Cloud (R2)" and fs:
    default_path = BUCKET
else:
    default_path = os.path.join(os.path.dirname(__file__), "../Local")
    default_path = os.path.abspath(default_path)

base_path = st.sidebar.text_input("Data Root (Path or Bucket)", value=default_path)

# (Rest of the sidebar and processing logic...)
# (Skipping identical lines to reach the TABS section correctly)
with st.sidebar.expander("🧬 Channel Configuration", expanded=False):
    tag_type = st.selectbox("Protein Tag Type", ["His (FITC/APC)", "FLAG (APC/PE)", "Custom"], index=0)

    # Placeholders for custom channels if tag_type is Custom
    custom_expr_placeholder = st.empty()
    custom_bind_placeholder = st.empty()

if tag_type == "His (FITC/APC)":
    def_expr_ch = "FITC-A"
    def_bind_ch = "APC-A"
elif tag_type == "FLAG (APC/PE)":
    def_expr_ch = "APC-A"
    def_bind_ch = "PE-A"
else:
    def_expr_ch = "FITC-A"
    def_bind_ch = "APC-A"

# --- GLOBAL QC & METRIC SETTINGS ---
with st.sidebar.expander("🛠️ Global QC & Metric Settings", expanded=False):
    pc_keywords = st.text_input("PC Keywords (Comma separated)", value="Positive Control, TIMP 3", help="Keywords to identify which rows are positive controls.", key="sidebar_pc_keywords")
    pc_pattern = "|".join([k.strip().upper() for k in pc_keywords.split(",") if k.strip()])
    pc_thresh = st.slider("QC Threshold: Min PC Double+ %", 0.0, 10.0, 2.0, step=0.1, help="Exclude trials where the identified positive control performed poorly.", key="sidebar_pc_thresh")
    
    st.divider()
    metric_options = {
        "Norm Median Ratio": "Norm Median Ratio",
        "Binding Efficiency": "Binding Efficiency",
        "Intensity-Weighted (IWB)": "Norm Intensity-Weighted Binding Index",
        "Norm Median of Expr+": "Norm Bind Med (Expr+)"
    }

    sel_metric_label = st.radio("Select Analysis Metric", list(metric_options.keys()), horizontal=True, index=0, key="sidebar_metric_radio")
    sel_metric_col = metric_options[sel_metric_label]

# --- MANUAL THRESHOLD SETTINGS ---
# Note: Values will be populated later once file/controls are loaded
thresh_expander = st.sidebar.expander("🎯 Manual Threshold Settings", expanded=False)
with thresh_expander:
    thresh_expr_placeholder = st.empty()
    thresh_bind_placeholder = st.empty()



with st.sidebar.expander("🛠️ Advanced Gating Settings"):
    g_min_fsc = st.number_input("Min FSC-A (Singlets)", value=500000, step=50000)
    g_min_ssc = st.number_input("Min SSC-A (Singlets)", value=20000, step=5000)
    g_upper_pct = st.slider("Upper Filter Percentile", 50, 100, 95)
    g_gate_frac = st.slider("Gate Target Fraction", 0.5, 1.0, 0.9)
    g_nc_pct = st.slider("NC Cutoff Percentile", 90.0, 100.0, 99.5, step=0.1)

is_cloud = fs and base_path.startswith(BUCKET)

if is_cloud or os.path.exists(base_path):
    folders = get_fcs_folders(base_path)
    if folders:
        def format_folder(path):
            try:
                if is_cloud:
                    if path == base_path: return "."
                    # Remove the base_path prefix to show relative path
                    return path[len(base_path):].strip("/")
                rel = os.path.relpath(path, base_path)
                return "." if rel == "." else rel
            except:
                return path
                
        selected_folder = st.sidebar.selectbox("Select Folder Context", folders, format_func=format_folder)
        search_path = selected_folder
        
        files = get_files(selected_folder)
        if files:
            selected_file = st.sidebar.selectbox("Select FCS File", files, format_func=os.path.basename)
        else:
            st.sidebar.warning("No .fcs files found in selected folder.")
            selected_file = None
    else:
        st.sidebar.warning("No folders with .fcs files found in Data Root.")
        selected_file = None
        search_path = base_path
else:
    st.sidebar.error("Data Root not found.")
    selected_file = None

st.sidebar.markdown("---")

if selected_file:
    # 1. Preliminary Load to get columns for Custom selection if needed
    meta, df = load_fcs(selected_file)
    
    if df is not None:
        expr_col = def_expr_ch if def_expr_ch in df.columns else df.columns[0]
        bind_col = def_bind_ch if def_bind_ch in df.columns else df.columns[1]
        
        if tag_type == "Custom":
            with custom_expr_placeholder:
                expr_col = st.selectbox("Expression Channel (Y-Axis)", df.columns, index=list(df.columns).index(expr_col) if expr_col in df.columns else 0)
            with custom_bind_placeholder:
                bind_col = st.selectbox("Binding Channel (X-Axis)", df.columns, index=list(df.columns).index(bind_col) if bind_col in df.columns else 1)

    with st.spinner("Analyzing folder controls..."):
        # CRITICAL: Pass the actual channels to ensure correct NC thresholds
        ctrls = analyze_folder_controls(search_path, 
                                      expr_col=expr_col, 
                                      bind_col=bind_col,
                                      nc_percentile=g_nc_pct, 
                                      min_fsc=g_min_fsc, 
                                      min_ssc=g_min_ssc, 
                                      upper_pct=g_upper_pct, 
                                      gate_fraction=g_gate_frac)
    
    st.sidebar.subheader("Gating Context")
    if ctrls["has_nc"]:
        st.sidebar.success("Global Negative Control applied.")
    else:
        st.sidebar.warning("Local Isolated Gating applied.")
        
    if ctrls["has_pc"]:
        st.sidebar.info("Global Positive Control data available.")
    
    # meta, df = load_fcs(selected_file) # Already loaded above
    
    if df is not None:
        # expr_col = "FITC-A" if "FITC-A" in df.columns else df.columns[0]
        # bind_col = "APC-A" if "APC-A" in df.columns else df.columns[1]
        
        if ctrls["has_nc"]:
            def_thresh_expr = float(ctrls["thresh_expr"])
            def_thresh_bind = float(ctrls["thresh_bind"])
        else:
            def_thresh_expr = float(np.percentile(df[expr_col], g_nc_pct))
            def_thresh_bind = float(np.percentile(df[bind_col], g_nc_pct))

        with thresh_expr_placeholder:
            thresh_expr_val = st.number_input(f"Expression Thresh ({expr_col})", value=def_thresh_expr, key=f"input_thresh_expr_{expr_col}_{g_nc_pct}_{search_path}")
        with thresh_bind_placeholder:
            thresh_bind_val = st.number_input(f"Binding Thresh ({bind_col})", value=def_thresh_bind, key=f"input_thresh_bind_{bind_col}_{g_nc_pct}_{search_path}")

        log_thresh_expr = np.log10(max(1, thresh_expr_val))
        log_thresh_bind = np.log10(max(1, thresh_bind_val))


st.sidebar.divider()
st.sidebar.subheader("📊 Data Selection")
if fs:
    st.sidebar.radio("Data Source", ["Cloud (R2)", "Local"], horizontal=True, key="data_mode")
else:
    st.sidebar.warning("☁️ R2 Credentials not found.")

if fs and st.session_state["data_mode"] == "Cloud (R2)":
    st.sidebar.success(f"☁️ Using Cloudflare R2 Data")
else:
    st.sidebar.info("📂 Using Local Data")

# --- AGGREGATE / SELECTIVITY SOURCE & TRIAL EXCLUSION ---
st.sidebar.divider()
st.sidebar.subheader("♻️ Aggregate & Selectivity")

agg_source_mode = st.sidebar.radio(
    "Data Source",
    ["Load from CSVs", "Recalculate (live)"],
    index=0,
    key="agg_source_mode",
    help=(
        "'Load from CSVs' reads the pre-computed sheets in Aggregate_FCS_Analysis. "
        "'Recalculate (live)' computes summary stats on the fly from the raw FCS "
        "files in the selected data source (Local folder or Cloud R2), using the "
        "*_Renamed trial folders — so trial exclusions apply without editing the "
        "saved sheets, and it works even where no *_Analysis folders exist (cloud)."
    ),
)
recalc_mode = (agg_source_mode == "Recalculate (live)")

# Trial folders live under the currently selected Data Root (Local path or bucket).
_all_trial_dirs = [os.path.basename(f.rstrip("/")) for f in discover_trial_folders(base_path)]
excluded_trials = st.sidebar.multiselect(
    "Exclude Trials (folders)",
    _all_trial_dirs,
    default=[],
    key="excluded_trials",
    help=(
        "Trial folders (e.g. *_Renamed) removed here are dropped from BOTH the "
        "aggregate and selectivity analyses. Use for days/trials you know were "
        "unreliable. Applies to the live recalculation; also honored as a "
        "Target+Date filter when loading from CSVs."
    ),
)
if recalc_mode:
    if not AGG_MODULE_OK:
        st.sidebar.warning(f"Live recalculation unavailable: {_AGG_IMPORT_ERROR}")
    elif not _all_trial_dirs:
        st.sidebar.warning("No trial folders (with .fcs) found under the current Data Root.")
    else:
        st.sidebar.caption(f"♻️ Recalculating from {len(_all_trial_dirs) - len(excluded_trials)} trial folder(s).")
        st.sidebar.caption(
            "Updates automatically (like Channel Config): **Advanced Gating**, "
            "**Channel Config**, and **PC Keywords** re-run the recalculation "
            "(gating/channel changes reload FCS — may take a moment; PC-keyword "
            "changes are quick). **QC Threshold** and **Analysis Metric** apply "
            "instantly with no recalculation."
        )

# --- R2 Migration Debug ---
with st.sidebar.expander("🛠️ Debug R2 Connection"):
    st.write(f"**FS Status**: {fs_status}")
    if fs and BUCKET:
        st.write(f"**Bucket**: `{BUCKET}`")
        try:
            test_ls = fs.ls(BUCKET, detail=False)
            st.write(f"Items found: {len(test_ls)}")
        except Exception as e:
            st.write(f"❌ Test Failed: {e}")
    else:
        st.write("❌ No active R2 connection.")

st.sidebar.divider()

# ---- MAIN CONTENT ----
if selected_file and df is not None:
    st.header(f"📂 {os.path.basename(selected_file)}")
    
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Total Events", f"{len(df):,}")
    
    fsc_col = "FSC-A" if "FSC-A" in df.columns else df.columns[0]
    ssc_col = "SSC-A" if "SSC-A" in df.columns else df.columns[1]
    
    if ctrls["has_nc"]:
        pent_path, pent_verts = ctrls["pentagon_path"], ctrls["pentagon_verts"]
    else:
        with st.spinner("Learning pentagon gate on current file..."):
            pent_path, pent_verts = learn_pentagon_gate(df, fsc_col, ssc_col, fraction=g_gate_frac, min_fsc=g_min_fsc, min_ssc=g_min_ssc, upper_pct=g_upper_pct)
            
    df_sing = apply_polygon_gate(df, pent_path, fsc_col, ssc_col, min_fsc=g_min_fsc, min_ssc=g_min_ssc)
    
    if len(df_sing) > 0 and expr_col in df_sing.columns and bind_col in df_sing.columns:
        df_sing = df_sing[(df_sing[expr_col] > 0) & (df_sing[bind_col] > 0)]
    
    c2.metric("Gated Events", f"{len(df_sing):,}")
    c3.metric("% Gated", f"{(len(df_sing)/len(df))*100:.1f}%" if len(df) > 0 else "0%")
    
    if len(df_sing) > 0:
        expr_pos = df_sing[expr_col] > thresh_expr_val
        bind_pos = df_sing[bind_col] > thresh_bind_val
        double_pos = expr_pos & bind_pos
        
        c4.metric("Double+ %", f"{double_pos.mean()*100:.1f}%")
        
        mfi_expr = df_sing[expr_col].median()
        mfi_bind = df_sing[bind_col].median()
        
        if ctrls["has_pc"] and ctrls["has_nc"]:
            fc_expr = mfi_expr / max(ctrls["neg_mfi_expr"], 1)
            fc_bind = mfi_bind / max(ctrls["neg_mfi_bind"], 1)
            fc_vs_pos = mfi_bind / ctrls["pos_mfi_bind"] if ctrls["pos_mfi_bind"] > 0 else 0
            pct_of_pos = (mfi_bind / ctrls["pos_mfi_bind"])*100 if ctrls["pos_mfi_bind"] > 0 else 0
            
            cc1, cc2, cc3 = st.columns(3)
            cc1.metric("Fold Change vs NC (Bind)", f"{fc_bind:.2f}x")
            cc2.metric("Ratio vs PC MFI", f"{fc_vs_pos:.2f}")
            cc3.metric("% of PC MFI", f"{pct_of_pos:.1f}%")
            st.markdown("---")
            
    # --- GLOBAL DATA LOADING (for tab_agg and tab_selectivity) ---
    if fs and BUCKET:
        global_agg_dir = os.path.join(BUCKET, "Aggregate_FCS_Analysis")
    else:
        global_agg_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../Local/Aggregate_FCS_Analysis"))
    
    g_agg_csv = os.path.join(global_agg_dir, "aggregate_summary.csv")
    df_global_trials = None
    if recalc_mode and AGG_MODULE_OK:
        _prog = st.progress(0.0, text="♻️ Recalculating summary stats from raw FCS...")

        def _prog_cb(i, total, name):
            frac = (i / total) if total else 1.0
            label = f"♻️ Recalculating ({i}/{total}) {name}" if name else "♻️ Aggregating..."
            _prog.progress(min(frac, 1.0), text=label)

        _pc_keyword_list = [k.strip() for k in pc_keywords.split(",") if k.strip()]
        df_global_trials = recalc_global_aggregate(
            base_path, tuple(sorted(excluded_trials)),
            expr_col, bind_col, g_nc_pct, g_min_fsc, g_min_ssc, g_upper_pct, g_gate_frac,
            pos_patterns=_pc_keyword_list, progress_cb=_prog_cb,
        )
        _prog.empty()
        if df_global_trials is not None and df_global_trials.empty:
            df_global_trials = None
            st.warning("Live recalculation produced no data (no valid trial folders with a "
                       "negative control in the selected data source).")
    else:
        if recalc_mode and not AGG_MODULE_OK:
            st.warning(f"Live recalculation unavailable: {_AGG_IMPORT_ERROR}. Loading CSVs instead.")
        if cloud_exists(g_agg_csv):
            df_global_trials = cloud_read_csv(g_agg_csv)
        # In CSV mode, still honor the trial-exclusion selection as a post-filter.
        if df_global_trials is not None and excluded_trials:
            excl_td = set()
            for _name in excluded_trials:
                _m = re.match(r'^([^_]+)_(\d{8})', _name)
                if _m:
                    excl_td.add((_m.group(1).upper(), _m.group(2)))
            if excl_td and {"Target", "Date"}.issubset(df_global_trials.columns):
                df_global_trials = df_global_trials[~df_global_trials.apply(
                    lambda r: (str(r["Target"]).upper(), str(r["Date"])) in excl_td, axis=1
                )].copy()

    # Ensure a manufacturer/source column exists (older CSVs predate it); derive
    # it from the sample name so both modes support the manufacturer breakdown.
    if df_global_trials is not None and "Source" not in df_global_trials.columns:
        if "Raw Name" in df_global_trials.columns and AGG_MODULE_OK:
            df_global_trials["Source"] = df_global_trials["Raw Name"].apply(aggregate_analysis.derive_source)
        else:
            df_global_trials["Source"] = "Unknown"

    # --- GLOBAL DATA FILTERING ---
    df_global_filtered = None
    df_pc_status = pd.DataFrame()
    if df_global_trials is not None:
        df_f = df_global_trials.copy()
        
        # Robust column renaming for Binding Efficiency
        rename_map = {}
        for col in df_f.columns:
            if "Binding Efficiency (DP/" in col:
                rename_map[col] = "Binding Efficiency"
        if rename_map:
            df_f = df_f.rename(columns=rename_map)
            
        # 1. Internal QC Columns (from individual folder analysis)
        if 'Trial Failed' in df_f.columns:
            df_f = df_f[df_f['Trial Failed'] == False]
        if 'Low Expression' in df_f.columns:
            df_f = df_f[df_f['Low Expression'] == False]
        if 'Low Events' in df_f.columns:
            df_f = df_f[df_f['Low Events'] == False]
            
        # 2. Dynamic PC QC (based on sidebar threshold)
        df_f["IsPC"] = df_f["Construct"].str.upper().str.contains(pc_pattern) | \
                       df_f["Raw Name"].str.upper().str.contains(pc_pattern)
        
        trial_groups = df_f.groupby(["Target", "Date"])
        valid_trial_ids = []
        pc_info = []
        for (tgt, dt), group in trial_groups:
            pcs = group[group["IsPC"]]
            if not pcs.empty:
                best_pc_dp = pcs["Double+ %"].max()
                pc_row = pcs[pcs["Double+ %"] == best_pc_dp].iloc[0]
                pc_name = pc_row["Raw Name"]
                if best_pc_dp >= pc_thresh:
                    valid_trial_ids.append(f"{tgt}_{dt}")
            else:
                best_pc_dp = 0
                pc_name = "None Found"
            pc_info.append({"Target": tgt, "Date": dt, "PC Name": pc_name, "PC Double+ %": best_pc_dp})
        
        df_pc_status = pd.DataFrame(pc_info)
        df_f["TrialID"] = df_f["Target"].astype(str) + "_" + df_f["Date"].astype(str)
        df_global_filtered = df_f[df_f["TrialID"].isin(valid_trial_ids)].copy()

    # TABS


    tab_main, tab_pos, tab_folder, tab_agg, tab_selectivity, tab_custom, tab_raw, tab_meta = st.tabs([
        "Main Analysis", "Pos Control Analytics", "Folder Analysis", "Aggregate Analysis", "Selectivity Analysis", "Custom Plots", "Raw Data", "Metadata"
    ])
    
    with tab_main:
        if len(df_sing) > 0:
            st.subheader("Analysis Dashboard")
            col_a, col_b = st.columns(2)
            
            SUB_N = 50000
            plot_df_full = df if len(df) <= SUB_N else df.sample(SUB_N, random_state=42)
            plot_df_sing = df_sing if len(df_sing) <= SUB_N else df_sing.sample(SUB_N, random_state=42)
            
            with col_a:
                z = calculate_density(plot_df_full[fsc_col].values, plot_df_full[ssc_col].values)
                fig1 = go.Figure()
                fig1.add_trace(go.Scattergl(
                    x=plot_df_full[fsc_col], 
                    y=plot_df_full[ssc_col],
                    mode='markers',
                    marker=dict(size=3, color=z, colorscale='Inferno', opacity=0.5),
                    name='Events'
                ))
                if pent_verts is not None:
                    # Create path string for SVG path
                    path_parts = [f"{p[0]},{p[1]}" for p in pent_verts]
                    path_str = "M " + " L ".join(path_parts) + " Z"
                    
                    fig1.add_shape(
                        type="path",
                        path=path_str,
                        line=dict(color='cyan', dash='dash', width=3),
                        layer='above'
                    )
                fig1.update_layout(
                    title=f"Original Gate ({len(df_sing)}/{len(df)})",
                    xaxis_title=fsc_col,
                    yaxis_title=ssc_col,
                    template='plotly_dark',
                    height=450, margin=dict(l=0,r=0,b=0)
                )
                st.plotly_chart(fig1, width='stretch', key="main_gate_fig")

            with col_b:
                t_expr = np.log10(np.clip(plot_df_sing[expr_col], 1, None))
                fig2 = px.histogram(
                    x=t_expr, 
                    nbins=100, 
                    histnorm='density',
                    color_discrete_sequence=['purple'],
                    template='plotly_dark'
                )
                fig2.add_vline(x=log_thresh_expr, line_dash="dash", line_color="red")
                fig2.update_layout(
                    title=f"Expression Distribution ({expr_col.replace('-A', '')})",

                    xaxis_title=f"Log10 {expr_col}",
                    yaxis_title="Density",
                    height=450, margin=dict(l=0,r=0,b=0),
                    showlegend=False
                )
                st.plotly_chart(fig2, width='stretch', key="main_expr_hist")
                
            col_c, col_d = st.columns(2)
            
            with col_c:
                t_expr_all = np.log10(np.clip(plot_df_sing[expr_col], 1, None)).values
                t_bind_all = np.log10(np.clip(plot_df_sing[bind_col], 1, None)).values
                
                z2 = calculate_density(t_bind_all, t_expr_all)
                fig3 = go.Figure()
                fig3.add_trace(go.Scattergl(
                    x=t_bind_all, 
                    y=t_expr_all,
                    mode='markers',
                    marker=dict(size=4, color=z2, colorscale='Jet', opacity=0.6),
                    name='Events'
                ))
                fig3.add_vline(x=log_thresh_bind, line_dash="dash", line_color="white")
                fig3.add_hline(y=log_thresh_expr, line_dash="dash", line_color="white")
                
                pct_ll = (~expr_pos & ~bind_pos).mean() * 100
                pct_lr = (~expr_pos & bind_pos).mean() * 100
                pct_ul = (expr_pos & ~bind_pos).mean() * 100
                pct_ur = (expr_pos & bind_pos).mean() * 100
                
                title_str = (
                    f"Q1 (UL): {pct_ul:.1f}% | Q2 (UR): {pct_ur:.1f}%<br>"
                    f"Q3 (LL): {pct_ll:.1f}% | Q4 (LR): {pct_lr:.1f}%"
                )
                
                fig3.update_layout(
                    title=title_str,
                    xaxis_title=f"Log10 {bind_col} (Binding)",
                    yaxis_title=f"Log10 {expr_col} (Expression)",
                    template='plotly_dark',
                    height=450, margin=dict(l=0,r=0,b=0),
                    showlegend=False
                )
                st.plotly_chart(fig3, width='stretch', key="main_quad_scatter")

            with col_d:
                t_bind = np.log10(np.clip(plot_df_sing[bind_col], 1, None))
                fig4 = px.histogram(
                    x=t_bind, 
                    nbins=100, 
                    histnorm='density',
                    color_discrete_sequence=['green'],
                    template='plotly_dark'
                )
                fig4.add_vline(x=log_thresh_bind, line_dash="dash", line_color="red")
                fig4.update_layout(
                    title=f"Binding Distribution ({bind_col.replace('-A', '')})",

                    xaxis_title=f"Log10 {bind_col}",
                    yaxis_title="Density",
                    height=450, margin=dict(l=0,r=0,b=0),
                    showlegend=False
                )
                st.plotly_chart(fig4, width='stretch', key="main_bind_hist")
                
            if ctrls.get("neg_sample") is not None:
                st.markdown("---")
                st.subheader("Comparison vs Controls")
                
                comp_a, comp_b = st.columns(2)
                
                # Combine data for Expression Plotly Histogram
                val_expr = []
                cond_expr = []
                
                neg_e_vals = np.log10(np.clip(ctrls["neg_sample"][expr_col], 1, None)).values
                val_expr.extend(neg_e_vals)
                cond_expr.extend(["Neg Ctrl"] * len(neg_e_vals))
                
                if ctrls["has_pc"] and ctrls.get("pos_sample") is not None:
                    pos_e_vals = np.log10(np.clip(ctrls["pos_sample"][expr_col], 1, None)).values
                    val_expr.extend(pos_e_vals)
                    cond_expr.extend(["Pos Ctrl"] * len(pos_e_vals))
                
                sam_e_vals = t_expr.values
                val_expr.extend(sam_e_vals)
                cond_expr.extend(["Sample"] * len(sam_e_vals))
                
                df_e_comp = pd.DataFrame({'Value': val_expr, 'Condition': cond_expr})
                fig_c2 = px.histogram(df_e_comp, x='Value', color='Condition', barmode='overlay', histnorm='density', template='plotly_dark')
                fig_c2.add_vline(x=log_thresh_expr, line_dash="dash", line_color="red")
                fig_c2.update_layout(title="Expression Comparison", xaxis_title=f"Log10 {expr_col}", yaxis_title="Density", height=400, margin=dict(l=0,r=0,b=0))
                
                with comp_a:
                    st.plotly_chart(fig_c2, width='stretch', key="main_expr_comp_hist")

                # Combine data for Binding Plotly Histogram
                val_bind = []
                cond_bind = []
                
                neg_b_vals = np.log10(np.clip(ctrls["neg_sample"][bind_col], 1, None)).values
                val_bind.extend(neg_b_vals)
                cond_bind.extend(["Neg Ctrl"] * len(neg_b_vals))
                
                if ctrls["has_pc"] and ctrls.get("pos_sample") is not None:
                    pos_b_vals = np.log10(np.clip(ctrls["pos_sample"][bind_col], 1, None)).values
                    val_bind.extend(pos_b_vals)
                    cond_bind.extend(["Pos Ctrl"] * len(pos_b_vals))
                
                sam_b_vals = t_bind.values
                val_bind.extend(sam_b_vals)
                cond_bind.extend(["Sample"] * len(sam_b_vals))
                
                df_b_comp = pd.DataFrame({'Value': val_bind, 'Condition': cond_bind})
                fig_c1 = px.histogram(df_b_comp, x='Value', color='Condition', barmode='overlay', histnorm='density', template='plotly_dark')
                fig_c1.add_vline(x=log_thresh_bind, line_dash="dash", line_color="red")
                fig_c1.update_layout(title="Binding Comparison", xaxis_title=f"Log10 {bind_col}", yaxis_title="Density", height=400, margin=dict(l=0,r=0,b=0))
                
                with comp_b:
                    st.plotly_chart(fig_c1, width='stretch', key="main_bind_comp_hist")
                    
        else:
            st.warning("No gated events found.")

    with tab_pos:
        st.subheader("Positive Event Distributions")
        
        expr_ch_lbl = expr_col.replace("-A", "")
        bind_ch_lbl = bind_col.replace("-A", "")
        
        pos_expr_evts = df_sing[df_sing[expr_col] > thresh_expr_val]
        pos_bind_evts = df_sing[df_sing[bind_col] > thresh_bind_val]

        
        # Function to build the styled histogram with overlay and markers
        def build_pos_hist(series, nc_series, thresh, title, x_label, color):
            fig = go.Figure()
            # Calculate metrics
            median_val = np.log10(np.clip(series.median(), 1, None))
            mean_val = np.log10(np.clip(series.mean(), 1, None))
            log_thresh = np.log10(np.clip(thresh, 1, None))
            
            # 1. Negative Control (NC) Overlay Trace
            if nc_series is not None:
                nc_vals = np.log10(np.clip(nc_series, 1, None))
                fig.add_trace(go.Histogram(
                    x=nc_vals, nbinsx=100, histnorm='density',
                    name='Neg Ctrl', marker_color='rgba(150, 150, 150, 0.4)',
                    showlegend=True
                ))
            
            # 2. Sample Trace
            fig.add_trace(go.Histogram(
                x=np.log10(np.clip(series, 1, None)), nbinsx=100, histnorm='density',
                name='Sample', marker_color=color,
                opacity=0.8, showlegend=True
            ))
            
            # 3. Vertical Markers
            fig.add_vline(x=log_thresh, line_dash="dash", line_color="red", annotation_text="Thresh", annotation_position="top left")
            fig.add_vline(x=median_val, line_dash="solid", line_color="white", annotation_text="Med", annotation_position="top left")
            fig.add_vline(x=mean_val, line_dash="dot", line_color="white", annotation_text="Mean", annotation_position="bottom right")
            
            fig.update_layout(
                title=title, xaxis_title=x_label, yaxis_title="Density",
                template='plotly_dark', barmode='overlay',
                height=400, margin=dict(l=0,r=0,b=20,t=40),
                legend=dict(yanchor="top", y=0.99, xanchor="right", x=0.99, bgcolor="rgba(0,0,0,0)")
            )
            return fig

        pc_c1, pc_c2 = st.columns(2)
        nc_sample = ctrls.get("neg_sample")
        
        if len(pos_expr_evts) > 0:
            fig_pe = build_pos_hist(pos_expr_evts[expr_col], nc_sample[expr_col] if nc_sample is not None else None, 
                                   thresh_expr_val, "Positive Expression Dist", f"Log10 {expr_col} (Expr+ ONLY)", "purple")
            pc_c1.plotly_chart(fig_pe, width='stretch', key="pos_expr_dist_hist")
            pe_med = pos_expr_evts[expr_col].median()
            pe_mean = pos_expr_evts[expr_col].mean()
            pc_c1.latex(rf"\text{{Med: }} {pe_med:.1f} \quad \text{{Mean: }} {pe_mean:.1f}")
        else:
            pc_c1.info("No Expr+ events.")
            
        if len(pos_bind_evts) > 0:
            fig_pb = build_pos_hist(pos_bind_evts[bind_col], nc_sample[bind_col] if nc_sample is not None else None, 
                                   thresh_bind_val, "Positive Binding Dist", f"Log10 {bind_col} (Bind+ ONLY)", "green")
            pc_c2.plotly_chart(fig_pb, width='stretch', key="pos_bind_dist_hist")
            pb_med = pos_bind_evts[bind_col].median()
            pb_mean = pos_bind_evts[bind_col].mean()
            pc_c2.latex(rf"\text{{Med: }} {pb_med:.1f} \quad \text{{Mean: }} {pb_mean:.1f}")
        else:
            pc_c2.info("No Bind+ events.")

        st.subheader("Filtered Histograms (Cross-Gated)")
        f1, f2 = st.columns(2)
        
        # Swapped Order: Expression on Left, Binding on Right
        if len(pos_bind_evts) > 0:
            fig_fe = build_pos_hist(pos_bind_evts[expr_col], nc_sample[expr_col] if nc_sample is not None else None, 
                                   thresh_expr_val, f"Filtered Expression ({expr_ch_lbl} for {bind_ch_lbl}+)", f"Log10 {expr_col} ({bind_ch_lbl}+ ONLY)", "orange")

            f1.plotly_chart(fig_fe, width='stretch', key="filt_expr_hist")
            exp_med = pos_bind_evts[expr_col].median()
            exp_mean = pos_bind_evts[expr_col].mean()
            f1.latex(rf"\text{{Med: }} {exp_med:.1f} \quad \text{{Mean: }} {exp_mean:.1f}")
        else:
            f1.info("No APC+ events.")

        if len(pos_expr_evts) > 0:
            fig_fb = build_pos_hist(pos_expr_evts[bind_col], nc_sample[bind_col] if nc_sample is not None else None, 
                                   thresh_bind_val, f"Filtered Binding ({bind_ch_lbl} for {expr_ch_lbl}+)", f"Log10 {bind_col} ({expr_ch_lbl}+ ONLY)", "cyan")

            f2.plotly_chart(fig_fb, width='stretch', key="filt_bind_hist")
            bin_med = pos_expr_evts[bind_col].median()
            bin_mean = pos_expr_evts[bind_col].mean()
            f2.latex(rf"\text{{Med: }} {bin_med:.1f} \quad \text{{Mean: }} {bin_mean:.1f}")
        else:
            f2.info(f"No {expr_ch_lbl}+ events.")
            
        st.subheader("Positive Metrics for this file")

        pos_med_expr_v = pos_expr_evts[expr_col].median() if len(pos_expr_evts) > 0 else 0
        pos_med_bind_v = pos_bind_evts[bind_col].median() if len(pos_bind_evts) > 0 else 0
        raw_pos_med_ratio = (pos_med_bind_v / pos_med_expr_v) if pos_med_expr_v > 0 else 0
        
        cmp1, cmp2 = st.columns(2)
        cmp1.metric(f"Raw Pos Median Ratio ({bind_ch_lbl} / {expr_ch_lbl})", f"{raw_pos_med_ratio:.3f}")

        
        if ctrls["has_pc"] and ctrls.get("pos_sample") is not None:
            pc_expr_p = ctrls["pos_sample"][expr_col] > thresh_expr_val
            pc_bind_p = ctrls["pos_sample"][bind_col] > thresh_bind_val
            pc_m_e = ctrls["pos_sample"][pc_expr_p][expr_col].median() if pc_expr_p.sum() > 0 else 1
            pc_m_b = ctrls["pos_sample"][pc_bind_p][bind_col].median() if pc_bind_p.sum() > 0 else 1
            pc_ratio = pc_m_b / pc_m_e if pc_m_e > 0 else 1
            
            norm_ratio = raw_pos_med_ratio / max(pc_ratio, 1e-9)
            cmp2.metric("Normalized Pos Median Ratio vs Pos Ctrl", f"{norm_ratio:.3f}")
            
            df_comp_ratio = pd.DataFrame({
                "Sample": [os.path.basename(selected_file), "Positive Control"],
                "Ratio": [raw_pos_med_ratio, pc_ratio]
            })
            fig_ratio = px.bar(df_comp_ratio, x="Sample", y="Ratio", color="Sample", template="plotly_dark", title="Pos Median Ratio Comparison")
            st.plotly_chart(fig_ratio, width='stretch', key="pos_ratio_comp_bar")
        else:
            st.info("No Positive Control in directory to compute normalized ratios.")

    with tab_custom:
        enable_comparison = st.checkbox("Enable Comparison Mode (Side-by-Side)", value=False)
        if enable_comparison:
            st.markdown("### Comparison View")
            c_a, c_b = st.columns(2)
            with c_a:
                st.caption("Plot A")
                render_scatter(df, "A")
            with c_b:
                st.caption("Plot B")
                render_scatter(df, "B")
        else:
            render_scatter(df, "Main")

    with tab_raw:
        st.dataframe(df, width='stretch', height=600)
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button("Download CSV", csv, "fcs_data.csv", "text/csv")
        
    with tab_meta:
        st.subheader("File Metadata")
        meta_items = []
        for k, v in meta.items():
            if not k.startswith("$P"):
                meta_items.append({"Parameter": k, "Value": str(v)})
        st.dataframe(pd.DataFrame(meta_items), width='stretch', hide_index=True)
        
        st.subheader("Channel Parameters")
        params = []
        i = 1
        while f"$P{i}N" in meta:
            params.append({
                "Index": i, 
                "Name": meta.get(f"$P{i}N"), 
                "Description": meta.get(f"$P{i}S", "-"), 
                "Range": meta.get(f"$P{i}R", "-")
            })
            i += 1
        if params:
            st.dataframe(pd.DataFrame(params), width='stretch', hide_index=True)

    with tab_folder:
        if search_path and cloud_exists(search_path):
            st.subheader(f"Folder Aggregates: {os.path.basename(search_path)}")
            agg_csv = os.path.join(search_path, "aggregate_summary.csv")
            cross_csv = os.path.join(search_path, "cross_target_summary.csv")
            sum_csv = os.path.join(search_path, "summary_stats.csv")
            
            has_agg = False
            if cloud_exists(cross_csv):
                st.markdown("**Cross-Target Summary**")
                st.dataframe(cloud_read_csv(cross_csv), width='stretch')
                has_agg = True
            elif cloud_exists(agg_csv):
                st.markdown("**Aggregate Summary**")
                st.dataframe(cloud_read_csv(agg_csv), width='stretch')
                has_agg = True
            elif cloud_exists(sum_csv):
                st.markdown("**Summary Stats**")
                st.dataframe(cloud_read_csv(sum_csv), width='stretch')
                has_agg = True
                
            if not has_agg:
                st.info("No static aggregate CSVs found in this folder. Generating dynamically...")
                with st.spinner("Generating Folder Aggregate Stats dynamically..."):
                    df_stats, df_ridge = get_folder_aggregate_stats(search_path, ctrls, fsc_col, ssc_col, expr_col, bind_col, thresh_expr_val, thresh_bind_val, min_fsc=g_min_fsc, min_ssc=g_min_ssc, upper_pct=g_upper_pct)
                    
                if not df_stats.empty:
                    hide_nc_folder = st.checkbox("Hide Negative Controls", value=True, key="hide_nc_fld")
                    
                    df_stats_plot = df_stats.copy()
                    if hide_nc_folder:
                        df_stats_plot = df_stats_plot[~df_stats_plot['Filename'].str.upper().str.contains("NC|NEGATIVE CONTROL")]
                        
                    df_stats_sorted = df_stats_plot.sort_values("Filename")
                    st.dataframe(df_stats_plot, width='stretch')
                    
                    row1_c1, row1_c2 = st.columns(2)
                    
                    # Robust color column selection
                    agg_color_col = "Norm Median Ratio" if "Norm Median Ratio" in df_stats_sorted.columns else ("Norm Pos Med Ratio" if "Norm Pos Med Ratio" in df_stats_sorted.columns else None)
                    
                    fig_dp = px.bar(df_stats_sorted, x="Filename", y="Double+ %", color=agg_color_col, template="plotly_dark", title="Double Positive %")
                    row1_c1.plotly_chart(fig_dp, width='stretch', key="agg_dp_bar")
                    
                    fig_mfi = px.scatter(df_stats, x="Bind MFI", y="Expr MFI", color="Filename", template="plotly_dark", title=f"MFI Scatter ({bind_ch_lbl} vs {expr_ch_lbl})")

                    row1_c2.plotly_chart(fig_mfi, width='stretch', key="agg_mfi_scatter")
                    
                    row2_c1, row2_c2 = st.columns(2)
                    
                    fig_fc = px.bar(df_stats_sorted, x="Filename", y="Bind Fold Change", color="Double+ %", template="plotly_dark", title="Binding Fold Change (vs NC)")
                    row2_c1.plotly_chart(fig_fc, width='stretch', key="agg_bind_fc_bar")

                    fig_efc = px.bar(df_stats_sorted, x="Filename", y="Expr Fold Change", color="Double+ %", template="plotly_dark", title="Expression Fold Change (vs NC)")
                    row2_c2.plotly_chart(fig_efc, width='stretch', key="agg_expr_fc_bar")
                    
                    df_ridge_sorted = df_ridge.sort_values("Sample", ascending=False)
                    if hide_nc_folder:
                        df_ridge_sorted = df_ridge_sorted[~df_ridge_sorted['Sample'].str.upper().str.contains("NC|NEGATIVE CONTROL")]

                    fig_ridge = px.violin(df_ridge_sorted, x="LogBinding", y="Sample", color="Sample", orientation="h", template="plotly_dark", title=f"Ridgeline: Log({bind_ch_lbl})", points=False)

                    fig_ridge.update_traces(side='positive', width=2.5)
                    fig_ridge.update_layout(xaxis_showgrid=False, xaxis_zeroline=False, showlegend=False)
                    st.plotly_chart(fig_ridge, width='stretch', key="agg_ridge_bind")
                    
                    fig_ridge_ex = px.violin(df_ridge_sorted, x="LogExpression", y="Sample", color="Sample", orientation="h", template="plotly_dark", title=f"Ridgeline: Log({expr_ch_lbl})", points=False)

                    fig_ridge_ex.update_traces(side='positive', width=2.5)
                    fig_ridge_ex.update_layout(xaxis_showgrid=False, xaxis_zeroline=False, showlegend=False)
                    st.plotly_chart(fig_ridge_ex, width='stretch', key="agg_ridge_expr")

                    if ctrls["has_pc"]:
                        st.markdown("---")
                        st.subheader("📁 Positive Control Aggregate Analytics")
                        
                        view_type = st.radio("Display Analytics Mode:", ["Normalized (vs Positive Control)", "Raw (MFI / Ratios)"], horizontal=True, key="an_mode_final")
                        is_norm = "Normalized" in view_type
                        pref = "Norm " if is_norm else ""

                        plots = [
                            ("Norm Median Ratio" if is_norm else "Pos Med Ratio", "Pos Median Ratio (Bind/Expr)"),
                            ("Norm Mean Ratio" if is_norm else "Pos Mean Ratio", "Pos Mean Ratio (Bind/Expr)"),
                            (f"{pref}Expr Med (Bind+)", f"Median {expr_ch_lbl} for Positive {bind_ch_lbl}"),
                            (f"{pref}Bind Med (Expr+)", f"Median {bind_ch_lbl} for Positive {expr_ch_lbl}"),

                            (f"{pref}Expr Mean (Bind+)", f"Mean {expr_ch_lbl} for Positive {bind_ch_lbl}"),
                            (f"{pref}Bind Mean (Expr+)", f"Mean {bind_ch_lbl} for Positive {expr_ch_lbl}")

                        ]
                        
                        for i in range(0, len(plots), 2):
                            rrc1, rrc2 = st.columns(2)
                            p1_col, p1_title = plots[i]
                            fig_p1 = px.bar(df_stats_sorted, x="Filename", y=p1_col, color="Double+ %", template="plotly_dark", title=f"{pref}{p1_title}" if is_norm else p1_title)
                            if is_norm: fig_p1.add_hline(y=1, line_dash="dash", line_color="white", annotation_text="Pos Ctrl")
                            rrc1.plotly_chart(fig_p1, width='stretch', key=f"agg_pos_ctrl_{p1_col}")
                            
                            if i+1 < len(plots):
                                p2_col, p2_title = plots[i+1]
                                fig_p2 = px.bar(df_stats_sorted, x="Filename", y=p2_col, color="Double+ %", template="plotly_dark", title=f"{pref}{p2_title}" if is_norm else p2_title)
                                if is_norm: fig_p2.add_hline(y=1, line_dash="dash", line_color="white", annotation_text="Pos Ctrl")
                                rrc2.plotly_chart(fig_p2, width='stretch', key=f"agg_pos_ctrl_{p2_col}")
                else:
                    st.warning("Could not compute aggregate stats on the fly.")

    with tab_agg:
        if df_global_filtered is not None:
            st.subheader("📊 Global Aggregate Analysis")
            if recalc_mode:
                _src_lbl = "Cloud R2" if (fs and BUCKET and str(base_path).startswith(BUCKET)) else "Local"
                _excl_note = f" ({len(excluded_trials)} trial(s) excluded)" if excluded_trials else ""
                st.caption(f"♻️ Live recalculation from raw FCS in {_src_lbl}{_excl_note} — computed on the fly, not read from saved CSVs.")

            g_cross_csv = os.path.join(global_agg_dir, "cross_target_summary.csv")
            
            found_g = True
            available_targets = sorted(df_global_filtered["Target"].unique())
            target_choice = st.selectbox("🎯 Filter by Target", ["All Targets"] + list(available_targets), key="global_target_sel")
            
            df_agg_f = df_global_filtered.copy()
            if target_choice != "All Targets":
                df_agg_f = df_agg_f[df_agg_f["Target"] == target_choice]
            
            # Separate Constructs from PCs for plotting
            df_constructs = df_agg_f[~df_agg_f["IsPC"]].copy()

            
            # --- NEW: Aggregated Summary View (Mean/SD) ---
            st.markdown("### 🏆 Combined Variant Performance (Mean/SD)")
            if not df_constructs.empty:
                # Ensure alphabetical sorting here
                df_constructs = df_constructs.sort_values("Construct")
                
                # Update grouping to use selected metric
                if sel_metric_col not in df_constructs.columns:
                    st.error(f"❌ Metric **{sel_metric_col}** not found in current dataset.")
                    if st.session_state.get("data_mode") == "Cloud (R2)":
                        st.info("You are in **Cloud (R2)** mode. This metric might require uploading the latest `aggregate_summary.csv` to R2. Try switching to **Local** mode in the sidebar if you just ran a local analysis.")
                    else:
                        st.info("You are in **Local** mode. Please ensure you have run `aggregate_analysis.py` to regenerate the summary data with the required metrics.")
                    st.stop()

                df_grouped = df_constructs.groupby("Construct").agg({
                    sel_metric_col: ["mean", "std", "count"],
                    "Double+ %": "mean"
                }).reset_index()
                df_grouped.columns = ["Construct", "Mean Val", "SD Val", "n", "Avg Double+ %"]
                
                # Calculate SEM and 95% CI
                df_grouped['SEM'] = df_grouped.apply(lambda row: row['SD Val'] / (row['n']**0.5) if row['n'] > 1 else 0, axis=1)
                df_grouped['95% CI'] = df_grouped['SEM'] * 1.96
                
                # User requested alphabetical order
                df_grouped = df_grouped.sort_values("Construct")
                
                fig_sum = px.bar(
                    df_grouped, x="Construct", y="Mean Val", error_y="95% CI",
                    color="Avg Double+ %", color_continuous_scale="Viridis",
                    template="plotly_dark", title=f"{sel_metric_label} by Construct (Error Bars: 95% CI)",
                    text_auto='.2f',
                    hover_data={"n": True, "95% CI": ":.2f"}
                )


                fig_sum.update_layout(height=550)
                st.plotly_chart(fig_sum, width='stretch', key="agg_mean_sd_bar")
                
                # DYNAMIC ANOVA: Across Constructs for this Target
                p_val_agg = run_anova_p(df_constructs, 'Construct', sel_metric_col)
                if not np.isnan(p_val_agg):
                    msg = f"📊 **Statistical Analysis (ANOVA)**: Comparing all constructs for **{target_choice}** using **{sel_metric_label}**.  \n**p-value**: `{p_val_agg:.4e}` ({'Significant' if p_val_agg < 0.05 else 'Not Significant'})"
                    if p_val_agg < 0.05:
                        st.success(msg)
                    else:
                        st.info(msg)
                    
                    if p_val_agg < 0.05:

                        sig_pairs = run_tukey_summary(df_constructs, 'Construct', sel_metric_col)
                        if sig_pairs is not None and not sig_pairs.empty:
                            with st.expander("🔍 View Significant Pairwise Differences (Tukey HSD)"):
                                st.caption(TUKEY_NOTE)
                                st.dataframe(sig_pairs[['group1', 'group2', 'meandiff', 'p-adj']], hide_index=True, width='stretch')
                else:
                    st.caption("ANOVA could not be calculated (requires multiple trials per construct).")


            else:
                st.info("No constructs match current filters or QC criteria.")

            with st.expander("📋 Review Identified Positive Controls"):
                st.dataframe(df_pc_status, hide_index=True, width='stretch')

            # --- NEW: COMPREHENSIVE STATISTICAL SUMMARY (ALL METRICS) ---
            st.divider()
            with st.expander("📊 Comprehensive Statistical Summary (All Metrics)", expanded=False):
                st.markdown(f"ANOVA-gated **Tukey HSD** pairwise comparisons for **{target_choice}** across all metrics.")
                st.caption(TUKEY_NOTE)
                all_sig_findings = []
                
                for m_label, m_col in metric_options.items():
                    p_val = run_anova_p(df_constructs, 'Construct', m_col)
                    if not np.isnan(p_val) and p_val < 0.05:
                        sig_pairs = run_tukey_summary(df_constructs, 'Construct', m_col)
                        if sig_pairs is not None and not sig_pairs.empty:
                            for _, row in sig_pairs.iterrows():
                                all_sig_findings.append({
                                    "Metric": m_label,
                                    "Pair": f"{row['group1']} vs {row['group2']}",
                                    "Diff": f"{row['meandiff']:.4f}",
                                    "p-adj": f"{row['p-adj']:.4f}"
                                })
                
                if all_sig_findings:
                    st.dataframe(pd.DataFrame(all_sig_findings), hide_index=True, width='stretch')
                else:
                    st.info("No significant differences found across any metrics for this target at the current QC threshold.")


            # Visualizations
            st.markdown(f"### 📅 Trial-by-Trial Analytics: {target_choice}")
            
            # Sort for alphabetical order in trial plots
            df_agg_f = df_agg_f.sort_values("Construct")
            
            v_col1, v_col2 = st.columns(2)
            
            with v_col1:
                # Adaptive Bar Chart based on selected metric
                fig_agg_ratio = px.bar(
                    df_agg_f, x="Construct", y=sel_metric_col, color="Date",
                    template="plotly_dark", title=f"{sel_metric_label} by Construct",
                    text=sel_metric_col
                )
                fig_agg_ratio.update_traces(texttemplate='%{text:.2f}', textposition='outside')
                fig_agg_ratio.add_hline(y=1, line_dash="dash", line_color="white", annotation_text="Pos Ctrl")
                fig_agg_ratio.update_layout(height=500)
                st.plotly_chart(fig_agg_ratio, width='stretch', key="agg_ratio_global_bar")
                
            with v_col2:
                # Double+ % Bar Chart
                fig_agg_dp = px.bar(
                    df_agg_f, x="Construct", y="Double+ %", color="Date",

                    template="plotly_dark", title="Double Positive Percentage (%)",
                    text="Double+ %"
                )
                fig_agg_dp.update_traces(texttemplate='%{text:.1f}%', textposition='outside')
                fig_agg_dp.update_layout(height=500)
                st.plotly_chart(fig_agg_dp, width='stretch', key="agg_dp_global_bar")
            
            st.markdown("**Global Aggregate Data Table**")
            st.dataframe(df_agg_f, width='stretch', hide_index=True)
            st.markdown("---")


            if recalc_mode:
                st.caption("ℹ️ The static Cross-Target Summary sheet is not shown in live "
                           "recalculation mode (it would not reflect your exclusions).")
            elif cloud_exists(g_cross_csv):
                found_g = True
                st.markdown("**Global Cross-Target Summary (Experimental Mapping)**")
                st.dataframe(cloud_read_csv(g_cross_csv), width='stretch', hide_index=True)

            if not recalc_mode and not found_g:
                st.info(f"Global directory found at '{global_agg_dir}', but no summary CSVs were detected.")
        else:
            st.info(f"Global aggregate directory not found at: {global_agg_dir}")

    with tab_selectivity:
        st.subheader("🎯 Selectivity Analysis")
        if st.button("🔄 Refresh Selectivity Data", key="refresh_selectivity"):
            st.rerun()

        # Metric options: raw column + normalized (vs Pos Ctrl) column where available.
        SEL_METRICS = {
            "Median Binding Ratio": {"raw": "Pos Med Ratio", "norm": "Norm Median Ratio"},
            "Median of Expr+": {"raw": "Bind Med (Expr+)", "norm": "Norm Bind Med (Expr+)"},
            "Binding Efficiency": {"raw": "Binding Efficiency", "norm": None},
            "IWB Index": {"raw": "Intensity-Weighted Binding Index", "norm": "Norm Intensity-Weighted Binding Index"},
        }

        if df_global_filtered is None or df_global_filtered.empty:
            if recalc_mode:
                st.info("No recalculated trials available for selectivity. Check that trial "
                        "folders exist in the selected data source and are not all excluded.")
            else:
                st.info("No aggregate data available. Ensure aggregate_summary.csv exists "
                        "(Load from CSVs) or switch the source to Recalculate (live).")
        else:
            _src_lbl = "♻️ Live recalculation" if recalc_mode else "📄 Loaded from CSV"
            st.caption(f"{_src_lbl} — selectivity computed from the current aggregate data "
                       "(honors data source, trial exclusions, gating & QC settings).")

            csel1, csel2 = st.columns([2, 2])
            with csel1:
                selected_label = st.selectbox("Metric", list(SEL_METRICS.keys()),
                                              key="selectivity_metric_sel")
            with csel2:
                view_mode = st.radio("Values", ["Raw", "Normalized (vs Pos Ctrl)"],
                                     horizontal=True, key="selectivity_norm_toggle")

            _m_info = SEL_METRICS[selected_label]
            _use_norm = view_mode.startswith("Normalized")
            metric_col = _m_info["norm"] if (_use_norm and _m_info["norm"]) else _m_info["raw"]
            _is_norm = _use_norm and bool(_m_info["norm"])
            if _use_norm and not _m_info["norm"]:
                st.caption(f"ℹ️ '{selected_label}' has no normalized form (it is already a "
                           "ratio); showing raw values.")

            if metric_col not in df_global_filtered.columns:
                st.error(f"Metric column '{metric_col}' not found in the current dataset.")
            else:
                y_label = f"{selected_label}" + (" (Norm)" if _is_norm else "")

                # Ensure a Source (manufacturer) column exists.
                df_src = df_global_filtered.copy()
                if "Source" not in df_src.columns:
                    if "Raw Name" in df_src.columns and AGG_MODULE_OK:
                        df_src["Source"] = df_src["Raw Name"].apply(aggregate_analysis.derive_source)
                    else:
                        df_src["Source"] = "Unknown"

                # Restrict to constructs tested against more than one target.
                _tc = df_src.groupby("Construct")["Target"].nunique()
                _multi = _tc[_tc > 1].index.tolist()
                df_multi = df_src[df_src["Construct"].isin(_multi)].copy()

                def _sel_summary(df_in, group_cols):
                    g = df_in.groupby(group_cols)[metric_col].agg(["mean", "std", "count"]).reset_index()
                    g.columns = list(group_cols) + ["Mean", "StdDev", "Count"]
                    g["SEM"] = g.apply(lambda r: r["StdDev"] / (r["Count"] ** 0.5) if r["Count"] > 1 else 0, axis=1)
                    g["95% CI"] = g["SEM"] * 1.96
                    return g

                if df_multi.empty:
                    st.info("No constructs were tested against more than one target, so there "
                            "is nothing to compare for selectivity.")
                else:
                    # ---- MAIN GRAPH: Construct x Target ----
                    df_by_target = _sel_summary(df_multi, ["Construct", "Target"]).sort_values(["Construct", "Target"])
                    st.markdown(f"### Selectivity by Target — {y_label}")
                    fig_sel = px.bar(
                        df_by_target, x="Construct", y="Mean", color="Target",
                        barmode="group", error_y="95% CI", template="plotly_dark",
                        title=f"{y_label} Across Targets (Error Bars: 95% CI)",
                        labels={"Mean": y_label}, hover_data={"Count": True, "95% CI": ":.2f"},
                    )
                    if _is_norm:
                        fig_sel.add_hline(y=1, line_dash="dash", line_color="white", annotation_text="Pos Ctrl")
                    fig_sel.update_layout(height=600, xaxis_tickangle=-45)
                    st.plotly_chart(fig_sel, width='stretch', key="selectivity_main_plotly")

                    # ---- SECOND GRAPH: Construct x Target (Manufacturer) ----
                    df_by_src = _sel_summary(df_multi, ["Construct", "Target", "Source"])
                    if AGG_MODULE_OK:
                        df_by_src["Provenance"] = df_by_src["Source"].apply(aggregate_analysis.source_provenance)
                    else:
                        df_by_src["Provenance"] = "Unknown"
                    df_by_src["Target / Source"] = df_by_src["Target"].astype(str) + " (" + df_by_src["Source"].astype(str) + ")"
                    df_by_src = df_by_src.sort_values(["Construct", "Target", "Source"])
                    st.markdown(f"### Selectivity by Target & Manufacturer — {y_label}")
                    st.caption("Each target split by manufacturer/vendor, canonicalized via "
                               "`Protein-Analysis/vendor_manifest.csv` (e.g. MMP9 (Sino, Purchased) "
                               "vs MMP9 (Masoud, In-house)). Source is parsed from the sample name, "
                               "falling back to the trial-folder token.")
                    fig_src = px.bar(
                        df_by_src, x="Construct", y="Mean", color="Target / Source",
                        barmode="group", error_y="95% CI", template="plotly_dark",
                        title=f"{y_label} Across Targets x Manufacturer (Error Bars: 95% CI)",
                        labels={"Mean": y_label},
                        hover_data={"Target": True, "Source": True, "Provenance": True, "Count": True, "95% CI": ":.2f"},
                    )
                    if _is_norm:
                        fig_src.add_hline(y=1, line_dash="dash", line_color="white", annotation_text="Pos Ctrl")
                    fig_src.update_layout(height=600, xaxis_tickangle=-45)
                    st.plotly_chart(fig_src, width='stretch', key="selectivity_source_plotly")

                    with st.expander("📊 View Summary Data Tables"):
                        st.markdown("**By Target**")
                        st.dataframe(df_by_target, hide_index=True, width='stretch')
                        st.markdown("**By Target & Manufacturer**")
                        st.dataframe(df_by_src.drop(columns=["Target / Source"]), hide_index=True, width='stretch')

                    # ---- MASTER SELECTIVITY SUMMARY (ALL VARIANTS & METRICS) ----
                    st.divider()
                    with st.expander("📜 Master Selectivity Summary — ANOVA + Tukey HSD (All Variants & Metrics)", expanded=False):
                        st.markdown("For every multi-target construct, a one-way **ANOVA** compares its targets; "
                                    "when significant (p < 0.05), the significant **Tukey HSD** target pairs are listed below.")
                        st.caption(TUKEY_NOTE)
                        master_sel_findings = []
                        for _construct in sorted(_multi):
                            c_trials = df_multi[df_multi["Construct"] == _construct]
                            for _lbl, _info in SEL_METRICS.items():
                                _col = _info["norm"] if (_is_norm and _info["norm"]) else _info["raw"]
                                if _col not in c_trials.columns:
                                    continue
                                p_val = run_anova_p(c_trials, "Target", _col)
                                if not np.isnan(p_val) and p_val < 0.05:
                                    sig_targets = run_tukey_summary(c_trials, "Target", _col)
                                    if sig_targets is not None and not sig_targets.empty:
                                        for _, row in sig_targets.iterrows():
                                            master_sel_findings.append({
                                                "Construct": _construct,
                                                "Metric": _lbl,
                                                "Targets": f"{row['group1']} vs {row['group2']}",
                                                "Diff": f"{row['meandiff']:.4f}",
                                                "p-adj": f"{row['p-adj']:.4f}",
                                            })
                        if master_sel_findings:
                            st.dataframe(pd.DataFrame(master_sel_findings), hide_index=True, width='stretch')
                        else:
                            st.info("No significant selectivity findings across any constructs/metrics.")

                    # ---- INDIVIDUAL CONSTRUCT ----
                    st.divider()
                    st.subheader("🔍 Individual Construct Selectivity")
                    selected_variant = st.selectbox("Select Construct", sorted(_multi),
                                                    key="selectivity_variant_sel")

                    df_var = df_by_target[df_by_target["Construct"] == selected_variant].sort_values("Target")
                    fig_var = px.bar(
                        df_var, x="Target", y="Mean", color="Target", error_y="95% CI",
                        template="plotly_dark", title=f"Selectivity Profile: {selected_variant}",
                        labels={"Mean": y_label}, hover_data={"Count": True, "95% CI": ":.2f"},
                    )
                    if _is_norm:
                        fig_var.add_hline(y=1, line_dash="dash", line_color="white", annotation_text="Pos Ctrl")
                    fig_var.update_layout(height=450, showlegend=False)
                    st.plotly_chart(fig_var, width='stretch', key="selectivity_variant_plotly")

                    # Manufacturer breakdown for this construct
                    df_var_src = df_by_src[df_by_src["Construct"] == selected_variant].sort_values(["Target", "Source"])
                    if not df_var_src.empty and df_var_src["Source"].nunique() > 1:
                        fig_var_src = px.bar(
                            df_var_src, x="Target", y="Mean", color="Source", barmode="group",
                            error_y="95% CI", template="plotly_dark",
                            title=f"{selected_variant}: by Manufacturer",
                            labels={"Mean": y_label}, hover_data={"Count": True, "95% CI": ":.2f"},
                        )
                        if _is_norm:
                            fig_var_src.add_hline(y=1, line_dash="dash", line_color="white", annotation_text="Pos Ctrl")
                        fig_var_src.update_layout(height=450)
                        st.plotly_chart(fig_var_src, width='stretch', key="selectivity_variant_src_plotly")

                    # ANOVA across targets for this construct
                    c_trials = df_multi[df_multi["Construct"] == selected_variant]
                    p_val_sel = run_anova_p(c_trials, "Target", metric_col)
                    if not np.isnan(p_val_sel):
                        msg_sel = (f"📊 **Selectivity ANOVA**: **{selected_variant}** across targets "
                                   f"using **{y_label}**.  \n**p-value**: `{p_val_sel:.4e}` "
                                   f"({'Significant' if p_val_sel < 0.05 else 'Not Significant'})")
                        if p_val_sel < 0.05:
                            st.success(msg_sel)
                        else:
                            st.info(msg_sel)
                        if p_val_sel < 0.05:
                            sig_targets = run_tukey_summary(c_trials, "Target", metric_col)
                            if sig_targets is not None and not sig_targets.empty:
                                with st.expander(f"🔍 Significant Target Differences for {selected_variant} (Tukey HSD)"):
                                    st.caption(TUKEY_NOTE)
                                    st.dataframe(sig_targets[['group1', 'group2', 'meandiff', 'p-adj']], hide_index=True, width='stretch')
                    else:
                        st.caption("ANOVA could not be calculated for this construct across targets "
                                   "(requires multiple trials).")


else:
    if selected_file is None:
        st.info("👈 Please select an FCS file.")
    elif df is None:
        st.error("Failed to load FCS file.")
