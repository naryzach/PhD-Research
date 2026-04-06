import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import fcsparser
import os
import glob
from scipy import stats
from scipy.stats import gaussian_kde
from scipy.spatial import ConvexHull
from matplotlib.path import Path

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

# ---- CONFIGURATION ----
MIN_FSC = 500000
MIN_SSC = 20000
UPPER_FILTER_PERCENTILE = 95
GATE_TARGET_FRACTION = 0.90
NC_CUTOFF_PERCENTILE = 99.5

# ---- HELPER FUNCTIONS ----

@st.cache_data
def load_fcs(file_path):
    """Loads an FCS file and returns metadata and dataframe."""
    try:
        meta, data = fcsparser.parse(file_path, reformat_meta=True)
        return meta, data
    except Exception as e:
        return None, None

def get_files(directory):
    """Recursively finds FCS files."""
    if not os.path.exists(directory):
        return []
    files = glob.glob(os.path.join(directory, "*.fcs"))
    return sorted(files)

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
def learn_pentagon_gate(df, fsc_col="FSC-A", ssc_col="SSC-A", fraction=GATE_TARGET_FRACTION):
    origin_filter = (df[fsc_col] > MIN_FSC) & (df[ssc_col] > MIN_SSC)
    core_df = df[origin_filter]
    
    if len(core_df) < 100:
        core_df = df
        
    x = core_df[fsc_col].values
    y = core_df[ssc_col].values
    
    fsc_upper = np.percentile(x, UPPER_FILTER_PERCENTILE)
    ssc_upper = np.percentile(y, UPPER_FILTER_PERCENTILE)
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

def apply_polygon_gate(df, path, fsc_col="FSC-A", ssc_col="SSC-A"):
    if len(df) == 0:
        return df
    points = np.vstack([df[fsc_col], df[ssc_col]]).T
    poly_mask = path.contains_points(points)
    origin_mask = (df[fsc_col] > MIN_FSC) & (df[ssc_col] > MIN_SSC)
    return df[poly_mask & origin_mask]

@st.cache_data(show_spinner=False)
def analyze_folder_controls(directory, fsc_col="FSC-A", ssc_col="SSC-A", expr_col="FITC-A", bind_col="APC-A", nc_percentile=NC_CUTOFF_PERCENTILE):
    all_files = glob.glob(os.path.join(directory, "*.fcs"))
    neg_files = [f for f in all_files if any(os.path.basename(f).startswith(p) for p in ["NC", "Negative Control"])]
    pos_files = [f for f in all_files if any(os.path.basename(f).startswith(p) for p in ["Positive Control", "TIMP"])]
    
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
    }
    
    if neg_files:
        _, d_neg_learn = load_fcs(neg_files[0])
        if d_neg_learn is not None and fsc_col in d_neg_learn.columns and ssc_col in d_neg_learn.columns:
            path, verts = learn_pentagon_gate(d_neg_learn, fsc_col, ssc_col)
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
                    results["has_nc"] = True

    if pos_files and results["has_nc"]:
        pos_dfs = []
        for f in pos_files:
            _, d = load_fcs(f)
            if d is not None:
                d_gated = apply_polygon_gate(d, results["pentagon_path"], fsc_col, ssc_col)
                if expr_col in d_gated.columns and bind_col in d_gated.columns:
                    d_gated = d_gated[(d_gated[expr_col] > 0) & (d_gated[bind_col] > 0)]
                    pos_dfs.append(d_gated)
                
        if pos_dfs:
            pos_concat = pd.concat(pos_dfs)
            if len(pos_concat) > 0:
                results["pos_mfi_expr"] = pos_concat[expr_col].median()
                results["pos_mfi_bind"] = pos_concat[bind_col].median()
                results["has_pc"] = True

    return results

@st.cache_data
def calculate_density(x, y, subsample=5000):
    if len(x) > subsample:
        idx = np.random.choice(len(x), subsample, replace=False)
        if hasattr(x, 'iloc'):
            x_sub, y_sub = x.iloc[idx], y.iloc[idx]
        else:
            x_sub, y_sub = x[idx], y[idx]
    else:
        x_sub, y_sub = x, y
        
    values = np.vstack([x_sub, y_sub])
    kernel = stats.gaussian_kde(values)
    z = kernel(np.vstack([x, y]))
    return z

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
        st.plotly_chart(fig, use_container_width=True)

# ---- SIDEBAR ----
st.sidebar.title("🧬 FCS Viewer")
st.sidebar.markdown("Explore raw flow cytometry data interactively.")

default_path = os.path.join(os.path.dirname(__file__), "../Local/TwistBio_FCS")
default_path = os.path.abspath(default_path)
search_path = st.sidebar.text_input("Data Directory", value=default_path)

if os.path.exists(search_path):
    files = get_files(search_path)
    if files:
        selected_file = st.sidebar.selectbox("Select FCS File", files, format_func=os.path.basename)
    else:
        st.sidebar.warning("No .fcs files found in directory.")
        selected_file = None
else:
    st.sidebar.error("Directory not found.")
    selected_file = None

st.sidebar.markdown("---")

if selected_file:
    with st.spinner("Analyzing folder controls..."):
        ctrls = analyze_folder_controls(search_path)
    
    st.sidebar.subheader("Gating Context")
    if ctrls["has_nc"]:
        st.sidebar.success("Global Negative Control applied.")
    else:
        st.sidebar.warning("Local Isolated Gating applied.")
        
    if ctrls["has_pc"]:
        st.sidebar.info("Global Positive Control data available.")
        
    st.sidebar.subheader("Threshold Settings")
    
    meta, df = load_fcs(selected_file)
    
    if df is not None:
        expr_col = "FITC-A" if "FITC-A" in df.columns else df.columns[0]
        bind_col = "APC-A" if "APC-A" in df.columns else df.columns[1]
        
        if ctrls["has_nc"]:
            def_thresh_expr = float(ctrls["thresh_expr"])
            def_thresh_bind = float(ctrls["thresh_bind"])
        else:
            def_thresh_expr = float(np.percentile(df[expr_col], NC_CUTOFF_PERCENTILE))
            def_thresh_bind = float(np.percentile(df[bind_col], NC_CUTOFF_PERCENTILE))

        thresh_expr_val = st.sidebar.number_input("Expression Thresh", value=def_thresh_expr)
        thresh_bind_val = st.sidebar.number_input("Binding Thresh", value=def_thresh_bind)

        log_thresh_expr = np.log10(max(1, thresh_expr_val))
        log_thresh_bind = np.log10(max(1, thresh_bind_val))

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
            pent_path, pent_verts = learn_pentagon_gate(df, fsc_col, ssc_col)
            
    df_sing = apply_polygon_gate(df, pent_path, fsc_col, ssc_col)
    
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
            
    # TABS
    tab_main, tab_custom, tab_raw, tab_meta = st.tabs(["📊 Main Analysis", "🛠️ Custom Plots", "📄 Raw Data", "ℹ️ Metadata"])
    
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
                    fig1.add_trace(go.Scatter(
                        x=pent_verts[:,0], y=pent_verts[:,1],
                        mode='lines',
                        line=dict(color='cyan', dash='dash', width=2),
                        name='Pentagon Gate'
                    ))
                fig1.update_layout(
                    title=f"Original Gate ({len(df_sing)}/{len(df)})",
                    xaxis_title=fsc_col,
                    yaxis_title=ssc_col,
                    template='plotly_dark',
                    height=450, margin=dict(l=0,r=0,b=0)
                )
                st.plotly_chart(fig1, use_container_width=True)

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
                    title="Expression Distribution (FITC)",
                    xaxis_title=f"Log10 {expr_col}",
                    yaxis_title="Density",
                    height=450, margin=dict(l=0,r=0,b=0),
                    showlegend=False
                )
                st.plotly_chart(fig2, use_container_width=True)
                
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
                st.plotly_chart(fig3, use_container_width=True)

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
                    title="Binding Distribution (APC)",
                    xaxis_title=f"Log10 {bind_col}",
                    yaxis_title="Density",
                    height=450, margin=dict(l=0,r=0,b=0),
                    showlegend=False
                )
                st.plotly_chart(fig4, use_container_width=True)
        else:
            st.warning("No gated events found.")

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
        st.dataframe(df, use_container_width=True, height=600)
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button("Download CSV", csv, "fcs_data.csv", "text/csv")
        
    with tab_meta:
        st.subheader("File Metadata")
        meta_items = []
        for k, v in meta.items():
            if not k.startswith("$P"):
                meta_items.append({"Parameter": k, "Value": str(v)})
        st.dataframe(pd.DataFrame(meta_items), use_container_width=True, hide_index=True)
        
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
            st.dataframe(pd.DataFrame(params), use_container_width=True, hide_index=True)

else:
    if selected_file is None:
        st.info("👈 Please select an FCS file.")
    elif df is None:
        st.error("Failed to load FCS file.")
