import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import fcsparser
import os
import glob
from scipy import stats

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
    # Try simple glob first
    files = glob.glob(os.path.join(directory, "*.fcs"))
    return sorted(files)

# ---- SIDEBAR ----
st.sidebar.title("🧬 FCS Viewer")
st.sidebar.markdown("Explore raw flow cytometry data interactively.")

# Default Search Path
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
st.sidebar.info(f"**Files Found:** {len(files) if 'files' in locals() else 0}")

# ---- MAIN CONTENT ----

@st.cache_data
def calculate_density(x, y, subsample=5000):
    """Calculates density using Gaussian KDE for color mapping."""
    # Subsample for speed if needed
    if len(x) > subsample:
        idx = np.random.choice(len(x), subsample, replace=False)
        x_sub, y_sub = x.iloc[idx], y.iloc[idx]
    else:
        x_sub, y_sub = x, y
        
    # Stack and calc KDE
    values = np.vstack([x_sub, y_sub])
    kernel = stats.gaussian_kde(values)
    
    # Eval on full (or subsampled) set? Evaluating on full set is slow.
    # Let's eval on the plot points. If plot is subsampled, we eval on those.
    # We will compute density on the exact points passed to this function.
    z = kernel(np.vstack([x, y]))
    return z

def render_scatter(df, key_prefix):
    """Renders a single scatter plot interface."""
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
        
        # Gating
        show_gate = st.checkbox("Show Gates", value=False, key=f"{key_prefix}_gate")
        if show_gate:
            gate_x = st.number_input(f"Gate X ({x_axis})", value=0.0, key=f"{key_prefix}_gx")
            gate_y = st.number_input(f"Gate Y ({y_axis})", value=0.0, key=f"{key_prefix}_gy")

        # Subsample First (to keep UI responsive, but meaningful density requires enough points)
        # However, if we filter first, we might have few points left.
        # Ideally: Filter -> Subsample -> Density
        
        # Manual Limits
        with st.expander("Manual Axis Limits (Recalculate Density)"):
            use_manual_limits = st.checkbox("Enable Manual Limits", key=f"{key_prefix}_man_lim")
            l_col1, l_col2 = st.columns(2)
            
            # Default to data range
            x_min_d, x_max_d = float(df[x_axis].min()), float(df[x_axis].max())
            y_min_d, y_max_d = float(df[y_axis].min()), float(df[y_axis].max())
            
            with l_col1:
                x_min = st.number_input(f"X Min ({x_axis})", value=x_min_d, key=f"{key_prefix}_xmin")
                x_max = st.number_input(f"X Max ({x_axis})", value=x_max_d, key=f"{key_prefix}_xmax")
            with l_col2:
                y_min = st.number_input(f"Y Min ({y_axis})", value=y_min_d, key=f"{key_prefix}_ymin")
                y_max = st.number_input(f"Y Max ({y_axis})", value=y_max_d, key=f"{key_prefix}_ymax")
                
    with col2:
        # 1. Filter Data (if manual limits)
        if use_manual_limits:
            stats_df = df[
                (df[x_axis] >= x_min) & (df[x_axis] <= x_max) &
                (df[y_axis] >= y_min) & (df[y_axis] <= y_max)
            ]
        else:
            stats_df = df
            
        # 2. Subsample for plotting
        if len(stats_df) > subsample_n:
            plot_df = stats_df.sample(n=subsample_n, random_state=42).copy()
        else:
            plot_df = stats_df.copy()
            
        if len(plot_df) == 0:
            st.warning("No events in current range.")
            return

        # 3. Transformations
        if log_x:
            # Avoid log(<=0)
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

        # 4. Colors (Density Recalculation on Filtered Data)
        if color_mode == "Density":
            with st.spinner("Recalculating relative density..."):
                # If we have very few points, KDE might fail or look weird
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

        # Add Gates
        if show_gate:
            # Stats (Quadrants) based on FILTERED view
            # If we want stats relative to the VIEW, we leverage stats_df or plot_df?
            # Usually users expect stats on what they see.
            
            # Count events in quadrants
            q1 = len(plot_df[(plot_df[x_axis] < gx_val) & (plot_df[y_axis] >= gy_val)]) # Top Left
            q2 = len(plot_df[(plot_df[x_axis] >= gx_val) & (plot_df[y_axis] >= gy_val)]) # Top Right
            q3 = len(plot_df[(plot_df[x_axis] >= gx_val) & (plot_df[y_axis] < gy_val)]) # Bottom Right
            q4 = len(plot_df[(plot_df[x_axis] < gx_val) & (plot_df[y_axis] < gy_val)]) # Bottom Left
            total = len(plot_df)
            
            # Add Lines
            fig.add_vline(x=gx_val, line_dash="dash", line_color="white")
            fig.add_hline(y=gy_val, line_dash="dash", line_color="white")
            
            # Annotation
            st.code(f"""
            Q1 (UL): {q1/total*100:.1f}% | Q2 (UR): {q2/total*100:.1f}%
            Q4 (LL): {q4/total*100:.1f}% | Q3 (LR): {q3/total*100:.1f}%
            (Relative to visible/subsampled events)
            """)
            
        fig.update_layout(
            autosize=True, 
            height=500,
            xaxis_title=plot_xlabel,
            yaxis_title=plot_ylabel,
            margin=dict(l=0, r=0, t=0, b=0)
        )
        st.plotly_chart(fig, use_container_width=True)


if selected_file:
    st.header(f"📂 {os.path.basename(selected_file)}")
    
    with st.spinner("Loading FCS file..."):
        meta, df = load_fcs(selected_file)
        
    if df is not None:
        # Top Metrics
        c1, c2, c3 = st.columns(3)
        c1.metric("Total Events", f"{len(df):,}")
        st.sidebar.markdown("---")
        
        # Tabs
        tab_scatter, tab_hist, tab_raw, tab_meta = st.tabs(["📊 Scatter Plot", "📈 Histogram", "📄 Raw Data", "ℹ️ Metadata"])
        
        # ---- 1. SCATTER PLOT ----
        with tab_scatter:
            
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
                
        # ---- 2. HISTOGRAM ----
        with tab_hist:
            col1, col2 = st.columns([1, 4])
            
            with col1:
                st.subheader("Settings")
                target_chs = st.multiselect("Channels", df.columns, default=["FITC-A", "APC-A"] if "FITC-A" in df.columns else [df.columns[0]])
                log_hist = st.checkbox("Log Scale (X-Axis)", value=True)
                bins = st.slider("Bins", 10, 200, 100)
                
            with col2:
                if target_chs:
                    st.subheader("Signal Distribution")
                    plot_data = df.copy()
                    long_df = plot_data[target_chs].melt(var_name="Channel", value_name="Intensity")
                    
                    if log_hist:
                        # Clip to a small positive value so they don't disappear on log scale
                        # Standard flow cytometry often treats < 1 as 1 for log viz
                        long_df["Intensity"] = long_df["Intensity"].clip(lower=1)
                    
                    fig_hist = px.histogram(long_df, x="Intensity", color="Channel", 
                                            barmode="overlay", nbins=bins, 
                                            log_x=log_hist,
                                            opacity=0.7,
                                            template="plotly_dark")
                    
                    fig_hist.update_layout(height=500, margin=dict(l=0, r=0, t=30, b=0), xaxis_title="Intensity (Log)" if log_hist else "Intensity")
                    st.plotly_chart(fig_hist, use_container_width=True)
                else:
                    st.info("Select channels to view histograms.")

        # ---- 3. RAW DATA ----
        with tab_raw:
            st.dataframe(df, use_container_width=True, height=600)
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("Download CSV", csv, "fcs_data.csv", "text/csv")
            
        # ---- 4. METADATA ----
        with tab_meta:
            st.subheader("File Metadata")
            
            # Convert to neat DataFrame
            meta_items = []
            for k, v in meta.items():
                # Filter out standard Pn parameters since we have a separate table for them
                if not k.startswith("$P"):
                    meta_items.append({"Parameter": k, "Value": str(v)})
            
            st.dataframe(pd.DataFrame(meta_items), use_container_width=True, hide_index=True)
            
            st.subheader("Channel Parameters")
            params = []
            i = 1
            while f"$P{i}N" in meta:
                p_name = meta.get(f"$P{i}N")
                p_desc = meta.get(f"$P{i}S", "-")
                p_range = meta.get(f"$P{i}R", "-")
                params.append({"Index": i, "Name": p_name, "Description": p_desc, "Range": p_range})
                i += 1
            
            if params:
                st.dataframe(pd.DataFrame(params), use_container_width=True, hide_index=True)

    else:
        st.error("Failed to load FCS file.")
else:
    st.info("👈 Please select an FCS file.")

