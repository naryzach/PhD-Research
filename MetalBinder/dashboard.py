import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os
import numpy as np
from stmol import showmol
import py3Dmol
import math

# Set page config - MUST be first streamlit command
st.set_page_config(page_title="Lanm Output Dashboard", layout="wide")

# --- Authentication ---
def check_password():
    """Returns True if the user has entered the correct password."""
    def password_entered():
        """Checks whether a password entered by the user is correct."""
        if st.session_state["password"] == st.secrets["admin"]["password"]:
            st.session_state["password_correct"] = True
            del st.session_state["password"]  # don't store password
        else:
            st.session_state["password_correct"] = False

    if st.session_state.get("password_correct", False):
        return True

    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        st.title("🔒 Password Protected")
        st.text_input(
            "Enter password to access the Lanm Output Dashboard", 
            type="password", 
            on_change=password_entered, 
            key="password"
        )
        if "password_correct" in st.session_state and not st.session_state["password_correct"]:
            st.error("😕 Password incorrect")
    return False

if not check_password():
    st.stop()

# --- Constants & Heuristics ---
def calculate_binding_probability(row):
    """
    Heuristic score (0-1) for metal binding quality.
    - Radius: Ideal ~2.45A (Gaussian decay)
    - CN: Bonus for 6-8 coordinating oxygens
    - Bidentate: Bonus for at least one bidentate ligand
    - Confidence: Scaled pLDDT
    """
    # 1. Radius Factor (Gaussian around 2.45A)
    # sigma=0.5 means at 3.0A (0.55 away) it's ~0.55
    r = row.get('binding_radius_A', 5.0)
    if pd.isna(r) or r <= 0: return 0.0
    r_score = math.exp(-((r - 2.45)**2) / (2 * 0.5**2))
    
    # 2. Coordination Number Factor
    cn = row.get('coordination_number', 0)
    if pd.isna(cn): cn = 0
    if 6 <= cn <= 8: cn_score = 1.0
    elif cn > 0: cn_score = 0.5
    else: cn_score = 0.0
    
    # 3. Bidentate Factor
    bid = row.get('bidentate_count', 0)
    bid_score = 1.0 if bid >= 1 else 0.0
    
    # 4. Confidence Factor (pLDDT)
    p = row.get('plddt', 0)
    p_score = p / 100.0 if p > 1.0 else p
    
    # Weighted average
    total = (r_score * 0.4) + (cn_score * 0.2) + (bid_score * 0.1) + (p_score * 0.3)
    return round(total, 3)

st.title("🧬 Lanm Output Dashboard")
st.markdown("Explore metal-binding protein designs with advanced structural and sequence analysis.")

# --- Categorization ---
ALKALI = ['LI', 'NA', 'K', 'RB', 'CS', 'FR']
ALKALINE_EARTH = ['BE', 'MG', 'CA', 'SR', 'BA', 'RA']
# Rare Earth: Sc, Y, Lanthanides (La-Lu), and Actinides (Ac-Lr)
RARE_EARTH = [
    'SC', 'Y', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU',
    'AC', 'TH', 'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR'
]
TRANSITION = [
    'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD',
    'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG'
]
POST_TRANSITION = ['AL', 'GA', 'IN', 'SN', 'TL', 'PB', 'BI', 'PO']

def get_category(ion):
    ion_upper = str(ion).upper().replace('2+', '').replace('3+', '')
    if ion_upper in ALKALI: return "Alkali"
    if ion_upper in ALKALINE_EARTH: return "Alkaline Earth"
    if ion_upper in RARE_EARTH: return "Rare Earth"
    if ion_upper in TRANSITION: return "Transition"
    if ion_upper in POST_TRANSITION: return "Post-transition"
    return "Other"

# --- Data Loading ---
@st.cache_data
def get_cif_data(cif_path):
    if not os.path.exists(cif_path):
        return None
    with open(cif_path, "r") as f:
        return f.read()

@st.cache_data
def load_data():
    # Get the directory where dashboard.py is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Potential data directories (constructed relative to this script)
    potential_paths = [
        os.path.join(script_dir, "..", "Local", "lanm_output"),
        os.path.join(script_dir, "lanm_data")
    ]
    base_path = None
    
    # Try to find a valid base path
    for path in potential_paths:
        if os.path.exists(os.path.join(path, "global_sequence_catalog.csv")):
            base_path = path
            break
            
    # Fallback to the first path if none found (to avoid base_path being None)
    if base_path is None:
        base_path = potential_paths[0]
        
    catalog_path = os.path.join(base_path, "global_sequence_catalog.csv")
    full_seq_path = os.path.join(base_path, "full_sequences_log.csv")
    
    df = pd.read_csv(catalog_path) if os.path.exists(catalog_path) else pd.DataFrame()
    full_seq_df = pd.read_csv(full_seq_path) if os.path.exists(full_seq_path) else pd.DataFrame()
    
    # Check for ions presence even if not in catalog
    ions = []
    if os.path.exists(base_path):
        ions = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d)) and d.isupper() and len(d) <= 2]
    
    if not df.empty:
        # Update Metal Category
        df['metal_category'] = df['metal_ion'].apply(get_category)
        # Add Binding Probability
        df['binding_probability'] = df.apply(calculate_binding_probability, axis=1)
    
    return df, full_seq_df, ions, base_path

@st.cache_data
def get_b_metrics(cif_path):
    """Parses B-factors from CIF file and determines if they represent pLDDT or Error."""
    if not os.path.exists(cif_path):
        return 0.0, 1.0, False
        
    with open(cif_path, "r") as f:
        lines = f.readlines()
        
    col_b = -1
    col_atom = -1
    col_count = 0
    b_factors = []
    header_section = True
    
    for line in lines:
        line = line.strip()
        if not line: continue
        if line.startswith('_atom_site.'):
            if '_atom_site.B_iso_or_equiv' in line:
                col_b = col_count
            if '_atom_site.label_atom_id' in line or '_atom_site.auth_atom_id' in line:
                col_atom = col_count
            col_count += 1
        elif line.startswith('ATOM') or line.startswith('HETATM'):
            header_section = False
            parts = line.split()
            if col_b != -1 and col_atom != -1 and len(parts) > max(col_b, col_atom):
                atom_name = parts[col_atom].strip('"').strip("'")
                if atom_name == "CA":
                    try:
                        b_factors.append(float(parts[col_b]))
                    except ValueError: continue
        elif not header_section and line.startswith('#'):
            break
            
    if not b_factors:
        return 0.0, 1.0, False # Default pLDDT False (roygb)
    
    b_min, b_max = min(b_factors), max(b_factors)
    avg_b = sum(b_factors) / len(b_factors)
    
    # Adaptive Inversion Logic:
    if b_max <= 1.1: # 0-1 scale
        invert = avg_b < 0.3 # Mostly small values = Error/Distance
    else: # 0-100 scale
        invert = avg_b < 30 # Likely RMSD or B-factor
        
    return b_min, b_max, invert

@st.cache_data
def load_cross_docking_data(base_path):
    cross_dock_path = os.path.join(base_path, "cross_docking_catalog.csv")
    if os.path.exists(cross_dock_path):
        return pd.read_csv(cross_dock_path)
    return pd.DataFrame()

df, full_seq_df, available_ions, base_path = load_data()

if df.empty:
    st.warning("⚠️ No data found in global_sequence_catalog.csv. Please ensure the pipeline has generated results.")
    st.stop()

# Add categories
df['metal_category'] = df['metal_ion'].apply(get_category)

# --- Helper Functions for Plots ---
# --- Cached Plotting Functions ---
@st.cache_resource
def get_residue_heatmap_fig(data, title, color_scale="YlOrRd"):
    if data.empty:
        return None
    
    max_len = data['loop_length'].max()
    aa_list = list("ACDEFGHIKLMNPQRSTVWY")
    freq_matrix = pd.DataFrame(0.0, index=aa_list, columns=range(max_len))
    
    for _, row in data.iterrows():
        seq = row['loop_sequence']
        for i, aa in enumerate(seq):
            if aa in freq_matrix.index and i < max_len:
                freq_matrix.loc[aa, i] += 1
    
    # Normalize
    freq_matrix = freq_matrix.div(freq_matrix.sum(axis=0), axis=1).fillna(0)
    
    fig = px.imshow(
        freq_matrix,
        labels=dict(x="Loop Position", y="Amino Acid", color="Frequency"),
        x=[f"Pos {i}" for i in range(max_len)],
        y=aa_list,
        color_continuous_scale=color_scale,
        title=title
    )
    return fig

@st.cache_resource
def get_eval_plot(df, selected_shape, selected_color, rad_range_vals, rmsd_range_vals):
    shape_options = {
        "Loop Index": "loop_index",
        "Bidentate Count": "bidentate_count",
        "Metal Ion": "metal_ion",
        "Configuration": "config_index",
        "Motif Match": "motif_match"
    }
    color_options = {
        "Metal Ion": "metal_ion",
        "Binding Probability": "binding_probability",
        "pLDDT": "plddt",
        "Loop Index": "loop_index"
    }
    
    # Filter for existing columns
    active_shape_options = {k: v for k, v in shape_options.items() if v in df.columns}
    symbol_col = active_shape_options.get(selected_shape, None)
    
    active_color_options = {k: v for k, v in color_options.items() if v in df.columns}
    color_col = active_color_options.get(selected_color, "metal_ion")
    
    fig = px.scatter(
        df, 
        x='loop_rmsd', 
        y='binding_radius_A', 
        color=color_col, 
        size='plddt' if 'plddt' in df.columns else None,
        symbol=symbol_col,
        hover_data=['design_id', 'loop_index', 'loop_sequence', 'binding_probability', 'coordination_number', 'net_charge'] if 'coordination_number' in df.columns else ['design_id', 'loop_index', 'loop_sequence', 'binding_probability'],
        color_continuous_scale="Viridis" if color_col != "metal_ion" else None,
        size_max=8, # Reduced from default 20
        labels={'loop_rmsd': 'Individual Loop RMSD (Å)', 'binding_radius_A': 'Binding Radius (Å)', 'binding_probability': 'Probability'},
        title="Candidate Evaluation (Thresholds: RMSD < 1.5, Radius 2.3-2.6)"
    )
    
    # If size is NOT mapped to pLDDT, set a small fixed size
    if 'plddt' not in df.columns:
        fig.update_traces(marker=dict(size=6))
    
    # Force axis range to match filters
    if rad_range_vals:
        fig.update_yaxes(range=[rad_range_vals[0]-0.1, rad_range_vals[1]+0.1])
    if rmsd_range_vals:
        fig.update_xaxes(range=[0, max(1.5, df['loop_rmsd'].max() + 0.2)])
        
    # Add horizontal region for ideal radius
    fig.add_hrect(y0=2.3, y1=2.6, line_width=0, fillcolor="green", opacity=0.1)
    # Add vertical line for RMSD threshold
    fig.add_vline(x=1.5, line_dash="dash", line_color="red", opacity=0.6)
    
    # Improve legend layout on the right side to prevent overlapping
    fig.update_layout(
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1.0,
            xanchor="left",
            x=1.15, # Offset discrete legend further right
            font=dict(size=10),
            itemsizing="constant"
        ),
        margin=dict(r=200) # Increased margin for two side-by-side legends
    )
    
    # Position colorbar specifically if it exists (for continuous scales)
    fig.update_coloraxes(
        colorbar=dict(
            thickness=15,
            x=1.02, # Keep colorbar closer to the plot
            y=0.5,
            len=0.75,
            title=dict(side="right")
        )
    )
    return fig

@st.cache_resource
def get_accuracy_plots(df):
    fig_violin_rmsd = px.violin(df, x='metal_ion', y='loop_rmsd', color='metal_ion', box=True, points="all", title="Loop RMSD Distribution")
    
    fig_box_rad = px.box(df, x='metal_ion', y='binding_radius_A', color='metal_ion', points="all", title="Binding Radius Distribution")
    fig_box_rad.add_hrect(y0=2.3, y1=2.6, line_width=0, fillcolor="green", opacity=0.1)
    
    # Success Rate
    p_max = df['plddt'].max()
    p_thresh = 0.8 if p_max <= 1.0 else 80.0
    success_df = df.copy()
    success_df['is_success'] = (success_df['plddt'] >= p_thresh) & (success_df['overall_rmsd'] <= 2.0)
    success_rate = success_df.groupby('metal_ion')['is_success'].mean().reset_index()
    success_rate['Success Rate (%)'] = success_rate['is_success'] * 100
    fig_success = px.bar(success_rate, x='metal_ion', y='Success Rate (%)', color='metal_ion', title="Binding Success Rate by Ion Type")
    
    return fig_violin_rmsd, fig_box_rad, fig_success

@st.cache_resource
def get_loop_dist_plots(df):
    fig_hist = px.histogram(
        df, x="loop_length", color="metal_ion", barmode="group",
        title="Loop Length Distributions",
        labels={"loop_length": "Loop Length", "count": "Frequency"}
    )
    fig_len_plddt = px.scatter(
        df, x='loop_length', y='plddt', color='metal_ion', trendline="ols", 
        title="Loop Length Impact on Confidence"
    )
    return fig_hist, fig_len_plddt

@st.cache_resource
def get_radar_chart(df, design_a, design_b, loop_a=None, loop_b=None):
    metrics = ['binding_probability', 'plddt', 'overall_rmsd', 'loop_rmsd', 'binding_radius_A', 'coordination_number']
    metrics = [m for m in metrics if m in df.columns]
    
    def normalize(val, col, global_df):
        v_min = global_df[col].min()
        v_max = global_df[col].max()
        if v_max == v_min: return 0.5
        return (val - v_min) / (v_max - v_min)
        
    fig = go.Figure()
    for d_id, l_idx, color in [(design_a, loop_a, 'blue'), (design_b, loop_b, 'red')]:
        sub = df[df['design_id'] == d_id]
        if l_idx is not None:
            sub = sub[sub['loop_index'] == l_idx]
        row = sub.iloc[0]
        
        name = f"{d_id}" + (f" (Loop {l_idx})" if l_idx is not None else "")
        r_values = [normalize(row[m], m, df) for m in metrics]
        r_values.append(r_values[0])
        theta_labels = metrics + [metrics[0]]
        fig.add_trace(go.Scatterpolar(r=r_values, theta=theta_labels, fill='toself', name=name, line_color=color))
        
    fig.update_layout(polar=dict(radialaxis=dict(visible=True, range=[0, 1])), showlegend=True, title="Relative Metric Strength")
    return fig

@st.cache_resource
def get_analytics_plots(df, numeric_cols):
    corr_fig = None
    if len(numeric_cols) > 1:
        corr = df[numeric_cols].corr()
        corr_fig = px.imshow(corr, text_auto=True, color_continuous_scale='RdBu_r', zmin=-1, zmax=1, title="Property Correlation")
        
    pareto_fig = px.scatter(df, x='overall_rmsd', y='plddt', color='metal_ion', hover_data=['design_id'], title="pLDDT vs. Overall RMSD")
    return corr_fig, pareto_fig

@st.cache_resource
def get_category_plots(df):
    fig_plddt = px.box(df, x='metal_category', y='plddt', color='metal_category', notched=True, points="all", title="Design Confidence by Category")
    fig_rmsd = px.box(df, x='metal_category', y='loop_rmsd', color='metal_category', notched=True, points="all", title="Loop Accuracy by Category")
    
    # Success per category
    unique_cats = df['metal_category'].unique()
    p_max = df['plddt'].max()
    p_thresh = 0.8 if p_max <= 1.0 else 80.0
    success_data = []
    for cat in unique_cats:
        sub = df[df['metal_category'] == cat]
        rate = ((sub['plddt'] >= p_thresh) & (sub['overall_rmsd'] <= 2.0)).mean() * 100
        success_data.append({'Category': cat, 'Success Rate (%)': rate})
    fig_success = px.bar(pd.DataFrame(success_data), x='Category', y='Success Rate (%)', color='Category', title="Success Rate per Category")
    
    fig_radius = px.violin(df, x='metal_category', y='binding_radius_A', color='metal_category', box=True, title="Binding Radius Variance")
    fig_radius.add_hrect(y0=2.3, y1=2.6, line_width=0, fillcolor="green", opacity=0.1)
    
    return fig_plddt, fig_rmsd, fig_success, fig_radius

@st.cache_resource
def get_res_stability_plot(df, ion):
    sub = df[df['metal_ion'] == ion]
    res_data = []
    for _, row in sub.iterrows():
        for res in set(row['loop_sequence']):
            res_data.append({'Residue': res, 'pLDDT': row['plddt'], 'Metal': row['metal_ion']})
    res_df = pd.DataFrame(res_data)
    fig = px.box(res_df, x='Residue', y='pLDDT', color='Metal', title=f"Residue Influence on Stability ({ion})")
    return fig

# --- Sidebar Filters ---
st.sidebar.header("Filters")

# Session State for Optimal Presets
if 'optimal_filters' not in st.session_state:
    st.session_state.optimal_filters = False

def set_optimal():
    st.session_state.optimal_filters = True
    # Reset sliders to optimal values with scale detection
    if 'plddt' in df.columns:
        p_max = df['plddt'].max()
        # Scale to 0.8-1.0 if data is 0-1, otherwise 80-100
        val = (80.0, 100.0) if p_max > 1.0 else (0.8, 1.0)
        # Clip to data bounds to avoid StreamlitValueAboveMaxError
        st.session_state.plddt_slider = (max(float(df['plddt'].min()), val[0]), min(float(df['plddt'].max()), val[1]))
        
    if 'overall_rmsd' in df.columns:
        st.session_state.overall_rmsd_slider = (max(float(df['overall_rmsd'].min()), 0.0), min(float(df['overall_rmsd'].max()), 2.0))
        
    if 'binding_radius_A' in df.columns:
        st.session_state.binding_radius_A_slider = (max(float(df['binding_radius_A'].min()), 2.1), min(float(df['binding_radius_A'].max()), 2.8))
        
    if 'coordination_number' in df.columns:
        st.session_state.coordination_number_slider = (max(float(df['coordination_number'].min()), 7.0), min(float(df['coordination_number'].max()), 9.0))
        
    if 'net_charge' in df.columns:
        st.session_state.net_charge_slider = (max(float(df['net_charge'].min()), -5.0), min(float(df['net_charge'].max()), -3.0))

if st.sidebar.button("Set Optimal Presets", on_click=set_optimal):
    st.sidebar.success("Optimal filters applied!")

# Metal Category Filter
st.sidebar.subheader("Metal Groups")
available_categories = sorted(df['metal_category'].unique())
category_choice = st.sidebar.multiselect("Metal Category", options=available_categories, default=available_categories)

# Ion selection based on category
possible_ions_df = df[df['metal_category'].isin(category_choice)]
possible_ions = possible_ions_df['metal_ion'].unique()
selected_ions = st.sidebar.multiselect("Select Target Ions", options=sorted(possible_ions), default=sorted(possible_ions))

# EF Loop selection
st.sidebar.subheader("Binding Site Selection")
available_loops = sorted(df['loop_index'].unique().tolist()) if 'loop_index' in df.columns else [1, 2, 3, 4]
selected_loops = st.sidebar.multiselect("Select EF Loops", options=available_loops, default=available_loops)

# Filter dataframe
filtered_df = df.copy()

if selected_ions:
    filtered_df = filtered_df[filtered_df['metal_ion'].isin(selected_ions)]
filtered_df = filtered_df[filtered_df['metal_category'].isin(category_choice)]

if 'loop_index' in filtered_df.columns:
    filtered_df = filtered_df[filtered_df['loop_index'].isin(selected_loops)]

# Metric filters (checking if columns exist)
st.sidebar.subheader("Quality Metrics")

def get_range_filter(df, col, label, key):
    if col in df.columns:
        m_min, m_max = float(df[col].min()), float(df[col].max())
        if m_min == m_max:
            return (m_min, m_max)
            
        # Initialize or FIX session state for slider
        if key not in st.session_state:
            st.session_state[key] = (m_min, m_max)
        else:
            # Persistent session state might have stale values (e.g. 80.0 when max is 0.87)
            # We MUST clip the existing state to current data bounds to prevent StreamlitValueAboveMaxError
            curr_val = st.session_state[key]
            new_val = (
                max(m_min, min(m_max, float(curr_val[0]))),
                max(m_min, min(m_max, float(curr_val[1])))
            )
            st.session_state[key] = new_val
            
        return st.sidebar.slider(label, m_min, m_max, key=key)
    return None

plddt_range = get_range_filter(df, 'plddt', "pLDDT Range", "plddt_slider")
rmsd_range = get_range_filter(df, 'overall_rmsd', "Overall RMSD Range", "overall_rmsd_slider")
rad_range = get_range_filter(df, 'binding_radius_A', "Binding Radius (A) Range", "binding_radius_A_slider")
cn_range = get_range_filter(df, 'coordination_number', "Coordination Number Range", "coordination_number_slider")
charge_range = get_range_filter(df, 'net_charge', "Net Charge Range", "net_charge_slider")
prob_range = get_range_filter(df, 'binding_probability', "Binding Probability Range", "prob_slider")

if selected_ions:
    filtered_df = filtered_df[filtered_df['metal_ion'].isin(selected_ions)]
filtered_df = filtered_df[filtered_df['metal_category'].isin(category_choice)]

if plddt_range:
    filtered_df = filtered_df[(filtered_df['plddt'] >= plddt_range[0]) & (filtered_df['plddt'] <= plddt_range[1])]
if rmsd_range:
    filtered_df = filtered_df[(filtered_df['overall_rmsd'] >= rmsd_range[0]) & (filtered_df['overall_rmsd'] <= rmsd_range[1])]
if rad_range:
    filtered_df = filtered_df[(filtered_df['binding_radius_A'] >= rad_range[0]) & (filtered_df['binding_radius_A'] <= rad_range[1])]
if cn_range:
    filtered_df = filtered_df[(filtered_df['coordination_number'] >= cn_range[0]) & (filtered_df['coordination_number'] <= cn_range[1])]
if charge_range:
    filtered_df = filtered_df[(filtered_df['net_charge'] >= charge_range[0]) & (filtered_df['net_charge'] <= charge_range[1])]
if prob_range:
    filtered_df = filtered_df[(filtered_df['binding_probability'] >= prob_range[0]) & (filtered_df['binding_probability'] <= prob_range[1])]

# --- Export (In Sidebar) ---
st.sidebar.divider()
st.sidebar.subheader("📥 Export Filtered Data")
csv_data = filtered_df.to_csv(index=False).encode('utf-8')
st.sidebar.download_button(
    label="Download Filtered CSV",
    data=csv_data,
    file_name='lanm_filtered_catalog.csv',
    mime='text/csv',
    width='stretch'
)

# --- Main Dashboard ---
metrics_cols = st.columns(6)
metrics_cols[0].metric("Total Designs", len(filtered_df['design_id'].unique()))
metrics_cols[1].metric("Mean pLDDT", f"{filtered_df['plddt'].mean():.2f}")
metrics_cols[2].metric("Mean RMSD", f"{filtered_df['overall_rmsd'].mean():.2f}")
if 'coordination_number' in filtered_df.columns:
    metrics_cols[3].metric("Mean CN", f"{filtered_df['coordination_number'].mean():.1f}")
if 'net_charge' in filtered_df.columns:
    metrics_cols[4].metric("Mean Net Charge", f"{filtered_df['net_charge'].mean():.1f}")
metrics_cols[5].metric("Mean Probability", f"{filtered_df['binding_probability'].mean():.2f}")

# --- Plots ---
st.header("📈 Comparative Analysis")
tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9 = st.tabs([
    "Candidate Evaluation", 
    "Structural Accuracy", 
    "Sequence Analysis", 
    "Loop Distributions", 
    "Pairwise Comparison",
    "Advanced Analytics",
    "Category Comparison",
    "Design Explorer",
    "Cross-Docking Analysis"
])

with tab1:
    col_eval_1, col_eval_2 = st.columns([1, 4])
    with col_eval_1:
        st.write("**Visual Settings**")
        shape_options = {
            "Loop Index": "loop_index",
            "Bidentate Count": "bidentate_count",
            "Metal Ion": "metal_ion",
            "Configuration": "config_index",
            "Motif Match": "motif_match"
        }
        # Filter for existing columns
        shape_options = {k: v for k, v in shape_options.items() if v in filtered_df.columns}
        selected_shape = st.selectbox("Point Shape Mapping", options=list(shape_options.keys()))
        
        selected_color = st.selectbox("Color By", options=["Metal Ion", "Binding Probability", "pLDDT", "Loop Index"])
        
    with col_eval_2:
        fig_eval = get_eval_plot(filtered_df, selected_shape, selected_color, rad_range, rmsd_range)
        st.plotly_chart(fig_eval, width='stretch')

with tab2:
    st.subheader("Structural Quality Distributions")
    fig_violin_rmsd, fig_box_rad, fig_success = get_accuracy_plots(filtered_df)
    
    col_v1, col_v2 = st.columns(2)
    with col_v1:
        st.plotly_chart(fig_violin_rmsd, width='stretch')
    with col_v2:
        st.plotly_chart(fig_box_rad, width='stretch')

    # Tertiary Analysis: Success Rate by Ion
    st.divider()
    st.subheader("Design Success Rate (pLDDT > 0.8 & RMSD < 2.0)")
    st.plotly_chart(fig_success, width='stretch')


with tab3:
    st.subheader("Residue Frequency Analysis")
    
    # Analyze specificity
    seq_to_ions = filtered_df.groupby('loop_sequence')['metal_ion'].unique()
    multi_ion_seqs = seq_to_ions[seq_to_ions.apply(len) > 1]
    
    st.write(f"Total unique loop sequences: **{len(seq_to_ions)}**")
    if not multi_ion_seqs.empty:
        st.warning(f"Found **{len(multi_ion_seqs)}** generalized sequences binding to multiple ions targets.")
    else:
        st.success("100% Sequence Specificity! Every generated sequence is unique to its target ion.")
    
    st.divider()
    col_heat1, col_heat2 = st.columns(2)
    
    with col_heat1:
        st.write("### Global Distribution")
        fig_global = get_residue_heatmap_fig(filtered_df, "Global Residue Frequency", "YlOrRd")
        if fig_global:
            st.plotly_chart(fig_global, width='stretch')
        else:
            st.info("No data available for global heatmap.")
            
    with col_heat2:
        st.write("### Ion-Specific Comparison")
        ion_for_heatmap = st.selectbox("Select Ion for Comparison", options=sorted(filtered_df['metal_ion'].unique().tolist()))
        
        heatmap_data = filtered_df[filtered_df['metal_ion'] == ion_for_heatmap]
        fig_ion = get_residue_heatmap_fig(heatmap_data, f"{ion_for_heatmap} Residue Frequency", "Blues")
        if fig_ion:
            st.plotly_chart(fig_ion, width='stretch')
        else:
            st.info(f"No data available for {ion_for_heatmap}.")
    
    if not filtered_df.empty:
        
        # Tertiary Analysis: pLDDT by Residue Type
        if 'loop_sequence' in filtered_df.columns and 'plddt' in filtered_df.columns:
            st.divider()
            st.subheader("Structural Confidence vs. Loop Composition")
            fig_res = get_res_stability_plot(filtered_df, ion_for_heatmap)
            st.plotly_chart(fig_res, width='stretch')
            st.caption("Shows if certain amino acids consistently correlate with higher structural confidence.")

with tab4:
    fig_hist, fig_len_plddt = get_loop_dist_plots(filtered_df)
    st.plotly_chart(fig_hist, width='stretch')
    
    # Tertiary Analysis: Loop Length vs structural Stability
    if 'plddt' in filtered_df.columns:
        st.divider()
        st.subheader("Loop Length vs. Structural Stability")
        st.plotly_chart(fig_len_plddt, width='stretch')
        st.caption("Investigates if longer or shorter loops tend to be more stable for specific metal ions.")

with tab5:
    st.subheader("Compare Two Molecules")
    col_a, col_b = st.columns(2)
    with col_a:
        design_a = st.selectbox("Select Molecule A", options=filtered_df['design_id'].unique(), key="mol_a")
    with col_b:
        design_b = st.selectbox("Select Molecule B", options=filtered_df['design_id'].unique(), index=min(1, len(filtered_df)-1), key="mol_b")
    
    if design_a and design_b:
        df_a_all = filtered_df[filtered_df['design_id'] == design_a]
        df_b_all = filtered_df[filtered_df['design_id'] == design_b]
        
        col_loop_a, col_loop_b = st.columns(2)
        with col_loop_a:
            loop_a = st.selectbox(f"Select Loop for {design_a}", options=sorted(df_a_all['loop_index'].unique()), key="loop_a")
        with col_loop_b:
            loop_b = st.selectbox(f"Select Loop for {design_b}", options=sorted(df_b_all['loop_index'].unique()), key="loop_b")
            
        df_a = df_a_all[df_a_all['loop_index'] == loop_a].iloc[0]
        df_b = df_b_all[df_b_all['loop_index'] == loop_b].iloc[0]
        
        comparison_df = pd.DataFrame([df_a, df_b])
        st.table(comparison_df[["design_id", "loop_index", "metal_ion", "binding_probability", "plddt", "overall_rmsd", "loop_rmsd", "binding_radius_A", "loop_sequence"]])

        # Radar Chart for Comparison
        st.subheader("Radar Chart Comparison")
        fig_radar = get_radar_chart(filtered_df, design_a, design_b, loop_a, loop_b)
        st.plotly_chart(fig_radar, width='stretch')

        # Population Context
        st.subheader("Population Context")
        fig_context = px.scatter(
            df, x='overall_rmsd', y='plddt', color='metal_ion', opacity=0.3,
            title="Position in Total Population (RMSD vs pLDDT)"
        )
        fig_context.add_trace(go.Scatter(
            x=[df_a['overall_rmsd']], y=[df_a['plddt']], 
            mode='markers+text', text=[f"A: {design_a}"], 
            marker=dict(color='blue', size=15, symbol='star'), name='Molecule A'
        ))
        fig_context.add_trace(go.Scatter(
            x=[df_b['overall_rmsd']], y=[df_b['plddt']], 
            mode='markers+text', text=[f"B: {design_b}"], 
            marker=dict(color='red', size=15, symbol='star'), name='Molecule B'
        ))
        st.plotly_chart(fig_context, width='stretch')
# --- Advanced Analytics Tab ---
with tab6:
    st.header("📈 Advanced Analytics")
    st.write("Deep dive into property correlations and feature influence.")
    
    # Select numeric columns for analysis
    numeric_cols = filtered_df.select_dtypes(include=[np.number]).columns.tolist()
    # Remove some less useful numeric columns if they exist
    for col in ['config_index', 'design_index', 'loop_index']:
        if col in numeric_cols: numeric_cols.remove(col)
    
    fig_corr, fig_pareto = get_analytics_plots(filtered_df, numeric_cols)
    
    col_adv1, col_adv2 = st.columns(2)
    
    with col_adv1:
        st.subheader("Property Correlation Heatmap")
        if fig_corr:
            st.plotly_chart(fig_corr, width='stretch')
        else:
            st.info("Not enough numeric columns for correlation analysis.")
            
    with col_adv2:
        st.subheader("Pareto Front (Confidence vs Accuracy)")
        if 'plddt' in filtered_df.columns and 'overall_rmsd' in filtered_df.columns:
            st.plotly_chart(fig_pareto, width='stretch')
            st.caption("Ideal designs are in the Top-Left: High Confidence, Low Error.")
        else:
            st.info("pLDDT and RMSD required for Pareto analysis.")

    st.divider()
    
    st.subheader("Feature Influence Study")
    influence_col = st.selectbox("Select Categorical Feature", ["metal_ion", "bidentate_count", "loop_index", "is_motif_match"] if "is_motif_match" in filtered_df.columns else ["metal_ion", "bidentate_count", "loop_index"])
    target_metric = st.selectbox("Select Metric to Compare", ["plddt", "overall_rmsd", "loop_rmsd", "binding_radius_A", "coordination_number", "net_charge"], index=1)
    
    if influence_col in filtered_df.columns and target_metric in filtered_df.columns:
        fig_infl = px.box(
            filtered_df, 
            x=influence_col, 
            y=target_metric, 
            color=influence_col,
            notched=True,
            points="all",
            title=f"Influence of {influence_col} on {target_metric}"
        )
        st.plotly_chart(fig_infl, width='stretch')
    
    st.divider()
    
    st.subheader("Multi-Dimensional Property Matrix")
    selected_metrics = st.multiselect("Select Metrics for Matrix", numeric_cols, default=numeric_cols[:4] if len(numeric_cols) >= 4 else numeric_cols)
    if len(selected_metrics) > 1:
        fig_matrix = px.scatter_matrix(
            filtered_df,
            dimensions=selected_metrics,
            color='metal_ion',
            hover_data=['design_id'],
            title="Multi-Metric Relationship Matrix"
        )
        # Scatter matrix can be huge, adjust height
        fig_matrix.update_layout(height=800)
        st.plotly_chart(fig_matrix, width='stretch')

# --- Category Comparison Tab ---
with tab7:
    st.header("⚖️ Group Performance Comparison")
    
    unique_cats = filtered_df['metal_category'].unique()
    
    if len(unique_cats) < 2:
        st.warning("Please select at least two metal categories in the sidebar to perform a comparative analysis.")
    else:
        st.write(f"Comparing performance across **{len(unique_cats)}** selected categories.")
        fig_cat_plddt, fig_cat_rmsd, fig_cat_success, fig_cat_radius = get_category_plots(filtered_df)
        
        col_cat1, col_cat2 = st.columns(2)
        with col_cat1:
            st.subheader("Structural Confidence (pLDDT)")
            st.plotly_chart(fig_cat_plddt, width='stretch')
        with col_cat2:
            st.subheader("Geometric Accuracy (Loop RMSD)")
            st.plotly_chart(fig_cat_rmsd, width='stretch')
            
        st.divider()
        st.subheader("Global Yield Rates")
        st.plotly_chart(fig_cat_success, width='stretch')
        
        st.divider()
        st.subheader("Binding Site Geometry")
        st.plotly_chart(fig_cat_radius, width='stretch')

# --- Design Explorer Tab ---
with tab8:
    st.subheader("Design Inspection & 3D Preview")
    selected_design = st.selectbox("Select a Design to Inspect", options=filtered_df['design_id'].unique())
    
    if selected_design:
        design_loops = filtered_df[filtered_df['design_id'] == selected_design]
        design_info = design_loops.iloc[0]
        
        st.subheader("Global Metrics")
        met1, met2, met3 = st.columns(3)
        met1.write(f"**Ion:** {design_info['metal_ion']}")
        met2.write(f"**pLDDT:** {design_info['plddt']:.4f}")
        met3.write(f"**Overall RMSD:** {design_info['overall_rmsd']:.4f}")
        
        st.subheader("Loop Summary")
        loop_table = design_loops[["loop_index", "binding_probability", "loop_rmsd", "binding_radius_A", "loop_length", "loop_sequence"]]
        st.dataframe(loop_table, width='stretch')
        
        st.divider()
        
        col_l, col_r = st.columns([1, 2])
        
        with col_l:
            st.subheader("Detailed Loop Info")
            selected_loop_idx = st.selectbox("Select Loop to Inspect", options=sorted(design_loops['loop_index'].unique()))
            loop_info = design_loops[design_loops['loop_index'] == selected_loop_idx].iloc[0]
            
            st.write(f"**Binding Probability:** {loop_info['binding_probability']:.3f}")
            if 'coordination_number' in loop_info:
                st.write(f"**Coordination Number:** {loop_info['coordination_number']}")
            if 'net_charge' in loop_info:
                st.write(f"**Net Charge:** {loop_info['net_charge']}")
            if 'bidentate_count' in loop_info:
                st.write(f"**Bidentate Count:** {loop_info['bidentate_count']}")
            
            st.code(loop_info['loop_sequence'], language=None)
            
        with col_r:
            st.subheader("3D Preview")
            cif_file = f"{base_path}/{design_info['metal_ion']}/rf3/{selected_design}_refolded.cif"
            cif_data = get_cif_data(cif_file)
            
            if cif_data:
                st.info(f"Structure: `{os.path.basename(cif_file)}`")
                st.download_button(
                    label="Download Structure (.cif)",
                    data=cif_data,
                    file_name=os.path.basename(cif_file),
                    mime="chemical/x-cif"
                )
                
                view = py3Dmol.view(width=800, height=500)
                view.addModel(cif_data, "cif")
                
                # Adaptive Coloring (Blue=Good, Red=Bad)
                b_min, b_max, invert = get_b_metrics(cif_file)
                
                # If the range is tiny, fall back to a standard 0-1 scale to avoid flickering
                if b_max - b_min < 0.0001:
                    b_min, b_max = 0.0, 1.0
                    invert = False
                    
                # Reverse the gradient by swapping min and max if needed
                gradient_settings = {'prop': 'b', 'gradient': 'roygb', 'min': b_max if invert else b_min, 'max': b_min if invert else b_max}
                view.setStyle({'cartoon': {'colorscheme': gradient_settings}})
                
                # Show metal
                view.addStyle({'resn': design_info['metal_ion']}, {'sphere': {'radius': 1.5, 'color': 'silver'}})
                
                # Highlight loop residues - we need to find them in the structure.
                # Usually the loop is what we redesigned.
                # We can't easily know the residue IDs from the catalog alone without indexing.
                # But the catalog has loop_sequence.
                
                view.zoomTo()
                showmol(view, height=500, width=800)
            else:
                st.warning(f"Structure file not found at {cif_file}")

# --- Cross-Docking Analysis Tab ---
with tab9:
    cross_dock_df = load_cross_docking_data(base_path)
    if not cross_dock_df.empty and selected_ions:
        cross_dock_df = cross_dock_df[
            cross_dock_df['original_ion'].isin(selected_ions) & 
            cross_dock_df['new_ion'].isin(selected_ions)
        ]
        
    if not cross_dock_df.empty:
        st.subheader("🔄 Cross-Docking Analysis")
        st.write("Analysis of how designs for one ion perform when modeled with a different ion.")
        
        st.subheader("Cross-Docking Metrics")
        
        # Base Metrics
        m1, m2, m3 = st.columns(3)
        m1.metric("Total Cross-Docking Runs", len(cross_dock_df))
        m2.metric("Mean RMSD Deviation", f"{cross_dock_df['structural_deviation_rmsd'].mean():.2f}")
        m3.metric("Mean Swap pLDDT", f"{cross_dock_df['plddt'].mean():.2f}")
        
        # Conditional New Metrics
        extra_metrics = [m for m in ['coordination_number', 'net_charge', 'bidentate_count'] if m in cross_dock_df.columns]
        if extra_metrics:
            extra_cols = st.columns(len(extra_metrics))
            for i, col_name in enumerate(extra_metrics):
                title = col_name.replace("_", " ").title()
                val = cross_dock_df[col_name].mean()
                extra_cols[i].metric(f"Mean {title}", f"{val:.2f}")
        
        col_cd1, col_cd2 = st.columns(2)
        
        with col_cd1:
            st.subheader("Structural Deviation by Original Ion")
            fig_cd_rmsd = px.box(cross_dock_df, x='original_ion', y='structural_deviation_rmsd', color='original_ion', 
                               points="all", title="Deviation When Swapping Metal Ion")
            st.plotly_chart(fig_cd_rmsd, width='stretch')
            
        with col_cd2:
            st.subheader("Confidence (pLDDT) After Swap")
            fig_cd_plddt = px.box(cross_dock_df, x='new_ion', y='plddt', color='new_ion', 
                                points="all", title="pLDDT for Swapped Targets")
            st.plotly_chart(fig_cd_plddt, width='stretch')
            
        # Optional Property Distributions
        if extra_metrics:
            st.divider()
            st.subheader("Chemical Property Distributions Post-Swap")
            prop_cols = st.columns(len(extra_metrics))
            
            for i, col_name in enumerate(extra_metrics):
                with prop_cols[i]:
                    title = col_name.replace("_", " ").title()
                    fig_prop = px.box(cross_dock_df, x='new_ion', y=col_name, color='new_ion',
                                      points="all", title=f"{title} by Target Ion")
                    st.plotly_chart(fig_prop, width='stretch')
            
        st.divider()
        st.subheader("Cross-Docking Matrix (Average Deviation)")
        
        # Create a pivot table for the heatmap
        pivot_df = cross_dock_df.pivot_table(values='structural_deviation_rmsd', index='original_ion', columns='new_ion', aggfunc='mean')
        
        fig_cd_heat = px.imshow(pivot_df, text_auto=".2f", aspect="auto", 
                              labels=dict(x="new_ion", y="original_ion", color="Mean RMSD"),
                              title="Average Structural Deviation (RMSD) Heatmap",
                              color_continuous_scale="Viridis_r")
        st.plotly_chart(fig_cd_heat, width='stretch')

        st.divider()
        st.subheader("Cross-Docking Explorer")
        st.dataframe(cross_dock_df, width='stretch')
        
    else:
        st.info("No cross-docking data available for the selected target ions, or the cross-docking workflow has not been run.")
