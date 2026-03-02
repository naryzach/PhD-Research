import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os
import numpy as np
from stmol import showmol
import py3Dmol

# Set page config
st.set_page_config(page_title="Lanm Output Dashboard", layout="wide")

st.title("🧬 Lanm Output Dashboard")
st.markdown("Explore metal-binding protein designs with advanced structural and sequence analysis.")

# --- Data Loading ---
@st.cache_data
def load_data():
    base_path = "/home/ryangustafson/Documents/GitHubProj/PhD-Research/Local/lanm_output"
    catalog_path = os.path.join(base_path, "global_sequence_catalog.csv")
    full_seq_path = os.path.join(base_path, "full_sequences_log.csv")
    
    df = pd.read_csv(catalog_path) if os.path.exists(catalog_path) else pd.DataFrame()
    full_seq_df = pd.read_csv(full_seq_path) if os.path.exists(full_seq_path) else pd.DataFrame()
    
    # Check for ions presence even if not in catalog
    ions = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d)) and d.isupper() and len(d) <= 2]
    
    return df, full_seq_df, ions, base_path

def get_b_factor_range(cif_data):
    """Parses B-factors from CIF focusing on CA atoms (amino acid scores)."""
    lines = cif_data.split('\n')
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
                # Filter for CA atoms to get "Amino Acid" level scores
                atom_name = parts[col_atom].strip('"').strip("'")
                if atom_name == "CA":
                    try:
                        b_factors.append(float(parts[col_b]))
                    except ValueError: continue
        elif not header_section and line.startswith('#'):
            break
            
    if b_factors:
        return min(b_factors), max(b_factors)
    return 0.0, 1.0 # Fallback

df, full_seq_df, available_ions, base_path = load_data()

if df.empty:
    st.warning("⚠️ No data found in global_sequence_catalog.csv. Please ensure the pipeline has generated results.")
    st.stop()

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

# Ion selection
selected_ions = st.sidebar.multiselect("Select Target Ions", options=sorted(available_ions), default=sorted(df['metal_ion'].unique()) if not df.empty else [])

# Filter dataframe
filtered_df = df.copy()

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

if selected_ions:
    filtered_df = filtered_df[filtered_df['metal_ion'].isin(selected_ions)]

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

# --- Main Dashboard ---
metrics_cols = st.columns(5)
metrics_cols[0].metric("Total Designs", len(filtered_df['design_id'].unique()))
metrics_cols[1].metric("Mean pLDDT", f"{filtered_df['plddt'].mean():.2f}")
metrics_cols[2].metric("Mean RMSD", f"{filtered_df['overall_rmsd'].mean():.2f}")
if 'coordination_number' in filtered_df.columns:
    metrics_cols[3].metric("Mean CN", f"{filtered_df['coordination_number'].mean():.1f}")
if 'net_charge' in filtered_df.columns:
    metrics_cols[4].metric("Mean Net Charge", f"{filtered_df['net_charge'].mean():.1f}")

# --- Plots ---
st.header("📈 Comparative Analysis")
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "Candidate Evaluation", 
    "Structural Accuracy", 
    "Sequence Analysis", 
    "Loop Distributions", 
    "Pairwise Comparison"
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
        
    with col_eval_2:
        fig_eval = px.scatter(
            filtered_df, 
            x='loop_rmsd', 
            y='binding_radius_A', 
            color='metal_ion', 
            size='plddt' if 'plddt' in filtered_df.columns else None,
            symbol=shape_options[selected_shape],
            hover_data=['design_id', 'loop_sequence', 'coordination_number', 'net_charge'] if 'coordination_number' in filtered_df.columns else ['design_id', 'loop_sequence'],
            labels={'loop_rmsd': 'Individual Loop RMSD (Å)', 'binding_radius_A': 'Binding Radius (Å)'},
            title="Candidate Evaluation (Thresholds: RMSD < 1.5, Radius 2.3-2.6)"
        )
        
        # Force axis range to match filters
        if rad_range:
            fig_eval.update_yaxes(range=[rad_range[0]-0.1, rad_range[1]+0.1])
        if rmsd_range:
            # We use loop_rmsd on X, so let's allow it to show a bit more
            fig_eval.update_xaxes(range=[0, max(1.5, filtered_df['loop_rmsd'].max() + 0.2)])
    
    # Add horizontal region for ideal radius
    fig_eval.add_hrect(y0=2.3, y1=2.6, line_width=0, fillcolor="green", opacity=0.1)
    # Add vertical line for RMSD threshold
    fig_eval.add_vline(x=1.5, line_dash="dash", line_color="red", opacity=0.6)
    
    st.plotly_chart(fig_eval, width='stretch')

with tab2:
    st.subheader("Structural Quality Distributions")
    col_v1, col_v2 = st.columns(2)
    with col_v1:
        fig_violin_rmsd = px.violin(filtered_df, x='metal_ion', y='loop_rmsd', color='metal_ion', box=True, points="all", title="Loop RMSD Distribution")
        st.plotly_chart(fig_violin_rmsd, width='stretch')
    with col_v2:
        fig_box_rad = px.box(filtered_df, x='metal_ion', y='binding_radius_A', color='metal_ion', points="all", title="Binding Radius Distribution")
        fig_box_rad.add_hrect(y0=2.3, y1=2.6, line_width=0, fillcolor="green", opacity=0.1)
        st.plotly_chart(fig_box_rad, width='stretch')

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
    
    # Heatmap Generation
    ion_for_heatmap = st.selectbox("Select Ion for Residue Heatmap", options=["Global"] + sorted(filtered_df['metal_ion'].unique().tolist()))
    
    if ion_for_heatmap == "Global":
        heatmap_data = filtered_df
    else:
        heatmap_data = filtered_df[filtered_df['metal_ion'] == ion_for_heatmap]
    
    if not heatmap_data.empty:
        max_len = heatmap_data['loop_length'].max()
        aa_list = list("ACDEFGHIKLMNPQRSTVWY")
        freq_matrix = pd.DataFrame(0.0, index=aa_list, columns=range(max_len))
        
        for _, row in heatmap_data.iterrows():
            seq = row['loop_sequence']
            for i, aa in enumerate(seq):
                if aa in freq_matrix.index and i < max_len:
                    freq_matrix.loc[aa, i] += 1
        
        # Normalize
        freq_matrix = freq_matrix.div(freq_matrix.sum(axis=0), axis=1).fillna(0)
        
        fig_heatmap = px.imshow(
            freq_matrix,
            labels=dict(x="Loop Position", y="Amino Acid", color="Frequency"),
            x=[f"Pos {i}" for i in range(max_len)],
            y=aa_list,
            color_continuous_scale="YlOrRd" if ion_for_heatmap == "Global" else "Blues",
            title=f"Residue Frequency Heatmap ({ion_for_heatmap})"
        )
        st.plotly_chart(fig_heatmap, width='stretch')

with tab4:
    fig_hist = px.histogram(
        filtered_df, 
        x="loop_length", 
        color="metal_ion", 
        barmode="group",
        title="Loop Length Distributions",
        labels={"loop_length": "Loop Length", "count": "Frequency"}
    )
    st.plotly_chart(fig_hist, width='stretch')

with tab5:
    st.subheader("Compare Two Molecules")
    col_a, col_b = st.columns(2)
    with col_a:
        design_a = st.selectbox("Select Molecule A", options=filtered_df['design_id'].unique(), key="mol_a")
    with col_b:
        design_b = st.selectbox("Select Molecule B", options=filtered_df['design_id'].unique(), index=min(1, len(filtered_df)-1), key="mol_b")
    
    if design_a and design_b:
        df_a = filtered_df[filtered_df['design_id'] == design_a].iloc[0]
        df_b = filtered_df[filtered_df['design_id'] == design_b].iloc[0]
        
        comparison_df = pd.DataFrame([df_a, df_b])
        st.table(comparison_df[["design_id", "metal_ion", "plddt", "overall_rmsd", "loop_rmsd", "binding_radius_A", "loop_sequence"]])

        # Radar Chart for Comparison
        st.subheader("Radar Chart Comparison")
        metrics = ['plddt', 'overall_rmsd', 'loop_rmsd', 'binding_radius_A', 'coordination_number']
        metrics = [m for m in metrics if m in filtered_df.columns]
        
        # Normalize metrics for radar chart (0-1 scale)
        def normalize(val, col):
            v_min = df[col].min()
            v_max = df[col].max()
            if v_max == v_min: return 0.5
            return (val - v_min) / (v_max - v_min)

        fig_radar = go.Figure()

        for d_id, color in [(design_a, 'blue'), (design_b, 'red')]:
            row = filtered_df[filtered_df['design_id'] == d_id].iloc[0]
            r_values = [normalize(row[m], m) for m in metrics]
            # Close the circle
            r_values.append(r_values[0])
            theta_labels = metrics + [metrics[0]]
            
            fig_radar.add_trace(go.Scatterpolar(
                r=r_values,
                theta=theta_labels,
                fill='toself',
                name=d_id,
                line_color=color
            ))

        fig_radar.update_layout(
            polar=dict(radialaxis=dict(visible=True, range=[0, 1])),
            showlegend=True,
            title="Relative Metric Strength (Normalized to Global Min/Max)"
        )
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

# --- Design Explorer ---
st.header("🔍 Design Explorer")
selected_design = st.selectbox("Select a Design to Inspect", options=filtered_df['design_id'].unique())

if selected_design:
    design_info = filtered_df[filtered_df['design_id'] == selected_design].iloc[0]
    
    col_l, col_r = st.columns([1, 2])
    
    with col_l:
        st.subheader("Metrics")
        st.write(f"**Ion:** {design_info['metal_ion']}")
        st.write(f"**pLDDT:** {design_info['plddt']:.4f}")
        st.write(f"**Overall RMSD:** {design_info['overall_rmsd']:.4f}")
        st.write(f"**Loop RMSD:** {design_info['loop_rmsd']:.4f}")
        st.write(f"**Binding Radius:** {design_info['binding_radius_A']:.4f} Å")
        
        if 'coordination_number' in design_info:
            st.write(f"**Coordination Number:** {design_info['coordination_number']}")
        if 'net_charge' in design_info:
            st.write(f"**Net Charge:** {design_info['net_charge']}")
        if 'bidentate_count' in design_info:
            st.write(f"**Bidentate Count:** {design_info['bidentate_count']}")
            
        st.write(f"**Loop Length:** {design_info['loop_length']}")
        st.code(design_info['loop_sequence'], language=None)
        
    with col_r:
        st.subheader("3D Preview")
        cif_file = f"{base_path}/{design_info['metal_ion']}/rf3/{selected_design}_refolded.cif"
        if os.path.exists(cif_file):
            st.info(f"Structure: `{os.path.basename(cif_file)}`")
            
            with open(cif_file, "r") as f:
                cif_data = f.read()
            
            view = py3Dmol.view(width=800, height=500)
            view.addModel(cif_data, "cif")
            
            # Relative pLDDT coloring (Blue=Local Max, Red=Local Min)
            # This ensures contrast even if absolute pLDDT is low.
            b_min, b_max = get_b_factor_range(cif_data)
            
            # If the range is tiny, fall back to a standard 0-1 scale to avoid flickering
            if b_max - b_min < 0.01:
                is_0_1 = b_max <= 1.0
                b_min, b_max = (0.0, 1.0) if is_0_1 else (0, 100)
            
            view.setStyle({'cartoon': {'colorscheme': {'prop': 'b', 'gradient': 'roygb', 'min': b_min, 'max': b_max}}})
            
            # Show metal
            view.addStyle({'resn': design_info['metal_ion']}, {'sphere': {'radius': 1.5, 'color': 'red'}})
            
            # Highlight loop residues - we need to find them in the structure.
            # Usually the loop is what we redesigned.
            # We can't easily know the residue IDs from the catalog alone without indexing.
            # But the catalog has loop_sequence.
            
            view.zoomTo()
            showmol(view, height=500, width=800)
        else:
            st.warning(f"Structure file not found at {cif_file}")

# --- Export ---
st.header("📥 Export Data")
csv = filtered_df.to_csv(index=False).encode('utf-8')
st.download_button(
    label="Download Filtered Catalog as CSV",
    data=csv,
    file_name='lanm_filtered_catalog.csv',
    mime='text/csv',
)
