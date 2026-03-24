import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os
import numpy as np

st.set_page_config(page_title="TIMP3 PPI Dashboard", layout="wide")

st.title("🧬 TIMP3 Protein-Protein Interaction Dashboard")
st.markdown("Explore generated TIMP3 variants and their AlphaFold prediction metrics against various MMP/ADAM targets.")

@st.cache_data
def load_data():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    potential_paths = [
        os.path.join(script_dir, "..", "Local", "TIMP-Dashboard_output", "global_sequence_catalog.csv"),
        os.path.join(script_dir, "global_sequence_catalog.csv")
    ]
    for p in potential_paths:
        if os.path.exists(p):
            return pd.read_csv(p)
    return pd.DataFrame()

df = load_data()

if df.empty:
    st.warning("⚠️ No data found in global_sequence_catalog.csv. Please run the generation pipeline first.")
    st.stop()

# --- Sidebar Filters ---
st.sidebar.header("Filters")

available_targets = sorted(df['target'].unique().tolist())
selected_targets = st.sidebar.multiselect("Select Targets", options=available_targets, default=available_targets)

available_combos = sorted(df['loop_combo'].unique().tolist())
selected_combos = st.sidebar.multiselect("Select Loop Combinations", options=available_combos, default=available_combos)

filtered_df = df[df['target'].isin(selected_targets) & df['loop_combo'].isin(selected_combos)].copy()

st.sidebar.subheader("Quality Metrics")
def get_range_filter(df, col, label, key):
    if col in df.columns:
        m_min, m_max = float(df[col].min()), float(df[col].max())
        if m_min == m_max: return (m_min, m_max)
        if key not in st.session_state:
            st.session_state[key] = (m_min, m_max)
        else:
            curr_val = st.session_state[key]
            st.session_state[key] = (max(m_min, min(m_max, float(curr_val[0]))), max(m_min, min(m_max, float(curr_val[1]))))
        return st.sidebar.slider(label, m_min, m_max, key=key)
    return None

score_range = get_range_filter(df, 'heuristic_score', "Heuristic Score", "score_slider")
plddt_range = get_range_filter(df, 'plddt', "pLDDT", "plddt_slider")
rmsd_range = get_range_filter(df, 'overall_rmsd', "Overall RMSD", "rmsd_slider")
contacts_range = get_range_filter(df, 'contacts', "Contacts", "contacts_slider")
ia_range = get_range_filter(df, 'interface_area', "Interface Area (Å²)", "ia_slider")
clashes_range = get_range_filter(df, 'clashes', "Clashes", "clash_slider")

if score_range: filtered_df = filtered_df[(filtered_df['heuristic_score'] >= score_range[0]) & (filtered_df['heuristic_score'] <= score_range[1])]
if plddt_range: filtered_df = filtered_df[(filtered_df['plddt'] >= plddt_range[0]) & (filtered_df['plddt'] <= plddt_range[1])]
if rmsd_range: filtered_df = filtered_df[(filtered_df['overall_rmsd'] >= rmsd_range[0]) & (filtered_df['overall_rmsd'] <= rmsd_range[1])]
if contacts_range: filtered_df = filtered_df[(filtered_df['contacts'] >= contacts_range[0]) & (filtered_df['contacts'] <= contacts_range[1])]
if ia_range: filtered_df = filtered_df[(filtered_df['interface_area'] >= ia_range[0]) & (filtered_df['interface_area'] <= ia_range[1])]
if clashes_range: filtered_df = filtered_df[(filtered_df['clashes'] >= clashes_range[0]) & (filtered_df['clashes'] <= clashes_range[1])]

# --- Top Metrics ---
m1, m2, m3, m4, m5 = st.columns(5)
m1.metric("Total Designs", len(filtered_df['design_id'].unique()))
if 'plddt' in filtered_df.columns:
    m2.metric("Mean pLDDT", f"{filtered_df['plddt'].mean():.2f}")
if 'interface_area' in filtered_df.columns:
    m3.metric("Mean Interface Area", f"{filtered_df['interface_area'].mean():.1f}")
if 'contacts' in filtered_df.columns:
    m4.metric("Mean Contacts", f"{filtered_df['contacts'].mean():.1f}")
if 'heuristic_score' in filtered_df.columns:
    m5.metric("Mean Score", f"{filtered_df['heuristic_score'].mean():.1f}")

st.header("📈 Comparative Analysis")

tab1, tab2, tab3, tab4 = st.tabs([
    "Candidate Evaluation", 
    "Loop Combination Impact", 
    "Structural Quality", 
    "Data Explorer"
])

with tab1:
    col1, col2 = st.columns([1, 4])
    with col1:
        color_col = st.selectbox("Color By", options=["target", "loop_combo", "heuristic_score", "plddt", "overall_rmsd"])
        size_col = st.selectbox("Size By", options=["None", "contacts", "interface_area", "heuristic_score"])
    
    with col2:
        hover_data = ["design_id", "target", "loop_combo", "plddt", "contacts", "clashes", "interface_area"]
        fig1 = px.scatter(
            filtered_df, 
            x="interface_area", 
            y="contacts", 
            color=color_col,
            size=size_col if size_col != "None" else None,
            hover_data=hover_data,
            title="Interface Area vs. Contact Count",
            labels={"interface_area": "Interface Area (Å²)", "contacts": "Contact Count"}
        )
        st.plotly_chart(fig1, use_container_width=True)

with tab2:
    st.subheader("Performance by Loop Combination")
    y_metric = st.selectbox("Metric", options=["heuristic_score", "interface_area", "contacts", "plddt", "overall_rmsd"], key="loop_metric")
    fig2 = px.box(filtered_df, x="loop_combo", y=y_metric, color="target", points="all", title=f"{y_metric} Across Loop Combinations")
    fig2.update_layout(xaxis={'categoryorder':'total descending'})
    st.plotly_chart(fig2, use_container_width=True)

with tab3:
    st.subheader("Structural Quality Distributions")
    c1, c2 = st.columns(2)
    with c1:
        fig3 = px.violin(filtered_df, x="target", y="plddt", color="target", points="all", title="pLDDT Distribution")
        st.plotly_chart(fig3, use_container_width=True)
    with c2:
        fig4 = px.violin(filtered_df, x="target", y="overall_rmsd", color="target", points="all", title="RMSD Distribution")
        st.plotly_chart(fig4, use_container_width=True)

with tab4:
    st.subheader("Raw Data Table")
    st.dataframe(filtered_df)
