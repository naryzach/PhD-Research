import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os
import numpy as np
import json
import py3Dmol
import streamlit.components.v1 as components

st.set_page_config(page_title="TIMP3 PPI Dashboard", layout="wide")

st.title("🧬 TIMP3 Protein-Protein Interaction Dashboard")
st.markdown("Explore generated TIMP3 variants and their RF3 prediction metrics against various MMP/ADAM targets.")

@st.cache_data
def load_data():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(script_dir, "..", "Local", "TIMP-Dashboard_output")
    
    df = pd.DataFrame()
    lead_df = pd.DataFrame()
    xd_df = pd.DataFrame()
    
    p_adv = os.path.join(out_dir, "advanced_metrics.csv")
    p_cat = os.path.join(out_dir, "global_sequence_catalog.csv")
    
    if os.path.exists(p_adv):
        df = pd.read_csv(p_adv)
    elif os.path.exists(p_cat):
        df = pd.read_csv(p_cat)
    elif os.path.exists(os.path.join(script_dir, "global_sequence_catalog.csv")):
        df = pd.read_csv(os.path.join(script_dir, "global_sequence_catalog.csv"))
        
    p_lead = os.path.join(out_dir, "consensus_ranking.csv")
    if os.path.exists(p_lead):
        lead_df = pd.read_csv(p_lead)
        
    p_xd = os.path.join(out_dir, "cross_docking_metrics.csv")
    if os.path.exists(p_xd):
        xd_df = pd.read_csv(p_xd)
        
    return df, lead_df, xd_df

df, lead_df, xd_df = load_data()

if df.empty:
    st.warning("⚠️ No data found in global_sequence_catalog.csv. Please run the generation pipeline first.")
    st.stop()

df['target_path'] = df['target']
df['target'] = df['target'].astype(str).str.replace('TIMP3_vs_', '').str.replace('_AF', '')
if 'heuristic_score' in df.columns and 'probability_of_binding_score' not in df.columns:
    df.rename(columns={'heuristic_score': 'probability_of_binding_score'}, inplace=True)

if not lead_df.empty:
    lead_df['target'] = lead_df['target'].astype(str).str.replace('TIMP3_vs_', '').str.replace('_AF', '')
    if 'heuristic_score' in lead_df.columns and 'probability_of_binding_score' not in lead_df.columns:
        lead_df.rename(columns={'heuristic_score': 'probability_of_binding_score'}, inplace=True)

def get_seq_str(row):
    loops = str(row['loop_combo']).split('_')
    seqs = []
    for l in loops:
        val = row.get(f"loop_{l}_seq")
        if pd.notna(val) and str(val).upper() != "MISSING":
            seqs.append(str(val))
    return "-".join(seqs) if seqs else row['design_id'].split('_')[-1]

df['sequence'] = df.apply(get_seq_str, axis=1)
if not lead_df.empty:
    lead_df['sequence'] = lead_df.apply(get_seq_str, axis=1)

if not xd_df.empty:
    # Map sequence from main df to xd_df based on design_id
    seq_map = df.set_index('design_id')['sequence'].to_dict()
    xd_df['sequence'] = xd_df['design_id'].map(seq_map)
    # Ensure all targets are cleaned up for consistency
    xd_df['intended_target'] = xd_df['intended_target'].astype(str).str.replace('TIMP3_vs_', '').str.replace('_AF', '')
    xd_df['folded_target'] = xd_df['folded_target'].astype(str).str.replace('TIMP3_vs_', '').str.replace('_AF', '')
    # Add a 'target' column for compatibility with general heatmap logic (using folded_target as the effective target)
    if 'target' not in xd_df.columns:
        xd_df['target'] = xd_df['folded_target']
    # Ensure loop_combo is available in xd_df by mapping from df
    if 'loop_combo' not in xd_df.columns:
        combo_map = df.set_index('design_id')['loop_combo'].to_dict()
        xd_df['loop_combo'] = xd_df['design_id'].map(combo_map)

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
available_all_metrics = [c for c in ["probability_of_binding_score", "plddt", "TIMP_pLDDT", "Target_pLDDT", "TIMP_pTM", "Target_pTM", "mean_loop_plddt", "overall_rmsd", "contacts", "bsa", "mean_loop_pae", "iptm", "h_bonds", "interface_area"] if c in df.columns]

score_range = get_range_filter(df, 'probability_of_binding_score', "Probability of Binding Score", "score_slider")
plddt_range = get_range_filter(df, 'plddt', "Overall pLDDT", "plddt_slider")
loop_plddt_range = get_range_filter(df, 'mean_loop_plddt', "Loop pLDDT", "loop_plddt_slider")
rmsd_range = get_range_filter(df, 'overall_rmsd', "Overall RMSD", "rmsd_slider")
contacts_range = get_range_filter(df, 'contacts', "Contacts", "contacts_slider")
ia_range = get_range_filter(df, 'bsa', "BSA (Å²)", "bsa_slider") if 'bsa' in df.columns else get_range_filter(df, 'interface_area', "Interface Area (Å²)", "ia_slider")
clashes_range = get_range_filter(df, 'clashes', "Clashes", "clash_slider")

if score_range: filtered_df = filtered_df[(filtered_df['probability_of_binding_score'] >= score_range[0]) & (filtered_df['probability_of_binding_score'] <= score_range[1])]
if plddt_range: filtered_df = filtered_df[(filtered_df['plddt'] >= plddt_range[0]) & (filtered_df['plddt'] <= plddt_range[1])]
if loop_plddt_range: filtered_df = filtered_df[(filtered_df['mean_loop_plddt'] >= loop_plddt_range[0]) & (filtered_df['mean_loop_plddt'] <= loop_plddt_range[1])]
if rmsd_range: filtered_df = filtered_df[(filtered_df['overall_rmsd'] >= rmsd_range[0]) & (filtered_df['overall_rmsd'] <= rmsd_range[1])]
if contacts_range: filtered_df = filtered_df[(filtered_df['contacts'] >= contacts_range[0]) & (filtered_df['contacts'] <= contacts_range[1])]
if ia_range and 'bsa' in df.columns: filtered_df = filtered_df[(filtered_df['bsa'] >= ia_range[0]) & (filtered_df['bsa'] <= ia_range[1])]
elif ia_range: filtered_df = filtered_df[(filtered_df['interface_area'] >= ia_range[0]) & (filtered_df['interface_area'] <= ia_range[1])]
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
if 'probability_of_binding_score' in filtered_df.columns:
    m5.metric("Mean Score", f"{filtered_df['probability_of_binding_score'].mean():.1f}")

st.header("📈 Comparative Analysis")

tab1, tab_adv, tab_comp, tab_heat, tab_leader, tab_spec, tab_viewer, tab_data = st.tabs([
    "Candidate Evaluation", 
    "Structural Metrics",
    "Comparative Metrics",
    "Heatmaps",
    "Leaderboard",
    "Specificity Analysis",
    "Variant Viewer",
    "Data Explorer"
])

with tab1:
    col1, col2 = st.columns([1, 4])
    with col1:
        color_col = st.selectbox("Color By", options=["target", "loop_combo", "probability_of_binding_score", "plddt", "overall_rmsd"])
        size_col = st.selectbox("Size By", options=["None", "contacts", "interface_area", "probability_of_binding_score"])
    
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
        st.plotly_chart(fig1, width='stretch')

with tab_adv:
    st.subheader("Structural Metrics")
    if 'probability_of_binding_score' not in filtered_df.columns:
        st.info("Advanced metrics not available. Run analyze_results.py first.")
    else:
        c1, c2 = st.columns(2)
        with c1:
            adv_m1 = st.selectbox("Heatmap Metric", ["iptm", "bsa", "mean_loop_plddt", "mean_loop_pae", "probability_of_binding_score", "contacts", "overall_rmsd"])
            fig_ip = px.density_heatmap(filtered_df, x="target", y="loop_combo", z=adv_m1, histfunc="avg", title=f"{adv_m1} Heatmap")
            st.plotly_chart(fig_ip, width='stretch', key='adv_heatmap')
        with c2:
            adv_m2 = st.selectbox("Violin Metric", ["h_bonds", "salt_bridges", "contacts", "interface_area", "plddt", "iptm", "probability_of_binding_score"])
            fig_hbond = px.violin(filtered_df, x="target", y=adv_m2, color="target", points="all", title=f"{adv_m2} Distribution")
            st.plotly_chart(fig_hbond, width='stretch', key='adv_violin')
            
        c3, c4 = st.columns(2)
        with c3:
            adv_m3_y = st.selectbox("Scatter Y", available_all_metrics, index=available_all_metrics.index("mean_loop_plddt") if "mean_loop_plddt" in available_all_metrics else 0)
            adv_m3_x = st.selectbox("Scatter X", available_all_metrics, index=available_all_metrics.index("plddt") if "plddt" in available_all_metrics else 0)
            fig_plddt = px.scatter(filtered_df, x=adv_m3_x, y=adv_m3_y, color="target", hover_data=["design_id"], title=f"{adv_m3_y} vs {adv_m3_x}")
            st.plotly_chart(fig_plddt, width='stretch', key='adv_scatter')
        with c4:
            adv_m4 = st.selectbox("Boxplot Metric", ["mean_loop_pae", "plddt", "iptm", "bsa", "probability_of_binding_score", "overall_rmsd", "contacts"])
            fig_pae = px.box(filtered_df, x="target", y=adv_m4, color="loop_combo", title=f"{adv_m4} vs Target")
            if adv_m4 in ["mean_loop_pae", "overall_rmsd", "clashes"]:
                fig_pae.update_layout(xaxis={'categoryorder':'total descending'}, yaxis={'autorange': 'reversed'}) # Lower is better
            st.plotly_chart(fig_pae, width='stretch', key='adv_box')

with tab_leader:
    st.subheader("🏆 Leaderboard")
    lead_c1, lead_c2, lead_c3 = st.columns(3)
    
    with lead_c1: lead_chart = st.selectbox("Leaderboard Format", ["Consensus Best Binders (Table)", "Variant Leaderboard (Bar)", "Raw Ranking Table"])
    with lead_c2: lead_metric = st.selectbox("Primary Ranking Metric", available_all_metrics, index=0, key="lead_metric")
    with lead_c3: lead_trans = st.selectbox("Score Transformation", ["Raw Scores", "T-Score (Overall)", "Local T-Score (Within Targets)"], key="lead_trans")

    if lead_chart == "Consensus Best Binders (Table)":
        if lead_df.empty:
            st.info("Leaderboard not generated yet. Run analyze_results.py.")
        else:
            filtered_lead = lead_df[lead_df['target'].isin(selected_targets) & lead_df['loop_combo'].isin(selected_combos)]
            cols_to_show = ["sequence", "target", "loop_combo", "probability_of_binding_score", "iptm", "bsa", "mean_loop_pae", "plddt", "mean_loop_plddt", "h_bonds", "salt_bridges", "overall_rmsd"]
            cols_to_show = [c for c in cols_to_show if c in filtered_lead.columns]
            st.dataframe(filtered_lead[cols_to_show])
    else:
        df_lead = filtered_df.copy()
        if lead_trans == "T-Score (Overall)":
            mean_val = df_lead[lead_metric].mean()
            std_val = df_lead[lead_metric].std()
            sem_val = std_val / np.sqrt(len(df_lead)) if len(df_lead) > 1 and std_val > 0 else 1
            df_lead["plot_val"] = (df_lead[lead_metric] - mean_val) / sem_val
        elif lead_trans == "Local T-Score (Within Targets)":
            df_lead["plot_val"] = 0.0
            for tgt, group in df_lead.groupby("target"):
                g_mean, g_std = group[lead_metric].mean(), group[lead_metric].std()
                g_sem = g_std / np.sqrt(len(group)) if len(group) > 1 and g_std > 0 else 1
                df_lead.loc[group.index, "plot_val"] = (group[lead_metric] - g_mean) / g_sem
        else:
            df_lead["plot_val"] = df_lead[lead_metric]
            
        if lead_chart == "Variant Leaderboard (Bar)":
            df_lead = df_lead.sort_values(by="plot_val", ascending=False).head(50)
            fig_lead = px.bar(df_lead, x="sequence", y="plot_val", color="target", hover_data=["loop_combo", "design_id"], title=f"Top 50 Variants by {lead_trans} ({lead_metric})")
            st.plotly_chart(fig_lead, width='stretch', key='leader_bar_chart')
            
        elif lead_chart == "Raw Ranking Table":
            st.dataframe(df_lead[["sequence", "target", "loop_combo", lead_metric, "plot_val", "design_id"]].sort_values(by="plot_val", ascending=False))

with tab_spec:
    st.header("🎯 Target-Specificity Analysis (Cross-Docking)")
    st.markdown("This analysis evaluates how well variants designed for one target bind to other 'alien' targets.")
    
    xd_csv = os.path.join("../Local/TIMP-Dashboard_output", "cross_docking_metrics.csv")
    if not os.path.exists(xd_csv):
        st.info("Cross-docking analysis results not found. Please run `cross_docking_analysis.py` to generate this data.")
    else:
        df_xd = pd.read_csv(xd_csv)
        
        # heatmap metric selection
        xd_metric = st.selectbox("Specificity Metric", ["plddt", "iptm", "ptm"], key="xd_metric_sel")
        
        # Pivot: Intended (Y) vs Folded (X)
        xd_pivot = df_xd.pivot_table(index="intended_target", columns="folded_target", values=xd_metric, aggfunc="mean")
        
        import plotly.express as px
        fig_xd = px.imshow(xd_pivot, 
                          labels=dict(x="Folded Against Target", y="Intended Target", color=xd_metric.upper()),
                          color_continuous_scale="Viridis",
                          text_auto=".3f")
        st.plotly_chart(fig_xd, width='stretch')
        
        # Specificity Leaderboard
        st.subheader("Binder Specificity Ranking")
        # Calculate specificity score: Score_intended - Mean(Score_alien)
        spec_results = []
        for d_id, group in df_xd.groupby("design_id"):
            intended = group[group['intended_target'] == group['folded_target']]
            alien = group[group['intended_target'] != group['folded_target']]
            
            if not intended.empty and not alien.empty:
                i_val = intended[xd_metric].values[0]
                a_avg = alien[xd_metric].mean()
                spec_score = i_val - a_avg
                spec_results.append({
                    "design_id": d_id,
                    "target": intended['intended_target'].values[0],
                    f"intended_{xd_metric}": i_val,
                    f"alien_avg_{xd_metric}": a_avg,
                    "specificity_score": spec_score
                })
        
        if spec_results:
            df_spec = pd.DataFrame(spec_results).sort_values("specificity_score", ascending=False)
            st.dataframe(df_spec, width='stretch')
            
with tab_viewer:
    st.subheader("🔬 Variant Viewer")
    st.markdown("Select a modeled variant to inspect its predicted 3D structure and PAE interaction matrix.")
    
    if filtered_df.empty:
        st.info("No variants match the current filter.")
    else:
        selected_design = st.selectbox("Select Variant to Inspect", options=filtered_df["design_id"].tolist())
        target_row = filtered_df[filtered_df["design_id"] == selected_design].iloc[0]
        
        cif_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "..", "Local", "TIMP-Dashboard_output", 
            target_row.get("target_path", target_row["target"]), target_row["loop_combo"], "rf3", f"{selected_design}_refolded.cif"
        )
        conf_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "..", "Local", "TIMP-Dashboard_output", 
            target_row.get("target_path", target_row["target"]), target_row["loop_combo"], "rf3", f"{selected_design}_confidences.json"
        )
        
        c1, c2 = st.columns([1, 1])
        with c1:
            if os.path.exists(cif_path):
                import json
                import tempfile
                import numpy as np
                import biotite.structure.io as strucio
                import biotite.structure as struc
                
                try:
                    atom_array = strucio.load_structure(cif_path)
                    if os.path.exists(conf_path):
                        with open(conf_path, "r") as f: conf = json.load(f)
                        if "plddt" in conf:
                            plddts = np.array(conf["plddt"]).flatten()
                            res_starts = struc.get_residue_starts(atom_array)
                            for i, start_idx in enumerate(res_starts):
                                end_idx = res_starts[i+1] if i+1 < len(res_starts) else len(atom_array)
                                if i < len(plddts):
                                    atom_array.b_factor[start_idx:end_idx] = plddts[i]
                                    
                    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
                        strucio.save_structure(tmp.name, atom_array)
                        with open(tmp.name, "r") as r: model_data = r.read()
                    os.remove(tmp.name)
                    fmt = "pdb"
                except Exception as e:
                    with open(cif_path, "r") as f: model_data = f.read()
                    fmt = "mmcif"
                    
                view = py3Dmol.view(width=600, height=600)
                view.addModel(model_data, fmt)
                
                # RF3 outputs pLDDT on a 0.0 to 1.0 logic scale, AlphaFold uses 0 to 100.
                view.setStyle({'chain': 'A'}, {'cartoon': {'colorscheme': {'prop': 'b', 'gradient': 'roygb', 'min': 0.4, 'max': 1.0}}})
                view.setStyle({'chain': 'B'}, {'sphere': {'color': '#cccccc'}})
                view.zoomTo()
                components.html(view._make_html(), height=600)
                st.caption(" TIMP3 (Chain A) as ribbons colored by pLDDT. Target (Chain B) as grey spheres.")
            else:
                st.error(f"Structure CIF file not found for {selected_design}")
                
        with c2:
            st.markdown(f"**ipTM**: {target_row.get('iptm', 'N/A')} | **Mean pLDDT**: {target_row.get('plddt', 'N/A')} | **Loop PAE**: {target_row.get('mean_loop_pae', 'N/A')}")
            if os.path.exists(conf_path):
                try:
                    with open(conf_path, "r") as f: conf = json.load(f)
                    if "pae" in conf:
                        pae_matrix = np.array(conf["pae"])
                        if len(pae_matrix.shape) == 3: pae_matrix = pae_matrix[0]
                        fig_pae = px.imshow(pae_matrix, color_continuous_scale="RdYlBu_r", zmin=0, zmax=30, title="Predicted Aligned Error (PAE)")
                        st.plotly_chart(fig_pae, width='stretch', key='viewer_pae')
                except Exception as e:
                    st.error(f"Error loading PAE matrix: {e}")
            else:
                st.warning("Confidences JSON (PAE matrix) is not available for this legacy variant.")

with tab_comp:
    st.subheader("🔍 Comparative Diagnostics")
    st.markdown("Statistically normalize metrics via T-Scores across all combinations or dynamically within distinct targets.")
    
    comp_c1, comp_c2, comp_c3 = st.columns(3)
    available_metrics = [c for c in ["probability_of_binding_score", "plddt", "mean_loop_plddt", "overall_rmsd", "contacts", "bsa", "mean_loop_pae", "iptm"] if c in df.columns]
    
    with comp_c1: chart_type = st.selectbox("Chart Format", [
        "Metric vs Loop Distributions (Boxplot)", 
        "Metric vs Metric (Scatter)",
        "Target Matrix Pivot (Loop vs Target)"
    ])
    with comp_c2: current_metric = st.selectbox("Primary Metric (Y-Axis)", available_all_metrics, index=0)
    with comp_c3: score_transformation = st.selectbox("Score Transformation", ["Raw Scores", "T-Score (Overall)", "Local T-Score (Within Targets)"])

    df_comp = filtered_df.copy()
    
    if score_transformation == "T-Score (Overall)":
        mean_val = df_comp[current_metric].mean()
        std_val = df_comp[current_metric].std()
        sem_val = std_val / np.sqrt(len(df_comp)) if len(df_comp) > 1 and std_val > 0 else 1
        df_comp["plot_val"] = (df_comp[current_metric] - mean_val) / sem_val
    elif score_transformation == "Local T-Score (Within Targets)":
        df_comp["plot_val"] = 0.0
        for tgt, group in df_comp.groupby("target"):
            g_mean, g_std = group[current_metric].mean(), group[current_metric].std()
            g_sem = g_std / np.sqrt(len(group)) if len(group) > 1 and g_std > 0 else 1
            df_comp.loc[group.index, "plot_val"] = (group[current_metric] - g_mean) / g_sem
    else:
        df_comp["plot_val"] = df_comp[current_metric]

    if chart_type == "Metric vs Loop Distributions (Boxplot)":
        fig_ts = px.box(df_comp, x="loop_combo", y="plot_val", color="target", points="all", title=f"{score_transformation} of {current_metric} across Loops")
        fig_ts.update_layout(xaxis={'categoryorder':'total descending'}, yaxis_title=f"{current_transformation if 'current_transformation' in locals() else score_transformation} {current_metric}")
        st.plotly_chart(fig_ts, width='stretch')
        
    elif chart_type == "Metric vs Metric (Scatter)":
        x_metric = st.selectbox("Secondary Metric (X-Axis)", [m for m in available_metrics if m != current_metric])
        fig_ts = px.scatter(df_comp, x=x_metric, y="plot_val", color="target", symbol="loop_combo", hover_data=["design_id"], title=f"{score_transformation} of {current_metric} vs {x_metric}")
        st.plotly_chart(fig_ts, width='stretch', key='comp_scatter')

    elif chart_type == "Target Matrix Pivot (Loop vs Target)":
        try:
            pivot_df = df_comp.pivot_table(index="loop_combo", columns="target", values="plot_val", aggfunc="mean")
            st.dataframe(pivot_df)
        except Exception as e:
            st.error(f"Cannot assemble pivot table with current data mappings.")
        
with tab_heat:
    st.subheader("🔥 Target Heatmaps")
    
    heat_mode = st.radio("Data Source", ["Primary Design Metrics", "Cross-Docking Sensitivity"], horizontal=True)
    
    heat_c1, heat_c2, heat_c3 = st.columns(3)
    
    with heat_c1: heat_chart = st.selectbox("Heatmap Format", [
        "Comparison Matrix (Sequence vs Target)",
        "Variant Ranking per Target (Styled Heatmap)",
        "Sequence vs Metrics Overview"
    ])
    
    if heat_mode == "Cross-Docking Sensitivity" and not xd_df.empty:
        available_heat_metrics = ["plddt", "iptm", "ptm"]
        df_heat_source = xd_df.copy()
        # Use folded_target as the primary categorical target for heatmaps
    else:
        available_heat_metrics = available_all_metrics
        df_heat_source = filtered_df.copy()

    with heat_c2: heat_metric = st.selectbox("Primary Metric (Color)", available_heat_metrics, index=0, key="heat_metric")
    with heat_c3: heat_trans = st.selectbox("Score Transformation", ["Raw Scores", "T-Score (Overall)", "Local T-Score (Within Targets)"], key="heat_trans")
    
    df_heat = df_heat_source.copy()
    
    if heat_trans == "T-Score (Overall)":
        mean_val = df_heat[heat_metric].mean()
        std_val = df_heat[heat_metric].std()
        sem_val = std_val / np.sqrt(len(df_heat)) if len(df_heat) > 1 and std_val > 0 else 1
        df_heat["plot_val"] = (df_heat[heat_metric] - mean_val) / sem_val
    elif heat_trans == "Local T-Score (Within Targets)":
        df_heat["plot_val"] = 0.0
        for tgt, group in df_heat.groupby("target"):
            g_mean, g_std = group[heat_metric].mean(), group[heat_metric].std()
            g_sem = g_std / np.sqrt(len(group)) if len(group) > 1 and g_std > 0 else 1
            df_heat.loc[group.index, "plot_val"] = (group[heat_metric] - g_mean) / g_sem
    else:
        df_heat["plot_val"] = df_heat[heat_metric]
        
    if heat_chart == "Variant Ranking per Target (Styled Heatmap)":
        text_data, numeric_data = {}, {}
        HIGH_IS_BETTER = ['plddt', 'TIMP_pLDDT', 'Target_pLDDT', 'TIMP_pTM', 'Target_pTM', 'mean_loop_plddt', 'iptm', 'contacts', 'bsa', 'probability_of_binding_score']
        LOW_IS_BETTER = ['overall_rmsd', 'mean_loop_pae', 'clashes']
        sort_ascending = True if heat_metric in LOW_IS_BETTER else False
        
        for target_name, group in df_heat.groupby('target'):
            sorted_group = group.sort_values(by=['loop_combo', 'plot_val'], ascending=[True, sort_ascending])
            
            def format_row(row):
                loops = str(row['loop_combo']).split('_')
                seqs = []
                for l in loops:
                    val = row.get(f"loop_{l}_seq")
                    if pd.notna(val) and str(val).upper() != "MISSING":
                        seqs.append(str(val))
                seq_str = "-".join(seqs) if seqs else row['design_id'].split('_')[-1]
                return f"{row['loop_combo']} - {seq_str} (T={row['plot_val']:.2f})"
                
            text_data[target_name] = [format_row(row) for _, row in sorted_group.iterrows()]
            numeric_data[target_name] = sorted_group['plot_val'].tolist()
            
        max_len = max(len(v) for v in text_data.values()) if text_data else 0
        for k in text_data.keys():
            text_data[k].extend([""] * (max_len - len(text_data[k])))
            numeric_data[k].extend([np.nan] * (max_len - len(numeric_data[k])))
            
        rank_df_text = pd.DataFrame(text_data)
        rank_df_numeric = pd.DataFrame(numeric_data)
        rank_df_text.index = rank_df_numeric.index = rank_df_text.index + 1
        rank_df_text.index.name = "Rank"
        
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        from matplotlib import colors
        
        cmap_name = 'RdYlGn_r' if heat_metric in LOW_IS_BETTER else 'RdYlGn'
        cmap = plt.get_cmap(cmap_name)
        g_min, g_max = rank_df_numeric.min().min(), rank_df_numeric.max().max()
        if pd.isna(g_min): g_min, g_max = 0, 1
        norm = Normalize(vmin=g_min, vmax=g_max)
        
        def color_cells(val):
            if pd.isna(val): return ''
            rgba = cmap(norm(val))
            hex_color = colors.to_hex(rgba)
            lum = 0.299*rgba[0] + 0.587*rgba[1] + 0.114*rgba[2]
            text_c = 'white' if lum < 0.5 else 'black'
            return f'background-color: {hex_color}; color: {text_c}'
            
        styled_df = rank_df_text.style.apply(lambda x: rank_df_numeric.map(color_cells), axis=None)
        st.dataframe(styled_df, width='stretch')

    elif heat_chart == "Comparison Matrix (Sequence vs Target)":
        try:
            pivot_df = df_heat.pivot_table(index=["loop_combo", "sequence"], columns="target", values="plot_val", aggfunc="mean")
            
            HIGH_IS_BETTER = ['plddt', 'TIMP_pLDDT', 'Target_pLDDT', 'TIMP_pTM', 'Target_pTM', 'mean_loop_plddt', 'iptm', 'contacts', 'bsa', 'probability_of_binding_score']
            LOW_IS_BETTER = ['overall_rmsd', 'mean_loop_pae', 'clashes']
            cmap_name = 'RdYlGn_r' if heat_metric in LOW_IS_BETTER else 'RdYlGn'
            
            styled_df = pivot_df.style.background_gradient(cmap=cmap_name, axis=None)
            st.dataframe(styled_df, width='stretch')
        except Exception as e:
            st.error("Cannot assemble pivot table with current data mappings.")
            
    elif heat_chart == "Sequence vs Metrics Overview":
        try:
            metrics_to_show = [m for m in available_all_metrics if m in df_heat.columns]
            df_plot_ms = df_heat.copy()
            
            if heat_trans == "T-Score (Overall)":
                for m in metrics_to_show:
                    m_mean, m_std = df_plot_ms[m].mean(), df_plot_ms[m].std()
                    m_sem = m_std / np.sqrt(len(df_plot_ms)) if len(df_plot_ms) > 1 and m_std > 0 else 1
                    df_plot_ms[m] = (df_plot_ms[m] - m_mean) / m_sem
            elif heat_trans == "Local T-Score (Within Targets)":
                for m in metrics_to_show:
                    for tgt, group in df_plot_ms.groupby("target"):
                        g_mean, g_std = group[m].mean(), group[m].std()
                        g_sem = g_std / np.sqrt(len(group)) if len(group) > 1 and g_std > 0 else 1
                        df_plot_ms.loc[group.index, m] = (group[m] - g_mean) / g_sem
            
            pivot_ms = df_plot_ms.pivot_table(index=["loop_combo", "sequence"], values=metrics_to_show, aggfunc="mean")
            
            # Use original metrics for sort order if possible
            pivot_ms = pivot_ms[metrics_to_show]
            
            HIGH_IS_BETTER = ['plddt', 'TIMP_pLDDT', 'Target_pLDDT', 'TIMP_pTM', 'Target_pTM', 'mean_loop_plddt', 'iptm', 'contacts', 'bsa', 'probability_of_binding_score']
            LOW_IS_BETTER = ['overall_rmsd', 'mean_loop_pae', 'clashes']
            
            styled_ms = pivot_ms.style
            for m in metrics_to_show:
                c_name = 'RdYlGn_r' if m in LOW_IS_BETTER else 'RdYlGn'
                styled_ms = styled_ms.background_gradient(cmap=c_name, subset=[m])
                
            styled_ms = styled_ms.format("{:.2f}")
            st.dataframe(styled_ms, width='stretch')
        except Exception as e:
            st.error(f"Cannot assemble sequence vs metrics table: {e}")

with tab_data:
    st.subheader("Raw Data Table")
    st.dataframe(filtered_df)
