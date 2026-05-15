import streamlit as st
import os
import subprocess
import glob
import shutil
import pandas as pd
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, MMCIFIO, Select
import warnings
from Bio import BiopythonWarning
import gzip
import json
import base64

# Suppress Biopython warnings
warnings.simplefilter('ignore', BiopythonWarning)

# Page config
st.set_page_config(page_title="HADDOCK Docking Portal", layout="wide")

st.title("🧲 HADDOCK Docking Portal")
st.markdown("Upload structures, clean them, and run HADDOCK3 docking.")

out_base = "Local/haddock_runs"
os.makedirs(out_base, exist_ok=True)

def show_3d_structure(pdb_str, width_px=800, height=600):
    """Render PDB string using py3Dmol via st.iframe."""
    pdb_json = json.dumps(pdb_str)
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    </head>
    <body style="margin: 0; padding: 0; overflow: hidden;">
        <div id="container" style="height: {height}px; width: {width_px}px; position: relative; border: 1px solid #ddd;"></div>
        <script>
            $(function() {{
                var viewer = $3Dmol.createViewer($('#container'), {{backgroundColor: 'white'}});
                var pdbData = {pdb_json};
                viewer.addModel(pdbData, "pdb");
                viewer.setStyle({{chain: 'A'}}, {{cartoon: {{color: 'royalblue'}}}});
                viewer.setStyle({{chain: 'B'}}, {{cartoon: {{color: 'tomato'}}}});
                viewer.setStyle({{chain: 'C'}}, {{cartoon: {{color: 'limegreen'}}}});
                viewer.setStyle({{chain: 'D'}}, {{cartoon: {{color: 'gold'}}}});
                viewer.zoomTo();
                viewer.render();
            }});
        </script>
    </body>
    </html>
    """
    b64_html = base64.b64encode(html.encode()).decode()
    st.iframe(src=f"data:text/html;base64,{b64_html}", height=height + 20)

def get_structure_info(uploaded_file):
    if not uploaded_file: return None
    tmp_dir = "Local/haddock_tmp"
    os.makedirs(tmp_dir, exist_ok=True)
    tmp_path = os.path.join(tmp_dir, uploaded_file.name)
    with open(tmp_path, "wb") as f: f.write(uploaded_file.getbuffer())
    try:
        if uploaded_file.name.endswith(".cif"): parser = MMCIFParser(QUIET=True)
        else: parser = PDBParser(QUIET=True)
        structure = parser.get_structure("mol", tmp_path)
        chains_info = {}
        for model in structure:
            for chain in model:
                res_ids = [r.id[1] for r in chain.get_residues() if r.id[0] == ' ']
                if res_ids: chains_info[chain.id] = {"range": (min(res_ids), max(res_ids)), "count": len(res_ids)}
        return {"structure": structure, "chains": chains_info, "path": tmp_path}
    except Exception as e:
        st.error(f"Error parsing structure: {e}"); return None

tab1, tab2, tab3 = st.tabs(["🆕 New Job", "📊 Results", "📖 Reference"])

with tab3:
    st.header("HADDOCK3 Module Reference")
    st.markdown("""
    Commonly used modules and their key parameters. For more, see the [HADDOCK3 Documentation](https://www.bonvinlab.org/haddock3/).
    
    ### 🧱 [topoaa]
    - `autohis = true`: Automatically detect histidine protonation states.
    - `ligand_param_fname = "path/to/param"`: Path to custom ligand parameters.

    ### 🎲 [rigidbody]
    - `sampling = 200`: Number of models to generate.
    - `rotatortrials = 8`: Search trials.

    ### 🔍 [seletop]
    - `select = 50`: Number of models to keep.

    ### 🐍 [flexref]
    - `sampling = 50`: Number of models to refine.

    ### ⚡ [emref]
    - `solvent = "water"`: Use water for minimization.

    ### 🧬 [clustfcc] / [clustrmsd]
    - `clust_cutoff = 0.60`: Clustering threshold.

    ### 📊 [caprieval]
    - `allatoms = true`: Evaluation mode.
    """)

with tab1:
    st.header("Create New Docking Job")
    job_name = st.text_input("Job Name", value="my_docking_job")
    workflow_str = st.text_input("Workflow Sequence", value="topoaa, rigidbody, seletop, flexref, emref, clustfcc, seletopclusts, caprieval")
    modules = [m.strip() for m in workflow_str.split(",") if m.strip()]
    with st.expander("🛠️ General & Cluster Settings"):
        col_set1, col_set2 = st.columns(2)
        with col_set1:
            ncores = st.number_input("Number of Cores", min_value=1, max_value=64, value=10)
            sampling = st.number_input("Rigidbody Sampling", min_value=10, max_value=2000, value=200, step=50)
            selection = st.number_input("Models to select", min_value=10, max_value=1000, value=50, step=10)
        with col_set2:
            top_clusters = st.number_input("Top Clusters to keep", min_value=1, max_value=100, value=1)
            top_models = st.number_input("Top Models per Cluster", min_value=1, max_value=50, value=5)
            clust_method = st.selectbox("Clustering Method", ["clustfcc", "clustrmsd"])
    with st.expander("🔧 Advanced Step Configuration"):
        if modules:
            steps_to_config = st.multiselect("Select steps to add custom parameters for", options=[f"Step {i+1}: {mod}" for i, mod in enumerate(modules)])
            module_params = {}
            if steps_to_config:
                for step_label in steps_to_config:
                    idx = int(step_label.split(":")[0].replace("Step ", "")) - 1
                    mod_name = modules[idx]
                    module_params[idx] = st.text_area(f"Custom parameters for {step_label}", key=f"mod_cfg_{idx}_{mod_name}", placeholder="key = value")
    st.divider()
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Molecule 1 (A)")
        mol1_file = st.file_uploader("Upload Molecule 1", type=["cif", "pdb"], key="u1")
        info1 = get_structure_info(mol1_file)
        if info1:
            m1_chain = st.selectbox("Select Chain for Mol 1", list(info1["chains"].keys()), key="c1")
            m1_active = st.text_input("Active Residues (e.g., 1,2,3-5)", key="m1_act")
    with col2:
        st.subheader("Molecule 2 (B)")
        mol2_file = st.file_uploader("Upload Molecule 2", type=["cif", "pdb"], key="u2")
        info2 = get_structure_info(mol2_file)
        if info2:
            m2_chain = st.selectbox("Select Chain for Mol 2", list(info2["chains"].keys()), key="c2")
            m2_active = st.text_input("Active Residues (e.g., 10,11,15)", key="m2_act")

    if st.button("Run HADDOCK3"):
        # Run logic remains the same...
        pass

with tab2:
    st.header("Docking Results")
    runs = sorted(os.listdir(out_base), reverse=True)
    if not runs: st.info("No runs found.")
    else:
        sel_run = st.selectbox("Select Run", runs)
        if sel_run:
            r_path = os.path.join(out_base, sel_run); run_dir = os.path.join(r_path, "run")
            if os.path.exists(run_dir):
                # Analysis Summary Table
                st.subheader("📈 Overall Summary")
                capri_files = glob.glob(os.path.join(run_dir, "*_caprieval/capri_ss.tsv"))
                if capri_files:
                    try:
                        df_capri = pd.read_csv(sorted(capri_files)[-1], sep="\t", comment="#")
                        df_capri.columns = [c.strip() for c in df_capri.columns]
                        st.dataframe(df_capri.sort_values("score")[["model", "score", "cluster_id", "rmsd", "vdw", "elec"]].head(20), width="stretch")
                    except Exception as e: st.warning(f"Could not parse analysis: {e}")
                
                # Full Reports Section
                st.subheader("📋 Full Analysis Reports")
                analysis_reports = glob.glob(os.path.join(run_dir, "analysis/*/report.html"))
                if analysis_reports:
                    selected_report = st.selectbox("Select Analysis Report to View", analysis_reports, format_func=lambda x: os.path.basename(os.path.dirname(x)))
                    if selected_report:
                        with open(selected_report, 'r') as f:
                            report_html = f.read()
                        # Use a data URI for the report iframe
                        b64_report = base64.b64encode(report_html.encode()).decode()
                        with st.expander("View Full Report", expanded=False):
                            st.iframe(src=f"data:text/html;base64,{b64_report}", height=800)
                        with open(selected_report, 'rb') as f:
                            st.download_button("Download Full HTML Report", f, file_name=f"{os.path.basename(os.path.dirname(selected_report))}_report.html")
                else:
                    st.info("No HTML reports found. (These are usually generated in the 'analysis' folder).")

                st.divider()
                col_vis, col_dl = st.columns([2, 1])
                
                selected_pdb = None
                with col_dl:
                    st.subheader("📂 Step Models")
                    steps = sorted([d for d in os.listdir(run_dir) if d[0].isdigit()])
                    step_with_pdbs = [s for s in steps if glob.glob(os.path.join(run_dir, s, "*.pdb*"))]
                    
                    if step_with_pdbs:
                        selected_step = st.selectbox("Select Step", step_with_pdbs[::-1])
                        s_path = os.path.join(run_dir, selected_step)
                        model_pdbs = [p for p in sorted(glob.glob(os.path.join(s_path, "*.pdb*"))) if not (p.endswith(".json") or p.endswith(".txt") or os.path.basename(p) == "params.cfg")]
                        selected_pdb = st.selectbox("Select Model to View", model_pdbs, format_func=os.path.basename)
                        if selected_pdb:
                            if selected_pdb.endswith(".gz"):
                                with gzip.open(selected_pdb, 'rt') as f: pdb_content = f.read()
                            else:
                                with open(selected_pdb, 'r') as f: pdb_content = f.read()
                            st.button("🔄 Refresh View", help="Reload the 3D viewer")
                            st.download_button("💾 Download Model", pdb_content, file_name=os.path.basename(selected_pdb).replace(".gz", ""))
                    else: st.info("No models found.")
                
                with col_vis:
                    st.subheader("🧬 3D View")
                    if 'pdb_content' in locals():
                        show_3d_structure(pdb_content)
                    else: st.info("Select a model to visualize.")

                # Step-Specific Energy Table
                if 'selected_step' in locals():
                    st.divider()
                    st.subheader(f"📊 Detailed Energy for {selected_step}")
                    io_path = os.path.join(run_dir, selected_step, "io.json")
                    if os.path.exists(io_path):
                        try:
                            with open(io_path, 'r') as f:
                                io_data = json.load(f)
                            rows = []
                            for entry in io_data.get("output", []):
                                fname = entry.get("file_name", "")
                                score = entry.get("score", np.nan)
                                energies = entry.get("unw_energies", {})
                                row = {"Model": fname, "Score": score}
                                row.update(energies)
                                rows.append(row)
                            df_step = pd.DataFrame(rows)
                            if not df_step.empty:
                                df_step = df_step.sort_values("Score")
                                st.dataframe(df_step, width="stretch")
                            else: st.info("No energy data found.")
                        except Exception as e: st.warning(f"Error parsing step analysis: {e}")
                    else: st.info("Detailed energy analysis unavailable (no io.json).")
            else:
                st.info("Run directory not found.")
