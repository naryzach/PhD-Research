import streamlit as st
import os
import json
import numpy as np
import pandas as pd
import plotly.express as px
import py3Dmol
from stmol import showmol
from rf3_utils import RF3Runner
import zipfile
import shutil
from biotite.structure.io.pdbx import CIFFile, get_structure
from biotite.structure import superimpose, rmsd

# Page config
st.set_page_config(page_title="RF3 Generation Portal", layout="wide")

def get_b_metrics(cif_data):
    """Parses B-factors from CIF file and determines the range and if they represent pLDDT or Error."""
    if not cif_data:
        return 0.0, 100.0, False
        
    lines = cif_data.splitlines()
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
        return 0.0, 100.0, False
    
    b_min, b_max = min(b_factors), max(b_factors)
    avg_b = sum(b_factors) / len(b_factors)
    
    # Adaptive Inversion Logic:
    if b_max <= 1.1: # 0-1 scale
        invert = avg_b < 0.3 # Mostly small values = Error/Distance
    else: # 0-100 scale
        invert = avg_b < 30 # Likely RMSD or B-factor
        
    return b_min, b_max, invert

def make_view(cif_file, overlay_file=None):
    """Creates a py3Dmol view with adaptive B-factor coloring and optional overlay."""
    view = py3Dmol.view(width=800, height=600)
    
    # Load primary structure
    with open(cif_file, 'r') as f:
        cif_str = f.read()
    view.addModel(cif_str, 'cif')
    
    # Adaptive Coloring
    b_min, b_max, invert = get_b_metrics(cif_str)
    
    # If the range is tiny, fall back to a standard 50-90 scale
    if b_max - b_min < 0.0001:
        if b_max > 1.1:
            b_min, b_max = 50.0, 90.0
        else:
            b_min, b_max = 0.5, 0.9
        invert = False
        
    gradient_settings = {
        'prop': 'b', 
        'gradient': 'roygb', 
        'min': b_max if invert else b_min, 
        'max': b_min if invert else b_max
    }
    view.setStyle({'model': 0}, {'cartoon': {'colorscheme': gradient_settings}})
    
    if overlay_file:
        try:
            from biotite.structure.io.pdbx import CIFFile, get_structure
            from biotite.structure import superimpose
            from biotite.structure.io.pdb import PDBFile
            
            struct1 = get_structure(CIFFile.read(cif_file))[0]
            struct2 = get_structure(CIFFile.read(overlay_file))[0]
            
            struct1_ca = struct1[struct1.atom_name == "CA"]
            struct2_ca = struct2[struct2.atom_name == "CA"]
            
            min_len = min(len(struct1_ca), len(struct2_ca))
            fitted, transformation = superimpose(struct1_ca[:min_len], struct2_ca[:min_len])
            
            struct2_fitted = struct2.copy()
            struct2_fitted.coord = transformation.apply(struct2_fitted.coord)
            
            tmp_pdb = PDBFile()
            tmp_pdb.set_structure(struct2_fitted)
            overlay_str = "\n".join(tmp_pdb.lines)
            
            view.addModel(overlay_str, 'pdb')
            view.setStyle({'model': 1}, {'cartoon': {'color': 'gray', 'opacity': 0.5}})
            st.success("Structures superimposed successfully!")
        except Exception as e:
            st.warning(f"Superimposition failed: {e}. Showing side-by-side.")
            with open(overlay_file, 'r') as f:
                overlay_str = f.read()
            view.addModel(overlay_str, 'cif')
            view.setStyle({'model': 1}, {'cartoon': {'color': 'gray', 'opacity': 0.5}})
    
    view.zoomTo()
    return view

st.title("🧬 RF3 Generation Portal")
st.markdown("Predict protein structures and complexes using RoseTTAFold 3.")

# Sidebar - Job History and Settings
with st.sidebar:
    st.header("⚙️ Settings")
    out_base = "Local/rf3_runs"
    os.makedirs(out_base, exist_ok=True)
    
    num_cycles = st.number_input("Number of Recycles", min_value=1, max_value=20, value=10)
    num_models = st.number_input("Number of Models", min_value=1, max_value=20, value=1)
    
    st.divider()
    st.header("📂 Previous Runs")
    runs = sorted(os.listdir(out_base), reverse=True)
    selected_run = st.selectbox("Select a previous run", ["New Run"] + runs)

# Main Portal logic
if selected_run == "New Run":
    st.header("🆕 Create New Job")
    job_name = st.text_input("Job Name", value="my_rf3_job")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("⛓️ Chains")
        num_chains = st.number_input("Number of Chains", min_value=1, max_value=10, value=1)
        
        chains = []
        for i in range(num_chains):
            chain_id = chr(65 + i)
            with st.expander(f"Chain {chain_id}", expanded=True):
                c_type = st.selectbox(f"Type for {chain_id}", ["protein", "ligand", "dna", "rna"], key=f"type_{i}")
                seq = st.text_area(f"Sequence (or SMILES) for {chain_id}", key=f"seq_{i}")
                chains.append({"id": chain_id, "type": c_type, "sequence": seq.strip()})
    
    with col2:
        st.subheader("🚀 Run")
        if st.button("Generate Structure"):
            if not job_name:
                st.error("Please provide a job name.")
            elif any(not c["sequence"] for c in chains):
                st.error("Please provide sequences for all chains.")
            else:
                out_dir = os.path.join(out_base, f"{job_name}_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}")
                
                with st.spinner("Running RF3 Inference... this may take several minutes."):
                    try:
                        runner = RF3Runner(n_recycles=num_cycles, diffusion_batch_size=num_models)
                        results = runner.run_from_sequences(job_name, chains, out_dir)
                        st.success(f"Job completed! Results saved to {out_dir}")
                        st.rerun()
                    except Exception as e:
                        st.error(f"Error running RF3: {e}")
                        st.exception(e)

else:
    run_path = os.path.join(out_base, selected_run)
    st.header(f"📊 Results: {selected_run}")
    
    cif_files = sorted([f for f in os.listdir(run_path) if f.endswith(".cif") and not f.endswith("_model.cif")])
    if not cif_files:
        cif_files = sorted([f for f in os.listdir(run_path) if f.endswith(".cif")])
        
    selected_model_file = None
    if cif_files:
        if len(cif_files) > 1:
            selected_model_file = st.selectbox("Select Model", cif_files)
        else:
            selected_model_file = cif_files[0]
            st.info(f"Model: {selected_model_file}")

    if selected_model_file:
        model_prefix = selected_model_file.replace(".cif", "")
        metrics_file = f"{model_prefix}_metrics.json"
        metrics_path = os.path.join(run_path, metrics_file)
        
        if os.path.exists(metrics_path):
            with open(metrics_path, 'r') as f:
                metrics = json.load(f)
            
            m_col1, m_col2, m_col3, m_col4 = st.columns(4)
            m_col1.metric("Average pLDDT", f"{metrics.get('overall_plddt', 0):.2f}")
            m_col2.metric("pTM Score", f"{metrics.get('ptm', 0):.3f}")
            iptm = metrics.get('iptm')
            if iptm is not None:
                m_col3.metric("ipTM Score", f"{iptm:.3f}")
                ranking = metrics.get('ranking_score')
                if ranking is not None:
                    m_col4.metric("Ranking Score", f"{ranking:.3f}")
            else:
                m_col3.info("ipTM not applicable")
    
    c1, c2 = st.columns([2, 1])
    
    if selected_model_file:
        cif_path = os.path.join(run_path, selected_model_file)
        
        with c1:
            st.subheader("🏗️ 3D Structure")
            # AF3 Upload comparison
            af3_zip = st.file_uploader("Upload AF3 Server Zip for Comparison (Overlay)", type="zip")
            overlay_cif = None
            if af3_zip:
                # Use a unique extract path per run to avoid collisions
                tmp_extract = os.path.join(run_path, "af3_overlay_tmp")
                if os.path.exists(tmp_extract):
                    shutil.rmtree(tmp_extract)
                os.makedirs(tmp_extract)
                
                with zipfile.ZipFile(af3_zip, 'r') as z:
                    z.extractall(tmp_extract)
                    
                # Look for CIF files recursively (some zips have subfolders)
                all_cifs = []
                for root, dirs, files in os.walk(tmp_extract):
                    for f in files:
                        if f.endswith(".cif"):
                            all_cifs.append(os.path.join(root, f))
                
                if all_cifs:
                    if len(all_cifs) > 1:
                        selected_overlay = st.selectbox(
                            "Select Model from ZIP to Overlay", 
                            all_cifs, 
                            format_func=lambda x: os.path.basename(x)
                        )
                        overlay_cif = selected_overlay
                    else:
                        overlay_cif = all_cifs[0]
                        st.info(f"Overlaying: {os.path.basename(overlay_cif)}")
                else:
                    st.error("No CIF files found in the uploaded ZIP.")
            
            view = make_view(cif_path, overlay_cif)
            showmol(view, height=600, width=800)
            st.caption("Coloring: Red (Low Confidence) → Blue (High Confidence)")
            
        with c2:
            st.subheader("📉 PAE Plot")
            pae_file = f"{model_prefix}_pae.npy"
            pae_path = os.path.join(run_path, pae_file)
            if os.path.exists(pae_path):
                pae_data = np.load(pae_path)
                fig = px.imshow(pae_data, color_continuous_scale="Viridis_r", 
                               labels={'color': 'Error (Å)'},
                               title="Predicted Aligned Error")
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.warning("PAE data not available for this model.")
    else:
        st.error("No CIF structures found in this run directory.")
