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

st.title("🧬 RF3 Generation Portal")
st.markdown("Predict protein structures and complexes using RoseTTAFold 3.")

# Sidebar - Job History and Settings
with st.sidebar:
    st.header("⚙️ Settings")
    out_base = "Local/rf3_runs"
    os.makedirs(out_base, exist_ok=True)
    
    num_cycles = st.number_input("Number of Recycles", min_value=1, max_value=20, value=10)
    
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
                        runner = RF3Runner()
                        results = runner.run_from_sequences(job_name, chains, out_dir, num_cycles=num_cycles)
                        st.success(f"Job completed! Results saved to {out_dir}")
                        st.rerun() # Refresh to show in previous runs
                    except Exception as e:
                        st.error(f"Error running RF3: {e}")
                        st.exception(e)

else:
    # Display results of selected run
    run_path = os.path.join(out_base, selected_run)
    st.header(f"📊 Results: {selected_run}")
    
    # Load metrics
    metrics_files = [f for f in os.listdir(run_path) if f.endswith("_metrics.json")]
    if metrics_files:
        with open(os.path.join(run_path, metrics_files[0]), 'r') as f:
            metrics = json.load(f)
        
        m_col1, m_col2, m_col3 = st.columns(3)
        m_col1.metric("Average pLDDT", f"{metrics.get('overall_plddt', 0):.2f}")
        m_col2.metric("pTM Score", f"{metrics.get('ptm', 0):.3f}")
        m_col3.metric("ipTM Score", f"{metrics.get('iptm', 0):.3f}")
    
    # Structure and PAE
    c1, c2 = st.columns([2, 1])
    
    cif_files = [f for f in os.listdir(run_path) if f.endswith(".cif")]
    if cif_files:
        cif_path = os.path.join(run_path, cif_files[0])
        
        with c1:
            st.subheader("🏗️ 3D Structure")
            
            # AF3 Upload comparison
            af3_zip = st.file_uploader("Upload AF3 Server Zip for Comparison (Overlay)", type="zip")
            overlay_cif = None
            if af3_zip:
                # Extract zip and find cif
                with zipfile.ZipFile(af3_zip, 'r') as z:
                    tmp_extract = os.path.join(run_path, "af3_tmp")
                    z.extractall(tmp_extract)
                    found_cifs = [f for f in os.listdir(tmp_extract) if f.endswith(".cif") and "_model_" in f]
                    if found_cifs:
                        overlay_cif = os.path.join(tmp_extract, found_cifs[0])
                        st.info(f"Comparing with: {found_cifs[0]}")
            
            # View structure
            def make_view(cif_file, overlay_file=None):
                view = py3Dmol.view(width=800, height=600)
                
                # Load primary structure
                with open(cif_file, 'r') as f:
                    cif_str = f.read()
                view.addModel(cif_str, 'cif')
                view.setStyle({'model': 0}, {'cartoon': {'colorscheme': {'prop': 'b', 'gradient': 'roygb', 'min': 50, 'max': 90}}})
                
                if overlay_file:
                    try:
                        # Superimpose logic using biotite
                        from biotite.structure.io.pdbx import CIFFile, get_structure
                        from biotite.structure import superimpose
                        
                        # Load structures
                        struct1 = get_structure(CIFFile.read(cif_file))[0]
                        struct2 = get_structure(CIFFile.read(overlay_file))[0]
                        
                        # Filter to CA atoms for alignment
                        struct1_ca = struct1[struct1.atom_name == "CA"]
                        struct2_ca = struct2[struct2.atom_name == "CA"]
                        
                        # Align struct2 to struct1
                        min_len = min(len(struct1_ca), len(struct2_ca))
                        fitted, transformation = superimpose(struct1_ca[:min_len], struct2_ca[:min_len])
                        
                        # Apply transformation to the whole struct2
                        struct2_fitted = struct2.copy()
                        struct2_fitted.coord = transformation.apply(struct2_fitted.coord)
                        
                        # Convert fitted struct2 to PDB string for py3Dmol
                        from biotite.structure.io.pdb import PDBFile
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

                
            view = make_view(cif_path, overlay_cif)
            showmol(view, height=600, width=800)
            st.caption("Coloring: Red (pLDDT < 50) → Blue (pLDDT > 90)")
            
        with c2:
            st.subheader("📉 PAE Plot")
            pae_files = [f for f in os.listdir(run_path) if f.endswith("_pae.npy")]
            if pae_files:
                pae_data = np.load(os.path.join(run_path, pae_files[0]))
                fig = px.imshow(pae_data, color_continuous_scale="Viridis_r", 
                               labels={'color': 'Error (Å)'},
                               title="Predicted Aligned Error")
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.warning("PAE data not available for this run.")
                
    else:
        st.error("No CIF structures found in this run directory.")
