import streamlit as st
import os
import subprocess
import glob
import shutil
import pandas as pd
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, MMCIFIO

# Page config
st.set_page_config(page_title="HADDOCK Docking Portal", layout="wide")

st.title("🧲 HADDOCK Docking Portal")
st.markdown("Upload structures and run HADDOCK3 docking.")

out_base = "Local/haddock_runs"
os.makedirs(out_base, exist_ok=True)

# Main layout
tab1, tab2 = st.tabs(["🆕 New Job", "📊 Results"])

with tab1:
    st.header("Create New Docking Job")
    job_name = st.text_input("Job Name", value="my_docking_job")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Molecule 1 (A)")
        mol1 = st.file_uploader("Upload Molecule 1 (CIF/PDB)", type=["cif", "pdb"], key="mol1")
        m1_active = st.text_input("Active Residues (e.g., 1,2,3-5)", key="m1_act")
        
    with col2:
        st.subheader("Molecule 2 (B)")
        mol2 = st.file_uploader("Upload Molecule 2 (CIF/PDB)", type=["cif", "pdb"], key="mol2")
        m2_active = st.text_input("Active Residues (e.g., 10,11,15)", key="m2_act")

    if st.button("Run HADDOCK3"):
        if not mol1 or not mol2:
            st.error("Please upload both molecules.")
        elif not job_name:
            st.error("Please provide a job name.")
        else:
            job_dir = os.path.join(out_base, f"{job_name}_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}")
            input_dir = os.path.join(job_dir, "inputs")
            run_dir = os.path.join(job_dir, "run")
            os.makedirs(input_dir, exist_ok=True)
            
            # Save files
            def save_and_preprocess(uploaded_file, name, chain):
                path = os.path.join(input_dir, f"{name}_{chain}.pdb")
                temp_path = os.path.join(input_dir, uploaded_file.name)
                with open(temp_path, "wb") as f:
                    f.write(uploaded_file.getbuffer())
                
                # Convert to PDB and change chain ID to given chain
                if uploaded_file.name.endswith(".cif"):
                    parser = MMCIFParser()
                    structure = parser.get_structure(name, temp_path)
                else:
                    parser = PDBParser()
                    structure = parser.get_structure(name, temp_path)
                
                # Force chain ID
                for model in structure:
                    for i, chain_obj in enumerate(model):
                        chain_obj.id = chain
                
                io = PDBIO()
                io.set_structure(structure)
                io.save(path)
                return path

            p1 = save_and_preprocess(mol1, "mol1", "A")
            p2 = save_and_preprocess(mol2, "mol2", "B")
            
            # Generate TBL
            def parse_residues(res_str):
                # Simple parser for "1,2,3-5"
                res = []
                for part in res_str.split(','):
                    if '-' in part:
                        start, end = part.split('-')
                        res.extend(range(int(start), int(end)+1))
                    else:
                        res.append(int(part))
                return res

            def generate_tbl(res1, res2, out_path):
                with open(out_path, 'w') as f:
                    # Assign res1 on A to all res2 on B
                    res2_str = " or ".join([f"resid {r}" for r in res2])
                    for r1 in res1:
                        f.write(f"assign (resid {r1} and segid A) (({res2_str}) and segid B) 2.0 2.0 0.0\n")
                    # Assign res2 on B to all res1 on A
                    res1_str = " or ".join([f"resid {r}" for r in res1])
                    for r2 in res2:
                        f.write(f"assign (resid {r2} and segid B) (({res1_str}) and segid A) 2.0 2.0 0.0\n")

            tbl_path = os.path.join(input_dir, "restraints.tbl")
            try:
                r1 = parse_residues(m1_active)
                r2 = parse_residues(m2_active)
                generate_tbl(r1, r2, tbl_path)
            except Exception as e:
                st.error(f"Error parsing residues: {e}")
                st.stop()

            # Generate CFG
            cfg_content = f"""
run_dir = "{run_dir}"
mode = "local"
ncores = 10
postprocess = true

molecules = [
    "{p1}",
    "{p2}"
]

[topoaa]
[rigidbody]
sampling = 200
ambig_fname = "{tbl_path}"
[seletop]
select = 50
[flexref]
[emref]
[clustfcc]
[seletopclusts]
top_clusters = 1
top_models = 5
[caprieval]
reference_fname = "{p1}"
"""
            cfg_path = os.path.join(input_dir, "haddock3.cfg")
            with open(cfg_path, 'w') as f:
                f.write(cfg_content)
            
            # Run HADDOCK3
            with st.spinner("Running HADDOCK3..."):
                log_path = os.path.join(job_dir, "haddock3.log")
                with open(log_path, 'w') as log_file:
                    res = subprocess.run(f"conda run -n haddock haddock3 {cfg_path}", shell=True, stdout=log_file, stderr=subprocess.STDOUT)
                
                if res.returncode == 0:
                    st.success("Docking complete!")
                else:
                    st.error("HADDOCK3 failed. Check log.")
                    st.text_area("Log Output", value=open(log_path).read(), height=300)

with tab2:
    st.header("Docking Results")
    runs = sorted(os.listdir(out_base), reverse=True)
    sel_run = st.selectbox("Select Run", runs)
    if sel_run:
        r_path = os.path.join(out_base, sel_run)
        # Find best models
        # HADDOCK3 usually puts them in run/output or similar
        # Based on previous script, they are in the last caprieval folder
        st.write(f"Results for {sel_run}")
        # Add viewing/downloading logic here
