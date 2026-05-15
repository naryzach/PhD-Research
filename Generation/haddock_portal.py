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

# Suppress Biopython warnings
warnings.simplefilter('ignore', BiopythonWarning)

# Page config
st.set_page_config(page_title="HADDOCK Docking Portal", layout="wide")

st.title("🧲 HADDOCK Docking Portal")
st.markdown("Upload structures, clean them, and run HADDOCK3 docking.")

out_base = "Local/haddock_runs"
os.makedirs(out_base, exist_ok=True)

class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id
    def accept_chain(self, chain):
        return chain.id == self.chain_id

def get_structure_info(uploaded_file):
    if not uploaded_file:
        return None
    
    tmp_dir = "Local/haddock_tmp"
    os.makedirs(tmp_dir, exist_ok=True)
    tmp_path = os.path.join(tmp_dir, uploaded_file.name)
    with open(tmp_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    
    try:
        if uploaded_file.name.endswith(".cif"):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
        
        structure = parser.get_structure("mol", tmp_path)
        
        chains_info = {}
        for model in structure:
            for chain in model:
                res_ids = [r.id[1] for r in chain.get_residues() if r.id[0] == ' ']
                if res_ids:
                    chains_info[chain.id] = {
                        "range": (min(res_ids), max(res_ids)),
                        "count": len(res_ids)
                    }
        return {"structure": structure, "chains": chains_info, "path": tmp_path}
    except Exception as e:
        st.error(f"Error parsing structure: {e}")
        return None

# Main layout
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
    
    workflow_str = st.text_input(
        "Workflow Sequence", 
        value="topoaa, rigidbody, seletop, flexref, emref, clustfcc, seletopclusts, caprieval"
    )
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
            # Dropdown to select which modules to configure
            steps_to_config = st.multiselect(
                "Select steps to add custom parameters for",
                options=[f"Step {i+1}: {mod}" for i, mod in enumerate(modules)]
            )
            
            module_params = {}
            if steps_to_config:
                for step_label in steps_to_config:
                    idx = int(step_label.split(":")[0].replace("Step ", "")) - 1
                    mod_name = modules[idx]
                    module_params[idx] = st.text_area(
                        f"Custom parameters for {step_label}", 
                        key=f"mod_cfg_{idx}_{mod_name}",
                        placeholder="key = value"
                    )
        else:
            st.warning("Please define a workflow sequence first.")

    st.divider()
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Molecule 1 (A)")
        mol1_file = st.file_uploader("Upload Molecule 1", type=["cif", "pdb"], key="u1")
        info1 = get_structure_info(mol1_file)
        
        m1_chain = None
        if info1:
            chain_ids = list(info1["chains"].keys())
            m1_chain = st.selectbox("Select Chain for Mol 1", chain_ids, key="c1")
            c_info = info1["chains"][m1_chain]
            st.info(f"Chain {m1_chain}: Residues {c_info['range'][0]}-{c_info['range'][1]} ({c_info['count']} residues)")
            m1_active = st.text_input("Active Residues (e.g., 1,2,3-5)", key="m1_act")
        
    with col2:
        st.subheader("Molecule 2 (B)")
        mol2_file = st.file_uploader("Upload Molecule 2", type=["cif", "pdb"], key="u2")
        info2 = get_structure_info(mol2_file)
        
        m2_chain = None
        if info2:
            chain_ids = list(info2["chains"].keys())
            m2_chain = st.selectbox("Select Chain for Mol 2", chain_ids, key="c2")
            c_info = info2["chains"][m2_chain]
            st.info(f"Chain {m2_chain}: Residues {c_info['range'][0]}-{c_info['range'][1]} ({c_info['count']} residues)")
            m2_active = st.text_input("Active Residues (e.g., 10,11,15)", key="m2_act")

    if st.button("Run HADDOCK3"):
        if not info1 or not info2:
            st.error("Please upload both molecules.")
        elif not job_name:
            st.error("Please provide a job name.")
        elif not modules:
            st.error("Workflow sequence cannot be empty.")
        else:
            job_dir = os.path.join(out_base, f"{job_name}_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}")
            input_dir = os.path.join(job_dir, "inputs")
            run_dir = os.path.join(job_dir, "run")
            os.makedirs(input_dir, exist_ok=True)
            
            def clean_and_save(info, selected_chain, target_chain, out_name):
                struct = info["structure"]
                model = struct[0]
                chain_obj = model[selected_chain]
                from Bio.PDB.Structure import Structure
                from Bio.PDB.Model import Model
                from Bio.PDB.Chain import Chain
                new_struct = Structure(out_name)
                new_model = Model(0)
                new_chain = Chain(target_chain)
                for i, residue in enumerate(chain_obj.get_residues(), start=1):
                    if residue.id[0] != ' ': continue
                    new_res = residue.copy()
                    new_res.id = (' ', i, ' ')
                    new_chain.add(new_res)
                new_model.add(new_chain)
                new_struct.add(new_model)
                out_path = os.path.join(input_dir, f"{out_name}.pdb")
                io = PDBIO()
                io.set_structure(new_struct)
                io.save(out_path)
                return out_path

            p1 = clean_and_save(info1, m1_chain, "A", "mol1")
            p2 = clean_and_save(info2, m2_chain, "B", "mol2")
            
            def parse_residues(res_str):
                res = []
                for part in res_str.split(','):
                    part = part.strip()
                    if not part: continue
                    if '-' in part:
                        start, end = part.split('-')
                        res.extend(range(int(start), int(end)+1))
                    else:
                        res.append(int(part))
                return res

            def generate_tbl(res1, res2, out_path):
                with open(out_path, 'w') as f:
                    res2_str = " or ".join([f"resid {r}" for r in res2])
                    for r1 in res1:
                        f.write(f"assign (resid {r1} and segid A) (({res2_str}) and segid B) 2.0 2.0 0.0\n")
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

            cfg_lines = [
                f'run_dir = "{run_dir}"',
                'mode = "local"',
                f'ncores = {ncores}',
                'postprocess = true',
                '',
                'molecules = [',
                f'    "{p1}",',
                f'    "{p2}"',
                ']',
                ''
            ]
            
            # Retrieve module_params logic needs to be robust to the multiselect
            # (In streamlit, the module_params dict is already populated by the loop above)
            
            for idx, mod in enumerate(modules):
                cfg_lines.append(f'[{mod}]')
                if mod == 'rigidbody':
                    cfg_lines.append(f'sampling = {sampling}')
                    cfg_lines.append(f'ambig_fname = "{tbl_path}"')
                elif mod == 'seletop':
                    cfg_lines.append(f'select = {selection}')
                elif mod == 'seletopclusts':
                    cfg_lines.append(f'top_clusters = {top_clusters}')
                    cfg_lines.append(f'top_models = {top_models}')
                elif mod == 'caprieval':
                    cfg_lines.append(f'reference_fname = "{p1}"')
                
                if idx in module_params and module_params[idx].strip():
                    cfg_lines.append(module_params[idx].strip())
                cfg_lines.append('')
            
            cfg_path = os.path.join(input_dir, "haddock3.cfg")
            with open(cfg_path, 'w') as f:
                f.write("\n".join(cfg_lines))
            
            with st.spinner("Running HADDOCK3..."):
                log_path = os.path.join(job_dir, "haddock3.log")
                haddock_cmd = f"haddock3 {cfg_path}"
                if not shutil.which("haddock3"):
                    haddock_cmd = f"conda run --no-capture-output -n haddock {haddock_cmd}"
                with open(log_path, 'w') as log_file:
                    res = subprocess.run(haddock_cmd, shell=True, stdout=log_file, stderr=subprocess.STDOUT)
                if res.returncode == 0:
                    st.success("Docking complete!")
                else:
                    st.error("HADDOCK3 failed. Check log.")
                    st.text_area("Log Output", value=open(log_path).read(), height=300)

with tab2:
    st.header("Docking Results")
    runs = sorted(os.listdir(out_base), reverse=True)
    if not runs:
        st.info("No runs found.")
    else:
        sel_run = st.selectbox("Select Run", runs)
        if sel_run:
            r_path = os.path.join(out_base, sel_run)
            st.write(f"Results for {sel_run}")
            run_dir = os.path.join(r_path, "run")
            if os.path.exists(run_dir):
                output_dirs = sorted([d for d in os.listdir(run_dir) if d[0].isdigit()], reverse=True)
                if output_dirs:
                    st.subheader("Model Files")
                    latest_dir = os.path.join(run_dir, output_dirs[0])
                    pdbs = glob.glob(os.path.join(latest_dir, "*.pdb"))
                    if pdbs:
                        for p in sorted(pdbs):
                            col_p1, col_p2 = st.columns([3, 1])
                            col_p1.text(os.path.basename(p))
                            with open(p, 'rb') as f:
                                col_p2.download_button("Download", f, file_name=os.path.basename(p), key=p)
                    else:
                        st.info("No PDB models found in the latest step.")
                else:
                    st.info("No output steps found. Docking might be in progress or failed.")
            else:
                st.info("Run directory not found.")
