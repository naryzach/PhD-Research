#!/bin/bash

# MetalBinder Dependency Installer
# This script installs all necessary Python packages for the MetalBinder pipeline.
# It assumes you have Conda installed and are in your desired environment (e.g., SE3nv).

echo "Installing Core Dependencies..."

# 1. BioPython: Essential for PDB/CIF parsing and structural manipulation.
# Used in: almost all scripts (download, catalog, scaffold, graft, swap).
conda install -y biopython

# 2. Pandas: Used for handling CSV catalogs (catalog_sites.py, finalize_designs.py).
conda install -y pandas

# 3. Requests: Used for querying the RCSB PDB API (download_datasets.py).
conda install -y requests

echo "Installing LigandMPNN/ProteinMPNN Dependencies..."

# 4. ProDy: Protein Dynamics analysis, required by LigandMPNN's helper scripts.
conda install -y -c conda-forge prody

# 5. ML Collections: Configuration library used by AlphaFold/OpenFold/LigandMPNN.
pip install ml_collections

# 6. DM-Tree: DeepMind Tree library, often a hidden dependency for OpenFold/LigandMPNN data structures.
pip install dm-tree

# 7. NumPy/SciPy: Standard scientific stack (usually pre-installed, but ensuring).
conda install -y numpy scipy

echo "Dependencies installed successfully."
