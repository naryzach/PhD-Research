# MetalBinder Pipeline

This repository contains the scripts for the MetalBinder pipeline, designed to engineer metal-binding sites into proteins using RFdiffusion, ProteinMPNN/LigandMPNN, and AllMetal3D.

## Pipeline Overview

The main entry point is `run_pipeline.py`, which orchestrates the following stages:

1.  **Scaffold Generation (RFdiffusion)**: Generates backbone designs.
2.  **Sequence Design (ProteinMPNN / LigandMPNN)**: Designs sequences for the backbones.
3.  **Validation (AllMetal3D)**: Predicts metal binding probabilities for the designed sequences.
4.  **Finalization**: Aggregates results and prepares inputs for AlphaFold validation.

## Usage

### Complete Pipeline Run
To run the full pipeline, including all standard steps:
```bash
python MetalBinder/run_pipeline.py --run_all
```

### Modular Execution
You can run specific stages of the pipeline using flags:

*   **RFdiffusion**: `--run_scaffold`
*   **ProteinMPNN**: `--run_mpnn`
*   **LigandMPNN**: `--run_ligandmpnn` (Generates packed PDBs with sidechains)
*   **AllMetal3D (Standard)**: `--run_allmetal3d` (Uses ProteinMPNN seqs threaded onto backbones)
*   **AllMetal3D (LigandMPNN)**: `--run_allmetal3d_ligand` (Uses LigandMPNN packed structures directly)
*   **Finalize**: `--finalize`

**Example:**
```bash
# Run LigandMPNN and then validate with AllMetal3D
python MetalBinder/run_pipeline.py --run_ligandmpnn
python MetalBinder/run_pipeline.py --run_allmetal3d_ligand
```

### Global Options
*   `--pred_dir <path>`: Directory containing predictions (Default: `../Local/Metal_Predictions`)
*   `--overwrite`: Overwrite existing outputs.
*   `--dry_run`: Print commands without executing them.
*   `--water`: Enable Water3D prediction in AllMetal3D steps.

## Script Descriptions

### Orchestrator
*   `run_pipeline.py`: The master script that calls other scripts in the correct order and manages environment switching (e.g., activating `allmetal3d` environment when needed).

### Generation & Design
*   `run_scaffold.py`: Runs RFdiffusion to generate protein backbones.
*   `run_mpnn.py`: Runs ProteinMPNN to design sequences for the backbones.
*   `run_ligandmpnn.py`: Runs LigandMPNN, specifically designed for ligand-context designs. It enables sidechain packing by default to produce full-atom structures.

### Validation
*   `run_allmetal3d.py`: Runs AllMetal3D on ProteinMPNN designs. It threads the designed sequence onto the backbone before prediction.
*   `run_allmetal3d_ligand.py`: Runs AllMetal3D on LigandMPNN designs (using packed structures). Outputs to `.../allmetal3d_ligand/`.

### Utilities
*   `finalize_designs.py`: Collects all results into a summary CSV (`design_summary.csv`) and generates JSON files for AlphaFold (`alphafold_inputs.json`).
*   `download_datasets.py`: Tools for downloading PDB datasets.
*   `metal_miner.py`: Mining metal binding sites from PDBs.
*   `catalog_sites.py`: Cataloging observed metal sites.
*   `run_graft.py` / `run_swap.py`: Specialized generation scripts for grafting or swapping motifs.
*   `thread_sequence.py`: Utility to thread a sequence onto a PDB backbone.
*   `utils_rfd2.py`: Helper functions for RFdiffusion integration.

## Environment Requirements

The pipeline expects specific Conda environments to be available:
*   `SE3nvNew`: For RFdiffusion.
*   `ligandmpnn`: For LigandMPNN.
*   `allmetal3d`: For AllMetal3D validation.
*   `base`: General purpose usage (run `run_pipeline.py` from here).

The `run_pipeline.py` script usually handles switching (via `conda run`), but for `allmetal3d` steps, you may need to verify the environment is set up correctly if running manually.
