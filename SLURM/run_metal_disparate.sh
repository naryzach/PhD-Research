#!/bin/bash

# Setup SLURM
source SLURM_credentials.sh
#SBATCH --job-name=RFDiffusion
#SBATCH --time=5-00:00:00         # Set a time limit (5 days)
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G # Amount of memory to reserve
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

# Open necessary virtual environment
cd ../MetalBinder
pwd
source ~/miniconda3/etc/profile.d/conda.sh

# Run RFdiffusion2
#conda activate rfd2_env
#python run_pipeline.py --catalog --scaffold --swap --graft --num_designs 5
#conda deactivate

# Run LigandMPNN
conda activate ligandmpnn
python run_pipeline.py --run_ligandmpnn
conda deactivate

# Run AllMetal3D
conda activate allmetal3d
python run_pipeline.py --run_allmetal3d_ligand
conda deactivate

conda activate
python run_pipeline.py --finalize
conda deactivate