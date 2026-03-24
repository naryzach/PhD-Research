#!/bin/bash

# Setup SLURM
source SLURM_credentials.sh
#SBATCH --job-name=TimpPipeline
#SBATCH --time=5-00:00:00         # Set a time limit (5 days)
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G # Amount of memory to reserve
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

# Open necessary virtual environment
cd ../TIMP-Dashboard
pwd
source ~/miniconda3/etc/profile.d/conda.sh

# Run Pipeline
conda activate foundry
export DISABLE_CUEQUIVARIANCE=1
python run_pipeline.py
python analyze_results.py
python cross_docking_test.py
conda deactivate