#!/bin/bash

# Setup SLURM
source SLURM_credentials.sh
#SBATCH --job-name=Debug
#SBATCH --time=5-00:00:00         # Set a time limit (5 days)
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G # Amount of memory to reserve
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

# Open necessary virtual environment
pwd
source ~/miniconda3/etc/profile.d/conda.sh

# Run Pipeline
conda activate foundry

pip install torchvision==0.21.0 --index-url https://download.pytorch.org/whl/cu124
python -c "import torch, torchvision; print(torch.__version__, torchvision.__version__)"