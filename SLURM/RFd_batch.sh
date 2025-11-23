#!/bin/bash

# Setup SLURM
#SBATCH --job-name=RFDiffusion
#SBATCH --time=5-00:00:00         # Set a time limit (5 days)
#SBATCH --gres=gpu:volta:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G # Amount of memory to reserve
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

# Open necessary virtual environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate SE3nv

# Run Script
cd ../Generation
ls
python RFd_Batch.py

source deactivate
