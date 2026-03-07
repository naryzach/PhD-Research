#!/bin/bash

# Setup SLURM
source SLURM/SLURM_credentials.sh
#SBATCH --job-name=LanMDashboard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G # Amount of memory to reserve
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

# Open necessary virtual environment
streamlit run MetalBinder/dashboard.py --server.port 8502