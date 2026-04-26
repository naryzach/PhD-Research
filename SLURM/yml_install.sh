#!/bin/bash

# Setup SLURM
source SLURM_credentials.sh
#SBATCH --job-name=YMLInst
#SBATCH --time=0-00:30:00
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G # Amount of memory to reserve
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

# Open necessary virtual environment
nvidia-smi
conda env create -f "../Environment Installers/combined_metal.yml" 
source activate MetalCombine

python -c 'import torch; print(torch.cuda.is_available()); print(torch.__version__)'