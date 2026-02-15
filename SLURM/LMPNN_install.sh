#!/bin/bash

# Setup SLURM
source SLURM_credentials.sh
#SBATCH --job-name=LMPNNInst
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
#conda create -n ligandmpnn python=3.11
source activate ligandmpnn

pip3 install -r ../Tools/LigandMPNN/requirements.txt

python -c 'import torch; print(torch.cuda.is_available()); print(torch.__version__)'