#!/bin/bash

# Setup SLURM
source SLURM_credentials.sh
#SBATCH --job-name=AM3DInst
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
conda create -n allmetal3d python=3.10
source activate allmetal3d

pip install torch==1.12.1+cu116 --extra-index-url https://download.pytorch.org/whl/cu116
pip install allmetal3d
pip install pdbfixer

python -c 'import torch; print(torch.cuda.is_available()); print(torch.__version__)'