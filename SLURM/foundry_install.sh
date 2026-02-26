#!/bin/bash

# Setup SLURM
source SLURM_credentials.sh
#SBATCH --job-name=RFd2Inst
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

conda create --name foundry python=3.12
source activate foundry

conda install -c conda-forge gcc_linux-64 gxx_linux-64
pip install numpy pandas
pip3 install torch torchvision --index-url https://download.pytorch.org/whl/cu128
pip install "rc-foundry[all]"
pip3 install torch torchvision --index-url https://download.pytorch.org/whl/cu128 # rerun

echo $CC
$CC --version
python -c 'import torch; print(torch.cuda.is_available()); print(torch.__version__)'

conda deactivate