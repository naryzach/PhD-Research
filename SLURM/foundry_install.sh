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
# Detect Compute Capability
CC_VERSION=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1)
echo "Detected GPU Compute Capability: $CC_VERSION"

# Pick the best CUDA wheel index
# 7.0 (V100) requires cu124 or lower. Newer cards can use cu128.
if [[ "$CC_VERSION" == "7.0" ]]; then
    CUDA_INDEX="cu124"
else
    CUDA_INDEX="cu128"
fi

pip3 install torch torchvision --index-url https://download.pytorch.org/whl/$CUDA_INDEX
pip install "rc-foundry[all]"
pip3 install torch torchvision --index-url https://download.pytorch.org/whl/$CUDA_INDEX # rerun

echo $CC
$CC --version
python -c 'import torch; print(torch.cuda.is_available()); print(torch.__version__)'

conda deactivate