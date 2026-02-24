#!/bin/bash
set -e

# Setup conda environment for Foundry
conda create -y -n foundry python=3.12
eval "$(conda shell.bash hook)"
conda activate foundry

echo "Installing PyTorch for RTX 5070 Ti (CUDA 12.8 compatible)..."
pip3 install --pre torch torchvision torchaudio --index-url https://download.pytorch.org/whl/nightly/cu128

echo "Installing Foundry..."
pip install "rc-foundry[all]"

echo "Checking installation..."
python -c "import torch; print('CUDA available:', torch.cuda.is_available()); print('PyTorch version:', torch.__version__)"
