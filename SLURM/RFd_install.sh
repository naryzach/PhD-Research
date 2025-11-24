#!/bin/bash

# Setup SLURM
source SLURM_credentials.sh
#SBATCH --job-name=SE3nvInst
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
cd ../Tools/RFdiffusion/

conda env create -f env/SE3nv.yml
source activate SE3nv
cd env/SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install
cd ../.. # change into the root directory of the repository
pip install -e . # install the rfdiffusion module from the root of the repository

# May need
pip install pandas
conda install pytorch==1.9 torchvision torchaudio cudatoolkit=11.1 -c pytorch -c nvidia

python -c 'import torch; print(torch.cuda.is_available()); ; print(torch.__version__)'