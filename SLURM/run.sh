#!/bin/bash

#SBATCH --job-name=TIMP3_ESM
#SBATCH --time=5-00:00:00         # Set a time limit (5 days)
#SBATCH --gres=gpu:volta:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G # Amount of memory to reserve
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

##cp -r esm_multitarget_out_650M_ALL esm_multitarget_out

. pyvenv/bin/activate

##export OPENBLAS_NUM_THREADS=16

pip install numpy pandas
pip3 install torch torchvision #--index-url https://download.pytorch.org/whl/cu128
pip install transformers datasets scikit-learn pandas seaborn matplotlib accelerate

# Run Script
#python Training_Layered.py
python Training_with_Gen.py

deactivate

#mv esm_gen_out esm_gen_out_C
