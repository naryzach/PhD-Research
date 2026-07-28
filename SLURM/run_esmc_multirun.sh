#!/bin/bash

# Setup SLURM
source SLURM_credentials.sh
#SBATCH --job-name=EsmcMultirun
#SBATCH --time=5-00:00:00         # Set a time limit (5 days) -- 5 variants x
                                   # (fine-tune + SHAP + up to two 64M enumeration sweeps)
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

# Open necessary virtual environment
cd ../"ESM-C Classification"
pwd
source ~/miniconda3/etc/profile.d/conda.sh
conda activate esmc

# Optional: run only a subset of variants (e.g. one per array-job/GPU) by
# passing e.g. --variants abloop_only,cloop_only as arguments to this script.
python multirun/run_all.py --log-every 100 "$@"

conda deactivate
