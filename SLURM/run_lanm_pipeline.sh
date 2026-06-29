#!/bin/bash

# Setup SLURM
source SLURM_credentials.sh
#SBATCH --job-name=LanMPipeline
#SBATCH --time=5-00:00:00         # Set a time limit (5 days)
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G # Amount of memory to reserve
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

# Open necessary virtual environment
cd ../MetalBinder
pwd
source ~/miniconda3/etc/profile.d/conda.sh

# Run Pipeline
conda activate foundry
export DISABLE_CUEQUIVARIANCE=1
python run_pipeline.py --ion ZN CA CO NI CU FE SC Y LA CE PR ND PM SM EU GD TB DY HO ER TM YB LU MN --num-configs 40 --batch-size 25 --out-dir ../Local/lanm_output_tmp
python cross_docking_test.py --catalog ../Local/lanm_output_tmp/global_sequence_catalog.csv --output-dir ../Local/lanm_output_tmp/cross_docking --top-k 10
python analyze_results.py --out-dir ../Local/lanm_output_tmp
conda deactivate