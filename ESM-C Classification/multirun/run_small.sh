#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate esmc
cd ~/Documents/GitHubProj/PhD-Research/"ESM-C Classification"
mkdir -p multirun/logs
python multirun/run_all.py --configs-dir configs_small --steps prepare,data,train,visualize,shap --log-every 100 \
  > multirun/logs/run_all_small.out 2>&1
