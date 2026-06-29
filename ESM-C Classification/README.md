# ESM-C Binder Classification

Fine-tunes **ESM-C** (via the HF-native [Synthyra ESM++](https://huggingface.co/Synthyra/ESMplusplus_small)
port — identical weights to EvolutionaryScale's ESM-C) to predict, for a TIMP-scaffold
loop-graft sequence, whether it **binds each protease target** (ADAM17 / MMP3 / MMP9).

This replaces the notebook pipeline in `../ESM2 Classification/`. See
`../../.claude/plans/hazy-swimming-wave.md` for the full design rationale.

## Why this is structured the way it is

The lab data (`Data/TIMP3_Binding_Results/TIMP_binder_all.csv`) has properties that drive
every design choice here:

| Observation | Consequence |
|---|---|
| Every sequence is the **same 188-aa scaffold** differing only in a **6-aa grafted loop** | Pool only the **loop tokens** — mean-pooling the whole sequence buries the signal. |
| Each row is a **single-target binary assay** (`Encoding` 1/0 = bound that target) | Model is **multi-task binary**: one backbone, one bind/no-bind head per target, **masked loss** (each example only supervises the target it was screened against). No artificial shared "Non-Binder" class. |
| Negatives are target-specific; binding shows little selectivity | Predict P(bind) for **all** targets so you can inspect (non-)selectivity directly. |
| Severe, target-specific imbalance (MMP3 92% pos, ADAM17 27%) | Per-target `pos_weight`; report **PR-AUC / ROC-AUC / MCC**, not accuracy. |
| `Count` = read/observation confidence; 162 conflicting labels | Count-weighted conflict resolution + optional `log(count)` loss weighting. |
| Constant scaffold ⇒ random splits leak | **Loop-group split** (optionally merging edit-distance≤1 loops) with one fixed held-out test set. |

## Install

### Python environment

Use **Python 3.11 or 3.12** in a dedicated **conda or venv** environment. That range has the
widest binary-wheel coverage for the whole stack (PyTorch, `numba`/`llvmlite` → `umap-learn`,
scikit-learn). 3.13 works for the core train/predict path but `umap-learn` needs `numba >= 0.61`
on 3.13 and is more likely to give install friction.

Avoid the **Microsoft Store** Python build (`...WindowsApps\PythonSoftwareFoundation...`): it is
sandboxed, which is what triggers the Hugging Face symlink-cache warning and can complicate native
packages like `numba`. A regular CPython / Miniconda install in a venv is much smoother. Example:

```bash
conda create -n esmc python=3.11 -y && conda activate esmc
pip install -r requirements.txt
```

`umap-learn` is only needed for `visualize.py`; if it is missing, that script automatically
falls back to PCA, so the rest of the pipeline runs on any of 3.11–3.13.

First run downloads the ESM-C weights (~1.2 GB for the 300M model) from Hugging Face.

## Pipeline

All scripts read `config.yaml` (paths are resolved relative to it).

```bash
# 1. Prepare data: resolve conflicts, pivot to per-sequence multi-task labels,
#    locate loops, leakage-safe split -> <output_dir>/data/{train,val,test}.parquet
python data_prep.py --config config.yaml

# 2. Fine-tune (staged: frozen heads -> unfreeze top-N layers, early stop on val).
#    Saves model + thresholds + held-out test report to <output_dir>/model/
python train.py --config config.yaml
#    --strict-test  also reports the novel-loop subset (no train near-neighbour)
#    --model-id ... try a different backbone without editing the config

# 3a. Predict a single sequence (give --loop if you know the grafted motif)
python predict_seq.py --config config.yaml --seq <SEQ> --loop AADSIY

# 3b. Predict a CSV of sequences
python predict_csv.py --config config.yaml --input path/to/seqs.csv \
    --seq-col "Full Seq" --loop-col "Residues"

# 4. (optional) Visualize the embedding space with UMAP
python visualize.py --config config.yaml --split test               # color = assayed target
python visualize.py --config config.yaml --split test --color-by binding

# 5. (optional) Compare pooling strategies head-to-head on the held-out test
python compare_pooling.py --config config.yaml --modes loop mean cls

# 6. (optional) Selectivity: does the model see binders as non-selective?
python selectivity.py --config config.yaml --split test
```

### Running on a cluster / bigger GPU

`config.yaml` -> `train:` controls scaling, so the same code runs on a laptop GPU or a
big cluster card:

- `precision: auto` — uses **bf16** on capable GPUs (A100/H100/Ampere+; no loss scaler,
  more stable), else fp16. Set explicitly to `bf16`/`fp16`/`fp32` if you prefer.
- `batch_size` — raise it on a big-VRAM card (e.g. 32–64); lower to 8 if a 12 GB GPU OOMs
  with the 600M model.
- `gradient_checkpointing: true` — trades compute for memory to fit the 600M model + large
  batches.
- `multi_gpu` — leave **false** for this model. `DataParallel` re-replicates the 600M model
  across GPUs every step, so it shows both cards "full" with little real speedup and a garbled
  progress bar. A single 32 GB GPU fits the 600M model fine. To use multiple cards productively,
  run independent jobs (`CUDA_VISIBLE_DEVICES=0/1 ...`) or a real DDP launcher (`torchrun`).
- `log_every` — `>0` prints a plain `step i/N loss=.. it/s` line every N steps; recommended on a
  cluster (tqdm auto-hides when output is redirected to a log file). Also `--log-every N`.
- `num_workers` / `pin_memory` — bump on many-core cluster nodes.

**CUDA build vs driver (common cluster gotcha):** `pip install torch` grabs the newest CUDA
build, which may be newer than the node's GPU driver supports — torch then prints "NVIDIA driver
too old" and silently runs on **CPU**. Check the driver's CUDA in `nvidia-smi` and install a
matching build, e.g. `pip install torch --index-url https://download.pytorch.org/whl/cu128`
(any `cu12x` ≤ the driver's CUDA works). `train.py` prints a loud warning if it ends up on CPU.

**Older GPUs (V100 / Volta, CC 7.0):** recent torch wheels dropped Volta kernels (symptom:
`no kernel image is available for execution on the device`). Install a torch that still ships
sm_70, e.g. `pip install torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121`
(verify with `python -c "import torch;x=torch.randn(64,64,device='cuda');print(float((x@x).sum()))"`).
That older torch is in turn missing some Dynamo config knobs the ESM++ remote code sets
(`recompile_limit`); `esmc_utils.patch_dynamo_config_compat()` neutralises those so the model
still imports. On Volta, `precision: auto` correctly selects fp16 (no native bf16).
If your cluster has an **Ampere+ partition (A100/H100)**, prefer it — the default modern torch
works as-is with bf16 and avoids this whole chain.

### Quick smoke run (tiny + fast, validates the whole chain)

```bash
python data_prep.py --config config.yaml --limit 1500
python train.py     --config config.yaml --smoke      # writes <output_dir>/model_smoke/
```

## Configuration highlights (`config.yaml`)

- `model.model_id` — `Synthyra/ESMplusplus_large` (600M, default) or `..._small` (300M).
- `model.pooling` — `loop` (default), `mean`, or `cls`. (Empirically near-identical on this
  data, since the constant scaffold is just absorbed as head bias.)
- `split.cluster_hamming1` — keep **false**: the 6-aa loops percolate into one giant
  edit-distance≤1 component (~58%), so clustering wrecks the split. `data_prep` instead reports
  a Hamming-1 exposure diagnostic and tags each eval row (`near_train_h1`) for the strict test.
- `preprocess.{min_count,count_weighting}` — use read-count as a confidence filter / loss weight.
- `train.precision / multi_gpu / gradient_checkpointing` — cluster scaling (see above).
- `train.phase1 / phase2` — epochs, LRs, and `unfreeze_layers` for the staged curriculum.

## Outputs

- `<output_dir>/data/` — `train/val/test.parquet` + `manifest.json` (class balance, `pos_weight`).
- `<output_dir>/model/` — `model_state.pt`, `model_meta.json` (targets, pooling, **per-target
  decision thresholds**), `test_report.{txt,json}`.
- `<output_dir>/model_<mode>/` — one per pooling mode from `compare_pooling.py`.
- `<output_dir>/predictions/` — timestamped prediction CSVs.
- `<output_dir>/visualizations/` — UMAP/PCA PNGs + 2D coordinate CSVs.
- `<output_dir>/selectivity/` — cross-target correlation / co-binding report + heatmaps.
- `<output_dir>/pooling_comparison_*.csv` — pooling head-to-head table.

## Notes / caveats

- The model emits an **independent probability per target**. Because the assay data shows
  binders are largely non-selective, expect correlated per-target probabilities — that is a
  finding, not a bug.
- Prediction accuracy depends on supplying the **loop motif** so loop pooling matches training.
  Generated/unknown sequences should carry a loop column (as the old `Run_Prediction` notebook
  did with `Residues`).
