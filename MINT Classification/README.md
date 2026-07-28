# MINT Binder Classification

Fine-tunes **MINT** ([VarunUllanat/mint](https://github.com/VarunUllanat/mint), an ESM2-650M-based
protein-protein interaction model trained on 96M real PPIs from STRING-DB) to predict whether a
TIMP3-scaffold loop-graft sequence **binds a given protease target**.

This is a sibling pipeline to `../ESM-C Classification/`, not a drop-in backbone swap -- MINT is
architecturally a different kind of model. Read "Why this is structured differently" before using it.

## Why this is structured differently from the ESM-C pipeline

| ESM-C pipeline | This pipeline |
|---|---|
| Input: **one** sequence (the TIMP3 variant) | Input: a **pair** -- (TIMP3 variant, target protease sequence) |
| Target = one of 3 output **heads** (multi-task, masked loss) | Target = **chain B of the input**, one shared head |
| Backbone never saw the target's sequence | Backbone is *specifically pretrained on interacting pairs* |
| AB-loop/C-loop target mismatch forced separate multirun variants | Target is just an input value -- no head-count/target-overlap issues |
| Loop-region pooling (only the graft embeds into the head) | Whole-pair mean pooling (MINT has no loop-aware pooling hook) |

The pair framing follows directly from what MINT actually is: it's trained to produce an
interaction-*aware* embedding of two sequences together, not a general-purpose single-sequence
encoder, so scoring "does A bind B" only really makes sense by giving it A and B jointly. Your data
already has a `Target Seq` column (the protease's full sequence) that the ESM-C pipeline doesn't use
at all -- this pipeline is built entirely around that column.

**Cost tradeoff:** because there's no shared forward pass across targets, scoring one candidate
sequence against all 3 targets costs 3 forward passes instead of 1. Keep that in mind if you later
want a MINT equivalent of the ESM-C pipeline's 64M-combination enumeration -- it would be ~3x the
compute of the ESM-C version, per target you enumerate against.

**Not yet ported from the ESM-C pipeline:** SHAP hotspot analysis, the 64M-combination enumeration,
UMAP visualization, and the 5-variant multirun wrapper. They all follow the same black-box
graft-and-score pattern used there (fix a scaffold window, substitute candidate loops, score), so
porting them is mechanical once this core pipeline is confirmed working -- ask if you want them.

**Unverified assumption (flagging honestly):** `model.py`'s `_backbone()` assumes
`MINTWrapper` exposes the underlying ESM2 model as `self.model` (inferred from
`mint/helpers/extract.py`'s forward method referencing `self.model.cls_idx`/`self.model.eos_idx`).
This has **not been tested against an actual installed `mint` package** -- there's no GPU/checkpoint
available in the environment this was written in. If `_backbone()` raises `AttributeError` on the
server, open the installed `mint.helpers.extract.MINTWrapper.__init__` and fix that one line.

## Install

### 1. The `mint` package itself (once, on the server -- not pip-installable from PyPI)

```bash
git clone https://github.com/VarunUllanat/mint.git
cd mint
conda env create --name mint --file=environment.yml
conda activate mint
pip install -e .

# checkpoint (~2.5GB, ESM2-650M-based) -- confirm the exact size/prompt before downloading
wget https://huggingface.co/varunullanat2012/mint/resolve/main/mint.ckpt
```

Note the absolute paths to `mint/data/esm2_t33_650M_UR50D.json` and the downloaded `mint.ckpt` --
`config.yaml`'s `mint.config_json` / `mint.checkpoint_path` need them.

### 2. This pipeline's own dependencies (same `mint` conda env, or a separate one -- either works
   since `mint` itself is only imported lazily by `model.py`)

```bash
pip install -r requirements.txt
```

### 3. Point `config.yaml` at your install

Edit `mint.config_json` and `mint.checkpoint_path` in `config.yaml` to the absolute paths from step 1.

## Pipeline

```bash
# 1. Prepare data: resolve conflicts, leakage-safe loop-group split ->
#    <output_dir>/data/{train,val,test}.parquet + manifest.json (incl. target_sequences)
python data_prep.py --config config.yaml

# 2. Fine-tune (staged: frozen backbone -> unfreeze via phase2.freeze_percent, early stop on val).
#    Saves model + per-target thresholds + held-out test report to <output_dir>/model/
python train.py --config config.yaml
#    --strict-test  also reports the novel-loop subset (no train near-neighbour)

# 3a. Predict a single sequence against every known target (or one, or an arbitrary target sequence)
python predict_seq.py --config config.yaml --seq <SEQ>
python predict_seq.py --config config.yaml --seq <SEQ> --target MMP9
python predict_seq.py --config config.yaml --seq <SEQ> --seqB <SOME_OTHER_PROTEIN_SEQ>

# 3b. Predict a CSV of sequences (see predict_csv.py docstring for the 3 modes)
python predict_csv.py --config config.yaml --input path/to/seqs.csv --seq-col "Full Seq"
```

### Quick smoke run (tiny + fast, validates the whole chain)

```bash
python data_prep.py --config config.yaml --limit 300
python train.py     --config config.yaml --smoke      # writes <output_dir>/model_smoke/
```

## Configuration highlights (`config.yaml`)

- `data.seqA_col` / `data.seqB_col` -- the TIMP3-variant column (`Full Seq`) and target column
  (`Target Seq`). No `targets:` list like the ESM-C pipeline -- targets are just whatever appears
  in `data.target_name_col`.
- `mint.sep_chains` -- `true` (default, recommended by MINT's own README) pools each chain
  separately and concatenates (2560-dim); `false` mean-pools both chains together (1280-dim).
- `mint.truncation_seq_length` -- CollateFn crop length; protease target sequences can run long.
- `train.metric` (`pr_auc`) -- threshold-free val selection metric across epochs, same rationale
  as the ESM-C pipeline (robust to per-target class imbalance).
- `train.threshold_metric` (`fbeta`) / `train.fbeta_beta` (`0.5`) -- picks each target's final
  decision threshold on val; beta<1 biases toward precision (fewer false "binder" calls).
- `train.phase2.freeze_percent` -- fraction of the backbone kept **frozen**, counted from the
  input side (0.5 unfreezes the final 50% of layers). Note the inverted convention vs the ESM-C
  pipeline's `unfreeze_layers` (an int count of trainable top layers) -- this one matches MINT's
  own vocabulary instead.

## Outputs

- `<output_dir>/data/` -- `train/val/test.parquet` + `manifest.json` (class balance, pos_weight,
  `target_sequences` lookup used by `predict_seq.py`/`predict_csv.py`).
- `<output_dir>/model/` -- `model_state.pt`, `model_meta.json` (targets, per-target decision
  thresholds, `target_sequences`), `test_report.{txt,json}`.
- `<output_dir>/predictions/` -- timestamped prediction CSVs.

## Notes / caveats

- Only rows with a target in `data.target_name_col` are usable; if you want to score against a
  target the model never saw, supply its sequence explicitly (`--seqB` / `--seqB-col`) -- MINT's
  pretraining means this is at least plausible, unlike the ESM-C pipeline where an unseen target
  has no head at all.
- Whole-pair mean pooling, not loop-region pooling -- the ESM-C pipeline found loop vs. mean
  pooling made little empirical difference there (constant scaffold, so pooling scope mostly
  became head bias), so this is a reasonable default, not untested guesswork, but it hasn't been
  validated on MINT specifically.
- Batch sizes here default smaller than the ESM-C pipeline (pairs are longer, and the backbone is
  650M) -- raise `train.batch_size` if your cluster GPU has room.
