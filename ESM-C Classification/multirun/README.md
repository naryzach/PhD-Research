# ESM-C multirun

Runs the fine-tune -> UMAP -> SHAP hotspot -> 64M-combination enumeration ->
figures pipeline across 5 dataset variants instead of one. See `run_all.py`'s
docstring for the full variant/step breakdown.

| variant | data | targets |
|---|---|---|
| `all3_original` | `TIMP_binder_all.csv` | ADAM17, MMP3, MMP9 (AB-loop + C-loop mixed) |
| `abloop_only` | `TIMP_binder_all.csv`, AB-loop rows only | MMP3, MMP9 |
| `cloop_only` | `TIMP_binder_all.csv`, C-loop rows only | ADAM17, MMP9 |
| `mmp9_other` | `TIMP_binder_Cloop_IHTE_CGLK.csv` (separate dataset) | MMP9 |
| `everything_combined` | both raw CSVs unioned | ADAM17, MMP3, MMP9 (AB-loop + C-loop mixed) |

The two mixed-loop variants (`all3_original`, `everything_combined`) run SHAP
and enumeration **twice** -- once for the AB-loop window, once for the C-loop
window -- since a single 6-position sweep only covers one graft window.

## Run

```bash
python multirun/prepare_variant_datasets.py   # one-time; run_all.py also does this automatically
python multirun/run_all.py --list             # see the plan without running anything
python multirun/run_all.py                    # everything, sequentially
```

**Enumeration is the long pole, not training.** Fine-tuning + SHAP + UMAP for all 5
variants together is a few hours. Each 64M-sequence enumeration sweep, by contrast, has
been observed at ~76-88 seq/s (batch-size 512, fp16 on a non-Ampere GPU) -> **~180-230
hours per sweep**, and there are 7 sweeps total (2 each for `all3_original` and
`everything_combined`, 1 each for the rest). Before committing a long job to this:
raise `--enum-batch-size` well past the default 512 if VRAM allows (inference-only, no
optimizer state, so headroom is much larger than during training), and prefer an
Ampere+ (A100/H100) partition if your cluster has one (bf16, much faster than the
fp16 path older GPUs fall back to).

To fan out across GPUs/nodes instead of one long sequential job, submit one
`sbatch` per variant (see `../../SLURM/run_esmc_multirun.sh`):

```bash
sbatch run_esmc_multirun.sh --variants abloop_only
sbatch run_esmc_multirun.sh --variants cloop_only
sbatch run_esmc_multirun.sh --variants mmp9_other
sbatch run_esmc_multirun.sh --variants all3_original
sbatch run_esmc_multirun.sh --variants everything_combined
```

Each step is skipped if its output already exists (`--force` to rerun anyway),
so a crashed job can just be resubmitted. **`enumerate` and `shap` both
checkpoint mid-run**, so a kill only costs the time since the last checkpoint,
not the whole step:

- `enumerate_cloop.py` saves progress (current index + the running top-K
  heaps) to `<enum_dir>/checkpoint.json` every `--enum-checkpoint-every`
  sequences (default 1,000,000; atomic write).
- `shap_hotspots.py` explains loops in chunks (`--shap-chunk-size`, default
  10) and appends each chunk's rows to `<shap_dir>/shap_values.csv` as it
  goes (flushed + fsynced per chunk), printing a progress line after every
  chunk -- SHAP's own `KernelExplainer` call gives zero output otherwise, so
  without this a multi-hour run looks identical whether it's healthy or
  hung. Each explained loop needs ~(2^6) x background_size forward passes,
  so 300 explained loops is genuinely GPU-hours, not seconds -- check the
  log's progress lines (or `nvidia-smi` utilization) before assuming a quiet
  run is stuck.

Re-running the *exact same* `run_all.py` invocation after a crash resumes
each sweep/explain-set from its last checkpoint instead of restarting from 0
-- no flag needed, this is automatic (matched against the run's saved
parameters, so it won't silently resume onto a mismatched config). `--force`
clears any checkpoint and starts clean instead.

**`mmp9_other` fix:** `enumerate_cloop.py` used to crash immediately on any
single-target model -- the "most selective" tracking needs a runner-up target
to compare against, which doesn't exist when there's only one head. Fixed:
single-target sweeps (just `mmp9_other` in this repo) now skip
`top_selective.csv` entirely (with a one-line note in the log) instead of
crashing; `analyze_enumeration.py` already handled that file being absent, so
nothing downstream needed to change. Caught this on a local CPU smoke test
before it could crash `mmp9_other`'s real sweep on the cluster.

**Persistent full-sweep record (`--enum-save-all` / `enumerate_cloop.py
--save-all`):** off by default. Normally only the bounded top-50K-per-target
+ top-50K-selective heaps ever get written -- every other one of the 64M
evaluated loops is scored once to maybe update those heaps, then discarded,
so there's no way to retroactively ask "what did the model say about loop X"
for anything outside those lists. `--save-all` persists **every** loop's
per-target probabilities to `<enum_dir>/all_probs/<start>_<end>.parquet`
shards, one shard per checkpoint interval (same cadence and same
crash-safety guarantee as the heap checkpoint -- each shard is a fully
closed, valid parquet file the moment it's written). Read them all at once
with `pd.read_parquet(enum_dir / "all_probs")` (pandas/pyarrow read a
directory of parquet files as one dataset natively, no manual concatenation
needed). Rough size: ~1-2.5 GB per full 64M sweep with this on -- opt in per
run with `--enum-save-all` when you actually want it.

Useful flags:

- `--steps data,train` / `--steps shap,enumerate,analyze` / `--steps visualize` -- run a subset
- `--smoke` -- tiny fast pass through every step (sanity-checks the whole chain)
- `--topk` / `--enum-batch-size` / `--enum-checkpoint-every` / `--enum-save-all` -- enumerate_cloop.py tuning
- `--n-explain` / `--background-size` / `--shap-chunk-size` -- shap_hotspots.py tuning
- `--dry-run` -- print the exact commands without running them

## Outputs

Everything lands under `../../Local/esmc_multirun/<variant>/`:

- `data/` -- train/val/test parquet + manifest (from `data_prep.py`)
- `model/` -- fine-tuned weights + test report (from `train.py`)
- `visualizations/target.png` / `binding.png` (+ matching `.csv` coords) -- UMAP of the
  test-split embedding space, colored by assayed target / true binding label
  (from `visualize.py`; runs once per variant, not per loop subtype)
- `shap/<ab|c>loop/` -- SHAP hotspot heatmaps + CSVs (from `shap_hotspots.py`)
- `enumeration/<ab|c>loop_full/` -- top-K per target + selective loops from the
  full 20^6 sweep (from `enumerate_cloop.py`; no `top_selective.csv` for
  single-target models), plus `analysis/` figures (from `analyze_enumeration.py`)
  and, only with `--enum-save-all`, `all_probs/*.parquet` (every loop, not just top-K)
- `logs/` -- per-step stdout, so a multi-day cluster run leaves a trail
