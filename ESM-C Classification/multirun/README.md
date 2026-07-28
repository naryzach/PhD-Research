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
python multirun/run_all.py                    # everything, sequentially, ~5 days on one GPU
```

To fan out across GPUs/nodes instead of one long sequential job, submit one
`sbatch` per variant (see `../../SLURM/run_esmc_multirun.sh`):

```bash
sbatch run_esmc_multirun.sh --variants abloop_only
sbatch run_esmc_multirun.sh --variants cloop_only
sbatch run_esmc_multirun.sh --variants mmp9_other
sbatch run_esmc_multirun.sh --variants all3_original
sbatch run_esmc_multirun.sh --variants everything_combined
```

Each step is skipped if its output already exists (`--force` to rerun
anyway), so a crashed job can just be resubmitted. Useful flags:

- `--steps data,train` / `--steps shap,enumerate,analyze` / `--steps visualize` -- run a subset
- `--smoke` -- tiny fast pass through every step (sanity-checks the whole chain)
- `--topk` / `--enum-batch-size` -- enumerate_cloop.py tuning
- `--n-explain` / `--background-size` -- shap_hotspots.py tuning
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
  full 20^6 sweep (from `enumerate_cloop.py`), plus `analysis/` figures (from
  `analyze_enumeration.py`)
- `logs/` -- per-step stdout, so a multi-day cluster run leaves a trail
