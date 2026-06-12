# Ground-truth calibration & the multi-term prediction recipe (v1)

Calibrates AlphaFold metrics against experimental FCS binding (the Dec-2025 batch:
12 constructs × 3 usable targets ADAM17/MMP2/MMP9) and turns the result into an
importable, multi-factor scoring recipe for the generation pipeline.

- `calibration.py` — runs the three stages below, writes the CSVs + plot. Default source
  is the AF-server matrices (== AF3); `--esm-scores <csv>` calibrates ESMFold2 instead.
- `../../Generation/binding_recipe.py` — the recipe (config weights + `score_design()`),
  in the tracked `Generation/` dir so the pipeline can import it (`Local/` is gitignored).
- `make_esmfold2_manifest.py` — builds the cluster input manifest for the ESMFold2 run.
- Re-run with `python calibration.py` whenever a new ordered batch's FCS results land
  (add the constructs to `CONSTRUCTS`, point `AGGF` at the updated aggregate).

> **This is v1 on n=12.** It is unsolved-problem territory: no metric is close to
> perfect. The recipe is deliberately multi-factor, every term is justified by both
> mechanism and the calibration below, and it is built to be re-fit as data grows.

## Stage 1 — where the experimental signal lives (`variance_decomposition.csv`)
Decomposing consensus binding (constructs × targets) into main effects + interaction:

| component | fraction of variance |
|-----------|---------------------:|
| **construct (stickiness / avidity)** | **49%** |
| target (assay baseline) | 17% |
| **interaction (target-specific binding)** | **33%** |

Half the signal is a non-specific, construct-level "stickiness" factor (a construct
that binds one protease tends to bind all — cross-target r≈0.7). Only a third is
target-specific. **You must calibrate against the target-specific part**, or you just
predict avidity. This is *the* reason raw-binding correlations looked hopeless earlier.

## Stage 2 — what AF actually predicts (`afmetric_calibration.csv`, `selectivity_calibration.csv`)
Once stickiness is removed (double-centred interaction), AF metrics predict
**target-specific binding far better than raw binding** (oriented Spearman):

| AF metric | predicts stickiness | **predicts target-specific binding** |
|-----------|--------------------:|-------------------------------------:|
| PAE | 0.15 | **0.50** |
| LpLDDT | 0.22 | **0.47** |
| BpTM | 0.30 | **0.44** |
| ipTM | −0.22 | 0.41 |
| pTM | −0.22 | 0.38 |
| pLDDT | 0.04 | 0.30 |

vs ≈0.15 on raw binding. The signal was there, masked by stickiness. Note ipTM/pTM
*anti*-track stickiness — consistent with ipTM behaving as an expression/foldability
signal rather than a binding one.

**Selectivity (on−off gap)** is the weak link and we report it honestly: the only
high gap-correlations land on flat metrics (ApTM, n=9 — suspect), while the real-range
ADAM17>MMP2 gap is ≈0. So selectivity is predictable in part (ADAM17>MMP9) but not
reliably (ADAM17>MMP2). → include it as a term, not as a load-bearing factor.

## Stage 3 — the recipe and how it scores (`recipe_validation.csv`, `weight_sensitivity.csv`)
`binding_recipe.py` combines four terms (weights = reasoned v1, normalized):

| term | weight | metric(s) | why |
|------|-------:|-----------|-----|
| **affinity** | 0.40 | LpLDDT & −PAE (mean) | only AF metrics with real dynamic range; backbone |
| **selectivity** | 0.30 | on − mean(off) gap on that confidence | the project goal is *selective* binders; the off-target term — present but not lone |
| **independent** | 0.12 | BpTM | most independent axis (ρ≈0.16 vs LpLDDT), binding-leaning, but flat → low weight, flagged |
| **expressibility** | 0.18 | ipTM | tracks expression not binding → a developability prior, kept small |

**Performance (Spearman vs experimental binding, ranking within a target):**

| scorer | ρ | p |
|--------|--:|--:|
| **RECIPE (multi-term v1)** | **0.31** | 0.07 |
| LpLDDT only | 0.21 | 0.23 |
| ipTM only (current pipeline objective) | **−0.11** | 0.53 |
| data-fitted blend (leave-1-construct-out) | 0.20 | 0.23 |

**Read this with the small-n in mind.** The recipe ρ=0.31 has a **95% bootstrap CI of
[−0.15, +0.71]** (resampling constructs, n=12) — it crosses zero. The recipe is
*directionally* better than every single metric and beats a data-fitted blend, but at
this n we **cannot** distinguish it from zero (or from quite strong). Treat it as a
reasoned prior, not a measured coefficient; the CI tightens only with more constructs.

Two things matter here:
1. The multi-term recipe beats every single metric — and beats **ipTM**, which the
   generation pipeline currently optimizes (ipTM is *worse than random* on this data).
2. The reasoned recipe (0.31) beats a **data-fitted** blend cross-validated
   leave-one-construct-out (0.20). At n=12, least-squares weights overfit and
   generalize *worse* — so weights are set by reasoning, not regression, on purpose.

Robustness (`weight_sensitivity.csv`): variants span 0.28–0.32, so the result is not
knife-edge. `v4_no_iptm` ≈ `v1`, i.e. the data is indifferent to ipTM for *binding* —
it is retained only as a developability prior. Selectivity: recipe gap predicts
ADAM17>MMP9 (ρ=0.52) but not ADAM17>MMP2 (≈0) — the known MMP2 anomaly.

## Using it in the generation pipeline
The recipe lives in **`Generation/binding_recipe.py`** (tracked; `Local/` is gitignored).
It is **robust to missing metrics** — any term whose input is absent is dropped and the
remaining weights renormalized, so the same call works on AF3 (no BpTM) and ESMFold2
(ipTM/pTM/pLDDT only).
```python
from binding_recipe import score_design
# on_metrics: AF metrics for the design vs its ON-target
# off_list:   one metrics-dict per OFF-target protease (use ALL of them, not just tested)
s = score_design(on_metrics, off_list)
```
It is already wired into **`select_binders_to_order.py`** as an optional `recipe_score`
column carried *alongside* `composite` (not replacing it). Note: the main ordering table
folds each design vs one target, so the **selectivity term only fires when cross-target
folds are present** (a `--specificity-scores` / crossfold run); otherwise it auto-drops.

**Predictor-scale caveat.** Calibrated on **AF-server** metrics. **AF-server and AF3 are
the same model run in different places → this calibration is directly valid for AF3.**
**ESMFold2 is a different scale** and does not natively emit loop-pLDDT or interface PAE,
so re-derive `NORM_RANGES` and re-run the calibration for ESMFold2 (next section) before
trusting it there.

## ESMFold2 calibration (the pipeline's mass-prediction ranker)
ESMFold2 is what will rank designs en masse, so it needs its own calibration against the
same FCS data. **It can't be run in this WSL env (no conda) — run it on the cluster.**

1. `python make_esmfold2_manifest.py` → `esmfold2_inputs/esmfold2_manifest.csv`
   (13 constructs × 6 targets: `design_id, target_name, binder_seq, target_seq`). Binder
   = the expressed twist construct; `design_id` is the construct name so the join is by
   name, not sequence.
2. Run ESMFold2 on that manifest on the cluster (see `Generation/score_with_esmfold2.py`).
   **Also emit `esm_lplddt` (loop pLDDT) and `esm_pae` (interface PAE)** — its default
   scorer only saves whole-binder pLDDT/ipTM/pTM, but `SAVE_ESMFOLD2_STRUCTURES=True`
   keeps the predicted CIF, so derive them with `loop_plddt()` / `loop_interface_pae()`
   from `select_binders_to_order.py`. Without those two, the recipe's affinity backbone
   runs in reduced (pLDDT-only) form.
3. Produce a scores CSV `Construct,Target,esm_iptm,esm_ptm,esm_plddt[,esm_lplddt,esm_pae]`
   and run `python calibration.py --esm-scores <csv>`. This re-runs Stages 2–3 on ESMFold2
   metrics and writes the same outputs suffixed `_esmfold2`, so the AF3 and ESMFold2
   calibrations sit side by side and can be compared directly.

## Honest limitations
- **n=12 constructs; the headline ρ=0.31 has a 95% CI of [−0.15, +0.71]** — directional,
  not a measured coefficient. Don't over-fit decisions to it.
- Per-target weights are **not** fit (hopeless at this n) — global weights, MMP2 flagged.
- Selectivity validation is inconclusive for ADAM17>MMP2.
- ~half the binding variance (stickiness) is not predicted by any AF metric; an explicit
  expression/developability predictor would help and is the natural next ingredient.
- Calibrated on AF-server/AF3 metrics; the ESMFold2 calibration is set up but not yet run
  (needs the cluster) — until it is, don't assume the AF3 weights transfer to ESMFold2.

## Re-fit checklist when new results arrive
1. Add new constructs+sequences to `CONSTRUCTS` (and `AGGF` if the aggregate moved).
2. `python calibration.py` → re-check the decomposition, Stage-2 table, and whether the
   recipe still beats baselines and the data-fitted blend.
3. Only once n is comfortably larger (≳30–40 constructs) consider letting the LOO
   data-fitted weights replace the reasoned ones, and revisit per-target weighting.
