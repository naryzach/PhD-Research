# Metric regression sweep — which predicted metric tracks which experimental readout

Every AlphaFold metric (ipTM, LpLDDT, PAE, ApTM, BpTM, pTM, pLDDT) regressed against
every experimental FCS readout (20 metrics: binding raw / binding normalized / binding
quality / expression / count), per target (ADAM17, MMP2, MMP9) and pooled
within-target-standardized. Built by `../regression_sweep.py`. Correlations are oriented
so **+ρ = predicts in the expected direction** (PAE is sign-flipped since lower = better).

n ≈ 12 constructs per target (36 pooled). This is **hypothesis-generating to inform the
blend**, not a definitive calibration — both sides have compressed dynamic range and the
per-target n is small.

## The result that reframes the ranking: dynamic range
A rank correlation is meaningless if the metric barely varies. `af_metric_dynamic_range.csv`:

| AF metric | within-target range | verdict |
|-----------|--------------------:|---------|
| **LpLDDT** | **48.5** (43–94) | real range |
| **PAE** | **18.7** (3–23) | real range |
| pLDDT | 4.0 (82–90) | modest |
| ipTM | 0.10 (0.77–0.90) | tight |
| pTM | 0.06 | ~flat |
| **ApTM** | 0.04 | **FLAT — rank = noise** |
| **BpTM** | 0.03 (0.84–0.93) | **FLAT — rank = noise** |

So the highest "correlations" in the naive ranking (BpTM/ApTM ρ≈0.40) are **rank-ordering
values that span 0.03–0.04** — almost certainly numerical noise, not signal. The only
metrics with real spread are **LpLDDT and PAE**, and they are near-perfectly redundant
(Spearman −0.97 → one axis). That one real-range axis carries only **weak** binding
signal (oriented ρ ≈ 0.15–0.21, not significant at n=36).

## Expression confound (`expression_vs_binding_confound.csv`)
| AF metric | Expr+ % | Pos Med Ratio (binding) | Norm Median Ratio (binding) |
|-----------|--------:|------------------------:|----------------------------:|
| **ipTM** | **+0.36** | −0.09 | −0.22 |
| pTM | +0.33 | 0.00 | −0.14 |
| BpTM | +0.17 | **+0.40** | +0.18 |
| LpLDDT | +0.14 | +0.20 | +0.15 |

**ipTM — the metric the pipeline optimizes — tracks expression, not binding** (and
slightly *anti*-tracks binding). Consistent across all 3 targets (ADAM17 +0.48, MMP2
+0.34, MMP9 +0.27). This is direct evidence for the foldability/expressibility-bias
concern: in this dataset ipTM is closer to "is this construct well-behaved / expressible"
than "does it bind." pTM behaves the same (it is 0.89-correlated with ipTM).

## What this says for an informed blend (toward BINDING)
- **Backbone: LpLDDT (≡ −PAE).** Only real-range metric, weak-but-positive and
  binding-leaning. Use one of the pair, not both (redundant).
- **Second, independent axis: BpTM.** Lowest redundancy with the LpLDDT/PAE cluster
  (Spearman 0.16 / −0.18), binding-leaning, and the only metric **positive across all 3
  targets** for `% of Pos Ctrl` / `Pos Med Ratio`. But it is FLAT → include at **low
  weight** and re-check on a wider-range dataset before trusting it.
- **De-weight for binding: ipTM, pTM, pLDDT** — expression-confounded and/or redundant
  with ipTM. (This is the opposite of the current ipTM-centric composite — worth an A/B.)

## Best experimental readout to calibrate against (`exp_readout_predictability.csv`)
Most predictable + consistent across targets: **`% of Pos Ctrl`, `Pos Med Ratio`,
`Bind FC vs Pos Ctrl`** (oriented ρ ≈ 0.40, positive in all 3 targets). These — binding
relative to a known positive control — are better ground-truth targets than
**`Norm Median Ratio`** (least predictable, best ρ only 0.18), which is what the earlier
analysis and the lab notebooks lead with. `% of Pos Ctrl` came from the per-folder
`summary_stats.csv` enrichment, so the enrichment paid off. `Binding Efficiency` and
`Bind Stain Index` are **sign-unstable** (invert on MMP2) → poor single targets.

## MMP2 anomaly (reproducible)
MMP2 repeatedly inverts sign vs ADAM17/MMP9 (binding efficiency, stain index, double+%).
This matches `filter_methods.md` §5.2's independent flag that ESMFold2-vs-AF3 inverted
for MMP2. **Two independent calibrations now flag MMP2** → real, not noise. Calibrate
per-target and treat MMP2 (assay gating? target prep? structure?) as a known special case.

## Files
| file | description |
|------|-------------|
| `sweep_long.csv` | every AF×Exp×scope correlation (Pearson + Spearman, oriented) |
| `af_predictor_ranking.csv` | AF metrics ranked, **with dynamic-range flag** |
| `exp_readout_predictability.csv` | each experimental readout's best AF predictor |
| `af_metric_dynamic_range.csv` | within-target spread per AF metric (the key caveat) |
| `expression_vs_binding_confound.csv` | binding vs expression signal per AF metric |
| `af_metric_intercorrelation.csv` | redundancy among AF metrics (for blend design) |
| `heatmap_pooledz.png` / `heatmap_ADAM17.png` / `heatmap_MMP2.png` | oriented-ρ heatmaps |

## Bottom line
No AlphaFold metric is a trustworthy *fine-ranker* of binding here: the real-range axis
(LpLDDT/PAE) is only weakly predictive, and the apparently-stronger metrics are flat.
For a blend, lead with LpLDDT/PAE, add BpTM as a low-weight independent term, de-weight
ipTM/pTM (expression-confounded), and calibrate per-target against `% of Pos Ctrl` /
`Pos Med Ratio`. Then redo this on a purpose-built ordered batch with wide predicted
range and an assay tuned for dynamic range — that is what turns this into a real
ground-truth calibration.
