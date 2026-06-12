# Construct Metric Summary

A comprehensive per-construct × per-target rollup of **every** FCS metric from the
`*_Renamed_Analysis` aggregation — the generalization of
`Aggregate_FCS_Analysis/Selectivity_Analysis/selectivity_summary.csv`, which only
rolled up the single selectivity metric.

## Source
Built from `Aggregate_FCS_Analysis/aggregate_summary.csv` (the consolidation of the
per-experiment `*_Renamed_Analysis/summary_stats.csv` files). Only **valid trials**
are included — rows with `Trial Failed == True` (failed positive control, low
expression, too few events) are dropped before aggregating.

Regenerate with `../Prediction_vs_Result_Analysis/analyze.py`.

## Files
| file | shape | description |
|------|-------|-------------|
| `construct_target_metric_summary.csv` | long | one row per (Construct, Target, Metric) with `Mean, Std, SEM, N` |
| `construct_target_metric_summary_wide.csv` | wide | one row per (Construct, Target) with the mean of each metric + `N_trials` |

Covers 15 constructs × 5 targets (ADAM10, ADAM17, MMP2, MMP3, MMP9) × 16 metrics.

## Caveats
- **ADAM10**: the only valid trials are the TIMP 3 positive control — there is **no
  usable ADAM10 data for any test construct**, so any ADAM10 / A17>A10 comparison is
  not yet verifiable.
- **N is small and uneven** (often 1–3) for ADAM17/MMP2 on some constructs; treat
  single-trial means as provisional.
- Inclusion differs slightly from `selectivity_summary.csv`, which used a different
  trial filter; this rollup keys strictly off the `Trial Failed` flag.
