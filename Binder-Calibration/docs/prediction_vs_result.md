# Prediction vs Result Analysis — Dec-2025 Twist order

Compares the AlphaFold-driven predictions that selected the December 2025 Twist
library (`../Twist_Order_Dec2025`, with metrics in `../AlphaFoldMetrics`) against the
experimental FCS binding results (`../Aggregate_FCS_Analysis`).

Regenerate everything with `python3 analyze.py`.

## What was predicted vs what is testable
The AlphaFold pipeline scored each loop variant against **adam10, adam17, mmp2, mmp9**
and tagged constructs with an expected behavior (`Final_Ordering_List_*.xlsx`
`Reasons` → the `Expected` column in `cross_target_summary.csv`):

`High` / `Low` overall affinity, `A17+` (ADAM17 binder), `A17>A10`, `M9+`, `M9>M2`.

Experimentally usable targets right now: **ADAM17, MMP2, MMP9**.
- **ADAM10** has zero valid test-construct trials (only the TIMP 3 positive control
  passed QC) → every `A17>A10` prediction is **unverifiable**.
- **MMP9** ANOVA across constructs is non-significant (p = 0.084) → no between-construct
  spread, so `M9+` magnitude calls are not testable and `M9>M2` is only a weak
  within-construct direction check.

## Headline result
**The AlphaFold interface scores did not predict measured binding strength in this
round.** Across the three usable targets, predicted interface confidence vs measured
binding (Norm Median Ratio) is essentially flat:

| Target | Pearson r (LpLDDT vs NMR) | p |
|--------|---------------------------|---|
| ADAM17 | 0.07 | 0.83 |
| MMP2   | 0.22 | 0.48 |
| MMP9   | 0.19 | 0.56 |

`ipTM` is uninformative by construction (range only 0.77–0.90). The one *significant*
relationship is **inverse**: for MMP2, **Binding Efficiency** falls as predicted
interface quality rises (ipTM r = −0.65, LpLDDT r = −0.61, ApTM r = −0.81; all five
interface metrics agree). This persists within the Exp-only subset (ipTM r = −0.79,
p = 0.01), so it is **not** an Exp-vs-PMPNN redesign batch artifact. Net: no robust,
consistent, positive numerical predictivity — depending on the readout the sign even
flips. See `correlation_summary.csv` and `plot_LpLDDT_vs_binding.png`.

Both sides have a **compressed dynamic range**: predicted ipTM is nearly constant, and
the assay readouts cluster between 0.6× and 1.26× the TIMP 3 positive control — there
is little signal to correlate against.

## The dominant signal is construct-level "stickiness", NOT target selectivity
Across the 12 constructs, binding to the three targets is **strongly correlated**
(ADAM17–MMP2 r=+0.72 p=0.009; ADAM17–MMP9 r=+0.69 p=0.013; MMP2–MMP9 r=+0.50). A
construct that binds one protease tends to bind all of them. This single construct-level
factor (general avidity / display level — likely tied to expression; see the
ipTM↔expression result in `Regression_Sweep/`) dwarfs any target-specific or
predicted-selectivity signal, and is stronger than any AF-metric-vs-binding correlation.

Consequences for the categorical calls (`prediction_scorecard.csv`,
`plot_overall_affinity_by_target.png`):
- **A17+ designs bind ADAM17 well but are NOT selective.** AB 4, AB 7, C 11, C 14 are
  among the strongest ADAM17 binders — but they are *also* the strongest **MMP2** binders
  (AB 4 and C 14 top the MMP2 ranking). ADAM17 only exceeds the best off-target by
  +0.02–0.08 NMR (within noise). So "binds ADAM17" is loosely true; "ADAM17-selective"
  is **not** supported. (An earlier version of this file wrongly stated these bind MMP2
  *less* than MMP-directed designs — that was a sign error; they bind MMP2 *more*.)
- **M9>M2 direction: 3/4 within-construct** (AB 1, AB 6, C 12 higher MMP9 than MMP2;
  C 15 not) — but AB 1/AB 6 are among the *weakest* binders overall, MMP9 spread is n.s.,
  and margins are tiny, so this is very weak support, not validated selectivity.

## Where they failed
- **Overall-affinity axis did not pan out.** The predicted **`Low`** designs bound as
  well as everything else — **C 13 ranked 4/12** and AB 5 ranked 6/12 by overall
  binding (a clear miss). The predicted **`High`** design **AB 3 was mid-pack (7/12)**;
  the other `High` design (C 16) failed QC, so no data.
- No usable validation for `M9+` (MMP9 flat) or any `A17>A10` (no ADAM10 data).

## Files
| file | description |
|------|-------------|
| `predicted_vs_actual.csv` | per construct: all AF metrics × 4 targets + experimental means (NMR, IWB, BindEff, Double+%) for ADAM17/MMP2/MMP9 |
| `correlation_summary.csv` | Pearson/Spearman for each AF-metric × exp-metric × target, plus pooled |
| `prediction_scorecard.csv` | per construct: expected label → Hit / Partial / Miss / Untestable, with the numbers |
| `plot_LpLDDT_vs_binding.png` | scatter of interface pLDDT vs Norm Median Ratio per target |
| `plot_overall_affinity_rank.png` | constructs ranked by measured binding, colored by predicted tier |
| `analyze.py` | reproducible script (also builds `../Construct_Metric_Summary`) |

## Caveats
- The ordered library mixes original **Exp** designs and **ProteinMPNN-redesigned**
  sequences (AB 1, AB 3, C 11, C 16). Scores are on the same AlphaFold scale and the
  two groups have near-identical mean binding, but the flavor is recorded in
  `AF_Source` so any batch effect can be checked.
- Per-construct trial counts are small (often n = 1–3 for ADAM17/MMP2); means are noisy.
- Primary experimental readout is **Norm Median Ratio** (binding:expression normalized
  to the TIMP 3 positive control); IWB index and binding efficiency are included for
  cross-checking and tell a similar (non-predictive) story.
