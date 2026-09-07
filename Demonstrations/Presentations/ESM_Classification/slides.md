---
theme: seriph
background: https://cover.sli.dev
title: "ESM-Driven Classification and Prediction of Known Binders"
info: |
  ## ESM Classification
  Experimental pipeline: ESM-driven classification and prediction of known binders.
  
  Author: Ryan Gustafson
class: text-center
drawings:
  persist: false
transition: slide-left
mdc: true
---

# ESM-Driven Classification and Prediction of Known Binders

Experimental Pipeline

<div class="abs-br m-6 flex gap-2">
  <span class="text-sm opacity-50">Ryan Gustafson</span>
</div>

---

# Introduction to Evolutionary Scale Modeling

- Complementary approach to *de novo* generation
- Exploits machine learning classification from established laboratory data
- ESM and ESM3 sequence representations

---

# Data Assimilation & Curation

## Preprocessing Known Laboratory Binders

- Collect and standardize existing binding data
- Encode sequences into ESM-compatible formats

### ESM-C Multi-Task Data Pipeline (`data_prep.py`)
- **Raw Data Curation:** Processes TIMP3 binding results (`TIMP_binder_all.csv`) containing variable 6-aa grafted loops on the constant 188-aa scaffold.
- **Conflict Resolution:** Resolves conflicting sequences (162 instances) using count-weighted majority vote, dropping exact ties.
- **Leakage-Safe Splits:** Implements exact loop-group splitting (15% test, 15% val, 70% train) to ensure evaluation data remains completely unseen.
- **Graph Percolation Guard:** Disabled Hamming-1 loop clustering (`cluster_hamming1: false`). Since edit-distance $\le 1$ loops form a giant component (~58% of loops), Hamming-1 clustering collapsed splits and distorted ratios.
- **Exposure Diagnostics:** Added Hamming-1 exposure metrics to monitor the fraction of validation/test loops within 1 edit of train data.

---

# Machine Learning Classification Framework

## Constructing Layered ESM Representations
- Extract per-residue embeddings from ESM models
- Build layered feature representations

## Implementation of Multi-Target Models
- Train classifiers on known binder/non-binder labels
- Multi-target prediction heads

### ESM-C Fine-Tuning Setup (`train.py`)
- **Base Architecture:** evolutionaryscale/ESM-C (via Hugging Face `Synthyra/ESMplusplus_small`, 300M params).
- <span class="text-amber-400">**Provenance note (2026-09-05):**</span> this specific run used the *small* ESM-C backbone — a real, since-flagged pilot-data risk (see `ESM_Classification.tex` §"ESM-C Fine-Tuning for TIMP3 Binder Classification" and `DATA_PROVENANCE.md`). The current production pipeline (`Local/esmc_multirun/`) retrains all targets on the confirmed **`Synthyra/ESMplusplus_large`** (600M) backbone across 5 data-slice variants — see updated results below.
- **Target-Specific Heads:** Independent binary heads for ADAM17, MMP3, and MMP9.
- **Sequence Pooling:** Loop-token pooling (pooling only the 6-aa variable loop tokens, preventing the constant scaffold from washing out the signal).
- **Loss Function:** Masked Binary Cross-Entropy (BCE) loss so each sample only supervises the target head it was screened against.
- **Loss Weighting:** Samples weighted by $log1p(count)$ to incorporate read-count confidence.
- **Staged Curriculum:**
  - **Phase 1 (3 epochs):** Random heads un-frozen, backbone frozen (lr = 1e-3).
  - **Phase 2 (15 epochs):** Top 6 transformer blocks un-frozen (lr = 1e-5).

---

# Predictive Binding Pipeline

## Benchmarking ML Sorting Against Empirical Datasets

- Cross-validate predictions against wet-lab results
- Compare ML-flagged candidates vs. generative candidates

### Held-Out Evaluation Performance
<span class="text-amber-400 text-[10px]">**Updated 2026-09-05:** the table below was an early, since-superseded snapshot of a single variant (`all3_original`) and used raw PR-AUC, which overstates MMP3 (see below). Replaced with the current, confirmed-`ESMplusplus_large`, best-per-target values across all 5 trained variants (source: `Local/Daily_Brief/2026-08-14/esmc_variant_model_selection_20260814_stats.csv`; full 15-row table in `ESM_Classification.tex` Table `tab:esmc_variants_current`).</span>
- **Validation Selection:** Evaluated on MCC and PR-lift (PR-AUC minus the positive-rate floor) rather than raw PR-AUC, which inflates MMP3 under its 93% positive rate.
- **Best-variant-per-target, current data:**

| Target Protease | Best Variant | Test Size (N) | Pos. Rate | ROC-AUC | MCC | PR-lift |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **MMP9** | `everything_combined` | 8,059 | 14.8% | **0.955** | **0.721** | **0.697** |
| **ADAM17** | `cloop_only` | 348 | 30.7% | **0.750** | **0.329** | **0.263** |
| **MMP3** | `everything_combined` | 5,610 | 92.4% | 0.804 | **0.252** | 0.055 |

- **Key Takeaways:** MMP9 is the one consistently strong target (MCC ≈0.72). ADAM17 is modest but real everywhere (MCC ≈0.27–0.33). **MMP3's headline raw PR-AUC (≈0.98, not shown above) is a base-rate mirage** driven by its ~93% positive rate — true PR-lift is only 0.05–0.06 and MCC 0.19–0.25 regardless of variant; MMP3 predictions from this classifier should not be trusted for ranking until that class is rebalanced.

---

# ESM-C Binder Enumeration & Motif Analysis

## Cross-Target Selectivity & Motif Architecture (June 29)

- **Target Overlap:** MMP3 and MMP9 share 38.1% of their top predicted binder sequences (Jaccard index 0.381). ADAM17 predictions are highly distinct (Jaccard overlap of 14.7% with MMP3 and 6.6% with MMP9).
- **Sticky MMP3 Subset:** Among the top 50,000 predicted binders, MMP9 is assigned 44,235 winning-target calls while MMP3 receives 0, indicating MMP3 hits represent a sticky subset of MMP9 hits.
- **Consensus Loop Motifs:**
  - **ADAM17:** `LSSDTT`
  - **MMP3:** `LSPDTT`
  - **MMP9:** `LSPTTL`
- **Mutational Insights:** Specificity is dictated by three positions on a shared `LS-x-x-T` core scaffold (positions 1, 2, 5). P3 separates ADAM17 (S) from MMP3/9 (P), and P4/P6 separate MMP9 (T, L) from MMP3 (D, T).
- **Novelty:** 100% of the top 300 deduped candidate loops are completely absent from the classifier training data.
- <span class="text-amber-400 text-[10px]">**Current re-run (`all3_original`, 2026-09-05) validates this:** Jaccard MMP9-MMP3 0.368, ADAM17-MMP3 0.143, ADAM17-MMP9 0.063 (near-identical); consensus motifs MMP3 `LSPDTT` and MMP9 `LSPTTL` unchanged, ADAM17 now `LPSDTT` (small shift from `LSSDTT`, position 2-3). The original June 29 source file no longer exists, so this current run is the traceable citation going forward.</span>

---

# Fine-Tuned ESM-C Classifier Evaluation

## Held-Out Test Performance & Validation Gaps (June 30)

- **MMP9 Test Performance:** Evaluation on a test set of 6,788 sequences (381 positives, 5.61% positive class balance) shows highly robust classification results:
  - **ROC-AUC:** `0.912` (high discriminative power)
  - **PR-AUC:** `0.746`
  - **MCC:** `0.722`
  - **F1-Score:** `0.734`
  - **Optimal Classification Threshold:** `0.872`
- **Validation Gaps & Limitations:**
  - **Single Target Focus:** Quantitative performance validation is currently restricted to MMP9 due to data availability.
  - **Strict Novel-Loop Challenge:** The out-of-distribution test set contains 0 positives (out of 152 sequences), leaving out-of-distribution generalization unvalidated.
- **Next Steps:** Future wet-lab binding screening is required to obtain labeled data for ADAM17, MMP3, and novel out-of-distribution loop designs.
- <span class="text-amber-400 text-[10px]">**Provenance note (2026-09-05):** this classifier is the small-model (`ESMplusplus_small`) pilot run flagged above — not an independent confirmation. The confirmed-large-model equivalent on the same source data (`mmp9_other` variant) scores ROC-AUC 0.920, MCC 0.716, PR-lift 0.695 on the same 6,788-sequence test set.</span>

---
layout: default
transition: fade-out
---

# ESM-C Classifier Confidence Distributions (July 4)

<div class="grid grid-cols-2 gap-8 mt-6 items-start">
  <div class="space-y-4">
    <div class="p-4 bg-white/5 rounded border border-white/10">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-3">Model Prediction Histograms</h4>
      <p class="text-[11px] leading-relaxed opacity-75">
        Probability distribution histograms comparing predicted binding confidence of the top 50,000 enumerated loop candidates across three targets.
      </p>
      <ul class="text-[10px] list-disc pl-4 space-y-1.5 opacity-70 mt-3">
        <li><b>Superseded 2026-09-05:</b> the original July 4 source file no longer exists and this "MMP9 saturates, others flat" claim could not be re-verified (see `DATA_PROVENANCE.md`).</li>
        <li><b>Current, real data (`all3_original`, confirmed `ESMplusplus_large`):</b> all three targets' retained top-50,000 sets cluster near 1.0 by construction, but at visibly different scales — ADAM17 lowest (0.994–0.997), MMP9 middle (0.996–0.998), MMP3 highest (0.998–1.000).</li>
      </ul>
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px]">
      <b class="text-amber-400">Calibration Insights (revised):</b> the different probability scales per target are evidence of per-target calibration skew (MMP3's classifier, trained on ~93% positive data, pushes probabilities closer to 1.0 regardless of true binding quality) — not an MMP9-specific saturation effect as originally claimed.
    </div>
  </div>
  <div class="space-y-3 flex flex-col items-center">
    <img src="../../Papers/Resources/ESMC_General/esmc_target_probability_distributions.png"
         alt="ESM-C Target Probability Distributions (current, 2026-09-05)"
         class="w-full rounded-lg border border-white/10 shadow-sm"
         style="max-height: 280px; object-fit: contain;">
    <p class="text-[9px] opacity-40 italic text-center">
      Current, traceable distribution (2026-09-05): confirmed-`ESMplusplus_large` `all3_original` classifier, top-50,000 C-loop designs per target.
    </p>
  </div>
</div>

---
layout: default
transition: fade-out
---

# ESM-C Classifier Evaluation (July 5)

<div class="grid grid-cols-2 gap-8 mt-6 items-start">
  <div class="space-y-4">
    <div class="p-4 bg-white/5 rounded border border-white/10">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-3">Held-Out Test Performance</h4>
      <p class="text-[11px] leading-relaxed opacity-75">
        Performance evaluation of the ESM-C classifier on held-out test data, showing metrics across targets (ROC-AUC, PR-AUC, MCC, and F1-score).
      </p>
      <ul class="text-[10px] list-disc pl-4 space-y-1.5 opacity-70 mt-3">
        <li><b>Updated 2026-09-05:</b> the ROC-AUC = 0.91 figure below was the single small-model (`ESMplusplus_small`) MMP9-only run, since superseded (see the Held-Out Evaluation Performance slide above). Current best MMP9 variant (confirmed `ESMplusplus_large`, `everything_combined`): ROC-AUC 0.955, MCC 0.721.</li>
        <li><b>Target Specificity:</b> MCC and PR-lift (not raw PR-AUC) are used to compare targets fairly across class-imbalance differences — see the current 5-variant chart at right.</li>
        <li><b>Pipeline Utility:</b> Deployment verdict: rank by the MMP9/`everything_combined` head; ADAM17 directionally useful but modest; do not trust MMP3 until rebalanced.</li>
      </ul>
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px]">
      <b class="text-amber-400">Next Steps:</b> Model-provenance check (confirmed all 5 current variants are `ESMplusplus_large`) and SHAP feature attribution are both complete — see `ESM_Classification.tex` §"Update: Model-Provenance Check, Final Performance, and Enumeration Refresh."
    </div>
  </div>
  <div class="space-y-3 flex flex-col items-center">
    <img src="../../Papers/Resources/ESMC_General/esmc_performance.png"
         alt="ESM-C Classifier Performance (current, 5-variant, 2026-09-04)"
         class="w-full rounded-lg border border-white/10 shadow-sm"
         style="max-height: 280px; object-fit: contain;">
    <p class="text-[9px] opacity-40 italic text-center">
      Current (2026-09-04) held-out performance (MCC, PR-lift) across all 5 confirmed-large-model variants, by target.
    </p>
  </div>
</div>


---
layout: default
transition: fade-out
---

# AlphaFold LpLDDT vs. MMP9 Binding Calibration (July 6)

<div class="grid grid-cols-2 gap-8 mt-6 items-start">
  <div class="space-y-4">
    <div class="p-4 bg-white/5 rounded border border-white/10">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-3">Structural Metric Calibration</h4>
      <p class="text-[11px] leading-relaxed opacity-75">
        Correlation analysis plotting AlphaFold structural confidence score (LpLDDT) against bench-measured flow cytometry binding (NMR) for MMP9.
      </p>
      <ul class="text-[10px] list-disc pl-4 space-y-1.5 opacity-70 mt-3">
        <li><b>Poor Correlation:</b> Correlation coefficient of $r = 0.19$, showing essentially no predictive relationship. <span class="text-amber-400">(Consistent with the companion De_Novo_Binder_Generation paper's current recompute: LpLDDT $\rho=0.20$, $p=0.23$, not statistically significant, $n=36$ — this slide's claim was not itself re-derived this pass but is not contradicted.)</span></li>
        <li><b>Implications for Design:</b> Optimizing sequence libraries solely for structural confidence does not yield higher binding affinities.</li>
        <li><b>Pipeline Complement:</b> Emphasizes the need for fine-tuned sequence classifiers (like ESM-C) to filter candidates.</li>
      </ul>
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px]">
      <b class="text-amber-400">Calibration Takeaway:</b> Structural confidence metrics from folding models are poor surrogates for actual binding activity on this target.
    </div>
  </div>
  <div class="space-y-3 flex flex-col items-center">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/af_calibration_LpLDDT_mmp9.png"
         alt="AlphaFold LpLDDT vs MMP9 Binding"
         class="w-full rounded-lg border border-white/10 shadow-sm"
         style="max-height: 280px; object-fit: contain;">
    <p class="text-[9px] opacity-40 italic text-center">
      AlphaFold LpLDDT score vs. measured MMP9 binding ratios (July 6).
    </p>
  </div>
</div>

---
layout: default
transition: fade-out
---

# ESM-C vs. FCS Binding Validation (July 8)

<div class="grid grid-cols-2 gap-8 mt-6 items-start">
  <div class="space-y-4">
    <div class="p-4 bg-white/5 rounded border border-white/10">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-3">Experimental Binding Validation</h4>
      <p class="text-[11px] leading-relaxed opacity-75">
        Evaluation of ESM-C predicted probabilities against actual flow-cytometry binding measurements for 11 matched constructs.
      </p>
      <ul class="text-[10px] list-disc pl-4 space-y-1.5 opacity-70 mt-3">
        <li><b>Not re-verified (2026-09-05):</b> ESM-C predicted probability was reported to correlate positively with raw binding (ADAM17 Pos Med Ratio: $\rho = 0.86^*$ on AB-loops; MMP9 Double+ %: $\rho = 0.62^*$ overall) — <b>but this run used the small ESM-C model (`ESMplusplus_small`)</b>, and the companion De_Novo_Binder_Generation paper's data-provenance audit separately flagged this exact claim as unconfirmed given the small $n$ (4–11) and how much the classifier landscape has changed since.</li>
        <li><b>Generalization claim:</b> reported at the time for AB-loop insertions ($\rho=0.86^*$ ADAM17, $\rho=0.75$ MMP9); not re-derived from current data this pass.</li>
        <li><b>Current basis for the deployment decision:</b> the confirmed-large-model held-out performance table (previous slides), not this unconfirmed wet-lab correlation.</li>
      </ul>
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px]">
      <b class="text-amber-400">Design Implications (revised):</b> This slide's historical correlation claim should not be cited as a current validated result. Structural confidence scores remain reasonable foldability/developability checks; the current ESM-C ranking decision rests on the confirmed-large-model MCC/PR-lift table, not on this figure.
    </div>
  </div>
  <div class="space-y-3 flex flex-col items-center">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/fig6_esmc_vs_fcs.png"
         alt="ESM-C vs FCS Binding Validation"
         class="w-full rounded-lg border border-white/10 shadow-sm"
         style="max-height: 280px; object-fit: contain;">
    <p class="text-[9px] opacity-40 italic text-center">
      ESM-C predicted probability vs. measured binding ratios (July 8).
    </p>
  </div>
</div>

---

# Conclusions

- ESM-based classification provides a data-driven complement to generative design
- Predictive models leverage evolutionary information for binding prediction
- Comparative framework validates ML categorizations against *de novo* outputs
- <span class="text-amber-400 text-[10px]">**Current status (2026-09-05):** 5-variant, confirmed-`ESMplusplus_large` comparison shows MMP9 strong (MCC≈0.72), ADAM17 modest (MCC≈0.33), and MMP3's high raw PR-AUC as a base-rate mirage (true PR-lift 0.05–0.06). Two open questions remain unresolved: an early small-vs-large model discrepancy, and a motif disagreement between two large-model MMP9 classifiers. The $\rho=0.86$ wet-lab correlation is unconfirmed and traces to the small model — see `ESM_Classification.tex` and `DATA_PROVENANCE.md`.</span>


