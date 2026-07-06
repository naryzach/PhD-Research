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
- **Validation Selection:** Evaluated on PR-AUC, ROC-AUC, MCC, and F1-score.
- **Test Metrics & Optimal Classification Thresholds:**

| Target Protease | Test Size (N) | Pos. Rate | PR-AUC | ROC-AUC | MCC | F1-Score | Optimal Threshold |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **ADAM17** | 347 | 31.4% | **0.517** | **0.712** | **0.266** | **0.538** | 0.246 |
| **MMP3** | 5,618 | 92.9% | **0.982** | **0.815** | **0.273** | **0.907** | 0.290 |
| **MMP9** | 1,231 | 66.1% | **0.869** | **0.791** | **0.442** | **0.824** | 0.341 |

- **Key Takeaways:** Excellent classification performance on MMP3 and MMP9 targets. Moderate predictive capability on ADAM17, which serves as a baseline for active de novo design selection.

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
        <li><b>MMP9 Domination:</b> Extremely high/tight distribution (mean = 0.999996, max = 1.0), defining a strict preference.</li>
        <li><b>MMP3 Modesty:</b> Flat and low-confidence distribution (mean = 0.5599, max = 0.6663).</li>
        <li><b>ADAM17 Flatness:</b> Extremely narrow/low distribution (mean = 0.5323, max = 0.5375).</li>
      </ul>
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px]">
      <b class="text-amber-400">Calibration Insights:</b> Because MMP9 confidence peaks near 1.0 while MMP3/ADAM17 remain near 0.5, selectivity metrics (like margin) are heavily skewed towards MMP9, requiring target-specific classification calibration.
    </div>
  </div>
  <div class="space-y-3 flex flex-col items-center">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/esmc_target_probability_distributions.png"
         alt="ESM-C Target Probability Distributions"
         class="w-full rounded-lg border border-white/10 shadow-sm"
         style="max-height: 280px; object-fit: contain;">
    <p class="text-[9px] opacity-40 italic text-center">
      Comparative predicted probability distributions for the top 50,000 designs per target (July 4).
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
        <li><b>High Classification Skill:</b> Excellent classification performance on MMP9 (ROC-AUC = 0.91), showing high discriminative power.</li>
        <li><b>Target Specificity:</b> Model metrics allow assessing predictive capability relative to class imbalances (PR-AUC vs. pos-rate).</li>
        <li><b>Pipeline Utility:</b> Validates using the fine-tuned classifier to sort and select TwistBio constructs.</li>
      </ul>
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px]">
      <b class="text-amber-400">Next Steps:</b> The baseline classification capability is established; feature attribution (SHAP plots) represents the next analytical extension.
    </div>
  </div>
  <div class="space-y-3 flex flex-col items-center">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/esmc_performance.png"
         alt="ESM-C Classifier Performance"
         class="w-full rounded-lg border border-white/10 shadow-sm"
         style="max-height: 280px; object-fit: contain;">
    <p class="text-[9px] opacity-40 italic text-center">
      Fine-tuned ESM-C classifier held-out test performance per target (July 5).
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
        <li><b>Poor Correlation:</b> Correlation coefficient of $r = 0.19$, showing essentially no predictive relationship.</li>
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

# Conclusions

- ESM-based classification provides a data-driven complement to generative design
- Predictive models leverage evolutionary information for binding prediction
- Comparative framework validates ML categorizations against *de novo* outputs


