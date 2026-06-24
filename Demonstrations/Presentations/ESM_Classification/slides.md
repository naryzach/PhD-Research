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

# Conclusions

- ESM-based classification provides a data-driven complement to generative design
- Predictive models leverage evolutionary information for binding prediction
- Comparative framework validates ML categorizations against *de novo* outputs
