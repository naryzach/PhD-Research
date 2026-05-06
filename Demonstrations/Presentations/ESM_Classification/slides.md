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

<!-- TODO: Populate with data pipeline details from ESM2 Classification/ -->

---

# Machine Learning Classification Framework

## Constructing Layered ESM Representations
- Extract per-residue embeddings from ESM models
- Build layered feature representations

## Implementation of Multi-Target Models
- Train classifiers on known binder/non-binder labels
- Multi-target prediction heads

<!-- TODO: Populate with model architecture details -->

---

# Predictive Binding Pipeline

## Benchmarking ML Sorting Against Empirical Datasets

- Cross-validate predictions against wet-lab results
- Compare ML-flagged candidates vs. generative candidates

<!-- TODO: Populate with benchmark results -->

---

# Conclusions

- ESM-based classification provides a data-driven complement to generative design
- Predictive models leverage evolutionary information for binding prediction
- Comparative framework validates ML categorizations against *de novo* outputs
