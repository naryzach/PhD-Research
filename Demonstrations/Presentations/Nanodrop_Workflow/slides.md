---
theme: seriph
background: https://cover.sli.dev
title: "Integration of Nanodrop Quantification Data"
info: |
  ## Nanodrop Workflow
  Standalone workflow: Integration of Nanodrop quantification data.
  
  Author: Ryan Gustafson
class: text-center
drawings:
  persist: false
transition: slide-left
mdc: true
---

# Integration of Nanodrop Quantification Data

Standalone Workflow

<div class="abs-br m-6 flex gap-2">
  <span class="text-sm opacity-50">Ryan Gustafson</span>
</div>

---

# Introduction and Baseline Quantification

- Nanodrop spectroscopic data validates DNA purity
- Standardizes nucleic acid concentrations early in sample preparation
- Ensures experimental consistency across workflows

---

# Automated Spectroscopic Data Ingestion

## Direct Parser Scripting
- Automated ingestion of Nanodrop numerical summaries
- Organizes absorption data and contamination estimates

<!-- TODO: Populate with parser details from DNA Analysis/nanodrop_parser.py -->

---

# Quality Control and Benchmarking

## Assessing Purity and Contamination Curves
- 260/280 and 260/230 ratio analysis
- Threshold-based quality gating

<!-- TODO: Populate with QC results -->

---

# Impact on Downstream Standardization

- Normalized concentrations feed directly into library preparation
- Mitigates human bottlenecking in sample QC workflows
