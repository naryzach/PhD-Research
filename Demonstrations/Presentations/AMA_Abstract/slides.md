---
theme: seriph
background: https://cover.sli.dev
title: "Generative Deep Learning Pipeline for Patient-Specific Therapeutic Binders"
info: |
  ## AMA Abstract Presentation
  Generative deep learning pipeline for the rapid design of patient-specific therapeutic binders against customized cancer markers.
  
  Author: Ryan Gustafson
class: text-center
drawings:
  persist: false
transition: slide-left
mdc: true
---

# Generative Deep Learning Pipeline for Patient-Specific Therapeutic Binders

Against Customized Cancer Markers

<div class="abs-br m-6 flex gap-2">
  <span class="text-sm opacity-50">Ryan Gustafson</span>
</div>

---
transition: fade-out
---

# Importance

- Precision medicine demands therapeutic agents tailored to individual genomic variations
- Current approaches remain computationally expensive and error-prone
- Urgent need for accelerated pipelines capable of rapid generation and validation

---

# Objective

Develop and validate an **end-to-end computational generative framework** that systematically synthesizes and screens customized protein binders optimized for individualized clinical targets.

---

# Design & Methods

1. **RFdiffusion** — Hallucinate *de novo* structural backbones as molecular scaffolds
2. **ProteinMPNN** — Populate poly-glycine coordinates with specific amino acid identities
3. **AlphaFold 3 / AF2-multimer** — High-fidelity structural validation
4. **Custom metric extraction** — Loop pLDDT and Loop PAE confidence scoring

---

# Main Outcomes & Measures

- Computational stratification of sequence candidates
- Distinct binding profiles evaluated:
  - **General thermodynamic affinity**
  - **Overall target specificity**
  - **Complex pairwise selectivity**

---

# Results

- Z-score standardization across localized sequence evaluations
- Successfully filtered theoretical designs into high-confidence consensus set
- Eliminated malformed configurations (pTM < 0.8)
- Isolated *de novo* variants outperforming wild-type controls

---

# Conclusions

- Combining deep-learning structural diffusion with AlphaFold verification yields a **robust, repeatable blueprint** for custom molecular design
- Permits rapid conception and verification of tailored binders against evolving patient-specific genomic targets
- Accelerates timeline from **tumor marker discovery → viable therapeutic candidate**
