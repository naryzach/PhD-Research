---
theme: seriph
title: "De Novo Protein Prediction to Binding Evaluation"
info: |
  ## De Novo Binder Generation Pipeline
  TIMP3 loop engineering via RFdiffusion → ProteinMPNN → AlphaFold → Flow Cytometry

  **Author:** Ryan Gustafson
class: text-center
drawings:
  persist: false
transition: slide-left
mdc: true
config:
  nav: true
background: "#0f172a"
layout: center
---

# De Novo Protein Prediction <br/> to Binding Evaluation

Engineering TIMP3 variants with generative AI for targeted metalloproteinase inhibition

<div class="abs-br m-6 flex gap-2">
  <span class="text-sm opacity-50">Ryan Gustafson · PhD Research</span>
</div>

<style>
h1 {
  background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%);
  -webkit-background-clip: text;
  -webkit-text-fill-color: transparent;
  font-size: 2.8rem !important;
  font-weight: 800;
}
</style>

---
transition: fade-out
layout: default
---

# Pipeline Overview

<PipelineFlow />

<div class="text-center text-[11px] opacity-40 mt-4">
  Hover each stage for details
</div>

---

# Introduction — TIMP3 & Target Engineering

<div class="grid grid-cols-2 gap-12 mt-8 text-[13px] items-start max-w-5xl mx-auto overflow-hidden">
<div>

## The Challenge
- **TIMP3** (Tissue Inhibitor of Metalloproteinases-3) is a natural regulator of MMPs and ADAMs
- Engineering variants with **altered specificity** enables selective therapeutic intervention
- Wet-lab screening of random mutants is slow and expensive

## Our Approach
- **Computational-first** pipeline: generate → verify → synthesize → test
- Loop-focused redesign: target binding interface directly
- Thousands of candidates → consensus filtering → top hits

</div>
<div>

## Targets
| Target | Type | Role |
|--------|------|------|
| **ADAM10** | Disintegrin | Notch signaling |
| **ADAM17** | Disintegrin | TNF-α shedding |
| **MMP2** | Gelatinase | ECM degradation |
| **MMP9** | Gelatinase | Tumor invasion |

<div class="mt-4 text-xs opacity-60">

All targets modeled as catalytic domain (cd) complexes with TIMP3

</div>

</div>
</div>

---

# TIMP3 Loop Architecture

<LoopConfig />

<div class="grid grid-cols-3 gap-4 mt-2 text-xs">
<div class="p-3 rounded-lg" style="background:rgba(79,172,254,0.08); border:1px solid rgba(79,172,254,0.2)">

**AB Loop** (pos 30–36) — Primary binding interface. 6 aa native, expandable to 15 aa. Flanks: `LVK...LVY`

</div>
<div class="p-3 rounded-lg" style="background:rgba(79,172,254,0.08); border:1px solid rgba(79,172,254,0.2)">

**C Loop** (pos 62–68) — Secondary contact. 6 aa native, expandable to 15 aa. Flanks: `HTE...GLK`

</div>
<div class="p-3 rounded-lg" style="background:rgba(79,172,254,0.08); border:1px solid rgba(79,172,254,0.2)">

**EF Loop** (pos 92–96) — Tertiary contact. 4 aa native, expandable to 10 aa. Flanks: `MYT...FVE`

</div>
</div>

<div class="text-center text-xs opacity-40 mt-2">
  Hover loops on the diagram for detailed parameters · Data from RFd_Batch.py loop_configs
</div>

---
layout: two-cols
---

# Stage 1: RFdiffusion

<div class="text-sm">

## Backbone Generation
- Hallucinate novel 3D backbone coordinates
- Input: TIMP3–target complex PDB
- Loop regions replaced with **poly-glycine** placeholders
- Diffuser steps: T=20 (short sequence mode)

## Generation Parameters
```python
loop_configs = {
  "AB": { normal: 6,  max: 15, pos: 30 },
  "C":  { normal: 6,  max: 15, pos: 62 },
  "EF": { normal: 4,  max: 10, pos: 92 },
}
num_designs = 25  # per target complex
```

## Infrastructure
- **A30 GPU** node on HPC cluster
- ~20 hours per loop expansion run

</div>

::right::

<div class="ml-4 text-sm">

## Contig String Assembly

The contig string tells RFdiffusion which regions to keep fixed vs. regenerate:

```
A1-30/6-15/A37-62/6-15/A69-92/4-10/
A97-121/0 B1-{target_len}
```

- Fixed scaffold: `A1-30` (N-terminus)
- Variable loop: `6-15` (AB loop expansion)  
- Fixed region: `A37-62`
- Variable loop: `6-15` (C loop expansion)
- Chain break: `/0`
- Fixed target: `B1-{len}`

## Output
Each run produces **25 PDB files** with novel poly-glycine loop geometries per target complex

</div>

---
layout: center
---

# Stage 2: ProteinMPNN

<div class="grid grid-cols-2 gap-10 mt-6 text-[13px] items-start max-w-6xl mx-auto overflow-hidden">
<div>

## Design Strategy
- Input: RFdiffusion backbone ensemble
- **Fixed-backbone** sequence design
- 8 sequences designed per backbone structure
- All scaffold residues frozen — only loop residues redesigned

## Key Parameters
- Sampling Temp: **0.1**
- Model: `v_48_020`
- Target: Multimer-aware packing

</div>
<div class="p-4 rounded-lg bg-white bg-opacity-5 border border-white border-opacity-10">

<div class="text-[11px] font-bold mb-2">Top Predicted Loops (ADAM10)</div>

<div class="text-[9px] leading-tight">

| Len | Top Loop | Score | Rec |
|---|---|---|---|
| 6 | `KGEIDE` | 0.857 | 16% |
| 7 | `KKEGYEE` | 1.095 | 14% |
| 8 | `KLPDGFEE` | 0.895 | 12% |
| 9 | `KDEKTGFET` | 0.885 | 11% |
| 10 | `VLDEKTGLKE` | 0.992 | 10% |
| 11 | `TVKTPGATFKE` | 1.041 | 9% |
| 12 | `KVKDWGGGEFEE` | 1.037 | 25% |
| 14 | `EVTSPEDPSVKFKE` | 1.067 | 0% |
| 15 | `TKTVTENGETAKVEE` | 1.024 | 6% |

</div>

<div class="mt-3 text-[9px] opacity-50 italic">
Lower score = better packing. Data: Results/best_loops_per_length_ADAM10.csv
</div>

</div>
</div>

---

# Stage 3: AlphaFold Verification

<div class="grid grid-cols-2 gap-6 text-sm">
<div>

## Blind Structural Prediction
- **AlphaFold 3** and **AF2-multimer** fold each designed TIMP3 variant against each target
- No spatial biases from the design phase — purely predictive
- Extract per-variant confidence metrics from JSON outputs

## Quality Gate
```python
# Discard misfolded configurations
if variant.pTM < 0.8:
    discard(variant)
```

## Extracted Metrics
- **pTM** — Overall predicted TM-score
- **ipTM** — Interface-specific predicted TM-score
- **Loop pLDDT** — Per-residue confidence for AB/C/EF loops
- **Interface PAE** — Predicted aligned error at binding interface

</div>
<div>

## Z-Score Normalization

Raw metrics are standardized across all candidates:

$$Z_i = \frac{x_i - \mu}{\sigma}$$

This identifies variants significantly outperforming the mean, independent of absolute metric scale.

### Key Insight
> Z-scores enable **cross-metric** consensus filtering: a variant must score well on **all** metrics simultaneously, not just one.

</div>
</div>

---

# AlphaFold Results — Z-Score Analysis

<ZScoreTable />

<div class="mt-4 text-sm">

### Key Findings
- **`wt`** and **`asesnc`** show consistently strong Z-scores across all metrics for ADAM10
- **`asesta`** and **`agestc`** show ADAM17-specific Loop pLDDT advantages
- **`agesna`** uniquely shows MMP9 specificity in ipTM

</div>

---

# Consensus Filtering & Variant Stratification

<TopVariants />

<div class="grid grid-cols-3 gap-4 mt-6 text-xs">
<div class="p-3 rounded-lg" style="background:rgba(34,197,94,0.08); border:1px solid rgba(34,197,94,0.2)">

### 🎯 General Affinity
Globally highest-scoring binders across **all** targets. The "safe bet" for broad-spectrum inhibition.

</div>
<div class="p-3 rounded-lg" style="background:rgba(250,204,21,0.08); border:1px solid rgba(250,204,21,0.2)">

### 🔬 Target Specificity
High relative affinity for a **single** target. Enables selective therapeutic intervention.

</div>
<div class="p-3 rounded-lg" style="background:rgba(245,87,108,0.08); border:1px solid rgba(245,87,108,0.2)">

### ⚖️ Pairwise Selectivity
Differentiates between **closely related** target pairs (e.g., ADAM17 vs ADAM10).

</div>
</div>

---
layout: two-cols
---

# Stage 4: Experimental Validation

<div class="text-[13px] pr-8 overflow-hidden">

## Synthesis Pipeline
1. **Codon optimization** for yeast
2. **DNA synthesis** (Twist)
3. **Cloning** into display vectors
4. **Expression** and surface display

## Analysis Pipeline
- **`analyze_fcs.py`** — Single-target
- **`aggregate_analysis.py`** — Cross-target
- **`fcs_viewer.py`** — Interactive viewer

<div class="mt-4 p-3 rounded-lg border border-white border-opacity-10 bg-white bg-opacity-5 text-[11px] opacity-80">
  Click the targets on the right to toggle between different target binding results.
</div>

</div>

::right::

<div class="w-full h-full flex items-center pl-4 overflow-hidden">
  <FcsChart />
</div>

---
layout: center
---

# Experimental Validation — Gating Strategy

<div class="grid grid-cols-2 gap-10 mt-6 items-center max-w-6xl mx-auto overflow-hidden">
<div class="text-[13px]">

## Population Stratification
- **FITC (X-axis)**: Represents protein expression levels on the yeast surface.
- **APC (Y-axis)**: Represents binding affinity to the fluorescently-labeled target.

## Control Calibration
- **Neg Ctrl**: Unlabeled yeast or target-only to establish background gating.
- **Pos Ctrl**: Native TIMP3 to define the "Binding" quadrant (LR).
- **Variant**: Test population. Dual-positive signals indicate successful de novo binders.

<div class="mt-8 p-4 rounded-lg bg-primary bg-opacity-5 border border-primary border-opacity-20 italic text-[11px]">
  The quadrant gating allows us to normalize binding signal relative to expression, eliminating false negatives from low-expressing variants.
</div>

</div>
<div>
  <FcsScatter />
</div>
</div>

---
layout: center
class: text-center
---

# Conclusions & Future Directions

<div class="grid grid-cols-2 gap-10 text-left max-w-4xl mx-auto text-[13px] mt-10 overflow-hidden">
<div>

## What We Demonstrated
- ✅ End-to-end **generative pipeline** from backbone hallucination to validation
- ✅ Z-score consensus filtering **eliminates >90%** of candidates
- ✅ Loop confidence metrics (pLDDT, PAE) **predict functional binding**

</div>
<div>

## What's Next
- 🔄 **Feedback loop**: Flow cytometry results → retrain selection thresholds
- 🧬 **Combined loops**: Test multi-loop variants (AB + C + EF simultaneously)
- 📊 **Expanded targets**: MMP7, MMP10 catalytic domains
- 🤖 **ESM integration**: Complement generative designs with evolutionary classification

</div>
</div>

<div class="mt-12 opacity-50 text-[11px]">
  Ryan Gustafson · PhD Research · University of Nevada, Reno
</div>
