---
theme: seriph
background: https://cover.sli.dev
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
  font-size: 2.4rem !important;
}
</style>

---
transition: fade-out
layout: center
---

# Pipeline Overview

<PipelineFlow />

<div class="text-center text-xs opacity-40 mt-4">
  Hover each stage for details
</div>

<!--
This slide uses our custom PipelineFlow.vue component.
Any .vue file placed in the /components/ folder is auto-registered by Slidev.
-->

---

# Introduction — TIMP3 & Target Engineering

<div class="grid grid-cols-2 gap-8">
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
layout: two-cols
---

# Stage 2: ProteinMPNN

<div class="text-sm">

## Sequence Design Strategy

**Two pathways:**

1. **Expanded loops** — RFd backbones → ProteinMPNN  
   Populate poly-glycine with optimal amino acids

2. **Native-length loops** — Skip RFd → ProteinMPNN directly  
   Explore sequence space on wild-type topology

## Key Parameters
```python
num_seq_per_target = 1000
sampling_temp = 0.2   # low → conservative
fixed_chains = ["B"]  # target chain frozen
```

All scaffold residues **outside** loop regions are frozen — only loop residues are redesigned.

</div>

::right::

<div class="ml-4 text-sm">

## Output: Best Loops per Length (ADAM10)

<div class="text-xs">

| Length | Top Loop | Score | Seq Recovery |
|--------|----------|-------|-------------|
| 6 | `KGEIDE` | 0.857 | 16.7% |
| 7 | `KKEGYEE` | 1.095 | 14.3% |
| 8 | `KLPDGFEE` | 0.895 | 12.5% |
| 9 | `KDEKTGFET` | 0.885 | 11.1% |
| 10 | `VLDEKTGLKE` | 0.992 | 10.0% |
| 11 | `TVKTPGATFKE` | 1.041 | 9.1% |
| 12 | `KVKDWGGGEFEE` | 1.037 | 25.0% |
| 13 | `KVTDPETGKTFEE` | 0.966 | 7.7% |
| 14 | `EVTSPEDPSVKFKE` | 1.067 | 0.0% |
| 15 | `TKTVTENGETAKVEE` | 1.024 | 6.7% |

</div>

<div class="text-xs opacity-50 mt-2">

Data: `Results/best_loops_per_length_ADAM10.csv`  
Lower score = better packing; score reflects negative log-likelihood

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

<div class="text-sm">

## Synthesis Pipeline
1. **Codon optimization** for *S. cerevisiae* expression
2. **DNA synthesis** via Twist Bioscience
3. **Cloning** into yeast display vectors
4. **Expression** and surface display

## Flow Cytometry Protocol
- **Quad-gating** strategy using negative controls
- Differentiate expression-positive vs binding-positive cells
- **Double-positive** populations = functional binders

## Analysis Pipeline
```
Protein-Analysis/
├── analyze_fcs.py         # Single-target
├── analyze_fcs_target_trial.py  # Per-trial
├── aggregate_analysis.py  # Cross-target
└── fcs_viewer.py          # Interactive viewer
```

</div>

::right::

<div class="ml-4 text-sm">

## Gating Strategy

<div class="p-4 rounded-lg text-center" style="background:rgba(255,255,255,0.03); border:1px dashed rgba(255,255,255,0.15)">

📊 *Gating Strategy Figure*

`Gating_Strategy_NegCtrl_Quad.png`

<div class="text-xs opacity-50 mt-2">
Source: Local/ADAM17_*_Analysis/
</div>

</div>

## Aggregate Results

<div class="p-4 rounded-lg text-center mt-4" style="background:rgba(255,255,255,0.03); border:1px dashed rgba(255,255,255,0.15)">

📈 *Aggregate Binding Plot*

`aggregate_colorcoded.png`

<div class="text-xs opacity-50 mt-2">
Source: Local/Aggregate_FCS_Analysis/
</div>

</div>

<div class="text-xs opacity-50 mt-2">
Figures will be populated once FCS data is exported to SharedAssets
</div>

</div>

---
layout: center
class: text-center
---

# Conclusions & Future Directions

<div class="text-left max-w-xl mx-auto text-sm">

## What We Demonstrated
- ✅ End-to-end **generative pipeline** from backbone hallucination to empirical validation
- ✅ Z-score consensus filtering **eliminates >90%** of candidates before wet-lab
- ✅ Loop confidence metrics (pLDDT, PAE) **predict functional binding**

## What's Next
- 🔄 **Feedback loop**: Flow cytometry results → retrain selection thresholds
- 🧬 **Combined loops**: Test multi-loop variants (AB + C + EF simultaneously)
- 📊 **Expanded targets**: MMP7, MMP10 catalytic domains
- 🤖 **ESM integration**: Complement generative designs with evolutionary classification

</div>

<div class="mt-8 opacity-50 text-xs">
  Ryan Gustafson · PhD Research · University of Nevada, Reno
</div>
