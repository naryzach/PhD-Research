---
theme: seriph
background: '#0f172a'
title: "De Novo Binder Generation"
info: |
  ## De Novo Protein Prediction to Binding Evaluation
  Engineering MMP9-selective TIMP3 loop variants via generative deep learning.
  
  Author: Ryan Gustafson | University of Nevada, Reno
class: text-center
drawings:
  persist: false
transition: slide-left
mdc: true
---

<div class="layout: center">

# De Novo Binder Generation
### RFdiffusion → ProteinMPNN → AlphaFold3 → Yeast Display

<div class="mt-10 flex justify-center gap-4">
  <div class="px-4 py-2 bg-white/5 rounded border border-white/10">
    <div class="text-[10px] uppercase opacity-40 tracking-widest">Validated Results</div>
    <div class="text-blue-400 font-bold uppercase text-xs">3/4 Predicted Selective Confirmed</div>
  </div>
  <div class="px-4 py-2 bg-white/5 rounded border border-white/10">
    <div class="text-[10px] uppercase opacity-40 tracking-widest">Primary Targets</div>
    <div class="text-blue-400 font-bold uppercase text-xs">MMP9 vs MMP2</div>
  </div>
  <div class="px-4 py-2 bg-white/5 rounded border border-white/10">
    <div class="text-[10px] uppercase opacity-40 tracking-widest">PPV</div>
    <div class="text-emerald-400 font-bold uppercase text-xs">75%</div>
  </div>
</div>

</div>

<style>
h1 {
  background: linear-gradient(135deg, #58a6ff 0%, #00f2fe 100%);
  -webkit-background-clip: text;
  -webkit-text-fill-color: transparent;
  font-weight: 800;
  letter-spacing: -0.02em;
}
</style>

---
layout: default
transition: fade-out
---

# The Problem: Near-Identical Active Sites, Opposite Clinical Roles

<div class="grid grid-cols-2 gap-8 mt-5 items-start">
  <div class="space-y-3">
    <div class="p-4 bg-red-500/10 rounded-lg border border-red-500/20">
      <h4 class="text-red-400 font-bold text-xs uppercase tracking-widest mb-2">MMP9 — Pathogenic</h4>
      <ul class="text-[11px] space-y-1 opacity-80 list-none p-0">
        <li>→ Basement membrane degradation → tumor intravasation</li>
        <li>→ Releases ECM-sequestered VEGF → angiogenesis</li>
        <li>→ Overexpressed in >15 solid tumor types</li>
        <li>→ Correlates with metastasis and poor prognosis</li>
      </ul>
    </div>
    <div class="p-4 bg-blue-500/10 rounded-lg border border-blue-500/20">
      <h4 class="text-blue-400 font-bold text-xs uppercase tracking-widest mb-2">MMP2 — Homeostatic</h4>
      <ul class="text-[11px] space-y-1 opacity-80 list-none p-0">
        <li>→ Normal connective tissue maintenance</li>
        <li>→ Vascular remodeling in healthy tissue</li>
        <li>→ Inhibition → musculoskeletal toxicity</li>
        <li>→ 65% catalytic domain identity with MMP9</li>
      </ul>
    </div>
  </div>
  <div class="space-y-3">
    <div class="p-4 bg-white/5 rounded border border-white/10">
      <h4 class="text-[10px] font-bold text-blue-300 uppercase tracking-widest mb-2">Why TIMP3?</h4>
      <p class="text-[11px] leading-relaxed opacity-75">
        TIMP3 naturally inhibits metalloproteinases via three N-terminal reactive loops (AB, C, EF).
        Its β-barrel scaffold tolerates loop insertions of 6–15 aa while retaining fold stability —
        making it ideal for loop-level precision engineering.
      </p>
    </div>
    <div class="p-4 bg-white/5 rounded border border-white/10">
      <h4 class="text-[10px] font-bold text-emerald-300 uppercase tracking-widest mb-2">Our Approach</h4>
      <p class="text-[11px] leading-relaxed opacity-75">
        Generative AI de novo redesigns the AB and C loops to amplify MMP9 preference over MMP2.
        Computationally selected from >1,000 candidates, then experimentally validated via
        yeast surface display and flow cytometry.
      </p>
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# TIMP3 Scaffold — Engineerable Loop Architecture
Hover over each loop to see its position, native length, and expansion range.

<div class="mt-6">
  <LoopConfig />
</div>

<div class="grid grid-cols-3 gap-4 mt-6">
  <div class="p-3 bg-blue-500/10 rounded border border-blue-500/20 text-[10px]">
    <b class="text-blue-400 block mb-1">AB Loop (res 31–36)</b>
    <span class="opacity-70">Primary MMP active-site contact. AB 6 (13 aa insertion) achieved strongest selectivity signal, p = 0.0002.</span>
  </div>
  <div class="p-3 bg-violet-500/10 rounded border border-violet-500/20 text-[10px]">
    <b class="text-violet-400 block mb-1">C Loop (res 63–68)</b>
    <span class="opacity-70">Second primary contact loop. C 12 and C 15 both confirmed MMP9-selective. Up to 13 aa insertions tolerated.</span>
  </div>
  <div class="p-3 bg-emerald-500/10 rounded border border-emerald-500/20 text-[10px]">
    <b class="text-emerald-400 block mb-1">EF Loop (res 93–96)</b>
    <span class="opacity-70">Secondary contact loop. 4 residues native. EF variants generated but selection concentrated on AB and C in this cycle.</span>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Full Pipeline Architecture
Six stages from backbone hallucination to experimentally validated selectivity.

<div class="grid grid-cols-6 gap-2 mt-8">
  <div class="flex flex-col items-center p-3 bg-blue-500/10 rounded border border-blue-500/30 text-center">
    <div class="text-blue-400 font-black text-xl mb-1">①</div>
    <div class="text-[11px] font-bold text-white mb-2">RFdiffusion</div>
    <div class="text-[9px] opacity-50 leading-relaxed">Loop backbone hallucination<br/>AB / C / EF loops<br/>6–24 aa expansions<br/>A30 GPU, ~20 hr/run</div>
  </div>
  <div class="flex items-center justify-center text-blue-400 text-xl">→</div>
  <div class="flex flex-col items-center p-3 bg-violet-500/10 rounded border border-violet-500/30 text-center">
    <div class="text-violet-400 font-black text-xl mb-1">②</div>
    <div class="text-[11px] font-bold text-white mb-2">ProteinMPNN</div>
    <div class="text-[9px] opacity-50 leading-relaxed">1000 seqs/target<br/>T=0.2, fixed scaffold<br/>Top 10 per loop<br/>by log-likelihood</div>
  </div>
  <div class="flex items-center justify-center text-blue-400 text-xl">→</div>
  <div class="flex flex-col items-center p-3 bg-cyan-500/10 rounded border border-cyan-500/30 text-center">
    <div class="text-cyan-400 font-black text-xl mb-1">③</div>
    <div class="text-[11px] font-bold text-white mb-2">AlphaFold3</div>
    <div class="text-[9px] opacity-50 leading-relaxed">Batch co-folding<br/>5 targets<br/>ipTM / PAE / pLDDT<br/>loop & interface</div>
  </div>
  <div class="flex items-center justify-center text-blue-400 text-xl">→</div>
  <div class="flex flex-col items-center p-3 bg-emerald-500/10 rounded border border-emerald-500/30 text-center">
    <div class="text-emerald-400 font-black text-xl mb-1">④</div>
    <div class="text-[11px] font-bold text-white mb-2">T-Score Filter</div>
    <div class="text-[9px] opacity-50 leading-relaxed">T > 2.0, ≥2 metrics<br/>ipTM ≥ 0.82<br/>13 of >1000 pass<br/>>98% reduction</div>
  </div>
  <div class="flex items-center justify-center text-blue-400 text-xl">→</div>
  <div class="flex flex-col items-center p-3 bg-amber-500/10 rounded border border-amber-500/30 text-center">
    <div class="text-amber-400 font-black text-xl mb-1">⑤</div>
    <div class="text-[11px] font-bold text-white mb-2">DNA Synthesis</div>
    <div class="text-[9px] opacity-50 leading-relaxed">RE site screening<br/>BsrGI, BamHI, BsaI<br/>13/13 passed<br/>Twist, Dec 2025</div>
  </div>
  <div class="flex items-center justify-center text-blue-400 text-xl">→</div>
  <div class="flex flex-col items-center p-3 bg-red-500/10 rounded border border-red-500/30 text-center">
    <div class="text-red-400 font-black text-xl mb-1">⑥</div>
    <div class="text-[11px] font-bold text-white mb-2">Yeast Display</div>
    <div class="text-[9px] opacity-50 leading-relaxed">FITC = expression<br/>APC = binding<br/>Bind Med (Expr+)<br/>ANOVA + Tukey-HSD</div>
  </div>
</div>

<div class="mt-8 p-3 bg-white/5 rounded border border-white/10 text-[10px] opacity-75">
  <span class="text-blue-400 font-bold">Scaffold constraint:</span> TIMP3 residues outside the three engineered loops are held fixed throughout stages ①–③. Only loop positions are hallucinated (①) or designed (②), ensuring the validated scaffold fold is preserved.
</div>

---
layout: default
transition: fade-out
---

# Pipeline Phase Detail
Computational filtering (stages ①–④) vs experimental validation (stages ⑤–⑥).

<div class="grid grid-cols-2 gap-6 mt-6">
  <div class="space-y-3">
    <h3 class="text-blue-400 font-bold uppercase text-xs tracking-widest border-b border-blue-500/20 pb-1">Computational Phase</h3>
    <div class="p-3 bg-blue-500/5 rounded border border-blue-500/15 space-y-2 text-[10px] leading-relaxed opacity-80">
      <p><b class="text-blue-300">HADDOCK templates:</b> Each target provides a docked starting structure. Scaffold residues are frozen in all subsequent stages — only loop backbone (①) and sequence (②) vary.</p>
      <p><b class="text-blue-300">T-score normalization:</b> Each variant is scored relative to its own mean across all targets — not the population. A mediocre global binder with one uniquely-high target score earns a high T. This explicitly selects preferential binders, not globally strong ones.</p>
      <p><b class="text-blue-300">Multi-metric consensus:</b> T > 2.0 required on ≥2 independent AF3 metrics (ipTM, loop pLDDT, interface PAE). Prevents noise-driven false positives from a single metric.</p>
      <p><b class="text-blue-300">Reduction:</b> >1,000 designed sequences → 13 synthesis candidates (>98% reduction).</p>
    </div>
  </div>
  <div class="space-y-3">
    <h3 class="text-cyan-400 font-bold uppercase text-xs tracking-widest border-b border-cyan-500/20 pb-1">Experimental Phase</h3>
    <div class="p-3 bg-cyan-500/5 rounded border border-cyan-500/15 space-y-2 text-[10px] leading-relaxed opacity-80">
      <p><b class="text-cyan-300">Yeast surface display:</b> Golden Gate cloning into pCT302; EBY100 yeast transformation. Surface expression detected by FITC-conjugated anti-Myc antibody. His-tagged target detected by APC-conjugated anti-His.</p>
      <p><b class="text-cyan-300">Gating:</b> Pentagon FSC/SSC singlet gate learned from NC via KDE. Quadrant thresholds at 99.5th percentile NC per channel — &lt;0.5% NC events above threshold in either channel.</p>
      <p><b class="text-cyan-300">Primary metric:</b> Bind Med (Expr+) — raw APC-A median for FITC+ cells. Cross-target comparable because MMP9 and MMP2 values share the same fluorescence scale.</p>
      <p><b class="text-cyan-300">Statistics:</b> One-way ANOVA per construct across targets; Tukey-HSD post-hoc for pairwise target comparisons.</p>
    </div>
  </div>
</div>

---
layout: section
transition: slide-up
---

# Phase 1
### Computational Candidate Generation and Selection

---
layout: default
transition: fade-out
---

# AlphaFold Confidence Metrics — What We Extract
Five independent structural confidence signals per co-folded complex.

<div class="grid grid-cols-5 gap-3 mt-8">
  <div class="p-4 bg-blue-500/10 rounded border border-blue-500/20 space-y-2 text-center">
    <div class="text-blue-400 font-black text-base">pTM</div>
    <div class="text-[10px] opacity-60 leading-relaxed">Global complex TM-score. Hard threshold ≥ 0.80. Confirms overall fold integrity of the TIMP3–target complex.</div>
  </div>
  <div class="p-4 bg-cyan-500/10 rounded border border-cyan-500/20 space-y-2 text-center">
    <div class="text-cyan-400 font-black text-base">ipTM</div>
    <div class="text-[10px] opacity-60 leading-relaxed">Interface TM-score. Hard threshold ≥ 0.82. Most sensitive to loop-mediated binding geometry. Primary selection gatekeeper.</div>
  </div>
  <div class="p-4 bg-violet-500/10 rounded border border-violet-500/20 space-y-2 text-center">
    <div class="text-violet-400 font-black text-base">Mean pLDDT</div>
    <div class="text-[10px] opacity-60 leading-relaxed">Per-residue confidence, all residues. ≥70 = confident global fold. Dominated by scaffold; less discriminating for loop design.</div>
  </div>
  <div class="p-4 bg-emerald-500/10 rounded border border-emerald-500/20 space-y-2 text-center">
    <div class="text-emerald-400 font-black text-base">Loop pLDDT</div>
    <div class="text-[10px] opacity-60 leading-relaxed">pLDDT restricted to engineered loop residues only. More discriminating than mean pLDDT — directly reflects predicted loop structure confidence.</div>
  </div>
  <div class="p-4 bg-amber-500/10 rounded border border-amber-500/20 space-y-2 text-center">
    <div class="text-amber-400 font-black text-base">PAE (iface)</div>
    <div class="text-[10px] opacity-60 leading-relaxed">Expected position error of loop residues relative to the target. Lower = more precise predicted binding geometry. Units: Å.</div>
  </div>
</div>

<div class="mt-8 grid grid-cols-2 gap-4">
  <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px] opacity-70 leading-relaxed">
    <b class="text-blue-400">Why five metrics?</b> Each captures a different aspect of the predicted interface: global fold (pTM), interface quality (ipTM), whole-structure confidence (mean pLDDT), loop-specific confidence (loop pLDDT), geometric precision (PAE). A binder that scores high on all five is more likely to be a true positive than one that excels on only one.
  </div>
  <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px] opacity-70 leading-relaxed">
    <b class="text-cyan-400">Selection logic:</b> Hard thresholds (ipTM ≥ 0.82, pTM ≥ 0.80) are applied first as gates. Remaining candidates are ranked by T-score — a self-normalized metric that compares each variant against its own cross-target distribution, not the population. See next slide for T-score details.
  </div>
</div>

---
layout: default
transition: fade-out
---

# T-Score: Selecting for Preference, Not Just Strength
Self-normalized — each variant is its own control across all targets.

<div class="mt-4">
  <TScoreFormula />
</div>

<div class="grid grid-cols-3 gap-3 mt-4">
  <div class="p-2 bg-emerald-500/10 rounded border border-emerald-500/20 text-center">
    <div class="text-emerald-400 font-black text-lg">T > 2.0</div>
    <div class="text-[8px] uppercase tracking-widest text-emerald-400 mt-1">97.5th %ile</div>
    <div class="text-[8px] opacity-50">Required on ≥2 metrics</div>
  </div>
  <div class="p-2 bg-blue-500/10 rounded border border-blue-500/20 text-center">
    <div class="text-blue-400 font-black text-lg">ipTM ≥ 0.82</div>
    <div class="text-[8px] uppercase tracking-widest text-blue-400 mt-1">Interface floor</div>
    <div class="text-[8px] opacity-50">Hard gate applied first</div>
  </div>
  <div class="p-2 bg-violet-500/10 rounded border border-violet-500/20 text-center">
    <div class="text-violet-400 font-black text-lg">13 / >1000</div>
    <div class="text-[8px] uppercase tracking-widest text-violet-400 mt-1">Pass filter</div>
    <div class="text-[8px] opacity-50">>98% candidate reduction</div>
  </div>
</div>

<div class="grid grid-cols-2 gap-4 mt-3">
  <div class="p-3 bg-white/4 rounded border border-white/10 text-[10px] leading-relaxed opacity-80">
    <b class="text-blue-400">T-score vs Z-score:</b> Z-score identifies variants that bind better <em>than the population</em> — a globally excellent binder earns a high Z regardless of selectivity.
    T-score identifies variants that bind <em>this target better than their own cross-target average</em>.
    A mediocre global binder with one uniquely high score earns a high T — exactly the selectivity-optimizing behavior we want.
  </div>
  <div class="p-3 bg-emerald-500/5 rounded border border-emerald-500/20 text-[10px] leading-relaxed opacity-80">
    <b class="text-emerald-400">Multi-metric consensus:</b> T > 2.0 required on ≥2 independent metrics (ipTM, loop pLDDT, interface PAE). Agreement across orthogonal signals prevents noise-driven false positives from any single metric.
  </div>
</div>

---
layout: default
transition: fade-out
---

# T-Score Results — Statistical Significance Matrix
Filter by loop type; scroll for all variants. Green = T > 2.0 (significant win on that target).

<div class="mt-4" style="height: calc(100vh - 140px);">
  <TScoreTable />
</div>

---
layout: default
transition: fade-out
---

# Final Ordered Library — 13 Constructs
All constructs passed RE site screening (BsrGI, BamHI, BsaI).

<div class="mt-6">
  <TopVariants />
</div>

---
layout: default
transition: fade-out
---

# Final Library — Loop Group Breakdown

<div class="grid grid-cols-3 gap-5 mt-8">
  <div class="p-4 bg-blue-500/10 rounded border border-blue-500/20 space-y-2">
    <h4 class="text-blue-300 font-bold text-[11px] uppercase tracking-widest">AB-Loop Variants (7)</h4>
    <p class="text-[10px] leading-relaxed opacity-70">
      Native residues 31–36; extensions 6–15 aa. Targeted at MMP9 active site
      via HADDOCK geometry. AB 6 showed the strongest statistical signal (p = 0.0002).
      AB 1 was predicted selective — directionally correct but ANOVA non-significant due to high trial-to-trial variance.
    </p>
  </div>
  <div class="p-4 bg-violet-500/10 rounded border border-violet-500/20 space-y-2">
    <h4 class="text-violet-300 font-bold text-[11px] uppercase tracking-widest">C-Loop Variants (5)</h4>
    <p class="text-[10px] leading-relaxed opacity-70">
      Native residues 63–68. C 12 (loop: ASGPITVNGETIW) and C 15 both confirmed
      MMP9-selective at p &lt; 0.02. Up to 13 aa insertions tolerated without
      loss of yeast surface display function.
    </p>
  </div>
  <div class="p-4 bg-amber-500/10 rounded border border-amber-500/20 space-y-2">
    <h4 class="text-amber-300 font-bold text-[11px] uppercase tracking-widest">ADAM17-Targeted (4)</h4>
    <p class="text-[10px] leading-relaxed opacity-70">
      AB 4, AB 7, C 11, C 14 designed for ADAM17. ADAM10 comparison pending
      positive control optimization. Included in within-target analysis only.
      Extend the pipeline to inflammatory signaling vs Alzheimer's disease target pair.
    </p>
  </div>
</div>

<div class="mt-6 p-3 bg-white/5 rounded border border-white/10 text-[10px] opacity-70">
  <b class="text-blue-400">Library design rationale:</b> MMP9-targeted constructs span both primary contact loops (AB and C) to probe whether selectivity is loop-position–specific or generalizable across the TIMP3 interface. The ADAM17 group tests pipeline generalizability beyond the MMP scaffold in the same experimental cycle.
</div>

---
layout: section
transition: slide-up
---

# Phase 2
### Experimental Validation — Flow Cytometry Metrics

---
layout: default
transition: fade-out
---

# Gating Strategy
Pentagon scatter gate + 99.5th percentile NC quadrant thresholds.

<div class="grid grid-cols-2 gap-8 mt-6 items-start">
  <div class="space-y-3">
    <div class="p-3 bg-blue-500/5 rounded border border-blue-500/20">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-2">Step 1 — Pentagon Scatter Gate</h4>
      <p class="text-[10px] leading-relaxed opacity-70">
        KDE on FSC-A vs SSC-A from negative control. Binary search for 90% core-event inclusion.
        Convex hull simplified to 5 vertices (Visvalingam-Whyatt). Applied identically to all
        samples in the trial — ensures no sample-specific selection bias.
      </p>
    </div>
    <div class="p-3 bg-blue-500/5 rounded border border-blue-500/20">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-2">Step 2 — Dynamic Quadrant Gate</h4>
      <p class="text-[10px] leading-relaxed opacity-70">
        Expression threshold (θ<sub>expr</sub>): 99.5th percentile of FITC-A in NC.<br/>
        Binding threshold (θ<sub>bind</sub>): 99.5th percentile of APC-A in NC.<br/>
        &lt;0.5% NC events above threshold per channel. Q2 (upper-right) = double positive (DP) — expressing <em>and</em> binding.
      </p>
    </div>
  </div>
  <div>
    <FcsViewer />
  </div>
</div>

---
layout: default
transition: fade-out
---

# Flow Cytometry Metrics — Hover for Detail
Six quantitative metrics extracted per sample. Hover any card for full definition, scope, and use.

<div class="mt-4">
  <MetricExplorer />
</div>

<div class="mt-4 grid grid-cols-3 gap-2 text-[9px] opacity-50">
  <div><span class="text-cyan-400 font-bold">Cross-target:</span> Bind Med (Expr+), Binding Efficiency</div>
  <div><span class="text-violet-400 font-bold">Within-target:</span> Norm Bind Med, Norm Median Ratio, Norm IWB</div>
  <div><span class="text-slate-400 font-bold">QC only:</span> Stain Index</div>
</div>

---
layout: default
transition: fade-out
---

# Primary Metric: Positive Median Ratio (APC/FITC)
APC/FITC median ratio in expressing cells — cross-target comparable, expression-corrected.

<div class="grid grid-cols-2 gap-8 mt-4 items-start">
  <div class="space-y-4">
    <div class="p-4 bg-white/5 rounded border border-white/10">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-3">Definition</h4>
      <div class="font-mono text-[11px] bg-black/30 rounded p-3 leading-relaxed">
        <span class="text-cyan-400">Pos Med Ratio</span> =<br/>
        &nbsp;&nbsp;median(APC-A) / median(FITC-A)&nbsp;&nbsp;[FITC+ gate]
        <br/><br/>
        <span class="text-emerald-400">Norm Median Ratio</span> =<br/>
        &nbsp;&nbsp;Pos Med Ratio(sample) / Pos Med Ratio(TIMP3-WT)
      </div>
      <p class="text-[10px] leading-relaxed opacity-60 mt-3">
        Restricted to FITC+ cells. Dividing APC by FITC corrects for surface-display
        variation between individual cells. <b>Cross-target comparable:</b> MMP9 and MMP2 share
        the same FITC detection antibody so the denominator is equivalent across targets.
        Error bars = 95% CI (t-distribution). Failed trials excluded.
      </p>
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px]">
      <b class="text-amber-400">vs Binding Efficiency:</b>
      Efficiency = <em>fraction</em> of expressors that cross the APC threshold.
      Pos Med Ratio = <em>continuous intensity</em> of binding per unit expression.
      Pos Med Ratio uses all FITC+ cells; Efficiency uses a binary gate. Both are cross-target comparable.
    </div>
  </div>
  <div class="space-y-3">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/fig_bind_med_expr_pos.png"
         alt="Positive Median Ratio"
         class="w-full rounded-lg border border-white/10"
         style="max-height: 300px; object-fit: contain;">
    <p class="text-[9px] opacity-40 italic text-center">
      Pos Med Ratio per construct — MMP9 (green) vs MMP2 (red), APC/FITC.
      Dashed lines = TIMP3-WT ratio per target. Error bars = 95% CI. Significance by Tukey-HSD.
    </p>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Within-Target Ranking: Norm Bind Med (Expr+)
How do variants compare to TIMP3-WT on the <em>same</em> target?

<div class="grid grid-cols-2 gap-6 mt-4 items-start">
  <div>
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/fig_bind_med_heatmap.png"
         alt="Norm Bind Med (Expr+) heatmap"
         class="w-full rounded-lg border border-white/10"
         style="max-height: 360px; object-fit: contain;">
    <p class="text-[9px] opacity-40 italic text-center mt-1">
      Each target normalized independently to TIMP3-WT (= 1.0).
      Do not compare MMP9 values to MMP2 values across rows.
    </p>
  </div>
  <div class="space-y-2">
    <div class="p-2 bg-white/5 rounded border border-white/10 text-[9.5px] leading-snug opacity-80">
      <b class="text-blue-400">Why normalize within-target:</b> Raw MMP9 and MMP2 baselines differ due to protein concentration, staining efficiency, and TIMP3 affinity per target. Dividing by TIMP3-WT for each target independently isolates construct-level effects from between-target scale differences.
    </div>
    <div class="p-2 bg-emerald-500/10 rounded border border-emerald-500/20 text-[9.5px]">
      <b class="text-emerald-400 block mb-1">Key constructs vs TIMP3-WT (Norm Bind Med):</b>
      <table class="w-full font-mono border-collapse text-[9px]">
        <thead>
          <tr class="border-b border-white/20">
            <th class="text-left py-0.5 opacity-50">Construct</th>
            <th class="text-right py-0.5 opacity-50">MMP9</th>
            <th class="text-right py-0.5 opacity-50">MMP2</th>
          </tr>
        </thead>
        <tbody>
          <tr class="border-b border-white/5"><td class="py-0.5 text-emerald-400 font-bold">C 12</td><td class="text-right text-emerald-400">0.99</td><td class="text-right opacity-60">0.88</td></tr>
          <tr class="border-b border-white/5"><td class="py-0.5 text-emerald-400 font-bold">C 15</td><td class="text-right text-emerald-400">0.92</td><td class="text-right opacity-60">0.89</td></tr>
          <tr class="border-b border-white/5"><td class="py-0.5 text-emerald-400 font-bold">AB 6</td><td class="text-right text-emerald-400">0.85</td><td class="text-right opacity-60">0.78</td></tr>
          <tr class="border-b border-white/5 opacity-50"><td class="py-0.5">C 13</td><td class="text-right">0.73</td><td class="text-right">1.11</td></tr>
          <tr><td class="py-0.5 text-blue-400">TIMP 3</td><td class="text-right text-blue-400">1.00</td><td class="text-right text-blue-400">1.00</td></tr>
        </tbody>
      </table>
    </div>
    <div class="p-2 bg-amber-500/10 rounded border border-amber-500/20 text-[9.5px] opacity-80">
      <b class="text-amber-400">Interpretation:</b> Values near 1.0 = binding similar to TIMP3-WT for that target. Selectivity signal comes from the raw cross-target comparison — this normalized view confirms variants are not simply super-binders on every target.
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Secondary Metric: Binding Efficiency (DP/FITC+)
Cross-target comparable fraction metric — primary selectivity proof.

<div class="grid grid-cols-2 gap-8 mt-4 items-start">
  <div class="space-y-4">
    <div class="p-4 bg-white/5 rounded border border-white/10">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-3">Definition</h4>
      <div class="font-mono text-[11px] bg-black/30 rounded p-3 leading-relaxed">
        <span class="text-cyan-400">Binding Efficiency</span> =<br/>
        &nbsp;&nbsp;%DP / %Expr+ = N(Q2) / N(Q1+Q2)
      </div>
      <p class="text-[10px] leading-relaxed opacity-60 mt-3">
        Of the cells successfully displaying the variant (FITC+), what fraction
        also show detectable target binding (APC+)? Denominator is expression-normalized,
        making it directly comparable across targets and trials.
      </p>
    </div>
    <div class="p-3 bg-emerald-500/10 rounded border border-emerald-500/20 text-[10px]">
      <b class="text-emerald-400">Cross-target comparable</b> because the FITC
      antibody (anti-Myc) is target-independent. The denominator is the same
      cell population regardless of which MMP is added. This makes MMP9 vs MMP2
      Binding Efficiency values directly comparable.
    </div>
  </div>
  <div>
    <div class="p-4 bg-white/5 rounded border border-white/10 space-y-3">
      <h3 class="text-xs font-bold uppercase tracking-[0.2em] text-blue-400">Selectivity Results</h3>
      <table class="w-full text-[10px] font-mono border-collapse">
        <thead>
          <tr class="border-b border-white/20">
            <th class="text-left py-1 opacity-50">Construct</th>
            <th class="text-right py-1 opacity-50">MMP2</th>
            <th class="text-right py-1 opacity-50 text-blue-400">MMP9</th>
            <th class="text-right py-1 opacity-50">Ratio</th>
            <th class="text-right py-1 opacity-50">Result</th>
          </tr>
        </thead>
        <tbody>
          <tr class="border-b border-white/5">
            <td class="py-1 text-emerald-400 font-bold">C 12</td>
            <td class="text-right py-1">32.1%</td>
            <td class="text-right py-1 text-emerald-400 font-bold">91.3%</td>
            <td class="text-right py-1 text-emerald-400">2.8×</td>
            <td class="text-right py-1 text-emerald-400">✓ **</td>
          </tr>
          <tr class="border-b border-white/5">
            <td class="py-1 text-emerald-400 font-bold">C 15</td>
            <td class="text-right py-1">33.8%</td>
            <td class="text-right py-1 text-emerald-400 font-bold">88.6%</td>
            <td class="text-right py-1 text-emerald-400">2.6×</td>
            <td class="text-right py-1 text-emerald-400">✓ *</td>
          </tr>
          <tr class="border-b border-white/5">
            <td class="py-1 text-emerald-400 font-bold">AB 6</td>
            <td class="text-right py-1">23.8%</td>
            <td class="text-right py-1 text-emerald-400 font-bold">89.9%</td>
            <td class="text-right py-1 text-emerald-400">3.8×</td>
            <td class="text-right py-1 text-emerald-400">✓ ***</td>
          </tr>
          <tr class="border-b border-white/10 opacity-50">
            <td class="py-1">C 13</td>
            <td class="text-right py-1">16.0%</td>
            <td class="text-right py-1">45.4%</td>
            <td class="text-right py-1">—</td>
            <td class="text-right py-1">n.s. ✓</td>
          </tr>
          <tr class="border-b border-white/10 opacity-50">
            <td class="py-1">AB 5</td>
            <td class="text-right py-1">11.8%</td>
            <td class="text-right py-1">37.6%</td>
            <td class="text-right py-1">—</td>
            <td class="text-right py-1">n.s. ✓</td>
          </tr>
          <tr>
            <td class="py-1 text-blue-400">TIMP 3</td>
            <td class="text-right py-1">21.9%</td>
            <td class="text-right py-1">54.5%</td>
            <td class="text-right py-1">2.5×</td>
            <td class="text-right py-1 text-blue-400">ref *</td>
          </tr>
        </tbody>
      </table>
      <div class="text-[8px] opacity-40 italic border-t border-white/10 pt-2">
        * p&lt;0.05, ** p&lt;0.01, *** p&lt;0.001 (Tukey-HSD, MMP9 vs MMP2 pair)
      </div>
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Positive Median Ratio — Cross-Target Comparison
APC/FITC ratio for FITC+ cells. Click group buttons to switch target pairs.

<div class="mt-4 flex justify-center">
  <FcsChart />
</div>

<div class="mt-4 grid grid-cols-4 gap-3 text-[10px]">
  <div class="p-2 bg-white/5 rounded border border-white/10 opacity-80">
    <b class="text-cyan-400">Bind Med (Expr+)</b><br/>
    Raw APC-A intensity in expressors. Cross-target comparable. Primary poster metric.
  </div>
  <div class="p-2 bg-white/5 rounded border border-white/10 opacity-80">
    <b class="text-emerald-400">Binding Efficiency</b><br/>
    Fraction of expressors that bind. Cross-target. Primary selectivity proof.
  </div>
  <div class="p-2 bg-white/5 rounded border border-white/10 opacity-80">
    <b class="text-violet-400">Norm Median Ratio</b><br/>
    APC/FITC ratio, WT-normed. Corrects for expression-level variation.
  </div>
  <div class="p-2 bg-white/5 rounded border border-white/10 opacity-80">
    <b class="text-amber-400">Norm IWB Index</b><br/>
    Per-cell binding quality in DP population. Sensitive to threshold-crossing cells only.
  </div>
</div>

---
layout: default
transition: fade-out
---

# Key Findings: 3 of 4 Predictions Confirmed
75% PPV — strong given the vast sequence landscape. Every control non-significant.

<div class="grid grid-cols-2 gap-8 mt-6">
  <div class="space-y-4">
    <h3 class="text-emerald-400 font-bold uppercase text-xs tracking-widest border-b border-emerald-500/20 pb-1">Positive Outcomes (3/4 Predicted Selective)</h3>
    <ul class="text-[11px] space-y-4 opacity-85 list-none p-0">
      <li class="flex gap-3 items-start">
        <div class="w-7 h-7 rounded-full bg-emerald-500/20 flex items-center justify-center text-[10px] font-black text-emerald-400 shrink-0">01</div>
        <span><b class="text-emerald-300">3 of 4 Selective Confirmed:</b> C 12 (p=0.012), C 15 (p=0.018), AB 6 (p=0.0002) confirmed by ANOVA + Tukey-HSD. Selectivity ratios 2.6–3.8× vs TIMP3-WT baseline (2.5×). AB 1 showed directional preference but ANOVA non-significant due to high trial-to-trial variance (σ = 43%).</span>
      </li>
      <li class="flex gap-3 items-start">
        <div class="w-7 h-7 rounded-full bg-blue-500/20 flex items-center justify-center text-[10px] font-black text-blue-400 shrink-0">02</div>
        <span><b class="text-blue-300">Zero False Positives:</b> All 4 non-selective control constructs (C 13, AB 5, AB 2, AB 3) showed non-significant ANOVA (p > 0.05) — validating the T-score filter's specificity, not just sensitivity.</span>
      </li>
      <li class="flex gap-3 items-start">
        <div class="w-7 h-7 rounded-full bg-violet-500/20 flex items-center justify-center text-[10px] font-black text-violet-400 shrink-0">03</div>
        <span><b class="text-violet-300">Scaffold Plasticity Confirmed:</b> 6–15 aa loop insertions tolerated structurally and biologically. WT TIMP3 itself shows MMP9 preference (p=0.024) — engineered variants amplify a native preference.</span>
      </li>
    </ul>
  </div>

  <div class="space-y-4">
    <h3 class="text-blue-400 font-bold uppercase text-xs tracking-widest border-b border-blue-500/20 pb-1">Pipeline Validation Logic</h3>
    <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px] leading-relaxed opacity-75 space-y-2">
      <p><b class="text-blue-300">Context for 75% PPV:</b> Binders are not 50/50 by chance — the background rate of a randomly chosen sequence being both functional on yeast display AND MMP9-selective is far below 75%. Three statistically confirmed predictions from four targeted designs is a strong result against that baseline.</p>
      <p><b class="text-blue-300">Bidirectional concordance:</b> 0 false positives in the non-selective set. A filter that only gets selective constructs right might be too lenient — confirmed true negatives validate the filter's specificity, not just sensitivity.</p>
      <p><b class="text-blue-300">Platform implication:</b> 75% PPV at N=4 with 0 FP at N=4 justifies SPR kinetic characterization and next-cycle loop redesign without experimental redesign of the gating or statistical framework.</p>
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Next Steps and Translational Vision

<div class="grid grid-cols-2 gap-8 mt-6">
  <div class="space-y-3">
    <h3 class="text-cyan-400 font-bold uppercase text-xs tracking-widest border-b border-cyan-500/20 pb-1">Near-Term Experiments</h3>
    <ul class="text-[11px] space-y-3 opacity-80 list-none p-0">
      <li class="flex items-start gap-3">
        <div class="w-1.5 h-1.5 rounded-full bg-cyan-400 mt-1.5 shrink-0"></div>
        <span><b>SPR Kinetics:</b> K<sub>D</sub>, k<sub>on</sub>, k<sub>off</sub> for C 12, C 15, AB 6 vs purified MMP9 and MMP2. Absolute affinities enable direct comparison to marimastat/prinomastat clinical benchmarks.</span>
      </li>
      <li class="flex items-start gap-3">
        <div class="w-1.5 h-1.5 rounded-full bg-cyan-400 mt-1.5 shrink-0"></div>
        <span><b>Reverse Selectivity (MMP2 > MMP9):</b> Demonstrates controllable pipeline steering. If the same pipeline produces MMP2-preferring variants on demand, selectivity is engineerable in either direction — not an artifact of the TIMP3 scaffold.</span>
      </li>
      <li class="flex items-start gap-3">
        <div class="w-1.5 h-1.5 rounded-full bg-cyan-400 mt-1.5 shrink-0"></div>
        <span><b>ADAM10 Protocol Optimization:</b> Titrate antibody concentration and induction time to achieve sufficient positive control signal; enables ADAM17/ADAM10 selectivity analysis with identical ANOVA framework.</span>
      </li>
      <li class="flex items-start gap-3">
        <div class="w-1.5 h-1.5 rounded-full bg-cyan-400 mt-1.5 shrink-0"></div>
        <span><b>ESM-2 Saturation Mutagenesis:</b> Per-residue fitness landscape for all positions in C 12, C 15, AB 6 loops — identify critical selectivity determinants and tolerant positions for second-round design.</span>
      </li>
    </ul>
  </div>

  <div class="space-y-3">
    <h3 class="text-violet-400 font-bold uppercase text-xs tracking-widest border-b border-violet-500/20 pb-1">Translational Context</h3>
    <div class="p-4 bg-violet-500/5 rounded border border-violet-500/20 text-[10px] leading-relaxed opacity-80 space-y-3">
      <p><b class="text-violet-300">MMP9 as oncology target:</b> Overexpressed in >15 solid tumor types; correlates with metastasis and poor prognosis. Prior pan-MMP inhibitors failed Phase III due to MMP2 off-target toxicity — not because MMP9 inhibition is wrong, but because selectivity was insufficient. An MMP9-specific TIMP3 variant addresses the actual failure mode.</p>
      <p><b class="text-violet-300">Target-agnostic framework:</b> Any protein–protein interaction with available structural data is addressable by this pipeline. MMP9/MMP2 is the validation case; the computational workflow is identical for kinase isoforms, cytokine receptors, or viral surface proteins.</p>
      <p><b class="text-violet-300">Precision CAR-T analogy:</b> Patient tumor genome → neoantigen identification → de novo binder generation → CAR construct with hyperspecific, genome-informed recognition. Combines precision genomics with AI-driven protein design to minimize on-target/off-tumor toxicity.</p>
    </div>
  </div>
</div>

---
layout: section
transition: slide-up
---

# Cross-Reactivity Panel
### TwistBio Constructs vs MMP1, MMP7, MMP8, MMP10

---
layout: default
transition: fade-out
---

# Cross-Reactivity Binding Profile
Flow Cytometry assessment of TwistBio constructs against a multi-MMP panel.

<div class="grid grid-cols-2 gap-6 mt-6 items-start">
  <div class="space-y-3">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/20260529_CrossReactivity/Aggregate_FoldChange_Binding.png"
         alt="Binding Fold Change"
         class="w-full rounded-lg border border-white/10"
         style="max-height: 300px; object-fit: contain;">
  </div>
  <div class="space-y-4">
    <div class="p-4 bg-red-500/10 rounded border border-red-500/20">
      <h4 class="text-red-400 font-bold text-[10px] uppercase tracking-widest mb-2">Observation: Weak Signal Across Panel</h4>
      <p class="text-[10px] leading-relaxed opacity-70">
        All constructs (AB 3, 4, 6 and C 12) showed weak/minimal binding across all tested targets (MMP1, 7, 8, 10). The average binding fold change was ~1.18x over background.
      </p>
    </div>
    <div class="p-4 bg-amber-500/10 rounded border border-amber-500/20">
      <h4 class="text-amber-400 font-bold text-[10px] uppercase tracking-widest mb-2">Control Validation Issue</h4>
      <p class="text-[10px] leading-relaxed opacity-70">
        The positive control (TIMP 3) also exhibited extremely low binding (1.17x fold change). This indicates a potential issue with target biotinylation/labeling or the secondary staining protocol, rather than an inherent failure of the novel binders. Assays should be repeated with fresh/validated target batches.
      </p>
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Yeast Display Target Trial — ADAM10 & MMP9 (June 10\textsuperscript{th}, 2026)
Flow Cytometry assessment of TwistBio constructs against commercial and in-house targets.

<div class="grid grid-cols-2 gap-6 mt-6 items-start">
  <div class="space-y-3">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/20260610_YeastDisplayTargetTrial/ADAM10_F6_R_Ridgeline_Binding.png"
         alt="ADAM10 Fraction 6 Ridgeline Binding"
         class="w-full rounded-lg border border-white/10"
         style="max-height: 250px; object-fit: contain;">
    <p class="text-[9px] opacity-40 italic text-center">
      ADAM10 Fraction 6 (Ryan's run _R, FITC-A Binding, APC-A Expression).
    </p>
  </div>
  <div class="space-y-4">
    <div class="p-3 bg-emerald-500/10 rounded border border-emerald-500/20 text-[10px]">
      <h4 class="text-emerald-400 font-bold text-[10px] uppercase tracking-widest mb-1">MMP9 Positives Validated</h4>
      <p class="leading-relaxed opacity-75">
        <b>Sino MMP9 binding</b> was strong: positive control <b>TIMP3</b> achieved 38.2% Double+ (491 MFI), and Twist variants <b>AB 7</b> and <b>AB 3</b> achieved 33.5% and 32.9% Double+ respectively, validating the display system.
      </p>
    </div>
    <div class="p-3 bg-blue-500/10 rounded border border-blue-500/20 text-[10px]">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-1">ADAM10 Binding Confirmed</h4>
      <p class="leading-relaxed opacity-75">
        Re-analysis with the correct channel configuration (FITC-A binding) revealed robust binding for <b>ADAM10</b> F6 and F9. <b>AB 7</b> achieved 13.29% Double+ (1.56x MFI vs TIMP 3 positive control), and <b>AB 1</b> reached 11.71% Double+.
      </p>
    </div>
    <div class="p-3 bg-violet-500/10 rounded border border-violet-500/20 text-[10px]">
      <h4 class="text-violet-400 font-bold text-[10px] uppercase tracking-widest mb-1">Commercial AbCam Validation</h4>
      <p class="leading-relaxed opacity-75">
        Commercial <b>AbCam ADAM10</b> confirmed target activity: <b>TIMP 3</b> positive control achieved 12.65% Double+, and Twist variants <b>AB 5</b> and <b>AB 3</b> showed strong binding (MFI 1708 and 1665.5).
      </p>
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Target Trial Analysis & Channel Correction
Resolving analysis channel mismatch and explaining detection reagent failures.

<div class="grid grid-cols-2 gap-6 mt-6 items-start">
  <div class="space-y-3">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/20260610_YeastDisplayTargetTrial/ADAM10_F6_R_FoldChange_Binding.png"
         alt="ADAM10 Fold Change Binding"
         class="w-full rounded-lg border border-white/10"
         style="max-height: 220px; object-fit: contain;">
    <div class="p-3 bg-blue-500/10 rounded border border-blue-500/20 text-[9px]">
      <h4 class="text-blue-400 font-bold text-[9px] uppercase tracking-widest mb-1">Channel Mismatch Resolved</h4>
      <p class="leading-relaxed opacity-75">
        The default pipeline assumed a PE-A binding channel (Gavin's protocol). Re-running the pipeline with the correct <b>FITC-A</b> channel configuration for Ryan's run (<code>_R</code>) corrected the false-negative results.
      </p>
    </div>
  </div>
  <div class="space-y-4">
    <div class="p-4 bg-white/5 rounded border border-white/10 text-[10px]">
      <h4 class="text-blue-300 font-bold text-[10px] uppercase tracking-widest mb-2">Detection Reagent Troubleshooting</h4>
      <ul class="list-disc pl-4 space-y-2 opacity-80 leading-relaxed">
        <li><b>PE-A Reagent Failure:</b> Gavin's run (<code>_G</code>) used 555-conjugated anti-FLAG and showed low signal (&lt;1% Double+) across all samples. This points to a failure in his PE detection reagents or sample preparation.</li>
        <li><b>FITC-A Validation:</b> Using AF488 anti-mouse (FITC) for anti-FLAG binding and AF647 anti-rabbit (APC) for anti-c-myc expression confirmed high target activity.</li>
        <li><b>Spillover QC:</b> Gating and spillover correction worked as expected (Alpha = 0.1638), eliminating false-positive bleed-through between channels.</li>
      </ul>
    </div>
    <div class="p-3 bg-emerald-500/10 rounded border border-emerald-500/20 text-[10px]">
      <h4 class="text-emerald-400 font-bold text-[10px] uppercase tracking-widest mb-1">Active Target Conclusion</h4>
      <p class="leading-relaxed opacity-75">
        Both in-house purified ADAM10 FPLC fractions (F6 and F9) are folded and active, showing clear target binding and making them suitable for sorting and selectivity validations.
      </p>
    </div>
  </div>
</div>

---
layout: section
transition: slide-up
---

# Appendix: ADAM10/17 Expression & Purification
### FPLC and SDS-PAGE Optimization

---
layout: default
transition: fade-out
---

# FPLC Purification — ADAM10 & 17
Optimization of Anion Exchange Chromatography (HiTrap Q FF).

<div class="grid grid-cols-2 gap-4 mt-6">
  <div class="space-y-2">
    <div class="p-3 bg-white/5 rounded border border-white/10">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-2">ADAM10 (Sam's sample)</h4>
      <p class="text-[10px] leading-relaxed opacity-70">
        Run on AEC SAMPLE HiTrap Q FF preset. Only UV peak is in sample application.
      </p>
      <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/ADAM_10_Sam_AEC_1.png" class="w-full rounded mt-2 border border-white/10" style="max-height: 150px; object-fit: contain;">
    </div>
  </div>
  <div class="space-y-2">
    <div class="p-3 bg-white/5 rounded border border-white/10">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-2">ADAM10-pMopac</h4>
      <p class="text-[10px] leading-relaxed opacity-70">
        Diluted 1:4 with Buffer A. Something is eluting during wash. May need a more aggressive NaCl concentration during wash.
      </p>
      <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/ADAM_10_AEC_T1.png" class="w-full rounded mt-2 border border-white/10" style="max-height: 150px; object-fit: contain;">
    </div>
  </div>
</div>

---
# SDS-PAGE & Silver Stain Validation (June 2\textsuperscript{nd}, 2026)
Successful detection and optimization of ADAM10 expression and FPLC fractions.

<div class="grid grid-cols-2 gap-6 mt-4">
  <div class="p-4 bg-white/5 rounded border border-white/10 text-center">
    <b class="text-blue-400 text-xs block mb-1">ADAM10 Expression (Coomassie)</b>
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/ADAM1017_SDS-Page_Gel_20260602.png" class="w-full rounded mt-2 border border-white/10" style="max-height: 220px; object-fit: contain;">
    <p class="text-[10px] leading-snug opacity-70 mt-2">
      Expression test part 2. Bands visible at ~60-70 kDa corresponding to catalytic domain.
      Comparison with Jake's samples provides an independent expression benchmark.
    </p>
  </div>
  <div class="p-4 bg-white/5 rounded border border-white/10 text-center">
    <b class="text-emerald-400 text-xs block mb-1">Low-Abundance (Silver Stain)</b>
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/ADAM1017_SDS-Page_Gel_20260602_Silver.png" class="w-full rounded mt-2 border border-white/10" style="max-height: 220px; object-fit: contain;">
    <p class="text-[10px] leading-snug opacity-70 mt-2">
      Pierce Silver Stain of FPLC fractions. Confirmed presence in low-concentration samples.
      Detects F6 leakage in filtrate. Ladder compromised, so band identity is approximate.
    </p>
  </div>
</div>

---
# Protein Concentration & Recovery (June 2\textsuperscript{nd}, 2026)
Evaluation of FPLC fractions following spin concentration with Amicon Ultra-4 (10 kDa MWCO).

<div class="grid grid-cols-2 gap-4 mt-6">
  <div>
    <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-2">NanoDrop A280 Protein Concentrations</h4>
    <table class="w-full text-[9px] text-left border-collapse border border-white/10 rounded overflow-hidden">
      <thead>
        <tr class="bg-white/10 text-white font-bold border-b border-white/20">
          <th class="p-1 border-r border-white/10">Sample Name</th>
          <th class="p-1 border-r border-white/10">Protein (mg/mL)</th>
          <th class="p-1 border-r border-white/10">A280</th>
          <th class="p-1">A260/A280</th>
        </tr>
      </thead>
      <tbody>
        <tr class="border-b border-white/10 bg-white/5">
          <td class="p-1 border-r border-white/10 font-medium">A10 FXN 6 Unconcentrated</td>
          <td class="p-1 border-r border-white/10">0.163</td>
          <td class="p-1 border-r border-white/10">0.163</td>
          <td class="p-1">0.963</td>
        </tr>
        <tr class="border-b border-white/10 bg-red-500/10 text-red-200">
          <td class="p-1 border-r border-white/10 font-medium">A10 FXN 6 Filtrate (Leakage)</td>
          <td class="p-1 border-r border-white/10">0.176</td>
          <td class="p-1 border-r border-white/10">0.176</td>
          <td class="p-1">0.981</td>
        </tr>
        <tr class="border-b border-white/10 bg-emerald-500/10 text-emerald-200">
          <td class="p-1 border-r border-white/10 font-medium">A10 FXN 6 Residue (3.2x)</td>
          <td class="p-1 border-r border-white/10">0.519</td>
          <td class="p-1 border-r border-white/10">0.519</td>
          <td class="p-1">0.779</td>
        </tr>
        <tr class="border-b border-white/10 bg-white/5">
          <td class="p-1 border-r border-white/10 font-medium">A10 FXN 9 Unconcentrated</td>
          <td class="p-1 border-r border-white/10">0.141</td>
          <td class="p-1 border-r border-white/10">0.141</td>
          <td class="p-1">1.309</td>
        </tr>
        <tr class="border-b border-white/10 bg-white/5">
          <td class="p-1 border-r border-white/10 font-medium">A10 FXN 9 Filtrate</td>
          <td class="p-1 border-r border-white/10">0.063</td>
          <td class="p-1 border-r border-white/10">0.063</td>
          <td class="p-1">1.358</td>
        </tr>
        <tr class="border-b border-white/10 bg-emerald-500/10 text-emerald-200">
          <td class="p-1 border-r border-white/10 font-medium">A10 FXN 9 Residue (6.1x)</td>
          <td class="p-1 border-r border-white/10">0.863</td>
          <td class="p-1 border-r border-white/10">0.863</td>
          <td class="p-1">0.803</td>
        </tr>
        <tr class="bg-white/5">
          <td class="p-1 border-r border-white/10 font-medium">ADAM 17 FXN 12</td>
          <td class="p-1 border-r border-white/10">0.141</td>
          <td class="p-1 border-r border-white/10">0.141</td>
          <td class="p-1">0.930</td>
        </tr>
      </tbody>
    </table>
  </div>
  <div class="space-y-2">
    <div class="p-2.5 bg-emerald-500/10 rounded border border-emerald-500/20 text-[9px]">
      <h4 class="text-emerald-400 font-bold text-[10px] uppercase tracking-widest mb-1">Key Accomplishments</h4>
      <p class="leading-relaxed opacity-75">
        <b>Successful Concentration:</b> Spin concentration is highly effective for concentrating FPLC fractions. ADAM10 FXN 9 yielded a 6.1x increase in concentration (up to 0.863 mg/mL) with minimal loss in filtrate (0.063 mg/mL).
      </p>
    </div>
    <div class="p-2.5 bg-red-500/10 rounded border border-red-500/20 text-[9px]">
      <h4 class="text-red-400 font-bold text-[10px] uppercase tracking-widest mb-1">Leakage Alert</h4>
      <p class="leading-relaxed opacity-75">
        <b>ADAM10 FXN 6 Membrane Failure:</b> A10 FXN 6 showed significant protein leakage in the filtrate (0.176 mg/mL vs 0.163 mg/mL unconcentrated), indicating a membrane seal failure in that specific centrifugal unit.
      </p>
    </div>
    <div class="p-2.5 bg-white/5 rounded border border-white/10 text-[9px]">
      <h4 class="text-blue-300 font-bold text-[10px] uppercase tracking-widest mb-1">Future Recommendations</h4>
      <p class="leading-relaxed opacity-75">
        Check membrane integrity before run. Check yields of ADAM17 FXN 12 which remains low (0.141 mg/mL) and may require additional optimization.
      </p>
    </div>
  </div>
</div>

---
# ADAM17 Temperature Expansion Trial (June 3\textsuperscript{rd}, 2026)
Coomassie Blue SDS-PAGE of expression across induction temperatures and timepoints.

<div class="grid grid-cols-2 gap-6 mt-6 items-start">
  <div class="space-y-2 text-center">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/ADAM1017_SDS-Page_Gel_20260603.png" class="w-full rounded-lg border border-white/10" style="max-height: 280px; object-fit: contain;">
    <p class="text-[9px] opacity-40 italic mt-1">
      Coomassie Blue stained 12% Tris-glycine gel of ADAM17 expression (15 lanes).
    </p>
  </div>
  <div class="space-y-3">
    <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px]">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-1.5">Overnight Expression Requirement</h4>
      <p class="leading-relaxed opacity-75">
        <b>Strong Overnight Expression:</b> The overnight lysate pellet (Lane 13) shows strong expression with a prominent band at $\approx$30 kDa, indicating that overnight expression is required to obtain detectable levels.
      </p>
    </div>
    <div class="p-3 bg-emerald-500/10 rounded border border-emerald-500/20 text-[10px]">
      <h4 class="text-emerald-400 font-bold text-[10px] uppercase tracking-widest mb-1.5">Localization & Leakage QC</h4>
      <p class="leading-relaxed opacity-75">
        <b>No Periplasmic Leakage:</b> Target protein is absent in the growth supernatants (Lanes 7, 8, 9). The overnight growth pellet (Lane 12) shows expression, while early growth pellets (Lanes 10, 11) show none.
      </p>
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px]">
      <h4 class="text-amber-400 font-bold text-[10px] uppercase tracking-widest mb-1.5">Temperature Optimization</h4>
      <p class="leading-relaxed opacity-75">
        <b>25°C Standardization:</b> Early induction points (Lanes 14 and 15) show no bands. 25°C overnight induction is recommended to minimize inclusion body formation in the cell pellet, aligning with ADAM10 observations.
      </p>
    </div>
  </div>
</div>

---
# SDS-PAGE & Silver Stain Rerun (June 10\textsuperscript{th}, 2026)
Resolving migration issues with updated running gel chemistry and fresh buffers.

<div class="grid grid-cols-2 gap-4 mt-6">
  <div class="p-3 bg-white/5 rounded border border-white/10 text-center text-[9px]">
    <b class="text-blue-400 text-xs block mb-1">ADAM10 Gel Rerun</b>
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/ADAM1017_SDS-Page_Gel_20260610_ADAM10.png" class="w-full rounded mt-2 border border-white/10" style="max-height: 220px; object-fit: contain;">
    <p class="opacity-70 mt-2">
      Successful resolution of FPLC Fraction 9 Retentate (Lane 14) showing a band in the 30-40 kDa range (expected catalytic domain size).
    </p>
  </div>
  <div class="p-3 bg-white/5 rounded border border-white/10 text-center text-[9px]">
    <b class="text-emerald-400 text-xs block mb-1">ADAM17 Gel Rerun</b>
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/ADAM1017_SDS-Page_Gel_20260610_ADAM17.png" class="w-full rounded mt-2 border border-white/10" style="max-height: 220px; object-fit: contain;">
    <p class="opacity-70 mt-2">
      Successful resolution of FPLC Fraction 6 Retentate (Lane 14) resolving candidate target bands at 30-40 kDa with clean migration.
    </p>
  </div>
</div>

---
# Western Blot & Spin Concentration (June 11\textsuperscript{th}-12\textsuperscript{th}, 2026)
Validating expression, secretion, and concentration of ADAM10/17-pMopac target preparations.

<div class="grid grid-cols-2 gap-6 mt-6 items-start">
  <div class="space-y-3">
    <div class="p-3.5 bg-emerald-500/10 rounded border border-emerald-500/20 text-[10px]">
      <h4 class="text-emerald-400 font-bold text-[10px] uppercase tracking-widest mb-1.5">Centrifugal Spin Concentration</h4>
      <p class="leading-relaxed opacity-75">
        <b>20-40x Concentration:</b> Culture supernatants from expression trials (25°C, 30°C, and overnight) were concentrated using Amicon Ultra-4 (10 kDa MWCO) with no membrane failures or seal leakage.
      </p>
    </div>
    <div class="p-3.5 bg-blue-500/10 rounded border border-blue-500/20 text-[10px]">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-1.5">Western Blot Protocol QC</h4>
      <p class="leading-relaxed opacity-75">
        <b>Trans-Blot Turbo:</b> Blotted onto PVDF membrane using a custom 1.3A 24V 10Min protocol. Blocked in 5% milk / TBST-Brij, then incubated with mouse anti-FLAG (1:1000) primary antibody.
      </p>
    </div>
  </div>
  <div class="space-y-3">
    <div class="p-4 bg-white/5 rounded border border-white/10 text-[10px]">
      <h4 class="text-blue-300 font-bold text-[10px] uppercase tracking-widest mb-2">Western Blot Parameters</h4>
      <ul class="list-disc pl-4 space-y-2 opacity-80 leading-relaxed">
        <li><b>Lanes Loaded:</b> Ladder, A10 lysis pellet, A10 25°C/30°C/ON supernatants, A10 FXN 6, A17 25°C/30°C/ON supernatants, MMP2-FLAG (positive control).</li>
        <li><b>Blocking/Wash:</b> Blocked in 5% milk / TBST-Brij, primary anti-FLAG incubation overnight, secondary HRP-conjugated anti-mouse.</li>
        <li><b>Imaging:</b> Chemiluminescent ECL imaging and Ponceau S membrane total protein transfer staining completed on June 12th.</li>
      </ul>
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Western Blot Results & Analysis (June 12\textsuperscript{th}, 2026)
Validating soluble target expression and membrane transfer quality.

<div class="grid grid-cols-2 gap-4 mt-6">
  <div class="p-3 bg-white/5 rounded border border-white/10 text-center text-[9px]">
    <b class="text-blue-400 text-xs block mb-1">Ponceau S Membrane Stain</b>
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/ponceau_labeled.png" class="w-full rounded mt-2 border border-white/10" style="max-height: 200px; object-fit: contain;">
    <p class="opacity-70 mt-2">
      Uniform total protein transfer showing high-fidelity transfer of ladder and major bands (20260612).
    </p>
  </div>
  <div class="p-3 bg-white/5 rounded border border-white/10 text-center text-[9px]">
    <b class="text-emerald-400 text-xs block mb-1">Western Blot (Anti-FLAG)</b>
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/western_blot_labeled.png" class="w-full rounded mt-2 border border-white/10" style="max-height: 200px; object-fit: contain;">
    <p class="opacity-70 mt-2">
      Lane 3 (ADAM10 25°C supernatant) shows a single soluble band near 60 kDa; lanes 7-9 (ADAM17) show no signal.
    </p>
  </div>
</div>

---
layout: default
transition: fade-out
---

# BCA Assay — ADAM10 Concentration Standardization (June 12\textsuperscript{th}, 2026)
Quantifying purified ADAM10 FPLC Fraction 6 for downstream binding trials.

<div class="grid grid-cols-2 gap-4 mt-4">
  <div class="space-y-3">
    <div class="p-3 bg-blue-500/10 rounded border border-blue-500/20 text-[10px]">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-1.5">SpectraMax i3x Run</h4>
      <p class="leading-relaxed opacity-75">
        <b>Quadratic Fit ($R^2 = 0.9895$):</b> Back-calculated neat, 1:2, and 1:5 dilutions yield a mean concentration of <b>205.46 μg/mL</b> (~0.21 mg/mL).
      </p>
      <div class="flex gap-2 mt-2">
        <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/20260612_2_standard_curve.png" class="w-1/2 rounded border border-white/5" style="max-height: 80px; object-fit: contain;">
        <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/20260612_2_heatmap.png" class="w-1/2 rounded border border-white/5" style="max-height: 80px; object-fit: contain;">
      </div>
    </div>
  </div>
  <div class="space-y-3">
    <div class="p-3 bg-emerald-500/10 rounded border border-emerald-500/20 text-[10px]">
      <h4 class="text-emerald-400 font-bold text-[10px] uppercase tracking-widest mb-1.5">SPECTRAmax M5 Run</h4>
      <p class="leading-relaxed opacity-75">
        <b>Linear Fit ($R^2 = 0.9930$):</b> Dilution adjusted values show a mean concentration of <b>193.89 μg/mL</b> (~0.19 mg/mL).
      </p>
      <div class="flex gap-2 mt-2">
        <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/BCA_template_2_standard_curve.png" class="w-1/2 rounded border border-white/5" style="max-height: 80px; object-fit: contain;">
        <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/BCA_template_2_heatmap.png" class="w-1/2 rounded border border-white/5" style="max-height: 80px; object-fit: contain;">
      </div>
    </div>
  </div>
</div>

<div class="mt-4 p-3 bg-white/5 rounded border border-white/10 text-[10px] text-center">
  <span class="text-emerald-400 font-bold">Standardized Concentration:</span> Purified ADAM10 FPLC Fraction 6 concentration is standardized to <b>0.20 mg/mL</b> for downstream display trials.
</div>

---
layout: center
class: text-center
---

# Questions
### De Novo Binder Generation — RFdiffusion to Flow Cytometry

<div class="mt-8 grid grid-cols-3 gap-6 max-w-2xl mx-auto text-center">
  <div>
    <div class="text-emerald-400 font-black text-2xl">75%</div>
    <div class="text-[9px] opacity-40 uppercase tracking-widest mt-1">Positive Predictive Value (3/4)</div>
  </div>
  <div>
    <div class="text-blue-400 font-black text-2xl">3.8×</div>
    <div class="text-[9px] opacity-40 uppercase tracking-widest mt-1">Peak Selectivity (AB 6)</div>
  </div>
  <div>
    <div class="text-violet-400 font-black text-2xl">>98%</div>
    <div class="text-[9px] opacity-40 uppercase tracking-widest mt-1">Candidate Reduction</div>
  </div>
</div>

<div class="mt-8 text-xs opacity-30 uppercase tracking-widest">
  Repository: PhD-Research / Demonstrations / De_Novo_Binder_Generation
</div>
