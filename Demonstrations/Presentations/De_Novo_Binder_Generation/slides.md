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
### Generative Design → Structural Verification → Yeast Display Validation

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

# Full Pipeline Architecture — Round 1
Six stages from backbone hallucination to experimentally validated selectivity. This is the pipeline that generated and validated the primary result of this talk; <b class="text-amber-300">Round 2</b> (an iterative-refinement successor, now 12 iterations deep) is covered later.

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
      <p><b class="text-cyan-300">Yeast surface display:</b> Ready-to-transform pCHA-TIMP3 plasmids received directly from Twist Bioscience (no in-lab cloning step); EBY100 yeast transformation. Surface expression detected by FITC-conjugated anti-Myc antibody. His-tagged target detected by APC-conjugated anti-His.</p>
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
Ready-to-transform pCHA-TIMP3 plasmids from Twist Bioscience, spanning three design intents.

<div class="mt-6">
  <TopVariants />
</div>

---
layout: default
transition: fade-out
---

# Final Library — Design-Intent Breakdown
Three parallel design axes in the same synthesis order — not just an ADAM17 side quest.

<div class="grid grid-cols-3 gap-5 mt-6">
  <div class="p-4 bg-emerald-500/10 rounded border border-emerald-500/20 space-y-2">
    <h4 class="text-emerald-300 font-bold text-[11px] uppercase tracking-widest">MMP9-Selective (5)</h4>
    <p class="text-[10px] leading-relaxed opacity-70">
      AB 1, AB 2, AB 6 (AB loop) and C 12, C 15 (C loop) — designed for M9&gt;M2 preference. AB 6, C 12, C 15 confirmed by ANOVA + Tukey-HSD; AB 1/AB 2 directionally correct but underpowered. The primary result of this campaign.
    </p>
  </div>
  <div class="p-4 bg-blue-500/10 rounded border border-blue-500/20 space-y-2">
    <h4 class="text-blue-300 font-bold text-[11px] uppercase tracking-widest">MMP9 High / Low Controls (3)</h4>
    <p class="text-[10px] leading-relaxed opacity-70">
      AB 3 ("High" — broad binder) and AB 5, C 13 ("Low" — designed non-binders) calibrate the assay's dynamic range on the MMP9 axis. Both Low constructs came back correctly non-significant — validating that the T-score filter isn't just calling everything a hit.
    </p>
  </div>
  <div class="p-4 bg-amber-500/10 rounded border border-amber-500/20 space-y-2">
    <h4 class="text-amber-300 font-bold text-[11px] uppercase tracking-widest">ADAM17-Targeted (5)</h4>
    <p class="text-[10px] leading-relaxed opacity-70">
      AB 4, AB 7, C 11, C 14, ABC 22 designed for ADAM17 (and A17&gt;A10 where noted). ADAM10 comparison was limited by positive-control activity in this campaign; AB 4/AB 7/C 11 later confirmed as reproducible MMP9 hits in a follow-up screen.
    </p>
  </div>
</div>

<div class="mt-6 p-3 bg-white/5 rounded border border-white/10 text-[10px] opacity-70">
  <b class="text-blue-400">Library design rationale:</b> The MMP9 axis (selective + High/Low controls) is the primary, statistically powered result of this campaign, spanning both contact loops (AB and C). The ADAM17 group runs in parallel to test pipeline generalizability beyond the MMP scaffold in the same synthesis order.
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

# Flow Cytometry Metrics
Six quantitative metrics extracted per sample.

<div class="mt-4">
  <MetricExplorer />
</div>

<div class="mt-4 grid grid-cols-3 gap-2 text-[9px] opacity-50">
  <div><span class="text-cyan-400 font-bold">Cross-target:</span> Bind Med (Expr+), Binding Efficiency</div>
  <div><span class="text-violet-400 font-bold">Within-target:</span> Norm Bind Med, Norm Median Ratio, Norm IWB</div>
  <div><span class="text-slate-400 font-bold">QC only:</span> Stain Index</div>
</div>

<div class="mt-2 text-[8px] opacity-30 italic text-center">Hover any card above for its full definition, scope, and use.</div>

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
  <div class="space-y-2">
    <BindMedBars />
    <p class="text-[9px] opacity-40 italic text-center">
      Bind Med (Expr+) per construct across all five targets — mean ± SEM, aggregated from real per-trial data. Hover a bar for exact values.
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
    <NormBindHeatmap />
  </div>
  <div class="space-y-2">
    <div class="p-2 bg-white/5 rounded border border-white/10 text-[9.5px] leading-snug opacity-80">
      <b class="text-blue-400">Why normalize within-target:</b> Raw MMP9 and MMP2 baselines differ due to protein concentration, staining efficiency, and TIMP3 affinity per target. Dividing by TIMP3-WT for each target independently isolates construct-level effects from between-target scale differences.
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
75% PPV — every non-selective control stayed non-significant.

<div class="mt-3">
  <SelectivityBars />
</div>

<div class="grid grid-cols-3 gap-4 mt-3 text-[10px]">
  <div class="p-2.5 bg-emerald-500/10 rounded border border-emerald-500/20">
    <b class="text-emerald-300">Confirmed selective:</b> C 12 (p=0.012), C 15 (p=0.018), AB 6 (p=0.0002) — ratios 2.6–3.8× vs. TIMP3-WT (2.5×).
  </div>
  <div class="p-2.5 bg-blue-500/10 rounded border border-blue-500/20">
    <b class="text-blue-300">Zero false positives:</b> All non-selective controls (C 13, AB 5, AB 2, AB 3) stayed non-significant (p&gt;0.05).
  </div>
  <div class="p-2.5 bg-violet-500/10 rounded border border-violet-500/20">
    <b class="text-violet-300">Scaffold plasticity:</b> 6–15 aa insertions tolerated; WT TIMP3 already shows mild MMP9 preference (p=0.024) — variants amplify it.
  </div>
</div>

<div class="mt-2 text-[9px] opacity-50 text-center">
  75% PPV at N=4 with 0 false positives justifies SPR kinetics next — this is the narrow, adequately-powered confirmatory test; see the full 13-construct scorecard ahead for the broader picture.
</div>

---
layout: default
transition: fade-out
---

# Design-Prediction Verdict Scorecard — All 13 Constructs

<div class="grid grid-cols-2 gap-8 mt-6 items-start">
  <div class="space-y-4">
    <div class="p-4 bg-white/5 rounded border border-white/10">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-3">Model vs. Bench Scorecard</h4>
      <p class="text-[11px] leading-relaxed opacity-75">
        A systematic tally comparing the design intent of all 13 engineered TIMP3 loop variants against aggregated flow-cytometry results — the full-library complement to the MMP9-vs-MMP2 primary result on the next slide.
      </p>
      <ul class="text-[10px] list-disc pl-4 space-y-1.5 opacity-70 mt-3">
        <li><b>Hits (6/13):</b> Binding trend matched computational design intent. Standout: <b>AB 4</b> (ADAM17-directed, normalized median ratio 1.26 — ranked 1st of 12 testable constructs).</li>
        <li><b>Partials (2/13):</b> Directionally correct but not statistically significant.</li>
        <li><b>Misses (2/13):</b> <b>C 13</b> predicted low-affinity but measured NMR 1.11; <b>C 15</b> predicted MMP9-preferential but measured MMP2-preferential on this broader (non-ANOVA) metric — flags scoring-function edge cases even though C 15 remains one of the three ANOVA-confirmed MMP9-selective winners (next slide).</li>
        <li><b>Untestable (3/13):</b> Limited by low target activity, degraded prep quality, or mismatch in screening thresholds (incl. C 16 / ABC 21, which failed to display on yeast at all).</li>
      </ul>
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px]">
      <b class="text-amber-400">Reading this against the headline number:</b> The 75% PPV / "3 of 4" statistic quoted throughout this talk is the <i>narrow, adequately-powered</i> confirmatory result — the MMP9-vs-MMP2 ANOVA on 4 constructs explicitly designed for that axis. This scorecard is the <i>broad</i> exploratory result — directional agreement across all 13 constructs on every design axis (MMP9, ADAM17, "Low") and every metric, most of which were underpowered (n=1–3 trials). Both are true; they answer different questions.
    </div>
  </div>
  <div class="space-y-3 flex flex-col items-center">
    <CategoryBars dataset="verdict" title="Verdict Tally — All 13 Constructs" />
    <p class="text-[9px] opacity-40 italic text-center">
      Tally of design campaign success rate across all 13 scored constructs.
    </p>
  </div>
</div>

---
layout: section
transition: slide-up
---

# Phase 3
### Second-Generation Design — Iterative Refinement & ESM-C (June–July 2026)

---
layout: default
transition: fade-out
---

# Round 2 Pipeline: Iterative Refinement
Round 1 validated the approach; Round 2 automates and accelerates it on a dedicated GPU cluster.

<div class="grid grid-cols-5 gap-2 mt-4">
  <div class="flex flex-col items-center p-2.5 bg-blue-500/10 rounded border border-blue-500/30 text-center">
    <div class="text-blue-400 font-black text-lg mb-1">①</div>
    <div class="text-[10px] font-bold text-white mb-1.5">RFd3</div>
    <div class="text-[8px] opacity-50 leading-tight">RoseTTAFold Diffusion 3<br/>20 backbones/target/iter</div>
  </div>
  <div class="flex items-center justify-center text-blue-400 text-lg">→</div>
  <div class="flex flex-col items-center p-2.5 bg-violet-500/10 rounded border border-violet-500/30 text-center">
    <div class="text-violet-400 font-black text-lg mb-1">②</div>
    <div class="text-[10px] font-bold text-white mb-1.5">LigandMPNN</div>
    <div class="text-[8px] opacity-50 leading-tight">Annealed T: 0.50→0.10<br/>0.85× decay/iteration</div>
  </div>
  <div class="flex items-center justify-center text-blue-400 text-lg">→</div>
  <div class="flex flex-col items-center p-2.5 bg-cyan-500/10 rounded border border-cyan-500/30 text-center">
    <div class="text-cyan-400 font-black text-lg mb-1">③</div>
    <div class="text-[10px] font-bold text-white mb-1.5">ESMFold2</div>
    <div class="text-[8px] opacity-50 leading-tight">MSA-free local ranker<br/>replaced Boltz-2</div>
  </div>
  <div class="flex items-center justify-center text-blue-400 text-lg">→</div>
  <div class="flex flex-col items-center p-2.5 bg-emerald-500/10 rounded border border-emerald-500/30 text-center">
    <div class="text-emerald-400 font-black text-lg mb-1">④</div>
    <div class="text-[10px] font-bold text-white mb-1.5">Best Binders</div>
    <div class="text-[8px] opacity-50 leading-tight">Top 75/target<br/>loop-length narrowing</div>
  </div>
  <div class="flex items-center justify-center text-blue-400 text-lg">→</div>
  <div class="flex flex-col items-center p-2.5 bg-amber-500/10 rounded border border-amber-500/30 text-center">
    <div class="text-amber-400 font-black text-lg mb-1">⑤</div>
    <div class="text-[10px] font-bold text-white mb-1.5">AF3 Gold-Standard</div>
    <div class="text-[8px] opacity-50 leading-tight">Every 2 iterations<br/>final verification</div>
  </div>
</div>

<div class="mt-4 grid grid-cols-2 gap-3">
  <div class="p-2.5 bg-white/5 rounded border border-white/10 text-[9px] opacity-75 leading-snug">
    <b class="text-blue-400">Why ESMFold2 over Boltz-2:</b> Predicts AF3 ipTM at ρ≈+0.34 (n=24) vs. Boltz-2's ρ≈+0.11 — a better cheap-ranker proxy for the expensive AF3 gold standard.
  </div>
  <div class="p-2.5 bg-red-500/10 rounded border border-red-500/20 text-[9px] opacity-75 leading-snug">
    <b class="text-red-400">RF3 disabled for scoring:</b> RFd3's own confidence scores are <i>anti-correlated</i> with AF3 ipTM (r≈-0.07 to -0.51, n=18) — contributes geometric features only.
  </div>
</div>

<div class="mt-2.5 p-2.5 bg-white/5 rounded border border-white/10 text-[9px] opacity-70 leading-snug">
  Compute: Caprine cluster, 2× V100 GPUs. Active targets: MMP2, MMP9, MMP3, MMP10, ADAM10, ADAM17. Selection script <code>select_binders_to_order.py</code> computes loop-pLDDT, loop-interface PAE, and framework RMSD vs. native TIMP3.
</div>

---
layout: default
transition: fade-out
---

# Per-Target Champions & Convergence Dynamics (Iteration 12)
Cumulative best composite score discovered per iteration, tracked across all four protease targets.

<div class="grid grid-cols-2 gap-6 mt-4 items-start">
  <div>
    <ConvergenceLines />
    <p class="text-[9px] opacity-40 italic text-center mt-1">
      Cumulative max composite score through iteration 12, by target. Hover a line for its champion score.
    </p>
  </div>
  <div class="space-y-2 text-[10px]">
    <table class="w-full font-mono border-collapse text-[9px]">
      <thead>
        <tr class="border-b border-white/20 text-blue-400 font-bold">
          <th class="text-left py-1">Target</th>
          <th class="text-left py-1">Champion</th>
          <th class="text-right py-1">Composite</th>
          <th class="text-right py-1">Iter</th>
        </tr>
      </thead>
      <tbody class="opacity-85">
        <tr class="border-b border-white/10"><td class="py-1">ADAM10</td><td>it7_d32_s0</td><td class="text-right">0.840</td><td class="text-right">7</td></tr>
        <tr class="border-b border-white/10"><td class="py-1 text-emerald-400">ADAM17</td><td>it8_d29_s0</td><td class="text-right text-emerald-400">0.831</td><td class="text-right">8</td></tr>
        <tr class="border-b border-white/10"><td class="py-1">MMP9</td><td>it6_d11_s1</td><td class="text-right">0.808</td><td class="text-right">6</td></tr>
        <tr><td class="py-1 text-amber-400">MMP2</td><td>it3_d18_s0</td><td class="text-right text-amber-400">0.762</td><td class="text-right">3</td></tr>
      </tbody>
    </table>
    <div class="p-2.5 bg-emerald-500/10 rounded border border-emerald-500/20 mt-2">
      <b class="text-emerald-400">ADAM17 still climbing:</b> steady late-stage gain (0.822→0.827→0.831 across iterations 1/4/8) — the most viable candidate for further cycles.
    </div>
    <div class="p-2.5 bg-amber-500/10 rounded border border-amber-500/20">
      <b class="text-amber-400">MMP2 plateaued:</b> peaked at iteration 3 with no improvement through iteration 12 — flags a scaffold or parameter rethink rather than "just run it longer."
    </div>
    <div class="p-2.5 bg-white/5 rounded border border-white/10">
      <b class="text-blue-300">MMP9 has the deepest pool:</b> median 0.731, worst-retained 0.675 among the top-75 Hall of Fame — all four targets emerged at iterations 3–8, confirming the loop consistently beats the seed pool.
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Raw ipTM Calibration Looks Broken (First Pass)

<div class="grid grid-cols-2 gap-8 mt-6 items-start">
  <div class="space-y-4">
    <div class="p-4 bg-white/5 rounded border border-white/10">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-3">ipTM vs. Yeast Display NMR — raw, per target</h4>
      <p class="text-[11px] leading-relaxed opacity-75">
        Linear correlation of computational ipTM against experimental yeast-display normalized median ratio (NMR), 12 constructs across three targets.
      </p>
      <ul class="text-[10px] list-disc pl-4 space-y-1.5 opacity-70 mt-3">
        <li><b>ADAM17:</b> Weak negative (r = -0.163, p = 0.613).</li>
        <li><b>MMP2:</b> Essentially zero (r = 0.029, p = 0.928).</li>
        <li><b>MMP9:</b> Moderate negative (r = -0.471, p = 0.123).</li>
      </ul>
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px]">
      <b class="text-amber-400">First-pass read (later revised, next slide):</b> Taken at face value, ipTM appears to carry no predictive signal for binding. That conclusion turns out to be an artifact of testing raw scores — see why on the next slide.
    </div>
  </div>
  <div class="space-y-3 flex flex-col items-center">
    <CalibrationScatter metric="ipTM" />
    <p class="text-[9px] opacity-40 italic text-center">
      ipTM predictions vs. experimental NMR, real per-construct data (n=12). Hover a point for its value.
    </p>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Resolving the Calibration Signal: An Avidity Confound
Quantitative calibration of AlphaFold3 metrics against binding, 12 constructs × 3 targets (ADAM17, MMP2, MMP9).

<div class="grid grid-cols-2 gap-6 mt-4 items-start">
  <div class="space-y-3">
    <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px] leading-relaxed opacity-80">
      <b class="text-blue-400">Binding is dominated by a non-specific "avidity" factor:</b> Binding correlates strongly <em>across targets</em> for a given construct (ADAM17–MMP2 r=0.72, ADAM17–MMP9 r=0.69) — a construct that binds one protease tends to bind them all. Two-way variance decomposition attributes <b>49%</b> of binding variance to this construct-level "stickiness," <b>17%</b> to target baseline, and only <b>33%</b> to genuine target-specific interaction.
    </div>
    <div class="p-3 bg-emerald-500/10 rounded border border-emerald-500/20 text-[10px] leading-relaxed opacity-80">
      <b class="text-emerald-400">Once avidity is regressed out, the signal reappears:</b> Per-target interface metrics predict the target-specific residual at oriented Spearman ρ ≈ 0.50 (interface PAE), 0.47 (loop pLDDT), 0.44 (BpTM) — far above their raw correlation. The categorical selectivity result (AB 6, C 12, C 15) was real all along; it was just masked in raw scores by stickiness.
    </div>
    <div class="p-3 bg-red-500/10 rounded border border-red-500/20 text-[10px] leading-relaxed opacity-80">
      <b class="text-red-400">ipTM is a developability signal, not an affinity signal:</b> ipTM correlates more with construct <em>expression</em> (ρ ≈ 0.36) than with binding (ρ ≈ -0.09), and is nearly saturated across these designs (0.77–0.90) — too compressed to rank fine differences. Only loop-pLDDT and interface PAE carry real dynamic range.
    </div>
  </div>
  <div class="space-y-2 flex flex-col items-center w-full">
    <CategoryBars dataset="variance" title="Variance Decomposition (%)" />
    <CategoryBars dataset="rho" title="Avidity-Removed Correlation (ρ)" />
    <p class="text-[9px] opacity-40 italic text-center">
      Top: 49% avidity / 17% target baseline / 33% target-specific. Bottom: once avidity is removed, PAE/LpLDDT/BpTM track binding; ipTM tracks expression, not binding.
    </p>
  </div>
</div>

---
layout: default
transition: fade-out
---

# A Multi-Term, Re-Fittable Selection Recipe
Combining the real-dynamic-range metrics into a transparent, directional scoring function for the next design cycle.

<div class="grid grid-cols-2 gap-6 mt-4 items-start">
  <div class="space-y-3">
    <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px] leading-relaxed opacity-80">
      <b class="text-blue-400">Recipe:</b> (i) an affinity term from loop-pLDDT + inverted interface PAE, (ii) an off-target selectivity term (on-target minus mean off-target confidence), (iii) BpTM as a low-weight term, (iv) ipTM as a small expressibility prior. Weights fixed by reasoning, not least-squares.
    </div>
    <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px] leading-relaxed opacity-80">
      <b class="text-cyan-400">Why not just fit the weights?</b> At n=12, a leave-one-construct-out regression-fit blend performed <em>worse</em> (ρ=0.20) than the reasoned weights (ρ=0.31) — fitting overfits at this scale. The recipe out-ranks every individual metric, including ipTM (currently the funnel's optimization target).
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px] leading-relaxed opacity-80">
      <b class="text-amber-400">Honest uncertainty:</b> 95% bootstrap CI on ρ=0.31 is [-0.15, +0.71] — it crosses zero. This is a reasoned prior that beats any single metric today, <em>not</em> a precisely estimated model. It is designed to be re-estimated as every new ordered batch returns data.
    </div>
  </div>
  <div class="space-y-3 flex flex-col items-center">
    <CategoryBars dataset="recipe" title="Ranking Power (ρ) vs. Individual Metrics" />
    <p class="text-[9px] opacity-40 italic text-center">
      Recipe score vs. consensus binding (ρ=0.31, wide bootstrap CI [-0.15, 0.71]) — out-ranks every individual metric and the overfit data-fitted blend.
    </p>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Are the Docked Poses Even Right? A Crystal-Free Check
The backbone-conditioning step docks TIMP3:target with HADDOCK before hallucination — but co-fold models may simply "remember" TIMP:MMP complexes from training. Two tests that consult no crystal structure probe whether that matters.

<div class="grid grid-cols-2 gap-6 mt-4 items-start">
  <div class="space-y-3">
    <div class="p-3 bg-emerald-500/10 rounded border border-emerald-500/20 text-[10px] leading-relaxed opacity-80">
      <b class="text-emerald-400">Test 1 — Mechanistic chelation geometry:</b> AF3/ESMFold2 co-folds reproduce the native crystal's edge-leading C1/T2 chelation geometry (d≈2.6–4.0 Å vs. crystal 2.9 Å, edge-leads fraction 1.0). HADDOCK tracks bury the zinc loop non-specifically (d≈5.6–7.7 Å, edge-leads only 0.39–0.61).
    </div>
    <div class="p-3 bg-emerald-500/10 rounded border border-emerald-500/20 text-[10px] leading-relaxed opacity-80">
      <b class="text-emerald-400">Test 2 — Cross-source self-consistency:</b> The two co-folds mutually agree to ≈5.7 Å construct Cα RMSD after target superposition. The four independent HADDOCK tracks disagree with <em>each other</em> by ≈14–22 Å — a method that can't reproduce its own answer across equivalent inputs is unreliable on its own terms.
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px] leading-relaxed opacity-80">
      <b class="text-amber-400">Practical takeaway:</b> Attempts to "fix" HADDOCK (catalytic Zn²⁺, conformer ensembles) stayed CAPRI-incorrect. Use co-folds for binding-mode assignment; reserve HADDOCK for energetics on poses already believed correct — structural models are a foldability/mechanistic filter, not a standalone ranking objective.
    </div>
  </div>
  <div class="space-y-3 flex flex-col items-center">
    <CategoryBars dataset="haddock" title="Construct Cα RMSD Across Sources (Å)" />
    <p class="text-[9px] opacity-40 italic text-center">
      Cross-source pose convergence: co-folds agree (~5.7 Å); HADDOCK tracks disagree with each other (~14–22 Å) — without invoking any ground-truth structure.
    </p>
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
    <CrossReactivityBars />
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

# Yeast Display Target Trial — MMP9 & ADAM10
Flow Cytometry assessment of TwistBio constructs against commercial and in-house targets.

<div class="grid grid-cols-2 gap-5 mt-4 items-start">
  <div class="space-y-2">
    <div class="p-2 bg-emerald-500/10 rounded border border-emerald-500/20 text-[9px] leading-snug">
      <h4 class="text-emerald-400 font-bold text-[9px] uppercase tracking-widest mb-0.5">MMP9 Positives Validated (Purchased Sino)</h4>
      <p class="opacity-75">
        Sino MMP9: positive control TIMP3 38.2% Double+ (491 MFI); AB 7 and AB 3 achieved 33.5%/32.9% Double+, validating the display system.
      </p>
    </div>
    <div class="p-2 bg-violet-500/10 rounded border border-violet-500/20 text-[9px] leading-snug">
      <h4 class="text-violet-400 font-bold text-[9px] uppercase tracking-widest mb-0.5">Commercial AbCam ADAM10 Validation</h4>
      <p class="opacity-75">
        AbCam ADAM10 (purchased): TIMP3 positive control 12.65% Double+; AB 5 and AB 3 showed strong binding (MFI 1708/1665.5).
      </p>
    </div>
  </div>
  <div class="space-y-2">
    <div class="p-2 bg-red-500/10 rounded border border-red-500/20 text-[9px] leading-snug">
      <h4 class="text-red-400 font-bold text-[9px] uppercase tracking-widest mb-0.5">In-House ADAM10 (F6/F9): Retracted</h4>
      <p class="opacity-75">
        This June 10 channel-corrected read (AB 7 13.29%, AB 1 11.71% Double+ vs. in-house Fraction 6/9) is <b>no longer trusted</b>. The construct-QC vendor manifest flags every in-house ADAM10 prep (Sam-FLAG, F6/F9) as out of calibration scope, and a dedicated July 24 control test (in-house ADAM10 vs. displayed WT TIMP1/TIMP3, proper NC) found <b>no genuine binding signal</b> — the earlier read is now understood as background/artifact, not confirmed binding.
      </p>
    </div>
    <div class="p-2 bg-amber-500/10 rounded border border-amber-500/20 text-[9px] leading-snug">
      <b class="text-amber-400">Status:</b> ADAM10 binding remains <b>unresolved</b>. The commercial AbCam reagent (left) is the only ADAM10 source currently in calibration scope; the in-house axis needs a validated positive control and a titrated, quantified target before any binding or selectivity claim can be made.
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Next Steps and Translational Vision

<div class="grid grid-cols-2 gap-6 mt-4">
  <div class="space-y-2">
    <h3 class="text-cyan-400 font-bold uppercase text-[11px] tracking-widest border-b border-cyan-500/20 pb-1">Near-Term Experiments</h3>
    <ul class="text-[10px] space-y-2 opacity-80 list-none p-0 leading-snug">
      <li class="flex items-start gap-2">
        <div class="w-1.5 h-1.5 rounded-full bg-cyan-400 mt-1 shrink-0"></div>
        <span><b>SPR Kinetics (queued, top priority):</b> K<sub>D</sub>/k<sub>on</sub>/k<sub>off</sub> for C 12, C 15, AB 6 vs. purified MMP9/MMP2 — direct comparison to marimastat/prinomastat.</span>
      </li>
      <li class="flex items-start gap-2">
        <div class="w-1.5 h-1.5 rounded-full bg-cyan-400 mt-1 shrink-0"></div>
        <span><b>Reverse Selectivity — stalled:</b> Round 2's MMP2 champion plateaued at iteration 3; needs a scaffold/parameter rethink, not more iterations.</span>
      </li>
      <li class="flex items-start gap-2">
        <div class="w-1.5 h-1.5 rounded-full bg-emerald-400 mt-1 shrink-0"></div>
        <span><b>ADAM10/17 — still unresolved:</b> In-house ADAM10 <i>production</i> is validated (0.20 mg/mL, Western+BCA), but its <i>binding</i> is not — the earlier "positive" read was retracted by a July 24 control test and the vendor manifest. Needs a validated positive control. ADAM17 construct separately came back an empty vector — needs re-cloning (details: companion <i>ADAM Target Production</i> deck).</span>
      </li>
      <li class="flex items-start gap-2">
        <div class="w-1.5 h-1.5 rounded-full bg-cyan-400 mt-1 shrink-0"></div>
        <span><b>ESM-C — exploratory only (appendix):</b> Consensus motifs found could nominate a future batch, but the calibrated multi-term recipe remains the active selection method for now.</span>
      </li>
    </ul>
  </div>

  <div class="space-y-2">
    <h3 class="text-violet-400 font-bold uppercase text-[11px] tracking-widest border-b border-violet-500/20 pb-1">Translational Context</h3>
    <div class="p-3 bg-violet-500/5 rounded border border-violet-500/20 text-[9px] leading-snug opacity-80 space-y-2">
      <p><b class="text-violet-300">MMP9 as oncology target:</b> Overexpressed in >15 tumor types. Prior pan-MMP inhibitors failed Phase III on MMP2 off-target toxicity, not a wrong hypothesis — insufficient selectivity. An MMP9-specific variant addresses that failure mode directly.</p>
      <p><b class="text-violet-300">Target-agnostic framework:</b> Any protein–protein interaction with structural data is addressable; MMP9/MMP2 is the validation case, but the workflow is identical for kinase isoforms, cytokine receptors, or viral surface proteins.</p>
      <p><b class="text-violet-300">Precision CAR-T analogy:</b> Tumor genome → neoantigen ID → de novo binder → CAR construct with genome-informed recognition, minimizing on-target/off-tumor toxicity.</p>
    </div>
  </div>
</div>

---
layout: section
transition: slide-up
---

# Appendix
### ESM-C Sequence Classifier — Exploratory, Not Yet a Production Ranker

---
layout: default
transition: fade-out
---

# ESM-C Sequence Classifier (Exploratory)
An early-stage, structure-free classifier explored as a possible future sequence-level pre-filter — not currently used to rank or select designs in the active pipeline.

<div class="grid grid-cols-2 gap-6 mt-6 items-start">
  <div class="space-y-3">
    <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px] leading-relaxed opacity-80">
      <b class="text-blue-400">Scale:</b> A fine-tuned 300M-parameter ESM-C classifier enumerated all ~64M possible 6-residue C-loop sequences; retained the top 50,000 predicted binders per target (P≥0.99) for ADAM17, MMP3, MMP9. MMP9 and MMP3 share 38.1% Jaccard overlap; ADAM17 is highly distinct (0.066).
    </div>
    <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px] leading-relaxed opacity-80">
      <b class="text-cyan-400">Consensus motifs found:</b> A shared core <code>L-S-x-x-T</code> (positions 1,2,5) across all three targets; selectivity is carried by positions 3/4/6 — ADAM17 <code>LSSDTT</code>, MMP3 <code>LSPDTT</code>, MMP9 <code>LSPTTL</code>. These are candidate leads for a future order, not yet acted on.
    </div>
  </div>
  <div class="flex flex-col items-center">
    <EsmcConfidenceBars />
  </div>
</div>

---
layout: default
transition: fade-out
---

# ESM-C: Early Signal, Not Yet Sufficient to Adopt
Held-out classifier performance, plus a direct check against flow-cytometry results — encouraging, but underpowered.

<div class="grid grid-cols-2 gap-6 mt-6 items-start">
  <div class="space-y-3">
    <div class="p-3 bg-emerald-500/10 rounded border border-emerald-500/20 text-[10px] leading-relaxed opacity-80">
      <b class="text-emerald-400">Held-out performance:</b> ROC-AUC 0.912, PR-AUC 0.746, MCC 0.722, F1 0.734 (n=6,788, 381 MMP9 positives). ~99.7% of top predictions are absent from training data — genuine generalization, not memorization.
    </div>
    <div class="p-3 bg-violet-500/10 rounded border border-violet-500/20 text-[10px] leading-relaxed opacity-80">
      <b class="text-violet-400">Encouraging early check against real FCS binding:</b> Spearman ρ=0.86 (ADAM17, p=0.014, n=7) and ρ=0.75 (MMP9 Double+%, p=0.052, n=7) on out-of-distribution AB-loop constructs — outperforming monomeric AlphaFold LpLDDT (r≈-0.26) on the same constructs.
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px] opacity-80">
      <b class="text-amber-400">Why it's not adopted yet:</b> n=7 is far too small to act on, and the strict novel-loop held-out subset contains 0 positive binders, so out-of-distribution generalization can't be quantified statistically. Wet-lab validation on a dedicated batch is required before ESM-C output could drive a synthesis order.
    </div>
  </div>
  <div class="flex flex-col items-center">
    <CategoryBars dataset="esmc" title="Held-Out Test Performance (MMP9)" />
    <p class="text-[9px] opacity-40 italic text-center mt-1">Held-out ROC/PR performance on MMP9.</p>
  </div>
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
