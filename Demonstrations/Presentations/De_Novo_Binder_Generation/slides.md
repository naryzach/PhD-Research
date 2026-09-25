---
theme: seriph
title: "De Novo Binder Generation"
info: |
  ## De Novo Protein Prediction to Binding Evaluation
  Engineering TIMP3 loop variants for MMP9 preference via generative deep learning.
  
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
    <div class="text-[10px] uppercase opacity-40 tracking-widest">Supported Selective (vendor-matched, one day)</div>
    <div class="text-blue-400 font-bold uppercase text-xs">AB 1, AB 2, AB 6, C 12, C 15 — MMP9 &gt; MMP2 (nominal p=0.008–0.044)</div>
  </div>
  <div class="px-4 py-2 bg-white/5 rounded border border-white/10">
    <div class="text-[10px] uppercase opacity-40 tracking-widest">Primary Targets</div>
    <div class="text-blue-400 font-bold uppercase text-xs">MMP9 vs MMP2</div>
  </div>
  <div class="px-4 py-2 bg-white/5 rounded border border-white/10">
    <div class="text-[10px] uppercase opacity-40 tracking-widest">Caveat</div>
    <div class="text-emerald-400 font-bold uppercase text-xs">n=2v2, one day; replication pending</div>
  </div>
</div>

</div>

<style>
h1 {
  color: #2563eb !important;
  font-weight: 800;
  letter-spacing: -0.02em;
}
html.dark h1 {
  color: #4facfe !important;
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
    <span class="opacity-70">Primary MMP active-site contact. AB 1, AB 2, and AB 6 bind MMP9 more strongly than MMP2 on vendor-matched constructs, as does wild-type TIMP3 (one day, replication pending; see Results).</span>
  </div>
  <div class="p-3 bg-violet-500/10 rounded border border-violet-500/20 text-[10px]">
    <b class="text-violet-400 block mb-1">C Loop (res 63–68)</b>
    <span class="opacity-70">Second primary contact loop. C 12 and C 15 bind MMP9 more strongly than MMP2 on vendor-matched data, as does wild-type TIMP3 (one day, replication pending; see Results). Up to 13 aa insertions tolerated.</span>
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
Six stages from backbone hallucination to experimentally validated selectivity. This is the pipeline that generated and validated the primary result of this talk; <b class="text-amber-300">Round 2</b> (an iterative-refinement successor, now converged at iteration 92) and <b class="text-red-300">Round 3</b> (a specificity-aware successor, in progress) are covered later.

<div class="flex items-stretch gap-1 mt-8">
  <div class="flex-1 min-w-0 flex flex-col items-center p-2 bg-blue-500/10 rounded border border-blue-500/30 text-center">
    <div class="text-blue-400 font-black text-xl mb-1">①</div>
    <div class="text-[11px] font-bold text-white mb-2">RFdiffusion</div>
    <div class="text-[9px] opacity-50 leading-relaxed">Loop backbone hallucination<br/>AB / C / EF loops<br/>6–24 aa expansions<br/>A30 GPU, ~20 hr/run</div>
  </div>
  <div class="flex-none w-4 flex items-center justify-center text-blue-400 text-xl">→</div>
  <div class="flex-1 min-w-0 flex flex-col items-center p-2 bg-violet-500/10 rounded border border-violet-500/30 text-center">
    <div class="text-violet-400 font-black text-xl mb-1">②</div>
    <div class="text-[11px] font-bold text-white mb-2">ProteinMPNN</div>
    <div class="text-[9px] opacity-50 leading-relaxed">1000 seqs/target<br/>T=0.2, fixed scaffold<br/>Top 10 per loop<br/>by log-likelihood</div>
  </div>
  <div class="flex-none w-4 flex items-center justify-center text-blue-400 text-xl">→</div>
  <div class="flex-1 min-w-0 flex flex-col items-center p-2 bg-cyan-500/10 rounded border border-cyan-500/30 text-center">
    <div class="text-cyan-400 font-black text-xl mb-1">③</div>
    <div class="text-[11px] font-bold text-white mb-2">AlphaFold3</div>
    <div class="text-[9px] opacity-50 leading-relaxed">Batch co-folding<br/>5 targets<br/>ipTM / PAE / pLDDT<br/>loop & interface</div>
  </div>
  <div class="flex-none w-4 flex items-center justify-center text-blue-400 text-xl">→</div>
  <div class="flex-1 min-w-0 flex flex-col items-center p-2 bg-emerald-500/10 rounded border border-emerald-500/30 text-center">
    <div class="text-emerald-400 font-black text-xl mb-1">④</div>
    <div class="text-[11px] font-bold text-white mb-2">T-Score Filter</div>
    <div class="text-[9px] opacity-50 leading-relaxed">T > 2.0, ≥2 metrics<br/>ipTM ≥ 0.82<br/>13 candidates selected<br/>for synthesis</div>
  </div>
  <div class="flex-none w-4 flex items-center justify-center text-blue-400 text-xl">→</div>
  <div class="flex-1 min-w-0 flex flex-col items-center p-2 bg-amber-500/10 rounded border border-amber-500/30 text-center">
    <div class="text-amber-400 font-black text-xl mb-1">⑤</div>
    <div class="text-[11px] font-bold text-white mb-2">DNA Synthesis</div>
    <div class="text-[9px] opacity-50 leading-relaxed">13 vetted constructs ordered<br/>ORF checked clean of cut sites<br/>Twist Bioscience, Dec 2025<br/>ready-to-transform plasmids</div>
  </div>
  <div class="flex-none w-4 flex items-center justify-center text-blue-400 text-xl">→</div>
  <div class="flex-1 min-w-0 flex flex-col items-center p-2 bg-red-500/10 rounded border border-red-500/30 text-center">
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
      <p><b class="text-blue-300">Generation scale:</b> >1,000 designed sequences narrowed to 13 synthesis candidates. This is a shotgun-generation strategy — the large starting pool buys shots on goal, not a precision funnel; the filter's job is picking a synthesizable batch, not proving itself efficient.</p>
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

# T-Score: Ranking for Preference, Not Just Strength
Self-normalized: each variant is compared with its own cross-target mean.

<div class="mt-4">
  <TScoreFormula />
</div>

<div class="grid grid-cols-4 gap-3 mt-4">
  <div class="p-2 bg-emerald-500/10 rounded border border-emerald-500/20 text-center">
    <div class="text-emerald-400 font-black text-lg">549</div>
    <div class="text-[8px] uppercase tracking-widest text-emerald-400 mt-1">Raw AF3 co-folds</div>
  </div>
  <div class="p-2 bg-blue-500/10 rounded border border-blue-500/20 text-center">
    <div class="text-blue-400 font-black text-lg">128</div>
    <div class="text-[8px] uppercase tracking-widest text-blue-400 mt-1">Unique variants</div>
    <div class="text-[8px] opacity-50">126 with all four targets</div>
  </div>
  <div class="p-2 bg-violet-500/10 rounded border border-violet-500/20 text-center">
    <div class="text-violet-400 font-black text-lg">39</div>
    <div class="text-[8px] uppercase tracking-widest text-violet-400 mt-1">Shortlisted by rank</div>
    <div class="text-[8px] opacity-50">ApTM ≥ 0.80 gate</div>
  </div>
  <div class="p-2 bg-amber-500/10 rounded border border-amber-500/20 text-center">
    <div class="text-amber-400 font-black text-lg">15</div>
    <div class="text-[8px] uppercase tracking-widest text-amber-400 mt-1">Ordered</div>
    <div class="text-[8px] opacity-50">Manual curation</div>
  </div>
</div>

<div class="grid grid-cols-2 gap-4 mt-3">
  <div class="p-3 bg-white/4 rounded border border-white/10 text-[10px] leading-relaxed opacity-80">
    <b class="text-blue-400">T-score vs Z-score:</b> Z-score identifies variants that bind better <em>than the population</em>, so a globally excellent binder earns a high Z regardless of selectivity.
    T-score identifies variants that score higher on <em>this target than their own cross-target average</em>, which is intended to favor preference over global strength.
  </div>
  <div class="p-3 bg-emerald-500/5 rounded border border-emerald-500/20 text-[10px] leading-relaxed opacity-80">
    <b class="text-emerald-400">Selection as implemented:</b> rank-based consensus (top 10 per metric; Metric_Count ≥ 3 flags consensus), then the top 3 per category and target (highest affinity, lowest affinity, most target-specific, pairwise gap). No hard T &gt; 2.0 cutoff is applied in code, and 39 to 15 was a manual step. The T-score was not tested against the wet-lab results.
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
    <h4 class="text-emerald-300 font-bold text-[11px] uppercase tracking-widest">MMP9 Preference (5)</h4>
    <p class="text-[10px] leading-relaxed opacity-70">
      AB 1, AB 2, AB 6 (AB loop) and C 12, C 15 (C loop), designed for MMP9 preference. On raw Pos Med Ratio with vendor-matched data (Enzo only, n=2v2, separate cultures measured on one day), all five have MMP9 above MMP2 in every replicate (nominal Welch p=0.044, 0.008, 0.035, 0.016, 0.019; uncorrected), and so does TIMP3-WT.
    </p>
  </div>
  <div class="p-4 bg-blue-500/10 rounded border border-blue-500/20 space-y-2">
    <h4 class="text-blue-300 font-bold text-[11px] uppercase tracking-widest">MMP9 High / Low Controls (3)</h4>
    <p class="text-[10px] leading-relaxed opacity-70">
      AB 3 ("High", broad binder) and AB 5, C 13 ("Low", designed non-binders) were included to span the assay's dynamic range. AB 5 and C 13 lack replicated Enzo MMP2 trials, so they cannot be tested on the MMP9 axis vendor-matched.
    </p>
  </div>
  <div class="p-4 bg-amber-500/10 rounded border border-amber-500/20 space-y-2">
    <h4 class="text-amber-300 font-bold text-[11px] uppercase tracking-widest">ADAM17-Targeted (5)</h4>
    <p class="text-[10px] leading-relaxed opacity-70">
      AB 4, AB 7, C 11, C 14, ABC 22 designed for ADAM17 (and A17&gt;A10 where noted). ADAM10 comparison was limited by positive-control activity in this campaign, and the ADAM17-versus-ADAM10 axis remains untested.
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
  <div><span class="text-cyan-400 font-bold">Primary (cross-target):</span> Pos Med Ratio. Comparison only: Bind Med (Expr+), Binding Efficiency</div>
  <div><span class="text-violet-400 font-bold">Within-target:</span> Norm Bind Med, Norm Median Ratio, Norm IWB</div>
  <div><span class="text-slate-400 font-bold">QC only:</span> Stain Index</div>
</div>

<div class="mt-2 text-[8px] opacity-30 italic text-center">Hover any card above for its full definition, scope, and use.</div>

---
layout: default
transition: fade-out
---

# Primary Metric: Pos Med Ratio
Median APC over APC+ cells divided by median FITC over FITC+ cells: binding intensity per unit display, comparable across targets.

<div class="grid grid-cols-2 gap-8 mt-4 items-start">
  <div class="space-y-4">
    <div class="p-4 bg-white/5 rounded border border-white/10">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-3">Definition</h4>
      <div class="font-mono text-[11px] bg-black/30 rounded p-3 leading-relaxed">
        <span class="text-cyan-400">Pos Med Ratio</span> =<br/>
        &nbsp;&nbsp;median(APC-A)&nbsp;[APC+] &divide; median(FITC-A)&nbsp;[FITC+]
      </div>
      <p class="text-[10px] leading-relaxed opacity-60 mt-3">
        Retains how strongly cells bind and corrects for how much protein each construct
        displays. It is raw (not WT-normalized), so it compares directly across targets.
        The WT-normalized version, <b>Norm Median Ratio</b>, is for within-target ranking.
      </p>
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px]">
      <b class="text-amber-400">Comparison metrics:</b>
      Bind Med (Expr+) is the median APC in FITC+ cells (no division by display).
      Binding Efficiency is the fraction of expressors above the APC threshold, a count with
      no intensity, so it is shown only for comparison and is not an endpoint.
    </div>
  </div>
  <div class="space-y-2">
    <img src="/SharedAssets/figures/De_Novo_Binder_Generation/fig_pos_med_heatmap.png" class="w-full rounded bg-white" />
    <p class="text-[9px] opacity-40 italic text-center">
      Pos Med Ratio per construct and target. MMP9/MMP2 are Enzo-only for AB 1, AB 2, AB 6, C 12, C 15 and TIMP3-WT; other cells are pooled across suppliers.
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
      <b class="text-amber-400">Interpretation:</b> Values near 1.0 = binding similar to TIMP3-WT for that target. Selectivity signal comes from the raw cross-target comparison — this normalized view indicates variants are not simply super-binders on every target.
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Comparison Only: Binding Efficiency (DP/FITC+)
Count-fraction metric with no intensity information; shown to compare with Pos Med Ratio, not used as an endpoint.

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
      <h3 class="text-xs font-bold uppercase tracking-[0.2em] text-blue-400">Binding Efficiency on the same trials</h3>
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
            <td class="py-1 text-emerald-400 font-bold">AB 2</td>
            <td class="text-right py-1">22.2%</td>
            <td class="text-right py-1 text-emerald-400 font-bold">92.3%</td>
            <td class="text-right py-1 text-emerald-400">4.2×</td>
            <td class="text-right py-1 text-emerald-400">Enzo p=0.0004</td>
          </tr>
          <tr class="border-b border-white/5">
            <td class="py-1 text-emerald-400 font-bold">AB 6</td>
            <td class="text-right py-1">26.9%</td>
            <td class="text-right py-1 text-emerald-400 font-bold">94.1%</td>
            <td class="text-right py-1 text-emerald-400">3.5×</td>
            <td class="text-right py-1 text-emerald-400">Enzo p=0.0011</td>
          </tr>
          <tr class="border-b border-white/5 opacity-75">
            <td class="py-1">AB 1</td>
            <td class="text-right py-1">33.7%</td>
            <td class="text-right py-1">95.7%</td>
            <td class="text-right py-1">2.8×</td>
            <td class="text-right py-1">Enzo p=0.074, borderline</td>
          </tr>
          <tr class="border-b border-white/5 opacity-75">
            <td class="py-1">C 15</td>
            <td class="text-right py-1">39.2%</td>
            <td class="text-right py-1">92.7%</td>
            <td class="text-right py-1">2.4×</td>
            <td class="text-right py-1">Enzo p=0.068, borderline</td>
          </tr>
          <tr class="border-b border-white/5 opacity-60">
            <td class="py-1">C 12</td>
            <td class="text-right py-1">36.5%</td>
            <td class="text-right py-1">94.1%</td>
            <td class="text-right py-1">2.6×</td>
            <td class="text-right py-1">directional, n.s. even matched</td>
          </tr>
          <tr class="border-b border-white/10 opacity-50">
            <td class="py-1">C 13</td>
            <td class="text-right py-1">17.2%</td>
            <td class="text-right py-1">28.1%</td>
            <td class="text-right py-1">—</td>
            <td class="text-right py-1">n.s. (pooled, no Enzo n)</td>
          </tr>
          <tr class="border-b border-white/10 opacity-50">
            <td class="py-1">AB 5</td>
            <td class="text-right py-1">13.6%</td>
            <td class="text-right py-1">25.6%</td>
            <td class="text-right py-1">—</td>
            <td class="text-right py-1">n.s. (pooled, no Enzo n)</td>
          </tr>
          <tr>
            <td class="py-1 text-blue-400">TIMP 3</td>
            <td class="text-right py-1">32.8%</td>
            <td class="text-right py-1">91.9%</td>
            <td class="text-right py-1">2.8×</td>
            <td class="text-right py-1 text-blue-400">ref, Enzo p=0.022</td>
          </tr>
        </tbody>
      </table>
      <div class="text-[8px] opacity-40 italic border-t border-white/10 pt-2">
        Welch t-test, MMP9 vs MMP2, vendor-matched (Enzo only, n=2–3 per target) except C 13/AB 5 (pooled — no Enzo MMP2 replication exists for these two). Enzo is the one vendor supplying same-species, catalytic-domain-matched constructs for both targets; other vendors' constructs are documented mismatches (see next slide).
      </div>
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Norm Median Ratio — Cross-Target Comparison
APC/FITC ratio for FITC+ cells, normalized to TIMP3-WT. Click group buttons to switch target pairs.

<div class="mt-4 flex justify-center">
  <FcsChart />
</div>

<div class="mt-4 grid grid-cols-4 gap-3 text-[10px]">
  <div class="p-2 bg-white/5 rounded border border-white/10 opacity-80">
    <b class="text-cyan-400">Bind Med (Expr+)</b><br/>
    Raw APC-A intensity in expressors. Comparison intensity metric.
  </div>
  <div class="p-2 bg-white/5 rounded border border-white/10 opacity-80">
    <b class="text-emerald-400">Binding Efficiency</b><br/>
    Fraction of expressors that bind. Comparison only; the primary metric is raw Pos Med Ratio (Key Findings).
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

<div class="mt-3 p-2.5 bg-amber-500/10 rounded border border-amber-500/20 text-[9px] opacity-80">
  <b class="text-amber-400">Metric choice:</b> the WT-normalized median ratio shown here differs from the raw per-cell Pos Med Ratio used for the headline test on the Key Findings slide, and from the threshold-crossing Binding Efficiency on the previous slide; the three metrics are sensitive to different features of the binding distribution, so each design is read on all of them. Current full dataset (293 QC-passing trials through 2026-07-01).
</div>

---
layout: default
transition: fade-out
---

# Key Findings: Five Designs Bind MMP9 More Strongly Than MMP2, as Does Wild Type
Enzo only (same species, same catalytic domain for both targets), raw per-cell Pos Med Ratio. MMP9 exceeds MMP2 in every replicate for all five designs, and TIMP3-WT shows the same direction, so the data do not separate the designs from the scaffold. Two separate cultures per group (different wells, separately mixed label preps), measured on one day (2026-04-24); replication on a separate day is pending.

<div class="mt-3 flex justify-center">
  <img src="/SharedAssets/figures/De_Novo_Binder_Generation/fig_mmp9_vs_mmp2_matched.png" class="h-[250px] rounded bg-white" />
</div>

<div class="grid grid-cols-3 gap-4 mt-3 text-[9px]">
  <div class="p-2 bg-emerald-500/10 rounded border border-emerald-500/20">
    <b class="text-emerald-300">Nominal Welch p (Enzo only, n=2v2, uncorrected):</b> AB 1 p=0.044 (3.1×), AB 2 p=0.008 (2.7×), AB 6 p=0.035 (2.4×), C 12 p=0.016 (2.8×), C 15 p=0.019 (2.2×). With five comparisons only AB 2 clears a Bonferroni threshold of 0.01, and a rank-based test cannot return p&lt;0.33 at n=2v2.
  </div>
  <div class="p-2 bg-blue-500/10 rounded border border-blue-500/20">
    <b class="text-blue-300">TIMP3-WT shows the same direction:</b> 2.4× on the same day (p=0.16, n=2v2) and 3.2× across three replicates (p=0.141). The designs' folds are 0.93 to 1.28 times WT's (target-by-construct interaction p=0.32 to 1.00). The data show the scaffolds prefer MMP9 in this assay, not that redesign changed the preference.
  </div>
  <div class="p-2 bg-amber-500/10 rounded border border-amber-500/20">
    <b class="text-amber-400">Why not pooled:</b> Sino's "MMP2" is mouse (missing the N-catalytic half TIMP3 binds) and its "MMP9" is full-length, so the size of the difference depends on the supplier (1.1 to 2.0× Sino, 1.9 to 3.9× Enzo). Binding Efficiency, for comparison only, is significant for AB 2 and AB 6 because it discards binding intensity.
  </div>
</div>

---
layout: default
transition: fade-out
---

# Wild Type and the Other Constructs
The preference for MMP9 is present in every construct measured against both targets, at both suppliers.

<div class="grid grid-cols-2 gap-6 mt-4 items-start">
  <div>
    <img src="/SharedAssets/figures/De_Novo_Binder_Generation/fig_fold_by_supplier.png" class="w-full rounded bg-white" />
    <p class="text-[9px] opacity-40 italic text-center mt-1">MMP9:MMP2 fold in Pos Med Ratio. Triangles: one trial per target, run on a different day.</p>
  </div>
  <div class="space-y-3 text-[10px] leading-relaxed opacity-80">
    <div class="p-3 bg-white/5 rounded border border-white/10">
      <b class="text-blue-400">Designs vs WT (same day, Enzo):</b> fold difference as a multiple of TIMP3-WT's (2.41): AB 1 1.28, AB 2 1.12, AB 6 0.99, C 12 1.18, C 15 0.93. No design has more than 7% higher MMP9 signal than WT. On Bind Med (Expr+) the multiples are 0.94 to 1.36, and no interaction p is below 0.30.
    </div>
    <div class="p-3 bg-white/5 rounded border border-white/10">
      <b class="text-emerald-400">Other Enzo constructs:</b> seven constructs not designed for MMP9 preference (single trials) have folds of 1.9 to 3.9 (median 3.1), against 2.2 to 3.1 (median 2.7) for the five designs. This is descriptive, and shows the designs are not separated from the rest of the library.
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20">
      <b class="text-amber-400">What it means for the campaign:</b> the data are compatible with the campaign producing MMP9-preferring binders, but do not yet show that redesign is the source of the preference. Two replicates per group limit every comparison. The T-score used to rank candidates was not tested by these data.
    </div>
  </div>
</div>

---
layout: section
transition: slide-up
---

# Phase 3
### Second-Generation Design — Iterative Refinement & ESM-C (June–September 2026)

---
layout: default
transition: fade-out
---

# Round 2 Pipeline: Iterative Refinement
Round 1 validated the approach; Round 2 automates and accelerates it on a dedicated GPU cluster, looping backbone → sequence → structure and re-seeding each cycle from its own best designs.

<div class="flex items-stretch gap-1 mt-4">
  <div class="flex-1 min-w-0 flex flex-col items-center p-2 bg-blue-500/10 rounded border border-blue-500/30 text-center">
    <div class="text-blue-400 font-black text-lg mb-1">①</div>
    <div class="text-[10px] font-bold text-white mb-1.5">RFd3</div>
    <div class="text-[8px] opacity-50 leading-tight">RoseTTAFold Diffusion 3<br/>backbone generation</div>
  </div>
  <div class="flex-none w-4 flex items-center justify-center text-blue-400 text-lg">→</div>
  <div class="flex-1 min-w-0 flex flex-col items-center p-2 bg-violet-500/10 rounded border border-violet-500/30 text-center">
    <div class="text-violet-400 font-black text-lg mb-1">②</div>
    <div class="text-[10px] font-bold text-white mb-1.5">LigandMPNN</div>
    <div class="text-[8px] opacity-50 leading-tight">Annealed sequence design<br/>hot → cold temperature</div>
  </div>
  <div class="flex-none w-4 flex items-center justify-center text-blue-400 text-lg">→</div>
  <div class="flex-1 min-w-0 flex flex-col items-center p-2 bg-cyan-500/10 rounded border border-cyan-500/30 text-center">
    <div class="text-cyan-400 font-black text-lg mb-1">③</div>
    <div class="text-[10px] font-bold text-white mb-1.5">ESMFold2</div>
    <div class="text-[8px] opacity-50 leading-tight">MSA-free local ranker<br/>replaced Boltz-2</div>
  </div>
  <div class="flex-none w-4 flex items-center justify-center text-blue-400 text-lg">→</div>
  <div class="flex-1 min-w-0 flex flex-col items-center p-2 bg-emerald-500/10 rounded border border-emerald-500/30 text-center">
    <div class="text-emerald-400 font-black text-lg mb-1">④</div>
    <div class="text-[10px] font-bold text-white mb-1.5">Best Binders</div>
    <div class="text-[8px] opacity-50 leading-tight">Per-target Hall of Fame<br/>seeds the next cycle</div>
  </div>
  <div class="flex-none w-4 flex items-center justify-center text-blue-400 text-lg">→</div>
  <div class="flex-1 min-w-0 flex flex-col items-center p-2 bg-amber-500/10 rounded border border-amber-500/30 text-center">
    <div class="text-amber-400 font-black text-lg mb-1">⑤</div>
    <div class="text-[10px] font-bold text-white mb-1.5">AF3 Gold-Standard</div>
    <div class="text-[8px] opacity-50 leading-tight">Periodic validation<br/>of the cheap ranker</div>
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
  Compute: dedicated GPU cluster. Active targets: MMP2, MMP9, MMP3, MMP10, ADAM10, ADAM17. Selection script <code>select_binders_to_order.py</code> computes loop-pLDDT, loop-interface PAE, and framework RMSD vs. native TIMP3.
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
        <li><b>ADAM17:</b> Weak negative (r = -0.238, p = 0.456).</li>
        <li><b>MMP2:</b> Weak negative (r = -0.207, p = 0.519).</li>
        <li><b>MMP9:</b> Essentially zero (r = -0.022, p = 0.947).</li>
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

# Some Binders Are Just Sticky
Raw AlphaFold3 confidence looks useless for ranking binding; a stickiness confound explains part of why, but correcting for it doesn't yet recover a significant signal on the full current dataset.

<div class="grid grid-cols-2 gap-6 mt-4 items-start">
  <div class="space-y-3">
    <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px] leading-relaxed opacity-80">
      <b class="text-blue-400">Some constructs are just "sticky":</b> binding correlates strongly <em>across targets</em> for a given construct — a construct that binds one protease tends to bind them all, regardless of target. This general stickiness accounts for roughly two-thirds of the variance in binding scores; well under a third reflects genuine target-specific interaction, with a small remainder from baseline differences between targets themselves.
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px] leading-relaxed opacity-80">
      <b class="text-amber-400">Correcting for stickiness narrows the gap, but doesn't clear significance:</b> after removing each construct's stickiness, loop-pLDDT's correlation with the target-specific residual moves in the right direction (ρ=0.20, p=0.23, n=36) but isn't significant; ipTM stays flat (ρ=0.09, p=0.58). AB 2/AB 6's supported, vendor-matched selectivity result rests on its own wet-lab statistics, not on this AF-metric correlation — that finding doesn't depend on this correlation holding up.
    </div>
    <div class="p-3 bg-red-500/10 rounded border border-red-500/20 text-[10px] leading-relaxed opacity-80">
      <b class="text-red-400">ipTM tracks expression, not binding:</b> ipTM correlates with how well a construct <em>expresses</em> (ρ=0.37, p=0.027) but not with target-specific binding (ρ=0.09, p=0.58), and is nearly saturated across these designs (0.77–0.90) — too compressed to rank fine differences.
    </div>
  </div>
  <div class="space-y-2 flex flex-col items-center w-full">
    <CategoryBars dataset="variance" title="Variance Decomposition (%)" />
    <CategoryBars dataset="rho" title="Stickiness-Corrected Correlation (ρ)" />
    <p class="text-[9px] opacity-40 italic text-center">
      Top: stickiness / target baseline / genuine target-specific signal, recomputed against the current full dataset. Bottom: even after removing stickiness, loop-pLDDT/ipTM correlation with binding is still weak and not statistically significant (n=36); ipTM instead tracks expression.
    </p>
  </div>
</div>

---
layout: default
transition: fade-out
---

# The Ranking Metric's Evolution
A reasoned first pass, then a directly validated replacement.

<div class="space-y-2 mt-3">
  <div class="p-2 bg-white/5 rounded border border-white/10 text-[9px] leading-snug opacity-80">
    <b class="text-cyan-400">V1 — A reasoned recipe:</b> loop-pLDDT + interface PAE + an off-target selectivity term + ipTM as a small expressibility prior, weights fixed by reasoning (a data-fitted blend performed <em>worse</em> at this sample size). Out-ranked every individual metric, but its own uncertainty interval crossed zero — a working prior, not a precise model.
  </div>
  <div class="p-2 bg-emerald-500/10 rounded border border-emerald-500/20 text-[9px] leading-snug opacity-80">
    <b class="text-emerald-400">V2 — <code>sv_pdockq</code>, directly validated:</b> a Structural-Validation interface-battery pDockQ score, prospectively validated against AlphaFold3 across three independent samples with no decay (ρ=+0.655 → +0.591 → +0.586). Replaced the reasoned recipe as the active ranker.
  </div>
  <div class="p-1.5 bg-amber-500/10 rounded border border-amber-500/20 text-[8.5px] opacity-80">
    <b class="text-amber-400">Caveat:</b> raw <code>composite_score</code> isn't comparable across pipeline eras. <code>sv_pdockq</code> is the apples-to-apples number for "did it get better."
  </div>
</div>

<div class="grid grid-cols-2 gap-4 mt-2">
  <div class="space-y-1">
    <table class="w-full text-[9px] font-mono border-collapse">
      <thead>
        <tr class="border-b border-white/20 text-blue-400 font-bold">
          <th class="text-left py-1">Target</th>
          <th class="text-right py-1">Earlier best</th>
          <th class="text-right py-1">Converged</th>
        </tr>
      </thead>
      <tbody class="opacity-85">
        <tr class="border-b border-white/10"><td class="py-1">MMP2</td><td class="text-right">0.617</td><td class="text-right text-emerald-400 font-bold">0.673</td></tr>
        <tr class="border-b border-white/10"><td class="py-1">MMP9</td><td class="text-right">0.538</td><td class="text-right text-emerald-400 font-bold">0.612</td></tr>
        <tr class="border-b border-white/10"><td class="py-1">ADAM10</td><td class="text-right">0.634</td><td class="text-right text-emerald-400 font-bold">0.650</td></tr>
        <tr><td class="py-1">ADAM17</td><td class="text-right">0.643</td><td class="text-right text-emerald-400 font-bold">0.667</td></tr>
      </tbody>
    </table>
    <p class="text-[8px] opacity-40 italic text-center">Best <code>sv_pdockq</code> per target under V2, converged.</p>
  </div>
  <div class="p-2 bg-white/5 rounded border border-white/10 text-[8.5px] leading-snug opacity-70">
    <b class="text-blue-300">Convergence, not stall:</b> parent-to-child transmission of quality kept rising while frontier gain per round collapsed toward zero on every target — the signature of a search that found its optimum. A previously-suspected ceiling turned out to be a misconfigured noise parameter, not a scaffold limit.
  </div>
</div>

---
layout: default
transition: fade-out
---

# Switching the Design Base: HADDOCK → AF3/ESMFold2 Co-Folds
The appeal of HADDOCK was exact retention of each target's real crystal structure. In practice, designs built against HADDOCK-docked poses didn't bind well enough — so the design base moved to co-folded AF3/ESMFold2 structures instead.

<div class="grid grid-cols-2 gap-6 mt-4 items-start">
  <div class="space-y-3">
    <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px] leading-relaxed opacity-80">
      <b class="text-blue-400">Why HADDOCK didn't work as a design base:</b> HADDOCK docks rigid, unbound monomers — it cannot reconstruct the induced-fit interface this system actually uses. Designing against that pose meant designing against a binding mode that wasn't real, which is consistent with the binding shortfall observed when this was tried.
    </div>
    <div class="p-3 bg-emerald-500/10 rounded border border-emerald-500/20 text-[10px] leading-relaxed opacity-80">
      <b class="text-emerald-400">Mechanistic check:</b> AF3/ESMFold2 co-folds reproduce the native crystal's real zinc-chelation geometry; the HADDOCK-docked poses bury the zinc loop non-specifically, well off that geometry.
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px] leading-relaxed opacity-80">
      <b class="text-amber-400">Self-consistency check:</b> both co-fold methods produce consistently good-quality poses (mean DockQ ≈0.47–0.49); all four HADDOCK docking tracks (varying which structure — AF3, ESMFold2, or crystal — supplied each monomer) score far lower and about equally poorly (DockQ 0.07–0.08), regardless of which inputs were docked. A method that can't produce a good pose across any combination of equivalent inputs isn't a reliable design base, independent of any binding data.
    </div>
  </div>
  <div class="space-y-3 flex flex-col items-center">
    <CategoryBars dataset="haddock" title="Mean DockQ Across Sources" />
    <p class="text-[9px] opacity-40 italic text-center">
      Cross-source pose quality (DockQ, higher = better): both co-fold methods land around 0.47–0.49; all four HADDOCK docking tracks sit far lower (0.07–0.08), each roughly as poor as the others.
    </p>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Second-Generation Order
A calibration ladder, not a top-N list — nineteen constructs designed to test whether <code>sv_pdockq</code> tracks measured binding, not just to bank the highest scores.

<div class="grid grid-cols-2 gap-6 mt-4 items-start">
  <div class="space-y-3">
    <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px] leading-relaxed opacity-80">
      <b class="text-blue-400">Two arms:</b> 15 cysteine-free constructs spanning <code>sv_pdockq</code> across four calibration bands and all four targets (ADAM10, ADAM17, MMP2, MMP9) — deliberately not just the top-ranked designs, so the score↔binding relationship can be <i>estimated</i> rather than assumed. Plus 4 cysteine-bearing frontier designs as a separate test arm.
    </div>
    <div class="p-3 bg-emerald-500/10 rounded border border-emerald-500/20 text-[10px] leading-relaxed opacity-80">
      <b class="text-emerald-400">All 19 clear AF3 gates:</b> pTM 0.83–0.91, interface pTM 0.76–0.88.
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px] leading-relaxed opacity-80">
      <b class="text-amber-400">Cysteine-liability flag:</b> mean loop cysteines rise monotonically with score band; ~70% of scored designs carry ≥1 free cysteine. TIMP3 already carries 6 structural disulfides, so an unpaired thiol in a grafted loop is a plausible folding/expression liability — recurring loop-terminal motifs (<code>TICD</code>, <code>SRCD</code>, <code>RLCD</code>, <code>SICD</code>, <code>SICY</code>) recur across independent designs and targets.
    </div>
  </div>
  <div class="space-y-3">
    <div class="p-3 bg-white/5 rounded border border-white/10 text-[9px] opacity-70 leading-relaxed">
      <b class="text-blue-300">Why span the range instead of ordering winners:</b> AF3 interface confidence saturates above <code>sv_pdockq</code>≈0.5 and did not predict measured binding in the July calibration. Ordering only top scores would answer nothing about whether the metric transfers; the calibration ladder can.
    </div>
    <p class="text-[9px] opacity-40 italic text-center">Source: <code>Local/Pipeline_Calibration_2026-08/data/FINAL_ORDER_2026-08-27.csv</code>; notebook entry "Campaign Convergence and Construct Selection," 2026-09-01.</p>
  </div>
</div>

---
layout: section
transition: slide-up
---

# Phase 4
### Round 3 — Specificity-Aware Refinement (September 2026, In Progress)

---
layout: default
transition: fade-out
---

# Round 2 Never Selected for Specificity
So it's no surprise that nothing in the second-generation set came out selective — that's the gap Round 3 exists to close.

<div class="grid grid-cols-2 gap-6 mt-6 items-start">
  <div class="space-y-3">
    <div class="p-4 bg-white/5 rounded border border-white/10">
      <h4 class="text-blue-400 font-bold text-[10px] uppercase tracking-widest mb-3">Selectivity Check</h4>
      <p class="text-[11px] leading-relaxed opacity-75">
        Each second-generation candidate was AF3-folded twice — once against its own target, once against its paired protease — and the interface-confidence gap computed.
      </p>
      <div class="font-mono text-[11px] bg-black/30 rounded p-3 leading-relaxed mt-3">
        <span class="text-cyan-400">Own-target − paired-protease ipTM</span><br/>
        &nbsp;&nbsp;median: <span class="text-amber-300 font-bold">+0.030</span> (small, roughly even)
      </div>
    </div>
  </div>
  <div class="space-y-3">
    <div class="p-4 bg-amber-500/10 rounded-lg border border-amber-500/30">
      <h4 class="text-amber-400 font-bold text-[10px] uppercase tracking-widest mb-2">Expected, Not a Surprise</h4>
      <p class="text-[10px] leading-relaxed opacity-80">
        No stage of the Round 2 design loop ever evaluated an off-target, so there was never a mechanism for selectivity to emerge — this result confirms that rather than revealing a new problem. Consistent with the bench-side finding on the first-generation library (design intent largely didn't translate to measured selectivity either).
      </p>
    </div>
    <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px] leading-relaxed opacity-80">
      Selection can observe selectivity after the fact; it can't manufacture it. Optimizing for it directly — which Round 3 does — is the fix.
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Round 3: Optimizing Selectivity Directly
Same engine, new objective — every design is folded against both its on-target and a paired off-target, and ranked on selectivity as well as affinity.

<div class="grid grid-cols-2 gap-6 mt-4 items-start">
  <div class="space-y-3">
    <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px] leading-relaxed opacity-80">
      <b class="text-blue-400">Mechanism:</b> subclasses Round 2's exact backbone → sequence → structure loop. Design still happens against the on-target only; the off-target is used purely for cross-scoring. Ranked on a composite of on-target quality and a selectivity term (normalized on/off interface contact-density gap), 50/50 weighted.
    </div>
    <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px] leading-relaxed opacity-80">
      <b class="text-cyan-400">Two pairs, code-enforced:</b> <b>MMP9-over-MMP2</b> and <b>ADAM17-over-ADAM10</b> — chosen to directly test the two axes this deck already flagged as unresolved (MMP2's plateau, ADAM17/ADAM10's weak data).
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px] leading-relaxed opacity-80">
      <b class="text-amber-400">Status:</b> an early run had the on/off-target pairing backwards — discovered, corrected, and that run discarded. Since the fix: an interim run refined the structural-interface scoring, and a live run — now with beam-search conformation refinement added — is still going as this deck is written.
    </div>
  </div>
  <div class="space-y-2">
    <table class="w-full text-[10px] font-mono border-collapse">
      <thead>
        <tr class="border-b border-white/20 text-blue-400 font-bold">
          <th class="text-left py-1">Pair</th>
          <th class="text-right py-1">Best selectivity so far</th>
          <th class="text-right py-1">Trend</th>
        </tr>
      </thead>
      <tbody class="opacity-85">
        <tr class="border-b border-white/10"><td class="py-1">MMP9&gt;MMP2</td><td class="text-right text-blue-300">0.573</td><td class="text-right">early peak, since noisy</td></tr>
        <tr><td class="py-1">ADAM17&gt;ADAM10</td><td class="text-right text-blue-300">0.615</td><td class="text-right">early peak, since noisy</td></tr>
      </tbody>
    </table>
    <p class="text-[9px] opacity-40 italic text-center mt-1">Selectivity score = on-target minus off-target ESMFold2 ipTM. A mid-flight snapshot — these numbers will move.</p>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Early Read: MMP9 Selectivity Looks Tractable; ADAM17 Does Not
Population-level signal from the most recent complete-enough window — hedged, and still early.

<div class="grid grid-cols-2 gap-6 mt-6 items-start">
  <div class="p-4 bg-emerald-500/10 rounded-lg border border-emerald-500/30 space-y-2">
    <h4 class="text-emerald-400 font-bold text-[11px] uppercase tracking-widest">MMP9-over-MMP2: Achievable</h4>
    <p class="text-[10px] leading-relaxed opacity-80">
      54% of scored designs prefer MMP9 over MMP2 (<code>selectivity_cd</code> &gt; 0). Of designs called selective, 79% are structurally corroborated by an independent ESMFold2 ipTM check. A clear majority-direction signal, even before the pipeline is optimizing for it directly.
    </p>
  </div>
  <div class="p-4 bg-red-500/10 rounded-lg border border-red-500/30 space-y-2">
    <h4 class="text-red-400 font-bold text-[11px] uppercase tracking-widest">ADAM17-over-ADAM10: Essentially Fails</h4>
    <p class="text-[10px] leading-relaxed opacity-80">
      Only 22% of scored designs prefer ADAM17 over ADAM10 — most grab the sticky ADAM10 off-target harder than the intended on-target. Structural corroboration is also weaker (55% vs. 79% for MMP9) — the least credible axis in this campaign.
    </p>
  </div>
</div>

<div class="mt-4 grid grid-cols-2 gap-4">
  <div class="p-2.5 bg-white/5 rounded border border-white/10 text-[9px] opacity-70 leading-snug">
    <b class="text-blue-400">Relation to Round 1:</b> the first-generation constructs designed for MMP9 preference all show MMP9 above MMP2 on vendor-matched data, but so does wild-type TIMP3, so that result does not yet show a redesign effect. The ADAM17&gt;ADAM10 claim (AB 4) rested on weak ADAM10 data, and the ADAM17/ADAM10 axis is the harder one in both generations.
  </div>
  <div class="p-2.5 bg-amber-500/10 rounded border border-amber-500/20 text-[9px] opacity-80">
    <b class="text-amber-400">Not yet actionable:</b> no ordering/synthesis output exists from this pipeline yet — it's still early, and these percentages will be re-measured before any candidate is proposed.
  </div>
</div>

<div class="mt-2 text-[8px] opacity-30 italic text-center">Source: <code>Local/specificity_refinement_20260901/</code>, <code>Local/specificity_refinement/</code>, and the 2026-09-02 daily brief.</div>

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
        <span><b>Replicate set and activity assays (top priority):</b> AB 1, AB 2, AB 6, C 12, C 15 and TIMP3-WT on every day, more than two replicates per group, to compare each design with wild type directly. Soluble expression of the variants is not yet available, so a whole-cell inhibition assay against WT-displaying cells at matched expression is the nearer-term proxy for activity; SPR (K<sub>D</sub>, k<sub>on</sub>, k<sub>off</sub>) follows once soluble protein exists.</span>
      </li>
      <li class="flex items-start gap-2">
        <div class="w-1.5 h-1.5 rounded-full bg-amber-400 mt-1 shrink-0"></div>
        <span><b>Second-generation order — awaiting decision:</b> 19-construct calibration ladder (AF3-gated, cysteine-controlled, spans the <code>sv_pdockq</code> score range) is ready but not yet placed with Twist. Go/no-go is the immediate next step, not more compute.</span>
      </li>
      <li class="flex items-start gap-2">
        <div class="w-1.5 h-1.5 rounded-full bg-cyan-400 mt-1 shrink-0"></div>
        <span><b>Reverse/orthogonal selectivity — now in progress, not stalled:</b> Round 3's specificity-aware refinement directly optimizes MMP9&gt;MMP2 and ADAM17&gt;ADAM10. Early read: MMP9&gt;MMP2 looks tractable (54% of designs, 79% structurally corroborated); ADAM17&gt;ADAM10 is hard (22%, 55% corroborated) — still early, pre-decision.</span>
      </li>
      <li class="flex items-start gap-2">
        <div class="w-1.5 h-1.5 rounded-full bg-emerald-400 mt-1 shrink-0"></div>
        <span><b>ADAM10/17 — still unresolved:</b> In-house ADAM10 <i>production</i> is validated (0.20 mg/mL, Western+BCA), but its <i>binding</i> is not — the earlier "positive" read was retracted by a July 24 control test and the vendor manifest. A first-ever ADAM10 ELISA pilot (Aug 30) is encouraging on the FLAG tag (~393× over background) but confirms tag presence, not catalytic identity. ADAM17 construct separately came back an empty vector — needs re-cloning (details: companion <i>ADAM Target Production</i> deck).</span>
      </li>
      <li class="flex items-start gap-2">
        <div class="w-1.5 h-1.5 rounded-full bg-cyan-400 mt-1 shrink-0"></div>
        <span><b>ESM-C — exploratory only (appendix):</b> Consensus motifs found could nominate a future batch, but <code>sv_pdockq</code> is now the validated, active selection method.</span>
      </li>
    </ul>
  </div>

  <div class="space-y-2">
    <h3 class="text-violet-400 font-bold uppercase text-[11px] tracking-widest border-b border-violet-500/20 pb-1">Translational Context</h3>
    <div class="p-3 bg-violet-500/5 rounded border border-violet-500/20 text-[9px] leading-snug opacity-80 space-y-2">
      <p><b class="text-violet-300">MMP9 as oncology target:</b> Overexpressed in many tumor types. Prior pan-MMP inhibitors failed in Phase III, and off-target inhibition with dose-limiting toxicity is one of the factors discussed. Variants that bind MMP9 in preference to MMP2 would address one of those factors; inhibition has not yet been tested.</p>
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
### ESM-C Sequence Classifier — Updated 2026-09-05: MMP9 Head Now Adopted

---
layout: default
transition: fade-out
---

# ESM-C Sequence Classifier
<b class="text-amber-300">Updated 2026-09-05</b> — the single-classifier numbers below are superseded by a 5-variant comparison; see paper §"Sequence-Based Ranking" for full detail. This slide now reports the current, traceable `all3_original` enumeration.

<div class="grid grid-cols-2 gap-6 mt-6 items-start">
  <div class="space-y-3">
    <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px] leading-relaxed opacity-80">
      <b class="text-blue-400">Scale:</b> A fine-tuned <code>Synthyra/ESMplusplus_large</code> classifier enumerated all ~64M possible 6-residue C-loop sequences; retained the top 50,000 predicted binders per target (P≥0.99) for ADAM17, MMP3, MMP9. MMP9/MMP3 share 36.8% Jaccard overlap; MMP9/ADAM17 are the most distinct pair (6.2%).
    </div>
    <div class="p-3 bg-white/5 rounded border border-white/10 text-[10px] leading-relaxed opacity-80">
      <b class="text-cyan-400">Consensus motifs found:</b> ADAM17 <code>LPSDTT</code>, MMP3 <code>LSPDTT</code>, MMP9 <code>LSPTTL</code>. Positions 1 (L) and 5 (T) are shared by all three; ADAM17 diverges at positions 2-3, MMP9 diverges at positions 4/6 — candidate leads for a future order, not yet acted on.
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

# ESM-C: Now the Adopted MMP9 Ranker
Held-out classifier performance across 5 training-data variants: MMP9 is strongest in the two pooled slices (MCC 0.72) and varies widely across slices (-0.01 to 0.72); ADAM17 modest; MMP3 a base-rate mirage. Deployment: rank by the MMP9/`everything_combined` head.

<div class="grid grid-cols-2 gap-6 mt-6 items-start">
  <div class="space-y-3">
    <div class="p-3 bg-emerald-500/10 rounded border border-emerald-500/20 text-[10px] leading-relaxed opacity-80">
      <b class="text-emerald-400">Best MMP9 variant (`everything_combined`):</b> ROC-AUC 0.955, PR-AUC 0.845, MCC 0.721, F1 0.753 (n=8,059). `mmp9_other`, trained on a much larger, more realistic MMP9-only dataset, is close behind (MCC 0.716). Both far outperform the small-data variants (MCC 0.43-0.45) and `cloop_only` (MCC ≈0, chance).
    </div>
    <div class="p-3 bg-amber-500/10 rounded border border-amber-500/20 text-[10px] leading-relaxed opacity-80">
    <b class="text-amber-400">Prior wet-lab correlation claim (ρ=0.86, n=7): withdrawn 2026-09-23.</b> It was the largest of 30 tests, the ADAM17 probabilities spanned only 0.506 to 0.520, the predictions came from the small model, and on the current data it is ρ=0.68 (p=0.094). Agreement with flow-cytometry outcomes has not been tested with the current models.
    </div>
    <div class="p-3 bg-blue-500/10 rounded border border-blue-500/20 text-[10px] opacity-80">
      <b class="text-blue-400">Open question, not yet resolved:</b> `mmp9_other`'s completed C-loop enumeration (2026-09-05) lands on a different top motif (<code>PYGSLC</code>) than `all3_original`'s (<code>LSPTTL</code>) for the same target — both confirmed current large-model runs. Not picked one way or the other pending investigation.
    </div>
  </div>
  <div class="flex flex-col items-center">
    <CategoryBars dataset="esmc" title="Held-Out Test Performance (MMP9, everything_combined)" />
    <p class="text-[9px] opacity-40 italic text-center mt-1">Held-out ROC/PR performance, current best MMP9 variant.</p>
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
    <div class="text-emerald-400 font-black text-2xl">5 of 5</div>
    <div class="text-[9px] opacity-40 uppercase tracking-widest mt-1">Designs with MMP9 &gt; MMP2 in every replicate (Enzo, n=2v2); TIMP3-WT shows the same direction</div>
  </div>
  <div>
    <div class="text-blue-400 font-black text-2xl">2.4×</div>
    <div class="text-[9px] opacity-40 uppercase tracking-widest mt-1">TIMP3-WT fold on the same day (designs: 2.2 to 3.1×)</div>
  </div>
  <div>
    <div class="text-violet-400 font-black text-2xl">3</div>
    <div class="text-[9px] opacity-40 uppercase tracking-widest mt-1">Design Generations (Rounds 1–3)</div>
  </div>
</div>

<div class="mt-8 text-xs opacity-30 uppercase tracking-widest">
  Repository: PhD-Research / Demonstrations / De_Novo_Binder_Generation
</div>
