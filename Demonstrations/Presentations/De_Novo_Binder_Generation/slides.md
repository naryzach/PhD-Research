---
layout: center
class: text-center
transition: fade-out
---

# De Novo Binder Generation
### RFdiffusion to Experimental Validation Pipeline

<div class="mt-10 flex justify-center gap-4">
  <div class="px-4 py-2 bg-white/5 rounded border border-white/10">
    <div class="text-[10px] uppercase opacity-40 tracking-widest">Thread Status</div>
    <div class="text-blue-400 font-bold uppercase text-xs">Validated Results</div>
  </div>
  <div class="px-4 py-2 bg-white/5 rounded border border-white/10">
    <div class="text-[10px] uppercase opacity-40 tracking-widest">Primary Targets</div>
    <div class="text-blue-400 font-bold uppercase text-xs">ADAM10/17 | MMP2/9</div>
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

# Workflow Architecture
The integration of generative modeling and high-throughput experimental screening.

<div class="mt-8 scale-90 origin-top">
  <PipelineFlow />
</div>

<div class="grid grid-cols-2 gap-8 mt-12">
  <div class="p-4 bg-white/5 rounded border border-white/10">
    <h3 class="text-blue-400 font-bold mb-2 uppercase text-xs tracking-widest">Computational Phase</h3>
    <ul class="text-[13px] space-y-2 opacity-80 list-none p-0">
      <li class="flex gap-2"><span>—</span> <span>Loop-centric RFdiffusion on catalytic domains</span></li>
      <li class="flex gap-2"><span>—</span> <span>Sequence optimization via ProteinMPNN</span></li>
      <li class="flex gap-2"><span>—</span> <span>Z-score filtering of AlphaFold ensemble predictions</span></li>
    </ul>
  </div>
  <div class="p-4 bg-white/5 rounded border border-white/10">
    <h3 class="text-cyan-400 font-bold mb-2 uppercase text-xs tracking-widest">Experimental Phase</h3>
    <ul class="text-[13px] space-y-2 opacity-80 list-none p-0">
      <li class="flex gap-2"><span>—</span> <span>Golden Gate assembly of design libraries</span></li>
      <li class="flex gap-2"><span>—</span> <span>Multiplexed Yeast Surface Display titration</span></li>
      <li class="flex gap-2"><span>—</span> <span>Automated FCS affinity and specificity profiling</span></li>
    </ul>
  </div>
</div>

---
layout: section
transition: slide-up
---

# Phase 1
### Computational Design and Validation

---
layout: default
transition: fade-out
---

# RFdiffusion Loop Scaffolding
Targeting variable surface loops on the TIMP3 scaffold.

<div class="grid grid-cols-2 gap-8 mt-8 items-center">
  <LoopConfig />
  <div class="space-y-6">
    <div class="text-sm leading-relaxed opacity-80">
      We utilized length-variable loop generation (8–24 residues) to optimize the interaction surface against ADAM10 and ADAM17 catalytic domains.
    </div>
    <div class="p-5 bg-blue-500/5 rounded-lg border border-blue-500/20">
      <h4 class="text-[10px] font-bold text-blue-400 mb-3 uppercase tracking-widest">Design Parameters</h4>
      <div class="grid grid-cols-2 gap-y-3 gap-x-4 text-[11px] font-mono">
        <span class="opacity-40 uppercase">Targets</span> <span class="text-blue-200">MMP / ADAM</span>
        <span class="opacity-40 uppercase">Lengths</span> <span class="text-blue-200">8, 12, 16, 20, 24</span>
        <span class="opacity-40 uppercase">Noise</span> <span class="text-blue-200">0.1 - 0.5</span>
        <span class="opacity-40 uppercase">Iterations</span> <span class="text-blue-200">50</span>
      </div>
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Sequence Optimization & Metrics
Ranking design confidence via ProteinMPNN and AlphaFold Z-Scores.

<div class="grid grid-cols-2 gap-8 mt-4">
  <div>
    <h3 class="text-xs font-bold mb-4 uppercase tracking-widest text-blue-400">Final Twist Ordered Library <span class="opacity-40">(Dec 2025)</span></h3>
    <TopVariants />
  </div>
  <div class="space-y-4">
    <div class="p-4 bg-white/5 rounded border border-white/10">
      <h4 class="text-[10px] font-bold uppercase tracking-widest mb-2 text-blue-300">Consensus Selection Logic</h4>
      <p class="text-[10px] leading-relaxed opacity-70 mb-2">
        To minimize false positives, we employed a <b>multi-dimensional consensus filter</b>. Candidates were ranked across independent metrics for interface confidence (iPTM), local loop stability (LpLDDT), and distance consistency (PAE).
      </p>
      <div class="text-[9px] opacity-40 italic border-t border-white/5 pt-2">
        • <b>Elite:</b> Candidates with exceptional interface confidence (iPTM > 0.86)<br/>
        • <b>Design:</b> Robust winners passing multi-metric consensus thresholds
      </div>
    </div>
    <div class="p-4 bg-white/5 rounded border border-white/10">
      <h4 class="text-[10px] font-bold uppercase tracking-widest mb-2 text-blue-300">T-Score Calibration</h4>
      <p class="text-[10px] leading-relaxed opacity-70">
        Structural scores are transformed into T-scores by normalizing each variant's performance against its own mean across all targets. This identifies <b>statistically significant outliers</b>, ensuring that selected binders represent true structural "wins" relative to the design's baseline.
      </p>
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# AlphaFold Statistical Significance
Variant-specific selectivity analysis using T-scores across targets.

<div class="mt-8 flex justify-center">
  <TScoreTable />
</div>

---
layout: default
transition: fade-out
---

# Methodology: T-Score Calibration
How we identify statistically significant binders.

<div class="mt-2 flex flex-col items-center text-center max-w-3xl mx-auto space-y-4">
  <TScoreFormula />

  <div class="grid grid-cols-2 gap-6 w-full">
    <div class="bg-emerald-500/10 p-4 rounded-lg border border-emerald-500/20">
      <div class="text-emerald-400 font-black text-xl mb-1">T > 2.0</div>
      <div class="text-[10px] uppercase tracking-widest opacity-60">Significant Win</div>
    </div>
    <div class="bg-rose-500/10 p-4 rounded-lg border border-rose-500/20">
      <div class="text-rose-400 font-black text-xl mb-1">T < -2.0</div>
      <div class="text-[10px] uppercase tracking-widest opacity-60">Significant Loss</div>
    </div>
  </div>

  <p class="text-xs opacity-40 italic">
    *T-scores effectively map the "winning" target for every candidate in the library by normalizing out the intrinsic scaffold confidence levels.
  </p>
</div>

---
layout: default
---

# Wild Type Structural Folds
Interactive AlphaFold models of the WT TIMP3-Target complexes.

<div class="mt-4 h-[400px]">
  <StructuralViewer />
</div>

---
layout: section
transition: slide-up
---

# Phase 2
### Experimental Validation and Selectivity

---
layout: default
transition: fade-out
---

# Interactive FCS Analysis
Detailed exploration of raw population distributions and gating metrics.

<div class="flex justify-center w-full">
  <FcsViewer />
</div>

---
layout: default
transition: fade-out
---

# Binding Efficiency and Specificity
Normalizing experimental results against TIMP3-WT baseline.

<div class="grid grid-cols-2 gap-8 mt-12 items-center">
  <div class="p-6 bg-white/5 rounded border border-white/10 space-y-6">
    <h3 class="text-sm font-bold uppercase tracking-widest text-blue-400">Analysis Methodology</h3>
    <div class="text-[11px] leading-relaxed opacity-80">
      Experimental results were normalized against <b>TIMP3-WT</b> to identify variants that outcompete the natural binder. 
      Using a <b>replicate-aware T-score filter</b>, we decoupled raw signal intensity from true binding efficiency. 
      This allows us to distinguish between <span class="text-blue-400 font-bold">promiscuous binders</span> (high signal across all targets) and <span class="text-emerald-400 font-bold">specific winners</span> (statistically significant outliers for a single target).
    </div>
    <div class="p-4 bg-blue-500/10 rounded-lg border border-blue-500/20 text-[10px] leading-relaxed">
      • <b>Normalization:</b> Median Intensity Ratio vs TIMP3-WT baseline<br/>
      • <b>Quality Control:</b> Filtering by Double+ % to ensure population stability<br/>
      • <b>Correlation:</b> AlphaFold ipTM T-scores vs Experimental Selectivity
    </div>
  </div>
  <div class="p-6 bg-white/5 rounded border border-white/10 space-y-4 flex flex-col justify-center h-full">
    <h3 class="text-xs font-bold uppercase tracking-[0.2em] opacity-40">Design Wins vs TIMP3-WT</h3>
    <div class="flex flex-col gap-4">
      <div class="flex items-center justify-between">
        <div>
          <div class="text-2xl font-black text-white">1.52x</div>
          <div class="text-[8px] uppercase tracking-widest text-blue-400 font-bold">AB 4 @ MMP3</div>
        </div>
        <div class="text-[10px] opacity-60 text-right"><b>Better than WT</b><br/>Enrichment Ratio</div>
      </div>
      <div class="flex items-center justify-between">
        <div>
          <div class="text-2xl font-black text-white">91.2%</div>
          <div class="text-[8px] uppercase tracking-widest text-emerald-400 font-bold">C 12 @ MMP9</div>
        </div>
        <div class="text-[10px] opacity-60 text-right"><b>Binding Efficiency</b><br/>(DP/FITC+)</div>
      </div>
      <div class="border-t border-white/5 pt-2">
        <div class="text-[9px] leading-tight opacity-50">
          <b>Key Specificity Switch:</b> Variant <code>C 12</code> outcompetes WT on <b>MMP3</b> (1.42x) but shows near-zero enrichment on <b>MMP2</b> (0.92x), achieving targeted discrimination between the twin targets.
        </div>
      </div>
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Cross-Target Binding Analysis
Statistical summary of normalized median intensity across replicates.

<div class="mt-2 flex justify-center">
  <FcsChart />
</div>

---
layout: default
transition: fade-out
---

# Conclusion and Strategic Roadmap
Summary of current progress and future research directions.

<div class="grid grid-cols-2 gap-12 mt-12">
  <div class="space-y-6">
    <h3 class="text-blue-400 font-bold uppercase text-xs tracking-widest">Key Findings</h3>
    <ul class="text-[12px] space-y-4 opacity-80 list-none p-0">
      <li class="flex gap-3">
        <div class="w-6 h-6 rounded-full bg-blue-500/20 flex items-center justify-center text-[10px] font-bold text-blue-400 shrink-0">01</div>
        <span><b>Superior Binding:</b> Identified de novo loops that outcompete TIMP3-WT (1.52x enrichment at MMP3).</span>
      </li>
      <li class="flex gap-3">
        <div class="w-6 h-6 rounded-full bg-blue-500/20 flex items-center justify-center text-[10px] font-bold text-blue-400 shrink-0">02</div>
        <span><b>Target Discrimination:</b> Validated selectivity switch between twin targets (MMP3 vs MMP2).</span>
      </li>
      <li class="flex gap-3">
        <div class="w-6 h-6 rounded-full bg-blue-500/20 flex items-center justify-center text-[10px] font-bold text-blue-400 shrink-0">03</div>
        <span><b>Predictive Success:</b> Established AlphaFold <b>T-scores</b> as a robust proxy for experimental selectivity.</span>
      </li>
    </ul>
  </div>
  <div class="grid grid-cols-2 gap-4">
    <div class="p-4 bg-white/5 rounded border border-white/10 flex flex-col justify-center items-center text-center hover:bg-blue-500/10 transition-colors">
      <div class="text-xs font-bold uppercase tracking-widest mb-1 text-blue-400">Kinetics</div>
      <div class="text-[9px] opacity-50">FRET Binding Assays</div>
    </div>
    <div class="p-4 bg-white/5 rounded border border-white/10 flex flex-col justify-center items-center text-center hover:bg-blue-500/10 transition-colors">
      <div class="text-xs font-bold uppercase tracking-widest mb-1 text-blue-400">ML Refining</div>
      <div class="text-[9px] opacity-50">ESM-2 Loop Optimization</div>
    </div>
    <div class="p-4 bg-white/5 rounded border border-white/10 flex flex-col justify-center items-center text-center hover:bg-blue-500/10 transition-colors">
      <div class="text-xs font-bold uppercase tracking-widest mb-1 text-blue-400">Structural</div>
      <div class="text-[9px] opacity-50">CD & Crystallography</div>
    </div>
    <div class="p-4 bg-white/5 rounded border border-white/10 flex flex-col justify-center items-center text-center hover:bg-blue-500/10 transition-colors">
      <div class="text-xs font-bold uppercase tracking-widest mb-1 text-blue-400">Expansion</div>
      <div class="text-[9px] opacity-50">TIMP-1/2 Cross-Targets</div>
    </div>
  </div>
</div>

---
layout: center
class: text-center
---

# Questions
### De Novo Binder Generation Research Thread

<div class="mt-8 text-xs opacity-40 uppercase tracking-widest">
  Repository: PhD-Research / Demonstrations / De_Novo_Binder_Generation
</div>
