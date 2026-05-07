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
    <h3 class="text-xs font-bold mb-4 uppercase tracking-widest opacity-60">Top Design Variants</h3>
    <TopVariants />
  </div>
  <div class="space-y-6">
    <div class="p-4 bg-white/5 rounded border border-white/10">
      <h3 class="text-xs font-bold mb-3 uppercase tracking-widest opacity-40">Z-Score Calibration</h3>
      <div class="text-[11px] leading-relaxed opacity-70">
        AlphaFold metrics (ipTM, pLDDT, PAE) were standardized into Z-scores relative to an ensemble of 500 random design iterations.
      </div>
    </div>
    <div class="p-4 bg-white/5 rounded border border-white/10">
      <h3 class="text-xs font-bold mb-2 uppercase tracking-widest opacity-40">Filtering Thresholds</h3>
      <ul class="text-[11px] space-y-1 opacity-70 list-none p-0">
        <li>• ipTM Z > 1.0 (High confidence interface)</li>
        <li>• Loop pLDDT Z > 0.5 (Structural stability)</li>
        <li>• Interface PAE Z < -0.5 (Distance consistency)</li>
      </ul>
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# AlphaFold Confidence Matrix
Standardized Z-Score analysis across target pairs.

<div class="mt-8 flex justify-center">
  <ZScoreTable />
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
    <div class="text-sm leading-relaxed opacity-80">
      Binders were evaluated for normalization against TIMP3-WT across multi-target trials. 
      The <span class="text-blue-400 font-bold">AB loops</span> showed a significant enrichment in binding efficiency (DP/FITC+) for the <b>MMP3 and MMP9 families</b>, while ADAM10 served as a negative control due to batch-specific target degradation.
    </div>
    <div class="p-4 bg-blue-500/10 rounded-lg border border-blue-500/20 text-xs leading-relaxed">
      • Target specificity switch: Successful discrimination between MMP and ADAM families<br/>
      • Correlation: AlphaFold ipTM Z-scores vs Experimental affinity
    </div>
  </div>
  <div class="p-6 bg-white/5 rounded border border-white/10 space-y-4 flex flex-col justify-center h-full">
    <h3 class="text-xs font-bold uppercase tracking-[0.2em] opacity-40">Experimental Results</h3>
    <div class="flex gap-8">
      <div>
        <div class="text-3xl font-black text-white">~1.7x</div>
        <div class="text-[9px] uppercase tracking-widest text-blue-400 font-bold">MMP3 Enrichment</div>
      </div>
      <div>
        <div class="text-3xl font-black text-white">93%</div>
        <div class="text-[9px] uppercase tracking-widest text-blue-400 font-bold">Binding Efficiency</div>
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
  <div class="space-y-4">
    <h3 class="text-blue-400 font-bold uppercase text-xs tracking-widest">Key Findings</h3>
    <ul class="text-[13px] space-y-3 opacity-80 list-none p-0">
      <li class="flex gap-3">
        <span class="text-blue-500 font-bold">01</span>
        <span>Successful de novo generation of high-affinity binders for ADAM10.</span>
      </li>
      <li class="flex gap-3">
        <span class="text-blue-500 font-bold">02</span>
        <span>Validated selectivity switch between ADAM10/17 and MMP2/9 families.</span>
      </li>
      <li class="flex gap-3">
        <span class="text-blue-500 font-bold">03</span>
        <span>Correlation established between AlphaFold Z-scores and binding efficiency.</span>
      </li>
    </ul>
  </div>
  <div class="grid grid-cols-2 gap-4">
    <div class="p-4 bg-white/5 rounded border border-white/10 flex flex-col justify-center items-center text-center">
      <div class="text-xs font-bold uppercase tracking-widest mb-1">In Vitro</div>
      <div class="text-[10px] opacity-50">FRET Kinetics</div>
    </div>
    <div class="p-4 bg-white/5 rounded border border-white/10 flex flex-col justify-center items-center text-center">
      <div class="text-xs font-bold uppercase tracking-widest mb-1">Structural</div>
      <div class="text-[10px] opacity-50">X-Ray Resolution</div>
    </div>
    <div class="p-4 bg-white/5 rounded border border-white/10 flex flex-col justify-center items-center text-center">
      <div class="text-xs font-bold uppercase tracking-widest mb-1">Optimization</div>
      <div class="text-[10px] opacity-50">ESM-2 Refinement</div>
    </div>
    <div class="p-4 bg-white/5 rounded border border-white/10 flex flex-col justify-center items-center text-center">
      <div class="text-xs font-bold uppercase tracking-widest mb-1">Expansion</div>
      <div class="text-[10px] opacity-50">Custom Markers</div>
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
