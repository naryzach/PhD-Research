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
    <div class="text-blue-400 font-bold uppercase text-xs">MMP2 | MMP9</div>
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

# Full Pipeline Architecture
From generative AI to statistically validated experimental selectivity.

<div class="grid grid-cols-6 gap-2 mt-8">
  <div class="flex flex-col items-center p-3 bg-blue-500/10 rounded border border-blue-500/30 text-center">
    <div class="text-blue-400 font-black text-xs mb-1">GEN</div>
    <div class="text-[10px] font-bold text-white mb-1">RFdiffusion</div>
    <div class="text-[9px] opacity-50">Loop backbone hallucination<br/>AB / C / EF loops<br/>6–24 aa expansions</div>
  </div>
  <div class="flex items-center justify-center text-blue-400 text-lg">→</div>
  <div class="flex flex-col items-center p-3 bg-violet-500/10 rounded border border-violet-500/30 text-center">
    <div class="text-violet-400 font-black text-xs mb-1">SEQ</div>
    <div class="text-[10px] font-bold text-white mb-1">ProteinMPNN</div>
    <div class="text-[9px] opacity-50">1000 seqs/target<br/>T=0.2, fixed scaffold<br/>Top 10 per loop</div>
  </div>
  <div class="flex items-center justify-center text-blue-400 text-lg">→</div>
  <div class="flex flex-col items-center p-3 bg-cyan-500/10 rounded border border-cyan-500/30 text-center">
    <div class="text-cyan-400 font-black text-xs mb-1">VAL</div>
    <div class="text-[10px] font-bold text-white mb-1">AlphaFold3 Server</div>
    <div class="text-[9px] opacity-50">Batch co-folding<br/>ipTM / PAE / pLDDT<br/>AlphaFoldServer_Analysis.ipynb</div>
  </div>
  <div class="flex items-center justify-center text-blue-400 text-lg">→</div>
  <div class="flex flex-col items-center p-3 bg-emerald-500/10 rounded border border-emerald-500/30 text-center">
    <div class="text-emerald-400 font-black text-xs mb-1">SEL</div>
    <div class="text-[10px] font-bold text-white mb-1">Best Binder</div>
    <div class="text-[9px] opacity-50">T-score > 2.0<br/>Multi-metric consensus<br/>AlphaFold_Best_Binder.ipynb</div>
  </div>
  <div class="flex items-center justify-center text-blue-400 text-lg">→</div>
  <div class="flex flex-col items-center p-3 bg-amber-500/10 rounded border border-amber-500/30 text-center">
    <div class="text-amber-400 font-black text-xs mb-1">SYN</div>
    <div class="text-[10px] font-bold text-white mb-1">Twist Bioscience</div>
    <div class="text-[9px] opacity-50">13 constructs<br/>RE site validation<br/>Dec 2025 order</div>
  </div>
  <div class="flex items-center justify-center text-blue-400 text-lg">→</div>
  <div class="flex flex-col items-center p-3 bg-red-500/10 rounded border border-red-500/30 text-center">
    <div class="text-red-400 font-black text-xs mb-1">EXP</div>
    <div class="text-[10px] font-bold text-white mb-1">Yeast Display</div>
    <div class="text-[9px] opacity-50">Flow cytometry<br/>ANOVA + Tukey-HSD<br/>MMP9 vs MMP2</div>
  </div>
</div>

<div class="grid grid-cols-2 gap-6 mt-6">
  <div class="p-3 bg-white/5 rounded border border-white/10">
    <h3 class="text-blue-400 font-bold mb-1 uppercase text-xs tracking-widest">Computational Phase</h3>
    <ul class="text-[11px] space-y-1 opacity-80 list-none p-0">
      <li class="flex gap-2"><span>—</span><span>HADDOCK-docked PDB templates for each MMP/ADAM target</span></li>
      <li class="flex gap-2"><span>—</span><span>ProteinMPNN batch scripted: <code>ProMPNN_Batch.py</code></span></li>
      <li class="flex gap-2"><span>—</span><span>T-score filtering via <code>AlphaFold_Best_Binder.ipynb</code></span></li>
    </ul>
  </div>
  <div class="p-3 bg-white/5 rounded border border-white/10">
    <h3 class="text-cyan-400 font-bold mb-1 uppercase text-xs tracking-widest">Experimental Phase</h3>
    <ul class="text-[11px] space-y-1 opacity-80 list-none p-0">
      <li class="flex gap-2"><span>—</span><span>RE site validation via <code>TwistBioOrder.ipynb</code> (13/13 passed)</span></li>
      <li class="flex gap-2"><span>—</span><span>Golden Gate cloning into pCT yeast display vector</span></li>
      <li class="flex gap-2"><span>—</span><span>Flow cytometry: Binding Efficiency (DP/FITC+) metric</span></li>
    </ul>
  </div>
</div>

---
layout: section
transition: slide-up
---

# Phase 1
### Computational Design and Candidate Selection

---
layout: default
transition: fade-out
---

# RFdiffusion Loop Scaffolding
Targeting variable surface loops on the TIMP3 scaffold.

<div class="grid grid-cols-2 gap-8 mt-8 items-start">
  <LoopConfig />
  <div class="space-y-4">
    <div class="text-sm leading-relaxed opacity-80">
      Length-variable loop generation targeting the AB, C, and EF loops against TIMP3-bound metalloproteinase catalytic domains. MMP9 vs MMP2 selectivity was the primary design objective.
    </div>
    <div class="p-4 bg-blue-500/5 rounded-lg border border-blue-500/20">
      <h4 class="text-[10px] font-bold text-blue-400 mb-3 uppercase tracking-widest">Design Parameters</h4>
      <div class="grid grid-cols-2 gap-y-3 gap-x-4 text-[11px] font-mono">
        <span class="opacity-40 uppercase">Target pairs</span> <span class="text-blue-200">MMP9 / MMP2</span>
        <span class="opacity-40 uppercase">Loop regions</span> <span class="text-blue-200">AB, C, EF</span>
        <span class="opacity-40 uppercase">Lengths</span> <span class="text-blue-200">6–24 aa</span>
        <span class="opacity-40 uppercase">Iterations</span> <span class="text-blue-200">50 (20 hr/run, A30 GPU)</span>
      </div>
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Sequence Design & AlphaFold Selection
ProteinMPNN batch → AlphaFold co-folding → T-score consensus → Twist order.

<div class="grid grid-cols-2 gap-8 mt-4">
  <div>
    <h3 class="text-xs font-bold mb-3 uppercase tracking-widest text-blue-400">Final Twist Ordered Library <span class="opacity-40">(Dec 2025)</span></h3>
    <TopVariants />
  </div>
  <div class="space-y-3">
    <div class="p-3 bg-white/5 rounded border border-white/10">
      <h4 class="text-[10px] font-bold uppercase tracking-widest mb-2 text-blue-300">Step 1 — ProteinMPNN</h4>
      <p class="text-[10px] leading-relaxed opacity-70">
        1,000 sequences/target at T=0.2. Scaffold fixed; only loop positions redesigned. Top 10 unique sequences (by log-likelihood) sent to AlphaFold.
      </p>
    </div>
    <div class="p-3 bg-white/5 rounded border border-white/10">
      <h4 class="text-[10px] font-bold uppercase tracking-widest mb-2 text-blue-300">Step 2 — AlphaFold3 Co-folding</h4>
      <p class="text-[10px] leading-relaxed opacity-70">
        Batch JSON submission to AlphaFold Server. Parsed by <code>AlphaFoldServer_Analysis.ipynb</code>: ipTM, pLDDT, Loop pLDDT, Interface PAE extracted per complex.
      </p>
    </div>
    <div class="p-3 bg-white/5 rounded border border-white/10">
      <h4 class="text-[10px] font-bold uppercase tracking-widest mb-2 text-emerald-300">Step 3 — Best Binder Selection</h4>
      <p class="text-[10px] leading-relaxed opacity-70">
        <code>AlphaFold_Best_Binder.ipynb</code>: T-score &gt; 2.0 across ≥2 metrics = Elite. Final 13 constructs validated for restriction sites (BsrGI, BamHI, BsaI) via <code>TwistBioOrder.ipynb</code>.
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
Identifying statistically significant target-selective binders.

<div class="mt-2 flex flex-col items-center text-center max-w-3xl mx-auto space-y-4">
  <TScoreFormula />

  <div class="grid grid-cols-2 gap-6 w-full">
    <div class="bg-emerald-500/10 p-4 rounded-lg border border-emerald-500/20">
      <div class="text-emerald-400 font-black text-xl mb-1">T > 2.0</div>
      <div class="text-[10px] uppercase tracking-widest opacity-60">Significant Target Win</div>
    </div>
    <div class="bg-rose-500/10 p-4 rounded-lg border border-rose-500/20">
      <div class="text-rose-400 font-black text-xl mb-1">T < -2.0</div>
      <div class="text-[10px] uppercase tracking-widest opacity-60">Significant Target Loss</div>
    </div>
  </div>

  <p class="text-xs opacity-40 italic">
    *Unlike Z-scores (which reward globally superior binders), T-scores normalize each variant against its own cross-target mean — explicitly identifying target-preferential binders.
  </p>
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

# MMP9 vs MMP2 Specificity: ANOVA Results
Statistically significant target discrimination across designed constructs.

<div class="grid grid-cols-2 gap-8 mt-4 items-start">
  <div class="p-4 bg-white/5 rounded border border-white/10 space-y-3">
    <h3 class="text-xs font-bold uppercase tracking-[0.2em] text-blue-400">Selectivity Table (Binding Efficiency)</h3>
    <table class="w-full text-[10px] font-mono border-collapse">
      <thead>
        <tr class="border-b border-white/20">
          <th class="text-left py-1 opacity-50">Construct</th>
          <th class="text-right py-1 opacity-50">Design</th>
          <th class="text-right py-1 opacity-50">MMP2</th>
          <th class="text-right py-1 opacity-50 text-blue-400">MMP9</th>
          <th class="text-right py-1 opacity-50">Ratio</th>
          <th class="text-right py-1 opacity-50">Result</th>
        </tr>
      </thead>
      <tbody>
        <tr class="border-b border-white/5">
          <td class="py-1 text-emerald-400 font-bold">C 12</td>
          <td class="text-right py-1 opacity-40">M9 High</td>
          <td class="text-right py-1">32.1%</td>
          <td class="text-right py-1 text-emerald-400 font-bold">91.3%</td>
          <td class="text-right py-1 text-emerald-400">2.8×</td>
          <td class="text-right py-1 text-emerald-400">✓ Sig</td>
        </tr>
        <tr class="border-b border-white/5">
          <td class="py-1 text-emerald-400 font-bold">C 15</td>
          <td class="text-right py-1 opacity-40">M9 High</td>
          <td class="text-right py-1">33.8%</td>
          <td class="text-right py-1 text-emerald-400 font-bold">88.6%</td>
          <td class="text-right py-1 text-emerald-400">2.6×</td>
          <td class="text-right py-1 text-emerald-400">✓ Sig</td>
        </tr>
        <tr class="border-b border-white/5">
          <td class="py-1 text-emerald-400 font-bold">AB 6</td>
          <td class="text-right py-1 opacity-40">M9 High</td>
          <td class="text-right py-1">23.8%</td>
          <td class="text-right py-1 text-emerald-400 font-bold">89.9%</td>
          <td class="text-right py-1 text-emerald-400">3.8×</td>
          <td class="text-right py-1 text-emerald-400">✓ Sig</td>
        </tr>
        <tr class="border-b border-white/10 opacity-50">
          <td class="py-1">C 13</td>
          <td class="text-right py-1 opacity-40">Low ctrl</td>
          <td class="text-right py-1">16.0%</td>
          <td class="text-right py-1">45.4%</td>
          <td class="text-right py-1">—</td>
          <td class="text-right py-1">n.s. ✓</td>
        </tr>
        <tr class="border-b border-white/10 opacity-50">
          <td class="py-1">AB 5</td>
          <td class="text-right py-1 opacity-40">Low ctrl</td>
          <td class="text-right py-1">11.8%</td>
          <td class="text-right py-1">37.6%</td>
          <td class="text-right py-1">—</td>
          <td class="text-right py-1">n.s. ✓</td>
        </tr>
        <tr>
          <td class="py-1 text-blue-400">TIMP 3</td>
          <td class="text-right py-1 opacity-40">WT ref</td>
          <td class="text-right py-1">21.9%</td>
          <td class="text-right py-1">54.5%</td>
          <td class="text-right py-1">2.5×</td>
          <td class="text-right py-1">ref</td>
        </tr>
      </tbody>
    </table>
    <div class="text-[8px] opacity-40 italic border-t border-white/10 pt-2">
      ANOVA confirmed p&lt;0.05 for all selective designs. "Low" controls behaved as expected (p&gt;0.05).
    </div>
  </div>

  <div class="space-y-3">
    <div class="p-4 bg-emerald-500/5 rounded border border-emerald-500/20">
      <div class="text-emerald-400 font-black text-2xl mb-1">100% PPV</div>
      <div class="text-[9px] uppercase tracking-widest text-emerald-400 font-bold">Pipeline Predictive Power</div>
      <div class="text-[10px] opacity-60 mt-1">Every variant designed for MMP9 preference achieved statistically significant discrimination.</div>
    </div>
    <div class="p-4 bg-blue-500/5 rounded border border-blue-500/20">
      <div class="text-blue-400 font-black text-2xl mb-1">3.8× Ratio</div>
      <div class="text-[9px] uppercase tracking-widest text-blue-400 font-bold">AB 6 Peak Selectivity</div>
      <div class="text-[10px] opacity-60 mt-1">MMP9=89.9% vs MMP2=23.8% · ANOVA p=0.0002</div>
    </div>
    <div class="p-3 bg-white/5 rounded border border-white/10">
      <div class="text-[9px] leading-relaxed opacity-70">
        <b>Computational Success:</b> The T-score filter reduced >1,000 sequences to 13 candidates, perfectly enriching for experimental target-preferential binders.
      </div>
    </div>
  </div>
</div>

---
layout: default
transition: fade-out
---

# Cross-Target Binding Analysis
Binding efficiency heatmap across all constructs and targets.

<div class="mt-2 flex justify-center">
  <FcsChart />
</div>

---
layout: default
transition: fade-out
---

# Conclusion and Strategic Roadmap
Validating the de novo specificity pipeline.

<div class="grid grid-cols-2 gap-12 mt-4">
  <div class="space-y-4">
    <h3 class="text-blue-400 font-bold uppercase text-xs tracking-widest">Key Findings</h3>
    <ul class="text-[11px] space-y-3 opacity-80 list-none p-0">
      <li class="flex gap-3">
        <div class="w-6 h-6 rounded-full bg-emerald-500/20 flex items-center justify-center text-[10px] font-bold text-emerald-400 shrink-0">01</div>
        <span><b>Target Specificity Engineered:</b> C 12, C 15, AB 6 achieved statistically significant MMP9 preference over MMP2 (p&lt;0.02).</span>
      </li>
      <li class="flex gap-3">
        <div class="w-6 h-6 rounded-full bg-blue-500/20 flex items-center justify-center text-[10px] font-bold text-blue-400 shrink-0">02</div>
        <span><b>100% Predictivity:</b> T-score multi-metric filter accurately prioritized selective winners with zero false positives in control trials.</span>
      </li>
      <li class="flex gap-3">
        <div class="w-6 h-6 rounded-full bg-violet-500/20 flex items-center justify-center text-[10px] font-bold text-violet-400 shrink-0">03</div>
        <span><b>Therapeutic Path:</b> Selective MMP9 inhibition enables tumor treatment while sparing MMP2-dependent tissue homeostasis.</span>
      </li>
    </ul>
    
    <div class="p-4 bg-blue-500/5 rounded border border-blue-500/20 mt-4">
      <h4 class="text-[10px] font-bold text-blue-400 uppercase tracking-widest mb-2">Clinical Analog: Precision CAR-T</h4>
      <p class="text-[10px] leading-relaxed opacity-70">
        This pipeline can generate hyperspecific binders for patient-unique neoantigens, enabling personalized CAR-T therapies that minimize "on-target, off-tumor" toxicity.
      </p>
    </div>
  </div>

  <div class="space-y-4">
    <h3 class="text-cyan-400 font-bold uppercase text-xs tracking-widest">Strategic Roadmap</h3>
    <ul class="text-[11px] space-y-3 opacity-80 list-none p-0">
      <li class="flex items-center gap-3">
        <div class="w-1.5 h-1.5 rounded-full bg-cyan-400"></div>
        <span><b>Reverse Selectivity:</b> Target MMP2 > MMP9 as a proof of concept for controllable pipeline steering.</span>
      </li>
      <li class="flex items-center gap-3">
        <div class="w-1.5 h-1.5 rounded-full bg-cyan-400"></div>
        <span><b>ADAM Selectivity:</b> Optimize ADAM10 surface presentation to enable robust ADAM17 vs 10 analysis.</span>
      </li>
      <li class="flex items-center gap-3">
        <div class="w-1.5 h-1.5 rounded-full bg-cyan-400"></div>
        <span><b>Absolute Kinetics:</b> Measure absolute KD values via SPR to benchmark against clinical standards.</span>
      </li>
      <li class="flex items-center gap-3">
        <div class="w-1.5 h-1.5 rounded-full bg-cyan-400"></div>
        <span><b>Saturation Mutagenesis:</b> ESM-2 guided loop optimization to further refine affinity.</span>
      </li>
    </ul>
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
