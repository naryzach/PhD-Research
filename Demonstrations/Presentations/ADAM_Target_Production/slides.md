---
theme: seriph
background: https://images.unsplash.com/photo-1576086213369-97a306d36557?q=80&w=1920&auto=format&fit=crop
class: text-white
drawings:
  persist: false
transition: slide-left
mdc: true
title: "Recombinant Production of ADAM10, ADAM17, and MMP7"
---

<div class="h-full flex flex-col justify-center items-center bg-black/60 p-8 rounded-xl backdrop-blur-sm border border-white/10">
  <h1 class="text-4xl font-extrabold text-transparent bg-clip-text bg-gradient-to-r from-teal-400 to-blue-500 text-center leading-tight mb-4">
    Recombinant Production of ADAM10, ADAM17, and MMP7
  </h1>
  <p class="text-xl text-gray-300 font-light text-center mb-8">
    Target Metalloproteinases for Yeast-Display Selectivity Screening
  </p>
  <div class="text-sm text-gray-400 flex flex-col items-center">
    <span>Ryan Gustafson</span>
    <span class="opacity-70 mt-1">Sarmazdeh Lab | PhD Research</span>
  </div>
</div>

---
layout: default
---

# Motivation & Context

<div class="grid grid-cols-2 gap-8 mt-6">
  <div>
    <h2 class="text-teal-600 font-semibold mb-2">Expansion of Selectivity Panels</h2>
    <p class="text-sm leading-relaxed text-gray-700 mb-4">
      The primary de novo TIMP3 binder campaign targets metalloproteinases. To expand this panel and perform robust selectivity assays (specifically ADAM10 vs. ADAM17), we require high-quality in-house recombinant target catalytic domains.
    </p>
    <h2 class="text-teal-600 font-semibold mb-2">Issues with Commercial Reagents</h2>
    <p class="text-sm leading-relaxed text-gray-700">
      Commercial ADAM10 material has yielded insufficient positive-control binding (<3% Double+ cells), creating a critical bottleneck in sorting pipelines.
    </p>
  </div>
  <div class="bg-gray-50 p-6 rounded-lg border border-gray-150 flex flex-col justify-between">
    <h3 class="font-bold text-gray-800 mb-2">Production Targets</h3>
    <ul class="space-y-3 text-sm text-gray-600">
      <li><strong>ADAM10 catalytic domain:</strong> Periplasmic expression (pMopac16) + FPLC purification.</li>
      <li><strong>ADAM17 catalytic domain:</strong> Co-transformed with helper DsbC plasmid.</li>
      <li><strong>MMP7 target:</strong> Stepwise refolding via dialysis.</li>
    </ul>
    <div class="mt-4 p-3 bg-teal-50 border-l-4 border-teal-500 rounded text-xs text-teal-800">
      <strong>Goal:</strong> Generate active, display-grade metalloproteinases to evaluate binder specificity.
    </div>
  </div>
</div>

---

# Periplasmic Expression Optimization

<div class="grid grid-cols-12 gap-6 items-center">
  <div class="col-span-5">
    <h2 class="text-teal-600 font-semibold mb-2">Dual-Plasmid Periplasmic System</h2>
    <p class="text-xs text-gray-700 mb-3 leading-relaxed">
      Induction of target catalytic domains (arabinose-inducible pMopac16) co-transformed with disulfide isomerase C (DsbC) helper plasmid to promote correct folding in <em>E. coli</em>.
    </p>
    <h2 class="text-teal-600 font-semibold mb-2">Temperature Sensitivity</h2>
    <p class="text-xs text-gray-700 mb-3 leading-relaxed">
      A comparison of 25°C and 30°C induction showed that 30°C induction drives the target predominantly into insoluble inclusion bodies (lysis pellet), whereas <strong>25°C induction reduces aggregation</strong> and favors soluble periplasmic accumulation.
    </p>
    <div class="p-2 bg-yellow-50 border border-yellow-200 rounded text-2xs text-yellow-800">
      <strong>Standard protocol:</strong> 25°C overnight induction is utilized as the baseline for both targets.
    </div>
  </div>
  <div class="col-span-7 flex flex-col items-center">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/ADAM1017_SDS-Page_Gel_20260603.png" class="h-68 object-contain rounded border shadow-sm" alt="ADAM17 Expression Gel"/>
    <span class="text-2xs text-gray-500 mt-2">Coomassie Gel (June 3): Overnight induction (Lane 13) confirms strong target expression.</span>
  </div>
</div>

---

# AEC FPLC Purification Profiles

<div class="grid grid-cols-2 gap-4 mt-4">
  <div>
    <h2 class="text-teal-600 font-semibold mb-2">Anion-Exchange (AEC)</h2>
    <p class="text-xs text-gray-700 leading-relaxed mb-3">
      Osmotic shock periplasmic extracts were loaded onto a Cytiva HiTrap Q FF 5mL column. 
      Buffer: 20 mM Tris-HCl, 5 mM CaCl₂, 3 µM ZnCl₂, pH 8.5. Gradient: 5-95% NaCl over 15 CV.
    </p>
    <h2 class="text-red-600 font-semibold mb-2">Low-Affinity Binding Caveat</h2>
    <p class="text-xs text-gray-700 leading-relaxed">
      Both ADAM10 and ADAM17 exhibit a large UV peak in the flow-through during sample loading, indicating minimal column retention at pH 8.5. Enriched fractions are recovered in the high-NaCl strip, suggesting low affinity.
    </p>
  </div>
  <div class="grid grid-rows-2 gap-2">
    <div class="flex items-center gap-2 border p-2 rounded bg-gray-50">
      <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/ADAM_10_AEC_T4.png" class="h-28 object-contain" alt="ADAM10 AEC"/>
      <div class="text-2xs text-gray-600"><strong>ADAM10 AEC profile:</strong> Large flow-through peak with secondary strip elution peak.</div>
    </div>
    <div class="flex items-center gap-2 border p-2 rounded bg-gray-50">
      <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/ADAM_17_AEC_T2.png" class="h-28 object-contain" alt="ADAM17 AEC"/>
      <div class="text-2xs text-gray-600"><strong>ADAM17 AEC profile:</strong> Identical wash-through behavior, prompting silver-stain validation.</div>
    </div>
  </div>
</div>

---

# Spin Concentration & Yield Quantification

<div class="mt-4">
  <p class="text-sm text-gray-700 mb-4">
    FPLC fractions were concentrated using Amicon Ultra-4 10 kDa MWCO filters (4000 rpm, 30 min, 4°C). Yields were quantified via NanoDrop A280 protein absorbance.
  </p>

  <table class="w-full text-xs text-left border border-collapse border-gray-200">
    <thead>
      <tr class="bg-teal-500 text-white">
        <th class="p-2 border">Sample Name</th>
        <th class="p-2 border">Protein (mg/mL)</th>
        <th class="p-2 border">A280 Absorbance</th>
        <th class="p-2 border">A260/A280 Ratio</th>
        <th class="p-2 border">Status / Outcome</th>
      </tr>
    </thead>
    <tbody class="text-gray-700">
      <tr class="bg-white">
        <td class="p-2 border font-medium">A10 FXN 6 Unconcentrated</td>
        <td class="p-2 border">0.163</td>
        <td class="p-2 border">0.163</td>
        <td class="p-2 border">0.963</td>
        <td class="p-2 border" rowspan="3" style="vertical-align: middle;"><strong>Membrane Leakage:</strong> Filtrate concentration exceeded residue, indicating filter seal failure.</td>
      </tr>
      <tr class="bg-gray-50">
        <td class="p-2 border font-medium">A10 FXN 6 Filtrate</td>
        <td class="p-2 border">0.176</td>
        <td class="p-2 border">0.176</td>
        <td class="p-2 border">0.981</td>
      </tr>
      <tr class="bg-white">
        <td class="p-2 border font-medium">A10 FXN 6 Residue (3.2×)</td>
        <td class="p-2 border">0.519</td>
        <td class="p-2 border">0.519</td>
        <td class="p-2 border">0.779</td>
      </tr>
      <tr class="bg-gray-50">
        <td class="p-2 border font-medium">A10 FXN 9 Unconcentrated</td>
        <td class="p-2 border">0.141</td>
        <td class="p-2 border">0.141</td>
        <td class="p-2 border">1.309</td>
        <td class="p-2 border" rowspan="3" style="vertical-align: middle;"><strong>Successful Enrichment:</strong> achieved 6.1× concentration with minimal loss in filtrate (0.063 mg/mL).</td>
      </tr>
      <tr class="bg-white">
        <td class="p-2 border font-medium">A10 FXN 9 Filtrate</td>
        <td class="p-2 border">0.063</td>
        <td class="p-2 border">0.063</td>
        <td class="p-2 border">1.358</td>
      </tr>
      <tr class="bg-gray-50">
        <td class="p-2 border font-medium">A10 FXN 9 Residue (6.1×)</td>
        <td class="p-2 border">0.863</td>
        <td class="p-2 border">0.863</td>
        <td class="p-2 border">0.803</td>
      </tr>
    </tbody>
  </table>
</div>

---

# Optimized Silver-Stain Validation (June 10)

<div class="grid grid-cols-2 gap-4 items-center">
  <div class="flex flex-col items-center">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/ADAM1017_SDS-Page_Gel_20260610_ADAM10.png" class="h-68 object-contain rounded border" alt="ADAM10 Rerun"/>
    <span class="text-2xs text-gray-500 mt-1">ADAM10 Silver Stain Rerun (FXN 9 in Lane 14)</span>
  </div>
  <div class="flex flex-col items-center">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/ADAM1017_SDS-Page_Gel_20260610_ADAM17.png" class="h-68 object-contain rounded border" alt="ADAM17 Rerun"/>
    <span class="text-2xs text-gray-500 mt-1">ADAM17 Silver Stain Rerun (FXN 6 in Lane 14)</span>
  </div>
</div>

<div class="text-2xs text-gray-600 mt-2 text-center">
  Rerun resolved previous ladder issues using PageRuler Plus and optimized gel chemistry (Tris base). Both gels show distinct candidate bands corresponding to target catalytic domains in the 30-40 kDa range.
</div>

---

# Flow Cytometry Binding Trials (June 10)

<div class="grid grid-cols-2 gap-8 mt-2">
  <div>
    <h2 class="text-teal-600 font-semibold mb-2">Display Target Gating & Spillover</h2>
    <p class="text-xs text-gray-700 leading-relaxed mb-4">
      TIMP3-display variants were incubated with target proteins. Batch processing applied a <strong>spillover coefficient (Alpha) of 0.1638</strong> to correct FITC bleed-through into the PE-FLAG channel.
    </p>
    <h2 class="text-teal-600 font-semibold mb-2">MMP9 Control Validation</h2>
    <p class="text-xs text-gray-700 leading-relaxed">
      Positive control <strong>TIMP 3-MMP9(Sino)</strong> showed strong binding (<strong>38.2% Double+</strong>, MFI: 491). Twist variants <strong>AB 7</strong> (33.5%) and <strong>AB 3</strong> (32.9%) also bound commercial Sino MMP9 strongly.
    </p>
  </div>
  <div class="bg-gray-50 p-4 rounded border flex flex-col justify-between">
    <div>
      <h3 class="font-bold text-gray-800 text-sm mb-2">Low ADAM10 Target Binding</h3>
      <p class="text-xs text-gray-600 leading-relaxed mb-2">
        Incubation with in-house ADAM10 fractions (F6 and F9) yielded negligible binding (<0.9% Double+ cells) for all display variants.
      </p>
      <p class="text-xs text-gray-600 leading-relaxed">
        Commercial AbCam ADAM10 also showed very low binding (<0.6% Double+), suggesting that the ADAM10 targets are inactive or require different fold/cofactor conditions.
      </p>
    </div>
    <div class="mt-4 p-2 bg-red-50 border-l-4 border-red-500 text-2xs text-red-800">
      <strong>Action Item:</strong> Optimize target folding, buffer cofactors (Zn/Ca), or evaluate tag accessibility.
    </div>
  </div>
</div>

---

# Cross-Lab Reference: Masoud's MMP9

<div class="grid grid-cols-12 gap-6 mt-4 items-center">
  <div class="col-span-5">
    <h2 class="text-teal-600 font-semibold mb-2">Benchmarking Reagent Sources</h2>
    <p class="text-xs text-gray-700 leading-relaxed mb-3">
      Masoud's in-house MMP9 target (His-tagged) was tested alongside commercial Sino MMP9.
    </p>
    <h2 class="text-teal-600 font-semibold mb-2">Reduced Binding Affinity</h2>
    <p class="text-xs text-gray-700 leading-relaxed">
      Twist variants showed positive but reduced binding to Masoud's MMP9: 
      <strong>AB 7-MMP9(Mas)</strong> achieved <strong>12.8% Double+</strong> (MFI: 319) vs. 33.5% against Sino MMP9. This indicates that while active, the local preparation has lower binding yield.
    </p>
  </div>
  <div class="col-span-7 bg-white p-4 rounded border shadow-sm">
    <div class="text-xs font-bold text-gray-800 mb-2">Top Binders Comparison (Double+ % / Bind MFI)</div>
    <div class="space-y-2 text-2xs">
      <div class="flex justify-between border-b pb-1">
        <span class="font-medium text-gray-700">TIMP 3 Control (Sino MMP9)</span>
        <span class="text-teal-600 font-semibold">38.2% / 491</span>
      </div>
      <div class="flex justify-between border-b pb-1">
        <span class="font-medium text-gray-700">AB 7 Variant (Sino MMP9)</span>
        <span class="text-teal-600 font-semibold">33.5% / 445</span>
      </div>
      <div class="flex justify-between border-b pb-1">
        <span class="font-medium text-gray-700">AB 3 Variant (Sino MMP9)</span>
        <span class="text-teal-600 font-semibold">32.9% / 456</span>
      </div>
      <div class="flex justify-between border-b pb-1 text-gray-500">
        <span>AB 7 Variant (Masoud's MMP9)</span>
        <span class="font-semibold">12.8% / 319</span>
      </div>
      <div class="flex justify-between border-b pb-1 text-gray-500">
        <span>C 11 Variant (Masoud's MMP9)</span>
        <span class="font-semibold">10.9% / 326</span>
      </div>
    </div>
  </div>
</div>

# Western Blot Validation (June 11-12)

<div class="grid grid-cols-2 gap-6 mt-6 items-center">
  <div class="space-y-3">
    <h2 class="text-teal-600 font-semibold mb-2">Soluble ADAM10 Expression Verified</h2>
    <p class="text-xs text-gray-700 leading-relaxed mb-3">
      Anti-FLAG Western blotting of concentrated culture supernatants confirmed <strong>soluble ADAM10 expression</strong>. A single, distinct soluble band was detected near <strong>60 kDa</strong> in the 25°C supernatant, indicating successful target production and folding.
    </p>
    <div class="p-2 bg-emerald-50 border border-emerald-200 rounded text-2xs text-emerald-800">
      <strong>Transfer Quality Check:</strong> Ponceau S staining of the PVDF membrane confirmed uniform total protein transfer across all lanes, validating blotting efficiency.
    </div>
  </div>
  <div class="flex flex-col items-center">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/western_blot_labeled.png" class="h-56 object-contain rounded border shadow-sm" alt="ADAM10 Western Blot"/>
    <span class="text-2xs text-gray-500 mt-2">Anti-FLAG Western Blot: Soluble band near 60 kDa in Lane 3.</span>
  </div>
</div>

---

# BCA Assay Target Standardization (June 12)

<div class="grid grid-cols-2 gap-6 mt-6 items-center">
  <div class="space-y-3">
    <h2 class="text-teal-600 font-semibold mb-2">Target Concentration Calibration</h2>
    <p class="text-xs text-gray-700 leading-relaxed mb-3">
      BCA assays were executed on the concentrated FPLC Fraction 6 of ADAM10 to standardize target inputs for downstream display and binding trials, resolving NanoDrop buffer interference.
    </p>
    <div class="p-3 bg-gray-50 rounded border border-gray-200 text-2xs">
      <ul class="space-y-1.5 text-gray-600">
        <li><strong>SpectraMax i3x:</strong> Quadratic fit ($R^2 = 0.9895$) yielded a dilution-adjusted mean of <strong>205.46 μg/mL</strong>.</li>
        <li><strong>SPECTRAmax M5:</strong> Linear fit ($R^2 = 0.9930$) yielded a dilution-adjusted mean of <strong>193.89 μg/mL</strong>.</li>
        <li><strong>Standardization Outcome:</strong> Purified target concentration was standardized to <strong>0.20 mg/mL</strong>.</li>
      </ul>
    </div>
  </div>
  <div class="flex items-center gap-2 border p-2 rounded bg-gray-50">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/20260612_2_standard_curve.png" class="h-52 object-contain" alt="BCA Standard Curve"/>
    <div class="text-2xs text-gray-600"><strong>i3x Calibration:</strong> Standard curve fitting for precise concentration calculations.</div>
  </div>
</div>

---

# Retry ADAM10 Expression Trials (June 16)

<div class="grid grid-cols-12 gap-6 mt-4 items-center">
  <div class="col-span-5">
    <h2 class="text-teal-600 font-semibold mb-2">Optimizing Soluble Yields</h2>
    <p class="text-xs text-gray-700 leading-relaxed mb-3">
      To address low binding levels and optimize soluble periplasmic yields, a secondary expression trial was initiated to evaluate helper plasmid variants and lysis methods.
    </p>
    <h2 class="text-teal-600 font-semibold mb-2">Trial Conditions & Status</h2>
    <p class="text-xs text-gray-700 leading-relaxed">
      Four conditions were set up evaluating combinations of <strong>BL21(DE3)</strong> host strains, <strong>pMopac-ADAM-10cd-Flg</strong>, and <strong>pBAD-helper</strong> vectors. Samples were collected at 1hr and 3hr post-induction. SDS-PAGE and Western blot validation runs are currently pending.
    </p>
  </div>
  <div class="col-span-7 bg-white p-4 rounded border shadow-sm">
    <div class="text-xs font-bold text-gray-800 mb-2">Expression Trial Design Matrix</div>
    <table class="w-full text-left text-2xs border-collapse">
      <thead>
        <tr class="border-b text-teal-600 font-bold">
          <th class="py-1">Condition</th>
          <th class="py-1">Strain / Plasmids</th>
          <th class="py-1">Induction Temp</th>
          <th class="py-1">Lysis / Prep</th>
        </tr>
      </thead>
      <tbody>
        <tr class="border-b">
          <td class="py-1">1</td>
          <td class="py-1">BL21(DE3) + pMopac-ADAM10</td>
          <td class="py-1">25°C, 3.5 hr</td>
          <td class="py-1">Osmotic Shock</td>
        </tr>
        <tr class="border-b">
          <td class="py-1">2</td>
          <td class="py-1">BL21(DE3) + pMopac-ADAM10 + pBAD-helper</td>
          <td class="py-1">25°C, 3.5 hr</td>
          <td class="py-1">Osmotic Shock + DsbC</td>
        </tr>
        <tr class="border-b">
          <td class="py-1">3</td>
          <td class="py-1">BL21(DE3) + pMopac-ADAM10</td>
          <td class="py-1">25°C, 3.5 hr</td>
          <td class="py-1">Total Lysate (Sonicated)</td>
        </tr>
        <tr>
          <td class="py-1">4</td>
          <td class="py-1">BL21(DE3) + pMopac-ADAM10 + pBAD-helper</td>
          <td class="py-1">25°C, 3.5 hr</td>
          <td class="py-1">Total Lysate + DsbC</td>
        </tr>
      </tbody>
    </table>
  </div>
</div>

---

# pMopac DNA Plasmid Audit (June 17)

<div class="grid grid-cols-12 gap-6 mt-4 items-center">
  <div class="col-span-5">
    <h2 class="text-teal-600 font-semibold mb-2">Vector Topology Verification</h2>
    <p class="text-xs text-gray-700 leading-relaxed mb-3">
      Audited DNA sequence files to resolve physical property discrepancies and verify cloning sites (SfiI restriction sites) and tag configurations.
    </p>
    <div class="p-2.5 bg-blue-50 border border-blue-200 rounded text-2xs text-blue-800">
      <strong>Purification Insights:</strong> Theoretical isoelectric points (pI 5.95 for ADAM17, pI 6.44 for ADAM10) suggest negative charges at pH 8.0, confirming suitability of HiTrap Q FF Anion Exchange chromatography.
    </div>
  </div>
  <div class="col-span-7 bg-white p-4 rounded border shadow-sm">
    <div class="text-xs font-bold text-gray-800 mb-2">Verified DNA Annotation Metadata</div>
    <div class="space-y-2 text-2xs">
      <div class="flex justify-between border-b pb-1">
        <span class="font-medium text-gray-700"><strong>(pET-)ADAM-17cd-HT</strong> (T7 / KanR / 6xHis)</span>
        <span class="text-teal-600 font-semibold">30.3 kDa / pI 5.95</span>
      </div>
      <div class="flex justify-between border-b pb-1">
        <span class="font-medium text-gray-700"><strong>pCHA-ADAM17cd</strong> (GAL1,10 / AmpR / Myc)</span>
        <span class="text-teal-600 font-semibold">42.4 kDa / pI 5.56</span>
      </div>
      <div class="flex justify-between border-b pb-1">
        <span class="font-medium text-gray-700"><strong>pMopac-ADAM-10cd-Flg</strong> (lac / CmR / FLAG)</span>
        <span class="text-teal-600 font-semibold">29.6 kDa / pI 6.44</span>
      </div>
    </div>
  </div>
</div>

---

# ADAM17 Plasmid Miniprep (June 18)

<div class="grid grid-cols-12 gap-6 mt-4 items-center">
  <div class="col-span-5">
    <h2 class="text-teal-600 font-semibold mb-2">Plasmid Isolation & Purity</h2>
    <p class="text-xs text-gray-700 leading-relaxed mb-3">
      Isolated the expression plasmid for <strong>(pET-)ADAM-17cd-HT</strong> using the Promega Wizard® Plus SV Miniprep kit.
    </p>
    <div class="p-2 bg-emerald-50 border border-emerald-200 rounded text-2xs text-emerald-800">
      <strong>Sequence Validation:</strong> Isolated plasmid DNA was submitted to Eurofins Genomics for whole plasmid nanopore sequencing to verify construct integrity.
    </div>
  </div>
  <div class="col-span-7 bg-white p-4 rounded border shadow-sm">
    <div class="text-xs font-bold text-gray-800 mb-2">Spectrophotometric DNA Yields & Ratios</div>
    <table class="w-full text-left text-2xs border-collapse">
      <thead>
        <tr class="border-b text-teal-600 font-bold">
          <th class="py-1">Parameter</th>
          <th class="py-1">Measured Value</th>
          <th class="py-1">Quality Assessment</th>
        </tr>
      </thead>
      <tbody>
        <tr class="border-b">
          <td class="py-1"><strong>DNA Concentration</strong></td>
          <td class="py-1">99.5 ng/μL</td>
          <td class="py-1 text-emerald-600 font-semibold">High Yield (9.95 μg total)</td>
        </tr>
        <tr class="border-b">
          <td class="py-1"><strong>A260 / A280 Ratio</strong></td>
          <td class="py-1">1.894</td>
          <td class="py-1 text-emerald-600 font-semibold">Excellent Purity (Target: ~1.8)</td>
        </tr>
        <tr>
          <td class="py-1"><strong>A260 / A230 Ratio</strong></td>
          <td class="py-1">1.785</td>
          <td class="py-1 text-emerald-600 font-semibold">Minimal salt/organic carryover</td>
        </tr>
      </tbody>
    </table>
  </div>
</div>

---

# Western Blot Trials (June 22)

<div class="grid grid-cols-12 gap-6 mt-4 items-center">
  <div class="col-span-5">
    <h2 class="text-teal-600 font-semibold mb-2">Total Protein Transfer Audit</h2>
    <p class="text-xs text-gray-700 leading-relaxed mb-3">
      Evaluated total protein transfer of ADAM10/17 periplasmic extracts onto nitrocellulose membranes.
    </p>
    <div class="p-2 bg-yellow-50 border border-yellow-200 rounded text-2xs text-yellow-800 mb-2">
      <strong>Key Issue:</strong> Target bands were run into the nitrocellulose divider sheet rather than the membrane itself.
    </div>
    <ul class="space-y-1 text-2xs text-gray-700 list-disc pl-4">
      <li>Resolved buffer dilution issue: running buffer labeled as 40x was actually 5x.</li>
      <li>Ponceau S staining performed before blocking to enable densitometry.</li>
      <li>Secondary antibody incubation and imaging scheduled for next steps.</li>
    </ul>
  </div>
  <div class="col-span-7 flex flex-col items-center border p-2 rounded bg-gray-50">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/20260624_Ponceau.jpg" class="h-52 object-contain" alt="Ponceau S Membrane Stain"/>
    <div class="text-2xs text-gray-600 mt-1"><strong>Ponceau S Stain (June 24):</strong> Verification of transfer stack lanes and ladder.</div>
  </div>
</div>

---

# FPLC Chromatogram Trace (June 25)

<div class="grid grid-cols-12 gap-6 mt-4 items-center">
  <div class="col-span-5">
    <h2 class="text-teal-600 font-semibold mb-2">AEC Purification of ADAM10</h2>
    <p class="text-xs text-gray-700 leading-relaxed mb-3">
      Purification of osmotic shock periplasmic extracts on AEC HiTrap Q FF 5mL column under salt gradient.
    </p>
    <div class="p-2 bg-emerald-50 border border-emerald-200 rounded text-2xs text-emerald-800 mb-2">
      <strong>Chromatogram Peak:</strong> Elution trace confirms target peak at low salt.
    </div>
    <ul class="space-y-1 text-2xs text-gray-700 list-disc pl-4">
      <li>Main UV peak reaches <strong>67 mAU</strong> at <strong>13.1 mL</strong>.</li>
      <li>Salt gradient %B (green curve) shows elution occurs at low conductivity.</li>
      <li>Enriched fractions pooled for yield optimization and downstream assays.</li>
    </ul>
  </div>
  <div class="col-span-7 flex flex-col items-center border p-2 rounded bg-gray-50">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/fplc_chromatogram_20260625.png" class="h-52 object-contain" alt="FPLC Purification Trace"/>
    <div class="text-2xs text-gray-600 mt-1"><strong>HiTrap Q AEC Trace (June 25):</strong> UV A280 (blue) and Salt Gradient %B (green).</div>
  </div>
</div>

---

# Western Blot Trials (June 26)

<div class="grid grid-cols-12 gap-6 mt-4 items-center">
  <div class="col-span-5">
    <h2 class="text-teal-600 font-semibold mb-2">Lysis Pellet and Supernatant Fractions</h2>
    <p class="text-xs text-gray-700 leading-relaxed mb-3">
      Evaluated total protein transfer of supernatant and lysis pellet fractions of ADAM10-pMopac expression trials.
    </p>
    <div class="p-2 bg-emerald-50 border border-emerald-200 rounded text-2xs text-emerald-800 mb-2">
      <strong>Transfer Verification:</strong> Successful protein transfer verified via Ponceau S staining.
    </div>
    <ul class="space-y-1 text-2xs text-gray-700 list-disc pl-4">
      <li>Ponceau S staining confirmed successful protein transfer and lane alignment for both pellet and supernatant membranes.</li>
      <li>Membranes are currently incubating in anti-FLAG primary antibody overnight at 4°C.</li>
      <li>Fluorescence detection is scheduled for tomorrow to verify soluble vs. insoluble target yields.</li>
    </ul>
  </div>
  <div class="col-span-7 grid grid-cols-2 gap-2 border p-2 rounded bg-gray-50">
    <div class="flex flex-col items-center">
      <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/20260626_Supernatant_Ponceau_labeled.png" class="h-44 object-contain" alt="Supernatant Ponceau"/>
      <div class="text-3xs text-gray-500 mt-1">Supernatant Ponceau Stain</div>
    </div>
    <div class="flex flex-col items-center">
      <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/20260626_Pellet_Ponceau_labeled.png" class="h-44 object-contain" alt="Pellet Ponceau"/>
      <div class="text-3xs text-gray-500 mt-1">Lysis Pellet Ponceau Stain</div>
    </div>
  </div>
</div>

---

# ADAM10 Silver-Stain Validation (June 26)

<div class="grid grid-cols-12 gap-6 mt-4 items-center">
  <div class="col-span-5">
    <h2 class="text-teal-600 font-semibold mb-2">30°C Induction Time Course</h2>
    <p class="text-xs text-gray-700 leading-relaxed mb-3">
      SDS-PAGE and high-sensitivity silver staining analysis of the 30°C ADAM10-pMopac expression trials.
    </p>
    <div class="p-2 bg-yellow-50 border border-yellow-200 rounded text-2xs text-yellow-800 mb-2">
      <strong>Expression Results:</strong> Strong target expression peaking at 20-22 hours.
    </div>
    <ul class="space-y-1 text-2xs text-gray-700 list-disc pl-4">
      <li>Visible target band at expected molecular weight (~60-70 kDa) starting around 3 hours post-induction.</li>
      <li>High abundance in lysis pellet confirms partitioning predominantly into insoluble inclusion bodies.</li>
      <li>Lysis supernatant shows faint band, suggesting presence of soluble fraction.</li>
    </ul>
  </div>
  <div class="col-span-7 flex flex-col items-center border p-2 rounded bg-gray-50">
    <img src="../../SharedAssets/figures/De_Novo_Binder_Generation/20260626_ADAM10_Silver_labeled.png" class="h-52 object-contain" alt="ADAM10 Silver Gel"/>
    <div class="text-2xs text-gray-600 mt-1"><strong>Silver-Stained Gel (June 26):</strong> Progressive 30°C induction timecourse.</div>
  </div>
</div>

---

# Summary & Next Steps

<div class="grid grid-cols-2 gap-8 mt-6">
  <div class="border-r pr-6">
    <h2 class="text-teal-600 font-bold mb-3 flex items-center gap-2">
      <span class="p-1 bg-teal-100 text-teal-800 rounded-full text-2xs">✓</span> Achieved Milestones
    </h2>
    <ul class="space-y-3 text-xs text-gray-700">
      <li><strong>Western Blot Ponceau (June 26):</strong> Confirmed transfer stack and lane alignment for pellet & supernatant samples.</li>
      <li><strong>Silver-Stained SDS-PAGE (June 26):</strong> Visualized 30°C induction of ADAM10 (target band at ~60-70 kDa peaking at 20-22 hr).</li>
      <li><strong>FPLC AEC Trace (June 25):</strong> Captured and analyzed elution profile of recombinant ADAM10 (peak at 13.1 mL).</li>
    </ul>
  </div>
  <div>
    <h2 class="text-teal-600 font-bold mb-3 flex items-center gap-2">
      <span class="p-1 bg-teal-100 text-teal-800 rounded-full text-2xs">→</span> Future Directions
    </h2>
    <ul class="space-y-3 text-xs text-gray-700">
      <li><strong>Western Blot Detection:</strong> Image anti-FLAG fluorophore-labeled membranes to determine soluble vs. insoluble yields.</li>
      <li><strong>ADAM17 Target Recovery:</strong> Re-isolate and sequence pMopac-ADAM17 after helper plasmid contamination was identified.</li>
      <li><strong>Functional Assays:</strong> Conduct Enzo MMP target assays to verify activity of purified ADAM10 fractions.</li>
    </ul>
  </div>
</div>
