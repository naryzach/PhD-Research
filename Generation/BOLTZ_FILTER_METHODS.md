# In-Silico Pre-Filtering for AF3: Methods & Decision Log

**Project:** Iterative TIMP3-scaffold binder design (`Generation/iterative_refinement.py`)
**Question this document answers:** *Why did we add Boltz-2, and how do we know whether it's worth the compute?*
**Last updated:** 2026-06-01

---

## 1. The constraint that drives everything

The AlphaFold3 web server is capped at **~30 jobs/day**. We generate far more designs
than that per iteration (currently 80/iteration across 4 targets). We therefore need an
**in-silico pre-filter** that runs locally, with no daily cap, to decide which ~30
designs are worth an AF3 slot. The whole point is **AF3 budget efficiency** — spend the
scarce gold-standard validation on designs most likely to succeed.

"Success" is operationalized as **AF3 interface pTM (ipTM) ≥ 0.55** for the
binder–target complex. ipTM is AF3's confidence that the two chains are docked in the
correct relative orientation; it is the field that best tracks real binding for de novo
binders.

---

## 2. Pipeline architecture

```
RFd3 (backbone)
   → LigandMPNN (sequence design on fixed scaffold, variable loops AB/C/EF)
   → RF3 complex prediction        ← cheap, single-sequence; geometric features only
   → Boltz-2 complex prediction    ← local AF3-class ranker (the pre-filter under test)
   → Hall of Fame (per-target, composite-ranked)
   → AF3 Server (top ~30, the gold-standard gate)
   → AF3 results imported back, overwrite composite with high-trust metrics
```

Temperature anneals 0.50 → 0.10 (LigandMPNN sampling) to move from exploration to
exploitation; loop-length ranges adaptively narrow toward HOF-successful lengths from
iteration 3 on.

### Source-aware composite score (`calc_composite`)

Designs are ranked by a composite that **depends on the highest-trust prediction
available**, because not all predictors are equally trustworthy (see §4):

| Source | Formula weights | Rationale |
|--------|-----------------|-----------|
| **AF3** (validated) | ipTM 0.45, pLDDT 0.25, pTM 0.15, PAE 0.10, RMSD 0.05 | full trust |
| **Boltz-2** | ipTM 0.45, pLDDT 0.25, pTM 0.15, PAE 0.15 | MSA/template-based, AF3-class |
| **RF3 only** | n_contacts 0.40, PAE 0.30, RMSD 0.30 — **capped at 0.50** | geometric only; RF3 confidence excluded |

RF3's ipTM/pTM/pLDDT are **deliberately excluded** from ranking (see §4.1). The RF3-only
branch is capped so any AF3- or Boltz-scored design out-ranks an unvalidated one.

---

## 3. Statistical methods

**Data sources.**
- *Local metrics:* every `Local/iterative_refinement/it_*/round_summary.csv`
  (one row per design, with RF3 and Boltz metrics).
- *AF3 ground truth:* the AF3 Server results ZIP. For each job we read
  `*_summary_confidences_0.json` (ipTM, pTM) and `*_full_data_0.json`
  (per-atom pLDDT → binder-chain mean; per-token PAE → mean of the two cross-chain
  blocks = interface PAE).

**Matching.** AF3 jobs are joined to local designs by **exact binder sequence**
(`full_seq`). This is unambiguous — sequences are unique per design.

**Statistics.**
- Pearson *r* (linear association) and Spearman *ρ* (rank association) between each
  local metric and AF3 ipTM.
- "Hit-rate": fraction of designs with AF3 ipTM ≥ 0.55.
- For ordinary batches, top-K-by-Boltz hit-rate vs. the random-selection expectation.

**Reproduce any time:** `python Generation/validate_boltz_filter.py <af3_results.zip>`
(auto-detects stratified vs. correlation mode).

---

## 4. Findings to date

### 4.1 RF3 single-sequence confidence is anti-correlated with AF3 (batch 1, n=18)

The first AF3 batch (RF3-only ranking era) showed RF3's own confidence metrics were
**worse than useless** as a ranker:

| RF3 metric vs AF3 ipTM | Pearson r | Spearman ρ |
|---|---:|---:|
| RF3 ipTM | −0.07 | −0.31 |
| RF3 pTM | −0.63 | −0.54 |
| RF3 pLDDT | −0.26 | −0.21 |

**Cause:** RF3 here runs single-sequence (no MSA, no template). For de novo binders with
no evolutionary signal, single-sequence confidence does not reflect true interface
quality. **Consequence:** RF3 confidence was removed from `calc_composite`; RF3 is kept
only for geometric features (interface contacts, backbone RMSD) and as a "did anything
fold at all" sanity check.

### 4.2 Boltz-2 batch (n=28): high absolute quality, but ranking signal unproven

Second AF3 batch, after Boltz-2 was added and used to rank:

| Metric vs AF3 ipTM | Pearson r | Spearman ρ |
|---|---:|---:|
| **Boltz ipTM** | **−0.22** | **−0.06** |
| Boltz pLDDT | +0.29 | **+0.34** ← best signal |
| Boltz pTM | −0.19 | −0.08 |
| RF3 ipTM | +0.03 | −0.16 |

**Absolute quality (the good news):**

| Batch | Selection method | AF3 ipTM ≥ 0.55 | mean AF3 ipTM |
|---|---|---:|---:|
| Batch 1 | RF3 composite (≈ random; §4.1) | 7/18 = 39% | 0.38 |
| Batch 2 | Boltz composite (top band) | **22/28 = 79%** | **0.65** |

Best designs reached AF3 ipTM 0.87 (ADAM10), 0.81 (ADAM17), 0.79 (MMP2, MMP9).

### 4.3 Why the near-zero Boltz correlation is NOT a clean verdict

Two confounds prevent concluding "Boltz doesn't work" from §4.2:

1. **Range restriction.** Boltz scores the full design pool across 0.15–0.88
   (std 0.18), but the composite-ranked batch sent **only the top ~18%**
   (Boltz ipTM 0.68–0.88) to AF3. Correlation measured on a truncated top-slice is
   mathematically attenuated. We have **no AF3 scores for the designs Boltz rejected**,
   so we cannot see whether it correctly filtered out bad ones.

2. **Statistical power.** At n=28 with restricted range, the critical |r| for p<0.05
   is ≈ 0.37. **None** of the §4.2 correlations are statistically significant. The only
   robust facts are the absolute hit-rate (79%) and that it beat the effectively-random
   batch 1 (39%).

The 39% → 79% jump is **circumstantial evidence that Boltz coarse-filters**
(separates good from bad) even though it does not **fine-rank** (good from better).
But because batch 1 and batch 2 are different runs at different iterations, this is
suggestive, not conclusive.

---

## 5. The decision: a stratified validation experiment

To break the range-restriction confound we need AF3 scores at **low, mid, and high**
Boltz ipTM, not just the top. So the next AF3 batch is deliberately stratified:

```bash
python Generation/iterative_refinement.py --targets MMP2 MMP9 ADAM10 ADAM17 \
    --stratified-export 30
```

This pools all designs from every `round_summary.csv`, bins each target's designs into
LO/MID/HI Boltz-ipTM bands (equal frequency), and samples evenly — writing
`af3_submission_stratified.json` plus `stratified_manifest.json` (records each job's
band). Verified composition on current data (24 jobs, balanced 6/target):

| Band | n | Boltz ipTM (min–mean–max) |
|------|---|---------------------------|
| LO | 8 | 0.13 – 0.32 – 0.62 |
| MID | 8 | 0.38 – 0.56 – 0.72 |
| HI | 8 | 0.53 – 0.75 – 0.89 |

After running these through AF3:

```bash
python Generation/validate_boltz_filter.py <af3_results.zip>
```

**Interpretation:**
- **AF3 ipTM rises LO→MID→HI (Δ ≥ 0.10):** Boltz is a useful coarse filter — keep it as
  the AF3 gate, continue current strategy.
- **Flat across bands (Δ < 0.05):** Boltz is not separating good from bad. Its compute
  cost is not buying filtering → drop it, or replace with the §4.2 best signal
  (Boltz pLDDT) or another model (§6).

### 5.1 Stratified results (24 designs, run 2026-06-01)

| Boltz band | n | Boltz ipTM | AF3 ipTM mean | AF3 ≥ 0.55 |
|------------|---|-----------:|--------------:|-----------:|
| LO  | 8 | 0.316 | 0.479 | 50% |
| MID | 8 | 0.559 | 0.379 | 38% |
| HI  | 8 | 0.745 | **0.599** | **62%** |

- HI vs rest contrast: **+0.170** AF3 ipTM. Full-range Spearman **+0.11** (weak, positive).
- **Not monotonic** — MID dips below LO. Per-target breakdown explains why:

  | Target | LO | MID | HI | behavior |
  |--------|----|----|----|----------|
  | MMP2   | 0.20 | 0.44 | 0.56 | clean rise |
  | MMP9   | 0.42 | 0.24 | 0.62 | top works |
  | ADAM17 | 0.64 | 0.16 | 0.67 | bimodal |
  | ADAM10 | 0.64 | 0.67 | 0.54 | **inverted** |

- Precision down the Boltz ranking: top-4 = 75% hits; top-6/8/12 = 50% (= pool baseline).

**Verdict.** Boltz is a **weak, top-concentrated, target-inconsistent** filter. It enriches
only the extreme top (~top 15%), is no better than random past that, has poor recall
(good binders sit in LO/MID — e.g. ADAM17_LO 0.74, ADAM10_LO 0.70), and is *anti*-correlated
for ADAM10. Usable as a soft "send the very top per target" heuristic, but **not** the
reliable AF3 surrogate we wanted. This motivates trialing ESMFold2 (§6) on the same 24.

### 5.2 ESMFold2 head-to-head results (same 24 designs, 2026-06-01)

ESMFold2 (ESMC-6B base, ~25 GB weights, MSA-free) scored the same 24 designs on the
cluster; compared against AF3 ipTM ground truth via `validate_boltz_filter.py`:

| Predictor | Pearson | Spearman | p | top-6 AF3 hit-rate |
|-----------|--------:|---------:|--:|-------------------:|
| Boltz ipTM | +0.148 | +0.109 | 0.61 | 50% (3/6) |
| ESMFold2 ipTM | +0.250 | +0.283 | 0.18 | 83% (5/6) |
| **ESMFold2 pLDDT** | **+0.340** | **+0.335** | 0.11 | **100% (6/6)** |
| ESMFold2 pLDDT×ipTM (rank-prod) | — | +0.346 | 0.10 | — |

**ESMFold2 clearly beats Boltz** on both correlation (ρ 0.34 vs 0.11) and precision
(top-6 100% vs 50%; top-8 88% vs 50%). The best single signal is **ESMFold2 binder
pLDDT** — consistent with §4.2, where Boltz's own best signal was also pLDDT. The most
transferable predictor of AF3 success for these de novo binders is *"does the binder fold
confidently in the complex context,"* not either model's interface ipTM.

**Honest caveats (do not over-trust yet):**
1. **n=24; nothing is significant.** Critical |ρ|≈0.34 at p<0.05; ESMFold2 pLDDT lands
   right at the edge (p=0.11). The top-K hit-rates are the trustworthy practical signal,
   not the p-values.
2. **Per-target signal is noisy (n=6 each) and inconsistent across models:**

   | Target | Boltz ρ | ESMFold2-pLDDT ρ |
   |--------|--------:|-----------------:|
   | MMP2   | +0.72 | **−0.84** (inverted!) |
   | MMP9   | +0.43 | +0.71 |
   | ADAM10 | −0.49 (Boltz inverted) | +0.66 (ESMFold2 fixes it) |
   | ADAM17 | +0.14 | +0.49 |

   ESMFold2 pLDDT **fixes the ADAM10 inversion** that broke Boltz, and improves MMP9/ADAM17
   — but looks **inverted for MMP2** (n=6, so fragile). Neither model is uniformly right.
3. **pLDDT magnitude is target-dependent.** Top-8 by ESMFold2 pLDDT contained 0 MMP2
   designs (ADAM10/ADAM17/MMP9 get systematically higher pLDDT). A *global* pLDDT ranker
   would starve MMP2 of AF3 slots → **selection must be per-target top-K**, as the AF3
   export already does.

**Verdict (per the pre-registered rule in the prior §5.2):** ESMFold2 beat Boltz on both
Spearman and top-K, so **adopt ESMFold2 (rank by binder pLDDT, per-target) as the pre-AF3
filter** — but confirm on one more batch before fully trusting, since n=24 is
underpowered and MMP2 is a flag to watch.

**Recommended next step.** Score the full design pool with ESMFold2
(`score_with_esmfold2.py --all`, hours on the cluster), select the next 30-job AF3 batch
as per-target top-K by ESMFold2 binder pLDDT, and re-run `validate_boltz_filter.py`. If
the batch hit-rate is ~80%+ (and MMP2 doesn't crater), the filter is real and confirmed at
n≈48; then wire ESMFold2 into the pipeline as the ranking stage (subprocess to the
`esmfold2` env, mirroring Boltz).

---

## 6. Alternative local predictors considered

| Tool | MSA? | Complex ipTM? | Footprint | Verdict |
|------|------|---------------|-----------|---------|
| **RF3** (in foundry) | no (MSA-*starved*) | yes | installed | Anti-correlated (§4.1); geometric features only. |
| **Boltz-2** | yes (ColabFold API) | yes | ~5 GB | Current pre-filter; under validation (§5). Affinity head is protein-ligand only, not usable here. |
| **ESMFold2** (Biohub/Rives, 2026) | **no (MSA-free *by design*)** | **yes — pLDDT/pTM/ipTM** | large (ESMC-6B base; VRAM TBD on V100) | **Strong candidate.** Language-model-based; explicitly best on no-MSA cases (antibodies) — the closest analog to our de novo binders. MIT, open weights on HuggingFace. See note below. |
| **Chai-1** | yes (MSA) | yes | ~5 GB | MSA-based Boltz/AF3 alternative; outputs ipTM. Secondary option if a MSA-based second opinion is wanted. |

**On ESMFold2 (corrected 2026-06-01).** Earlier this row dismissed "ESMFold" as
single-sequence and therefore likely to share RF3's de novo blindness. **That was wrong
and conflated two different things.** RF3 fails here because it is a *structure* model run
**without the MSA it was trained to expect** (MSA-*starved*). ESMFold2 (Chan Zuckerberg
Biohub, Alex Rives's team) is built on the ESMC protein **language model** (trained on
billions of sequences) and is MSA-free **by design** — the evolutionary signal is baked
into the weights, not supplied per-query. Crucially, Biohub reports it is *strongest* on
targets that lack MSAs, such as antibodies — which is the closest published analog to our
engineered TIMP3-loop binders (no natural homologs → no informative MSA). It predicts
protein–protein and antibody–antigen complexes and **emits pLDDT, pTM, and ipTM** — the
exact metrics our pipeline consumes. Open source (MIT), weights on HuggingFace,
`pip install esm@git+https://github.com/Biohub/esm.git`.

This makes ESMFold2 a **better-motivated pre-filter than Boltz for this specific problem**:
Boltz's ColabFold MSA is informative for the protease target but near-empty for the
designed binder, so Boltz is partially MSA-starved on the side we care about; ESMFold2
sidesteps that entirely and skips the ~26 s/design MSA fetch.

Open practical questions before adopting:
1. **VRAM on V100.** ESMC-6B is ~12 GB of weights at bf16 before activations; may not fit
   a 16 GB V100 (tight on 32 GB). Must test on the actual node; a smaller ESMC variant may
   be available for the folding trunk.
2. **Independent benchmarks.** "Rivals/surpasses AF3" is from Biohub's own release + Nature
   news, not yet independent peer review on *de novo designed* binders specifically.
3. **Env isolation.** Like Boltz, install in a separate conda env and invoke via subprocess
   to avoid the numpy/cublas/lightning conflicts that hit the foundry env.

**Zero-cost evaluation path:** the 24 stratified designs (§5) become a labeled benchmark
once AF3 returns. Run ESMFold2 on those same 24 and compare Boltz-vs-AF3 and
ESMFold2-vs-AF3 on identical sequences — a clean three-way comparison that spends **no
extra AF3 budget**. Keep whichever local model best separates AF3 hits from non-hits.

Sources: Nature news (d41586-026-01686-3); Biohub release (biohub.org/news/world-model-of-protein-biology);
Modal ESMFold2 example (modal.com/docs/examples/esmfold2); Latent Space interview (latent.space/p/esmfold2).

---

## 7. Change log (key hyperparameters & fixes)

- Data source switched to `HADDOCK_Outputs/` (chains flipped vs. old `HADDOCK_PDB/`):
  binder = chain B in source PDB, but RFd3 **output** convention is binder = chain A.
  `DESIGN_BINDER_CHAIN="A"` / `DESIGN_TARGET_CHAIN="B"` used for all post-RFd3 reads.
- pLDDT normalized to 0–100 everywhere (`_normalize_plddt`).
- `BACKBONES_PER_TARGET` 10 → 20; `AF3_EXPORT_EVERY_N` 5 → 2; `AF3_TOP_N` → 30 (daily cap).
- AF3 export uses equal per-target quota (was top-N pooled → one target monopolized slots).
- RF3 template input attempted (AF3-style `useStructureTemplate`) but **rejected by this
  RF3 build** → reverted to single-sequence RF3. Correct template API still TBD.
- `promising` flag now uses Boltz ipTM ≥ 0.55 (was RF3 ipTM, which never crosses ~0.21).
- Boltz-2 runs via subprocess (`BOLTZ_EXECUTABLE`); 5-min timeout; affinity head omitted
  (protein-ligand only).
