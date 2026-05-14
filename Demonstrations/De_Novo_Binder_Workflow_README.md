# De Novo Binder Generation — Workflow README

> **Project Goal:** Engineer MMP9-selective TIMP3 loop variants using a generative AI pipeline  
> (RFdiffusion → ProteinMPNN → AlphaFold3 → Yeast Display Flow Cytometry).  
> Primary selectivity target: **MMP9 vs MMP2** (structurally homologous, clinically distinct).  
> ADAM17 is included in normalized binding comparisons to TIMP3 but not in selectivity analyses.  
> MMP3 data was used for internal calibration trials only and is excluded from publications.
> **Clinical Objective:** Selective MMP9 inhibition for oncology applications (tumor metastasis) without off-target musculoskeletal toxicity from MMP2 inhibition.
> **Strategic Roadmap:** Includes "Reverse Selectivity" (MMP2 > MMP9) as a controllability proof of concept, ADAM10 presentation optimization, and absolute kinetic measurement (SPR).

---

## Table of Contents
1. [Reference Data](#1-reference-data)
2. [Computational Generation](#2-computational-generation--stage-1-rfdiffusion)
3. [Sequence Design](#3-sequence-design--stage-2-proteinmpnn)
4. [Structural Verification](#4-structural-verification--stage-3-alphafold-server)
5. [Candidate Selection](#5-candidate-selection--stage-4-best-binder)
6. [Library Synthesis Validation](#6-library-synthesis-validation--stage-5-twistbioorder)
7. [Experimental Flow Cytometry Data](#7-experimental-flow-cytometry-data--stage-6)
8. [Aggregate Statistical Analysis](#8-aggregate-statistical-analysis)
9. [Dissemination Materials](#9-dissemination-materials)
10. [Target Summary](#target-experiment-status-summary)

---

## 1. Reference Data

| File / Folder | Description |
|---|---|
| `Data/3CKI_Sarmazdeh.pdb` | Native TIMP3 crystal structure (Sarmazdeh et al.) — scaffold reference |
| `Data/8FNS.cif` | TIMP3-MMP9 complex CIF structure |
| `Data/TIMP3_vs_ADAM17_X_ray.pdb` | TIMP3–ADAM17 complex X-ray structure |
| `Data/Target_Crystal_Structures/` | Crystal PDBs for all MMP/ADAM targets |
| `Data/TIMP_Complexes/HADDOCK_PDB/` | HADDOCK-docked TIMP3:target PDB files (primary ProteinMPNN input) |
| `Data/TIMP_Complexes/AlphaFold_PDB/` | AlphaFold-predicted TIMP3:target PDBs |
| `Data/TIMP_Complexes/AlphaFold_CIF/` | AlphaFold CIF outputs for structural visualization |
| `Data/HADDOCK_Outputs/` | Raw HADDOCK docking outputs |
| `Data/Original_Naive_library_Backbone.xlsx` | Original TIMP3 loop geometry measurements |
| `Data/Sequence_Catalog/` | Curated sequence reference library |
| `Protein-Analysis/TIMP3_variants_names.csv` | **Master construct list** with expected binding outcomes per variant |

---

## 2. Computational Generation — Stage 1: RFdiffusion

| File | Role | Output |
|---|---|---|
| `Generation/RFDiffusion_Gen.ipynb` | Interactive RFdiffusion runner | Poly-glycine backbone PDBs |
| `Generation/RFd_Batch.py` | Batch RFdiffusion over all loop types | Poly-glycine PDBs in `Generation/outputs/` |
| `Generation/RFd3_Gen.ipynb` | RF3-based generation notebook | — |
| `Generation/RFd3_batch.py` | RF3 batch runner | — |
| `Generation/clean_renum_struct.py` | Cleans and renumbers PDB outputs before PMPNN input | Processed PDBs |

**Input:** `Data/TIMP_Complexes/HADDOCK_PDB/` (HADDOCK-docked TIMP3:target complexes)  
**Output:** `Generation/outputs/` — poly-glycine PDB backbone proposals per loop/expansion

> Expansion lengths: 6–24 aa. Targets: AB loop (res. 31–36), C loop (63–68), EF loop (93–96).  
> Compute: A30 GPU, ~20 hr/run.

---

## 3. Sequence Design — Stage 2: ProteinMPNN

| File | Role | Output |
|---|---|---|
| `Generation/ProMPNN_Batch.py` | **Primary script** — batch ProteinMPNN over all loop/target PDBs | Top sequences + batch JSON for AF submission |
| `Generation/ProteinMPNN_Gen.ipynb` | Interactive single-run notebook | Single loop sequences |
| `AlphaFold/PMPNN_Out_Loop_Fetch.py` | Extracts loop sequences from PMPNN output files | Loop-only FASTA sequences |
| `AlphaFold/Top_loops_AF_gen.py` | Formats top-N sequences into AF Server batch JSON | `AlphaFold/`-ready JSON files |

**Input:** Poly-glycine PDBs from `Generation/outputs/`  
**Parameters:** 1000 seqs/target, T=0.2, AB/C/EF loop positions free, all scaffold residues fixed  
**Output:** Ranked sequence FASTA files + batch JSON for AlphaFold Server submission

---

## 4. Structural Verification — Stage 3: AlphaFold Server

| File | Role | Output |
|---|---|---|
| `AlphaFold/AlphaFoldServer.ipynb` | Batch job submission to AlphaFold Server | ZIP archives submitted to AF Server |
| `AlphaFold/AF_batch_gen.py` | Auto-generates AF Server JSON job files | JSON job files |
| `AlphaFold/AF_rename.py` | Renames AF Server download ZIPs consistently | Organized ZIP files |
| `AlphaFold/AF_swap_seqs.py` | Swaps sequences within AF JSON jobs | Modified JSON jobs |
| `AlphaFold/AlphaFoldServer_Analysis.ipynb` | **Primary parser** — extracts confidence metrics from AF JSON outputs | Summary tables + comparison matrices |
| `AlphaFold/AlphaFoldServer_MDA.ipynb` | Multi-dimensional AlphaFold analysis | Extended metric analyses |
| `AlphaFold/AlphaFoldServer_SA.ipynb` | Statistical analysis of AlphaFold results | Statistical summaries |
| `AlphaFold/ReadAlphaFold2-multimer.ipynb` | AF2-multimer result parser | — |

**Input:** Top-10 sequences per loop/target (from ProteinMPNN batch JSON)  
**Output folders:**  
- `Local/AlphaFoldRun/` — Raw AF Server ZIP archives per target/loop combination
  - Naming: `fold_timp3_variant_{target}_{loop}_wt.zip`
  - Targets present: `mmp2`, `mmp9`, `adam17`, `adam10` (ADAM10 unusable — low pos ctrl)
- `Local/AlphaFoldRun/analysis_output/` — Per-run parsed metrics
- `Local/AlphaFoldMetrics/` — Cross-run comparison matrices (`.xlsx` and `.tex`)
  - `AlphaFold_Summary.xlsx` / `.tex` — Master summary across all targets/loops
  - `Comparison_Matrix_ipTM.xlsx` — ipTM comparison matrix (T-score basis)
  - `Comparison_Matrix_Loop_pLDDT.xlsx` — Loop structural confidence
  - `Comparison_Matrix_Interface_PAE_(Loop).xlsx` — Interface precision metric
  - `Comparison_Matrix_Mean_pLDDT.xlsx`, `pTM.xlsx`, `Chain_A_pTM.xlsx`, `Chain_B_pTM.xlsx`

**Metrics extracted per complex:** pTM, ipTM, Mean pLDDT, Loop pLDDT, Interface PAE (Loop)

---

## 5. Candidate Selection — Stage 4: Best Binder

| File | Role | Output |
|---|---|---|
| `AlphaFold/AlphaFold_Best_Binder.ipynb` | **Selection engine** — applies T-score filter to AF metrics | Final 13 construct candidate list |

**Input:** `Local/AlphaFoldMetrics/` comparison matrices  
**Selection criteria:**
- T-score > 2.0 on ≥1 target across ≥2 independent metrics (multi-metric consensus)
- Global pTM ≥ 0.80 and ipTM ≥ 0.82 for intended target
- Predicted selectivity aligned with design intent (MMP9 > MMP2 for MMP9-targeted loops)

**Output:** 13 final constructs — 7 AB-loop, 6 C-loop, 3 dual ABC-combination variants

---

## 6. Library Synthesis Validation — Stage 5: TwistBioOrder

| File | Role | Output |
|---|---|---|
| `Generation/TwistBioOrder.ipynb` | RE site validation + codon optimization (Saccharomyces cerevisiae) | Ordered library files |

**Restriction sites screened:** BsrGI (TGTACA), BamHI (GGATCC), BsaI (GGTCTC)  
**Result:** 13/13 constructs passed (0 failures)  
**Output files (in `Local/Twist_Order_Dec2025/`):**

| File | Description |
|---|---|
| `twist_library.csv` | Final construct DNA + amino acid sequences |
| `twist_library.fasta` | FASTA format of ordered sequences |
| `twist_library_w_Combo.csv` | Library with dual-loop combination constructs included |
| `validation_report.txt` | Per-construct RE site validation results (all PASS) |

**Order placed:** December 2025, Twist Bioscience  
**Format:** Linear dsDNA with Golden Gate assembly overhangs

---

## 7. Experimental Flow Cytometry Data — Stage 6

### Raw Data Organization (all in `Local/`)

Each experiment follows the naming convention:  
`{TARGET}_{YYYYMMDD}[_variant]/` → `{TARGET}_{YYYYMMDD}[_variant]_Renamed/` → `{TARGET}_{YYYYMMDD}[_variant]_Renamed_Analysis/`

| Target | Trial Dates | Notes |
|---|---|---|
| **MMP9** | 20260210, 20260227, 20260304, 20260310, 20260319, 20260331, 20260424 (A+B), 20260509 (Enzo, Sino, SinoThaw) | Primary selectivity target |
| **MMP2** | 20260424 (A+B), 20260509 (Arly-FLAG, Enzo, Sino, SinoThaw) | Primary comparison target |
| **ADAM17** | 20260403, 20260407, 20260416, 20260509 (Enzo, Sam, Sino, SinoThaw) | Secondary target — included in normalized TIMP3 comparisons |
| **ADAM10** | 20260403, 20260407, 20260509 (Sam-FLAG) | ⚠️ Excluded from analysis — pos ctrl Double+% below 3% threshold in all runs |
| **MMP3** | 20260127, 20260304, 20260330, 20260331, 20260407 (G), 20260416, 20260424 (A+B) | ⚠️ Internal calibration trials only — excluded from publications |
| **NoTarget** | 20260418 | Negative control (no target antigen) |

### Per-Trial Analysis Scripts

| Script | Role | Output per trial |
|---|---|---|
| `Protein-Analysis/rename_fcs.py` | Renames raw FCS files to standard naming convention | `*_Renamed/` folder |
| `Protein-Analysis/analyze_fcs.py` | **Primary single-trial FCS analysis** — gating, metrics, plots | `*_Renamed_Analysis/` folder |
| `Protein-Analysis/analyze_fcs_target_trial.py` | Target-trial-level analysis for cross-date summaries | Per-target summary files |
| `Protein-Analysis/fcs_viewer.py` | Interactive Streamlit FCS viewer dashboard | Web UI |
| `Protein-Analysis/run_pipeline.py` | Runs full rename → analyze pipeline | Analysis outputs |

---

## 8. Aggregate Statistical Analysis

| Script | Role | Output |
|---|---|---|
| `Protein-Analysis/aggregate_analysis.py` | **Primary aggregate script** — pools all trials, computes metrics, ANOVA, Tukey-HSD | All outputs in `Local/Aggregate_FCS_Analysis/` |

### Output Structure (`Local/Aggregate_FCS_Analysis/`)

| Subfolder / File | Metric | Valid for Cross-Target? |
|---|---|---|
| `Binding_Efficiency/` | DP/FITC+ fraction | ✅ **YES — use for MMP9 vs MMP2 selectivity** |
| `Norm_Median_Ratio/` | PE-MFI normalized to TIMP3-WT | ❌ Within-target only (MMP9 vs TIMP3, MMP2 vs TIMP3, ADAM17 vs TIMP3) |
| `Norm_Bind_Med_Expr_Positive/` | Binding MFI in expression-positive gate, normalized | ❌ Within-target only |
| `Norm_IWB_Index/` | Integrated weighted binding index | ❌ Within-target only |
| `Selectivity_Analysis/` | Binding efficiency-based selectivity ratios (MMP9/MMP2) | ✅ Cross-target ratio |
| `aggregate_summary.csv` | Full pooled per-trial summary | — |
| `cross_target_summary.csv` | Per-construct cross-target binding efficiency summary | — |
| `significant_tukey_summary.txt` | All constructs with ANOVA p<0.05 + Tukey-HSD pairs | — |

### Per-Metric Subfolders contain:
- `anova_results_{TARGET}.txt` — one-way ANOVA results across all constructs for that target
- `plot_{TARGET}.png` — aggregate comparison bar chart
- `aggregate_colorcoded.png` — color-coded heatmap across all constructs

### Key Statistical Results (from `significant_tukey_summary.txt`)

| Construct | Expected | ANOVA p | MMP9 Eff. | MMP2 Eff. | Tukey MMP2 vs MMP9 |
|---|---|---|---|---|---|
| **C 12** | M9+, M9>M2 | 0.0115 | 91.3% | 32.1% | p=0.015 ✅ |
| **C 15** | M9+, M9>M2 | 0.0177 | 88.6% | 33.8% | p=0.043 ✅ |
| **AB 6** | M9>M2 | 0.0002 | 89.9% | 23.8% | p=0.0004 ✅ |
| C 13 | Low | 0.096 | — | — | n.s. ✓ expected |
| AB 5 | Low | 0.601 | — | — | n.s. ✓ expected |

> **ADAM17:** Included in normalized TIMP3 comparisons (`Norm_Median_Ratio/`, within-target).  
> Not used for selectivity since ADAM10 data is unavailable.

---

## 9. Dissemination Materials

### Paper
| File | Description |
|---|---|
| `Demonstrations/Papers/De_Novo_Binder_Generation/De_Novo_Binder_Generation.tex` | LaTeX source |
| `Demonstrations/Papers/De_Novo_Binder_Generation/De_Novo_Binder_Generation.pdf` | Compiled PDF (17 pages) |

### Presentation (Slidev)
| File | Description |
|---|---|
| `Demonstrations/Presentations/De_Novo_Binder_Generation/slides.md` | Slidev Markdown source |

### AMA Abstract Poster
| File | Description |
|---|---|
| `Demonstrations/Presentations/AMA_Abstract_Poster/poster_data.yaml` | Content data (edit this to update content) |
| `Demonstrations/Presentations/AMA_Abstract_Poster/template.html` | Jinja2 HTML template |
| `Demonstrations/Presentations/AMA_Abstract_Poster/style.css` | Poster CSS (dimensions defined here) |
| `Demonstrations/Presentations/AMA_Abstract_Poster/build_poster.py` | Build script → renders output.html + PDF |
| `Demonstrations/Presentations/AMA_Abstract_Poster/output.html` | ✅ Rendered poster (open in browser to view) |
| `Demonstrations/Presentations/AMA_Abstract_Poster/generate_poster_figures.py` | Generates all publication figures |

> **⚠️ POSTER DIMENSIONS: 48" W × 36" H (landscape), aspect ratio 4:3 (≈ 19:6 approximate)**  
> Enforced in `style.css` (`.poster-container`) and `build_poster.py` (Playwright PDF export).  
> Do NOT change these dimensions for AMA submission.

### Shared Figure Assets
| Folder | Description |
|---|---|
| `Demonstrations/SharedAssets/figures/De_Novo_Binder_Generation/` | All poster/paper figures |
| → `fig_mmp9_vs_mmp2.png` | MMP9 vs MMP2 binding efficiency bar chart (ANOVA annotated) |
| → `fig_specificity_ratio.png` | MMP9/MMP2 selectivity ratio per construct |
| → `fig_binding_heatmap.png` | Binding efficiency heatmap (all constructs × targets) |
| → `fig_pipeline.png` | 6-stage pipeline overview diagram |
| `Demonstrations/SharedAssets/figures/UNR Logo.png` | University logo |

---

## Target Experiment Status Summary

| Target | Used In | Notes |
|---|---|---|
| **MMP9** | All documents — primary selectivity target | Full analysis, ANOVA confirmed |
| **MMP2** | All documents — primary comparison target | Full analysis, ANOVA confirmed |
| **ADAM17** | Within-target TIMP3 normalized comparisons only | Insufficient ADAM10 to show ADAM selectivity |
| ADAM10 | ❌ Excluded | Pos ctrl Double+% < 3% in all runs |
| MMP3 | ❌ Excluded from publications | Internal calibration only (Jan–Apr 2026) |

---

## Workflow Quick-Reference

```
Data/TIMP_Complexes/HADDOCK_PDB/          ← HADDOCK templates (ProteinMPNN input)
    │
    ▼ Generation/RFd_Batch.py (optional for expansions)
Generation/outputs/                        ← Poly-glycine backbones
    │
    ▼ Generation/ProMPNN_Batch.py
AlphaFold/Top_loops_AF_gen.py              ← Batch JSON for AF Server
    │
    ▼ AlphaFold/AlphaFoldServer.ipynb (submission)
Local/AlphaFoldRun/                        ← Raw AF Server ZIP downloads
    │
    ▼ AlphaFold/AlphaFoldServer_Analysis.ipynb
Local/AlphaFoldMetrics/                    ← Comparison matrices (ipTM, pLDDT, PAE)
    │
    ▼ AlphaFold/AlphaFold_Best_Binder.ipynb (T-score > 2.0 filter)
    │  → 13 final constructs selected
    │
    ▼ Generation/TwistBioOrder.ipynb (RE validation)
Local/Twist_Order_Dec2025/                 ← Final ordered library
    │
    ▼ [Lab: yeast transformation + flow cytometry]
Local/{TARGET}_{DATE}/                     ← Raw FCS data
Local/{TARGET}_{DATE}_Renamed/             ← Standardized filenames
Local/{TARGET}_{DATE}_Renamed_Analysis/    ← Per-trial metrics
    │
    ▼ Protein-Analysis/aggregate_analysis.py
Local/Aggregate_FCS_Analysis/              ← Pooled metrics, ANOVA, Tukey-HSD
    │
    ▼ Demonstrations/Presentations/AMA_Abstract_Poster/generate_poster_figures.py
Demonstrations/SharedAssets/figures/       ← Publication figures
    │
    ▼ Demonstrations/Presentations/AMA_Abstract_Poster/build_poster.py
Demonstrations/Presentations/AMA_Abstract_Poster/output.html  ← Final poster
Demonstrations/Papers/De_Novo_Binder_Generation/De_Novo_Binder_Generation.pdf
```
