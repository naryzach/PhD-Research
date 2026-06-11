# Demonstrations — AI Authoring Guide

> **This document is the authoritative source of truth for all research demonstration threads.**
> If you are an AI agent working on this project, read this file first. It defines what threads exist, how Papers and Presentations are structured, where data comes from, and how to maintain both outputs simultaneously.

> **For humans:** See [SETUP.md](SETUP.md) for step-by-step instructions on how to install dependencies, run presentations, build papers, and create Vue components — no prior Slidev or LaTeX experience needed.

---

## Table of Contents

1. [Role & Responsibilities](#role--responsibilities)
2. [Thread Registry](#thread-registry)
3. [Directory Architecture](#directory-architecture)
4. [Architectural Conventions](#architectural-conventions)
5. [Data-to-Demonstration Pipeline](#data-to-demonstration-pipeline)
6. [Thread Management Protocol](#thread-management-protocol)
7. [LaTeX Paper Conventions](#latex-paper-conventions)
8. [Slidev Presentation Conventions](#slidev-presentation-conventions)
9. [SharedAssets Conventions](#sharedassets-conventions)
10. [Codebase Map](#codebase-map)

---

## Role & Responsibilities

You are the **author and maintainer** of all materials inside `Demonstrations/`. Your responsibilities are:

1. **Papers** (LaTeX → PDF): Draft, expand, and refine academic manuscripts for each research thread.
2. **Presentations** (Slidev/Vue → browser slides): Create and maintain slide decks that visually communicate each thread's key findings.
3. **SharedAssets**: Produce, clean, and organize data exports, figures, and reusable components that are consumed by both Papers and Presentations.
4. **Thread Registry**: Keep this README's thread registry up to date as threads are added, modified, or deprecated.

**Critical rule:** This README is the **single source of truth** for which threads exist. Both `Papers/` and `Presentations/` are derived from the thread definitions here. If a thread is not listed in the registry below, it should not have a folder in either `Papers/` or `Presentations/`.

---

## Data Availability & Local Storage

**Understanding where data lives is critical to your work.** Research data is distributed across several locations, and not all of it is available at all times.

### Where Data Can Be Found

Data relevant to a thread may exist in **any** of the following locations:

| Priority | Location | Availability | Contents |
|----------|----------|-------------|----------|
| **Primary** | `Local/` | ⚠️ Lab computer only | The vast majority of experimental outputs — AlphaFold results, flow cytometry data, RFdiffusion outputs, ProteinMPNN sequences, pipeline intermediates, gating strategy PNGs, aggregate analysis plots |
| **Secondary** | `Data/` | ✅ In repository | Reference structures (PDBs, CIFs), sequence catalogs, crystal structure coordinates, binding result spreadsheets |
| **Secondary** | `Results/` | ✅ In repository | Curated summary files (e.g., `best_loops_per_length_*.csv`) |
| **Tertiary** | Thread-specific dirs | ✅ In repository | Code-adjacent data within `AlphaFold/`, `Generation/`, `Protein-Analysis/`, `MetalBinder/`, `TIMP-Dashboard/`, `ESM2 Classification/`, `FASTQ/`, `DNA Analysis/` |

### Key Facts About `Local/`

1. **`Local/` is gitignored.** It contains too much data to store in GitHub. It exists only on the **lab computer**.
2. **Data in `Local/` may come and go.** Due to storage constraints on the lab computer, datasets are sometimes archived or moved off the machine. A directory that existed last week may not be there today.
3. **`Local/` is the primary source for SharedAssets.** Most figures and cleaned data that feed into Papers and Presentations originate from `Local/` subdirectories.
4. **This machine (personal computer) may not have `Local/` populated at all.** When working on this computer, `Local/` may be empty or contain only a subset of data.

### What To Do When Data Is Missing

**It is acceptable — and encouraged — to ask the user if a specific piece of data can be added to `Local/`.** For example:

- *"I need the gating strategy PNGs from `Local/ADAM17_20260407_Renamed_Analysis/` to populate the flow cytometry slide. Can you make sure those files are available in `Local/` on this machine?"*
- *"The aggregate analysis CSVs at `Local/Aggregate_FCS_Analysis/` would allow me to generate the binding comparison figure. Are those available, or can they be copied over?"*
- *"I found a reference to `Local/lanm_output/global_sequence_catalog.csv` in the MetalBinder thread. Is that dataset currently on this machine?"*

Do not silently skip data or generate placeholder content when real data might be obtainable. **Always ask first.**

When data **is** available, clean and export it to `SharedAssets/` so that both Papers and Presentations can reference a stable, committed copy rather than depending on the volatile `Local/` directory.

---

## Thread Registry

Every research thread is listed below. Each thread MUST have a corresponding folder in both `Papers/<ThreadID>/` and `Presentations/<ThreadID>/`.

### Thread 1: `AMA_Abstract`

| Property | Value |
|----------|-------|
| **Title** | Generative Deep Learning Pipeline for Patient-Specific Therapeutic Binders Against Customized Cancer Markers |
| **Type** | Conference Abstract |
| **Paper Status** | ✅ Complete |
| **Presentation Status** | 🟡 Scaffold |
| **Description** | A structured AMA-format abstract summarizing the full generative pipeline. Covers the motivation for precision medicine binder design, the RFdiffusion → ProteinMPNN → AlphaFold methodology, Z-score-based consensus filtering, and the downstream implications for rapid therapeutic candidate generation. |
| **Codebase Sources** | General overview — references all pipeline components |
| **Key Data Paths** | None (self-contained abstract) |

---

### Thread 2: `De_Novo_Binder_Generation`

| Property | Value |
|----------|-------|
| **Title** | Experimental Pipeline: De Novo Protein Prediction to Binding Evaluation |
| **Type** | Full Research Paper |
| **Paper Status** | 🟢 Advanced Draft (figures, tables, 116 lines) |
| **Presentation Status** | 🟡 Scaffold |
| **Description** | The primary experimental thread. Details the bottom-up generation of novel TIMP3 variants targeting MMP2, MMP9, ADAM10, and ADAM17. Covers RFdiffusion backbone generation (AB/C/EF loop expansions), ProteinMPNN sequence optimization, AlphaFold conformation validation, Z-score standardization of loop confidence metrics, consensus filtering into affinity/specificity/selectivity categories, and downstream yeast display flow cytometry validation. |
| **Codebase Sources** | `Generation/` (RFd_Batch.py, ProMPNN_Batch.py, RFDiffusion_Gen.ipynb), `AlphaFold/` (AlphaFoldServer_Analysis.ipynb, Top_loops_AF_gen.py), `Protein-Analysis/` (analyze_fcs.py, aggregate_analysis.py, fcs_viewer.py), `FASTA Mod/` (make_fasta.py) |
| **Key Data Paths** | `Local/Aggregate_FCS_Analysis/` (flow cytometry aggregates), `Local/ADAM17_*_Analysis/` (gating strategy PNGs), `Local/proteinmpnn_output/`, `Local/rfdiffusion_output/`, `Results/best_loops_per_length_*.csv` |
| **Key Figures** | `Gating_Strategy_NegCtrl_Quad.png`, `aggregate_colorcoded.png` |

---

### Thread 3: `ESM_Classification`

| Property | Value |
|----------|-------|
| **Title** | Experimental Pipeline: ESM-Driven Classification and Prediction of Known Binders |
| **Type** | Full Research Paper |
| **Paper Status** | 🔴 Skeleton (abstract + empty sections) |
| **Presentation Status** | 🟡 Scaffold |
| **Description** | A complementary approach to *de novo* generation. Uses ESM and ESM3 sequence representations to classify known laboratory binders into machine-learned functional classes. These representations predict new target-binding properties from unseen sequences. Candidates flagged by the predictive model are evaluated using the same downstream protocols as generated constructs. |
| **Codebase Sources** | `ESM2 Classification/` (Training_Layered.ipynb, Training_Multi_target.ipynb, Run_Prediction.ipynb) |
| **Key Data Paths** | Training data from `Data/TIMP3_Binding_Results/`, prediction outputs from ESM notebooks |

---

### Thread 4: `FASTQ_Workflow`

| Property | Value |
|----------|-------|
| **Title** | Standalone Workflow: Automation of FASTQ Sequence Parsing |
| **Type** | Technical Workflow Report |
| **Paper Status** | 🔴 Skeleton (abstract + empty sections) |
| **Presentation Status** | 🟡 Scaffold |
| **Description** | An independent, streamlined pipeline for rapidly parsing sequencing outputs. Identifies specific motif hits or sequence matches and assembles read-level insights from broad sequencing datasets. Provides structured, reproducible sorting architectures for FASTQ configurations. |
| **Codebase Sources** | `FASTQ/` (read_fastq.py, find_hits.py, CDRH_unique.py, ID_sequences.py, compare_similarity.py, NGS_SOP.ipynb) |
| **Key Data Paths** | FASTQ input files (external), processed outputs from FASTQ scripts |

---

### Thread 5: `Metal_Binder_Generation`

| Property | Value |
|----------|-------|
| **Title** | Experimental Pipeline: Scaffold-Based Generation of Specific Metal Binders |
| **Type** | Full Research Paper |
| **Paper Status** | 🔴 Skeleton (abstract + empty sections) |
| **Presentation Status** | 🟡 Scaffold |
| **Description** | Adapts the generative deep-learning pipeline to synthesize novel metal-binding metalloproteins. Uses a scaffold pre-disposed to rare earth elements and redesigns loop regions and coordination sites to selectively uptake distinct rare earth or heavy metal ions. Explores how discrete conformational changes dictate selective metal coordination. |
| **Codebase Sources** | `MetalBinder/` (run_pipeline.py, analyze_results.py, dashboard.py, cross_docking_test.py, generate_summary_tables.py), `Generation/` (RFd3_batch.py, run_batch_allmetal3d.py) |
| **Key Data Paths** | `Local/lanm_output/` (global_sequence_catalog.csv, cross_docking_catalog.csv), `Local/Metal_Catalog/`, `Local/Metal_PDBs/`, `Local/Metal_Predictions/` |
| **Key Metrics** | overall_rmsd, loop_rmsd, structural_deviation_rmsd, plddt, ptm, binding_radius_A |

---

### Thread 6: `Nanodrop_Workflow`

| Property | Value |
|----------|-------|
| **Title** | Standalone Workflow: Integration of Nanodrop Quantification Data |
| **Type** | Technical Workflow Report |
| **Paper Status** | 🔴 Skeleton (abstract + empty sections) |
| **Presentation Status** | 🟡 Scaffold |
| **Description** | Covers the computational handling of Nanodrop spectroscopic data. Validates DNA purity and standardizes nucleic acid concentrations early in sample preparation. An independent parser digests numerical summaries from instrumentation, organizing absorption data and contamination estimates into analyzable downstream arrays. |
| **Codebase Sources** | `DNA Analysis/` (nanodrop_parser.py) |
| **Key Data Paths** | Nanodrop instrument exports (external), parsed outputs from nanodrop_parser.py |

---

### Thread 7: `TIMP_Dashboard_Pipeline`

| Property | Value |
|----------|-------|
| **Title** | Experimental Pipeline: Utilizing Advanced RoseTTAFold Architectures for TIMP Analysis |
| **Type** | Full Research Paper |
| **Paper Status** | 🔴 Skeleton (abstract + empty sections) |
| **Presentation Status** | 🟡 Scaffold |
| **Description** | An upgraded computational pipeline employing RoseTTAFold Diffusion 3 (RFd3) and LigandMPNN (LMPNN) for refined sequence hallucination, verified with RoseTTAFold 3 (RF3). The multi-dimensional binding and confidence metrics are systematically ingested into an interactive TIMP-Dashboard. Provides centralized, dynamic exploration of cross-docking performance, loop variations, and binding specificities. |
| **Codebase Sources** | `TIMP-Dashboard/` (dashboard.py, run_pipeline.py, analyze_results.py, cross_docking_analysis.py), `Generation/` (RFd3_batch.py, test_lmpnn_fixed.py, clean_renum_struct.py) |
| **Key Data Paths** | `Local/TIMP-Dashboard_output/`, `Local/rfd3_output/`, `Local/ligandmpnn_output/`, `Local/rf3_output/` |

---

### Thread 8: `ADAM_Target_Production`

| Property | Value |
|----------|-------|
| **Title** | Recombinant Production of ADAM10, ADAM17, and MMP7 Target Metalloproteinases for Yeast-Display Selectivity Screening |
| **Type** | Technical Paper / Working Record |
| **Paper Status** | 🟢 Advanced Draft (figures, tables, updated with June 10th results) |
| **Presentation Status** | 🟢 Draft |
| **Description** | A working record for the expression, purification, and validation of target metalloproteinases (ADAM10, ADAM17, MMP7) used for yeast-display selectivity screening. Tracks SDS-PAGE runs, silver stain enrichment QC, spin concentration yields, step-wise refolding, and binding performance in target flow cytometry trials. |
| **Codebase Sources** | `Protein-Analysis/` (gel_annotator.py, analyze_fcs_target_trial.py) |
| **Key Data Paths** | `Local/20260610_170522_Analysis/` (June 10 flow trial), `SharedAssets/figures/De_Novo_Binder_Generation/` (gel images) |
| **Key Figures** | `ADAM1017_SDS-Page_Gel_20260610_ADAM10.png`, `ADAM1017_SDS-Page_Gel_20260610_ADAM17.png` |

---

## Directory Architecture

```
PhD-Research/
├── Demonstrations/                 ← YOU ARE HERE
│   ├── README.md                   ← This file (source of truth)
│   ├── Papers/                     ← LaTeX manuscripts
│   │   └── <ThreadID>/
│   │       ├── <ThreadID>.tex      ← Main LaTeX source
│   │       ├── <ThreadID>.pdf      ← Compiled output
│   │       └── *.aux, *.log        ← Build artifacts
│   ├── Presentations/              ← Slidev slide decks
│   │   └── <ThreadID>/
│   │       ├── slides.md           ← Main Slidev markdown
│   │       └── components/         ← Optional Vue SFCs for this thread
│   └── SharedAssets/               ← Shared data, figures, and components
│       ├── figures/                ← Exported PNGs/SVGs used by both
│       ├── data/                   ← Cleaned CSVs/JSONs ready for rendering
│       └── components/             ← Shared Vue components or LaTeX macros
│
├── AlphaFold/                      ← AlphaFold analysis notebooks
├── Data/                           ← Raw/reference data (PDBs, FASTAs, etc.)
├── DNA Analysis/                   ← Nanodrop parsing scripts
├── ESM2 Classification/            ← ESM2/ESM3 classification notebooks
├── FASTA Mod/                      ← FASTA manipulation utilities
├── FASTQ/                          ← FASTQ parsing and analysis
├── Generation/                     ← RFdiffusion, ProteinMPNN, RFd3, HADDOCK
├── Local/                          ← Generated outputs (.gitignored)
├── MetalBinder/                    ← Metal binder pipeline and dashboard
├── Protein-Analysis/               ← Flow cytometry analysis pipeline
├── Results/                        ← Curated result files
├── SLURM/                          ← HPC job scripts
├── TIMP-Dashboard/                 ← TIMP interactive dashboard
└── Tools/                          ← External tool repositories (.gitignored)
```

---

## Architectural Conventions

### Papers (LaTeX)
- **One folder per thread**: `Papers/<ThreadID>/<ThreadID>.tex`
- **Document class**: `\documentclass{article}` with `geometry` (letterpaper, 1in margins)
- **Encoding**: `\usepackage[utf8]{inputenc}`
- **Graphics**: `\usepackage{graphicx}` — reference figures from `../SharedAssets/figures/` or `../../../Local/` when SharedAssets copy is not yet available
- **Compilation**: `pdflatex <ThreadID>.tex` from within the thread folder

### Presentations (Slidev)
- **One folder per thread**: `Presentations/<ThreadID>/slides.md`
- **Framework**: [Slidev](https://sli.dev/) with Vue SFC support
- **Theme**: `seriph` (default; can be overridden per thread)
- **Running**: `npx slidev Presentations/<ThreadID>/slides.md` from `Demonstrations/` root
- **Assets**: Reference shared figures via relative paths `../../SharedAssets/figures/`

### SharedAssets
- **Figures**: Publication-quality PNGs or SVGs at consistent DPI (300+ for print, 150 for slides)
- **Data**: Cleaned CSVs or JSONs with documented schemas; no raw instrument outputs
- **Components**: Reusable Vue SFCs (e.g., protein structure viewers, metric charts) or LaTeX macro files

---

## Data-to-Demonstration Pipeline

When you need to add data, figures, or tables to a Paper or Presentation, follow this pipeline:

### Step 1: Identify Raw Data Sources
Locate the relevant raw data in the codebase. Check `Local/` first (the primary data source), then `Data/`, `Results/`, and thread-specific directories. **Remember: `Local/` is only on the lab computer and its contents may change.** If a file you need is missing, ask the user whether it can be made available (see [Data Availability](#data-availability--local-storage) above).

Common locations:

| Data Type | Location | Example Files | Availability |
|-----------|----------|---------------|-------------|
| AlphaFold scores | `Local/` subdirectories | `*_scores_rank_*.json`, `*_predicted_aligned_error_v1.json` | Lab computer only |
| Flow cytometry | `Local/Aggregate_FCS_Analysis/` | `.fcs` files, parsed CSVs | Lab computer only |
| MetalBinder results | `Local/lanm_output/` | `global_sequence_catalog.csv`, `cross_docking_catalog.csv` | Lab computer only |
| TIMP Dashboard | `Local/TIMP-Dashboard_output/` | Pipeline outputs | Lab computer only |
| RFdiffusion outputs | `Local/rfdiffusion_output/` | `.pdb` files | Lab computer only |
| ProteinMPNN outputs | `Local/proteinmpnn_output/` | `.fa` files | Lab computer only |
| Best loop rankings | `Results/` | `best_loops_per_length_*.csv` | ✅ In repository |
| Gating figures | `Local/*_Analysis/` | `Gating_Strategy_*.png` | Lab computer only |
| Reference structures | `Data/` | `.pdb`, `.cif`, `.xlsx` files | ✅ In repository |
| Sequence catalogs | `Data/Sequence_Catalog/` | `.csv`, `.fa` files | ✅ In repository |

### Step 2: Clean and Export to SharedAssets
Write or invoke a script that:
1. Reads the raw data from `Local/` or other source directories
2. Cleans, filters, and formats it for presentation
3. Outputs cleaned data to `SharedAssets/data/<ThreadID>/`
4. Generates figures to `SharedAssets/figures/<ThreadID>/`

**Naming convention for exports:**
```
SharedAssets/
├── data/
│   └── De_Novo_Binder_Generation/
│       ├── zscore_summary.csv
│       └── top_variants_raw.csv
└── figures/
    └── De_Novo_Binder_Generation/
        ├── gating_strategy.png
        └── aggregate_binding.png
```

### Step 3: Reference from Papers and Presentations

**In LaTeX:**
```latex
\includegraphics[width=0.8\textwidth]{../SharedAssets/figures/De_Novo_Binder_Generation/gating_strategy.png}
```

**In Slidev:**
```markdown
![Gating Strategy](../../SharedAssets/figures/De_Novo_Binder_Generation/gating_strategy.png)
```

**In Vue components:**
```vue
<template>
  <img :src="'../../SharedAssets/figures/' + threadId + '/gating_strategy.png'" />
</template>
```

### Step 4: Document the Export
When you create a new data export or figure, add its path and description to the relevant thread's entry in the Thread Registry above.

---

## Thread Management Protocol

### Adding a New Thread
1. **Add to this README**: Create a new entry in the Thread Registry above with all required fields
2. **Create Paper folder**: `Papers/<NewThreadID>/<NewThreadID>.tex` — at minimum, include `\documentclass`, title, author, date, and an abstract
3. **Create Presentation folder**: `Presentations/<NewThreadID>/slides.md` — at minimum, include frontmatter, title slide, and section outlines
4. **Create SharedAssets subfolders**: `SharedAssets/figures/<NewThreadID>/` and `SharedAssets/data/<NewThreadID>/` if needed

### Modifying a Thread
1. **Update the registry** in this README first (description, status, data paths)
2. **Then update** the Paper and Presentation to reflect the changes

### Deprecating a Thread
1. **Mark as deprecated** in the registry (add `⛔ Deprecated` status)
2. **Do not delete** the folders — they serve as historical records
3. **Add a deprecation note** in both the `.tex` and `slides.md` files

### Status Labels
| Label | Meaning |
|-------|---------|
| ✅ Complete | Finished and publication-ready |
| 🟢 Advanced Draft | Substantial content, figures, tables; needs polish |
| 🟡 Scaffold | Structure and section headers exist; content is placeholder or TODO |
| 🔴 Skeleton | Only abstract and empty section headers |
| 🚧 In Progress | Actively being worked on |
| ⛔ Deprecated | No longer maintained |

---

## LaTeX Paper Conventions

### Standard Preamble
Every paper should start with:
```latex
\documentclass{article}
\usepackage[utf8]{inputenc}
\usepackage{geometry}
\usepackage{graphicx}
\geometry{letterpaper, margin=1in}

\title{<Paper Title>}
\author{Ryan Gustafson}
\date{\today}

\begin{document}
\maketitle
```

### Required Sections
At minimum, each paper must contain:
1. **Abstract** — Self-contained summary
2. **Introduction** — Context, motivation, scope
3. **Methods/Pipeline** — Technical approach (subsections for each tool/phase)
4. **Results** — Data, tables, figures with captions
5. **Conclusions** — Summary and future work

### Figure and Table Conventions
- All figures must have `\caption{}` and `\label{fig:<name>}`
- All tables must have `\caption{}` and `\label{tab:<name>}`
- Reference figures from SharedAssets when available: `../SharedAssets/figures/<ThreadID>/`
- Fallback to repository-relative paths for figures not yet in SharedAssets: `../../../Local/...`

### Building Papers
From within a thread's folder:
```bash
cd Demonstrations/Papers/<ThreadID>/
pdflatex <ThreadID>.tex
```

---

## Slidev Presentation Conventions

### Frontmatter Template
Every presentation's `slides.md` should begin with:
```yaml
---
theme: seriph
background: https://cover.sli.dev
title: "<Presentation Title>"
info: |
  ## <Thread Name>
  <One-line description>
  
  Author: Ryan Gustafson
class: text-center
drawings:
  persist: false
transition: slide-left
mdc: true
---
```

### Slide Structure
1. **Title Slide** — Thread title, subtitle, author attribution
2. **Overview/Outline** — Brief roadmap of the presentation
3. **Content Slides** — One slide per major section in the paper, with:
   - Key bullet points (not paragraph text)
   - Embedded figures from SharedAssets
   - Code snippets or diagrams where applicable
4. **Results Slides** — Data visualizations, tables, key metrics
5. **Conclusions Slide** — Summary and future directions

### Slidev Features to Use
- **Layouts**: `default`, `two-cols`, `image-right`, `center`, `fact`
- **Transitions**: `slide-left`, `fade-out`, `slide-up`
- **Mermaid diagrams**: For pipeline/workflow visualizations
- **Vue components**: For interactive elements (protein viewers, charts)
- **Code blocks**: For showing pipeline commands or scripts

### Running Presentations
```bash
# From the Demonstrations directory:
npx @slidev/cli Presentations/<ThreadID>/slides.md

# Or install globally:
npm i -g @slidev/cli
slidev Presentations/<ThreadID>/slides.md
```

### TODO Markers
When content needs to be populated from the codebase, use HTML comments:
```markdown
<!-- TODO: Insert <description> from <source path> -->
```
These markers tell you exactly what data to pull and from where.

---

## SharedAssets Conventions

### Figure Requirements
| Context | Format | DPI | Max Width |
|---------|--------|-----|-----------|
| Paper (print) | PNG or PDF | 300+ | Full column width |
| Presentation (screen) | PNG or SVG | 150 | 1920px |
| Shared (both) | PNG | 300 | 1920px |

When a figure is used by both a Paper and a Presentation, export at 300 DPI / 1920px width to satisfy both requirements.

### Data File Requirements
- Format: CSV (preferred) or JSON
- Include a header row with descriptive column names
- Include a `_schema.md` file alongside each dataset documenting:
  - Column descriptions
  - Units
  - Source script/notebook that produced the data
  - Date of last generation

### Component Organization
```
SharedAssets/components/
├── vue/                    ← Shared Vue SFCs for Slidev
│   ├── ProteinViewer.vue
│   └── MetricChart.vue
└── latex/                  ← Shared LaTeX macros
    └── common_macros.tex
```

---

## Codebase Map

Quick reference for finding scripts and data relevant to each thread.

### Pipeline Tools (Generation/)
| Script | Purpose |
|--------|---------|
| `RFd_Batch.py` | Batch RFdiffusion backbone generation |
| `RFd3_batch.py` | Batch RFd3 (v3) backbone generation |
| `ProMPNN_Batch.py` | Batch ProteinMPNN sequence design |
| `test_lmpnn_fixed.py` | LigandMPNN testing |
| `RFDiffusion_Gen.ipynb` | Interactive RFdiffusion generation |
| `RFd3_Gen.ipynb` | Interactive RFd3 generation |
| `ProteinMPNN_Gen.ipynb` | Interactive ProteinMPNN generation |
| `TwistBioOrder.ipynb` | Codon optimization and synthesis ordering |
| `clean_renum_struct.py` | Structure cleaning and renumbering |
| `run_haddock_timp_mp.py` | HADDOCK3 docking pipeline |
| `run_batch_allmetal3d.py` | AllMetal3D batch metal prediction |

### Analysis Tools (AlphaFold/)
| Script | Purpose |
|--------|---------|
| `AlphaFoldServer_Analysis.ipynb` | Primary AF score extraction and Z-score analysis |
| `AlphaFoldServer_MDA.ipynb` | Multi-dimensional AF analysis |
| `Top_loops_AF_gen.py` | Top loop ranking by AF metrics |
| `AF_batch_gen.py` | Batch AF job generation |

### Protein Analysis (Protein-Analysis/)
| Script | Purpose |
|--------|---------|
| `analyze_fcs.py` | Single-target FCS flow cytometry analysis |
| `analyze_fcs_target_trial.py` | Per-target per-trial FCS analysis |
| `aggregate_analysis.py` | Cross-target aggregate analysis |
| `fcs_viewer.py` | Interactive FCS data viewer (Streamlit) |
| `rename_fcs.py` | FCS file renaming utility |

### ESM Classification (ESM2 Classification/)
| Script | Purpose |
|--------|---------|
| `Training_Layered.ipynb` | Layered ESM embedding training |
| `Training_Multi_target.ipynb` | Multi-target classifier training |
| `Run_Prediction.ipynb` | Prediction on unseen sequences |

### Sequence Utilities (FASTQ/, FASTA Mod/, DNA Analysis/)
| Script | Purpose |
|--------|---------|
| `FASTQ/read_fastq.py` | FASTQ file parsing |
| `FASTQ/find_hits.py` | Motif hit detection in reads |
| `FASTQ/CDRH_unique.py` | CDR-H3 unique sequence extraction |
| `FASTA Mod/make_fasta.py` | FASTA file generation |
| `DNA Analysis/nanodrop_parser.py` | Nanodrop data parsing |

### Dashboards
| Script | Purpose |
|--------|---------|
| `TIMP-Dashboard/dashboard.py` | TIMP variant exploration dashboard |
| `MetalBinder/dashboard.py` | Metal binder exploration dashboard |
| `Protein-Analysis/fcs_viewer.py` | Flow cytometry data viewer |

---

## Quick Start for AI Agents

1. **Read this README completely** before making any changes to `Demonstrations/`.
2. **Check the Thread Registry** to understand what exists and its current status.
3. **Follow the Data-to-Demonstration Pipeline** when adding new data, figures, or tables.
4. **Update the Thread Registry** status labels after making substantive changes.
5. **Use TODO markers** in `slides.md` files to flag content that needs codebase data.
6. **Build both outputs**: When you add content to a Paper, also add the corresponding content to its Presentation (and vice versa).
7. **Keep SharedAssets current**: Any figure or data file used by either a Paper or Presentation should live in `SharedAssets/` with proper naming and schema documentation.
