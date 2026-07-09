# TIMP3 Construct × Target — Structural Validation

In-silico structural characterization of every engineered TIMP3 construct against
its protease targets: fold each entity with **AF3** and **ESMFold2**, dock every
construct×target pair with **HADDOCK3**, score the models and complexes on a broad
panel of structural metrics, and benchmark the predicted structures against known
**crystal structures**. Flow-cytometry correlation is deliberately **out of scope
for now** (held until the FCS data is finalized).

- **Code** lives here (`Structural-Validation/`).
- **All outputs** go to `Local/TIMP3_Structural_Validation_2026-07/` — nothing is
  written inside the code folder.

## Locked decisions (2026-07-08)
| Decision | Choice | Why |
|---|---|---|
| TIMP3 modeling unit | **N-terminal binding domain** (~121 aa, `CT…RYHLGCN`) | the MMP/ADAM-binding unit; matches prior AF3/HADDOCK runs and the FCS-calibrated metrics |
| Target panel | **MMP2, MMP9, ADAM17, MMP3, ADAM10** | all purchased + in-house targets; each has a CD sequence and a crystal on disk |
| Compute model | this repo **prepares** inputs/scripts; **you run** AF3/ESMFold2/HADDOCK | AF3 is a manual web server (~30 jobs/day), ESMFold2 needs the GPU cluster, HADDOCK needs the `haddock` env |

## Panel (19 entities)
- **14 constructs** (N-domain): `TIMP3_WT` + `AB_1..7` + `C_11..16`
- **5 targets** (catalytic domain): `MMP2 MMP9 ADAM17 MMP3 ADAM10`
- Sequences are assembled from repo data (no external fetch except WT TIMP3,
  UniProt P35625) into `sequences/registry.csv`.

## Layout
```
Structural-Validation/
  config.py                 single source of truth (panel, paths, thresholds)
  run_all.py                orchestrator: prep | compute | analysis | all
  sequences/build_sequences.py      -> registry.csv, all_sequences.fasta, fasta/
  folding/make_af3_json.py          -> AF3 monomer + chunked complex JSON
  folding/make_esmfold_inputs.py    -> ESMFold2 FASTA + SLURM job
  folding/run_esmfold2.py           ESMFold2 runner (cluster, esmfold2 env)
  crystal/catalog_structures.py     inventory + map crystal references
  crystal/extract_targets.py        CD-only clean crystal targets (consistent docking)
  docking/run_haddock_all.py        construct×target HADDOCK3 (motif-driven AIRs),
                                    multi-track (--matrix compare|crystal|full)
  analysis/metrics.py               self-contained metric library
  analysis/structure_io helpers in  utils/structure_io.py
  analysis/model_registry.py        tolerant output path discovery
  analysis/model_vs_crystal.py      monomer quality + AF3↔ESMFold2 agreement
  analysis/complex_metrics.py       interface battery + HADDOCK energies + DockQ
  analysis/aggregate.py             tidy tables, matrices, correlations, dict
  analysis/visualize.py             figures + reports/dashboard.html
```

## How to run
```bash
# 1. Prepare everything (runs anywhere; needs only pandas/openpyxl/biopython)
python Structural-Validation/run_all.py prep

# 2. Compute (your side) — see printed instructions
python Structural-Validation/run_all.py compute
#   AF3:      upload af3/*.json at https://alphafoldserver.com, unzip to af3/results/
#   ESMFold2: sbatch …/esmfold2/run_esmfold2.slurm   (folds monomers AND complexes)
#   HADDOCK:  python Structural-Validation/docking/run_haddock_all.py --matrix full --jobs 8

# 3. Analyse (re-runnable any time; fills in as outputs arrive)
python Structural-Validation/run_all.py analysis
```
The analysis stage is safe to run before compute finishes — every table has a
`status` column and unfinished cells read `pending_*`.

## Metrics computed (see `analysis/data_dictionary.csv` for all 35)
**Monomer vs crystal:** CA-RMSD, TM-score (over shared region), GDT-TS, GDT-HA,
lDDT, radius of gyration, clashscore, Ramachandran-outlier fraction, pLDDT, pTM.
**AF3 ↔ ESMFold2 agreement:** cross RMSD and TM-score.
**Complex / interface:** buried surface area, interface-residue counts, 5 Å
contacts, H-bonds, salt bridges, hydrophobic contacts, contact density, min
CA-CA from the TIMP3 reactive edge to the target HExxHxxGxxH zinc motif, HADDOCK
score + vdW/elec/desolv/AIR components + reported BSA + AIR violations, AF3
ipTM/pTM/PAE (co-folds), and **DockQ** (fnat, iRMS, LRMS, CAPRI class) against a
native complex where one exists.

The metric library is pure Python/NumPy/BioPython (validated against self- and
crystal comparisons: self-RMSD 0 / TM 1 / DockQ 1; AF ADAM17 CD vs crystal
RMSD 1.56 Å / TM 0.98). No PyMOL/TMalign/DockQ/freesasa binaries required.

## Complex model sources (the `source` column in master_complex_metrics.csv)
Every construct×target pair is scored across up to six independent complex sources,
so you can see which modeler/route gives the best information:
| source | how the complex is built |
|---|---|
| `HADDOCK:AF3xAF3` | dock AF3 construct + AF3 target |
| `HADDOCK:ESMFold2xESMFold2` | dock ESMFold2 construct + ESMFold2 target |
| `HADDOCK:AF3xCrystal` | dock AF3 construct + **extracted crystal** target |
| `HADDOCK:ESMFold2xCrystal` | dock ESMFold2 construct + **extracted crystal** target |
| `AF3_cofold` | AF3 co-fold of the two chains together |
| `ESMFold2_cofold` | ESMFold2 co-fold (esm_iptm / esm_pae / loop-pLDDT) |

## Caveats
- **DockQ/CAPRI need a bound-complex co-crystal — only ADAM17 has one** (`3CKI` /
  `TIMP3_vs_ADAM17`). Extracting + redocking targets gives input consistency but
  does **not** create a native complex reference; the other four targets get the
  full interface battery (BSA, contacts, H-bonds, salt bridges, zinc-loop geometry,
  energies) without a DockQ benchmark. Ask if you want homologous TIMP:MMP
  co-crystals used as approximate references.
- **Crystal targets are extracted CD-only** (`crystal/clean_targets/`): MMP2's
  fibronectin-II inserts and ADAM10's ectodomain are dropped so the docking target
  matches the purchased catalytic-domain construct. Model-vs-crystal TM is
  normalized over the aligned region.
- **MMP/ADAM crystals are often the catalytic-Glu→Ala inactive mutant** (e.g. MMP2
  E404A); the zinc-motif locator tolerates this (`H[EAQ]xxHxxGxxH`).
- Ramachandran-outlier fraction uses coarse basins — treat as a relative flag.
- **FCS correlation is intentionally not included yet.**
