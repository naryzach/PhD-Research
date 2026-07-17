"""
Central configuration for the TIMP3 structural-validation project.

Single source of truth for: the construct/target panel, sequence boundaries,
paths, crystal-structure map, and analysis parameters. Every other script
imports from here so the pipeline stays internally consistent.

Decisions locked 2026-07-08 (see README.md):
  * TIMP3 modeling unit = N-terminal binding domain (~121 aa, CT...RYHLGCN)
  * Target panel       = MMP2, MMP9, ADAM17, MMP3, ADAM10  (purchased + in-house)
  * Compute model      = this repo PREPARES inputs/scripts; the heavy folding
                         (AF3 web server, ESMFold2 on the cluster) and docking
                         (HADDOCK3 in the `haddock` conda env) run on your side.
"""
from __future__ import annotations

from pathlib import Path

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = REPO_ROOT / "Structural-Validation"

DATA_DIR = REPO_ROOT / "Data"
SEQ_CATALOG_DIR = DATA_DIR / "Sequence_Catalog"
CRYSTAL_DIR = DATA_DIR / "Target_Crystal_Structures"
TARGETS_XLSX = DATA_DIR / "Targets_w_CDs.xlsx"
CONSTRUCTS_CSV = (
    REPO_ROOT
    / "Local/TIMP3_Redesign_2026-07/data/esmc_predict_input_fcs_constructs.csv"
)

# All generated artifacts land here (never inside the code folder).
OUT_DIR = REPO_ROOT / "Local/TIMP3_Structural_Validation_2026-07"
OUT_SEQ = OUT_DIR / "sequences"
OUT_AF3 = OUT_DIR / "af3"
OUT_ESM = OUT_DIR / "esmfold2"
OUT_CRYSTAL = OUT_DIR / "crystal"
OUT_DOCK = OUT_DIR / "docking"
OUT_CLEAN_TARGETS = OUT_CRYSTAL / "clean_targets"   # CD-only single-chain crystals
OUT_ANALYSIS = OUT_DIR / "analysis"
OUT_FIG = OUT_DIR / "figures"
OUT_REPORT = OUT_DIR / "reports"

ALL_OUT_SUBDIRS = [
    OUT_SEQ, OUT_AF3, OUT_ESM, OUT_CRYSTAL, OUT_CLEAN_TARGETS,
    OUT_DOCK, OUT_ANALYSIS, OUT_FIG, OUT_REPORT,
]

# --------------------------------------------------------------------------
# TIMP3 constructs
# --------------------------------------------------------------------------
# The FCS construct CSV records sequences starting at the 2nd Cys ("CSPSHPQD..")
# and running through the full two-domain protein. The mature protein and every
# prior AF3/HADDOCK job used the N-terminal *binding* domain with the native
# "CT" dipeptide restored and the chain truncated at the N-domain terminus.
TIMP3_NDOMAIN_PREFIX = "CT"          # restores CTCSP... (mature N-terminus)
TIMP3_NDOMAIN_END_MOTIF = "RYHLGCN"  # last residues of the N-terminal domain

# Wild-type mature TIMP3 (UniProt P35625, signal/pro 1-23 removed).
TIMP3_WT_MATURE = (
    "CTCSPSHPQDAFCNSDIVIRAKVVGKKLVKEGPFGTLVYTIKQMKMYRGFTKMPHVQYIHTE"
    "ASESLCGLKLEVNKYQYLLTGRVYDGKMYTGLCNFVERWDQLTLSQRKGLNYRYHLGCNCKI"
    "KSCYYLPCFVTSKNECLWTDMLSNFGYPGYQSKHYACIRQKGGYCSWYRGWAPPDKSIINATDP"
)

# --------------------------------------------------------------------------
# Target panel (Core 5). Sequences are read from Targets_w_CDs.xlsx at runtime;
# this table carries the metadata that isn't in the spreadsheet.
# --------------------------------------------------------------------------
TARGETS = {
    "MMP2":   {"vendor": "Sino/Enzo",    "source": "purchased",
               "uniprot": "P08253", "crystal": "MMP2_Xray.pdb",
               "af_cif": "MMP2cd_AF.cif"},
    "MMP9":   {"vendor": "Masoud + Sino/Enzo", "source": "in-house+purchased",
               "uniprot": "P14780", "crystal": "MMP9_Xray.pdb",
               "af_cif": "MMP9cd_AF.cif"},
    "ADAM17": {"vendor": "Enzo/Sino",    "source": "purchased",
               "uniprot": "P78536", "crystal": "ADAM17_Xray.pdb",
               "af_cif": "ADAM17cd_AF.cif"},
    "MMP3":   {"vendor": "Masoud (MMP3-HT)", "source": "in-house",
               "uniprot": "P08254", "crystal": "MMP3cd_X_ray.pdb",
               "af_cif": "MMP3cd_AF.cif"},
    "ADAM10": {"vendor": "Enzo",        "source": "purchased",
               "uniprot": "O14672", "crystal": "ADAM10_Xray.pdb",
               "af_cif": "ADAM10cd_AF.cif"},
}
TARGET_ORDER = ["MMP2", "MMP9", "ADAM17", "MMP3", "ADAM10"]

# Reference *complex* crystal structures for docking QA / DockQ.
#   native      = an actual TIMP3:<this target> co-crystal (exact benchmark)
#   approximate = a homologous TIMP:metalloprotease co-crystal (conserved binding
#                 mode) used only where no native complex exists. DockQ against
#                 these is APPROXIMATE and flagged as such in the output.
HOMOLOG_DIR = OUT_CRYSTAL / "references"   # downloaded by crystal/fetch_homologous_refs.py
COMPLEX_REFERENCES = {
    "ADAM17": {"native": [DATA_DIR / "TIMP3_vs_ADAM17_X_ray.pdb",
                          DATA_DIR / "3CKI_Sarmazdeh.pdb"]},
    # ADAM10: sister ADAM, same TIMP3 ligand -> reuse the TIMP3:ADAM17 complex.
    "ADAM10": {"approximate": [DATA_DIR / "TIMP3_vs_ADAM17_X_ray.pdb"]},
    # MMPs: TIMP-1:MMP-3 catalytic complex (1UEA). MMP3 gets an exact receptor.
    "MMP3": {"approximate": [HOMOLOG_DIR / "1UEA.pdb"]},
    "MMP2": {"approximate": [HOMOLOG_DIR / "1UEA.pdb"]},
    "MMP9": {"approximate": [HOMOLOG_DIR / "1UEA.pdb"]},
}

TIMP3_CRYSTAL = CRYSTAL_DIR / "TIMP3_Xray.pdb"

# --------------------------------------------------------------------------
# Interface / restraint definition
# --------------------------------------------------------------------------
# TIMP3 residues that coordinate the catalytic zinc (N-terminal ridge, mature
# numbering with the CT prefix present): C1-T2-C3-S4-P5 is the reactive edge.
TIMP3_ACTIVE = (1, 5)
TIMP3_PASSIVE = (6, 10)

# Each target's zinc-binding motif (HExxHxxGxxH) is taken from the spreadsheet
# "Active_site" column and located in the CD sequence at runtime, so HADDOCK
# restraints and interface metrics don't depend on hand-typed residue ranges.

# --------------------------------------------------------------------------
# Fold modelers under comparison
# --------------------------------------------------------------------------
FOLDERS = ["AF3", "ESMFold2"]

# Sources of docking *input* structures. Construct always comes from a fold
# model (no construct crystal exists); the target may also come from the
# extracted crystal, so docked poses can use the ground-truth target.
CONSTRUCT_SOURCES = ["AF3", "ESMFold2"]
TARGET_SOURCES = ["AF3", "ESMFold2", "Crystal"]

# HADDOCK docking matrix presets: (construct_source, target_source) pairs.
DOCK_MATRIX = {
    # literal "try HADDOCK with both AF3 and ESMFold2 and compare"
    "compare": [("AF3", "AF3"), ("ESMFold2", "ESMFold2")],
    # "extract the crystal target and redock each construct" (consistent target)
    "crystal": [("AF3", "Crystal"), ("ESMFold2", "Crystal")],
    # everything
    "full": [("AF3", "AF3"), ("ESMFold2", "ESMFold2"),
             ("AF3", "Crystal"), ("ESMFold2", "Crystal")],
}
DOCK_MATRIX_DEFAULT = "full"

# Complex model sources scored by the complex analysis (co-folds + docked tracks
# are discovered dynamically; these are the two native co-fold predictors).
COFOLD_SOURCES = ["AF3_cofold", "ESMFold2_cofold"]

# --------------------------------------------------------------------------
# HADDOCK QA thresholds (carried over from Generation/run_haddock_timp_mp.py)
# --------------------------------------------------------------------------
HADDOCK_SCORE_GOOD = -80.0
HADDOCK_BSA_GOOD = 1200.0
HADDOCK_DIST_GOOD = 9.0
HADDOCK_VIOL_MAX = 3

# --------------------------------------------------------------------------
# Composite interface score (transparent, weights exposed — not buried)
# --------------------------------------------------------------------------
# Re-implements the AlphaFoldServer_SA.ipynb TIMP-tuned heuristic, but the raw
# metrics are always reported alongside it so the composite is never a black box.
# Confidence terms (ipTM, PAE) are only available for co-folds; for docked
# (HADDOCK) poses those terms are dropped and the remaining weights renormalised.
COMPOSITE_WEIGHTS = {
    "iptm": 0.35,        # AF3/ESM interface confidence
    "bsa": 0.25,         # buried surface area (critical for TIMP interfaces)
    "pae": 0.20,         # interface PAE (low = good)
    "hbond": 0.20,       # H-bond count
    "bsa_cap": 1000.0,   # BSA saturates the score at this many Å²
    "hbond_cap": 15.0,   # H-bonds saturate here (TIMP interfaces are large)
    "pae_norm": 15.0,    # PAE normalised 0..this (Å)
    "salt_bonus_per": 0.02, "salt_bonus_max": 0.10,
    "hydro_bonus_per": 0.01, "hydro_bonus_max": 0.05,
}


def ensure_out_dirs() -> None:
    """Create every output subdirectory (idempotent)."""
    for d in ALL_OUT_SUBDIRS:
        d.mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":
    ensure_out_dirs()
    print(f"Repo root : {REPO_ROOT}")
    print(f"Output dir: {OUT_DIR}")
    for d in ALL_OUT_SUBDIRS:
        print(f"  created  {d.relative_to(REPO_ROOT)}")
