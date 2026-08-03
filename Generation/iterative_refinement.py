"""
iterative_refinement.py

Iterative binder design with an in-silico funnel before any AF3 call:

  RFd3 (backbone) → LigandMPNN (sequence) → ESMFold2 (local AF3-class ranker) →
  AF3 Server (gold-standard validation, capped at 30/day)

ESMFold2 replaced Boltz-2 as the ranker after AF3 calibration (see
filter_methods.md §5): MSA-free, faster, and a better AF3 predictor.
RF3 and Boltz-2 remain in the code but are OFF by default (RF3_ENABLE / BOLTZ_ENABLE)
— RF3's confidence was anti-correlated with AF3, so it only added cost. Enable RF3
with --enable-rf3 to log its geometric features for curiosity.

TIMP3-scaffold binders for MMP2, MMP9, MMP3, MMP10, ADAM10, ADAM17 are produced by
redesigning loops AB / C / EF (GH optional).  Temperature anneals 0.5 → 0.1 to shift
LigandMPNN from exploration to exploitation; loop-length ranges are adaptively narrowed
toward what the HOF actually used (kicks in at iteration 3).

Scoring is source-aware (calc_composite):
  AF3-validated   →  full-trust composite using all AF3 metrics
  Boltz-scored    →  Boltz ipTM/pTM/pLDDT/PAE + binding-affinity head
  RF3-only        →  geometric features only (n_contacts, iface PAE, RMSD);
                     RF3 ipTM/pTM/pLDDT are excluded because they are anti-correlated
                     with AF3 ipTM for de novo binders (empirically calibrated).
                     This branch is capped at RF3_COMPOSITE_CEILING so any Boltz or
                     AF3 score automatically out-ranks an unvalidated RF3 entry.

Output layout:
  Local/iterative_refinement/
    refinement_state.json             persistent state (HOF, temperature, iteration)
    target_templates/                 target-only CIFs used as RF3 templates
    it_N/
      rfd3/                           RFd3 backbone CIFs
      lmpnn/                          LigandMPNN sequence CIFs
      rf3/                            RF3 prediction CIFs + metrics JSONs
      boltz/                          Boltz-2 prediction outputs (per design)
      round_summary.csv               all scored designs this iteration
    hof_structures/<target>/          best RF3 CIFs per target (for seeding)
    af3_submission_itN.json           AF3 Server input batch (top AF3_TOP_N)
    hof_summary.csv                   all-time best per target (updated each round)
"""

import os
import subprocess

# --- PORTABILITY: GPU-Aware Environment Setup ---
# We must detect the GPU and set DISABLE_CUEQUIVARIANCE before importing heavy ML libraries.
# cuEquivariance checks this at import time, so inside main() is too late.
try:
    smi_out = subprocess.check_output(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"]).decode()
    if "V100" in smi_out:
        os.environ["DISABLE_CUEQUIVARIANCE"] = "1"
        print(f"Detected V100 GPU. Automatically setting DISABLE_CUEQUIVARIANCE=1 for compatibility.")
except Exception:
    pass

import re
import sys
import time
import json
import logging
import dataclasses
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from biotite.structure import superimpose, rmsd as bio_rmsd
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile
from biotite.sequence import ProteinSequence
from atomworks.io.utils.io_utils import to_cif_file
from atomworks.constants import PROTEIN_BACKBONE_ATOM_NAMES

from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine
from rfd3.inference.input_parsing import DesignInputSpecification
from mpnn.inference_engines.mpnn import MPNNInferenceEngine
from rf3.inference_engines.rf3 import RF3InferenceEngine
from rf3.utils.inference import InferenceInput

# Calibrated in-silico -> binding priors (2026-07 exact-sequence, purchased-only
# calibration; see calibrated_scoring.py). Lives alongside this file in
# Generation/; make it importable whether run as a script or imported by
# specificity_refinement.
sys.path.insert(0, str(Path(__file__).parent.resolve()))
import calibrated_scoring as cs
import sv_bridge as svb

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)
for _noisy in ("transforms", "atomworks.io", "atomworks.ml", "foundry", "lightning"):
    logging.getLogger(_noisy).setLevel(logging.ERROR)

torch.set_float32_matmul_precision("medium")

# ── Paths ─────────────────────────────────────────────────────────────────────
_HERE     = Path(__file__).parent.resolve()
# AF3 co-fold complexes are the RFd3 design template (chain A = FULL-LENGTH 188-aa
# TIMP3, chain B = target, in the AF3-predicted binding pose). Two fixes vs the
# original:
#  (1) replaced the HADDOCK complexes, which placed TIMP3 ~30 A off-native (LRMS
#      ~30 A / DockQ ~0.05 vs AF3 co-fold LRMS ~2.3 A / DockQ ~0.76) — RFd3 was
#      inpainting loops into a wrong binding orientation.
#  (2) FULL LENGTH (188), not the N-domain (121): the pipeline orders + tests the
#      full protein (N-domain design + native C-domain), and the C-domain is needed
#      for ADAM10 binding (literature) — so design must see it. The C-domain
#      (122-188) is fixed scaffold in the contig; scaffold_len=188 in TARGETS.
# Regenerate templates with Generation/prep_af3_templates.py. Revert: point back to
# "HADDOCK_Outputs", set binder_chain="B"/target_chain="A", scaffold_len=121.
DATA_DIR  = _HERE / ".." / "Data" / "TIMP_Complexes" / "AF3_Templates"
# Output root is overridable via REFINE_OUT_BASE so an unrelated run (e.g. a fresh
# anneal) can target a separate directory without clobbering a preserved/ salvaged
# pool. Defaults to Local/iterative_refinement. rescore_pool.py and
# select_binders_to_order.py honor the same variable.
OUT_BASE  = Path(os.environ.get("REFINE_OUT_BASE") or (_HERE / ".." / "Local" / "iterative_refinement"))
CKPT_DIR  = _HERE / ".." / "Tools" / "foundry_checkpoints"

# ── Target definitions ────────────────────────────────────────────────────────
# binder_chain / target_chain refer to chain IDs in the SOURCE PDB only.
# They are used to read target_seqs and to construct the RFd3 contig string.
# After RFd3 runs, all output structures (RFd3, LMPNN, RF3) follow a fixed
# convention regardless of input — see DESIGN_BINDER_CHAIN below.
# pdb = AF3 co-fold template (chain A = full 188-aa TIMP3, chain B = target).
# scaffold_len = 188: RFd3 designs the N-domain loops (AB/C/EF, and GH at pos 127)
# with the native C-domain (122-188) held as fixed scaffold context.
TARGETS = {
    "MMP2":   {"pdb": "MMP2_TIMP3_AF3.pdb",   "binder_chain": "A", "target_chain": "B", "scaffold_len": 188},
    "MMP9":   {"pdb": "MMP9_TIMP3_AF3.pdb",   "binder_chain": "A", "target_chain": "B", "scaffold_len": 188},
    "MMP3":   {"pdb": "MMP3_TIMP3_AF3.pdb",   "binder_chain": "A", "target_chain": "B", "scaffold_len": 188},
    "MMP10":  {"pdb": "MMP10_TIMP3_AF3.pdb",  "binder_chain": "A", "target_chain": "B", "scaffold_len": 188},
    "ADAM10": {"pdb": "ADAM10_TIMP3_AF3.pdb",  "binder_chain": "A", "target_chain": "B", "scaffold_len": 188},
    "ADAM17": {"pdb": "ADAM17_TIMP3_AF3.pdb",  "binder_chain": "A", "target_chain": "B", "scaffold_len": 188},
}

# RFd3 emits chains in contig order: the first contig segment becomes chain A,
# anything after a "/0" chain break becomes B, etc.  Since our contig always
# puts the binder template first and the fixed target second, the convention
# downstream is always:
#   chain A = designed binder (TIMP3-derived)
#   chain B = fixed target (protease)
# These IDs are used everywhere we read from RFd3 / LMPNN / RF3 output arrays,
# and as the component IDs we hand to RF3.  They are independent of the source
# PDB's chain assignments above.
DESIGN_BINDER_CHAIN = "A"
DESIGN_TARGET_CHAIN = "B"

# ── Loop definitions ──────────────────────────────────────────────────────────
# pos: 1-indexed last fixed residue before the loop in the native scaffold.
# normal/max: native and maximum loop length for the contig string.
# left/right: flanking tripeptides for regex-based loop extraction from sequences.
LOOP_CONFIGS = {
    "AB":  {"normal": 6,  "max": 15, "pos": 30,  "left": "LVK", "right": "LVY"},
    "C":   {"normal": 6,  "max": 15, "pos": 62,  "left": "HTE", "right": "GLK"},
    "EF":  {"normal": 4,  "max": 10, "pos": 92,  "left": "MYT", "right": "FVE"},
    "GH":  {"normal": 10, "max": 20, "pos": 127, "left": "KSC", "right": "NEC"},
    # MTL = the C-terminal Multiple-Turn Loop (native DMLSNFGYPG). Lives in the
    # C-terminal domain (pos 143), so it requires the full-length TIMP3 construct
    # (not the 121-aa N-terminal design construct). Same loop the RFd_Batch /
    # RFd3_batch scripts label "Multi".
    "MTL": {"normal": 10, "max": 20, "pos": 143, "left": "LWT", "right": "YQS"},
}

# ── Hyperparameters ───────────────────────────────────────────────────────────
# Throughput per iteration = BACKBONES_PER_TARGET × LMPNN_SEQS_PER_BACKBONE × n_targets.
# RFd3 (backbones) is the slow stage; LMPNN sequences-per-backbone is nearly free
# and each one is one more (cheap) ESMFold2 prediction. Both are CLI-overridable
# (--backbones-per-target / --seqs-per-backbone).
BACKBONES_PER_TARGET   = 40    # RFd3 backbones per target per iteration (structural diversity; was 20)
LMPNN_SEQS_PER_BACKBONE = 1    # LMPNN sequences sampled per backbone. Raise (e.g. 2-3) via
                               # --seqs-per-backbone for cheap sequence diversity + more ESMFold2
                               # predictions without extra RFd3 cost.
INIT_TEMPERATURE       = 0.50  # LMPNN sampling temperature (exploration)
MIN_TEMPERATURE        = 0.10  # LMPNN temperature floor (exploitation)
TEMP_DECAY             = 0.85  # Temperature multiplier applied each iteration
HOF_SIZE_PER_TARGET    = 75    # Max Hall of Fame entries per target
ADAPTIVE_BIAS_START    = 3     # First iteration to apply adaptive loop-length bias
ADAPTIVE_BIAS_PCT      = 25    # Use top-N% HOF lengths to define new contig range
# AF3_EXPORT_EVERY_N: with RF3 + Boltz-2 both running, each iteration is roughly
# 6–10 h on an A100.  N=2 targets ~12–20 h between submissions, matching the
# 30 AF3 jobs/day Google-server cap.
AF3_EXPORT_EVERY_N     = 2     # Export AF3 JSON every N full iterations
AF3_TOP_N              = 30    # Designs to include per AF3 submission (max per day)
IPTM_PROMISING         = 0.55  # Boltz ipTM threshold to flag a design as "promising"
RMSD_CLIP              = 5.0   # Å beyond which RMSD contribution scores 0
PAE_MAX                = 30.0  # Å² beyond which interface PAE contribution scores 0
N_CONTACTS_NORM        = 60.0  # n_contacts above this saturates to 1.0

# Composite score weights — source-aware.  Higher-trust sources contribute more
# and use a richer formula; lower-trust sources use geometric features only.
# Trust order:  AF3  >  ESMFold2  >  Boltz-2  >  RF3-geometric.
#
# Calibration history (all vs AF3 ipTM ground truth):
#   - RF3 single-sequence ipTM/pTM were ANTI-correlated (r ≈ -0.07 / -0.51, n=18):
#     excluded from ranking entirely; RF3 contributes geometric features only.
#   - Boltz-2: weak, target-inconsistent (Spearman +0.11, top-6 50%, n=24).
#   - ESMFold2: best local predictor (binder pLDDT rho +0.34 / ipTM rho +0.28,
#     top-6 100%, n=24).  Now the primary pre-AF3 ranker.
#
# ESMFold2 weighting rationale: pLDDT correlated slightly better than ipTM, but
# both are noisy at n=24 and both are mechanistically informative (foldability +
# interface quality).  We weight them EQUALLY on purpose — leaning harder on
# pLDDT would overfit a low-n ranking difference and bias toward "foldable" over
# "binding".  Per-target selection (AF3 export) handles pLDDT's target-dependent
# magnitude, so no per-target special-casing is needed (avoids systematic bias).
#
# Each formula returns a value in [0, 1].  The RF3-only fallback is capped at
# RF3_COMPOSITE_CEILING so any model-scored entry outranks an unvalidated one.
# NOTE (2026-07 recalibration): the AF3 and ESMFold2 composites below are
# SUPERSEDED by calibrated_scoring.py and kept only as a legacy fallback for
# entries missing the calibrated inputs. calc_composite now routes:
#   AF3      -> cs.af3_binding_prior      (binder pLDDT + ApTM + interface size +
#                                          interface pLDDT + small ipTM; BpTM/loop-PAE
#                                          dropped — they don't predict binding)
#   ESMFold2 -> cs.esmfold2_stage_score   (foldability FLOOR + interface contact
#                                          density; esm_iptm NOT rewarded, to avoid
#                                          the negative-correlation selection bias)
# See calibrated_scoring.py for the evidence and the selection-bias rationale.
COMPOSITE_AF3 = {        # legacy fallback only (used if calibrated terms absent)
    "iptm":  0.45,
    "plddt": 0.25,
    "ptm":   0.15,
    "rmsd":  0.05,
    "pae":   0.10,
}
COMPOSITE_ESMFOLD2 = {   # legacy fallback only (superseded by cs.esmfold2_stage_score)
    "iptm":  0.50,
    "plddt": 0.50,
}
COMPOSITE_BOLTZ = {      # Boltz-2 scored: retained but no longer the default ranker
    "iptm":   0.45,
    "plddt":  0.25,
    "ptm":    0.15,
    "pae":    0.15,
}
COMPOSITE_RF3_GEOM = {   # RF3-only fallback: geometric features only (no ipTM/pTM/pLDDT)
    "n_contacts": 0.40,  # interface size
    "pae":        0.30,  # weak +signal in our calibration
    "rmsd":       0.30,  # backbone fidelity to RFd3 design
}
RF3_COMPOSITE_CEILING = 0.50  # cap RF3-only composite so model-scored entries dominate HOF

# ── ESMFold2 configuration (primary pre-AF3 ranker) ──────────────────────────
# ESMFold2 (Chan Zuckerberg Biohub, Rives lab) is MSA-free and faster than Boltz,
# and beat Boltz in our AF3 calibration (§5.2 of filter_methods.md).  It is
# the default local ranker.  It cannot share the foundry env (torch/CUDA pins),
# so we invoke score_with_esmfold2.py in a separate `esmfold2` conda env via
# subprocess.  Configure the env's python + the scorer path:
#     export ESMFOLD2_PYTHON=/x/capa/<user>/miniconda3/envs/esmfold2/bin/python
ESMFOLD2_ENABLE     = True
ESMFOLD2_PYTHON     = os.environ.get(
    "ESMFOLD2_PYTHON",
    str(Path.home() / "miniconda3" / "envs" / "esmfold2" / "bin" / "python"),
)
ESMFOLD2_SCRIPT     = str(_HERE / "score_with_esmfold2.py")
ESMFOLD2_TIMEOUT_S  = 3600   # whole-batch timeout (model load + N designs)
# Save the ESMFold2-predicted complex per design. The structure is already
# computed during folding, so this is ~free (disk I/O only, no extra GPU time);
# it gives an inspectable predicted complex from the actual ranker, and feeds
# hof_structures/. Off → scores only.
SAVE_ESMFOLD2_STRUCTURES = True

# Structural-Validation battery (sv_bridge). LOG-ONLY: sv_* columns are written to
# round_summary.csv for transparency + future wet-lab correlation; they never enter
# calc_composite. The occlusion filter is a mechanistic sanity gate (does the
# reactive edge reach the catalytic cleft), analogous to the foldability floor —
# not a ranking term. All opt-in; off by default (SASA/shape cost per candidate).
SV_BATTERY = False
SV_OCCLUSION_FILTER = False
SV_OCCLUSION_MIN = svb.DEFAULT_MIN_OCCLUSION
# Multi-GPU for the ESMFold2 stage (the dominant cost). Designs are embarrassingly
# parallel, so we data-parallel-shard them across GPUs, one scorer process each.
#   1      → single GPU (default; no behavior change)
#   N      → use up to N *free* GPUs (skips GPUs already busy with other jobs)
#   "auto" → use all free GPUs, capped at ESMFOLD2_MAX_GPUS
# Override per run with --esmfold2-gpus. RFd3/LMPNN stay single-GPU.
ESMFOLD2_GPUS       = 1
ESMFOLD2_MAX_GPUS   = 4
ESMFOLD2_GPU_FREE_MB = 2000   # a GPU counts as "free" if used memory < this

# ── Boltz-2 configuration (retained, OFF by default) ─────────────────────────
# Superseded by ESMFold2 as the ranker but kept available for A/B comparison.
# Boltz pins numpy<2.0 / cublas<12.5 / older lightning; install in a separate env
# and set BOLTZ_EXECUTABLE if you re-enable it.
BOLTZ_ENABLE            = False
BOLTZ_EXECUTABLE        = os.environ.get("BOLTZ_EXECUTABLE", "boltz")
BOLTZ_USE_MSA_SERVER    = True
BOLTZ_DIFFUSION_SAMPLES = 1
BOLTZ_TIMEOUT_S         = 300

# ── Pipeline stage toggles ───────────────────────────────────────────────────
# RF3 is OFF by default: ESMFold2 replaced it as the ranker (RF3's own ipTM/pTM/
# pLDDT were anti-correlated with AF3, §4.1), so running it is a full extra model
# inference per design that no longer informs selection — pure slowdown. Kept for
# curiosity: enable with --enable-rf3 (or RF3_ENABLE=True) to log its geometric
# features (n_contacts, RMSD-to-RFd3, interface PAE) into round_summary.csv.
RF3_ENABLE = False
# LMPNN already returns the designed sequence (all we need downstream). Writing a
# CIF per design is cheap but clutters a long run with thousands of files we never
# read. Off by default; the RFd3 backbone CIFs and (if enabled) RF3/AF3 structures
# remain the structural record.
SAVE_LMPNN_STRUCTURES = False


# ── Utility functions ─────────────────────────────────────────────────────────

def setup_env() -> None:
    OUT_BASE.mkdir(parents=True, exist_ok=True)
    CKPT_DIR.mkdir(parents=True, exist_ok=True)
    os.environ["FOUNDRY_CHECKPOINT_DIRS"] = str(CKPT_DIR)
    os.environ["DGLBACKEND"] = "pytorch"


def get_free_gpus(max_n: int = None, mem_threshold_mb: int = ESMFOLD2_GPU_FREE_MB) -> list:
    """
    Return indices of GPUs whose used memory is below `mem_threshold_mb` (i.e.
    free / not in use by another job), via nvidia-smi. Empty list if none/no smi
    (callers fall back to the default device). Respects other users on a shared node.
    """
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=index,memory.used", "--format=csv,noheader,nounits"]
        ).decode()
    except Exception:
        return []
    free = []
    for line in out.strip().splitlines():
        try:
            idx, used = (x.strip() for x in line.split(","))
            if int(used) < mem_threshold_mb:
                free.append(int(idx))
        except ValueError:
            continue
    return free[:max_n] if max_n else free


def resolve_esmfold2_gpus() -> list:
    """
    Turn the ESMFOLD2_GPUS setting into a concrete list of GPU ids to shard over.
    Returns [None] for the single-GPU default (no CUDA_VISIBLE_DEVICES pinning —
    unchanged behavior); otherwise a list of free GPU indices.
    """
    spec = ESMFOLD2_GPUS
    if isinstance(spec, str) and spec.lower() == "auto":
        return get_free_gpus(max_n=ESMFOLD2_MAX_GPUS) or [None]
    try:
        n = int(spec)
    except (TypeError, ValueError):
        n = 1
    if n <= 1:
        return [None]
    free = get_free_gpus(max_n=n)
    return free if len(free) > 1 else [None]


def get_seq(atom_array, chain_id: str = "A") -> str:
    """Extract one-letter sequence from an AtomArray (CA atoms only)."""
    mask = (atom_array.chain_id == chain_id) & (atom_array.atom_name == "CA")
    ca = atom_array[mask]
    if len(ca) == 0:
        return ""
    ca = ca[np.argsort(ca.res_id)]
    letters = []
    for rn in ca.res_name:
        try:
            letters.append(ProteinSequence.convert_letter_3to1(rn))
        except Exception:
            letters.append("X")
    return "".join(letters)


def renumber(atom_array):
    """Renumber residue IDs contiguously from 1 per chain (mutates in place)."""
    new_ids = np.zeros(len(atom_array), dtype=int)
    for ch in np.unique(atom_array.chain_id):
        mask = atom_array.chain_id == ch
        ids = atom_array.res_id[mask]
        uniq = list(dict.fromkeys(ids.tolist()))
        mapping = {old: new for new, old in enumerate(uniq, start=1)}
        new_ids[mask] = [mapping[r] for r in ids]
    atom_array.res_id = new_ids
    return atom_array


def pdb_chain_seq(pdb_path: str, chain_id: str) -> str:
    """Return one-letter sequence for a chain from a PDB file."""
    try:
        arr = PDBFile.read(pdb_path).get_structure()[0]
        return get_seq(arr, chain_id)
    except Exception as exc:
        logger.warning(f"Could not parse {pdb_path} for chain {chain_id}: {exc}")
        return ""


def extract_loops(sequence: str, selected_loops: list) -> dict:
    """
    Extract loop sub-sequences using flank tripeptides.
    Returns dict with loop_<name>_seq and loop_<name>_len entries.
    """
    result = {}
    cursor = 0
    for lc in selected_loops:
        name = lc["name"]
        pattern = re.compile(f"{lc['left']}([A-Z]*?){lc['right']}")
        m = pattern.search(sequence[cursor:])
        if m:
            seq = m.group(1)
            result[f"loop_{name}_seq"] = seq
            result[f"loop_{name}_len"] = len(seq)
            cursor += m.end() - len(lc["right"])
        else:
            result[f"loop_{name}_seq"] = "MISSING"
            result[f"loop_{name}_len"] = -1
    return result


def get_fixed_residues(atom_array, selected_loops: list, binder_chain: str, target_chain: str) -> list:
    """
    Identify fixed residues for LMPNN: all scaffold positions outside the loops
    (determined via flank tripeptide search) plus every residue on the target chain.
    """
    seq = get_seq(atom_array, binder_chain)
    fixed = []
    current_start = 1
    cursor = 0

    for lc in selected_loops:
        pattern = re.compile(f"{lc['left']}([A-Z]*?){lc['right']}")
        m = pattern.search(seq[cursor:])
        if m:
            abs_start = cursor + m.start()
            loop_1idx = abs_start + len(lc["left"]) + 1
            loop_len = len(m.group(1))
            fixed.extend(range(current_start, loop_1idx))
            current_start = loop_1idx + loop_len
            cursor = cursor + m.end() - len(lc["right"])

    fixed.extend(range(current_start, len(seq) + 1))
    fixed_residues = [f"{binder_chain}{r}" for r in fixed]

    # Pin entire target chain
    tgt_mask = (atom_array.chain_id == target_chain) & (atom_array.atom_name == "CA")
    for r in np.unique(atom_array.res_id[tgt_mask]):
        fixed_residues.append(f"{target_chain}{r}")

    return fixed_residues


def count_interface_contacts(atom_array, chain_a: str = "A", chain_b: str = "B",
                              cutoff: float = 5.0) -> int:
    """
    Count unique residues in chain_a with any heavy atom within cutoff Å of chain_b.
    Proxy for buried interface area.
    """
    arr_a = atom_array[atom_array.chain_id == chain_a]
    arr_b = atom_array[atom_array.chain_id == chain_b]
    if len(arr_a) == 0 or len(arr_b) == 0:
        return 0
    diff = arr_a.coord[:, None, :] - arr_b.coord[None, :, :]   # (N, M, 3)
    dist = np.sqrt((diff ** 2).sum(axis=-1))                    # (N, M)
    contact_mask = (dist < cutoff).any(axis=1)
    return int(len(np.unique(arr_a.res_id[contact_mask])))


def _esm_iface_feats(cif_path) -> dict:
    """
    Interface geometry from an ESMFold2 predicted-complex CIF, keyed with the
    esm_iface_* names calc_composite reads. Returns {} if the CIF is missing/
    unreadable (composite then falls back to the pLDDT foldability floor only).
    Uses the SAME extractor as the 2026-07 calibration (calibrated_scoring).
    """
    if not cif_path or not Path(cif_path).exists():
        return {}
    try:
        f = cs.interface_features_from_cif(Path(cif_path).read_text(),
                                           binder_chain=DESIGN_BINDER_CHAIN,
                                           target_chain=DESIGN_TARGET_CHAIN)
    except Exception:
        return {}
    return {"esm_iface_n_iface_res":     f.get("n_iface_res"),
            "esm_iface_contact_density": f.get("contact_density"),
            "esm_iface_iface_plddt":     f.get("iface_plddt")}


def _safe(x, default=0.0):
    """Coerce to float, treating None / NaN / missing as `default`."""
    try:
        v = float(x)
        return default if np.isnan(v) else v
    except (TypeError, ValueError):
        return default


def _rmsd_score(rmsd_val):
    if rmsd_val is None or np.isnan(_safe(rmsd_val, np.nan)):
        return 0.0
    return max(0.0, 1.0 - _safe(rmsd_val) / RMSD_CLIP)


def _pae_score(iface_pae):
    if iface_pae is None or np.isnan(_safe(iface_pae, np.nan)):
        return 0.0
    return max(0.0, 1.0 - _safe(iface_pae) / PAE_MAX)


def _contacts_score(n_contacts):
    return min(1.0, _safe(n_contacts) / N_CONTACTS_NORM)


def _normalize_plddt(val) -> float:
    """Coerce pLDDT to the 0–100 scale (Boltz/RF3 return 0–1, AF3 returns 0–100)."""
    v = _safe(val, np.nan)
    if np.isnan(v):
        return 0.0
    return v * 100.0 if 0.0 < v <= 1.0 else v


def calc_composite(entry: dict) -> float:
    """
    Source-aware composite score in [0, 1].  Selects metrics based on the
    highest-trust prediction available for the entry:

        AF3-validated   →  COMPOSITE_AF3      (full trust)
        ESMFold2-scored →  COMPOSITE_ESMFOLD2 (primary local ranker: ipTM + pLDDT)
        Boltz-scored    →  COMPOSITE_BOLTZ    (retained; off by default)
        RF3-only        →  COMPOSITE_RF3_GEOM (geometric features only, capped)

    The RF3-only branch deliberately excludes RF3's iptm/pTM/pLDDT because they
    are anti-correlated with AF3 ipTM for de novo binders (calibrated empirically).
    It is capped at RF3_COMPOSITE_CEILING so any model-scored entry dominates.
    """
    src = entry.get("source", "RF3")

    # ── AF3-validated → calibrated binding prior ──
    # Positive foldability/interface terms only (binder pLDDT, ApTM, interface
    # size + pLDDT, small ipTM). BpTM and loop-PAE are NOT used (calibration:
    # BpTM is target-determined & didn't replicate; loop-PAE ~0). Falls back to
    # the legacy blend only if no calibrated term is present.
    if src == "AF3":
        prior = cs.af3_binding_prior({
            "af3_plddt":               _normalize_plddt(entry.get("plddt")),  # binder-chain pLDDT (0-100)
            "af3_aptm":                entry.get("af3_aptm"),
            "af3_iptm":                entry.get("iptm"),
            "af3_iface_n_iface_res":   entry.get("af3_iface_n_iface_res"),
            "af3_iface_iface_plddt":   entry.get("af3_iface_iface_plddt"),
        })
        if prior == prior:   # not NaN
            return prior
        w = COMPOSITE_AF3    # legacy fallback
        return (
            w["iptm"]  * _safe(entry.get("iptm"))
            + w["plddt"] * (_normalize_plddt(entry.get("plddt")) / 100.0)
            + w["ptm"]   * _safe(entry.get("ptm"))
            + w["rmsd"]  * _rmsd_score(entry.get("rmsd_to_rfd3"))
            + w["pae"]   * _pae_score(entry.get("interface_pae"))
        )

    # ── ESMFold2-scored → foldability FILTER + contact-density prior ──
    # esm_iptm/esm_lplddt are NOT rewarded: their negative in-sample correlation
    # with binding is a selection-bias artifact (we never observe the non-folders
    # that would score even lower). We gate on esm_plddt (binding-neutral) and rank
    # by interface contact density (the one positive, mechanistic ESMFold2 feature).
    if entry.get("esm_iptm") is not None or entry.get("esm_plddt") is not None:
        return cs.esmfold2_stage_score({
            "esm_plddt":                 _normalize_plddt(entry.get("esm_plddt")),
            "esm_iface_contact_density": entry.get("esm_iface_contact_density"),
            "esm_iface_n_iface_res":     entry.get("esm_iface_n_iface_res"),
        })

    # ── Boltz-2-scored ──
    if entry.get("boltz_iptm") is not None:
        w = COMPOSITE_BOLTZ
        return (
            w["iptm"]  * _safe(entry.get("boltz_iptm"))
            + w["plddt"] * (_normalize_plddt(entry.get("boltz_plddt")) / 100.0)
            + w["ptm"]   * _safe(entry.get("boltz_ptm"))
            + w["pae"]   * _pae_score(entry.get("boltz_iface_pae"))
        )

    # ── RF3-only fallback (geometric features, capped) ──
    w = COMPOSITE_RF3_GEOM
    raw = (
        w["n_contacts"] * _contacts_score(entry.get("n_contacts"))
        + w["pae"]        * _pae_score(entry.get("interface_pae"))
        + w["rmsd"]       * _rmsd_score(entry.get("rmsd_to_rfd3"))
    )
    return min(raw, RF3_COMPOSITE_CEILING)


def _parse_boltz_output(run_dir: Path, job_name: str, binder_len: int) -> dict:
    """
    Pull metrics from a Boltz-2 prediction directory.

    Boltz 2.2.x writes to <run_dir>/boltz_results_<job>/predictions/<job>/ with
    confidence_<job>_model_0.json, pae_<job>_model_0.npz, plddt_<job>_model_0.npz,
    and <job>_model_0.cif.  Returns {iptm, ptm, plddt, iface_pae, cif_path};
    missing files are silently skipped.
    """
    out: dict = {}
    job_dir = run_dir / f"boltz_results_{job_name}" / "predictions" / job_name
    if not job_dir.exists():
        return out

    # ── Confidence JSON (iptm, ptm, summary plddt) ──
    conf_path = job_dir / f"confidence_{job_name}_model_0.json"
    if conf_path.exists():
        try:
            with open(conf_path) as f:
                conf = json.load(f)
            out["iptm"]  = conf.get("iptm", conf.get("complex_iptm", conf.get("protein_iptm")))
            out["ptm"]   = conf.get("ptm",  conf.get("complex_ptm"))
            out["plddt"] = _normalize_plddt(conf.get("complex_plddt", conf.get("plddt")))
        except Exception as e:
            logger.warning(f"Failed to read Boltz confidence: {e}")

    # ── PAE matrix → interface PAE (binder ↔ target sub-blocks) ──
    pae_path = job_dir / f"pae_{job_name}_model_0.npz"
    if pae_path.exists():
        try:
            arr = np.load(pae_path, allow_pickle=False)
            pae_mat = arr["pae"] if "pae" in arr.files else arr[arr.files[0]]
            if pae_mat.ndim == 2 and pae_mat.shape[0] > binder_len:
                out["iface_pae"] = float(
                    (pae_mat[:binder_len, binder_len:].mean()
                     + pae_mat[binder_len:, :binder_len].mean()) / 2
                )
        except Exception as e:
            logger.warning(f"Failed to read Boltz PAE: {e}")

    # ── Binder-only pLDDT (overrides complex pLDDT from confidence JSON) ──
    # plddt array is per-residue with binder residues first (binder is chain A in
    # the YAML we wrote), so the first binder_len entries are the binder.
    plddt_path = job_dir / f"plddt_{job_name}_model_0.npz"
    if plddt_path.exists():
        try:
            arr = np.load(plddt_path, allow_pickle=False)
            plddt_arr = (arr["plddt"] if "plddt" in arr.files else arr[arr.files[0]]).flatten()
            if len(plddt_arr) >= binder_len:
                out["plddt"] = _normalize_plddt(float(plddt_arr[:binder_len].mean()))
        except Exception as e:
            logger.warning(f"Failed to read Boltz pLDDT: {e}")

    # ── Structure path ──
    cif_path = job_dir / f"{job_name}_model_0.cif"
    if cif_path.exists():
        out["cif_path"] = str(cif_path)

    return out


def backbone_rmsd(ref_array, query_array) -> float:
    """RMSD between backbone atoms of two atom arrays after superimposition."""
    bb_names = set(PROTEIN_BACKBONE_ATOM_NAMES) - {"OXT"}
    ref_bb  = ref_array[np.isin(ref_array.atom_name, list(bb_names))]
    qry_bb  = query_array[np.isin(query_array.atom_name, list(bb_names))]
    n = min(len(ref_bb), len(qry_bb))
    if n == 0:
        return float("nan")
    qry_fit, _ = superimpose(ref_bb[:n], qry_bb[:n])
    return float(bio_rmsd(ref_bb[:n], qry_fit))


# ── Main class ────────────────────────────────────────────────────────────────

class IterativeRefiner:
    """
    Runs the RFd3 → LigandMPNN → RF3 (complex) iterative loop for TIMP3 binders.

    Parameters
    ----------
    active_targets : list[str]
        Subset of TARGETS keys to include.
    active_loops : list[str]
        Subset of LOOP_CONFIGS keys to redesign (e.g. ["AB", "C", "EF"]).
    state_path : Path | None
        Where to persist JSON state.  Defaults to OUT_BASE/refinement_state.json.
    """

    def __init__(self, active_targets: list, active_loops: list, state_path: Path = None):
        setup_env()
        self.active_targets = active_targets
        self.active_loops   = active_loops
        self.selected_loops = sorted(
            [{**LOOP_CONFIGS[n], "name": n} for n in active_loops],
            key=lambda x: x["pos"]
        )
        self.state_path = Path(state_path) if state_path else OUT_BASE / "refinement_state.json"
        self.hof_struct_dir = OUT_BASE / "hof_structures"

        # Pre-load target sequences from crystal PDB files
        self.target_seqs: dict[str, str] = {}
        for tname, tcfg in TARGETS.items():
            if tname not in active_targets:
                continue
            pdb_path = DATA_DIR / tcfg["pdb"]
            if pdb_path.exists():
                self.target_seqs[tname] = pdb_chain_seq(str(pdb_path), tcfg["target_chain"])
                logger.info(f"Loaded {tname} target sequence: {len(self.target_seqs[tname])} residues")
            else:
                logger.warning(f"PDB not found for {tname}: {pdb_path}")

        self._load_state()

    # ── State management ──────────────────────────────────────────────────────

    def _load_state(self) -> None:
        if self.state_path.exists():
            try:
                with open(self.state_path) as f:
                    self.state = json.load(f)
                logger.info(f"Resumed from iteration {self.state['iteration']} "
                            f"(T={self.state['temperature']:.3f})")
                return
            except Exception as exc:
                logger.error(f"State load failed ({exc}), starting fresh.")
        self._init_state()

    def _init_state(self) -> None:
        self.state = {
            "iteration":   0,
            "temperature": INIT_TEMPERATURE,
            "hof":         {t: [] for t in self.active_targets},  # per-target HOF
            "last_af3_it": -1,
        }
        self._save_state()

    def _save_state(self) -> None:
        self.state_path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.state_path, "w") as f:
            json.dump(self.state, f, indent=2)

    # ── Adaptive contig bias ──────────────────────────────────────────────────

    def _adaptive_loop_ranges(self) -> dict:
        """
        From the top ADAPTIVE_BIAS_PCT% of each target's HOF, compute the
        length range (min, max) actually used by high-scoring designs.
        Returns dict {loop_name: (min_len, max_len)} or empty if not enough data.
        """
        it = self.state["iteration"]
        if it < ADAPTIVE_BIAS_START:
            return {}

        all_hof = []
        for t in self.active_targets:
            all_hof.extend(self.state["hof"].get(t, []))

        if len(all_hof) < 20:
            return {}

        all_hof_sorted = sorted(all_hof, key=lambda x: x.get("composite_score", 0), reverse=True)
        top_n = max(5, len(all_hof_sorted) * ADAPTIVE_BIAS_PCT // 100)
        top = all_hof_sorted[:top_n]

        ranges = {}
        for lc in self.selected_loops:
            name = lc["name"]
            key  = f"loop_{name}_len"
            lens = [e[key] for e in top if e.get(key, -1) > 0]
            if len(lens) >= 3:
                lo = max(lc["normal"], min(lens))
                hi = min(lc["max"],    max(lens))
                if lo <= hi:
                    ranges[name] = (lo, hi)
                    logger.info(f"Adaptive bias: loop {name} → [{lo}, {hi}] "
                                f"(from {len(lens)} top designs)")
        return ranges

    # NOTE: a target-only template CIF extractor (_ensure_target_template) lived
    # here briefly to feed RF3 a template.  This RF3 build's SequenceComponent
    # doesn't accept the AF3-style `useStructureTemplate` field, so the helper
    # was removed.  When the correct RF3 template API is identified, re-add an
    # extractor using PDBFile.read + chain filter + renumber + to_cif_file.

    # ── Boltz-2 (local AF3-class ranker) ──────────────────────────────────────

    def run_boltz(self, target_name: str, candidates: list, out_dir: Path) -> list:
        """
        Score each candidate with Boltz-2.  Boltz-2 is an open-source AF3-class
        model used here as the local pre-AF3 ranker — much better correlated
        with AF3 ipTM than RF3 single-sequence (calibrated empirically).

        Adds boltz_iptm / boltz_ptm / boltz_plddt / boltz_iface_pae fields to
        each candidate and recomputes composite_score (the Boltz-2 affinity
        head only applies to protein-ligand and is not invoked here).
        Returns the (in-place mutated) candidates list.

        If `boltz` is not installed on the system, logs a warning and returns
        the candidates unchanged.  Set BOLTZ_ENABLE = False to skip explicitly.
        """
        if not BOLTZ_ENABLE or not candidates:
            return candidates

        out_dir.mkdir(parents=True, exist_ok=True)
        bc = DESIGN_BINDER_CHAIN
        fc = DESIGN_TARGET_CHAIN
        t0 = time.time()
        n_done = 0

        for cand in candidates:
            did  = cand["design_id"]
            bseq = cand.get("full_seq", "")
            tseq = cand.get("target_seq") or self.target_seqs.get(target_name, "")
            if not bseq or not tseq:
                continue

            # Boltz YAML input (one job per design).
            # NOTE: Boltz-2's affinity head is protein-LIGAND only (verified
            # 2.2.1 raises "Chain A is not a ligand!" on a protein binder).
            # We rely on its structure-confidence metrics (ipTM, pTM, pLDDT,
            # PAE), which are still AF3-class for ranking.
            yaml_path = out_dir / f"{did}.yaml"
            yaml_path.write_text(
                "version: 1\n"
                "sequences:\n"
                f"  - protein:\n"
                f"      id: {bc}\n"
                f"      sequence: {bseq}\n"
                f"  - protein:\n"
                f"      id: {fc}\n"
                f"      sequence: {tseq}\n"
            )

            run_dir = out_dir / f"{did}_run"
            cmd = [
                BOLTZ_EXECUTABLE, "predict", str(yaml_path),
                "--out_dir", str(run_dir),
                "--diffusion_samples", str(BOLTZ_DIFFUSION_SAMPLES),
                "--output_format", "mmcif",
            ]
            if BOLTZ_USE_MSA_SERVER:
                cmd.append("--use_msa_server")

            try:
                # Boltz runs in its own conda env (different numpy/cublas/lightning).
                # Wipe PYTHONPATH so the subprocess sees only its own site-packages.
                child_env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
                proc = subprocess.run(
                    cmd, capture_output=True, text=True,
                    timeout=BOLTZ_TIMEOUT_S, env=child_env,
                )
                if proc.returncode != 0:
                    logger.warning(f"Boltz failed on {did}: {proc.stderr[-300:].strip()}")
                    continue
            except FileNotFoundError:
                logger.error(
                    f"Boltz binary not found at: {BOLTZ_EXECUTABLE}\n"
                    "  Install Boltz in a separate conda env (NOT foundry — version conflicts):\n"
                    "    conda create -n boltz python=3.12 -y\n"
                    "    conda activate boltz && pip install boltz\n"
                    "  Then set BOLTZ_EXECUTABLE in iterative_refinement.py to the path of\n"
                    "  that env's `boltz` binary (e.g. ~/miniconda3/envs/boltz/bin/boltz),\n"
                    "  or export BOLTZ_EXECUTABLE=... before running.\n"
                    "  Or set BOLTZ_ENABLE=False to skip Boltz entirely."
                )
                return candidates
            except subprocess.TimeoutExpired:
                logger.warning(f"Boltz timed out on {did} (>{BOLTZ_TIMEOUT_S}s)")
                continue

            metrics = _parse_boltz_output(run_dir, did, len(bseq))
            if metrics:
                boltz_iptm = metrics.get("iptm")
                cand.update({
                    "boltz_iptm":      boltz_iptm,
                    "boltz_ptm":       metrics.get("ptm"),
                    "boltz_plddt":     metrics.get("plddt"),
                    "boltz_iface_pae": metrics.get("iface_pae"),
                    "boltz_cif":       metrics.get("cif_path"),
                    "source":          "Boltz",  # promote from "RF3"
                    "promising":       _safe(boltz_iptm) >= IPTM_PROMISING,
                })
                cand["composite_score"] = calc_composite(cand)
                n_done += 1

        logger.info(
            f"[{target_name}] Boltz-2 done in {(time.time()-t0)/60:.1f} min "
            f"({n_done}/{len(candidates)} scored)"
        )
        return candidates

    # ── ESMFold2 (primary local ranker) ──────────────────────────────────────

    def _esmfold2_one_shard(self, rows: list, out_dir: Path, tag: str, gpu) -> Path:
        """
        Launch one ESMFold2 scorer subprocess over `rows`. If `gpu` is an int,
        pin it via CUDA_VISIBLE_DEVICES (the process then sees it as cuda:0).
        Returns the Popen handle + its output CSV path (started, not awaited).
        """
        in_csv  = out_dir / f"esm_input{tag}.csv"
        out_csv = out_dir / f"esm_scores{tag}.csv"
        pd.DataFrame(rows).to_csv(in_csv, index=False)

        cmd = [ESMFOLD2_PYTHON, ESMFOLD2_SCRIPT, "--input", str(in_csv), "--out", str(out_csv)]
        if SAVE_ESMFOLD2_STRUCTURES:
            cmd += ["--cif-dir", str(out_dir / "structures")]

        child_env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
        if gpu is not None:
            child_env["CUDA_VISIBLE_DEVICES"] = str(gpu)
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                text=True, env=child_env)
        return proc, out_csv

    def _score_sequences_esmfold2(self, rows: list, out_dir: Path) -> dict:
        """
        Core ESMFold2 invocation, reused by run_esmfold2 and by the specificity
        subclass (which scores each binder against on- AND off-target).

        `rows`: list of {design_id, target_name, full_seq, target_seq}; design_id
        must be unique per row (specificity suffixes ::on / ::off).
        Data-parallel across free GPUs when ESMFOLD2_GPUS > 1 (designs sharded,
        one scorer process per GPU). Returns {design_id: {esm_iptm, esm_ptm,
        esm_plddt[, esm_cif]}} (empty on total failure).
        """
        if not ESMFOLD2_ENABLE or not rows:
            return {}
        out_dir.mkdir(parents=True, exist_ok=True)

        gpus = resolve_esmfold2_gpus()
        n_shards = min(len(gpus), len(rows))
        # Round-robin rows into shards (keeps each shard balanced)
        shards = [rows[i::n_shards] for i in range(n_shards)]
        if n_shards > 1:
            logger.info(f"ESMFold2: sharding {len(rows)} jobs across GPUs {gpus[:n_shards]}")

        procs = []
        for i, shard in enumerate(shards):
            tag = "" if n_shards == 1 else f"_g{i}"
            procs.append(self._esmfold2_one_shard(shard, out_dir, tag, gpus[i]))

        merged: dict = {}
        for proc, out_csv in procs:
            try:
                _, stderr = proc.communicate(timeout=ESMFOLD2_TIMEOUT_S)
                if proc.returncode != 0:
                    logger.warning(f"ESMFold2 shard failed: {stderr[-400:].strip()}")
                    continue
            except subprocess.TimeoutExpired:
                proc.kill()
                logger.warning(f"ESMFold2 shard timed out (>{ESMFOLD2_TIMEOUT_S}s)")
                continue
            except FileNotFoundError:
                logger.error(f"ESMFold2 python not found at: {ESMFOLD2_PYTHON}\n"
                             "  Set ESMFOLD2_PYTHON or ESMFOLD2_ENABLE=False.")
                return {}
            if not out_csv.exists():
                continue
            df = pd.read_csv(out_csv)
            for _, r in df.iterrows():
                rec = {"esm_iptm": float(r.get("esm_iptm", float("nan"))),
                       "esm_ptm":   float(r.get("esm_ptm", float("nan"))),
                       "esm_plddt": float(r.get("esm_plddt", float("nan")))}
                if "esm_cif" in r and isinstance(r["esm_cif"], str):
                    rec["esm_cif"] = r["esm_cif"]
                merged[r["design_id"]] = rec
        if not merged:
            logger.warning("ESMFold2 produced no scores.")
        else:
            # Consolidate the per-GPU shards (esm_scores_g*.csv) into one
            # authoritative CSV per stage dir; the merged values also flow into
            # round_summary.csv at end of iteration.
            pd.DataFrame([{"design_id": did, **m} for did, m in merged.items()]).to_csv(
                out_dir / "esm_scores.csv", index=False)
        return merged

    def run_esmfold2(self, target_name: str, candidates: list, out_dir: Path) -> list:
        """
        Score candidates with ESMFold2 against their design target, merge
        esm_iptm / esm_plddt / esm_ptm, set source="ESMFold2", recompute the
        composite, and flag `promising` on esm_iptm.  No-op if ESMFOLD2_ENABLE
        is False or the env is missing (pipeline keeps RF3 geometric composites).
        """
        if not ESMFOLD2_ENABLE or not candidates:
            return candidates

        rows = []
        for c in candidates:
            tseq = c.get("target_seq") or self.target_seqs.get(target_name, "")
            if c.get("full_seq") and tseq:
                rows.append({"design_id": c["design_id"], "target_name": target_name,
                             "full_seq": c["full_seq"], "target_seq": tseq})

        t0 = time.time()
        scores = self._score_sequences_esmfold2(rows, out_dir)
        n_done = 0
        for c in candidates:
            m = scores.get(c["design_id"])
            if not m:
                continue
            c.update({**m, "source": "ESMFold2",
                      "promising": _safe(m["esm_iptm"]) >= IPTM_PROMISING})
            # Interface geometry from the ESMFold2 predicted complex (the one
            # positive ESMFold2 binding-ish feature is contact density). Free when
            # SAVE_ESMFOLD2_STRUCTURES wrote the CIF; skipped otherwise.
            c.update(_esm_iface_feats(m.get("esm_cif")))
            # Full Structural-Validation battery — LOG-ONLY, never enters the
            # composite (calibration: these don't predict binding). Optional
            # occlusion gate can drop designs that miss the catalytic cleft.
            if SV_BATTERY:
                sv = svb.sv_battery(m.get("esm_cif"), DESIGN_BINDER_CHAIN,
                                    DESIGN_TARGET_CHAIN, target_id=target_name)
                c.update(sv)
                if SV_OCCLUSION_FILTER:
                    c["sv_occlusion_pass"] = svb.occlusion_pass(sv, SV_OCCLUSION_MIN)
            c["composite_score"] = calc_composite(c)
            if (SV_BATTERY and SV_OCCLUSION_FILTER
                    and c.get("sv_occlusion_pass") is False):
                c["composite_score"] = 0.0        # sanity gate: exclude from HOF/ranking
                c["sv_filtered"] = "occlusion"
            n_done += 1

        logger.info(f"[{target_name}] ESMFold2 done in {(time.time()-t0)/60:.1f} min "
                    f"({n_done}/{len(candidates)} scored)")
        return candidates

    # ── RFd3 ─────────────────────────────────────────────────────────────────

    def _build_contig(self, tcfg: dict, fix_chain_len: int,
                      adaptive_ranges: dict = None,
                      scaffold_len: int = None) -> tuple[str, str]:
        """
        Build contig string and length range for RFd3.
        Returns (contig_string, length_range_string).

        scaffold_len is the binder length, derived from the template at runtime
        (run_rfd3) so it tracks whatever TIMP3 form the template holds — 121 (N-domain)
        or 188 (full length). Falls back to tcfg["scaffold_len"] if not supplied.
        """
        bc     = tcfg["binder_chain"]
        fc     = tcfg["target_chain"]
        total  = scaffold_len if scaffold_len is not None else tcfg.get("scaffold_len")
        ar     = adaptive_ranges or {}

        parts = []
        current = 1
        for lc in self.selected_loops:
            if current <= lc["pos"]:
                parts.append(f"{bc}{current}-{lc['pos']}")
            name = lc["name"]
            lo = ar.get(name, (lc["normal"], lc["max"]))[0] if name in ar else lc["normal"]
            hi = ar.get(name, (lc["normal"], lc["max"]))[1] if name in ar else lc["max"]
            parts.append(f"{lo}-{hi}")
            current = lc["pos"] + lc["normal"] + 1

        if current <= total:
            parts.append(f"{bc}{current}-{total}")

        binder_contig = ",".join(parts)
        full_contig   = f"{binder_contig},/0,{fc}1-{fix_chain_len}"

        # Compute output length range
        base = total + fix_chain_len - sum(lc["normal"] for lc in self.selected_loops)
        min_len = base + sum(
            (ar.get(lc["name"], (lc["normal"], lc["max"]))[0] if lc["name"] in ar else lc["normal"])
            for lc in self.selected_loops
        )
        max_len = base + sum(
            (ar.get(lc["name"], (lc["normal"], lc["max"]))[1] if lc["name"] in ar else lc["max"])
            for lc in self.selected_loops
        )
        return full_contig, f"{min_len}-{max_len}"

    def run_rfd3(self, target_name: str, pdb_path: str, n_designs: int,
                 adaptive_ranges: dict = None) -> list:
        """Generate backbone atom arrays using RFd3. Returns list of AtomArrays."""
        tcfg     = TARGETS[target_name]
        pdb_arr  = PDBFile.read(pdb_path).get_structure()[0]
        fix_len  = len(np.unique(pdb_arr.res_id[pdb_arr.chain_id == tcfg["target_chain"]]))
        # Derive the scaffold (binder) length from the template so 121 vs 188 is
        # never hard-coded; warn if it disagrees with the TARGETS hint.
        scaffold_len = int(len(np.unique(pdb_arr.res_id[pdb_arr.chain_id == tcfg["binder_chain"]])))
        hint = tcfg.get("scaffold_len")
        if hint is not None and hint != scaffold_len:
            logger.warning(f"[{target_name}] template binder length {scaffold_len} "
                           f"!= TARGETS scaffold_len {hint}; using template length.")

        contig, length_range = self._build_contig(tcfg, fix_len, adaptive_ranges,
                                                  scaffold_len=scaffold_len)
        logger.info(f"[{target_name}] RFd3 scaffold_len (from template): {scaffold_len}")
        logger.info(f"[{target_name}] RFd3 contig: {contig}")
        logger.info(f"[{target_name}] Length range: {length_range}")

        rfd3_cfg = RFD3InferenceConfig(
            diffusion_batch_size=min(10, n_designs),
            low_memory_mode=False,
            specification={"length": length_range, "contig": contig, "extra": {}},
        )
        engine = RFD3InferenceEngine(**dataclasses.asdict(rfd3_cfg))
        spec   = DesignInputSpecification(
            input=pdb_path, contig=contig, length=length_range, extra={}
        )
        n_batches = (n_designs + rfd3_cfg.diffusion_batch_size - 1) // rfd3_cfg.diffusion_batch_size

        t0 = time.time()
        outputs = engine.run(inputs=spec, n_batches=n_batches, out_dir=None)
        logger.info(f"[{target_name}] RFd3 done in {(time.time()-t0)/60:.1f} min")

        arrays = []
        if outputs:
            idx = 0
            for key, out_list in outputs.items():
                if not key.startswith("backbone"):
                    continue
                for out in out_list:
                    if idx >= n_designs:
                        break
                    arrays.append(renumber(out.atom_array))
                    idx += 1

        del engine
        torch.cuda.empty_cache()
        return arrays

    # ── LigandMPNN ────────────────────────────────────────────────────────────

    def run_lmpnn(self, target_name: str, backbones: list, out_dir: Path,
                  temperature: float) -> list:
        """
        Design sequences for each backbone using LigandMPNN.
        Returns list of dicts with keys: design_id, array, rfd3_array, full_seq,
        target_seq, seq_recovery, loop_data.

        Reads chains by the RFd3 OUTPUT convention (A = designed binder,
        B = fixed target), NOT the source-PDB chain IDs in TARGETS.
        """
        bc     = DESIGN_BINDER_CHAIN
        fc     = DESIGN_TARGET_CHAIN
        out_dir.mkdir(parents=True, exist_ok=True)

        engine  = MPNNInferenceEngine(
            model_type="ligand_mpnn", is_legacy_weights=True,
            write_structures=False, write_fasta=True, out_directory=str(out_dir)
        )
        results = []
        t0 = time.time()

        for idx, rfd3_array in enumerate(backbones):
            design_id = f"{target_name}_it{self.state['iteration']}_d{idx}"
            try:
                fixed_res = get_fixed_residues(rfd3_array, self.selected_loops, bc, fc)
                mpnn_in = {
                    "name":           design_id,
                    "batch_size":     LMPNN_SEQS_PER_BACKBONE,  # sequences sampled per backbone
                    "remove_waters":  True,
                    "seed":           42 + idx,
                    "fixed_residues": fixed_res,
                    "sampling_temp":  temperature,
                }
                mpnn_outs = engine.run(input_dicts=[mpnn_in], atom_arrays=[rfd3_array])

                for si, mp_out in enumerate(mpnn_outs):
                    valid = ~np.isnan(mp_out.atom_array.coord[:, 0])
                    arr   = renumber(mp_out.atom_array[valid])
                    bseq  = get_seq(arr, bc)
                    tseq  = get_seq(arr, fc)
                    loops = extract_loops(bseq, self.selected_loops)
                    seq_rec = float(
                        getattr(mp_out, "output_dict", {}).get("sequence_recovery", 0.0)
                        or getattr(mp_out, "seq_recovery", 0.0)
                    )
                    if SAVE_LMPNN_STRUCTURES:
                        to_cif_file(arr, str(out_dir / f"{design_id}_s{si}.cif"), file_type="cif")
                    results.append({
                        "design_id":   f"{design_id}_s{si}",
                        "target_name": target_name,
                        # rfd3_array kept only for RF3 backbone-RMSD; "array" (the
                        # LMPNN structure) is never read downstream, so not stored.
                        "rfd3_array":  rfd3_array,
                        "full_seq":    bseq,
                        "target_seq":  tseq,
                        "seq_recovery": seq_rec,
                        **loops,
                    })
            except Exception as exc:
                logger.error(f"LMPNN error on {design_id}: {exc}")

        logger.info(f"[{target_name}] LMPNN done in {(time.time()-t0)/60:.1f} min "
                    f"({len(results)} sequences)")
        del engine
        torch.cuda.empty_cache()
        return results

    # ── Candidate finalization ────────────────────────────────────────────────

    def _finalize_record(self, cand: dict, extra: dict = None) -> dict:
        """
        Convert an in-memory candidate (carrying heavy AtomArrays) into a clean,
        JSON-serializable scored record for the HOF/state. Strips the arrays
        (otherwise json.dump on the state file would fail), stamps iteration and
        temperature, and computes the composite. `source`/`promising` default to
        an unscored placeholder and are overwritten by whichever ranker runs.
        """
        rec = {k: v for k, v in cand.items() if k not in ("array", "rfd3_array")}
        if extra:
            rec.update(extra)
        rec["iteration"]   = self.state["iteration"]
        rec["temperature"] = self.state["temperature"]
        rec.setdefault("source", "unscored")
        rec.setdefault("promising", False)
        rec["composite_score"] = calc_composite(rec)
        return rec

    # ── RF3 complex scoring ───────────────────────────────────────────────────

    def run_rf3_complex(self, target_name: str, candidates: list,
                        out_dir: Path) -> list:
        """
        Fold each candidate as a two-chain complex (binder + target) using RF3
        with sequence-only input.  Extracts pLDDT, ptm, ipTM, interface PAE,
        backbone RMSD, and interface contacts.  Saves CIF + metrics JSON.
        Returns list of scored metadata dicts (no atom arrays — too large to cache).

        OFF by default (RF3_ENABLE): RF3 no longer informs ranking (ESMFold2 does),
        so when disabled we just finalize the candidates (strip arrays, stamp
        bookkeeping) and let ESMFold2 score them.
        """
        if not candidates:
            return []

        if not RF3_ENABLE:
            return [self._finalize_record(c) for c in candidates]

        # Use the RFd3 OUTPUT chain convention (binder=A, target=B); the source-PDB
        # chain IDs in TARGETS apply only to reading target_seqs and building contigs.
        bc   = DESIGN_BINDER_CHAIN
        fc   = DESIGN_TARGET_CHAIN
        out_dir.mkdir(parents=True, exist_ok=True)

        target_seq = self.target_seqs.get(target_name, "")
        if not target_seq:
            logger.warning(f"No target sequence for {target_name}; RF3 complex skipped.")
            return []

        # NOTE: RF3 runs single-sequence here.  An earlier attempt to pass
        # target structure as an AF3-style template was rejected by this RF3
        # build (`SequenceComponent` doesn't accept `useStructureTemplate`).
        # Boltz-2 is the trusted local ranker; RF3 contributes only geometric
        # features (n_contacts, RMSD, iface PAE) until the right template API
        # is identified.
        precision = "bf16-mixed"
        if torch.cuda.is_available():
            device_name = torch.cuda.get_device_name(0)
            if "V100" in device_name:
                precision = "16-mixed"
        engine  = RF3InferenceEngine(ckpt_path="rf3", verbose=False, trainer_overrides={"precision": precision})
        scored  = []
        t0      = time.time()

        for cand in candidates:
            did   = cand["design_id"]
            bseq  = cand["full_seq"]
            if not bseq:
                continue
            try:
                rf3_data = {
                    "name": did,
                    "components": [
                        {"id": bc, "sequence": bseq},
                        {"id": fc, "sequence": target_seq},
                    ],
                }
                rf3_input = InferenceInput.from_json_dict(rf3_data)
                rf3_outs  = engine.run(inputs=rf3_input, annotate_b_factor_with_plddt=True)

                rf3_key  = next(iter(rf3_outs))
                rf3_out  = rf3_outs[rf3_key][0]
                rf3_arr  = renumber(rf3_out.atom_array)

                # Confidence metrics.  Normalize pLDDT to 0–100 so it lives on
                # the same scale as AF3's atom_plddts (RF3 returns 0–1 here).
                conf  = rf3_out.summary_confidences or {}
                plddt = float(conf.get("overall_plddt", conf.get("plddt", 0.0)))
                if 0.0 < plddt <= 1.0:
                    plddt *= 100.0
                ptm   = conf.get("ptm", 0.0)
                iptm  = (conf.get("iptm")
                         or conf.get("ipTM")
                         or conf.get("interface_ptm")
                         or conf.get("complex_pTM")
                         or 0.0)

                # Interface PAE (binder ↔ target sub-blocks)
                iface_pae = float("nan")
                if hasattr(rf3_out, "confidences") and rf3_out.confidences:
                    raw_pae = rf3_out.confidences.get("pae")
                    if raw_pae is not None:
                        pae_mat = np.array(raw_pae)
                        nb = len(bseq)
                        if pae_mat.shape[0] > nb:
                            iface_pae = float(
                                (np.mean(pae_mat[:nb, nb:]) + np.mean(pae_mat[nb:, :nb])) / 2
                            )

                # Backbone RMSD: RF3 prediction vs LMPNN input structure
                rmsd_val = backbone_rmsd(cand["rfd3_array"], rf3_arr)

                # Interface contacts from RF3 predicted complex
                n_contacts = count_interface_contacts(rf3_arr, bc, fc, cutoff=5.0)

                # Persist structure
                cif_path = out_dir / f"{did}_rf3.cif"
                to_cif_file(rf3_arr, str(cif_path), file_type="cif")

                metrics_rec = {
                    "plddt":          plddt,
                    "ptm":            ptm,
                    "iptm":           iptm,
                    "interface_pae":  iface_pae,
                    "rmsd_to_rfd3":   rmsd_val,
                    "n_contacts":     n_contacts,
                    "source":         "RF3",
                }
                # Source-aware composite (RF3-only branch uses geometric features only)
                metrics_rec["composite_score"] = calc_composite(metrics_rec)
                with open(out_dir / f"{did}_metrics.json", "w") as mf:
                    json.dump(metrics_rec, mf, indent=2)

                scored.append({
                    **{k: v for k, v in cand.items() if k not in ("array", "rfd3_array")},
                    **metrics_rec,
                    "iteration":    self.state["iteration"],
                    "temperature":  self.state["temperature"],
                    "rf3_cif":      str(cif_path),
                    "promising":    False,  # set by run_boltz once Boltz ipTM is in hand
                })

            except Exception as exc:
                logger.error(f"RF3 error on {did}: {exc}")

        logger.info(f"[{target_name}] RF3 done in {(time.time()-t0)/60:.1f} min "
                    f"({len(scored)}/{len(candidates)} scored)")
        del engine
        torch.cuda.empty_cache()
        return scored

    # ── Hall of Fame ──────────────────────────────────────────────────────────

    @staticmethod
    def _consolidate_hof(entries: list, size: int) -> list:
        """
        Dedup a HOF by design_id, protect AF3-validated entries, trim to `size`.

        - Dedup (idempotent across resume): per design_id keep the AF3 version if
          present, else the higher composite. Stops an interrupted-then-resumed
          iteration from double-counting designs.
        - AF3 entries are never trimmed (scarce ground truth must not be displaced
          by an un-validated local score); other slots fill by composite.

        Used for both the per-target HOF and the specificity HOF.
        """
        best_by_id: dict = {}
        for e in entries:
            did = e.get("design_id")
            cur = best_by_id.get(did)
            if cur is None:
                best_by_id[did] = e
                continue
            e_af3, cur_af3 = e.get("source") == "AF3", cur.get("source") == "AF3"
            if (e_af3 and not cur_af3) or (
                e_af3 == cur_af3
                and e.get("composite_score", 0) > cur.get("composite_score", 0)
            ):
                best_by_id[did] = e
        uniq = list(best_by_id.values())

        af3   = [e for e in uniq if e.get("source") == "AF3"]
        other = [e for e in uniq if e.get("source") != "AF3"]
        other.sort(key=lambda x: x.get("composite_score", 0), reverse=True)
        merged = af3 + other[: max(0, size - len(af3))]
        return sorted(merged, key=lambda x: x.get("composite_score", 0), reverse=True)

    def update_hof(self, new_scored: list) -> None:
        """
        Add new scored designs to per-target HOF; keep top HOF_SIZE_PER_TARGET.

        AF3-validated entries are NEVER trimmed: AF3 results are the scarce
        ground truth (30/day) and the local composite is only a proxy, so an
        un-validated ESMFold2 entry must not displace a validated one. We keep
        all AF3 entries, then fill the remaining slots with the best others.
        """
        for entry in new_scored:
            tname = entry.get("target_name")
            if tname not in self.state["hof"]:
                self.state["hof"][tname] = []
            self.state["hof"][tname].append(entry)

        for tname in self.active_targets:
            self.state["hof"][tname] = self._consolidate_hof(
                self.state["hof"].get(tname, []), HOF_SIZE_PER_TARGET
            )

        # Save best HOF structures for later seeding/inspection
        self._save_hof_structures()
        self._write_hof_summary()

    def _save_hof_structures(self) -> None:
        """
        Copy the top-3 predicted complexes per target into hof_structures/ for
        quick inspection. Prefers the ESMFold2 structure (our ranker's output),
        falls back to the RF3 CIF if RF3 was enabled. No-op if neither exists.
        """
        import shutil
        for tname in self.active_targets:
            tdir = self.hof_struct_dir / tname
            tdir.mkdir(parents=True, exist_ok=True)
            for rank, entry in enumerate(self.state["hof"].get(tname, [])[:3], start=1):
                src = entry.get("esm_cif") or entry.get("rf3_cif")
                if src and Path(src).exists():
                    dest = tdir / f"hof_rank{rank:02d}_{entry['design_id']}{Path(src).suffix}"
                    if not dest.exists():
                        shutil.copy2(src, dest)

    def _write_hof_summary(self) -> None:
        """Write a combined CSV of all HOF entries across targets."""
        rows = []
        for tname in self.active_targets:
            for rank, entry in enumerate(self.state["hof"].get(tname, []), start=1):
                rows.append({"hof_rank": rank, **{k: v for k, v in entry.items()
                                                   if k not in ("array", "rfd3_array")}})
        if rows:
            pd.DataFrame(rows).to_csv(OUT_BASE / "hof_summary.csv", index=False)

    # ── Per-iteration summary ─────────────────────────────────────────────────

    def _write_round_summary(self, all_scored: list, it_dir: Path) -> None:
        """Write CSV of all designs scored this iteration."""
        rows = [{k: v for k, v in e.items() if k not in ("array", "rfd3_array")}
                for e in all_scored]
        if rows:
            df = pd.DataFrame(rows)
            df.to_csv(it_dir / "round_summary.csv", index=False)
            n_promising = int(df["promising"].sum()) if "promising" in df.columns else 0
            best = df.loc[df["composite_score"].idxmax()] if "composite_score" in df.columns else None
            logger.info(f"Round summary: {len(df)} designs, {n_promising} promising.")
            if best is not None:
                logger.info(
                    f"Best this round: {best.get('design_id')} | "
                    f"target={best.get('target_name')} | "
                    f"ipTM={best.get('iptm', 0):.3f} | "
                    f"pLDDT={best.get('plddt', 0):.1f} | "
                    f"composite={best.get('composite_score', 0):.3f}"
                )

    # ── AF3 export / import ───────────────────────────────────────────────────

    def export_for_af3(self, force: bool = False) -> None:
        """
        Export AF3_TOP_N designs to an AF3 Server JSON file, allocated as an
        equal per-target quota so every target gets AF3 validation each cycle.
        Within each target we take the highest-composite unique sequences.
        """
        it = self.state["iteration"]
        if not force and (it - self.state.get("last_af3_it", -1)) < AF3_EXPORT_EVERY_N:
            return

        n_targets = max(1, len(self.active_targets))
        per_target_quota = max(1, AF3_TOP_N // n_targets)

        seen_seqs = set()
        unique_entries = []
        for tname in self.active_targets:
            target_hof = sorted(
                self.state["hof"].get(tname, []),
                key=lambda x: x.get("composite_score", 0),
                reverse=True,
            )
            picked = 0
            for e in target_hof:
                seq = e.get("full_seq", "")
                if seq and seq not in seen_seqs:
                    seen_seqs.add(seq)
                    unique_entries.append(e)
                    picked += 1
                if picked >= per_target_quota:
                    break
            logger.info(f"AF3 export: {tname} contributing {picked} design(s)")

        if not unique_entries:
            logger.warning("HOF empty; skipping AF3 export.")
            return

        # Build AF3 Server JSON
        target_seqs_by_entry = {
            e["design_id"]: self.target_seqs.get(e.get("target_name", ""), "")
            for e in unique_entries
        }
        jobs = []
        for i, e in enumerate(unique_entries):
            tseq = target_seqs_by_entry.get(e["design_id"], e.get("target_seq", ""))
            jobs.append({
                "name":        f"refine_it{it}_{e['target_name']}_{i:02d}",
                "modelSeeds":  [42],
                "sequences": [
                    {"proteinChain": {"sequence": e.get("full_seq", ""),  "count": 1}},
                    {"proteinChain": {"sequence": tseq,                   "count": 1}},
                ],
            })

        export_path = OUT_BASE / f"af3_submission_it{it}.json"
        with open(export_path, "w") as f:
            json.dump(jobs, f, indent=2)

        self.state["last_af3_it"] = it
        self._save_state()

        logger.info(f"AF3 export: {len(jobs)} designs → {export_path}")
        print("\n" + "=" * 60)
        print(f"  AlphaFold3 submission ready: {export_path}")
        print(f"  Upload this JSON to the AF3 Server, then place the results")
        print(f"  (parsed into af3_results_it{it}.json) in {OUT_BASE}/")
        print(f"  Expected format: list of dicts with keys: design_id, iptm,")
        print(f"  plddt, ptm, interface_pae, full_seq, target_name")
        print("=" * 60 + "\n")

    def export_for_af3_stratified(self, n_total: int = 30, n_bands: int = 3,
                                  strat_metric: str = None) -> None:
        """
        One-time VALIDATION / CALIBRATION-TRANCHE export: instead of sending the
        top-scoring designs (which all sit in a narrow high band and cause range
        restriction), sample evenly across LOCAL-SCORE bands so the returned AF3
        scores span predicted-strong → predicted-weak. This is the deliberate way
        to break the dynamic-range wall (docs/09): the next wet-lab calibration
        finally gets variance to fit.

        Stratifies on `strat_metric` (default: composite_score, the live ESMFold2
        ranker; falls back to boltz_iptm if composite_score is absent). Pools ALL
        designs from every it_*/round_summary.csv (the unbiased pool, not the
        range-restricted HOF), bins each target into `n_bands` equal-frequency
        bands (LO/MID/HI) and samples evenly. Writes the AF3 submission JSON plus a
        sidecar manifest (stratified_manifest.json) recording each job's band and
        score, so the follow-up analysis can compute per-band AF3 metrics without
        fragile job-name parsing.
        """
        # Gather the full design pool from round summaries (unique by sequence)
        pool_frames = []
        for it_dir in sorted(OUT_BASE.glob("it_*")):
            csv = it_dir / "round_summary.csv"
            if csv.exists():
                pool_frames.append(pd.read_csv(csv))
        if not pool_frames:
            logger.warning("No round_summary.csv found; cannot stratify.")
            return

        pool = pd.concat(pool_frames, ignore_index=True)
        # Pick the stratification metric: explicit > composite_score > boltz_iptm.
        sm = strat_metric
        if sm is None:
            sm = "composite_score" if pool.get("composite_score") is not None \
                 and pd.to_numeric(pool["composite_score"], errors="coerce").notna().any() else "boltz_iptm"
        if pool.get(sm) is None:
            logger.warning(f"Stratify metric '{sm}' not in round_summary; cannot stratify.")
            return
        pool[sm] = pd.to_numeric(pool[sm], errors="coerce")
        pool = pool.dropna(subset=[sm]).drop_duplicates(subset=["full_seq"])
        # Drop UNSCORED designs: an ESMFold2 failure (e.g. OOM) leaves esm_iptm NaN
        # and composite_score == 0. Including them would fill the LO band with junk
        # rather than genuinely predicted-weak designs. Require a real ESMFold2 score.
        if "esm_iptm" in pool.columns:
            pool = pool[pd.to_numeric(pool["esm_iptm"], errors="coerce").notna()]
        if sm == "composite_score":
            pool = pool[pool[sm] > 0]
        logger.info(f"Stratifying AF3 tranche on '{sm}' over {len(pool)} SCORED pooled designs.")

        n_targets       = max(1, len(self.active_targets))
        per_target      = max(n_bands, n_total // n_targets)
        per_band        = max(1, per_target // n_bands)

        jobs, manifest = [], []
        for tname in self.active_targets:
            tpool = pool[pool["target_name"] == tname].sort_values(sm)
            if len(tpool) < n_bands:
                logger.warning(f"[{tname}] too few designs ({len(tpool)}) to stratify; skipping.")
                continue

            # Equal-frequency bands (low → high local score)
            band_labels = ["LO", "MID", "HI"][:n_bands]
            bands = np.array_split(tpool, n_bands)   # ascending, so [0]=low ... [-1]=high
            for label, band_df in zip(band_labels, bands):
                # Even spread within the band rather than all from one end
                take = band_df.iloc[np.linspace(0, len(band_df) - 1, min(per_band, len(band_df))).astype(int)]
                for _, e in take.iterrows():
                    idx  = len(jobs)
                    tseq = self.target_seqs.get(tname, "") or e.get("target_seq", "")
                    name = f"strat_{tname}_{label}_{idx:02d}"
                    jobs.append({
                        "name":       name,
                        "modelSeeds": [42],
                        "sequences": [
                            {"proteinChain": {"sequence": e.get("full_seq", ""), "count": 1}},
                            {"proteinChain": {"sequence": tseq,                  "count": 1}},
                        ],
                    })
                    manifest.append({
                        "job_name":    name,
                        "design_id":   e.get("design_id"),
                        "target":      tname,
                        "band":        label,
                        "strat_metric": sm,
                        "strat_score": float(e[sm]),
                        "full_seq":    e.get("full_seq", ""),
                    })
            logger.info(f"[{tname}] stratified: "
                        + ", ".join(f"{lbl}={sum(1 for mm in manifest if mm['target']==tname and mm['band']==lbl)}"
                                    for lbl in band_labels))

        if not jobs:
            logger.warning("Stratified export produced no jobs.")
            return

        sub_path = OUT_BASE / "af3_submission_stratified.json"
        man_path = OUT_BASE / "stratified_manifest.json"
        with open(sub_path, "w") as f:
            json.dump(jobs, f, indent=2)
        with open(man_path, "w") as f:
            json.dump(manifest, f, indent=2)

        logger.info(f"Stratified AF3 export: {len(jobs)} designs → {sub_path}")
        print("\n" + "=" * 64)
        print(f"  STRATIFIED validation submission ready: {sub_path}")
        print(f"  Manifest (band assignments): {man_path}")
        print(f"  Upload to AF3, then join the results back to the manifest by band")
        print(f"  to see AF3 hit-rate per band (low/mid/high local score).")
        print(f"  See filter_methods.md for the calibration method.")
        print("=" * 64 + "\n")

    def import_af3_results(self, results_path: str) -> None:
        """
        Import AF3 Server results to update HOF with higher-quality metrics.
        Expected JSON: list of dicts with design_id, iptm, plddt, ptm,
                       interface_pae, full_seq, target_name.
        """
        with open(results_path) as f:
            results = json.load(f)

        logger.info(f"Importing {len(results)} AF3 results from {results_path}")
        for res in results:
            res["source"] = "AF3"
            loops = extract_loops(res.get("full_seq", ""), self.selected_loops)
            res.update(loops)
            tname = res.get("target_name", "unknown")
            if tname in self.state["hof"]:
                # Replace any existing entry with same design_id, else append
                existing_ids = {e["design_id"] for e in self.state["hof"][tname]}
                if res.get("design_id") in existing_ids:
                    self.state["hof"][tname] = [
                        res if e["design_id"] == res["design_id"] else e
                        for e in self.state["hof"][tname]
                    ]
                else:
                    self.state["hof"][tname].append(res)

        self.update_hof([])  # Re-sort and trim
        self._save_state()
        logger.info("AF3 import complete.")

    def import_af3_zip(self, zip_path: str) -> None:
        """
        Parse a batch AF3 Server ZIP (one subdirectory per job) and update the HOF.

        Expected layout inside the ZIP:
          <af3_job_name>/
            <af3_job_name>_job_request.json          original name + sequences
            <af3_job_name>_summary_confidences_0.json  top-ranked model metrics
            <af3_job_name>_full_data_0.json            PAE matrix, per-atom pLDDT

        AF3 lowercases and optionally prepends "fold_" to the submitted job name.
        The original name (e.g. "refine_it4_MMP9_01") is recovered from job_request.json.
        Designs are matched back to HOF entries by binder sequence identity.
        """
        import zipfile

        # One-time index of binder sequences across all HOFs so we can resolve
        # (target_name, design_id) in O(1) per job instead of scanning each HOF.
        seq_index: dict[str, tuple[str, str]] = {}
        for t in self.active_targets:
            for e in self.state["hof"].get(t, []):
                seq = e.get("full_seq")
                if seq and seq not in seq_index:
                    seq_index[seq] = (t, e["design_id"])

        with zipfile.ZipFile(zip_path) as zf:
            all_names = set(zf.namelist())
            job_requests = sorted(n for n in all_names if n.endswith("_job_request.json"))

            if not job_requests:
                logger.error(f"No *_job_request.json found in {zip_path}")
                return

            logger.info(f"Found {len(job_requests)} AF3 jobs in {zip_path}")
            results = []

            for jrf in job_requests:
                try:
                    raw = json.loads(zf.read(jrf))
                    job_data = raw[0] if isinstance(raw, list) else raw

                    job_name = job_data.get("name", "")
                    seqs = job_data.get("sequences", [])
                    if len(seqs) < 2:
                        logger.warning(f"Skipping {jrf}: fewer than 2 sequences")
                        continue
                    binder_seq = seqs[0]["proteinChain"]["sequence"]

                    # Recover target name from the job name pattern refine_it{N}_{TARGET}_{NN}
                    m = re.search(r"(?:fold_)?refine_it\d+_([A-Za-z0-9]+)_\d{2}$", job_name, re.I)
                    target_name = m.group(1).upper() if m else None

                    # Fallback: sequence-match against the precomputed seq_index
                    if target_name not in self.active_targets:
                        target_name = seq_index.get(binder_seq, (None, None))[0]

                    if not target_name:
                        logger.warning(f"Cannot determine target for job '{job_name}'; skipping.")
                        continue

                    # Sibling files share the same path prefix
                    prefix = jrf[: -len("_job_request.json")]

                    sc_name = prefix + "_summary_confidences_0.json"
                    if sc_name not in all_names:
                        logger.warning(f"Missing summary_confidences for {job_name}")
                        continue
                    sc    = json.loads(zf.read(sc_name))
                    iptm  = float(sc.get("iptm", 0.0))
                    ptm   = float(sc.get("ptm",  0.0))
                    has_clash = float(sc.get("has_clash", 0.0))
                    # ApTM = binder (chain A) pTM — a calibrated positive binding
                    # prior. chain_ptm = [ApTM, BpTM]; we take [0] and DROP BpTM.
                    cp    = sc.get("chain_ptm") or []
                    aptm  = float(cp[0]) if len(cp) > 0 else float("nan")

                    # Interface geometry from the AF3 complex CIF (same extractor
                    # as calibration): interface size + interface-residue pLDDT are
                    # calibrated positive priors.
                    af3_iface = {}
                    cif_name = prefix + "_model_0.cif"
                    if cif_name in all_names:
                        try:
                            af3_iface = cs.interface_features_from_cif(zf.read(cif_name).decode())
                        except Exception as exc:
                            logger.debug(f"AF3 interface feats failed for {job_name}: {exc}")

                    plddt     = 0.0
                    iface_pae = float("nan")

                    fd_name = prefix + "_full_data_0.json"
                    if fd_name in all_names:
                        fd = json.loads(zf.read(fd_name))

                        # Binder pLDDT: chain A = first submitted chain = binder
                        atom_plddts = np.array(fd.get("atom_plddts", []))
                        atom_chains = fd.get("atom_chain_ids", [])
                        if atom_plddts.size and atom_chains:
                            a_mask = np.array([c == "A" for c in atom_chains])
                            if a_mask.any():
                                plddt = float(atom_plddts[a_mask].mean())

                        # Interface PAE: average of both cross-chain PAE blocks
                        raw_pae      = fd.get("pae")
                        token_chains = fd.get("token_chain_ids", [])
                        if raw_pae and token_chains:
                            pae_mat = np.array(raw_pae)
                            a_idx   = np.array([i for i, c in enumerate(token_chains) if c == "A"])
                            b_idx   = np.array([i for i, c in enumerate(token_chains) if c == "B"])
                            if a_idx.size and b_idx.size:
                                iface_pae = float(
                                    (pae_mat[np.ix_(a_idx, b_idx)].mean()
                                     + pae_mat[np.ix_(b_idx, a_idx)].mean()) / 2
                                )

                    # Resolve design_id via the precomputed seq_index (O(1));
                    # fall back to the AF3 job name if no HOF match found.
                    design_id = seq_index.get(binder_seq, (None, job_name))[1]

                    loops = extract_loops(binder_seq, self.selected_loops)
                    entry = {
                        "design_id":     design_id,
                        "target_name":   target_name,
                        "full_seq":      binder_seq,
                        "iptm":          iptm,
                        "ptm":           ptm,
                        "plddt":         plddt,
                        "interface_pae": iface_pae,
                        "rmsd_to_rfd3":  float("nan"),
                        "has_clash":     has_clash,
                        # calibrated positive priors (BpTM intentionally omitted)
                        "af3_aptm":                  aptm,
                        "af3_iface_n_iface_res":     af3_iface.get("n_iface_res"),
                        "af3_iface_iface_plddt":     af3_iface.get("iface_plddt"),
                        "af3_iface_contact_density": af3_iface.get("contact_density"),
                        "source":        "AF3",
                        **loops,
                    }
                    entry["composite_score"] = calc_composite(entry)
                    results.append(entry)
                    logger.info(
                        f"  {job_name} | {target_name} | "
                        f"ipTM={iptm:.3f}  pLDDT={plddt:.1f}  iface_PAE={iface_pae:.2f}"
                        + ("  [CLASH]" if has_clash > 0.5 else "")
                    )

                except Exception as exc:
                    logger.error(f"Error parsing {jrf}: {exc}", exc_info=True)

        logger.info(f"Parsed {len(results)} AF3 results from ZIP.")

        # Merge into HOF: replace existing entry if design_id matches, else append
        for res in results:
            tname = res["target_name"]
            if tname not in self.state["hof"]:
                self.state["hof"][tname] = []
            existing_ids = {e["design_id"] for e in self.state["hof"][tname]}
            if res["design_id"] in existing_ids:
                self.state["hof"][tname] = [
                    res if e["design_id"] == res["design_id"] else e
                    for e in self.state["hof"][tname]
                ]
            else:
                self.state["hof"][tname].append(res)

        self.update_hof([])
        self._save_state()
        logger.info("AF3 ZIP import complete.")

    # ── Core iteration ────────────────────────────────────────────────────────

    def run_iteration(self) -> None:
        """Execute one full RFd3 → LigandMPNN → RF3 round for all active targets."""
        it   = self.state["iteration"]
        temp = self.state["temperature"]
        it_dir = OUT_BASE / f"it_{it}"
        it_dir.mkdir(parents=True, exist_ok=True)

        logger.info(
            f"\n{'='*60}\n  Iteration {it} | T={temp:.3f} | "
            f"Targets: {', '.join(self.active_targets)}\n{'='*60}"
        )

        adaptive_ranges = self._adaptive_loop_ranges()
        all_scored: list = []

        for tname in self.active_targets:
            # Per-target isolation: a failure on one target (OOM, CUDA hiccup,
            # MSA-server blip) must not abort the whole unattended run. Log it,
            # keep whatever scored, and move on.
            try:
                tcfg     = TARGETS[tname]
                pdb_path = str(DATA_DIR / tcfg["pdb"])
                if not Path(pdb_path).exists():
                    logger.warning(f"Skipping {tname}: PDB not found at {pdb_path}")
                    continue

                # 1. Backbone generation
                rfd3_dir = it_dir / "rfd3" / tname
                rfd3_dir.mkdir(parents=True, exist_ok=True)
                backbones = self.run_rfd3(tname, pdb_path, BACKBONES_PER_TARGET, adaptive_ranges)

                # Save backbone CIFs (the structural record of each design)
                for bi, arr in enumerate(backbones):
                    to_cif_file(arr, str(rfd3_dir / f"bb_{bi}.cif"), file_type="cif")

                if not backbones:
                    logger.warning(f"[{tname}] No backbones generated; skipping.")
                    continue

                # 2. Sequence design
                lmpnn_dir  = it_dir / "lmpnn" / tname
                candidates = self.run_lmpnn(tname, backbones, lmpnn_dir, temp)

                if not candidates:
                    logger.warning(f"[{tname}] No LMPNN sequences; skipping.")
                    continue

                # 3. RF3 (OFF by default; finalizes candidates + optional geometric features)
                scored = self.run_rf3_complex(tname, candidates, it_dir / "rf3" / tname)

                # 4. ESMFold2 ranker (default). Boltz retained but off by default.
                # Each is a no-op if disabled or its env is missing.
                scored = self.run_esmfold2(tname, scored, it_dir / "esmfold2" / tname)
                scored = self.run_boltz(tname, scored, it_dir / "boltz" / tname)
                all_scored.extend(scored)
            except Exception as exc:
                logger.error(f"[{tname}] iteration {it} failed: {exc}", exc_info=True)
                torch.cuda.empty_cache()
                continue

        # Update HOF and write summaries (even if some targets failed)
        self.update_hof(all_scored)
        self._write_round_summary(all_scored, it_dir)
        self._save_state()

    # ── Main loop ─────────────────────────────────────────────────────────────

    def main_loop(self, max_iterations: int = None) -> None:
        """
        Run the iterative design loop indefinitely (or up to max_iterations).
        Exports AF3 submission every AF3_EXPORT_EVERY_N iterations.
        """
        logger.info(
            f"Starting iterative refinement.\n"
            f"  Targets: {self.active_targets}\n"
            f"  Loops:   {self.active_loops}\n"
            f"  Output:  {OUT_BASE}\n"
        )

        it = 0
        while max_iterations is None or it < max_iterations:
            self.run_iteration()

            # AF3 export check
            self.export_for_af3()

            # Temperature annealing
            self.state["temperature"] = max(
                MIN_TEMPERATURE, self.state["temperature"] * TEMP_DECAY
            )
            self.state["iteration"] += 1
            self._save_state()

            it += 1
            logger.info(
                f"Iteration complete. Next temperature: {self.state['temperature']:.3f}\n"
            )

        logger.info("Reached max_iterations. Exiting.")
        self.export_for_af3(force=True)


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    import argparse
    # CLI flags can override these module-level settings — declare up front so the
    # argparse help strings (which read the defaults) don't trip Python's
    # "name used prior to global declaration" rule.
    global RF3_ENABLE, ESMFOLD2_GPUS, INIT_TEMPERATURE, MIN_TEMPERATURE, TEMP_DECAY, ADAPTIVE_BIAS_START
    global SV_BATTERY, SV_OCCLUSION_FILTER, SV_OCCLUSION_MIN
    global BACKBONES_PER_TARGET, LMPNN_SEQS_PER_BACKBONE

    parser = argparse.ArgumentParser(description="Iterative TIMP3 binder design.")
    parser.add_argument(
        "--targets", nargs="+", default=["MMP2", "MMP9", "ADAM10", "ADAM17"],
        choices=list(TARGETS.keys()),  # also available: MMP3, MMP10
        help="Protease targets to include.",
    )
    parser.add_argument(
        "--loops", nargs="+", default=["AB", "C", "EF"],
        choices=list(LOOP_CONFIGS.keys()),
        help="Loop regions to redesign.",
    )
    parser.add_argument(
        "--max-iterations", type=int, default=None,
        help="Stop after this many iterations (default: run indefinitely).",
    )
    parser.add_argument(
        "--stratified-export", type=int, default=None, metavar="N",
        help=(
            "Write a one-time STRATIFIED AF3 submission of N designs sampled "
            "evenly across local-score bands (to calibrate the filter against AF3), "
            "then exit without running iterations.  Use ~30 to match the daily "
            "AF3 cap.  Analyze the returned results per band (see filter_methods.md)."
        ),
    )
    parser.add_argument(
        "--import-af3", type=str, default=None,
        help=(
            "Path to AF3 results to import before continuing.  "
            "Accepts either a batch ZIP downloaded directly from the AF3 Server "
            "(one subdirectory per job) or a hand-crafted JSON list of dicts with keys: "
            "design_id, iptm, plddt, ptm, interface_pae, full_seq, target_name."
        ),
    )
    parser.add_argument(
        "--enable-rf3", action="store_true",
        help="Run RF3 too (OFF by default). Adds its geometric features to "
             "round_summary.csv; does not change ranking (ESMFold2 ranks). Slower.",
    )
    parser.add_argument(
        "--sv-battery", action="store_true",
        help="Log the full Structural-Validation interface battery (BSA, H-bonds, "
             "shape complementarity, catalytic occlusion, pDockQ, DockQ-vs-WT, ...) "
             "as sv_* columns in round_summary.csv. LOG-ONLY: does NOT change "
             "ranking (the FCS calibration showed these don't predict binding). "
             "Slower (SASA/shape per candidate).",
    )
    parser.add_argument(
        "--sv-occlusion-filter", action="store_true",
        help="With --sv-battery: drop designs whose reactive edge does not reach "
             "the catalytic zinc cleft — a mechanistic sanity gate like the "
             "foldability floor, NOT a metric in the composite.",
    )
    parser.add_argument(
        "--sv-occlusion-min", type=float, default=None, metavar="FRAC",
        help=f"Min buried fraction of the zinc loop for the occlusion gate "
             f"(default {SV_OCCLUSION_MIN}).",
    )
    parser.add_argument(
        "--esmfold2-gpus", default=None, metavar="N|auto",
        help="Data-parallel the ESMFold2 stage across GPUs (default 1). An integer "
             "uses up to N *free* GPUs; 'auto' uses all free GPUs (cap "
             f"{ESMFOLD2_MAX_GPUS}). Skips GPUs busy with other jobs.",
    )
    # ── Diversity / exploration controls ──
    parser.add_argument("--init-temperature", type=float, default=None,
                        help=f"Starting LMPNN temperature for a fresh run (default {INIT_TEMPERATURE}).")
    parser.add_argument("--min-temperature", type=float, default=None,
                        help=f"Temperature floor (default {MIN_TEMPERATURE}). Raise to ~0.30 to keep "
                             "exploring (sustained sequence diversity) instead of converging.")
    parser.add_argument("--temp-decay", type=float, default=None,
                        help=f"Per-iteration temperature multiplier (default {TEMP_DECAY}). Raise "
                             "toward ~0.95 to stay explorative for more iterations.")
    parser.add_argument("--no-adaptive-bias", action="store_true",
                        help="Disable adaptive loop-length narrowing — keeps the full loop-length "
                             "range every iteration (more structural diversity).")
    # ── Throughput controls ──
    parser.add_argument("--backbones-per-target", type=int, default=None,
                        help=f"RFd3 backbones per target per iteration (default {BACKBONES_PER_TARGET}). "
                             "The slow stage — scales RFd3 time.")
    parser.add_argument("--seqs-per-backbone", type=int, default=None,
                        help=f"LMPNN sequences per backbone (default {LMPNN_SEQS_PER_BACKBONE}). Cheap "
                             "way to add sequence diversity + ESMFold2 predictions without more RFd3.")
    args = parser.parse_args()

    if args.enable_rf3:
        RF3_ENABLE = True
    if args.sv_battery:
        SV_BATTERY = True
    if args.sv_occlusion_filter:
        SV_OCCLUSION_FILTER = True
    if args.sv_occlusion_min is not None:
        SV_OCCLUSION_MIN = args.sv_occlusion_min
    if args.esmfold2_gpus is not None:
        ESMFOLD2_GPUS = args.esmfold2_gpus  # int-like string or "auto"; resolved at call time
    if args.init_temperature is not None:
        INIT_TEMPERATURE = args.init_temperature
    if args.min_temperature is not None:
        MIN_TEMPERATURE = args.min_temperature
    if args.temp_decay is not None:
        TEMP_DECAY = args.temp_decay
    if args.no_adaptive_bias:
        ADAPTIVE_BIAS_START = 10**9   # effectively never
    if args.backbones_per_target is not None:
        BACKBONES_PER_TARGET = args.backbones_per_target
    if args.seqs_per_backbone is not None:
        LMPNN_SEQS_PER_BACKBONE = args.seqs_per_backbone

    refiner = IterativeRefiner(
        active_targets=args.targets,
        active_loops=args.loops,
    )

    if args.import_af3:
        if args.import_af3.lower().endswith(".zip"):
            refiner.import_af3_zip(args.import_af3)
        else:
            refiner.import_af3_results(args.import_af3)

    if args.stratified_export:
        # Validation mode: write the stratified batch and stop (no iterations).
        refiner.export_for_af3_stratified(n_total=args.stratified_export)
        return

    refiner.main_loop(max_iterations=args.max_iterations)


if __name__ == "__main__":
    main()
