"""
iterative_refinement.py

Iterative binder design with a four-stage in-silico funnel before any AF3 call:

  RFd3 (backbone) → LigandMPNN (sequence) → RF3 with target template (sanity) →
  Boltz-2 (local AF3-class ranker) → AF3 Server (gold-standard validation, capped at 30/day)

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
DATA_DIR  = _HERE / ".." / "Data" / "TIMP_Complexes" / "HADDOCK_Outputs"
OUT_BASE  = _HERE / ".." / "Local" / "iterative_refinement"
CKPT_DIR  = _HERE / ".." / "Tools" / "foundry_checkpoints"

# ── Target definitions ────────────────────────────────────────────────────────
# binder_chain / target_chain refer to chain IDs in the SOURCE PDB only.
# They are used to read target_seqs and to construct the RFd3 contig string.
# After RFd3 runs, all output structures (RFd3, LMPNN, RF3) follow a fixed
# convention regardless of input — see DESIGN_BINDER_CHAIN below.
TARGETS = {
    "MMP2":   {"pdb": "MMP2_TIMP3_HADDOCK.pdb",   "binder_chain": "B", "target_chain": "A", "scaffold_len": 121},
    "MMP9":   {"pdb": "MMP9_TIMP3_HADDOCK.pdb",   "binder_chain": "B", "target_chain": "A", "scaffold_len": 121},
    "MMP3":   {"pdb": "MMP3_TIMP3_HADDOCK.pdb",   "binder_chain": "B", "target_chain": "A", "scaffold_len": 121},
    "MMP10":  {"pdb": "MMP10_TIMP3_HADDOCK.pdb",  "binder_chain": "B", "target_chain": "A", "scaffold_len": 121},
    "ADAM10": {"pdb": "ADAM10_TIMP3_HADDOCK.pdb",  "binder_chain": "B", "target_chain": "A", "scaffold_len": 121},
    "ADAM17": {"pdb": "ADAM17_TIMP3_HADDOCK.pdb",  "binder_chain": "B", "target_chain": "A", "scaffold_len": 121},
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
    "AB": {"normal": 6,  "max": 15, "pos": 30,  "left": "LVK", "right": "LVY"},
    "C":  {"normal": 6,  "max": 15, "pos": 62,  "left": "HTE", "right": "GLK"},
    "EF": {"normal": 4,  "max": 10, "pos": 92,  "left": "MYT", "right": "FVE"},
    "GH": {"normal": 10, "max": 20, "pos": 127, "left": "KSC", "right": "NEC"},
}

# ── Hyperparameters ───────────────────────────────────────────────────────────
BACKBONES_PER_TARGET   = 20    # RFd3 designs per target per iteration
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
#
# Why split: in our calibration against AF3 (n=18), RF3 single-sequence ipTM/pTM
# were anti-correlated with AF3 ipTM (r = -0.07 / -0.51).  Including those metrics
# in a unified ranker introduces noise.  Boltz-2 and AF3, which use MSA/templates,
# should track each other well — we trust them more.
#
# Each formula returns a value in [0, 1].  The fallback RF3 composite is capped
# at RF3_COMPOSITE_CEILING so any AF3- or Boltz-scored entry naturally ranks
# above an RF3-only entry of similar geometric quality.
COMPOSITE_AF3 = {        # AF3-validated entries: full trust
    "iptm":  0.45,
    "plddt": 0.25,
    "ptm":   0.15,
    "rmsd":  0.05,
    "pae":   0.10,
}
COMPOSITE_BOLTZ = {      # Boltz-2 scored: AF3-class metrics, no affinity head
    "iptm":   0.45,      # Boltz-2's affinity head only runs for protein-ligand
    "plddt":  0.25,      # complexes (verified Boltz 2.2.1 — rejects protein-protein
    "ptm":    0.15,      # in `properties.affinity`), so we rely on the structure-
    "pae":    0.15,      # confidence metrics, which mirror AF3-class scoring.
}
COMPOSITE_RF3_GEOM = {   # RF3-only fallback: geometric features only (no ipTM/pTM/pLDDT)
    "n_contacts": 0.40,  # interface size
    "pae":        0.30,  # weak +signal in our calibration
    "rmsd":       0.30,  # backbone fidelity to RFd3 design
}
RF3_COMPOSITE_CEILING = 0.50  # cap RF3-only composite so AF3/Boltz dominate HOF

# ── Boltz-2 configuration ────────────────────────────────────────────────────
# Boltz-2 is an open-source AF3-class model used here as the pre-AF3 local ranker.
# Disk footprint ≈ 3 GB weights + ~2 GB ColabFold MSA cache (per-user, in ~/.boltz/).
#
# DEPLOYMENT NOTE: Boltz pins numpy<2.0 / cublas<12.5 / older lightning.  When
# installed into a foundry env that previously had newer versions, those get
# downgraded — fine on V100 (DISABLE_CUEQUIVARIANCE is auto-set there), but
# breaks RF3 on A100/H100 which requires cuEquivariance.  If you move this job
# to A100/H100 hardware, install Boltz in a separate conda env and set
# BOLTZ_EXECUTABLE to point at it:
#     export BOLTZ_EXECUTABLE=/x/capa/<user>/miniconda3/envs/boltz/bin/boltz
BOLTZ_ENABLE            = True
BOLTZ_EXECUTABLE        = os.environ.get("BOLTZ_EXECUTABLE", "boltz")
BOLTZ_USE_MSA_SERVER    = True   # use ColabFold's remote MMseqs2 API (no local DB)
BOLTZ_DIFFUSION_SAMPLES = 1      # samples per design (1 is fastest, good enough for ranking)
BOLTZ_TIMEOUT_S         = 300    # 5 min hard timeout (real runs are ~30s post-cache-warm)


# ── Utility functions ─────────────────────────────────────────────────────────

def setup_env() -> None:
    OUT_BASE.mkdir(parents=True, exist_ok=True)
    CKPT_DIR.mkdir(parents=True, exist_ok=True)
    os.environ["FOUNDRY_CHECKPOINT_DIRS"] = str(CKPT_DIR)
    os.environ["DGLBACKEND"] = "pytorch"


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

        AF3-validated  →  COMPOSITE_AF3      (full trust)
        Boltz-scored   →  COMPOSITE_BOLTZ    (structure-confidence metrics)
        RF3-only       →  COMPOSITE_RF3_GEOM (geometric features only, capped)

    The RF3-only branch deliberately excludes RF3's iptm/pTM/pLDDT because they
    are anti-correlated with AF3 ipTM for de novo binders (calibrated empirically).
    It is also capped at RF3_COMPOSITE_CEILING so any AF3 or Boltz score dominates.
    """
    src = entry.get("source", "RF3")

    # ── AF3-validated ──
    if src == "AF3":
        w = COMPOSITE_AF3
        return (
            w["iptm"]  * _safe(entry.get("iptm"))
            + w["plddt"] * (_normalize_plddt(entry.get("plddt")) / 100.0)
            + w["ptm"]   * _safe(entry.get("ptm"))
            + w["rmsd"]  * _rmsd_score(entry.get("rmsd_to_rfd3"))
            + w["pae"]   * _pae_score(entry.get("interface_pae"))
        )

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


def _parse_boltz_output(run_dir: Path, binder_chain: str, target_chain: str,
                        binder_len: int) -> dict:
    """
    Pull metrics from a Boltz-2 prediction directory.

    Boltz 2.2.x layout (verified on caprine, 2026-05):
      <run_dir>/
        boltz_results_<job>/
          predictions/
            <job>/
              confidence_<job>_model_0.json   summary scores (iptm, ptm, complex_plddt, ...)
              pae_<job>_model_0.npz           PAE matrix (per-token)
              plddt_<job>_model_0.npz         per-residue pLDDT
              pde_<job>_model_0.npz           predicted distance error (unused)
              <job>_model_0.cif               predicted structure
              affinity_<job>_model_0.json     ONLY if affinity head was requested
                                              via `properties: [{affinity: ...}]` in YAML
    """
    out: dict = {}

    # Step into boltz_results_<job>/predictions/<job>/
    results_dirs = sorted(run_dir.glob("boltz_results_*"))
    if not results_dirs:
        return out
    pred_root = results_dirs[0] / "predictions"
    if not pred_root.exists():
        return out
    job_dirs = [p for p in pred_root.iterdir() if p.is_dir()]
    if not job_dirs:
        return out
    job_dir  = job_dirs[0]
    job_name = job_dir.name

    # ── Confidence JSON (iptm, ptm, summary plddt) ──
    conf_path = job_dir / f"confidence_{job_name}_model_0.json"
    if conf_path.exists():
        try:
            with open(conf_path) as f:
                conf = json.load(f)
            out["iptm"] = conf.get("iptm", conf.get("complex_iptm", conf.get("protein_iptm")))
            out["ptm"]  = conf.get("ptm",  conf.get("complex_ptm"))
            plddt = conf.get("complex_plddt", conf.get("plddt"))
            if plddt is not None:
                plddt = float(plddt)
                if 0.0 < plddt <= 1.0:
                    plddt *= 100.0
                out["plddt"] = plddt
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

    # ── Binder-only pLDDT (more informative than complex pLDDT) ──
    # Prefer per-residue plddt npz; fall back to complex_plddt from confidence JSON
    plddt_path = job_dir / f"plddt_{job_name}_model_0.npz"
    cif_path   = job_dir / f"{job_name}_model_0.cif"
    if plddt_path.exists():
        try:
            arr = np.load(plddt_path, allow_pickle=False)
            plddt_arr = (arr["plddt"] if "plddt" in arr.files else arr[arr.files[0]]).flatten()

            n_binder_res = binder_len  # default assumption
            if cif_path.exists():
                try:
                    from biotite.structure.io.pdbx import CIFFile, get_structure
                    struct = get_structure(CIFFile.read(str(cif_path)), model=1)
                    ca = struct[struct.atom_name == "CA"]
                    n_binder_res = int((ca.chain_id == binder_chain).sum()) or binder_len
                except Exception:
                    pass

            # plddt may be per-residue (len == total residues) or per-atom (longer)
            if len(plddt_arr) >= n_binder_res:
                binder_mean = float(plddt_arr[:n_binder_res].mean())
                if 0.0 < binder_mean <= 1.0:
                    binder_mean *= 100.0
                out["plddt"] = binder_mean   # overrides complex_plddt with binder-only
        except Exception as e:
            logger.warning(f"Failed to read Boltz pLDDT: {e}")

    # ── Affinity head (Boltz-2 only; requires explicit request in YAML `properties`) ──
    aff_path = job_dir / f"affinity_{job_name}_model_0.json"
    if aff_path.exists():
        try:
            with open(aff_path) as f:
                aff = json.load(f)
            out["affinity"] = aff.get(
                "affinity_pred", aff.get("affinity", aff.get("delta_g",
                                          aff.get("affinity_pred_value")))
            )
        except Exception as e:
            logger.warning(f"Failed to read Boltz affinity: {e}")

    # ── Structure path ──
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

    # ── Target template (target-chain-only CIF, cached) ───────────────────────

    def _ensure_target_template(self, target_name: str) -> Path:
        """
        Extract just the target chain from the source HADDOCK PDB and save as a
        CIF that RF3 can consume as a template.  Cached per target.

        We deliberately do NOT template the binder chain — that would bias RF3
        toward the native TIMP3 conformation, defeating the point of redesign.
        """
        tcfg     = TARGETS[target_name]
        tpl_dir  = OUT_BASE / "target_templates"
        tpl_dir.mkdir(parents=True, exist_ok=True)
        cif_path = tpl_dir / f"{target_name}_target_only.cif"
        if cif_path.exists():
            return cif_path

        src_pdb = DATA_DIR / tcfg["pdb"]
        arr     = PDBFile.read(str(src_pdb)).get_structure()[0]
        keep    = arr[arr.chain_id == tcfg["target_chain"]]
        keep    = renumber(keep)
        to_cif_file(keep, str(cif_path), file_type="cif")
        logger.info(f"[{target_name}] target template cached: {cif_path}")
        return cif_path

    # ── Boltz-2 (local AF3-class ranker) ──────────────────────────────────────

    def run_boltz(self, target_name: str, candidates: list, out_dir: Path) -> list:
        """
        Score each candidate with Boltz-2.  Boltz-2 is an open-source AF3-class
        model used here as the local pre-AF3 ranker — much better correlated
        with AF3 ipTM than RF3 single-sequence (calibrated empirically).

        Adds boltz_iptm / boltz_ptm / boltz_plddt / boltz_iface_pae /
        boltz_affinity fields to each candidate and recomputes composite_score.
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

            metrics = _parse_boltz_output(run_dir, bc, fc, len(bseq))
            if metrics:
                cand.update({
                    "boltz_iptm":      metrics.get("iptm"),
                    "boltz_ptm":       metrics.get("ptm"),
                    "boltz_plddt":     metrics.get("plddt"),
                    "boltz_iface_pae": metrics.get("iface_pae"),
                    "boltz_affinity":  metrics.get("affinity"),
                    "boltz_cif":       metrics.get("cif_path"),
                    "source":          "Boltz",  # promote from "RF3"
                })
                cand["composite_score"] = calc_composite(cand)
                n_done += 1

        logger.info(
            f"[{target_name}] Boltz-2 done in {(time.time()-t0)/60:.1f} min "
            f"({n_done}/{len(candidates)} scored)"
        )
        return candidates

    # ── RFd3 ─────────────────────────────────────────────────────────────────

    def _build_contig(self, tcfg: dict, fix_chain_len: int,
                      adaptive_ranges: dict = None) -> tuple[str, str]:
        """
        Build contig string and length range for RFd3.
        Returns (contig_string, length_range_string).
        """
        bc     = tcfg["binder_chain"]
        fc     = tcfg["target_chain"]
        total  = tcfg["scaffold_len"]
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

        contig, length_range = self._build_contig(tcfg, fix_len, adaptive_ranges)
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
                    "batch_size":     1,
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
                    to_cif_file(arr, str(out_dir / f"{design_id}_s{si}.cif"), file_type="cif")
                    results.append({
                        "design_id":   f"{design_id}_s{si}",
                        "target_name": target_name,
                        "array":       arr,
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

    # ── RF3 complex scoring ───────────────────────────────────────────────────

    def run_rf3_complex(self, target_name: str, candidates: list,
                        out_dir: Path) -> list:
        """
        Fold each candidate as a two-chain complex (binder + target) using RF3
        with sequence-only input.  Extracts pLDDT, ptm, ipTM, interface PAE,
        backbone RMSD, and interface contacts.  Saves CIF + metrics JSON.
        Returns list of scored metadata dicts (no atom arrays — too large to cache).
        """
        if not candidates:
            return []

        # Use the RFd3 OUTPUT chain convention (binder=A, target=B); the source-PDB
        # chain IDs in TARGETS apply only to reading target_seqs and building contigs.
        bc   = DESIGN_BINDER_CHAIN
        fc   = DESIGN_TARGET_CHAIN
        out_dir.mkdir(parents=True, exist_ok=True)

        target_seq = self.target_seqs.get(target_name, "")
        if not target_seq:
            logger.warning(f"No target sequence for {target_name}; RF3 complex skipped.")
            return []

        # NOTE: an earlier attempt passed AF3-style `useStructureTemplate` +
        # `templates` fields per chain, but this RF3 build's SequenceComponent
        # rejects them with `unexpected keyword argument 'useStructureTemplate'`.
        # Until we identify the correct RF3 template API (likely a kwarg on
        # engine.run() or a separate templates= parameter on InferenceInput),
        # RF3 runs single-sequence.  Boltz-2 remains the trusted local ranker;
        # RF3 contributes geometric features (n_contacts, RMSD, iface PAE) only.
        # _ensure_target_template still caches the target-only CIF on disk so
        # it's ready the moment we wire up the right API.
        self._ensure_target_template(target_name)

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
                    **{k: v for k, v in cand.items() if k not in ("array",)},
                    **metrics_rec,
                    "iteration":    self.state["iteration"],
                    "temperature":  self.state["temperature"],
                    "rf3_cif":      str(cif_path),
                    "promising":    iptm >= IPTM_PROMISING,
                })

            except Exception as exc:
                logger.error(f"RF3 error on {did}: {exc}")

        logger.info(f"[{target_name}] RF3 done in {(time.time()-t0)/60:.1f} min "
                    f"({len(scored)}/{len(candidates)} scored)")
        del engine
        torch.cuda.empty_cache()
        return scored

    # ── Hall of Fame ──────────────────────────────────────────────────────────

    def update_hof(self, new_scored: list) -> None:
        """Add new scored designs to per-target HOF; keep top HOF_SIZE_PER_TARGET."""
        for entry in new_scored:
            tname = entry.get("target_name")
            if tname not in self.state["hof"]:
                self.state["hof"][tname] = []
            self.state["hof"][tname].append(entry)

        for tname in self.active_targets:
            hof = self.state["hof"].get(tname, [])
            self.state["hof"][tname] = sorted(
                hof, key=lambda x: x.get("composite_score", 0), reverse=True
            )[:HOF_SIZE_PER_TARGET]

        # Save best HOF structures for later seeding/inspection
        self._save_hof_structures()
        self._write_hof_summary()

    def _save_hof_structures(self) -> None:
        """
        Copy top-3 RF3 CIFs per target into hof_structures/ if not already there.
        These can be used as visual reference or seeded into future RFd3 runs.
        """
        for tname in self.active_targets:
            tdir = self.hof_struct_dir / tname
            tdir.mkdir(parents=True, exist_ok=True)
            for rank, entry in enumerate(self.state["hof"].get(tname, [])[:3], start=1):
                src = entry.get("rf3_cif")
                if src and Path(src).exists():
                    dest = tdir / f"hof_rank{rank:02d}_{entry['design_id']}.cif"
                    if not dest.exists():
                        import shutil
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

                    # Fallback: find target by sequence match across all active HOFs
                    if target_name not in self.active_targets:
                        target_name = next(
                            (t for t in self.active_targets
                             if any(e.get("full_seq") == binder_seq
                                    for e in self.state["hof"].get(t, []))),
                            None,
                        )

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

                    # Match to HOF entry by sequence; fall back to job name as design_id
                    design_id = job_name
                    for e in self.state["hof"].get(target_name, []):
                        if e.get("full_seq") == binder_seq:
                            design_id = e["design_id"]
                            break

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
            tcfg     = TARGETS[tname]
            pdb_path = str(DATA_DIR / tcfg["pdb"])
            if not Path(pdb_path).exists():
                logger.warning(f"Skipping {tname}: PDB not found at {pdb_path}")
                continue

            # 1. Backbone generation
            rfd3_dir = it_dir / "rfd3" / tname
            rfd3_dir.mkdir(parents=True, exist_ok=True)
            backbones = self.run_rfd3(tname, pdb_path, BACKBONES_PER_TARGET, adaptive_ranges)

            # Save backbone CIFs
            for bi, arr in enumerate(backbones):
                to_cif_file(arr, str(rfd3_dir / f"bb_{bi}.cif"), file_type="cif")

            if not backbones:
                logger.warning(f"[{tname}] No backbones generated; skipping.")
                continue

            # 2. Sequence design
            lmpnn_dir  = it_dir / "lmpnn" / tname
            candidates = self.run_lmpnn(tname, backbones, lmpnn_dir, temp)

            if not candidates:
                logger.warning(f"[{tname}] No LMPNN sequences; skipping RF3.")
                continue

            # 3. RF3 complex scoring (templated, sanity check + geometric features)
            rf3_dir = it_dir / "rf3" / tname
            scored  = self.run_rf3_complex(tname, candidates, rf3_dir)

            # 4. Boltz-2 local AF3-class scoring — the real pre-AF3 ranker.
            # Adds boltz_* fields and promotes source from "RF3" to "Boltz".
            # If BOLTZ_ENABLE is False or boltz CLI is missing, this is a no-op
            # and the pipeline continues with RF3-only geometric composites.
            boltz_dir = it_dir / "boltz" / tname
            scored    = self.run_boltz(tname, scored, boltz_dir)

            # Drop the heavy atom arrays before HOF storage
            for s in scored:
                s.pop("rfd3_array", None)
                s.pop("array", None)
            all_scored.extend(scored)

        # Update HOF and write summaries
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
        "--import-af3", type=str, default=None,
        help=(
            "Path to AF3 results to import before continuing.  "
            "Accepts either a batch ZIP downloaded directly from the AF3 Server "
            "(one subdirectory per job) or a hand-crafted JSON list of dicts with keys: "
            "design_id, iptm, plddt, ptm, interface_pae, full_seq, target_name."
        ),
    )
    args = parser.parse_args()

    refiner = IterativeRefiner(
        active_targets=args.targets,
        active_loops=args.loops,
    )

    if args.import_af3:
        if args.import_af3.lower().endswith(".zip"):
            refiner.import_af3_zip(args.import_af3)
        else:
            refiner.import_af3_results(args.import_af3)

    refiner.main_loop(max_iterations=args.max_iterations)


if __name__ == "__main__":
    main()
