"""
iterative_refinement.py

Iterative binder design: RFd3 → LigandMPNN → RF3 (complex) with temperature annealing.

Generates and refines TIMP3-scaffold binders for MMP2, MMP9, ADAM10, and ADAM17 by
redesigning loop regions (AB, C, EF). Each iteration scores all designs against the
target using RF3 complex-mode prediction to obtain ipTM, pLDDT, interface PAE, and
backbone RMSD. The best performers are collected into a per-target Hall of Fame; loop
length ranges for subsequent iterations are narrowed toward successful lengths (adaptive
contig bias). Temperature is annealed from 0.5 → 0.1 to shift from exploration to
exploitation. Every AF3_EXPORT_EVERY_N iterations the top AF3_TOP_N designs are written
to an AF3 Server JSON for manual submission; AF3 results can be imported to update the
Hall of Fame with higher-quality metrics.

Output layout:
  Local/iterative_refinement/
    refinement_state.json             persistent state (HOF, temperature, iteration)
    it_N/
      rfd3/                           RFd3 backbone CIFs
      lmpnn/                          LigandMPNN sequence CIFs
      rf3/                            RF3 complex prediction CIFs + metrics JSONs
      round_summary.csv               all scored designs this iteration
    hof_structures/<target>/          RF3 CIFs for HOF designs (for seeding)
    af3_submission_itN.json           AF3 Server input batch
    hof_summary.csv                   all-time best per target (updated each round)
"""

import os
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
DATA_DIR  = _HERE / ".." / "Data" / "TIMP_Complexes" / "HADDOCK_PDB"
OUT_BASE  = _HERE / ".." / "Local" / "iterative_refinement"
CKPT_DIR  = _HERE / ".." / "Tools" / "foundry_checkpoints"

# ── Target definitions ────────────────────────────────────────────────────────
# binder_chain: chain being redesigned (TIMP3); target_chain: fixed protease chain.
TARGETS = {
    "MMP2":   {"pdb": "TIMP3_vs_MMP2_HADDOCK_Xray.pdb",   "binder_chain": "A", "target_chain": "B", "scaffold_len": 121},
    "MMP9":   {"pdb": "TIMP3_vs_MMP9_HADDOCK_Xray.pdb",   "binder_chain": "A", "target_chain": "B", "scaffold_len": 121},
    "ADAM10": {"pdb": "TIMP3_vs_ADAM10_HADDOCK_Xray.pdb",  "binder_chain": "A", "target_chain": "B", "scaffold_len": 121},
    "ADAM17": {"pdb": "TIMP3_vs_ADAM17_HADDOCK_Xray.pdb",  "binder_chain": "A", "target_chain": "B", "scaffold_len": 121},
}

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
BACKBONES_PER_TARGET   = 10    # RFd3 designs per target per iteration
INIT_TEMPERATURE       = 0.50  # LMPNN sampling temperature (exploration)
MIN_TEMPERATURE        = 0.10  # LMPNN temperature floor (exploitation)
TEMP_DECAY             = 0.85  # Temperature multiplier applied each iteration
HOF_SIZE_PER_TARGET    = 75    # Max Hall of Fame entries per target
ADAPTIVE_BIAS_START    = 3     # First iteration to apply adaptive loop-length bias
ADAPTIVE_BIAS_PCT      = 25    # Use top-N% HOF lengths to define new contig range
AF3_EXPORT_EVERY_N     = 5     # Export AF3 JSON every N full iterations
AF3_TOP_N              = 20    # Designs to include per AF3 submission
IPTM_PROMISING         = 0.55  # ipTM threshold to flag a design as "promising"
RMSD_CLIP              = 5.0   # Å beyond which RMSD contribution scores 0
PAE_MAX                = 30.0  # Å² beyond which interface PAE contribution scores 0

# Composite score weights (must sum to 1.0)
COMPOSITE_WEIGHTS = {
    "iptm":  0.40,  # Interface quality — primary binding predictor
    "plddt": 0.20,  # Binder fold confidence
    "ptm":   0.15,  # Global structure confidence
    "rmsd":  0.15,  # Backbone fidelity (lower RMSD → higher score)
    "pae":   0.10,  # Interface PAE (lower → higher score)
}


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


def calc_composite(iptm: float, plddt: float, ptm: float,
                   rmsd_val: float, iface_pae: float) -> float:
    """Weighted composite score ∈ [0, 1]. Higher is better."""
    w = COMPOSITE_WEIGHTS
    rmsd_s = max(0.0, 1.0 - rmsd_val / RMSD_CLIP)
    pae_s  = max(0.0, 1.0 - iface_pae / PAE_MAX) if not np.isnan(iface_pae) else 0.0
    return (
        w["iptm"]  * float(iptm)
        + w["plddt"] * (float(plddt) / 100.0)
        + w["ptm"]   * float(ptm)
        + w["rmsd"]  * rmsd_s
        + w["pae"]   * pae_s
    )


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
        """
        tcfg   = TARGETS[target_name]
        bc     = tcfg["binder_chain"]
        fc     = tcfg["target_chain"]
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

        tcfg = TARGETS[target_name]
        bc   = tcfg["binder_chain"]
        fc   = tcfg["target_chain"]
        out_dir.mkdir(parents=True, exist_ok=True)

        target_seq = self.target_seqs.get(target_name, "")
        if not target_seq:
            logger.warning(f"No target sequence for {target_name}; RF3 complex skipped.")
            return []

        engine  = RF3InferenceEngine(ckpt_path="rf3", verbose=False)
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

                # Confidence metrics
                conf  = rf3_out.summary_confidences or {}
                plddt = conf.get("overall_plddt", conf.get("plddt", 0.0))
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

                comp = calc_composite(iptm, plddt, ptm, rmsd_val, iface_pae)

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
                    "composite_score": comp,
                }
                with open(out_dir / f"{did}_metrics.json", "w") as mf:
                    json.dump(metrics_rec, mf, indent=2)

                scored.append({
                    **{k: v for k, v in cand.items() if k not in ("array", "rfd3_array")},
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
        Export the top AF3_TOP_N designs (across all targets) to an AF3 Server
        JSON file.  Designs with AF3 source already in HOF are de-prioritised
        so new candidates always surface.
        """
        it = self.state["iteration"]
        if not force and (it - self.state.get("last_af3_it", -1)) < AF3_EXPORT_EVERY_N:
            return

        # Gather all HOF entries, prefer RF3-validated over AF3-already-seen
        all_entries = []
        for tname in self.active_targets:
            for e in self.state["hof"].get(tname, []):
                all_entries.append(e)

        if not all_entries:
            logger.warning("HOF empty; skipping AF3 export.")
            return

        seen_seqs = set()
        unique_entries = []
        for e in sorted(all_entries, key=lambda x: x.get("composite_score", 0), reverse=True):
            seq = e.get("full_seq", "")
            if seq and seq not in seen_seqs:
                seen_seqs.add(seq)
                unique_entries.append(e)
            if len(unique_entries) >= AF3_TOP_N:
                break

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

            # 3. RF3 complex scoring
            rf3_dir = it_dir / "rf3" / tname
            scored  = self.run_rf3_complex(tname, candidates, rf3_dir)
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
        choices=list(TARGETS.keys()),
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
        help="Path to an AF3 results JSON to import before continuing.",
    )
    args = parser.parse_args()

    refiner = IterativeRefiner(
        active_targets=args.targets,
        active_loops=args.loops,
    )

    if args.import_af3:
        refiner.import_af3_results(args.import_af3)

    refiner.main_loop(max_iterations=args.max_iterations)


if __name__ == "__main__":
    main()
