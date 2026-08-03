"""
specificity_refinement.py

Specificity-focused iterative binder design for TIMP3-scaffold loops.

Extends iterative_refinement by scoring every designed binder against BOTH the
on-target and off-target protease in a specificity pair:
  - MMP2 (on-target)   vs MMP9 (off-target)
  - ADAM10 (on-target) vs ADAM17 (off-target)

Scoring uses ESMFold2 (the validated local ranker — see filter_methods.md),
NOT RF3 (RF3 confidence was anti-correlated with AF3). Each binder is folded
against the on- and off-target; we record:
  selectivity_score = esm_iptm_on - esm_iptm_off   (>0 prefers on-target)
  composite         = reward on-target binding (ipTM_on, pLDDT_on) + selectivity

Inherits all the main-pipeline machinery: RFd3 → LigandMPNN → (RF3 off) →
ESMFold2, per-target failure isolation, AF3-protected/deduped HOF, multi-GPU
ESMFold2, and AF3 ZIP import. RFd3/LigandMPNN design against the ON-target only;
off-target sequences are used purely for cross-scoring.

Usage:
  python specificity_refinement.py --pair MMP --loops AB C EF
  python specificity_refinement.py --pair MMP ADAM --esmfold2-gpus auto
  python specificity_refinement.py --pair MMP --import-af3 <results.zip>

Output: Local/specificity_refinement/
"""

import os
import subprocess

# --- PORTABILITY: GPU-Aware Environment Setup ---
# Detect the GPU and set DISABLE_CUEQUIVARIANCE before importing heavy ML libs
# (RFd3/LigandMPNN still run in-process in the foundry env).
try:
    smi_out = subprocess.check_output(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"]).decode()
    if "V100" in smi_out:
        os.environ["DISABLE_CUEQUIVARIANCE"] = "1"
        print("Detected V100 GPU. Automatically setting DISABLE_CUEQUIVARIANCE=1 for compatibility.")
except Exception:
    pass

import re
import json
import logging
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

import iterative_refinement as ir  # for setting module globals (ESMFOLD2_GPUS)
from iterative_refinement import (
    IterativeRefiner,
    TARGETS,
    LOOP_CONFIGS,
    DATA_DIR,
    AF3_EXPORT_EVERY_N,
    AF3_TOP_N,
    IPTM_PROMISING,
    ESMFOLD2_ENABLE,
    pdb_chain_seq,
    extract_loops,
    _esm_iface_feats,
    _safe,
    _normalize_plddt,
)
import calibrated_scoring as cs  # calibrated priors + interface geometry (2026-07)

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)
for _noisy in ("transforms", "atomworks.io", "atomworks.ml", "foundry", "lightning"):
    logging.getLogger(_noisy).setLevel(logging.ERROR)

# ── Paths ─────────────────────────────────────────────────────────────────────
_HERE    = Path(__file__).parent.resolve()
OUT_BASE = Path(os.environ.get("SPEC_OUT_BASE") or (_HERE / ".." / "Local" / "specificity_refinement"))

# CRITICAL ISOLATION FIX: the inherited IterativeRefiner methods (run_iteration,
# _write_round_summary, update_hof, export_for_af3, ...) reference iterative_
# refinement's MODULE-GLOBAL OUT_BASE — NOT this one. Without redirecting it,
# specificity's it_N / round_summary.csv / hof_summary.csv land in
# Local/iterative_refinement and COLLIDE with a concurrent generation run (this
# happened 07-28: the iterative dir became a generation+specificity mix). Point the
# parent module's OUT_BASE at ours so ALL specificity output stays in one place.
ir.OUT_BASE = OUT_BASE

# ── Specificity pair definitions ──────────────────────────────────────────────
SPECIFICITY_PAIRS = {
    "MMP":  {"on_target": "MMP2",  "off_target": "MMP9"},
    "ADAM": {"on_target": "ADAM10", "off_target": "ADAM17"},
}

# Specificity composite (2026-07 recalibration; see calibrated_scoring.py).
#   on_quality  = the calibrated on-target score (ESMFold2 stage score, or the AF3
#                 binding prior) — foldability floor + interface geometry. It does
#                 NOT reward esm_iptm, whose negative in-sample correlation with
#                 binding is a selection-bias artifact.
#   selectivity = the normalized interface CONTACT-DENSITY gap (on - off). Contact
#                 density is the one ESMFold2 feature that tracked binding, and a
#                 within-design on/off contrast is the least bias-exposed way to
#                 read selectivity. Neutral (0.5) when the geometry isn't available.
SPECIFICITY_COMPOSITE_WEIGHTS = {"on_quality": 0.50, "selectivity": 0.50}
CONTACT_DENSITY_GAP_SCALE = 20.0   # a contact-density gap of this magnitude ~ saturates selectivity
HOF_SIZE = 75   # per specificity pair


def _cd_gap_norm(cd_on, cd_off) -> float:
    """Normalize an on-minus-off contact-density gap to [0,1]; 0.5 = neutral/unknown."""
    if not (cs._isnum(cd_on) and cs._isnum(cd_off)):
        return 0.5
    g = (float(cd_on) - float(cd_off)) / CONTACT_DENSITY_GAP_SCALE
    return float(min(1.0, max(0.0, (g + 1.0) / 2.0)))


def calc_specificity_composite(on_metrics: dict, off_metrics: dict = None,
                               model: str = "esm") -> float:
    """
    Composite in [0,1] rewarding on-target QUALITY (calibrated prior) AND
    SELECTIVITY (contact-density gap). `on_metrics`/`off_metrics` are metric dicts
    keyed with af3_*/esm_* names (see calibrated_scoring). `model` picks the
    on-quality scorer: "esm" -> esmfold2_stage_score, "af3" -> af3_binding_prior.
    """
    w = SPECIFICITY_COMPOSITE_WEIGHTS
    off_metrics = off_metrics or {}
    if model == "af3":
        on_q = cs.af3_binding_prior(on_metrics)
        cd_on, cd_off = on_metrics.get("af3_iface_contact_density"), off_metrics.get("af3_iface_contact_density")
    else:
        on_q = cs.esmfold2_stage_score(on_metrics)
        cd_on, cd_off = on_metrics.get("esm_iface_contact_density"), off_metrics.get("esm_iface_contact_density")
    if not (on_q == on_q):   # NaN guard
        on_q = 0.0
    return float(w["on_quality"] * on_q + w["selectivity"] * _cd_gap_norm(cd_on, cd_off))


class SpecificityRefiner(IterativeRefiner):
    """
    Score each design against both the on- and off-target (ESMFold2) and rank by
    a selectivity-aware composite. Design happens against the on-target only.
    """

    def __init__(self, pair_keys: list, active_loops: list):
        self.pair_keys = pair_keys
        self.pairs     = {pk: SPECIFICITY_PAIRS[pk] for pk in pair_keys}
        # Design (and the per-target HOF) operate on the ON-targets only.
        on_targets = list({self.pairs[pk]["on_target"] for pk in pair_keys})
        self._pair_for_on_target = {self.pairs[pk]["on_target"]: pk for pk in pair_keys}

        super().__init__(
            active_targets=on_targets,
            active_loops=active_loops,
            state_path=OUT_BASE / "specificity_state.json",
        )

        # Per-pair specificity HOF (separate from the inherited per-target HOF)
        if "specificity_hof" not in self.state:
            self.state["specificity_hof"] = {pk: [] for pk in pair_keys}

        # Load off-target sequences for cross-scoring (not designed against)
        for pk, pair in self.pairs.items():
            off = pair["off_target"]
            if off not in self.target_seqs:
                pdb_path = DATA_DIR / TARGETS[off]["pdb"]
                if pdb_path.exists():
                    self.target_seqs[off] = pdb_chain_seq(str(pdb_path), TARGETS[off]["target_chain"])

    # ── ESMFold2 dual (on + off) scoring — replaces the old RF3 override ───────

    def run_esmfold2(self, target_name: str, candidates: list, out_dir: Path) -> list:
        """
        Score each candidate against the on-target AND off-target with ESMFold2
        (one combined, GPU-sharded batch via the inherited helper), then rank by
        the selectivity-aware composite and update the specificity HOF.
        """
        if not ESMFOLD2_ENABLE or not candidates:
            return candidates

        pk = self._pair_for_on_target.get(target_name)
        if pk is None:   # not an on-target (shouldn't happen) → standard scoring
            return super().run_esmfold2(target_name, candidates, out_dir)

        on_t, off_t = self.pairs[pk]["on_target"], self.pairs[pk]["off_target"]
        on_seq  = self.target_seqs.get(on_t, "")
        off_seq = self.target_seqs.get(off_t, "")

        # One batch with two rows per design (::on / ::off); helper shards across GPUs.
        rows = []
        for c in candidates:
            bseq = c.get("full_seq", "")
            if not bseq:
                continue
            if on_seq:
                rows.append({"design_id": f"{c['design_id']}::on",  "target_name": on_t,
                             "full_seq": bseq, "target_seq": on_seq})
            if off_seq:
                rows.append({"design_id": f"{c['design_id']}::off", "target_name": off_t,
                             "full_seq": bseq, "target_seq": off_seq})

        t0 = time.time()
        scores = self._score_sequences_esmfold2(rows, out_dir)

        n_done = 0
        for c in candidates:
            on_m  = scores.get(f"{c['design_id']}::on")
            off_m = scores.get(f"{c['design_id']}::off")
            if not on_m:
                continue
            iptm_on  = on_m["esm_iptm"]
            plddt_on = on_m["esm_plddt"]
            iptm_off = off_m["esm_iptm"] if off_m else float("nan")
            sel = _safe(iptm_on) - (_safe(iptm_off) if not np.isnan(_safe(iptm_off, np.nan)) else 0.0)

            # Interface geometry from the on- and off-target ESMFold2 complexes:
            # contact density drives the (recalibrated) selectivity term.
            on_feats  = _esm_iface_feats(on_m.get("esm_cif"))
            off_feats = _esm_iface_feats(off_m.get("esm_cif")) if off_m else {}
            cd_on  = on_feats.get("esm_iface_contact_density")
            cd_off = off_feats.get("esm_iface_contact_density")

            c.update({
                "esm_iptm_on":   iptm_on,
                "esm_plddt_on":  plddt_on,
                "esm_ptm_on":    on_m.get("esm_ptm"),
                "esm_iptm_off":  iptm_off,
                "esm_plddt_off": off_m.get("esm_plddt") if off_m else float("nan"),
                "selectivity_score": sel,                      # legacy ipTM gap (display only)
                "cd_on": cd_on, "cd_off": cd_off,
                "selectivity_cd": ((cd_on - cd_off) if (cs._isnum(cd_on) and cs._isnum(cd_off)) else float("nan")),
                # esm_iptm/esm_plddt aliases = on-target (so inherited code/CSVs read sensibly)
                "esm_iptm":   iptm_on,
                "esm_plddt":  plddt_on,
                "esm_iface_contact_density": cd_on,
                "esm_iface_n_iface_res":     on_feats.get("esm_iface_n_iface_res"),
                "esm_cif":    on_m.get("esm_cif"),
                "specificity_pair": pk,
                "on_target": on_t, "off_target": off_t,
                "source": "ESMFold2",
                # "promising" now means: on-target folds+docks AND geometry prefers on-target
                "promising": (cs.esm_passes_fold_gate({"esm_plddt": plddt_on,
                                                       "esm_iface_n_iface_res": on_feats.get("esm_iface_n_iface_res")})
                              and (_safe(c.get("selectivity_cd"), 0.0) > 0)),
            })
            c["composite_score"] = calc_specificity_composite(
                {"esm_plddt": plddt_on, "esm_iface_contact_density": cd_on,
                 "esm_iface_n_iface_res": on_feats.get("esm_iface_n_iface_res")},
                {"esm_iface_contact_density": cd_off}, model="esm")
            self.state["specificity_hof"].setdefault(pk, []).append(c)
            n_done += 1

        # Consolidate specificity HOF with the same dedup + AF3-protection as the parent
        for pk2 in self.pair_keys:
            self.state["specificity_hof"][pk2] = self._consolidate_hof(
                self.state["specificity_hof"].get(pk2, []), HOF_SIZE
            )

        logger.info(f"[{on_t}↔{off_t}] ESMFold2 specificity scoring done in "
                    f"{(time.time()-t0)/60:.1f} min ({n_done}/{len(candidates)} scored)")
        return candidates

    # ── Specificity HOF summary ───────────────────────────────────────────────

    def _write_hof_summary(self) -> None:
        super()._write_hof_summary()
        rows = []
        for pk in self.pair_keys:
            for rank, e in enumerate(self.state["specificity_hof"].get(pk, []), start=1):
                rows.append({"pair": pk, "rank": rank,
                             **{k: v for k, v in e.items() if k not in ("array", "rfd3_array")}})
        if rows:
            pd.DataFrame(rows).to_csv(OUT_BASE / "specificity_hof_summary.csv", index=False)
            for pk in self.pair_keys:
                pkrows = [r for r in rows if r["pair"] == pk]
                if pkrows:
                    b = pkrows[0]
                    logger.info(
                        f"Specificity HOF best [{pk}]: {b.get('design_id')} | "
                        f"ipTM_on={_safe(b.get('esm_iptm_on')):.3f} | "
                        f"ipTM_off={_safe(b.get('esm_iptm_off')):.3f} | "
                        f"selectivity={_safe(b.get('selectivity_score')):.3f} | "
                        f"src={b.get('source')}"
                    )

    # ── AF3 export (on-target AND off-target jobs for selectivity validation) ──

    def export_for_af3(self, force: bool = False) -> None:
        it = self.state["iteration"]
        if not force and (it - self.state.get("last_af3_it", -1)) < AF3_EXPORT_EVERY_N:
            return

        # Each design needs an on AND an off job, so it consumes TWO of the daily
        # 30 slots — budget designs accordingly.
        designs_per_pair = max(1, AF3_TOP_N // (2 * max(1, len(self.pair_keys))))

        jobs, seen = [], set()
        for pk in self.pair_keys:
            on_t, off_t = self.pairs[pk]["on_target"], self.pairs[pk]["off_target"]
            on_seq, off_seq = self.target_seqs.get(on_t, ""), self.target_seqs.get(off_t, "")
            picked = 0
            for e in self.state["specificity_hof"].get(pk, []):
                bseq = e.get("full_seq", "")
                if not bseq or bseq in seen:
                    continue
                seen.add(bseq)
                idx = len(jobs)
                jobs.append({"name": f"spec_it{it}_{pk}_on_{idx:02d}", "modelSeeds": [42],
                             "sequences": [{"proteinChain": {"sequence": bseq,   "count": 1}},
                                           {"proteinChain": {"sequence": on_seq, "count": 1}}]})
                jobs.append({"name": f"spec_it{it}_{pk}_off_{idx:02d}", "modelSeeds": [42],
                             "sequences": [{"proteinChain": {"sequence": bseq,    "count": 1}},
                                           {"proteinChain": {"sequence": off_seq, "count": 1}}]})
                picked += 1
                if picked >= designs_per_pair:
                    break
            logger.info(f"AF3 specificity export: {pk} contributing {picked} design(s)")

        if not jobs:
            return
        export_path = OUT_BASE / f"af3_specificity_it{it}.json"
        export_path.write_text(json.dumps(jobs, indent=2))
        self.state["last_af3_it"] = it
        self._save_state()

        logger.info(f"AF3 specificity export: {len(jobs)} jobs → {export_path}")
        print("\n" + "=" * 60)
        print(f"  AF3 specificity submission: {export_path}")
        print(f"  {len(jobs)} jobs (on + off per design) for selectivity validation.")
        print(f"  Import results via --import-af3 <results.zip>")
        print("=" * 60 + "\n")

    # ── AF3 ZIP import (pairs on/off jobs, recomputes selectivity from AF3) ────

    def import_af3_zip(self, zip_path: str) -> None:
        """
        Parse a specificity AF3 ZIP whose jobs are named spec_it{N}_{PK}_{on|off}_{NN}.
        Pairs the on/off jobs for each binder, recomputes the selectivity composite
        from AF3 ipTM, and updates the specificity HOF (source="AF3", never trimmed).
        """
        # group[(pk, binder_seq)] = {"on": {...}, "off": {...}}
        group: dict = {}
        with zipfile.ZipFile(zip_path) as zf:
            names = set(zf.namelist())
            for jrf in sorted(n for n in names if n.endswith("_job_request.json")):
                try:
                    raw = json.loads(zf.read(jrf))
                    jd  = raw[0] if isinstance(raw, list) else raw
                    name = jd.get("name", "")
                    m = re.search(r"(?:fold_)?spec_it\d+_([A-Za-z0-9]+)_(on|off)_\d+", name, re.I)
                    if not m:
                        continue
                    pk, side = m.group(1).upper(), m.group(2).lower()
                    if pk not in self.pairs:
                        continue
                    seqs = jd.get("sequences", [])
                    if len(seqs) < 2:
                        continue
                    binder_seq = seqs[0]["proteinChain"]["sequence"]

                    prefix  = jrf[: -len("_job_request.json")]
                    sc_name = prefix + "_summary_confidences_0.json"
                    if sc_name not in names:
                        continue
                    sc   = json.loads(zf.read(sc_name))
                    iptm = float(sc.get("iptm", 0.0))
                    ptm  = float(sc.get("ptm", 0.0))
                    cp   = sc.get("chain_ptm") or []
                    aptm = float(cp[0]) if len(cp) > 0 else float("nan")  # binder-chain pTM (drop BpTM)

                    plddt = float("nan")
                    fd_name = prefix + "_full_data_0.json"
                    if fd_name in names:
                        fd = json.loads(zf.read(fd_name))
                        ap = np.array(fd.get("atom_plddts", []))
                        ac = fd.get("atom_chain_ids", [])
                        if ap.size and ac:
                            mask = np.array([c == "A" for c in ac])   # chain A = binder
                            if mask.any():
                                plddt = float(ap[mask].mean())

                    # Interface geometry (same extractor as calibration)
                    af3_iface = {}
                    cif_name = prefix + "_model_0.cif"
                    if cif_name in names:
                        try:
                            af3_iface = cs.interface_features_from_cif(zf.read(cif_name).decode())
                        except Exception:
                            af3_iface = {}

                    group.setdefault((pk, binder_seq), {})[side] = {
                        "iptm": iptm, "ptm": ptm, "plddt": plddt, "aptm": aptm,
                        "n_iface_res":     af3_iface.get("n_iface_res"),
                        "iface_plddt":     af3_iface.get("iface_plddt"),
                        "contact_density": af3_iface.get("contact_density"),
                    }
                except Exception as exc:
                    logger.error(f"Error parsing {jrf}: {exc}")

        n = 0
        for (pk, bseq), sides in group.items():
            on_m = sides.get("on")
            if not on_m:
                continue   # need at least the on-target result
            off_m = sides.get("off", {})
            iptm_on  = on_m["iptm"]
            plddt_on = on_m["plddt"]
            iptm_off = off_m.get("iptm", float("nan"))
            sel  = _safe(iptm_on) - (_safe(iptm_off) if not np.isnan(_safe(iptm_off, np.nan)) else 0.0)
            cd_on, cd_off = on_m.get("contact_density"), off_m.get("contact_density")
            on_metrics = {
                "af3_plddt": plddt_on, "af3_aptm": on_m.get("aptm"), "af3_iptm": iptm_on,
                "af3_iface_n_iface_res": on_m.get("n_iface_res"),
                "af3_iface_iface_plddt": on_m.get("iface_plddt"),
                "af3_iface_contact_density": cd_on,
            }
            comp = calc_specificity_composite(on_metrics,
                                              {"af3_iface_contact_density": cd_off}, model="af3")
            sel_cd = (cd_on - cd_off) if (cs._isnum(cd_on) and cs._isnum(cd_off)) else float("nan")
            loops = extract_loops(bseq, self.selected_loops)
            entry = {
                "design_id": f"AF3_{pk}_{abs(hash(bseq)) % 10**8}",
                "target_name": self.pairs[pk]["on_target"],
                "full_seq": bseq,
                "esm_iptm_on": iptm_on, "esm_iptm_off": iptm_off, "esm_plddt_on": plddt_on,
                "af3_aptm": on_m.get("aptm"),
                "af3_iface_contact_density": cd_on, "af3_iface_n_iface_res": on_m.get("n_iface_res"),
                "selectivity_score": sel, "selectivity_cd": sel_cd,
                "specificity_pair": pk,
                "on_target": self.pairs[pk]["on_target"], "off_target": self.pairs[pk]["off_target"],
                "composite_score": comp, "source": "AF3",
                "promising": (comp >= 0.5 and _safe(sel_cd, 0.0) > 0),
                **loops,
            }
            # Replace any existing specificity-HOF entry with the same sequence
            hof = self.state["specificity_hof"].setdefault(pk, [])
            hof = [e for e in hof if e.get("full_seq") != bseq]
            # Preserve the original design_id if we had this sequence before
            for e in self.state["specificity_hof"].get(pk, []):
                if e.get("full_seq") == bseq:
                    entry["design_id"] = e.get("design_id", entry["design_id"])
                    break
            hof.append(entry)
            self.state["specificity_hof"][pk] = hof
            n += 1
            logger.info(f"  [{pk}] AF3: ipTM_on={iptm_on:.3f} ipTM_off={iptm_off:.3f} "
                        f"selectivity={sel:.3f}")

        for pk in self.pair_keys:
            self.state["specificity_hof"][pk] = self._consolidate_hof(
                self.state["specificity_hof"].get(pk, []), HOF_SIZE
            )
        self._save_state()
        logger.info(f"AF3 specificity import complete: {n} designs updated.")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Specificity-focused TIMP3 binder design.")
    parser.add_argument("--pair", nargs="+", default=["MMP", "ADAM"],
                        choices=list(SPECIFICITY_PAIRS.keys()), help="Specificity pairs.")
    parser.add_argument("--loops", nargs="+", default=["AB", "C", "EF"],
                        choices=list(LOOP_CONFIGS.keys()), help="Loop regions to redesign.")
    parser.add_argument("--max-iterations", type=int, default=None, help="Stop after N iterations.")
    parser.add_argument("--import-af3", type=str, default=None,
                        help="AF3 results to import (.zip from the server, or legacy .json).")
    parser.add_argument("--esmfold2-gpus", default=None, metavar="N|auto",
                        help="Data-parallel ESMFold2 across GPUs (default 1; N free GPUs, or 'auto').")
    # Anneal / throughput controls — mirror iterative_refinement so the specificity
    # campaign can run the SAME hot->cold schedule. These set the inherited
    # iterative_refinement module globals (the parent's methods read them at call time).
    parser.add_argument("--backbones-per-target", type=int, default=None,
                        help=f"RFd3 backbones per target per iteration (default {ir.BACKBONES_PER_TARGET}).")
    parser.add_argument("--seqs-per-backbone", type=int, default=None,
                        help=f"LMPNN sequences per backbone (default {ir.LMPNN_SEQS_PER_BACKBONE}).")
    parser.add_argument("--init-temperature", type=float, default=None,
                        help=f"Starting LMPNN temperature for a fresh run (default {ir.INIT_TEMPERATURE}).")
    parser.add_argument("--min-temperature", type=float, default=None,
                        help=f"Temperature floor (default {ir.MIN_TEMPERATURE}).")
    parser.add_argument("--temp-decay", type=float, default=None,
                        help=f"Per-iteration temperature multiplier (default {ir.TEMP_DECAY}).")
    parser.add_argument("--no-adaptive-bias", action="store_true",
                        help="Disable adaptive loop-length narrowing (keeps full diversity).")
    args = parser.parse_args()

    if args.esmfold2_gpus is not None:
        ir.ESMFOLD2_GPUS = args.esmfold2_gpus   # resolved at call time by the parent
    if args.backbones_per_target is not None:
        ir.BACKBONES_PER_TARGET = args.backbones_per_target
    if args.seqs_per_backbone is not None:
        ir.LMPNN_SEQS_PER_BACKBONE = args.seqs_per_backbone
    if args.init_temperature is not None:
        ir.INIT_TEMPERATURE = args.init_temperature
    if args.min_temperature is not None:
        ir.MIN_TEMPERATURE = args.min_temperature
    if args.temp_decay is not None:
        ir.TEMP_DECAY = args.temp_decay
    if args.no_adaptive_bias:
        ir.ADAPTIVE_BIAS_START = 10**9

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    refiner = SpecificityRefiner(pair_keys=args.pair, active_loops=args.loops)

    if args.import_af3:
        if args.import_af3.lower().endswith(".zip"):
            refiner.import_af3_zip(args.import_af3)
        else:
            refiner.import_af3_results(args.import_af3)

    refiner.main_loop(max_iterations=args.max_iterations)


if __name__ == "__main__":
    main()
