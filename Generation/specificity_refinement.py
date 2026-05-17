"""
specificity_refinement.py

Specificity-focused iterative binder design for TIMP3-scaffold loops.

Extends the iterative_refinement pipeline by scoring every designed binder
against BOTH the on-target and off-target protease in a specificity pair:
  - MMP2 (on-target) vs MMP9 (off-target)
  - ADAM10 (on-target) vs ADAM17 (off-target)

Additional metrics computed per design:
  selectivity_score  = ipTM_on - ipTM_off
      > 0 : prefers on-target
      < 0 : prefers off-target (discard)
  selectivity_ratio  = ipTM_on / max(ipTM_off, 1e-3)

The Hall of Fame is ranked by a composite that balances absolute binding
affinity (composite_score) with selectivity (selectivity_score).

Usage:
  python specificity_refinement.py --pair MMP --loops AB C EF
  python specificity_refinement.py --pair ADAM --loops AB C EF
  python specificity_refinement.py --pair MMP ADAM --loops AB C EF

Output: Local/specificity_refinement/
"""

import os
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from iterative_refinement import (
    IterativeRefiner,
    TARGETS,
    LOOP_CONFIGS,
    DATA_DIR,
    OUT_BASE as _IR_OUT_BASE,
    BACKBONES_PER_TARGET,
    HOF_SIZE_PER_TARGET,
    AF3_EXPORT_EVERY_N,
    AF3_TOP_N,
    IPTM_PROMISING,
    MIN_TEMPERATURE,
    TEMP_DECAY,
    INIT_TEMPERATURE,
    setup_env,
    get_seq,
    renumber,
    pdb_chain_seq,
    extract_loops,
    get_fixed_residues,
    calc_composite,
    backbone_rmsd,
    count_interface_contacts,
)
from atomworks.io.utils.io_utils import to_cif_file
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
_HERE    = Path(__file__).parent.resolve()
OUT_BASE = _HERE / ".." / "Local" / "specificity_refinement"

# ── Specificity pair definitions ──────────────────────────────────────────────
# on_target: the protease we want to bind
# off_target: the protease to avoid
SPECIFICITY_PAIRS = {
    "MMP":  {"on_target": "MMP2",  "off_target": "MMP9"},
    "ADAM": {"on_target": "ADAM10", "off_target": "ADAM17"},
}

# Composite score weights for the combined specificity HOF ranking.
# Selectivity is given a dedicated weight alongside the standard affinity metrics.
SPECIFICITY_COMPOSITE_WEIGHTS = {
    "iptm_on":      0.35,  # Binding to on-target
    "selectivity":  0.25,  # ipTM_on - ipTM_off (selectivity)
    "plddt":        0.15,  # Fold confidence
    "ptm":          0.10,  # Global structure quality
    "rmsd":         0.10,  # Backbone fidelity
    "pae_on":       0.05,  # Interface PAE on-target
}
RMSD_CLIP  = 5.0
PAE_MAX    = 30.0
HOF_SIZE   = 75           # Per specificity pair


def calc_specificity_composite(iptm_on: float, iptm_off: float, plddt: float,
                                ptm: float, rmsd_val: float,
                                iface_pae_on: float) -> float:
    """Composite score that rewards both on-target binding and selectivity."""
    w   = SPECIFICITY_COMPOSITE_WEIGHTS
    sel = float(iptm_on) - float(iptm_off)          # in [-1, 1]
    sel_norm = (sel + 1.0) / 2.0                    # rescale to [0, 1]
    rmsd_s   = max(0.0, 1.0 - rmsd_val / RMSD_CLIP)
    pae_s    = max(0.0, 1.0 - iface_pae_on / PAE_MAX) if not np.isnan(iface_pae_on) else 0.0
    return (
        w["iptm_on"]      * float(iptm_on)
        + w["selectivity"]  * sel_norm
        + w["plddt"]        * (float(plddt) / 100.0)
        + w["ptm"]          * float(ptm)
        + w["rmsd"]         * rmsd_s
        + w["pae_on"]       * pae_s
    )


class SpecificityRefiner(IterativeRefiner):
    """
    Extends IterativeRefiner to score each design against both the on-target
    and off-target protease and compute a selectivity score.

    Parameters
    ----------
    pair_keys : list[str]
        Keys from SPECIFICITY_PAIRS to run (e.g. ["MMP", "ADAM"]).
    active_loops : list[str]
        Subset of LOOP_CONFIGS keys.
    """

    def __init__(self, pair_keys: list, active_loops: list):
        # Collect all targets needed across pairs
        all_targets = list({
            t
            for pk in pair_keys
            for t in (SPECIFICITY_PAIRS[pk]["on_target"], SPECIFICITY_PAIRS[pk]["off_target"])
        })
        self.pair_keys  = pair_keys
        self.pairs      = {pk: SPECIFICITY_PAIRS[pk] for pk in pair_keys}

        super().__init__(
            active_targets=all_targets,
            active_loops=active_loops,
            state_path=OUT_BASE / "specificity_state.json",
        )

        # Specificity HOF: per pair
        if "specificity_hof" not in self.state:
            self.state["specificity_hof"] = {pk: [] for pk in pair_keys}

        # Ensure off-target sequences are loaded
        for pk, pair in self.pairs.items():
            off = pair["off_target"]
            if off not in self.target_seqs:
                pdb_path = DATA_DIR / TARGETS[off]["pdb"]
                if pdb_path.exists():
                    self.target_seqs[off] = pdb_chain_seq(
                        str(pdb_path), TARGETS[off]["target_chain"]
                    )

    # ── RF3 on a specific target (not necessarily the design target) ──────────

    def _score_against_target(self, binder_seq: str, target_name: str,
                               rf3_engine, design_id: str,
                               out_dir: Path, suffix: str = "") -> dict:
        """
        Run RF3 complex prediction for binder_seq vs target_name.
        Returns dict with iptm, plddt, ptm, interface_pae.
        """
        target_seq = self.target_seqs.get(target_name, "")
        if not target_seq:
            return {"iptm": 0.0, "plddt": 0.0, "ptm": 0.0, "interface_pae": float("nan")}

        tcfg = TARGETS[target_name]
        bc   = tcfg["binder_chain"]
        fc   = tcfg["target_chain"]

        rf3_data = {
            "name": f"{design_id}_{suffix or target_name}",
            "components": [
                {"id": bc, "sequence": binder_seq},
                {"id": fc, "sequence": target_seq},
            ],
        }
        rf3_input = InferenceInput.from_json_dict(rf3_data)
        rf3_outs  = rf3_engine.run(inputs=rf3_input, annotate_b_factor_with_plddt=True)
        rf3_key   = next(iter(rf3_outs))
        rf3_out   = rf3_outs[rf3_key][0]
        rf3_arr   = renumber(rf3_out.atom_array)

        conf  = rf3_out.summary_confidences or {}
        plddt = conf.get("overall_plddt", conf.get("plddt", 0.0))
        ptm   = conf.get("ptm", 0.0)
        iptm  = (conf.get("iptm") or conf.get("ipTM") or
                 conf.get("interface_ptm") or conf.get("complex_pTM") or 0.0)

        iface_pae = float("nan")
        if hasattr(rf3_out, "confidences") and rf3_out.confidences:
            raw_pae = rf3_out.confidences.get("pae")
            if raw_pae is not None:
                pae_mat = np.array(raw_pae)
                nb = len(binder_seq)
                if pae_mat.shape[0] > nb:
                    iface_pae = float(
                        (np.mean(pae_mat[:nb, nb:]) + np.mean(pae_mat[nb:, :nb])) / 2
                    )

        # Save CIF for inspection
        cif_out = out_dir / f"{design_id}_{suffix or target_name}_rf3.cif"
        to_cif_file(rf3_arr, str(cif_out), file_type="cif")

        return {
            "iptm":          float(iptm),
            "plddt":         float(plddt),
            "ptm":           float(ptm),
            "interface_pae": iface_pae,
            "cif":           str(cif_out),
        }

    # ── Override RF3 scoring to include off-target cross-scoring ─────────────

    def run_rf3_complex(self, target_name: str, candidates: list,
                        out_dir: Path) -> list:
        """
        For each candidate, score against its designed target AND any off-targets
        in which it participates as an on-target (via SPECIFICITY_PAIRS).

        Falls back to parent implementation for targets not in any specificity pair.
        """
        # Determine if this target is an on-target for any pair
        pair_for_target = {
            v["on_target"]: pk for pk, v in self.pairs.items()
        }

        if target_name not in pair_for_target:
            # Not an on-target; use standard scoring
            return super().run_rf3_complex(target_name, candidates, out_dir)

        pk        = pair_for_target[target_name]
        on_tname  = self.pairs[pk]["on_target"]
        off_tname = self.pairs[pk]["off_target"]

        if not candidates:
            return []

        out_dir.mkdir(parents=True, exist_ok=True)
        engine = RF3InferenceEngine(ckpt_path="rf3", verbose=False)
        scored = []
        t0     = time.time()

        for cand in candidates:
            did  = cand["design_id"]
            bseq = cand.get("full_seq", "")
            if not bseq:
                continue
            try:
                # Score against on-target
                on_metrics = self._score_against_target(
                    bseq, on_tname, engine, did, out_dir, suffix="on"
                )
                # Score against off-target
                off_metrics = self._score_against_target(
                    bseq, off_tname, engine, did, out_dir, suffix="off"
                )

                sel_score = on_metrics["iptm"] - off_metrics["iptm"]
                sel_ratio = on_metrics["iptm"] / max(off_metrics["iptm"], 1e-3)

                # RMSD: RFd3 backbone vs on-target RF3 prediction
                # Load the on-target CIF to get the atom array for RMSD
                rmsd_val = float("nan")
                try:
                    from biotite.structure.io.pdbx import CIFFile as _CIF, get_structure as _gs
                    on_cif_path = on_metrics.get("cif", "")
                    if on_cif_path and Path(on_cif_path).exists():
                        on_arr = renumber(_gs(_CIF.read(on_cif_path), model=1))
                        rmsd_val = backbone_rmsd(cand["rfd3_array"], on_arr)
                except Exception:
                    pass

                comp = calc_specificity_composite(
                    iptm_on=on_metrics["iptm"],
                    iptm_off=off_metrics["iptm"],
                    plddt=on_metrics["plddt"],
                    ptm=on_metrics["ptm"],
                    rmsd_val=cand.get("rmsd_to_rfd3", 2.5),
                    iface_pae_on=on_metrics["interface_pae"],
                )

                entry = {
                    **{k: v for k, v in cand.items() if k not in ("array", "rfd3_array")},
                    # On-target metrics
                    "iptm_on":           on_metrics["iptm"],
                    "plddt_on":          on_metrics["plddt"],
                    "ptm_on":            on_metrics["ptm"],
                    "interface_pae_on":  on_metrics["interface_pae"],
                    # Off-target metrics
                    "iptm_off":          off_metrics["iptm"],
                    "plddt_off":         off_metrics["plddt"],
                    "ptm_off":           off_metrics["ptm"],
                    "interface_pae_off": off_metrics["interface_pae"],
                    # Selectivity
                    "selectivity_score": sel_score,
                    "selectivity_ratio": sel_ratio,
                    # Also alias iptm/plddt/ptm to the on-target values for HOF compat
                    "iptm":              on_metrics["iptm"],
                    "plddt":             on_metrics["plddt"],
                    "ptm":               on_metrics["ptm"],
                    "interface_pae":     on_metrics["interface_pae"],
                    "rmsd_to_rfd3":      rmsd_val,
                    "composite_score":   comp,
                    "specificity_pair":  pk,
                    "on_target":         on_tname,
                    "off_target":        off_tname,
                    "promising":         (on_metrics["iptm"] >= IPTM_PROMISING and sel_score > 0),
                    "rf3_cif":           on_metrics["cif"],
                    "iteration":         self.state["iteration"],
                    "temperature":       self.state["temperature"],
                }
                scored.append(entry)

                # Update specificity HOF
                self.state["specificity_hof"].setdefault(pk, []).append(entry)

            except Exception as exc:
                logger.error(f"Specificity RF3 error on {did}: {exc}")

        # Trim specificity HOF
        for pk in self.pair_keys:
            self.state["specificity_hof"][pk] = sorted(
                self.state["specificity_hof"].get(pk, []),
                key=lambda x: x.get("composite_score", 0),
                reverse=True,
            )[:HOF_SIZE]

        logger.info(
            f"[{target_name}→{off_tname}] RF3 specificity scoring done in "
            f"{(time.time()-t0)/60:.1f} min ({len(scored)}/{len(candidates)} scored)"
        )
        del engine
        torch.cuda.empty_cache()
        return scored

    # ── Specificity HOF summary ───────────────────────────────────────────────

    def _write_hof_summary(self) -> None:
        """Write both the standard HOF summary and a selectivity-focused table."""
        super()._write_hof_summary()

        # Specificity HOF summary
        rows = []
        for pk in self.pair_keys:
            for rank, entry in enumerate(self.state["specificity_hof"].get(pk, []), start=1):
                rows.append({"pair": pk, "rank": rank,
                             **{k: v for k, v in entry.items()
                                if k not in ("array", "rfd3_array")}})
        if rows:
            df = pd.DataFrame(rows)
            df.to_csv(OUT_BASE / "specificity_hof_summary.csv", index=False)
            # Log top entry per pair
            for pk in self.pair_keys:
                sub = df[df["pair"] == pk]
                if not sub.empty:
                    best = sub.iloc[0]
                    logger.info(
                        f"Specificity HOF best [{pk}]: {best.get('design_id')} | "
                        f"ipTM_on={best.get('iptm_on', 0):.3f} | "
                        f"ipTM_off={best.get('iptm_off', 0):.3f} | "
                        f"selectivity={best.get('selectivity_score', 0):.3f}"
                    )

    # ── AF3 export (specificity-aware) ────────────────────────────────────────

    def export_for_af3(self, force: bool = False) -> None:
        """
        Export AF3 submissions with both on-target AND off-target chains,
        so that AF3 Server can validate the selectivity predictions.
        """
        it = self.state["iteration"]
        if not force and (it - self.state.get("last_af3_it", -1)) < AF3_EXPORT_EVERY_N:
            return

        jobs = []
        seen_seqs: set = set()

        for pk in self.pair_keys:
            pair  = self.pairs[pk]
            on_t  = pair["on_target"]
            off_t = pair["off_target"]
            on_seq  = self.target_seqs.get(on_t, "")
            off_seq = self.target_seqs.get(off_t, "")

            hof_entries = self.state["specificity_hof"].get(pk, [])[:AF3_TOP_N // len(self.pair_keys)]
            for i, e in enumerate(hof_entries):
                bseq = e.get("full_seq", "")
                if not bseq or bseq in seen_seqs:
                    continue
                seen_seqs.add(bseq)
                # On-target job
                jobs.append({
                    "name":       f"spec_it{it}_{pk}_on_{i:02d}",
                    "modelSeeds": [42],
                    "sequences":  [
                        {"proteinChain": {"sequence": bseq,   "count": 1}},
                        {"proteinChain": {"sequence": on_seq,  "count": 1}},
                    ],
                })
                # Off-target job (same binder, different protease — for selectivity check)
                jobs.append({
                    "name":       f"spec_it{it}_{pk}_off_{i:02d}",
                    "modelSeeds": [42],
                    "sequences":  [
                        {"proteinChain": {"sequence": bseq,    "count": 1}},
                        {"proteinChain": {"sequence": off_seq,  "count": 1}},
                    ],
                })

        if not jobs:
            return

        export_path = OUT_BASE / f"af3_specificity_it{it}.json"
        with open(export_path, "w") as f:
            json.dump(jobs, f, indent=2)

        self.state["last_af3_it"] = it
        self._save_state()

        logger.info(f"AF3 specificity export: {len(jobs)} jobs → {export_path}")
        print("\n" + "=" * 60)
        print(f"  AF3 specificity submission: {export_path}")
        print(f"  Contains on-target AND off-target jobs for selectivity validation.")
        print(f"  Import results via --import-af3 <results.json>")
        print("=" * 60 + "\n")

    # ── Main loop ─────────────────────────────────────────────────────────────

    def main_loop(self, max_iterations: int = None) -> None:
        """Run the specificity-aware iterative loop."""
        logger.info(
            f"Starting specificity refinement.\n"
            f"  Pairs:   {self.pair_keys}\n"
            f"  Loops:   {self.active_loops}\n"
            f"  Output:  {OUT_BASE}\n"
        )
        it = 0
        while max_iterations is None or it < max_iterations:
            self.run_iteration()
            self.export_for_af3()
            self.state["temperature"] = max(
                MIN_TEMPERATURE, self.state["temperature"] * TEMP_DECAY
            )
            self.state["iteration"] += 1
            self._save_state()
            it += 1


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Specificity-focused iterative TIMP3 binder design."
    )
    parser.add_argument(
        "--pair", nargs="+", default=["MMP", "ADAM"],
        choices=list(SPECIFICITY_PAIRS.keys()),
        help="Specificity pairs to optimize.",
    )
    parser.add_argument(
        "--loops", nargs="+", default=["AB", "C", "EF"],
        choices=list(LOOP_CONFIGS.keys()),
        help="Loop regions to redesign.",
    )
    parser.add_argument(
        "--max-iterations", type=int, default=None,
        help="Stop after N iterations.",
    )
    parser.add_argument(
        "--import-af3", type=str, default=None,
        help="Path to AF3 results JSON to import.",
    )
    args = parser.parse_args()

    OUT_BASE.mkdir(parents=True, exist_ok=True)

    refiner = SpecificityRefiner(
        pair_keys=args.pair,
        active_loops=args.loops,
    )

    if args.import_af3:
        refiner.import_af3_results(args.import_af3)

    refiner.main_loop(max_iterations=args.max_iterations)


if __name__ == "__main__":
    main()
