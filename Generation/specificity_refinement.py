"""
specificity_refinement.py

Specificity-focused iterative binder design for TIMP3-scaffold loops.

Extends iterative_refinement by scoring every designed binder against BOTH the
on-target and off-target protease in a specificity pair:
  - MMP2 (on-target)   vs MMP9 (off-target)
  - ADAM10 (on-target) vs ADAM17 (off-target)

Scoring uses ESMFold2 (the validated local ranker — see BOLTZ_FILTER_METHODS.md),
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
    _safe,
    _normalize_plddt,
)

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
OUT_BASE = _HERE / ".." / "Local" / "specificity_refinement"

# ── Specificity pair definitions ──────────────────────────────────────────────
SPECIFICITY_PAIRS = {
    "MMP":  {"on_target": "MMP2",  "off_target": "MMP9"},
    "ADAM": {"on_target": "ADAM10", "off_target": "ADAM17"},
}

# Specificity composite (ESMFold2 signals only).  Mirrors the main pipeline's
# 0.5·ipTM + 0.5·pLDDT philosophy but reallocates weight to selectivity — the
# whole point here.  Equal-ish on-target binding vs selectivity, no per-target
# special-casing (avoids overfitting / systematic bias):
#   binding (iptm_on 0.35 + plddt_on 0.25 = 0.60)   selectivity 0.40
SPECIFICITY_COMPOSITE_WEIGHTS = {
    "iptm_on":     0.35,   # on-target interface quality (ESMFold2 ipTM)
    "plddt_on":    0.25,   # binder foldability in the on-target complex
    "selectivity": 0.40,   # (ipTM_on - ipTM_off), rescaled to [0, 1]
}
HOF_SIZE = 75   # per specificity pair


def calc_specificity_composite(iptm_on, iptm_off, plddt_on) -> float:
    """Composite in [0, 1] rewarding on-target binding AND selectivity (ESMFold2)."""
    w   = SPECIFICITY_COMPOSITE_WEIGHTS
    io  = _safe(iptm_on)
    off = _safe(iptm_off, 0.0)
    if np.isnan(off):
        off = 0.0
    sel_norm = (io - off + 1.0) / 2.0                 # [-1,1] → [0,1]
    pl = _normalize_plddt(plddt_on) / 100.0
    return w["iptm_on"] * io + w["plddt_on"] * pl + w["selectivity"] * sel_norm


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

            c.update({
                "esm_iptm_on":   iptm_on,
                "esm_plddt_on":  plddt_on,
                "esm_ptm_on":    on_m.get("esm_ptm"),
                "esm_iptm_off":  iptm_off,
                "esm_plddt_off": off_m.get("esm_plddt") if off_m else float("nan"),
                "selectivity_score": sel,
                "selectivity_ratio": (_safe(iptm_on) / max(_safe(iptm_off), 1e-3)
                                      if not np.isnan(_safe(iptm_off, np.nan)) else float("inf")),
                # esm_iptm/esm_plddt aliases = on-target (so inherited code/CSVs read sensibly)
                "esm_iptm":   iptm_on,
                "esm_plddt":  plddt_on,
                "esm_cif":    on_m.get("esm_cif"),
                "specificity_pair": pk,
                "on_target": on_t, "off_target": off_t,
                "source": "ESMFold2",
                "promising": (_safe(iptm_on) >= IPTM_PROMISING and sel > 0),
            })
            c["composite_score"] = calc_specificity_composite(iptm_on, iptm_off, plddt_on)
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

                    group.setdefault((pk, binder_seq), {})[side] = {
                        "iptm": iptm, "ptm": ptm, "plddt": plddt,
                    }
                except Exception as exc:
                    logger.error(f"Error parsing {jrf}: {exc}")

        n = 0
        for (pk, bseq), sides in group.items():
            on_m = sides.get("on")
            if not on_m:
                continue   # need at least the on-target result
            iptm_on  = on_m["iptm"]
            plddt_on = on_m["plddt"]
            iptm_off = sides.get("off", {}).get("iptm", float("nan"))
            sel  = _safe(iptm_on) - (_safe(iptm_off) if not np.isnan(_safe(iptm_off, np.nan)) else 0.0)
            comp = calc_specificity_composite(iptm_on, iptm_off, plddt_on)
            loops = extract_loops(bseq, self.selected_loops)
            entry = {
                "design_id": f"AF3_{pk}_{abs(hash(bseq)) % 10**8}",
                "target_name": self.pairs[pk]["on_target"],
                "full_seq": bseq,
                "esm_iptm_on": iptm_on, "esm_iptm_off": iptm_off, "esm_plddt_on": plddt_on,
                "selectivity_score": sel,
                "specificity_pair": pk,
                "on_target": self.pairs[pk]["on_target"], "off_target": self.pairs[pk]["off_target"],
                "composite_score": comp, "source": "AF3",
                "promising": (_safe(iptm_on) >= IPTM_PROMISING and sel > 0),
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
    args = parser.parse_args()

    if args.esmfold2_gpus is not None:
        ir.ESMFOLD2_GPUS = args.esmfold2_gpus   # resolved at call time by the parent

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
