"""
loop_probe.py

Probe the generative binder pipeline's *sequence preferences* per loop position.

For a chosen target and set of TIMP3 loops, this runs the same two generative
stages used in production —

    RFd3 (build a loop backbone of a FIXED user-specified length)
      -> LigandMPNN (design the loop sequence on that fixed scaffold)

— many times, then tallies which amino acids (and which biochemical groups)
LigandMPNN places at each loop position.  There is **no structure-validation
stage** (no AF3 / RF3 / ESMFold2): we only care about the emitted sequences, so
the funnel's scoring machinery is intentionally skipped.  The output is a set of
per-position frequency heatmaps (raw 20-AA + charge / size / type /
hydrophobicity / polarity / aromaticity) written per loop.

Starting structures are the AlphaFold complexes in
    Data/TIMP_Complexes/AlphaFold_CIF/TIMP3_vs_<TARGET>_AF.cif
where chain A is full-length TIMP3 (188 aa) and chain B is the target.  We trim
TIMP3 to the N-terminal design construct (residues 1..scaffold_len, default 121)
and relabel to the pipeline's canonical binder=A / target=B convention.

This module reuses the engine wrappers and helpers from `iterative_refinement`
(so it must run in the same GPU `foundry` conda env with FOUNDRY_CHECKPOINT_DIRS
set).  All heatmap/counting logic lives in the GPU-free `loop_probe_analysis`
module, so figures can be rebuilt on any machine from the CSVs written here.

Example
-------
    conda activate foundry
    python Generation/Loop-Probe/loop_probe.py --target MMP2 --loops AB C EF \
        --n-backbones 40 --seqs-per-backbone 3 --temperature 0.3

Native loop lengths (default) come from LOOP_CONFIGS: AB=6, C=6, EF=4, GH=10.
"""

from __future__ import annotations

import sys
import json
import time
import logging
import argparse
import dataclasses
from pathlib import Path

import numpy as np
import pandas as pd

# ── Reuse the production pipeline's engines + helpers ──────────────────────────
# iterative_refinement performs GPU-aware env setup at import time and pulls in
# rfd3 / mpnn / atomworks; importing it here means this script only runs in the
# foundry env (by design).
_HERE = Path(__file__).parent.resolve()
# This module lives in Generation/Loop-Probe/; add its own dir (for the sibling
# loop_probe_analysis) and the parent Generation/ (for iterative_refinement).
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent))

import torch  # noqa: E402
import biotite.structure.io.pdbx as pdbx  # noqa: E402
from biotite.structure.io.pdb import PDBFile  # noqa: E402

from iterative_refinement import (  # noqa: E402
    LOOP_CONFIGS,
    DESIGN_BINDER_CHAIN,
    DESIGN_TARGET_CHAIN,
    setup_env,
    renumber,
    get_seq,
    get_fixed_residues,
    extract_loops,
)
from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine  # noqa: E402
from rfd3.inference.input_parsing import DesignInputSpecification  # noqa: E402
from mpnn.inference_engines.mpnn import MPNNInferenceEngine  # noqa: E402

import loop_probe_analysis as lpa  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger("loop_probe")
for _noisy in ("transforms", "atomworks.io", "atomworks.ml", "foundry", "lightning"):
    logging.getLogger(_noisy).setLevel(logging.ERROR)

# ── Paths / target inputs ──────────────────────────────────────────────────────
# _HERE is Generation/Loop-Probe/, so the repo root is two levels up.
_REPO    = _HERE.parent.parent
AF_DIR   = _REPO / "Data" / "TIMP_Complexes" / "AlphaFold_CIF"
OUT_BASE = _REPO / "Local" / "loop_probe"

# Target -> AlphaFold complex CIF (chain A = TIMP3 binder, chain B = target).
# These six clean "TIMP3_vs_<T>_AF.cif" folds cover the pipeline's target set.
PROBE_TARGETS = {
    "MMP2":   "TIMP3_vs_MMP2_AF.cif",
    "MMP3":   "TIMP3_vs_MMP3_AF.cif",
    "MMP9":   "TIMP3_vs_MMP9_AF.cif",
    "MMP10":  "TIMP3_vs_MMP10_AF.cif",
    "ADAM10": "TIMP3_vs_ADAM10_AF.cif",
    "ADAM17": "TIMP3_vs_ADAM17_AF.cif",
}

# TIMP3 chain in the AF CIF is full-length (188 aa); the design construct is the
# N-terminal domain.  1..121 is the pipeline's scaffold_len; extended only if a
# C-terminal loop (GH) is requested.
DEFAULT_SCAFFOLD_LEN = 121
FULL_TIMP3_LEN       = 188
AF_BINDER_CHAIN      = "A"   # TIMP3 in the AlphaFold CIF
AF_TARGET_CHAIN      = "B"   # protease in the AlphaFold CIF


# ── Loop-length resolution ──────────────────────────────────────────────────────
def resolve_lengths(active_loops: list[str], overrides: dict[str, int] | None) -> dict[str, int]:
    """Native (LOOP_CONFIGS['normal']) length per loop, overridden where given."""
    overrides = overrides or {}
    lengths = {}
    for name in active_loops:
        if name not in LOOP_CONFIGS:
            raise ValueError(f"unknown loop {name!r}; known: {list(LOOP_CONFIGS)}")
        L = int(overrides.get(name, LOOP_CONFIGS[name]["normal"]))
        if L < 1:
            raise ValueError(f"loop {name} length must be >= 1 (got {L})")
        lengths[name] = L
    return lengths


def required_scaffold_len(active_loops: list[str]) -> int:
    """Smallest N-TIMP3 length that still contains every selected loop + flanks."""
    need = DEFAULT_SCAFFOLD_LEN
    for name in active_loops:
        lc = LOOP_CONFIGS[name]
        need = max(need, lc["pos"] + lc["normal"] + len(lc["right"]) + 1)
    return min(need, FULL_TIMP3_LEN)


# ── Input structure preparation (AF CIF -> trimmed design PDB) ──────────────────
def prepare_input_pdb(target: str, af_dir: Path, out_pdb: Path,
                      scaffold_len: int) -> tuple[Path, int]:
    """
    Read the AlphaFold complex CIF, keep TIMP3 chain A residues 1..scaffold_len
    (the design construct) + the whole target chain B, relabel to the pipeline's
    binder=A / target=B convention, renumber each chain from 1, and write a PDB
    RFd3 can consume.  Cached: skipped if `out_pdb` already exists.

    Returns (pdb_path, target_len).
    """
    cif_path = af_dir / PROBE_TARGETS[target]
    if not cif_path.exists():
        raise FileNotFoundError(f"AF CIF not found for {target}: {cif_path}")

    cif = pdbx.CIFFile.read(str(cif_path))
    arr = pdbx.get_structure(cif, model=1)

    binder = arr[(arr.chain_id == AF_BINDER_CHAIN) & (arr.res_id <= scaffold_len)]
    target_arr = arr[arr.chain_id == AF_TARGET_CHAIN]
    if len(binder) == 0 or len(target_arr) == 0:
        raise RuntimeError(
            f"{target}: expected chains {AF_BINDER_CHAIN}/{AF_TARGET_CHAIN} in "
            f"{cif_path.name}; got chains {sorted(set(arr.chain_id))}")

    binder.chain_id[:] = DESIGN_BINDER_CHAIN     # "A"
    target_arr.chain_id[:] = DESIGN_TARGET_CHAIN  # "B"
    combined = renumber(binder + target_arr)

    target_len = len(np.unique(combined.res_id[combined.chain_id == DESIGN_TARGET_CHAIN]))
    binder_len = len(np.unique(combined.res_id[combined.chain_id == DESIGN_BINDER_CHAIN]))
    logger.info(f"[{target}] prepared construct: binder(N-TIMP3)={binder_len} aa, "
                f"target={target_len} aa -> {out_pdb.name}")

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    pf = PDBFile()
    pf.set_structure(combined)
    pf.write(str(out_pdb))
    return out_pdb, target_len


# ── Fixed-length contig (loop lengths pinned to lengths[name]) ──────────────────
def build_fixed_contig(selected_loops: list[dict], scaffold_len: int,
                       target_len: int, lengths: dict[str, int]) -> tuple[str, str]:
    """
    RFd3 contig with each designed loop pinned to a single length.

    Mirrors iterative_refinement._build_contig but with lo==hi==lengths[name].
    The fixed-segment cursor advances by the loop's *native* length (`normal`)
    because those are the residue indices in the INPUT construct, independent of
    the designed length.
    """
    bc = DESIGN_BINDER_CHAIN
    fc = DESIGN_TARGET_CHAIN
    parts, current = [], 1
    for lc in selected_loops:
        if current <= lc["pos"]:
            parts.append(f"{bc}{current}-{lc['pos']}")
        L = lengths[lc["name"]]
        parts.append(f"{L}-{L}")
        current = lc["pos"] + lc["normal"] + 1
    if current <= scaffold_len:
        parts.append(f"{bc}{current}-{scaffold_len}")

    full_contig = ",".join(parts) + f",/0,{fc}1-{target_len}"
    total = (scaffold_len - sum(lc["normal"] for lc in selected_loops)
             + sum(lengths.values()) + target_len)
    return full_contig, f"{total}-{total}"


# ── Generation ──────────────────────────────────────────────────────────────────
def run_rfd3_fixed(input_pdb: Path, contig: str, length_range: str,
                   n_designs: int) -> list:
    """Sample `n_designs` fixed-length backbones with RFd3. Returns AtomArrays."""
    cfg = RFD3InferenceConfig(
        diffusion_batch_size=min(10, n_designs),
        low_memory_mode=False,
        specification={"length": length_range, "contig": contig, "extra": {}},
    )
    engine = RFD3InferenceEngine(**dataclasses.asdict(cfg))
    spec = DesignInputSpecification(input=str(input_pdb), contig=contig,
                                    length=length_range, extra={})
    n_batches = (n_designs + cfg.diffusion_batch_size - 1) // cfg.diffusion_batch_size

    t0 = time.time()
    outputs = engine.run(inputs=spec, n_batches=n_batches, out_dir=None)
    logger.info(f"RFd3: {n_designs} backbones in {(time.time()-t0)/60:.1f} min")

    arrays, idx = [], 0
    if outputs:
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


def run_lmpnn_fixed(backbones: list, selected_loops: list[dict], out_dir: Path,
                    temperature: float, seqs_per_backbone: int,
                    tag: str, seed: int) -> list[dict]:
    """
    Design loop sequences for each backbone with LigandMPNN (everything outside
    the selected loops — scaffold + flanks + entire target — held fixed).
    Returns one record per emitted sequence with the extracted loop sub-seqs.
    """
    bc, fc = DESIGN_BINDER_CHAIN, DESIGN_TARGET_CHAIN
    out_dir.mkdir(parents=True, exist_ok=True)
    engine = MPNNInferenceEngine(
        model_type="ligand_mpnn", is_legacy_weights=True,
        write_structures=False, write_fasta=False, out_directory=str(out_dir),
    )
    results: list[dict] = []
    t0 = time.time()
    for idx, bb in enumerate(backbones):
        design_id = f"{tag}_d{idx}"
        try:
            fixed_res = get_fixed_residues(bb, selected_loops, bc, fc)
            mpnn_in = {
                "name": design_id,
                "batch_size": seqs_per_backbone,
                "remove_waters": True,
                "seed": seed + idx,
                "fixed_residues": fixed_res,
                "sampling_temp": temperature,
            }
            outs = engine.run(input_dicts=[mpnn_in], atom_arrays=[bb])
            for si, mp in enumerate(outs):
                valid = ~np.isnan(mp.atom_array.coord[:, 0])
                arr = renumber(mp.atom_array[valid])
                bseq = get_seq(arr, bc)
                loops = extract_loops(bseq, selected_loops)
                results.append({"design_id": f"{design_id}_s{si}",
                                "full_seq": bseq, **loops})
        except Exception as exc:  # one backbone failing must not kill the run
            logger.error(f"LMPNN error on {design_id}: {exc}")
    logger.info(f"LigandMPNN: {len(results)} sequences in {(time.time()-t0)/60:.1f} min")
    del engine
    torch.cuda.empty_cache()
    return results


# ── Top-level probe (one target, one length configuration) ─────────────────────
def run_probe(target: str, active_loops: list[str], lengths: dict[str, int],
              n_backbones: int, seqs_per_backbone: int, temperature: float,
              out_dir: Path, af_dir: Path = AF_DIR, scaffold_len: int | None = None,
              seed: int = 42, make_plots: bool = True) -> dict:
    """
    Run RFd3->LigandMPNN at fixed loop lengths for one target and write, per
    loop: sequences.csv, position-count/frequency CSVs, and heatmaps.
    Returns a summary dict (also written as summary.json).
    """
    if target not in PROBE_TARGETS:
        raise ValueError(f"unknown target {target!r}; known: {list(PROBE_TARGETS)}")

    selected_loops = sorted(
        [{**LOOP_CONFIGS[n], "name": n} for n in active_loops],
        key=lambda x: x["pos"])
    scaffold_len = scaffold_len or required_scaffold_len(active_loops)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # Prepared design constructs are cached in one shared location, keyed by
    # target + construct length, so a whole sweep reuses a single trim per target.
    input_pdb = OUT_BASE / "inputs" / f"{target}_Nterm{scaffold_len}.pdb"
    if not input_pdb.exists():
        _, target_len = prepare_input_pdb(target, af_dir, input_pdb, scaffold_len)
    else:
        arr0 = PDBFile.read(str(input_pdb)).get_structure()[0]
        target_len = len(np.unique(arr0.res_id[arr0.chain_id == DESIGN_TARGET_CHAIN]))

    contig, length_range = build_fixed_contig(selected_loops, scaffold_len,
                                              target_len, lengths)
    length_sig = "_".join(f"{n}{lengths[n]}" for n in [lc["name"] for lc in selected_loops])
    logger.info(f"[{target}] loops={active_loops} lengths={lengths} "
                f"T={temperature} backbones={n_backbones}x{seqs_per_backbone}")
    logger.info(f"[{target}] contig: {contig}")

    backbones = run_rfd3_fixed(input_pdb, contig, length_range, n_backbones)
    tag = f"{target}_{length_sig}_T{temperature}"
    designs = run_lmpnn_fixed(backbones, selected_loops, out_dir / "lmpnn",
                              temperature, seqs_per_backbone, tag, seed)

    # ── Persist sequences ──────────────────────────────────────────────────────
    seq_rows = [{"design_id": d["design_id"], "target": target,
                 **{k: v for k, v in d.items() if k.startswith("loop_")},
                 "full_seq": d["full_seq"]} for d in designs]
    seq_df = pd.DataFrame(seq_rows)
    seq_df.to_csv(out_dir / "sequences.csv", index=False)
    with open(out_dir / "loops.fasta", "w") as fh:
        for d in designs:
            for lc in selected_loops:
                s = d.get(f"loop_{lc['name']}_seq", "MISSING")
                fh.write(f">{d['design_id']}|{lc['name']}|L{lengths[lc['name']]}\n{s}\n")

    # ── Count + heatmap per loop ────────────────────────────────────────────────
    per_loop = {}
    for lc in selected_loops:
        name, L = lc["name"], lengths[lc["name"]]
        col = f"loop_{name}_seq"
        raw = [s for s in seq_df.get(col, pd.Series([], dtype=str)).astype(str).tolist()
               if s and s != "MISSING"]
        good = [s for s in raw if len(s) == L and set(s) <= lpa.AA_SET]
        counts = lpa.count_positions(good, length=L)
        counts.attrs["n_sequences"] = len(good)
        per_loop[name] = {"n_designs": len(raw), "n_usable": len(good),
                          "n_parse_fail": len(raw) - len(good), "length": L}
        if make_plots:
            lpa.render_all_heatmaps(counts, out_dir, stem=name,
                                    title_prefix=f"{target} | {name} L{L} | ")
        else:
            lpa.to_frequency(counts).to_csv(out_dir / f"position_freq_{name}.csv")
            counts.to_csv(out_dir / f"position_counts_{name}.csv")

    summary = {
        "target": target, "loops": active_loops, "lengths": lengths,
        "scaffold_len": scaffold_len, "temperature": temperature,
        "n_backbones": n_backbones, "seqs_per_backbone": seqs_per_backbone,
        "n_designs_total": len(designs), "contig": contig,
        "length_range": length_range, "seed": seed, "per_loop": per_loop,
    }
    with open(out_dir / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    logger.info(f"[{target}] done -> {out_dir}  "
                + " ".join(f"{k}:{v['n_usable']}/{v['n_designs']}" for k, v in per_loop.items()))
    return summary


# ── CLI ─────────────────────────────────────────────────────────────────────────
def _parse_lengths(spec: str | None) -> dict[str, int]:
    """Parse '--lengths AB=8,C=6' into {'AB':8,'C':6}."""
    if not spec:
        return {}
    out = {}
    for tok in spec.replace(" ", "").split(","):
        if not tok:
            continue
        name, _, val = tok.partition("=")
        out[name] = int(val)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", required=True, choices=list(PROBE_TARGETS),
                    help="which TIMP3-vs-target AF complex to probe")
    ap.add_argument("--loops", nargs="+", default=["AB", "C", "EF"],
                    choices=list(LOOP_CONFIGS), help="loops to design (default AB C EF)")
    ap.add_argument("--lengths", default=None,
                    help="fix loop lengths, e.g. 'AB=8,C=6'; unspecified loops use native")
    ap.add_argument("--n-backbones", type=int, default=40,
                    help="RFd3 backbones to sample (default 40)")
    ap.add_argument("--seqs-per-backbone", type=int, default=3,
                    help="LigandMPNN sequences per backbone (default 3)")
    ap.add_argument("--temperature", type=float, default=0.3,
                    help="LigandMPNN sampling temperature (default 0.3)")
    ap.add_argument("--scaffold-len", type=int, default=None,
                    help="N-TIMP3 construct length (default auto: 121, extended for GH)")
    ap.add_argument("--out-dir", default=None,
                    help="output dir (default Local/loop_probe/<target>_<lengthsig>)")
    ap.add_argument("--af-dir", default=str(AF_DIR), help="AlphaFold CIF directory")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--no-plots", action="store_true",
                    help="write CSVs only (rebuild figures later with loop_probe_analysis)")
    args = ap.parse_args()

    setup_env()
    lengths = resolve_lengths(args.loops, _parse_lengths(args.lengths))
    length_sig = "_".join(f"{n}{lengths[n]}" for n in args.loops)
    out_dir = Path(args.out_dir) if args.out_dir else OUT_BASE / f"{args.target}_{length_sig}"

    run_probe(
        target=args.target, active_loops=args.loops, lengths=lengths,
        n_backbones=args.n_backbones, seqs_per_backbone=args.seqs_per_backbone,
        temperature=args.temperature, out_dir=out_dir, af_dir=Path(args.af_dir),
        scaffold_len=args.scaffold_len, seed=args.seed, make_plots=not args.no_plots,
    )


if __name__ == "__main__":
    main()
