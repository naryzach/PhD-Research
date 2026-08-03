"""
loop_probe.py

Probe the generative binder pipeline's *sequence preferences* per loop position.

For a chosen target and set of TIMP3 loops, this runs the same two generative
stages used in production —

    RFd3 (build loop backbones of FIXED user-specified lengths)
      -> LigandMPNN (design the loop sequences on that fixed scaffold)

— many times, then tallies which amino acids (and which biochemical groups)
LigandMPNN places at each loop position.  There is **no structure-validation
stage** (no AF3 / RF3 / ESMFold2): we only care about the emitted sequences, so
the funnel's scoring machinery is intentionally skipped.  Output is a set of
per-position frequency heatmaps (raw 20-AA + charge / size / type /
hydrophobicity / polarity / aromaticity) written per loop.

All selected loops are designed **simultaneously** in one LigandMPNN pass, so
loop-loop interactions are captured.  (The `loop_probe_sweep.py` wrapper adds
joint *length* sweeping for loops grouped as e.g. `AB_C`.)

Templates
---------
Any complex structure with a TIMP3-like chain and a target chain works — CIF or
PDB.  **Chain order is auto-detected** by sequence (see `identify_chains`), so
files whose binder/target chains are swapped need no special handling:

    af3 (default)  Data/TIMP_Complexes/AF3_Templates/<T>_TIMP3_AF3.pdb    (binder = chain A, 188 aa)
    alphafold      Data/TIMP_Complexes/AlphaFold_CIF/TIMP3_vs_<T>_AF.cif  (binder = chain A, 188 aa)
    haddock        Data/TIMP_Complexes/HADDOCK_Outputs/<T>_TIMP3_HADDOCK.pdb (binder = chain B, 121 aa)

By default the **full-length TIMP3 chain** is kept (188 aa — the N-terminal
domain that holds the loops plus the C-terminal domain as fixed structural
context).  Pass `--scaffold-len 121` to trim to just the N-terminal design
construct.  TIMP3 is relabelled to the pipeline's canonical binder=A / target=B
convention.  Loop positions and native lengths are **derived from the template's
own sequence** by locating the flanking tripeptides, so alternative templates and
numbering offsets work without editing LOOP_CONFIGS.

This module reuses the engine wrappers and helpers from `iterative_refinement`
(so it must run in the same GPU `foundry` conda env).  All heatmap/counting
logic lives in the GPU-free `loop_probe_analysis` module.

Example
-------
    conda activate foundry
    python Generation/Loop-Probe/loop_probe.py --target MMP2 --loops AB C EF
    python Generation/Loop-Probe/loop_probe.py --config my_run.yaml --target MMP2

Native loop lengths (default) come from LOOP_CONFIGS: AB=6, C=6, EF=4, GH=10,
MTL=10.  (GH and MTL are C-terminal-domain loops, so they need the full-length
TIMP3 construct — the default; they are unavailable under `--scaffold-len 121`.)
"""

from __future__ import annotations

import re
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

# ── Paths ──────────────────────────────────────────────────────────────────────
# _HERE is Generation/Loop-Probe/, so the repo root is two levels up.
_REPO    = _HERE.parent.parent
OUT_BASE = _REPO / "Local" / "loop_probe"

TARGET_NAMES = ["MMP2", "MMP3", "MMP9", "MMP10", "ADAM10", "ADAM17"]

# Built-in template sets: {name: (directory, filename pattern)}.  Anything else
# can be pointed at with --template-dir / --template-map.
TEMPLATE_SETS = {
    "af3":       (_REPO / "Data" / "TIMP_Complexes" / "AF3_Templates",
                  "{target}_TIMP3_AF3.pdb"),
    "alphafold": (_REPO / "Data" / "TIMP_Complexes" / "AlphaFold_CIF",
                  "TIMP3_vs_{target}_AF.cif"),
    "haddock":   (_REPO / "Data" / "TIMP_Complexes" / "HADDOCK_Outputs",
                  "{target}_TIMP3_HADDOCK.pdb"),
}
DEFAULT_TEMPLATE_SET = "af3"          # full-length TIMP3 (188 aa) AF3 complexes

# Mature human TIMP3 (188 aa) — the reference used to decide which chain of a
# template is the binder.  Designed/variant TIMP3s still score far above any
# protease chain, so detection is robust to loop mutations.
TIMP3_REF = (
    "CTCSPSHPQDAFCNSDIVIRAKVVGKKLVKEGPFGTLVYTIKQMKMYRGFTKMPHVQYIHTEASESLCGLK"
    "LEVNKYQYLLTGRVYDGKMYTGLCNFVERWDQLTLSQRKGLNYRYHLGCNCKIKSCYYLPCFVTSKNECLW"
    "TDMLSNFGYPGYQSKHYACIRQKGGYCSWYRGWAPPDKSIINATDP"
)
# Every loop flank tripeptide — a TIMP3 chain contains most of them, a protease
# essentially none.  (A 121-aa N-TIMP3 construct hits 6/8: GH's flanks lie
# beyond residue 121.)
FLANK_MOTIFS = sorted({m for lc in LOOP_CONFIGS.values()
                       for m in (lc["left"], lc["right"])})

# Construct length = how many N-terminal TIMP3 residues to keep.  The DEFAULT is
# now full length: scaffold_len=None keeps the entire binder chain (188 aa for
# the AF3/AlphaFold templates).  Pass an int (e.g. 121) to trim to the
# N-terminal domain instead.  NTERM_DOMAIN_LEN is the sensible floor / the classic
# N-TIMP3 construct used by the rest of the pipeline.
NTERM_DOMAIN_LEN = 121
FULL_TIMP3_LEN   = 188
BINDER_SCORE_MIN     = 0.50  # below this we refuse to guess the binder chain
BINDER_SCORE_MARGIN  = 0.20  # binder must beat runner-up by at least this

# ── Run defaults (tuned for multi-day sweeps; all CLI/config overridable) ───────
DEFAULT_N_BACKBONES       = 100   # RFd3 backbones per configuration
DEFAULT_SEQS_PER_BACKBONE = 5     # LigandMPNN sequences per backbone
DEFAULT_TEMPERATURE       = 0.5   # matches production's hot end (INIT_TEMPERATURE)
DEFAULT_SEED              = 42


# ── Config file ────────────────────────────────────────────────────────────────
def load_config(path: str | Path | None) -> dict:
    """
    Load a JSON or YAML run config.  Every CLI flag has a same-named key (dashes
    or underscores); explicit CLI flags always win over the config file.
    Returns {} when `path` is None.
    """
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"config not found: {p}")
    text = p.read_text()
    if p.suffix.lower() in (".yaml", ".yml"):
        try:
            import yaml
        except ImportError as exc:
            raise SystemExit(
                f"{p.name} is YAML but PyYAML is not installed — "
                f"`pip install pyyaml` or use a .json config") from exc
        cfg = yaml.safe_load(text) or {}
    else:
        cfg = json.loads(text)
    if not isinstance(cfg, dict):
        raise ValueError(f"{p}: config must be a mapping at the top level")
    return {str(k).replace("-", "_"): v for k, v in cfg.items()}


def cfg_get(cli_value, cfg: dict, key: str, default):
    """CLI value (if given) > config file > built-in default."""
    if cli_value is not None:
        return cli_value
    return cfg.get(key, default)


# ── Loop selection / grouping ──────────────────────────────────────────────────
def parse_loop_tokens(tokens: list[str]) -> tuple[list[list[str]], list[str]]:
    """
    Parse loop tokens into (groups, design_loops).

    A token joined by '_' is a GROUP whose loops are swept jointly by
    loop_probe_sweep (e.g. 'AB_C' sweeps AB and C lengths together, capturing
    their interaction).  Within a single probe run every listed loop is designed
    simultaneously regardless of grouping.

        ['AB_C', 'EF'] -> ([['AB','C'], ['EF']], ['AB','C','EF'])
    """
    groups: list[list[str]] = []
    design: list[str] = []
    for tok in tokens:
        grp = [x for x in str(tok).split("_") if x]
        if not grp:
            continue
        for name in grp:
            if name not in LOOP_CONFIGS:
                raise ValueError(f"unknown loop {name!r}; known: {list(LOOP_CONFIGS)}")
            if name not in design:
                design.append(name)
        groups.append(grp)
    if not design:
        raise ValueError("no loops selected")
    return groups, design


def resolve_lengths(active_loops: list[str], overrides: dict[str, int] | None,
                    geometry: list[dict] | None = None) -> dict[str, int]:
    """
    Native length per loop, overridden where given.  Native comes from the
    template-derived geometry when available, else LOOP_CONFIGS['normal'].
    """
    overrides = {k: int(v) for k, v in (overrides or {}).items()}
    native = {g["name"]: g["normal"] for g in (geometry or [])}
    lengths = {}
    for name in active_loops:
        if name not in LOOP_CONFIGS:
            raise ValueError(f"unknown loop {name!r}; known: {list(LOOP_CONFIGS)}")
        L = overrides.get(name, native.get(name, LOOP_CONFIGS[name]["normal"]))
        if int(L) < 1:
            raise ValueError(f"loop {name} length must be >= 1 (got {L})")
        lengths[name] = int(L)
    return lengths


def required_scaffold_len(active_loops: list[str]) -> int:
    """Smallest TIMP3 construct length that still contains every selected loop."""
    need = NTERM_DOMAIN_LEN
    for name in active_loops:
        lc = LOOP_CONFIGS[name]
        need = max(need, lc["pos"] + lc["normal"] + len(lc["right"]) + 1)
    return min(need, FULL_TIMP3_LEN)


# ── Template resolution + chain identification ─────────────────────────────────
def resolve_template(target: str, template_set: str = DEFAULT_TEMPLATE_SET,
                     template_dir: str | Path | None = None,
                     template_map: dict | None = None) -> Path:
    """
    Find the template structure for `target`.

    Priority: explicit template_map[target] > <template_dir or set dir>/<pattern>
    > any .cif/.pdb in that directory whose filename contains the target name.
    """
    if template_map and target in template_map:
        p = Path(template_map[target])
        return p if p.is_absolute() else (_REPO / p)

    if template_set not in TEMPLATE_SETS:
        raise ValueError(f"unknown template set {template_set!r}; "
                         f"known: {list(TEMPLATE_SETS)}")
    set_dir, pattern = TEMPLATE_SETS[template_set]
    directory = Path(template_dir) if template_dir else set_dir
    if not directory.is_absolute():
        directory = _REPO / directory

    cand = directory / pattern.format(target=target)
    if cand.exists():
        return cand
    hits = [p for p in sorted(directory.glob("*"))
            if p.suffix.lower() in (".cif", ".pdb") and target.lower() in p.name.lower()]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        raise ValueError(
            f"{target}: {len(hits)} candidate templates in {directory} "
            f"({[p.name for p in hits]}); disambiguate with --template-map")
    raise FileNotFoundError(f"{target}: no template found in {directory} "
                            f"(tried {cand.name})")


def load_structure(path: str | Path):
    """Read a .cif or .pdb into a biotite AtomArray (first model)."""
    path = Path(path)
    if path.suffix.lower() in (".cif", ".mmcif"):
        return pdbx.get_structure(pdbx.CIFFile.read(str(path)), model=1)
    return PDBFile.read(str(path)).get_structure()[0]


def timp3_likeness(seq: str) -> tuple[float, int, float]:
    """
    Score how TIMP3-like a chain sequence is.
    Combines loop-flank motif hits (robust to loop redesign / numbering offsets)
    with N-terminal identity to mature TIMP3.  Returns (score, n_flanks, identity).
    """
    if not seq:
        return 0.0, 0, 0.0
    hits = sum(1 for m in FLANK_MOTIFS if m in seq)
    n = min(len(seq), len(TIMP3_REF))
    ident = sum(1 for a, b in zip(seq[:n], TIMP3_REF[:n]) if a == b) / n if n else 0.0
    return 0.7 * (hits / len(FLANK_MOTIFS)) + 0.3 * ident, hits, ident


def identify_chains(arr, binder_chain: str | None = None,
                    target_chain: str | None = None,
                    label: str = "") -> tuple[str, str]:
    """
    Decide which chain is the TIMP3 binder and which is the target.

    Explicit `binder_chain` / `target_chain` are honoured (and validated).
    Otherwise the binder is the highest TIMP3-likeness chain and the target is
    the largest remaining chain.  Raises if the call is ambiguous, so a bad
    template fails loudly instead of silently designing the wrong chain.
    """
    chains = list(dict.fromkeys(arr.chain_id[arr.atom_name == "CA"].tolist()))
    if not chains:
        raise RuntimeError(f"{label}: no CA atoms / chains found")

    if binder_chain and target_chain:
        for ch, role in ((binder_chain, "binder"), (target_chain, "target")):
            if ch not in chains:
                raise ValueError(f"{label}: {role} chain {ch!r} not in {chains}")
        return binder_chain, target_chain

    scored = []
    for ch in chains:
        s, hits, ident = timp3_likeness(get_seq(arr, ch))
        scored.append((s, hits, ident, ch, int((arr.chain_id == ch).sum())))
    scored.sort(key=lambda r: -r[0])
    best = scored[0]

    if binder_chain:
        bc = binder_chain
    else:
        if best[0] < BINDER_SCORE_MIN:
            raise RuntimeError(
                f"{label}: no chain looks like TIMP3 (best {best[3]!r} score "
                f"{best[0]:.2f}). Pass --binder-chain/--target-chain explicitly.")
        if len(scored) > 1 and (best[0] - scored[1][0]) < BINDER_SCORE_MARGIN:
            raise RuntimeError(
                f"{label}: binder chain ambiguous ({best[3]!r} {best[0]:.2f} vs "
                f"{scored[1][3]!r} {scored[1][0]:.2f}). Pass --binder-chain.")
        bc = best[3]

    if target_chain:
        tc = target_chain
    else:
        others = [r for r in scored if r[3] != bc]
        if not others:
            raise RuntimeError(f"{label}: only one chain ({bc}); need a target chain")
        # largest remaining chain by atom count
        tc = max(others, key=lambda r: r[4])[3]
        if len(others) > 1:
            logger.warning(f"{label}: {len(chains)} chains present; using {tc!r} "
                           f"as target (largest non-binder)")

    bs = next(r for r in scored if r[3] == bc)
    logger.info(f"{label}: binder=chain {bc} (TIMP3 score {bs[0]:.2f}, "
                f"flanks {bs[1]}/{len(FLANK_MOTIFS)}, ident {bs[2]:.2f}), "
                f"target=chain {tc}")
    return bc, tc


# ── Template-derived loop geometry ─────────────────────────────────────────────
def derive_loop_geometry(binder_seq: str, active_loops: list[str]) -> list[dict]:
    """
    Locate each selected loop in THIS template's binder sequence via its flanking
    tripeptides, returning loop dicts with template-accurate `pos` (1-indexed
    last fixed residue before the loop) and `normal` (native loop length).

    This makes the probe independent of LOOP_CONFIGS' hard-coded numbering, so
    alternative templates / constructs / numbering offsets work unchanged.
    """
    geometry, cursor = [], 0
    for name in sorted(active_loops, key=lambda n: LOOP_CONFIGS[n]["pos"]):
        lc = LOOP_CONFIGS[name]
        m = re.compile(f"{lc['left']}([A-Z]*?){lc['right']}").search(binder_seq[cursor:])
        if not m:
            raise RuntimeError(
                f"loop {name}: flanks {lc['left']}...{lc['right']} not found in the "
                f"construct (len {len(binder_seq)}). The template may be truncated "
                f"(try a larger --scaffold-len) or not TIMP3-like.")
        abs_start = cursor + m.start()
        pos = abs_start + len(lc["left"])      # 1-indexed last fixed residue
        normal = len(m.group(1))
        cursor = cursor + m.end() - len(lc["right"])
        if pos != lc["pos"] or normal != lc["normal"]:
            logger.info(f"loop {name}: template geometry pos={pos} normal={normal} "
                        f"(LOOP_CONFIGS says pos={lc['pos']} normal={lc['normal']})")
        geometry.append({**lc, "name": name, "pos": pos, "normal": normal})
    return sorted(geometry, key=lambda g: g["pos"])


# ── Input structure preparation (template -> design PDB) ───────────────────────
def prepare_input_pdb(target: str, template: Path, out_pdb: Path,
                      scaffold_len: int | None = None,
                      binder_chain: str | None = None,
                      target_chain: str | None = None) -> tuple[Path, int]:
    """
    Read a template complex (CIF or PDB), auto-detect which chain is the TIMP3
    binder, keep the binder (full length when `scaffold_len` is None, else its
    N-terminal `scaffold_len` residues) plus the whole target chain, relabel to
    binder=A / target=B, renumber each chain from 1, and write a PDB RFd3 can
    consume.  Returns (pdb_path, target_len).
    """
    template = Path(template)
    if not template.exists():
        raise FileNotFoundError(f"{target}: template not found: {template}")
    arr = load_structure(template)
    bc, tc = identify_chains(arr, binder_chain, target_chain,
                             label=f"[{target}] {template.name}")

    binder_mask = (arr.chain_id == bc)
    if scaffold_len is not None:                 # trim to N-terminal domain
        binder_mask &= (arr.res_id <= scaffold_len)
    binder = arr[binder_mask]
    target_arr = arr[arr.chain_id == tc]
    if len(binder) == 0 or len(target_arr) == 0:
        raise RuntimeError(f"{target}: empty binder or target selection from "
                           f"{template.name} (chains {bc}/{tc})")

    binder.chain_id[:] = DESIGN_BINDER_CHAIN      # "A"
    target_arr.chain_id[:] = DESIGN_TARGET_CHAIN  # "B"
    combined = renumber(binder + target_arr)

    target_len = len(np.unique(combined.res_id[combined.chain_id == DESIGN_TARGET_CHAIN]))
    binder_len = len(np.unique(combined.res_id[combined.chain_id == DESIGN_BINDER_CHAIN]))
    logger.info(f"[{target}] construct from {template.name}: binder={binder_len} aa, "
                f"target={target_len} aa -> {out_pdb.name}")

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    pf = PDBFile()
    pf.set_structure(combined)
    pf.write(str(out_pdb))
    return out_pdb, target_len


# ── Fixed-length contig (loop lengths pinned to lengths[name]) ─────────────────
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


# ── Generation ─────────────────────────────────────────────────────────────────
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
    the selected loops — scaffold + flanks + entire target — held fixed).  All
    selected loops are designed together, so their interactions are modelled.

    The engine writes no files (write_structures/write_fasta are False); the
    sequences come back in memory and are persisted by run_probe.
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


# ── Construct preparation + geometry (shared by run_probe and the sweep) ───────
def prepare_construct(target: str, active_loops: list[str],
                      template_set: str = DEFAULT_TEMPLATE_SET,
                      template_dir: str | Path | None = None,
                      template_map: dict | None = None,
                      binder_chain: str | None = None,
                      target_chain: str | None = None,
                      scaffold_len: int | None = None) -> dict:
    """
    Resolve the template, build (or reuse) the design construct, and derive this
    template's loop geometry.  `scaffold_len=None` (the default) keeps the FULL
    TIMP3 chain (188 aa); an int trims to that many N-terminal residues.

    Returns {input_pdb, target_len, selected_loops, scaffold_len, template,
    binder_seq}.  Cheap and idempotent — the prepared PDB is cached on disk, so
    the sweep can call this once per target to learn the native loop lengths
    before running anything.
    """
    if scaffold_len is not None:
        need = required_scaffold_len(active_loops)
        if scaffold_len < need:
            raise ValueError(f"--scaffold-len {scaffold_len} is too short for loops "
                             f"{active_loops}: need >= {need} to contain them "
                             f"(or leave it unset to keep full-length TIMP3).")

    template = resolve_template(target, template_set, template_dir, template_map)
    tag = f"N{scaffold_len}" if scaffold_len is not None else "full"
    input_pdb = OUT_BASE / "inputs" / f"{target}_{template.stem}_{tag}.pdb"
    if not input_pdb.exists():
        prepare_input_pdb(target, template, input_pdb, scaffold_len,
                          binder_chain, target_chain)

    prepared = PDBFile.read(str(input_pdb)).get_structure()[0]
    binder_seq = get_seq(prepared, DESIGN_BINDER_CHAIN)
    target_len = len(np.unique(
        prepared.res_id[prepared.chain_id == DESIGN_TARGET_CHAIN]))
    selected_loops = derive_loop_geometry(binder_seq, active_loops)
    # The construct length actually realised (full binder chain, or the trim).
    actual_scaffold_len = len(binder_seq)
    return {"input_pdb": input_pdb, "target_len": target_len,
            "selected_loops": selected_loops, "scaffold_len": actual_scaffold_len,
            "template": template, "binder_seq": binder_seq}


def template_native_lengths(target: str, active_loops: list[str], **kwargs) -> dict[str, int]:
    """Native loop lengths as they appear in this target's template construct."""
    con = prepare_construct(target, active_loops, **kwargs)
    return {g["name"]: g["normal"] for g in con["selected_loops"]}


# ── Top-level probe (one target, one length configuration) ─────────────────────
def run_probe(target: str, active_loops: list[str],
              lengths: dict[str, int] | None = None,
              n_backbones: int = DEFAULT_N_BACKBONES,
              seqs_per_backbone: int = DEFAULT_SEQS_PER_BACKBONE,
              temperature: float = DEFAULT_TEMPERATURE,
              out_dir: Path = None,
              template_set: str = DEFAULT_TEMPLATE_SET,
              template_dir: str | Path | None = None,
              template_map: dict | None = None,
              binder_chain: str | None = None, target_chain: str | None = None,
              scaffold_len: int | None = None, seed: int = DEFAULT_SEED,
              make_plots: bool = True) -> dict:
    """
    Run RFd3->LigandMPNN at fixed loop lengths for one target and write, per
    loop: sequences.csv, position-count/frequency CSVs, and heatmaps.
    `lengths` defaults to each loop's template-native length.
    Returns a summary dict (also written as summary.json).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Prepared constructs are cached per (target, template, construct length) so
    # a whole sweep reuses a single trim per target.  Loop positions and native
    # lengths come from THIS template's own sequence.
    con = prepare_construct(target, active_loops, template_set, template_dir,
                            template_map, binder_chain, target_chain, scaffold_len)
    input_pdb      = con["input_pdb"]
    target_len     = con["target_len"]
    selected_loops = con["selected_loops"]
    scaffold_len   = con["scaffold_len"]
    template       = con["template"]
    lengths = resolve_lengths(active_loops, lengths, selected_loops)

    contig, length_range = build_fixed_contig(selected_loops, scaffold_len,
                                              target_len, lengths)
    ordered = [lc["name"] for lc in selected_loops]
    length_sig = "_".join(f"{n}{lengths[n]}" for n in ordered)
    logger.info(f"[{target}] loops={ordered} lengths={lengths} T={temperature} "
                f"backbones={n_backbones}x{seqs_per_backbone} template={template.name}")
    logger.info(f"[{target}] contig: {contig}")

    backbones = run_rfd3_fixed(input_pdb, contig, length_range, n_backbones)
    tag = f"{target}_{length_sig}_T{temperature}"
    designs = run_lmpnn_fixed(backbones, selected_loops, out_dir,
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

    # ── Count + heatmap per loop ───────────────────────────────────────────────
    per_loop = {}
    for lc in selected_loops:
        name, L = lc["name"], lengths[lc["name"]]
        col = f"loop_{name}_seq"
        raw = [s for s in seq_df.get(col, pd.Series([], dtype=str)).astype(str).tolist()
               if s and s != "MISSING"]
        good = [s for s in raw if len(s) == L and set(s) <= lpa.AA_SET]
        counts = lpa.count_positions(good, length=L)
        counts.attrs["n_sequences"] = len(good)
        n_uniq = len(set(good))
        per_loop[name] = {"n_designs": len(raw), "n_usable": len(good),
                          "n_parse_fail": len(raw) - len(good),
                          "n_unique": n_uniq, "length": L}
        if make_plots:
            lpa.render_all_heatmaps(counts, out_dir, stem=name,
                                    title_prefix=f"{target} | {name} L{L} | ")
        else:
            lpa.to_frequency(counts).to_csv(out_dir / f"position_freq_{name}.csv")
            counts.to_csv(out_dir / f"position_counts_{name}.csv")

    summary = {
        "target": target, "loops": ordered, "lengths": lengths,
        "template": str(template), "template_set": template_set,
        "scaffold_len": scaffold_len, "temperature": temperature,
        "n_backbones": n_backbones, "seqs_per_backbone": seqs_per_backbone,
        "n_designs_total": len(designs), "contig": contig,
        "length_range": length_range, "seed": seed,
        "loop_geometry": {g["name"]: {"pos": g["pos"], "normal": g["normal"]}
                          for g in selected_loops},
        "per_loop": per_loop,
    }
    with open(out_dir / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    logger.info(f"[{target}] done -> {out_dir}  " + " ".join(
        f"{k}:{v['n_usable']}/{v['n_designs']}(uniq {v['n_unique']})"
        for k, v in per_loop.items()))
    return summary


# ── CLI ────────────────────────────────────────────────────────────────────────
def parse_kv_ints(spec: str | dict | None) -> dict[str, int]:
    """Parse 'AB=8,C=6' (or pass through a dict from a config file)."""
    if not spec:
        return {}
    if isinstance(spec, dict):
        return {k: int(v) for k, v in spec.items()}
    out = {}
    for tok in str(spec).replace(" ", "").split(","):
        if not tok:
            continue
        name, _, val = tok.partition("=")
        out[name] = int(val)
    return out


def add_common_args(ap: argparse.ArgumentParser) -> None:
    """Arguments shared by loop_probe and loop_probe_sweep (defaults stay None
    so the config file can supply them; see cfg_get)."""
    ap.add_argument("--config", default=None, help="JSON/YAML run config")
    ap.add_argument("--loops", nargs="+", default=None,
                    help="loops to design; join with '_' to sweep jointly "
                         "(e.g. AB_C EF). Default: AB C EF")
    ap.add_argument("--n-backbones", type=int, default=None,
                    help=f"RFd3 backbones per config (default {DEFAULT_N_BACKBONES})")
    ap.add_argument("--seqs-per-backbone", type=int, default=None,
                    help=f"LigandMPNN sequences per backbone (default {DEFAULT_SEQS_PER_BACKBONE})")
    ap.add_argument("--temperature", type=float, default=None,
                    help=f"LigandMPNN sampling temperature (default {DEFAULT_TEMPERATURE})")
    ap.add_argument("--template-set", default=None, choices=list(TEMPLATE_SETS),
                    help=f"built-in template set (default {DEFAULT_TEMPLATE_SET}: "
                         f"full-length TIMP3 AF3 complexes)")
    ap.add_argument("--template-dir", default=None,
                    help="directory of template structures (overrides the set's dir)")
    ap.add_argument("--binder-chain", default=None,
                    help="force the TIMP3/binder chain id (default: auto-detect)")
    ap.add_argument("--target-chain", default=None,
                    help="force the target chain id (default: auto-detect)")
    ap.add_argument("--scaffold-len", type=int, default=None,
                    help="N-terminal TIMP3 residues to keep (default: full length, "
                         "188 aa; e.g. 121 trims to the N-terminal domain)")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--no-plots", action="store_true",
                    help="write CSVs only (rebuild figures later with loop_probe_analysis)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", default=None, choices=TARGET_NAMES,
                    help="which target complex to probe")
    ap.add_argument("--lengths", default=None,
                    help="fix loop lengths, e.g. 'AB=8,C=6'; unspecified use native")
    add_common_args(ap)
    args = ap.parse_args()

    cfg = load_config(args.config)
    target = cfg_get(args.target, cfg, "target", None)
    if not target:
        ap.error("--target is required (or set 'target' in the config)")

    loop_tokens = cfg_get(args.loops, cfg, "loops", ["AB", "C", "EF"])
    if isinstance(loop_tokens, str):
        loop_tokens = loop_tokens.split()
    _, design_loops = parse_loop_tokens(loop_tokens)
    lengths = parse_kv_ints(cfg_get(args.lengths, cfg, "lengths", None))

    setup_env()
    length_sig = "_".join(f"{n}{lengths[n]}" for n in design_loops if n in lengths) or "native"
    out_dir = cfg_get(args.out_dir, cfg, "out_dir", None)
    out_dir = Path(out_dir) if out_dir else OUT_BASE / f"{target}_{length_sig}"

    run_probe(
        target=target, active_loops=design_loops, lengths=lengths, out_dir=out_dir,
        n_backbones=cfg_get(args.n_backbones, cfg, "n_backbones", DEFAULT_N_BACKBONES),
        seqs_per_backbone=cfg_get(args.seqs_per_backbone, cfg, "seqs_per_backbone",
                                  DEFAULT_SEQS_PER_BACKBONE),
        temperature=cfg_get(args.temperature, cfg, "temperature", DEFAULT_TEMPERATURE),
        template_set=cfg_get(args.template_set, cfg, "template_set", DEFAULT_TEMPLATE_SET),
        template_dir=cfg_get(args.template_dir, cfg, "template_dir", None),
        template_map=cfg.get("template_map"),
        binder_chain=cfg_get(args.binder_chain, cfg, "binder_chain", None),
        target_chain=cfg_get(args.target_chain, cfg, "target_chain", None),
        scaffold_len=cfg_get(args.scaffold_len, cfg, "scaffold_len", None),
        seed=cfg_get(args.seed, cfg, "seed", DEFAULT_SEED),
        make_plots=not (args.no_plots or cfg.get("no_plots", False)),
    )


if __name__ == "__main__":
    main()
