"""
score_with_esmfold2.py

Score a fixed set of binder-target designs with ESMFold2 (Chan Zuckerberg Biohub,
Rives lab) so its ipTM/pTM/pLDDT can be compared head-to-head against Boltz-2 and
AF3 on identical sequences — at zero extra AF3 cost.

WHY ESMFold2 here: it is a protein-language-model folder, MSA-free *by design*
(no ColabFold dependency). Biohub reports it is strongest on no-MSA targets like
antibodies — the closest analog to our engineered TIMP3-loop binders, which have no
natural homologs. Our calibration showed Boltz is only a weak, target-inconsistent
filter (see filter_methods.md §4); ESMFold2 may separate AF3 hits better.

RUN THIS ON THE CLUSTER (needs a GPU; ESMC-6B base is large — verify it fits VRAM).
Install in a SEPARATE env to avoid the numpy/cublas/lightning conflicts that hit the
foundry env. ESMFold2 requires Python >= 3.12:
    conda create -n esmfold2 python=3.12 -y
    conda activate esmfold2
    pip install "esm @ git+https://github.com/Biohub/esm.git@main"
    pip install transformers torch pandas

Then, from the repo root (so it can find Local/iterative_refinement/):
    python Generation/score_with_esmfold2.py            # scores the stratified manifest
    python Generation/score_with_esmfold2.py --all       # scores every design in round summaries

Output: Local/iterative_refinement/esmfold2_scores.csv
        columns: design_id, target_name, full_seq, esm_iptm, esm_ptm, esm_plddt,
                 esm_lplddt, esm_pae
        esm_lplddt = pLDDT over ONLY the redesigned-loop residues; esm_pae = mean
        interface PAE between those loops and the target (NaN if the model emits no
        PAE matrix). These are the real-dynamic-range inputs the binding recipe
        ranks on, so the SAME multi-term recipe runs on ESMFold2 as on AF3 (not
        just per-metric). Loops are located by flank tripeptides; pass a
        `design_loops` input column (e.g. "AB" or "AB,C,EF") to restrict which
        loops (default: all active loops).
To validate against AF3: parse an AF3 results zip, join its jobs to this CSV by
binder sequence, and correlate each local metric (esm_iptm, esm_plddt) with AF3
ipTM (Spearman + top-K hit-rate, pooled across batches). See filter_methods.md
for the calibration method and results.

API NOTE: the ESMFold2 call below follows the published API
(github.com/atong01/esmfold2 README, 2026-06). If the installed package's module
paths differ, only `load_model()` and `predict_complex()` need editing — they are
isolated for exactly this reason.
"""

import re
import sys
import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

OUT_BASE = (Path(__file__).parent / ".." / "Local" / "iterative_refinement").resolve()

# Designs are submitted binder-first, so chain A = binder, chain B = target.
BINDER_CHAIN = "A"
TARGET_CHAIN = "B"


# ── Loop-focused metrics: LpLDDT + interface PAE ──────────────────────────────
# The binding recipe (Generation/binding_recipe.py) ranks on loop pLDDT and
# interface PAE — the only AF metrics with real dynamic range. ESMFold2 gives a
# per-residue pLDDT array and (if available) a PAE matrix; we reduce those to the
# same two numbers over the SAME flank-tripeptide loop definitions as the rest of
# the pipeline, so the ESMFold2 recipe score is comparable to the AF3 one.
try:
    from select_binders_to_order import loop_residue_positions, ACTIVE_LOOPS, LOOP_CONFIGS
except Exception:   # keep this scorer runnable in isolation (no sibling import)
    LOOP_CONFIGS = {
        "AB": {"left": "LVK", "right": "LVY"}, "C": {"left": "HTE", "right": "GLK"},
        "EF": {"left": "MYT", "right": "FVE"}, "GH": {"left": "KSC", "right": "NEC"}}
    ACTIVE_LOOPS = ["AB", "C", "EF"]

    def loop_residue_positions(binder_seq, loops=ACTIVE_LOOPS):
        out, cursor = {}, 0
        for name in loops:
            lc = LOOP_CONFIGS[name]
            m = re.compile(f"{lc['left']}([A-Z]*?){lc['right']}").search(binder_seq[cursor:])
            if m:
                start = cursor + m.start() + len(lc["left"]) + 1   # 1-indexed loop start
                out[name] = set(range(start, start + len(m.group(1))))
                cursor = cursor + m.end() - len(lc["right"])
            else:
                out[name] = set()
        return out

_LOOP_PAE_WARNED = {"done": False}


def _parse_design_loops(design_loops):
    """Normalize 'AB' / 'AB,C,EF' / ['AB','C'] -> known-loop list (default: all active)."""
    if design_loops is None or (isinstance(design_loops, float) and np.isnan(design_loops)):
        return list(ACTIVE_LOOPS)
    toks = re.split(r"[^A-Za-z]+", design_loops) if isinstance(design_loops, str) else list(design_loops)
    loops = [str(t).upper() for t in toks if t and str(t).upper() in LOOP_CONFIGS]
    return loops or list(ACTIVE_LOOPS)


def _loop_positions_0indexed(binder_seq, design_loops):
    """0-indexed binder positions of the redesigned loop(s), located by flanks."""
    pos1 = loop_residue_positions(binder_seq, _parse_design_loops(design_loops))
    allpos = sorted(set().union(*pos1.values())) if pos1 else []
    return [p - 1 for p in allpos]


def _binder_plddt_0_100(plddt_arr, binder_len):
    """Binder per-residue pLDDT on a 0-100 scale (None if unavailable)."""
    if plddt_arr is None:
        return None
    arr = np.asarray(plddt_arr, dtype=float).flatten()
    if arr.size == 0:
        return None
    arr = arr[:binder_len] if arr.size >= binder_len else arr
    finite = arr[np.isfinite(arr)]
    if finite.size and np.nanmax(finite) <= 1.0:        # 0-1 -> 0-100
        arr = arr * 100.0
    return arr


def _loop_plddt(plddt_res, loop_pos0):
    """Mean pLDDT (0-100) over loop residues; NaN if unavailable."""
    if plddt_res is None or not loop_pos0:
        return float("nan")
    vals = [plddt_res[i] for i in loop_pos0 if 0 <= i < plddt_res.size and np.isfinite(plddt_res[i])]
    return float(np.mean(vals)) if vals else float("nan")


def _extract_pae(res):
    """Best-effort [N, N] PAE matrix from the fold result, or None (API-defensive)."""
    cand = None
    for owner in (res, getattr(res, "complex", None), getattr(res, "structure", None)):
        if owner is None:
            continue
        for name in ("pae", "predicted_aligned_error", "pae_matrix", "aligned_error"):
            v = getattr(owner, name, None)
            if v is not None:
                cand = v
                break
        if cand is not None:
            break
    if cand is None:
        return None
    try:
        mat = np.asarray(cand, dtype=float)
    except Exception:
        try:
            mat = cand.detach().cpu().numpy().astype(float)     # torch tensor
        except Exception:
            return None
    mat = np.squeeze(mat)
    return mat if mat.ndim == 2 else None


def _loop_interface_pae(pae, loop_pos0, binder_len, target_len):
    """Mean PAE (both directions) between binder loop residues and the target chain."""
    if pae is None or not loop_pos0:
        return float("nan")
    total = binder_len + target_len
    if pae.shape != (total, total):
        # token count != residue count (e.g. chain-break tokens). Don't guess the
        # offset -- emit NaN; the raw matrix is saved as .npy for offline recompute.
        if not _LOOP_PAE_WARNED["done"]:
            print(f"  note: PAE shape {tuple(pae.shape)} != "
                  f"binder+target residues ({total}); skipping interface PAE.")
            _LOOP_PAE_WARNED["done"] = True
        return float("nan")
    a_idx = np.array([i for i in loop_pos0 if 0 <= i < binder_len], dtype=int)
    b_idx = np.arange(binder_len, total)
    if a_idx.size == 0 or b_idx.size == 0:
        return float("nan")
    return float((pae[np.ix_(a_idx, b_idx)].mean() + pae[np.ix_(b_idx, a_idx)].mean()) / 2.0)


# ── ESMFold2 interface (isolated; edit here if the installed API differs) ─────

def load_model(device: str = "cuda"):
    """Load ESMFold2 once. Returns (model, input_builder)."""
    from esm.models.esmfold2 import ESMFold2InputBuilder           # noqa: F401
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
    model = ESMFold2Model.from_pretrained("biohub/ESMFold2").to(device).eval()
    return model, ESMFold2InputBuilder()


def _save_structure(res, stub: str):
    """
    Write the predicted complex to disk if a structure is available. The folded
    structure is already computed during fold() — saving it is just disk I/O,
    no extra GPU time. Returns the written path or None. Defensive about the
    exact serialization API (verify against your installed `esm` if it differs).
    """
    comp = getattr(res, "complex", None) or getattr(res, "structure", None)
    if comp is None:
        return None
    # Try mmCIF text, then PDB text, then a direct save method.
    for meth, ext in (("to_mmcif", "cif"), ("to_mmcif_string", "cif"),
                      ("to_pdb", "pdb"), ("to_pdb_string", "pdb")):
        fn = getattr(comp, meth, None)
        if callable(fn):
            try:
                text = fn()
                if isinstance(text, str) and text.strip():
                    path = f"{stub}.{ext}"
                    with open(path, "w") as fh:
                        fh.write(text)
                    return path
            except Exception:
                continue
    for meth, ext in (("save_mmcif", "cif"), ("save_pdb", "pdb"), ("save", "cif")):
        fn = getattr(comp, meth, None)
        if callable(fn):
            try:
                path = f"{stub}.{ext}"
                fn(path)
                return path
            except Exception:
                continue
    return None


def predict_complex(model, builder, binder_seq: str, target_seq: str,
                    binder_len: int, seed: int = 0, cif_stub: str = None,
                    design_loops=None) -> dict:
    """
    Fold a two-chain complex and return
        {esm_iptm, esm_ptm, esm_plddt, esm_lplddt, esm_pae[, esm_cif]}.

      esm_plddt   binder-chain mean pLDDT (0-100).
      esm_lplddt  pLDDT over ONLY the redesigned-loop residues (0-100) -- the
                  loop-focused confidence the binding recipe ranks on.
                  `design_loops` restricts which loops (default: all active loops).
      esm_pae     mean interface PAE (A) between those loop residues and the target
                  chain, both directions. NaN if the model emits no PAE matrix.

    Field access is defensive so minor API drift doesn't crash the run. If
    `cif_stub` is given, the predicted structure (and PAE matrix, if any) is written
    to disk (free -- already computed by fold()).
    """
    from esm.models.esmfold2 import ProteinInput, StructurePredictionInput

    spi = StructurePredictionInput(sequences=[
        ProteinInput(id=BINDER_CHAIN, sequence=binder_seq),
        ProteinInput(id=TARGET_CHAIN, sequence=target_seq),
    ])
    res = builder.fold(model, spi, num_loops=3, num_sampling_steps=50,
                       num_diffusion_samples=1, seed=seed)

    def _get(obj, *names):
        for n in names:
            v = getattr(obj, n, None)
            if v is not None:
                return v
        return None

    iptm = _get(res, "iptm", "interface_ptm")
    ptm  = _get(res, "ptm")

    # Per-residue pLDDT (binder chain first) -> whole-binder mean + loop-only mean.
    plddt_res = _binder_plddt_0_100(_get(res, "plddt", "plddts"), binder_len)
    plddt = float(np.nanmean(plddt_res)) if plddt_res is not None and plddt_res.size else np.nan
    loop_pos0 = _loop_positions_0indexed(binder_seq, design_loops)
    lplddt = _loop_plddt(plddt_res, loop_pos0)

    # Interface PAE over the loop residues (needs a PAE matrix; NaN if unavailable).
    pae_mat = _extract_pae(res)
    loop_pae = _loop_interface_pae(pae_mat, loop_pos0, binder_len, len(target_seq))

    def _f(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return np.nan

    out = {"esm_iptm": _f(iptm), "esm_ptm": _f(ptm), "esm_plddt": plddt,
           "esm_lplddt": lplddt, "esm_pae": loop_pae}
    if cif_stub:
        saved = _save_structure(res, cif_stub)
        if saved:
            out["esm_cif"] = saved
        if pae_mat is not None:                 # keep the raw matrix for offline recompute
            try:
                np.save(f"{cif_stub}_pae.npy", pae_mat)
            except Exception:
                pass
    return out


# ── Design pool ───────────────────────────────────────────────────────────────

def load_pool() -> pd.DataFrame:
    """All designs from round summaries (full_seq, target_seq, design_id, target)."""
    frames = []
    for it_dir in sorted(OUT_BASE.glob("it_*")):
        csv = it_dir / "round_summary.csv"
        if csv.exists():
            frames.append(pd.read_csv(csv))
    if not frames:
        sys.exit(f"No round_summary.csv under {OUT_BASE}")
    df = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["full_seq"])
    return df[["design_id", "target_name", "full_seq", "target_seq"]]


def select_designs(pool: pd.DataFrame, score_all: bool) -> pd.DataFrame:
    """Default: the stratified-manifest designs (those with AF3 ground truth)."""
    if score_all:
        return pool
    man_path = OUT_BASE / "stratified_manifest.json"
    if not man_path.exists():
        print(f"No {man_path.name}; falling back to scoring all designs.")
        return pool
    man = pd.DataFrame(json.loads(man_path.read_text()))
    sel = pool[pool["full_seq"].isin(set(man["full_seq"]))]
    print(f"Scoring {len(sel)} designs from stratified manifest.")
    return sel


def main():
    ap = argparse.ArgumentParser(description="Score designs with ESMFold2.")
    ap.add_argument("--all", action="store_true", help="Score every design, not just the manifest.")
    ap.add_argument("--input", default=None,
                    help="Explicit CSV of designs to score (cols: design_id, full_seq, "
                         "target_seq[, target_name]). Used by the pipeline's run_esmfold2.")
    ap.add_argument("--device", default="cuda", help="torch device (cuda / cuda:0 / cpu).")
    ap.add_argument("--out", default=str(OUT_BASE / "esmfold2_scores.csv"))
    ap.add_argument("--cif-dir", default=None,
                    help="If set, write each predicted complex CIF/PDB here "
                         "(free — structure is already computed). Path recorded as esm_cif.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--resume", action="store_true",
                    help="Skip design_ids already present in --out (the CSV is appended to "
                         "as each fold completes, so a killed run restarts where it stopped).")
    args = ap.parse_args()

    # Fail fast on a broken CUDA/torch setup BEFORE the slow model download.
    if args.device.startswith("cuda"):
        import torch
        if not torch.cuda.is_available():
            sys.exit(
                "torch.cuda.is_available() is False — the torch build does not match "
                "this machine's NVIDIA driver.\n"
                f"  torch={torch.__version__}, compiled for CUDA {torch.version.cuda}\n"
                "  Fix: install a torch built for the cluster's CUDA (match the foundry env), e.g.\n"
                "    pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cu124\n"
                "  Or pass --device cpu (very slow)."
            )
        try:
            torch.zeros(1).to(args.device)   # surfaces driver-too-old before the download
        except RuntimeError as exc:
            sys.exit(f"CUDA init failed on {args.device}: {exc}\n"
                     "  Reinstall torch to match the cluster driver (see above).")
        print(f"CUDA OK: {torch.cuda.get_device_name(0)} | torch {torch.__version__} (CUDA {torch.version.cuda})")

    if args.input:
        designs = pd.read_csv(args.input)
        if "target_name" not in designs.columns:
            designs["target_name"] = ""
        print(f"Scoring {len(designs)} designs from {args.input}")
    else:
        designs = select_designs(load_pool(), args.all)
    if designs.empty:
        sys.exit("No designs selected to score.")

    if args.resume and Path(args.out).exists():
        try:
            done = set(pd.read_csv(args.out, usecols=["design_id"])["design_id"].astype(str))
        except Exception:
            done = set()
        if done:
            before = len(designs)
            designs = designs[~designs["design_id"].astype(str).isin(done)]
            print(f"--resume: {len(done)} already scored in {args.out}; "
                  f"{len(designs)}/{before} left.")
        if designs.empty:
            print("Nothing left to score (all design_ids already in --out).")
            return

    print(f"Loading ESMFold2 on {args.device} ... (large model; first load is slow)")
    model, builder = load_model(args.device)

    cif_dir = None
    if args.cif_dir:
        cif_dir = Path(args.cif_dir)
        cif_dir.mkdir(parents=True, exist_ok=True)

    has_loops = "design_loops" in designs.columns
    out_path = Path(args.out)
    n_ok = 0
    for i, d in enumerate(designs.itertuples(index=False), 1):
        bseq, tseq = d.full_seq, d.target_seq
        if not isinstance(bseq, str) or not isinstance(tseq, str) or not bseq or not tseq:
            continue
        dloops = getattr(d, "design_loops", None) if has_loops else None
        try:
            stub = str(cif_dir / str(d.design_id)) if cif_dir else None
            metrics = predict_complex(model, builder, bseq, tseq, len(bseq),
                                      seed=args.seed, cif_stub=stub, design_loops=dloops)
            row = {"design_id": d.design_id, "target_name": d.target_name,
                   "full_seq": bseq, "esm_iptm": np.nan, "esm_ptm": np.nan,
                   "esm_plddt": np.nan, "esm_lplddt": np.nan, "esm_pae": np.nan,
                   "esm_cif": "", **metrics}
            # Append-as-we-go: a hard process kill (CUDA OOM / cgroup OOM-killer)
            # then costs one design instead of the whole batch. Pair with --resume.
            # The key order above is fixed so appended rows always match the header
            # even when a structure fails to save (no esm_cif in `metrics`).
            pd.DataFrame([row]).to_csv(out_path, mode="a", index=False,
                                       header=not out_path.exists())
            n_ok += 1
            print(f"  [{i}/{len(designs)}] {d.design_id}: "
                  f"ipTM={metrics['esm_iptm']:.3f}  pLDDT={metrics['esm_plddt']:.1f}  "
                  f"LpLDDT={metrics['esm_lplddt']:.1f}  loopPAE={metrics['esm_pae']:.1f}",
                  flush=True)
        except Exception as exc:
            print(f"  [{i}/{len(designs)}] {d.design_id}: ERROR {exc}", flush=True)

    if out_path.exists():
        print(f"\nWrote {n_ok} ESMFold2 scores -> {out_path}")
        print("Validate against AF3 by joining an AF3 results zip to this CSV by "
              "binder sequence (Spearman + top-K vs AF3 ipTM); see filter_methods.md.")
    else:
        print("\nNo scores produced. Check the ESMFold2 install / API (see load_model).")


if __name__ == "__main__":
    main()
