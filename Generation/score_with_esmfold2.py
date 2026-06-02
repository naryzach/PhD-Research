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
        columns: design_id, target_name, full_seq, esm_iptm, esm_ptm, esm_plddt
To validate against AF3: parse an AF3 results zip, join its jobs to this CSV by
binder sequence, and correlate each local metric (esm_iptm, esm_plddt) with AF3
ipTM (Spearman + top-K hit-rate, pooled across batches). See filter_methods.md
for the calibration method and results.

API NOTE: the ESMFold2 call below follows the published API
(github.com/atong01/esmfold2 README, 2026-06). If the installed package's module
paths differ, only `load_model()` and `predict_complex()` need editing — they are
isolated for exactly this reason.
"""

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
                    binder_len: int, seed: int = 0, cif_stub: str = None) -> dict:
    """
    Fold a two-chain complex and return {esm_iptm, esm_ptm, esm_plddt[, esm_cif]}.
    pLDDT is the binder-chain mean (0-100). Field access is defensive so minor
    API drift doesn't crash the run. If `cif_stub` is given, the predicted
    structure is written to disk (free — it's already computed by fold()).
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

    # pLDDT: per-residue array → binder-chain mean → 0-100
    plddt_arr = _get(res, "plddt", "plddts")
    plddt = np.nan
    if plddt_arr is not None:
        arr = np.asarray(plddt_arr, dtype=float).flatten()
        if arr.size >= binder_len:
            arr = arr[:binder_len]
        if arr.size:
            plddt = float(arr.mean())
            if 0.0 < plddt <= 1.0:
                plddt *= 100.0

    def _f(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return np.nan

    out = {"esm_iptm": _f(iptm), "esm_ptm": _f(ptm), "esm_plddt": plddt}
    if cif_stub:
        saved = _save_structure(res, cif_stub)
        if saved:
            out["esm_cif"] = saved
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

    print(f"Loading ESMFold2 on {args.device} ... (large model; first load is slow)")
    model, builder = load_model(args.device)

    cif_dir = None
    if args.cif_dir:
        cif_dir = Path(args.cif_dir)
        cif_dir.mkdir(parents=True, exist_ok=True)

    rows, n_ok = [], 0
    for i, d in enumerate(designs.itertuples(index=False), 1):
        bseq, tseq = d.full_seq, d.target_seq
        if not isinstance(bseq, str) or not isinstance(tseq, str) or not bseq or not tseq:
            continue
        try:
            stub = str(cif_dir / str(d.design_id)) if cif_dir else None
            metrics = predict_complex(model, builder, bseq, tseq, len(bseq),
                                      seed=args.seed, cif_stub=stub)
            rows.append({"design_id": d.design_id, "target_name": d.target_name,
                         "full_seq": bseq, **metrics})
            n_ok += 1
            print(f"  [{i}/{len(designs)}] {d.design_id}: "
                  f"ipTM={metrics['esm_iptm']:.3f}  pLDDT={metrics['esm_plddt']:.1f}")
        except Exception as exc:
            print(f"  [{i}/{len(designs)}] {d.design_id}: ERROR {exc}")

    if rows:
        pd.DataFrame(rows).to_csv(args.out, index=False)
        print(f"\nWrote {n_ok} ESMFold2 scores -> {args.out}")
        print("Validate against AF3 by joining an AF3 results zip to this CSV by "
              "binder sequence (Spearman + top-K vs AF3 ipTM); see filter_methods.md.")
    else:
        print("\nNo scores produced. Check the ESMFold2 install / API (see load_model).")


if __name__ == "__main__":
    main()
