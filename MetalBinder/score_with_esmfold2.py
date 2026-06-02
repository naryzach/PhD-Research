"""
score_with_esmfold2.py  (MetalBinder)

Score metal-binding design sequences with ESMFold2 (Chan Zuckerberg Biohub / Rives lab)
using single-chain (monomer) folding.  Called as a subprocess by run_pipeline.py because
ESMFold2 requires Python >= 3.12 and different GPU library versions than the foundry env.

Install in a SEPARATE conda env:
    conda create -n esmfold2 python=3.12 -y
    conda activate esmfold2
    pip install "esm @ git+https://github.com/Biohub/esm.git@main"
    pip install transformers torch pandas

Usage (invoked via subprocess from run_pipeline.py):
    python MetalBinder/score_with_esmfold2.py \\
        --input /path/to/esm_input.csv \\
        --out   /path/to/esm_scores.csv \\
        [--cif-dir /path/to/structures] \\
        [--device cuda]

Input CSV columns : design_id, sequence
Output CSV columns: design_id, esm_plddt, esm_ptm[, esm_cif]
"""

import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


# ── ESMFold2 interface (isolated; edit here if the installed API differs) ──────

def load_model(device: str = "cuda"):
    """Load ESMFold2 once. Returns (model, input_builder)."""
    from esm.models.esmfold2 import ESMFold2InputBuilder
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
    model = ESMFold2Model.from_pretrained("biohub/ESMFold2").to(device).eval()
    return model, ESMFold2InputBuilder()


def _save_structure(res, stub: str):
    """Write predicted structure to disk. Returns path or None."""
    comp = getattr(res, "complex", None) or getattr(res, "structure", None)
    if comp is None:
        return None
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


def predict_monomer(model, builder, sequence: str, seed: int = 0,
                    cif_stub: str = None) -> dict:
    """
    Fold a single-chain protein and return {esm_plddt, esm_ptm[, esm_cif]}.
    pLDDT is the per-residue mean (0–100 scale). If cif_stub is given, the
    predicted structure is written to {cif_stub}.cif or .pdb at no extra GPU cost.
    """
    from esm.models.esmfold2 import ProteinInput, StructurePredictionInput

    spi = StructurePredictionInput(sequences=[
        ProteinInput(id="A", sequence=sequence),
    ])
    res = builder.fold(model, spi, num_loops=3, num_sampling_steps=50,
                       num_diffusion_samples=1, seed=seed)

    def _get(obj, *names):
        for n in names:
            v = getattr(obj, n, None)
            if v is not None:
                return v
        return None

    ptm = _get(res, "ptm")

    plddt_arr = _get(res, "plddt", "plddts")
    plddt = np.nan
    if plddt_arr is not None:
        arr = np.asarray(plddt_arr, dtype=float).flatten()
        if arr.size:
            plddt = float(arr.mean())
            if 0.0 < plddt <= 1.0:
                plddt *= 100.0

    def _f(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return np.nan

    out = {"esm_plddt": plddt, "esm_ptm": _f(ptm)}
    if cif_stub:
        saved = _save_structure(res, cif_stub)
        if saved:
            out["esm_cif"] = saved
    return out


def main():
    ap = argparse.ArgumentParser(
        description="Score metal-binding designs with ESMFold2 (monomer folding).")
    ap.add_argument("--input", required=True,
                    help="Input CSV: columns design_id, sequence")
    ap.add_argument("--out", required=True,
                    help="Output CSV: design_id, esm_plddt, esm_ptm[, esm_cif]")
    ap.add_argument("--device", default="cuda",
                    help="Torch device (cuda / cuda:0 / cpu)")
    ap.add_argument("--cif-dir", default=None,
                    help="Write predicted structure CIFs here (free — already computed by fold())")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    if args.device.startswith("cuda"):
        import torch
        if not torch.cuda.is_available():
            sys.exit(
                "torch.cuda.is_available() is False.\n"
                "Fix: install a torch build matching the cluster CUDA version, or pass --device cpu."
            )
        print(f"CUDA OK: {torch.cuda.get_device_name(0)} | torch {torch.__version__} "
              f"(CUDA {torch.version.cuda})", flush=True)

    designs = pd.read_csv(args.input)
    if designs.empty:
        sys.exit("No designs in input CSV.")
    print(f"Scoring {len(designs)} designs with ESMFold2 ({args.device}) ...", flush=True)

    print("Loading ESMFold2 model (first load downloads weights — can be slow) ...", flush=True)
    model, builder = load_model(args.device)

    cif_dir = None
    if args.cif_dir:
        cif_dir = Path(args.cif_dir)
        cif_dir.mkdir(parents=True, exist_ok=True)

    rows, n_ok = [], 0
    for i, d in enumerate(designs.itertuples(index=False), 1):
        seq = getattr(d, "sequence", None)
        if not isinstance(seq, str) or not seq:
            continue
        try:
            stub = str(cif_dir / str(d.design_id)) if cif_dir else None
            metrics = predict_monomer(model, builder, seq, seed=args.seed, cif_stub=stub)
            rows.append({"design_id": d.design_id, **metrics})
            n_ok += 1
            print(f"  [{i}/{len(designs)}] {d.design_id}: "
                  f"pLDDT={metrics['esm_plddt']:.1f}  pTM={metrics['esm_ptm']:.3f}",
                  flush=True)
        except Exception as exc:
            print(f"  [{i}/{len(designs)}] {d.design_id}: ERROR {exc}", flush=True)

    if rows:
        pd.DataFrame(rows).to_csv(args.out, index=False)
        print(f"\nWrote {n_ok} ESMFold2 scores to {args.out}", flush=True)
    else:
        print("No scores produced. Check the ESMFold2 install / API (see load_model).",
              flush=True)


if __name__ == "__main__":
    main()
