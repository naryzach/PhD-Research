"""
ESMFold2 folding runner (runs on the GPU cluster in the `esmfold2` env).

Uses the SAME ESMFold2 API as Generation/score_with_esmfold2.py
(`biohub/ESMFold2` via the `esm` package: ESMFold2InputBuilder / ProteinInput /
StructurePredictionInput / builder.fold), so results are consistent with the
lab's existing ESMFold2 scores. The model call is isolated in load_model() /
_fold() — edit only there if your installed API differs.

Modes:
  monomer  fold each sequence alone  -> results/<id>.pdb + esmfold2_metrics.csv
  complex  co-fold construct+target  -> results_complex/<c>__<t>.cif|pdb
                                        + esmfold2_complex_metrics.csv
           (esm_iptm, esm_ptm, esm_plddt, esm_lplddt, esm_pae)

Usage:
  python run_esmfold2.py --mode both \\
      --fasta esmfold2_input.fasta --pairs esmfold2_complex_pairs.csv \\
      --outdir results --complex-outdir results_complex \\
      --metrics esmfold2_metrics.csv --complex-metrics esmfold2_complex_metrics.csv
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

BINDER_CHAIN, TARGET_CHAIN = "A", "B"


# ── data ------------------------------------------------------------------
def read_fasta(path: Path):
    recs, name, buf = [], None, []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


# ── ESMFold2 backend (isolated; edit here if the installed API differs) ---
def load_model(device: str = "cuda"):
    from esm.models.esmfold2 import ESMFold2InputBuilder            # noqa: F401
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
    model = ESMFold2Model.from_pretrained("biohub/ESMFold2").to(device).eval()
    return model, ESMFold2InputBuilder()


def _fold(model, builder, chains, seed=0):
    """chains = [(chain_id, seq), ...]; returns the raw fold result object."""
    from esm.models.esmfold2 import ProteinInput, StructurePredictionInput
    spi = StructurePredictionInput(
        sequences=[ProteinInput(id=cid, sequence=seq) for cid, seq in chains])
    return builder.fold(model, spi, num_loops=3, num_sampling_steps=50,
                        num_diffusion_samples=1, seed=seed)


def _get(obj, *names):
    for n in names:
        v = getattr(obj, n, None)
        if v is not None:
            return v
    return None


def _save_structure(res, stub: str):
    comp = _get(res, "complex", "structure")
    if comp is None:
        return None
    for meth, ext in (("to_pdb", "pdb"), ("to_pdb_string", "pdb"),
                      ("to_mmcif", "cif"), ("to_mmcif_string", "cif")):
        fn = getattr(comp, meth, None)
        if callable(fn):
            try:
                text = fn()
                if isinstance(text, str) and text.strip():
                    Path(f"{stub}.{ext}").write_text(text)
                    return f"{stub}.{ext}"
            except Exception:
                continue
    for meth, ext in (("save_pdb", "pdb"), ("save_mmcif", "cif"), ("save", "cif")):
        fn = getattr(comp, meth, None)
        if callable(fn):
            try:
                fn(f"{stub}.{ext}"); return f"{stub}.{ext}"
            except Exception:
                continue
    return None


def _plddt_0_100(arr, n):
    if arr is None:
        return None
    a = np.asarray(arr, dtype=float).flatten()
    if a.size == 0:
        return None
    a = a[:n] if a.size >= n else a
    fin = a[np.isfinite(a)]
    if fin.size and np.nanmax(fin) <= 1.0:
        a = a * 100.0
    return a


def _loop_positions(construct_seq: str, loop: str):
    """0-indexed residues of the grafted loop (exact substring; the graft is known)."""
    idx = construct_seq.find(loop) if loop else -1
    return list(range(idx, idx + len(loop))) if idx >= 0 else []


def _interface_pae(res, binder_len, target_len):
    pae = _get(res, "pae", "predicted_aligned_error", "pae_matrix")
    if pae is None:
        return float("nan")
    try:
        m = np.squeeze(np.asarray(pae, dtype=float))
    except Exception:
        return float("nan")
    total = binder_len + target_len
    if m.ndim != 2 or m.shape != (total, total):
        return float("nan")
    a, b = np.arange(binder_len), np.arange(binder_len, total)
    return float((m[np.ix_(a, b)].mean() + m[np.ix_(b, a)].mean()) / 2.0)


# ── drivers ---------------------------------------------------------------
def run_monomers(model, builder, fasta, outdir, metrics_csv, device):
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, seq in read_fasta(fasta):
        print(f"[monomer] {name} ({len(seq)} aa)", flush=True)
        try:
            res = _fold(model, builder, [(BINDER_CHAIN, seq)])
            saved = _save_structure(res, str(outdir / name))
            plddt = _plddt_0_100(_get(res, "plddt", "plddts"), len(seq))
            rows.append({"id": name, "length": len(seq),
                         "plddt_mean": round(float(np.nanmean(plddt)), 3) if plddt is not None else "",
                         "ptm": round(float(_get(res, "ptm") or float("nan")), 4),
                         "file": Path(saved).name if saved else ""})
        except Exception as e:
            rows.append({"id": name, "length": len(seq), "error": str(e)})
    _write(metrics_csv, rows, ["id", "length", "plddt_mean", "ptm", "file", "error"])


def run_complexes(model, builder, pairs_csv, outdir, metrics_csv, device):
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    with open(pairs_csv, newline="") as fh:
        pairs = list(csv.DictReader(fh))
    for p in pairs:
        cid, tid = p["construct_id"], p["target_id"]
        cseq, tseq, loop = p["construct_seq"], p["target_seq"], p.get("loop", "")
        print(f"[complex] {cid}__{tid}", flush=True)
        try:
            res = _fold(model, builder, [(BINDER_CHAIN, cseq), (TARGET_CHAIN, tseq)])
            saved = _save_structure(res, str(outdir / f"{cid}__{tid}"))
            plddt = _plddt_0_100(_get(res, "plddt", "plddts"), len(cseq))
            loop_pos = _loop_positions(cseq, loop)
            lplddt = (round(float(np.nanmean([plddt[i] for i in loop_pos
                       if 0 <= i < plddt.size])), 3)
                      if plddt is not None and loop_pos else "")
            rows.append({
                "construct_id": cid, "target_id": tid,
                "esm_iptm": round(float(_get(res, "iptm", "interface_ptm") or float("nan")), 4),
                "esm_ptm": round(float(_get(res, "ptm") or float("nan")), 4),
                "esm_plddt": round(float(np.nanmean(plddt)), 3) if plddt is not None else "",
                "esm_lplddt": lplddt,
                "esm_pae": round(_interface_pae(res, len(cseq), len(tseq)), 3),
                "file": Path(saved).name if saved else ""})
        except Exception as e:
            rows.append({"construct_id": cid, "target_id": tid, "error": str(e)})
    _write(metrics_csv, rows, ["construct_id", "target_id", "esm_iptm", "esm_ptm",
                               "esm_plddt", "esm_lplddt", "esm_pae", "file", "error"])


def _write(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})
    print(f"  wrote {len(rows)} rows -> {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["monomer", "complex", "both"], default="both")
    ap.add_argument("--fasta", type=Path)
    ap.add_argument("--pairs", type=Path)
    ap.add_argument("--outdir", type=Path, default=Path("results"))
    ap.add_argument("--complex-outdir", type=Path, default=Path("results_complex"))
    ap.add_argument("--metrics", type=Path, default=Path("esmfold2_metrics.csv"))
    ap.add_argument("--complex-metrics", type=Path,
                    default=Path("esmfold2_complex_metrics.csv"))
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()

    model, builder = load_model(args.device)
    if args.mode in ("monomer", "both") and args.fasta:
        run_monomers(model, builder, args.fasta, args.outdir, args.metrics, args.device)
    if args.mode in ("complex", "both") and args.pairs:
        run_complexes(model, builder, args.pairs, args.complex_outdir,
                      args.complex_metrics, args.device)


if __name__ == "__main__":
    main()
