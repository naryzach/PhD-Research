"""
score_with_chai1.py  (MetalBinder)

Score metal-binding designs with Chai-1 (Chai Discovery), an AF3-class model that
accepts arbitrary chemical entities — including rare-earth / transition-metal ions —
as separate input entities alongside the protein sequence.  Unlike ESMFold2, Chai-1
predicts where the metal sits relative to the protein, so the pLDDT and pTM it reports
reflect the metal-bound conformation.  This is essential for EF-hand / lanthanide
designs where the binding loops are intrinsically disordered without the metal.

Each design is scored as one protein chain + one metal ion (specified as SMILES, which
is more reliable than CCD codes for lanthanides whose CCD entries are less consistent
across databases).

Install in a SEPARATE conda env (chai-lab pins its own torch/numpy):
    conda create -n chai1 python=3.10 -y
    conda activate chai1
    pip install chai-lab

Usage (invoked via subprocess from run_pipeline.py):
    python MetalBinder/score_with_chai1.py \\
        --input  /path/to/chai_input.csv \\
        --out    /path/to/chai_scores.csv \\
        [--cif-dir /path/to/structures] \\
        [--device cuda] [--num-recycles 3] [--num-diffn-steps 200]

Input CSV columns : design_id, sequence, metal_ccd
Output CSV columns: design_id, chai_plddt, chai_ptm[, chai_score, chai_cif]

NOTE on Chai-1 FASTA format and API:
  This script targets chai-lab ≥ 0.4.  If run_inference raises TypeError on keyword
  arguments, it falls back to a bare call (older API).  Verify the FASTA header schema
  against your installed version — the metal ligand section uses SMILES in the sequence
  field, which is stable across ≥ 0.4.  See: github.com/chaidiscovery/chai-lab
"""

import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


# ── Metal ion SMILES lookup ───────────────────────────────────────────────────
# SMILES is used instead of raw CCD codes to avoid lookup ambiguity for rare
# earths (e.g. TB conflicts with nucleotide codes in some databases).
# Charge states follow the dominant oxidation state in biological coordination.
METAL_SMILES: dict[str, str] = {
    # Lanthanides / rare earths (primary MetalBinder targets)
    "EU":  "[Eu+3]",   # Europium(III)
    "TB":  "[Tb+3]",   # Terbium(III)
    "DY":  "[Dy+3]",   # Dysprosium(III)
    "ND":  "[Nd+3]",   # Neodymium(III)   — original 8FNS metal
    "GD":  "[Gd+3]",   # Gadolinium(III)
    "ER":  "[Er+3]",   # Erbium(III)
    "YB":  "[Yb+3]",   # Ytterbium(III)
    "LA":  "[La+3]",   # Lanthanum(III)
    "CE":  "[Ce+3]",   # Cerium(III)
    "SM":  "[Sm+3]",   # Samarium(III)
    "HO":  "[Ho+3]",   # Holmium(III)
    "LU":  "[Lu+3]",   # Lutetium(III)
    "PR":  "[Pr+3]",   # Praseodymium(III)
    # Alkaline earth / common biological metals
    "CA":  "[Ca+2]",   # Calcium(II)
    "MG":  "[Mg+2]",   # Magnesium(II)
    # Transition metals
    "ZN":  "[Zn+2]",   # Zinc(II)
    "FE":  "[Fe+2]",   # Iron(II)
    "FE3": "[Fe+3]",   # Iron(III)
    "MN":  "[Mn+2]",   # Manganese(II)
    "CU":  "[Cu+2]",   # Copper(II)
    "CO":  "[Co+2]",   # Cobalt(II)
    "NI":  "[Ni+2]",   # Nickel(II)
    "MO":  "[Mo+2]",   # Molybdenum(II)
}


def _metal_smiles(metal_ccd: str) -> str:
    """Return SMILES for a known ion; fall through CCD code for unknown ions."""
    return METAL_SMILES.get(metal_ccd.upper(), metal_ccd)


# ── Chai-1 FASTA builder ──────────────────────────────────────────────────────

def _build_fasta(protein_seq: str, metal_smiles: str, chain_name: str = "design") -> str:
    """
    Chai-1 FASTA format (chai-lab >= 0.4):

        >protein|name=<chain_name>
        <amino-acid sequence>
        >ligand|name=metal
        <SMILES or CCD code>

    The SMILES string goes in the sequence field for ligand entities; Chai-1
    passes it through RDKit which covers the full periodic table including lanthanides.
    """
    return (
        f">protein|name={chain_name}\n"
        f"{protein_seq}\n"
        f">ligand|name=metal\n"
        f"{metal_smiles}\n"
    )


# ── Output parsing ────────────────────────────────────────────────────────────

def _parse_scores(pred_dir: Path, protein_len: int) -> dict:
    """
    Extract confidence metrics from Chai-1 output files.

    Chai-1 writes per-model NPZ files (scores.model_idx_N.npz) and CIF structures
    (pred.model_idx_N.cif) into the prediction directory.  We take index 0 (the
    only model by default with num_diffn_timesteps as the single sample count).

    Key names tried in order to accommodate chai-lab version differences.
    pLDDT is sliced to the first protein_len entries to exclude the metal token.
    """
    out: dict = {}

    # ── Confidence NPZ ────────────────────────────────────────────────────────
    score_files = sorted(pred_dir.glob("scores.model_idx_*.npz"))
    if score_files:
        try:
            s = np.load(score_files[0])

            # pLDDT — per-residue; protein residues come first, then the metal token
            for key in ("plddt", "per_residue_plddt", "plddt_scores"):
                if key in s.files:
                    arr = np.asarray(s[key], dtype=float).flatten()
                    # Slice to protein only — metal token at the end is excluded
                    prot = arr[:protein_len] if len(arr) >= protein_len else arr
                    if prot.size:
                        v = float(prot.mean())
                        if 0.0 < v <= 1.0:   # normalise 0–1 → 0–100
                            v *= 100.0
                        out["chai_plddt"] = v
                    break

            # pTM (scalar or 1-element array)
            for key in ("ptm", "pTM", "per_chain_ptm"):
                if key in s.files:
                    out["chai_ptm"] = float(
                        np.asarray(s[key], dtype=float).flatten()[0])
                    break

            # Aggregate ranking score (Chai-1 specific composite; useful for HOF)
            for key in ("aggregate_score", "ranking_score", "score"):
                if key in s.files:
                    out["chai_score"] = float(
                        np.asarray(s[key], dtype=float).flatten()[0])
                    break

        except Exception as exc:
            print(f"  Warning: could not parse {score_files[0].name}: {exc}",
                  flush=True)

    # ── Predicted structure CIF ───────────────────────────────────────────────
    cif_files = sorted(pred_dir.glob("pred.model_idx_*.cif"))
    if not cif_files:
        cif_files = sorted(pred_dir.glob("*.cif"))   # fallback glob
    if cif_files:
        out["chai_cif"] = str(cif_files[0])

    return out


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Score metal-binding designs with Chai-1 (protein + metal ion).")
    ap.add_argument("--input", required=True,
                    help="Input CSV: design_id, sequence, metal_ccd")
    ap.add_argument("--out", required=True,
                    help="Output CSV: design_id, chai_plddt, chai_ptm[, chai_score, chai_cif]")
    ap.add_argument("--cif-dir", default=None,
                    help="Parent directory for per-design prediction subdirs "
                         "(default: chai1_structures/ alongside --out)")
    ap.add_argument("--device", default="cuda",
                    help="Torch device (cuda / cuda:0 / cpu)")
    ap.add_argument("--num-recycles", type=int, default=3,
                    help="Chai-1 trunk recycling steps (default 3)")
    ap.add_argument("--num-diffn-steps", type=int, default=200,
                    help="Chai-1 diffusion timesteps (default 200; lower = faster, noisier)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    # ── CUDA check ────────────────────────────────────────────────────────────
    if args.device.startswith("cuda"):
        import torch
        if not torch.cuda.is_available():
            sys.exit(
                "torch.cuda.is_available() is False.\n"
                "Fix: install a torch build matching the cluster CUDA version, "
                "or pass --device cpu."
            )
        print(f"CUDA OK: {torch.cuda.get_device_name(0)} | "
              f"torch {torch.__version__} (CUDA {torch.version.cuda})", flush=True)

    # ── Import chai-lab ───────────────────────────────────────────────────────
    try:
        import torch
        from chai_lab.chai1 import run_inference
    except ImportError as exc:
        sys.exit(
            f"chai-lab not installed in this env: {exc}\n"
            "  conda create -n chai1 python=3.10 -y\n"
            "  conda activate chai1 && pip install chai-lab"
        )

    device = torch.device(args.device)

    # ── Load designs ──────────────────────────────────────────────────────────
    designs = pd.read_csv(args.input)
    if designs.empty:
        sys.exit("No designs in input CSV.")

    out_path = Path(args.out)
    cif_dir  = Path(args.cif_dir) if args.cif_dir else out_path.parent / "chai1_structures"
    cif_dir.mkdir(parents=True, exist_ok=True)

    print(f"Scoring {len(designs)} designs with Chai-1 on {args.device} ...", flush=True)
    print("Model weights load on first call — this may take a minute.", flush=True)

    rows, n_ok = [], 0

    for i, d in enumerate(designs.itertuples(index=False), 1):
        seq       = getattr(d, "sequence",  None)
        metal_ccd = getattr(d, "metal_ccd", "CA")
        did       = str(d.design_id)

        if not isinstance(seq, str) or not seq:
            print(f"  [{i}/{len(designs)}] {did}: skipped (empty sequence)", flush=True)
            continue

        smiles        = _metal_smiles(metal_ccd)
        fasta_content = _build_fasta(seq, smiles, chain_name="design")

        # Chai-1 requires its own output dir per call (files are named by model index)
        pred_dir = cif_dir / did
        pred_dir.mkdir(parents=True, exist_ok=True)
        fasta_path = pred_dir / "input.fasta"
        fasta_path.write_text(fasta_content)

        try:
            # Full kwarg call (chai-lab >= 0.4)
            try:
                run_inference(
                    fasta_file=fasta_path,
                    output_dir=pred_dir,
                    num_trunk_recycles=args.num_recycles,
                    num_diffn_timesteps=args.num_diffn_steps,
                    seed=args.seed,
                    device=device,
                    use_esm_embeddings=True,
                )
            except TypeError:
                # Fallback: older chai-lab version with different signature
                run_inference(
                    fasta_file=fasta_path,
                    output_dir=pred_dir,
                )
        except Exception as exc:
            print(f"  [{i}/{len(designs)}] {did}: ERROR — {exc}", flush=True)
            continue

        metrics = _parse_scores(pred_dir, len(seq))
        if metrics:
            rows.append({"design_id": did, "metal_ccd": metal_ccd, **metrics})
            n_ok += 1
            plddt = metrics.get("chai_plddt", float("nan"))
            ptm   = metrics.get("chai_ptm",   float("nan"))
            print(f"  [{i}/{len(designs)}] {did} ({metal_ccd}, SMILES={smiles}): "
                  f"pLDDT={plddt:.1f}  pTM={ptm:.3f}", flush=True)
        else:
            print(f"  [{i}/{len(designs)}] {did}: no scores parsed — check {pred_dir}",
                  flush=True)

    if rows:
        pd.DataFrame(rows).to_csv(args.out, index=False)
        print(f"\nWrote {n_ok}/{len(designs)} Chai-1 scores to {args.out}", flush=True)
    else:
        print(
            "No scores produced.\n"
            "Checklist:\n"
            "  1. Is chai-lab installed in this Python env?\n"
            "  2. Does the FASTA format match your chai-lab version? "
            "(see github.com/chaidiscovery/chai-lab)\n"
            "  3. Are the SMILES strings valid? "
            "(verify with: python -c \"from rdkit import Chem; "
            "print(Chem.MolFromSmiles('[Eu+3]'))\")",
            flush=True,
        )


if __name__ == "__main__":
    main()
