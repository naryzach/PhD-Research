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
    # On a Blackwell GPU (RTX 50-series, sm_120 — e.g. RTX 5070 Ti) the torch
    # 2.6.0+cu124 wheel that chai-lab pulls in has no kernels for that
    # architecture and fails at inference time with "CUDA error: no kernel
    # image is available for execution on the device". Fix by overriding
    # torch past chai-lab's declared <2.7 ceiling (verified working with
    # chai-lab 0.6.1 despite the pin):
    #   pip install torch==2.7.1 --index-url https://download.pytorch.org/whl/cu128
    # See Chai1Env.yml (this directory) for a full known-working env export.

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
import shutil
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
    "Y":   "[Y+3]",    # Yttrium(III)      — in test_metals.py rare-earth list
    "TM":  "[Tm+3]",   # Thulium(III)      — in test_metals.py rare-earth list
    "PM":  "[Pm+3]",   # Promethium(III)   — in test_metals.py rare-earth list
    "SC":  "[Sc+3]",   # Scandium(III)     — in test_metals.py rare-earth list
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


class UnknownMetalError(ValueError):
    pass


def _metal_smiles(metal_ccd: str) -> str:
    """Return SMILES for a known ion.

    This used to fall through to the raw CCD code string for unknown ions —
    e.g. bare "SC" (Scandium's CCD code) is NOT the SMILES for the Sc3+ ion,
    it's valid SMILES for an unrelated S-C (thiol-ish) two-atom fragment, and
    RDKit parses it without complaint. Chai-1 would then silently score a
    design against the wrong ligand entirely: no error anywhere, just a
    confident-looking pLDDT/pTM for the wrong chemistry. Some other unlisted
    codes (bare element symbols outside SMILES' organic subset, e.g. "TM",
    "PM") do at least fail to parse — but that's inconsistent luck, not a
    safeguard. Failing loudly here for any ion missing from METAL_SMILES is
    the fix; add the ion's SMILES to the table above instead of relying on
    the fallback.
    """
    smiles = METAL_SMILES.get(metal_ccd.upper())
    if smiles is None:
        raise UnknownMetalError(
            f"No SMILES entry for metal_ccd={metal_ccd!r} in METAL_SMILES. "
            f"Add it to score_with_chai1.py rather than guessing — the CCD "
            f"code itself is not valid/correct SMILES in general."
        )
    return smiles


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
    top-ranked / first diffusion sample).

    IMPORTANT (verified against chai-lab 0.6.1 source, chai_lab/ranking/rank.py
    get_scores()): the scores NPZ does NOT contain pLDDT under any key — only
    "aggregate_score", "ptm", "iptm", "per_chain_ptm", "per_chain_pair_iptm",
    "has_inter_chain_clashes", "chain_chain_clashes". pLDDT is only persisted as
    the per-atom B-factor column of the predicted CIF (already scaled 0-100 by
    chai-lab itself, see chai_lab/chai1.py run_folding_on_context), so it must be
    read from there instead. A key-based npz lookup for pLDDT would silently find
    nothing and leave chai_plddt unset without raising — exactly the "Chai-1 runs
    but downstream gets bad/missing data" failure mode this pipeline hit before.
    """
    out: dict = {}

    # ── Confidence NPZ (pTM / aggregate score only — no pLDDT here) ───────────
    score_files = sorted(pred_dir.glob("scores.model_idx_*.npz"))
    if score_files:
        try:
            s = np.load(score_files[0])

            for key in ("ptm", "pTM", "per_chain_ptm"):
                if key in s.files:
                    out["chai_ptm"] = float(
                        np.asarray(s[key], dtype=float).flatten()[0])
                    break

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

        # ── pLDDT from CIF B-factors ───────────────────────────────────────────
        # Chai-1 writes per-atom pLDDT (0-100 scale) into the B-factor column.
        # Average over protein-chain atoms only; the metal ligand is written as
        # its own hetero residue and would otherwise skew the mean.
        try:
            plddt_vals = _plddt_from_cif(cif_files[0])
            if plddt_vals.size:
                out["chai_plddt"] = float(plddt_vals.mean())
        except Exception as exc:
            print(f"  Warning: could not parse pLDDT from {cif_files[0].name}: {exc}",
                  flush=True)

    return out


def _plddt_from_cif(cif_path: Path) -> np.ndarray:
    """Read per-atom B-factors (== pLDDT, 0-100) for the protein chain of a
    Chai-1 output CIF, using gemmi (already a chai-lab dependency, avoids
    requiring biotite in this env)."""
    import gemmi
    structure = gemmi.read_structure(str(cif_path))
    model = structure[0]
    vals = []
    for chain in model:
        for residue in chain:
            # The metal ligand is written as its own single-atom HETATM residue
            # named "LIG" (chai_lab.data.residue_constants.new_ligand_residue_name,
            # written with het=True in chai_lab.data.io.cif_utils.save_to_cif);
            # skip it so it doesn't pull the protein pLDDT average off.
            if residue.het_flag == "H" or residue.name == "LIG":
                continue
            for atom in residue:
                vals.append(atom.b_iso)
    return np.asarray(vals, dtype=float)


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

        try:
            smiles = _metal_smiles(metal_ccd)
        except UnknownMetalError as exc:
            print(f"  [{i}/{len(designs)}] {did}: ERROR — {exc}", flush=True)
            continue
        fasta_content = _build_fasta(seq, smiles, chain_name="design")

        # run_inference() asserts that output_dir is EITHER absent OR empty
        # (chai_lab/chai1.py: "assert not any(output_dir.iterdir())"). Writing
        # the input FASTA into pred_dir before calling run_inference — as this
        # script used to do — makes that assertion fail on every single design,
        # every time: run_inference raises immediately, before any model call,
        # the bare `except Exception` below swallows it, and the design is
        # silently skipped (chai_plddt/chai_ptm end up NaN downstream with no
        # loud failure). So the FASTA goes in its own directory, and pred_dir
        # itself is left absent/empty for run_inference to populate. Any
        # leftover pred_dir from a previous run is cleared first so reruns
        # don't hit the same assertion.
        fasta_dir = cif_dir / "_fasta_inputs"
        fasta_dir.mkdir(parents=True, exist_ok=True)
        fasta_path = fasta_dir / f"{did}.fasta"
        fasta_path.write_text(fasta_content)

        pred_dir = cif_dir / did
        if pred_dir.exists():
            shutil.rmtree(pred_dir)

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
        missing = [k for k in ("chai_plddt", "chai_ptm") if k not in metrics]
        if metrics and not missing:
            rows.append({"design_id": did, "metal_ccd": metal_ccd, **metrics})
            n_ok += 1
            plddt = metrics.get("chai_plddt", float("nan"))
            ptm   = metrics.get("chai_ptm",   float("nan"))
            print(f"  [{i}/{len(designs)}] {did} ({metal_ccd}, SMILES={smiles}): "
                  f"pLDDT={plddt:.1f}  pTM={ptm:.3f}", flush=True)
        elif metrics:
            # Chai-1 "succeeded" (files were written) but one or more expected
            # metrics couldn't be parsed out of them. Still write the row (with
            # NaN for the missing field) so it's visible in the CSV, but flag it
            # loudly here rather than let it look like a clean success — this is
            # the exact "ran but produced bad/partial data silently" failure mode.
            rows.append({"design_id": did, "metal_ccd": metal_ccd, **metrics})
            print(f"  [{i}/{len(designs)}] {did}: *** PARTIAL — missing {missing} "
                  f"(check {pred_dir}) ***", flush=True)
        else:
            print(f"  [{i}/{len(designs)}] {did}: no scores parsed — check {pred_dir}",
                  flush=True)

    if rows:
        pd.DataFrame(rows).to_csv(args.out, index=False)
        print(f"\nWrote {len(rows)}/{len(designs)} Chai-1 rows to {args.out} "
              f"({n_ok} complete, {len(rows) - n_ok} partial)", flush=True)
        if n_ok < len(designs):
            # Non-zero exit even though a CSV was written, so the caller
            # (run_pipeline.py, which only inspects stderr/returncode on
            # failure) can't mistake "some rows are missing/partial" for a
            # clean run just because the subprocess "succeeded".
            print(f"*** {len(designs) - n_ok} / {len(designs)} designs incomplete "
                  f"or missing — see per-design lines above ***", file=sys.stderr,
                  flush=True)
            sys.exit(2)
    else:
        msg = (
            "No scores produced.\n"
            "Checklist:\n"
            "  1. Is chai-lab installed in this Python env?\n"
            "  2. Does the FASTA format match your chai-lab version? "
            "(see github.com/chaidiscovery/chai-lab)\n"
            "  3. Are the SMILES strings valid? "
            "(verify with: python -c \"from rdkit import Chem; "
            "print(Chem.MolFromSmiles('[Eu+3]'))\")"
        )
        print(msg, flush=True)
        print(msg, file=sys.stderr, flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
