"""
cross_docking_analysis.py

Target-specificity (cross-docking) analysis for TIMP3 binders, using ESMFold2.

For each top design we fold its binder sequence against EVERY target (its intended
target and all "alien" targets) as a two-chain complex with ESMFold2, recording
ipTM / pTM / pLDDT per pair.  analyze_results.py turns these into specificity gaps
(intended ipTM - mean alien ipTM), and the dashboard's Specificity tab visualizes
the intended-vs-folded matrix.

ESMFold2 replaced RF3 here (see run_pipeline.py / filter_methods.md): it folds
directly from sequence, so cross-docking no longer needs to splice atom arrays — it
just pairs the design's binder sequence with each target's chain-B sequence.

Output: <OUT_BASE_DIR>/cross_docking_metrics.csv
        columns: design_id, intended_target, folded_target, plddt, ptm, iptm
Per-pair scores are also cached to cross_docking/<target>/<id>_scores.json so a
re-run resumes instead of re-folding, and predicted complex CIFs are written for
inspection when SAVE_ESM_STRUCTS is True.
"""

import os
import json
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm
import biotite.structure.io.pdbx as pdbx
from biotite.sequence import ProteinSequence

# --- CONFIGURATION ---
OUT_BASE_DIR = "../Local/TIMP-Dashboard_output"
ADVANCED_METRICS_PATH = os.path.join(OUT_BASE_DIR, "advanced_metrics.csv")
CROSS_DOCK_OUT_DIR = os.path.join(OUT_BASE_DIR, "cross_docking")
DATA_DIR = "../Data/TIMP_Complexes/AlphaFold_CIF"

TARGET_PDBS = {
    "ADAM10": "TIMP3_vs_ADAM10_AF.cif",
    "ADAM17": "TIMP3_vs_ADAM17_AF.cif",
    "MMP2":   "TIMP3_vs_MMP2_AF.cif",
    "MMP3":   "TIMP3_vs_MMP3_AF.cif",
    "MMP9":   "TIMP3_vs_MMP9_AF.cif",
    "MMP10":  "TIMP3_vs_MMP10_AF.cif",
}

TOP_K = 10                # Number of top designs PER Target-Combo to cross-dock
SAVE_ESM_STRUCTS = True   # write predicted cross-docked complex CIFs (≈ free)

# In our pipeline and the reference CIFs: chain A = TIMP3 binder, chain B = target.
BINDER_CHAIN = "A"
TARGET_CHAIN = "B"


# ── ESMFold2 helpers (inlined; mirror run_pipeline.py) ───────────────────────
def _esm_load_model(device: str = "cuda"):
    """Load ESMFold2 once and return (model, input_builder)."""
    from esm.models.esmfold2 import ESMFold2InputBuilder
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
    model = ESMFold2Model.from_pretrained("biohub/ESMFold2").to(device).eval()
    return model, ESMFold2InputBuilder()


def _esm_save_structure(res, stub: str):
    """Write the predicted complex to {stub}.cif (or .pdb). Returns path or None."""
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


def _esm_predict_complex(model, builder, binder_seq: str, target_seq: str,
                         binder_len: int, seed: int = 0, cif_stub: str = None) -> dict:
    """
    Fold a two-chain complex and return {esm_iptm, esm_ptm, esm_plddt[, esm_cif]}.
    esm_plddt is the binder-chain mean on the 0-1 scale used across this dashboard
    pipeline.  If cif_stub is given the predicted structure is written at no extra
    GPU cost.
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

    plddt = np.nan
    plddt_arr = _get(res, "plddt", "plddts")
    if plddt_arr is not None:
        a = np.asarray(plddt_arr, dtype=float).flatten()
        if a.size:
            if np.nanmax(a) > 1.0:           # model emitted a 0-100 scale
                a = a / 100.0
            binder_part = a[:binder_len] if a.size >= binder_len else a
            plddt = float(np.mean(binder_part))

    def _f(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return np.nan

    out = {"esm_iptm": _f(iptm), "esm_ptm": _f(ptm), "esm_plddt": plddt}
    if cif_stub:
        saved = _esm_save_structure(res, cif_stub)
        if saved:
            out["esm_cif"] = saved
    return out


def get_sequence_from_array(atom_array, chain_id: str = "A") -> str:
    """One-letter sequence from the CA atoms of a chain."""
    mask = (atom_array.chain_id == chain_id) & (atom_array.atom_name == "CA")
    ca = atom_array[mask]
    if len(ca) == 0:
        return ""
    ca = ca[np.argsort(ca.res_id)]
    letters = []
    for rn in ca.res_name:
        try:
            letters.append(ProteinSequence.convert_letter_3to1(rn))
        except Exception:
            letters.append("X")
    return "".join(letters)


def load_target_sequence(target_name: str):
    """Return the target (protease) chain-B sequence from the reference CIF."""
    pdb_path = os.path.join(DATA_DIR, TARGET_PDBS[target_name])
    if not os.path.exists(pdb_path):
        return None
    arr = pdbx.get_structure(pdbx.CIFFile.read(pdb_path), model=1)
    return get_sequence_from_array(arr, chain_id=TARGET_CHAIN)


def run_cross_docking():
    os.makedirs(CROSS_DOCK_OUT_DIR, exist_ok=True)

    if not os.path.exists(ADVANCED_METRICS_PATH):
        print(f"Metrics not found at {ADVANCED_METRICS_PATH}")
        return

    df = pd.read_csv(ADVANCED_METRICS_PATH)

    # Pre-load each target's chain-B sequence once (cheap; avoids re-reading per design).
    target_seqs = {}
    for tname in TARGET_PDBS:
        seq = load_target_sequence(tname)
        if seq:
            target_seqs[tname] = seq
        else:
            print(f"Warning: could not load target sequence for {tname}; skipping it.")
    if not target_seqs:
        print("No target sequences available; aborting cross-docking.")
        return

    # Selection: top designs per (target, loop_combo) by binding probability.
    top_designs = []
    for (target, combo), group in df.groupby(['target', 'loop_combo']):
        clean_target = target.replace("TIMP3_vs_", "").replace("_AF", "").upper()
        top = group.sort_values("probability_of_binding_score", ascending=False).head(TOP_K)
        for _, row in top.iterrows():
            top_designs.append({
                "clean_target": clean_target,
                "full_seq": row['full_seq'],
                "design_id": row['design_id'],
            })

    print(f"Selected {len(top_designs)} designs for cross-docking against "
          f"{len(target_seqs)} targets.")

    # Load ESMFold2 once for the whole sweep.
    device = "cuda" if torch.cuda.is_available() else "cpu"
    try:
        print(f"Loading ESMFold2 on {device} ...", flush=True)
        model, builder = _esm_load_model(device)
    except ImportError as e:
        print(
            f"ESMFold2 not available ({e}).\n"
            "  pip install 'esm @ git+https://github.com/Biohub/esm.git@main' transformers\n"
            "  Cross-docking specificity metrics cannot be computed without it."
        )
        return

    cross_dock_results = []
    try:
        for design in tqdm(top_designs, desc="Designs"):
            orig_id    = design['design_id']
            intended   = design['clean_target']
            binder_seq = design['full_seq']
            if not isinstance(binder_seq, str) or not binder_seq:
                continue

            # Fold against every target (including the intended one) so the
            # specificity gap is an apples-to-apples ESMFold2 comparison.
            for target_name, target_seq in tqdm(target_seqs.items(),
                                                desc=f"Docking {orig_id}", leave=False):
                new_id     = f"{orig_id}_vs_{target_name}_XD"
                out_folder = os.path.join(CROSS_DOCK_OUT_DIR, target_name)
                os.makedirs(out_folder, exist_ok=True)
                scores_path = os.path.join(out_folder, f"{new_id}_scores.json")

                # Resume: reuse a previously computed score if present.
                if os.path.exists(scores_path):
                    try:
                        with open(scores_path) as f:
                            cross_dock_results.append(json.load(f))
                        continue
                    except Exception:
                        pass  # corrupt/old cache — fall through and recompute

                try:
                    stub = os.path.join(out_folder, new_id) if SAVE_ESM_STRUCTS else None
                    metrics = _esm_predict_complex(
                        model, builder, binder_seq, target_seq, len(binder_seq),
                        cif_stub=stub,
                    )
                    rec = {
                        "design_id":       orig_id,
                        "intended_target": intended,
                        "folded_target":   target_name,
                        "plddt":           metrics["esm_plddt"],
                        "ptm":             metrics["esm_ptm"],
                        "iptm":            metrics["esm_iptm"],
                    }
                    with open(scores_path, "w") as f:
                        json.dump(rec, f)
                    cross_dock_results.append(rec)
                except Exception as e:
                    print(f"Error folding {new_id}: {e}")
                finally:
                    torch.cuda.empty_cache()
    finally:
        del model, builder
        torch.cuda.empty_cache()

    # Save results
    if cross_dock_results:
        xd_df = pd.DataFrame(cross_dock_results)
        out_csv = os.path.join(OUT_BASE_DIR, "cross_docking_metrics.csv")
        xd_df.to_csv(out_csv, index=False)
        print(f"Cross-docking complete. Saved to {out_csv}")
    else:
        print("Cross-docking produced no results.")


if __name__ == "__main__":
    run_cross_docking()
