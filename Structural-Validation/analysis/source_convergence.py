"""
TEST 2 — cross-source pose convergence (no ground truth required).

For every construct x target pair and every pair of model sources: superpose on the
TARGET chain (sequence-matched CA), apply that transform to the CONSTRUCT chain, and
measure construct CA RMSD. Averaged over all pairs this gives a 6x6 "do these two
methods put the binder in the same place?" matrix.

The key question is whether the FOUR HADDOCK tracks agree with EACH OTHER. If they
do not, their poses are arbitrary rather than a consistent alternative binding mode.

Caveat: AF3 and ESMFold2 agreeing with each other is only weak evidence of realism —
they share PDB training data and may simply be reproducing the same memorized
complex. Disagreement among HADDOCK tracks, however, needs no such assumption: a
method that cannot reproduce its own answer across equivalent inputs is unreliable
regardless of what any other method says.

Outputs:
  analysis/source_convergence_rmsd.csv
  figures/source_convergence.png

Run:  python Structural-Validation/analysis/source_convergence.py
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
import config as C  # noqa: E402
import metrics as M  # noqa: E402
import model_registry as MR  # noqa: E402
import structure_io as sio  # noqa: E402
from complex_metrics import _resolve, _source_specs, assign_chains  # noqa: E402

SRC_LABEL = {"AF3_cofold": "AF3 co-fold", "ESMFold2_cofold": "ESM co-fold",
             "HADDOCK:AF3xAF3": "HDK AF3×AF3",
             "HADDOCK:AF3xCrystal": "HDK AF3×Xtal",
             "HADDOCK:ESMFold2xCrystal": "HDK ESM×Xtal",
             "HADDOCK:ESMFold2xESMFold2": "HDK ESM×ESM"}


def pose_rmsd(ci, ti, cj, tj) -> float:
    """Construct CA RMSD after superposing source i onto source j via the target."""
    mt, rt = M._matched_ca(ti, tj)
    if len(mt) < 3:
        return float("nan")
    R, tr, _ = M.kabsch(mt, rt)
    mc, rc = M._matched_ca(ci, cj)
    if len(mc) < 1:
        return float("nan")
    moved = M.apply_transform(mc, R, tr)
    return float(np.sqrt(((moved - rc) ** 2).sum(1).mean()))


def main() -> None:
    C.OUT_ANALYSIS.mkdir(parents=True, exist_ok=True)
    C.OUT_FIG.mkdir(parents=True, exist_ok=True)
    reg = {r["id"]: r for r in MR.load_registry()}
    constructs = [r for r in reg.values() if r["kind"] == "construct"]
    targets = [r for r in reg.values() if r["kind"] == "target"]
    specs = _source_specs()
    sources = [s["source"] for s in specs]

    acc: dict[tuple[str, str], list[float]] = {}
    n_pairs = 0
    for c in constructs:
        for t in targets:
            models = {}
            for spec in specs:
                mp = _resolve(spec, c["id"], t["id"])
                if mp is None:
                    continue
                try:
                    chains = sio.get_chains(mp)
                    if len(chains) < 2:
                        continue
                    models[spec["source"]] = assign_chains(
                        chains, c["sequence"], t["sequence"])
                except Exception:
                    continue
            if len(models) < 2:
                continue
            n_pairs += 1
            for a, b in combinations(sources, 2):
                if a not in models or b not in models:
                    continue
                ca_, ta_ = models[a]
                cb_, tb_ = models[b]
                try:
                    r = pose_rmsd(ca_, ta_, cb_, tb_)
                except Exception:
                    r = float("nan")
                if r == r:
                    acc.setdefault((a, b), []).append(r)

    mat = pd.DataFrame(np.nan, index=sources, columns=sources, dtype=float)
    cnt = pd.DataFrame(0, index=sources, columns=sources, dtype=int)
    for s in sources:
        mat.loc[s, s] = 0.0
    for (a, b), vals in acc.items():
        m = float(np.mean(vals))
        mat.loc[a, b] = mat.loc[b, a] = round(m, 2)
        cnt.loc[a, b] = cnt.loc[b, a] = len(vals)

    out_csv = C.OUT_ANALYSIS / "source_convergence_rmsd.csv"
    mat.to_csv(out_csv)

    # ── figure ────────────────────────────────────────────────────────────
    lab = [SRC_LABEL.get(s, s) for s in sources]
    fig, ax = plt.subplots(figsize=(8.2, 6.6))
    im = ax.imshow(mat.values, cmap="magma_r", vmin=0, vmax=np.nanmax(mat.values))
    ax.set_xticks(range(len(lab))); ax.set_xticklabels(lab, rotation=35, ha="right")
    ax.set_yticks(range(len(lab))); ax.set_yticklabels(lab)
    for i in range(len(lab)):
        for j in range(len(lab)):
            v = mat.values[i, j]
            if v == v:
                ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=9,
                        color="white" if v > np.nanmax(mat.values) * .55 else "#222")
    ax.set_title("Cross-source pose convergence\nconstruct CA RMSD (Å) after "
                 "superposing on the target — lower = same binding mode",
                 fontweight="bold", fontsize=11)
    fig.colorbar(im, ax=ax, fraction=.046, label="mean construct RMSD (Å)")
    ax.grid(False)
    fig.tight_layout()
    out_png = C.OUT_FIG / "source_convergence.png"
    fig.savefig(out_png, dpi=140)
    plt.close(fig)

    # ── summary ───────────────────────────────────────────────────────────
    hdk = [s for s in sources if s.startswith("HADDOCK")]
    cof = [s for s in sources if s.endswith("_cofold")]

    def block(rows, cols, exclude_diag=True):
        vals = []
        for a in rows:
            for b in cols:
                if exclude_diag and a == b:
                    continue
                v = mat.loc[a, b]
                if v == v:
                    vals.append(v)
        return float(np.mean(vals)) if vals else float("nan")

    print(f"source_convergence: {n_pairs} construct×target pairs compared\n")
    print(mat.to_string())
    print("\nBlock means (construct CA RMSD, Å):")
    print(f"  co-fold  vs co-fold   : {block(cof, cof):.2f}   "
          f"(AF3 vs ESMFold2 — weak evidence, shared training data)")
    print(f"  HADDOCK  vs HADDOCK   : {block(hdk, hdk):.2f}   "
          f"(do the 4 tracks agree with each other?)")
    print(f"  HADDOCK  vs co-fold   : {block(hdk, cof, False):.2f}")
    print(f"\n-> {out_csv.relative_to(C.REPO_ROOT)}")
    print(f"-> {out_png.relative_to(C.REPO_ROOT)}")


if __name__ == "__main__":
    main()
