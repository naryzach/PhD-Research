"""
Per-pair residue-residue contact matrices (deliverable G).

For a selected set of construct x target pairs, compute the matrix of minimum
heavy-atom distances between the TIMP3 reactive N-terminal edge (the ridge that
inserts into the catalytic cleft) and the target's catalytic zinc motif
(HExxHxxGxxH). This is the AlphaFoldServer_SA.ipynb "residue interaction matrix",
generalised across the co-fold and docking sources.

By default it renders, for every target, the best-scoring construct's AF3 co-fold
(the native-recovering source) as one panel of a gallery figure, and writes the
raw matrices as .csv alongside so any pair can be inspected numerically.

Run:  python Structural-Validation/analysis/contact_matrices.py
      python Structural-Validation/analysis/contact_matrices.py --source AF3_cofold --top 1
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
import config as C  # noqa: E402
import model_registry as MR  # noqa: E402
import structure_io as sio  # noqa: E402
from complex_metrics import assign_chains  # noqa: E402

EDGE_N = 14  # TIMP3 reactive N-terminal edge length to profile


def _min_dist_matrix(A, B, a_resids, b_resids):
    a_res = [r for r in A.residues if r.id[1] in set(a_resids)]
    b_res = [r for r in B.residues if r.id[1] in set(b_resids)]
    a_res.sort(key=lambda r: r.id[1])
    b_res.sort(key=lambda r: r.id[1])
    mat = np.full((len(a_res), len(b_res)), np.nan)
    for i, ra in enumerate(a_res):
        aa = np.array([at.coord for at in ra if at.element != "H"])
        for j, rb in enumerate(b_res):
            bb = np.array([at.coord for at in rb if at.element != "H"])
            if len(aa) and len(bb):
                mat[i, j] = np.linalg.norm(aa[:, None] - bb[None], axis=2).min()
    ylab = [f"{sio._one_letter(r.resname)}{r.id[1]}" for r in a_res]
    xlab = [f"{sio._one_letter(r.resname)}{r.id[1]}" for r in b_res]
    return mat, ylab, xlab


def matrix_for(model_path, construct_seq, target_seq):
    chains = sio.get_chains(model_path)
    if len(chains) < 2:
        return None
    cc, tc = assign_chains(chains, construct_seq, target_seq)
    edge = cc.resids[:EDGE_N]
    motif = sio.zinc_motif_resids(tc)
    if not motif:
        return None
    return _min_dist_matrix(cc, tc, edge, motif)


def _panel(ax, mat, ylab, xlab, title):
    im = ax.imshow(mat, cmap="viridis_r", vmin=2.0, vmax=12.0, aspect="auto")
    ax.set_xticks(range(len(xlab)))
    ax.set_xticklabels(xlab, rotation=90, fontsize=6)
    ax.set_yticks(range(len(ylab)))
    ax.set_yticklabels(ylab, fontsize=6)
    ax.set_title(title, fontsize=9, fontweight="bold")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if not np.isnan(mat[i, j]) and mat[i, j] < 5.0:
                ax.text(j, i, f"{mat[i, j]:.1f}", ha="center", va="center",
                        fontsize=5, color="white")
    return im


def select_pairs(cx, source, top):
    """Best `top` construct(s) per target for `source`, ranked by DockQ then composite."""
    ok = cx[(cx.status == "ok") & (cx.source == source)].copy()
    ok["dockq"] = pd.to_numeric(ok["dockq"], errors="coerce")
    ok["iface_composite"] = pd.to_numeric(ok["iface_composite"], errors="coerce")
    ok["rank"] = ok["dockq"].fillna(ok["iface_composite"] / 100.0)
    pairs = []
    for t in C.TARGET_ORDER:
        sub = ok[ok.target_id == t].sort_values("rank", ascending=False)
        for _, r in sub.head(top).iterrows():
            pairs.append((r.construct_id, r.target_id))
    return pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default="AF3_cofold")
    ap.add_argument("--top", type=int, default=1)
    args = ap.parse_args()

    cx = pd.read_csv(C.OUT_ANALYSIS / "master_complex_metrics.csv")
    reg = {r["id"]: r for r in MR.load_registry()}
    pairs = select_pairs(cx, args.source, args.top)

    out_dir = C.OUT_FIG / "contact_matrices"
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_dir = C.OUT_ANALYSIS / "contact_matrices"
    csv_dir.mkdir(parents=True, exist_ok=True)

    panels = []
    for cid, tid in pairs:
        mp = MR.af3_complex(cid, tid) if args.source == "AF3_cofold" else \
            MR.esmfold_complex(cid, tid)
        if mp is None:
            continue
        res = matrix_for(mp, reg[cid]["sequence"], reg[tid]["sequence"])
        if res is None:
            continue
        mat, ylab, xlab = res
        pd.DataFrame(mat, index=ylab, columns=xlab).to_csv(
            csv_dir / f"{cid}__{tid}_{args.source}.csv")
        panels.append((cid, tid, mat, ylab, xlab))

    if not panels:
        print("No contact matrices produced.")
        return

    n = len(panels)
    ncol = min(3, n)
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.6 * nrow),
                             squeeze=False)
    im = None
    for k, (cid, tid, mat, ylab, xlab) in enumerate(panels):
        ax = axes[k // ncol][k % ncol]
        im = _panel(ax, mat, ylab, xlab, f"{cid} × {tid}")
    for k in range(len(panels), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.suptitle(f"TIMP3 reactive edge × target zinc motif — min-distance (Å), "
                 f"{args.source}", fontsize=12, fontweight="bold")
    if im is not None:
        fig.colorbar(im, ax=axes, fraction=0.02, pad=0.02,
                     label="min heavy-atom distance (Å)")
    out_png = out_dir / f"contact_matrix_gallery_{args.source}.png"
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"{len(panels)} contact matrices -> {out_png.relative_to(C.REPO_ROOT)}")
    print(f"raw CSVs -> {csv_dir.relative_to(C.REPO_ROOT)}")


if __name__ == "__main__":
    main()
