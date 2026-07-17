"""
loop_probe_analysis.py

GPU-free analysis + plotting for the generative-pipeline loop probe.

This module contains everything needed to turn a bag of designed loop
sub-sequences (all the same length, one per generated design) into
per-position amino-acid-frequency heatmaps — both the raw 20-AA view and
several biochemical-group views (charge, size, chemical type, hydrophobicity,
polarity, aromaticity).

It deliberately has NO dependency on torch / rfd3 / mpnn / atomworks, so the
heatmaps can be (re)generated on any laptop from the CSVs written by the GPU
run.  `loop_probe.py` imports this module after generation; you can also run it
standalone to rebuild figures:

    python loop_probe_analysis.py --counts path/to/position_counts_AB.csv
    python loop_probe_analysis.py --sequences path/to/sequences.csv --loop AB

Author: probe tooling for the TIMP3-scaffold binder pipeline.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ── Amino-acid vocabulary ──────────────────────────────────────────────────────
# Fixed row order used everywhere so heatmaps from different runs line up.
AA20 = list("ACDEFGHIKLMNPQRSTVWY")
AA_SET = set(AA20)

# ── Biochemical grouping schemes ────────────────────────────────────────────────
# Each scheme is an ORDERED dict {group_label: "AAs"}.  Every scheme partitions
# the 20 amino acids (each AA appears in exactly one group) so the grouped
# frequencies at a position sum to the same total as the raw frequencies.  The
# order of groups is chosen so the heatmap reads along a meaningful axis
# (e.g. negative → neutral → positive for charge).
PROPERTY_SCHEMES: dict[str, dict[str, str]] = {
    # Net side-chain charge at physiological pH (H treated as weakly positive).
    "charge": {
        "Negative (D,E)":       "DE",
        "Polar/Neutral":        "STNQCGAPVILMFWY",
        "Positive (K,R,H)":     "KRH",
    },
    # Kyte–Doolittle-style tripartite split.
    "hydrophobicity": {
        "Hydrophobic":          "AVLIMFWC",
        "Neutral":              "GSTPNQYH",
        "Hydrophilic (charged)": "DEKR",
    },
    # Side-chain volume (van der Waals) binned into three sizes.
    "size": {
        "Small (G,A,S,C,D,P,N,T)": "GASCDPNT",
        "Medium (E,V,Q,H)":        "EVQH",
        "Large (M,I,L,K,R,F,Y,W)": "MILKRFYW",
    },
    # Chemical class of the side chain.
    "type": {
        "Aliphatic (G,A,V,L,I)": "GAVLI",
        "Hydroxyl (S,T)":        "ST",
        "Sulfur (C,M)":          "CM",
        "Aromatic (F,W,Y)":      "FWY",
        "Amide (N,Q)":           "NQ",
        "Acidic (D,E)":          "DE",
        "Basic (K,R,H)":         "KRH",
        "Imino (P)":             "P",
    },
    # Polar vs non-polar (classic two-way + charged split).
    "polarity": {
        "Nonpolar":             "AVLIMFWGP",
        "Polar uncharged":      "STCNQY",
        "Charged":              "DEKRH",
    },
    # Aromatic ring present or not (H counted as aromatic).
    "aromaticity": {
        "Aromatic (F,W,Y,H)":   "FWYH",
        "Non-aromatic":         "ACDEGIKLMNPQRSTV",
    },
}


def _validate_schemes() -> None:
    """Fail loudly if any scheme does not partition all 20 amino acids."""
    for name, groups in PROPERTY_SCHEMES.items():
        seen: dict[str, str] = {}
        for label, letters in groups.items():
            for a in letters:
                if a not in AA_SET:
                    raise ValueError(f"scheme {name!r}: '{a}' is not a standard AA")
                if a in seen:
                    raise ValueError(
                        f"scheme {name!r}: '{a}' in both {seen[a]!r} and {label!r}")
                seen[a] = label
        missing = AA_SET - set(seen)
        if missing:
            raise ValueError(f"scheme {name!r} omits {sorted(missing)}")


_validate_schemes()


# ── Counting ────────────────────────────────────────────────────────────────────
def count_positions(loop_seqs: list[str], length: int | None = None) -> pd.DataFrame:
    """
    Tally amino acids per position across a set of equal-length loop sequences.

    Parameters
    ----------
    loop_seqs : list[str]
        Designed loop sub-sequences (flanks already stripped).  Sequences whose
        length does not match `length` are dropped (and counted as parse fails
        by the caller if it cares).
    length : int | None
        Expected loop length.  If None, inferred as the modal sequence length.

    Returns
    -------
    DataFrame indexed by the 20 amino acids (AA20 order), columns
    ``pos01 … posLL`` (1-indexed), values = integer counts.  All-zero rows for
    unused amino acids are retained so every heatmap has the full 20-AA axis.
    """
    seqs = [s for s in loop_seqs if s and set(s) <= AA_SET]
    if length is None:
        if not seqs:
            length = 0
        else:
            length = int(pd.Series([len(s) for s in seqs]).mode().iloc[0])
    seqs = [s for s in seqs if len(s) == length]

    cols = [f"pos{p:02d}" for p in range(1, length + 1)]
    counts = pd.DataFrame(0, index=AA20, columns=cols, dtype=int)
    for s in seqs:
        for i, a in enumerate(s):
            counts.iat[AA20.index(a), i] += 1
    counts.attrs["n_sequences"] = len(seqs)
    return counts


def to_frequency(counts: pd.DataFrame) -> pd.DataFrame:
    """Convert a count matrix to per-position frequency (columns sum to 1)."""
    col_tot = counts.sum(axis=0).replace(0, np.nan)
    freq = counts.div(col_tot, axis=1).fillna(0.0)
    freq.attrs.update(counts.attrs)
    return freq


def group_counts(counts: pd.DataFrame, scheme: str) -> pd.DataFrame:
    """Collapse a 20-AA count/frequency matrix into biochemical-group rows."""
    groups = PROPERTY_SCHEMES[scheme]
    rows = {label: counts.loc[[a for a in letters], :].sum(axis=0)
            for label, letters in groups.items()}
    out = pd.DataFrame(rows).T
    out.attrs.update(counts.attrs)
    return out


# ── Plotting ────────────────────────────────────────────────────────────────────
_CMAP = LinearSegmentedColormap.from_list(
    "probe_blues", ["#f7fbff", "#6baed6", "#08306b"])


def heatmap(matrix: pd.DataFrame, title: str, out_path: Path,
            cbar_label: str = "frequency", annotate: bool = True,
            vmax: float | None = None) -> Path:
    """
    Render a rows×positions heatmap to `out_path` (PNG).

    Works for both the 20-AA matrix and the collapsed group matrices — the row
    labels come straight from the DataFrame index.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_rows, n_cols = matrix.shape
    fig_w = max(4.0, 0.55 * n_cols + 2.2)
    fig_h = max(2.4, 0.42 * n_rows + 1.4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=150)

    data = matrix.values.astype(float)
    im = ax.imshow(data, aspect="auto", cmap=_CMAP,
                   vmin=0.0, vmax=vmax if vmax is not None else max(data.max(), 1e-9))

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([c.replace("pos", "") for c in matrix.columns], fontsize=8)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(matrix.index, fontsize=8, family="monospace")
    ax.set_xlabel("loop position", fontsize=9)

    # Annotate cells when the grid is small enough to stay readable.
    if annotate and n_rows * n_cols <= 20 * 22:
        thresh = (vmax if vmax is not None else data.max()) * 0.55
        for i in range(n_rows):
            for j in range(n_cols):
                v = data[i, j]
                if v <= 0:
                    continue
                ax.text(j, i, f"{v:.2f}" if v < 1 else f"{int(v)}",
                        ha="center", va="center", fontsize=6,
                        color="white" if v > thresh else "#222")

    n = matrix.attrs.get("n_sequences")
    subtitle = f"  (n={n} designs)" if n is not None else ""
    ax.set_title(title + subtitle, fontsize=10)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label(cbar_label, fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


def render_all_heatmaps(counts: pd.DataFrame, out_dir: Path, stem: str,
                        title_prefix: str = "") -> list[Path]:
    """
    From one loop's count matrix, write the full set of figures + CSVs:
      * raw 20-AA frequency heatmap
      * one grouped heatmap per PROPERTY_SCHEMES entry
    Returns the list of PNG paths written.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    freq = to_frequency(counts)
    counts.to_csv(out_dir / f"position_counts_{stem}.csv")
    freq.to_csv(out_dir / f"position_freq_{stem}.csv")

    written.append(heatmap(
        freq, f"{title_prefix}{stem}: per-position AA frequency",
        out_dir / f"heatmap_{stem}_AA.png", cbar_label="frequency",
        annotate=True, vmax=1.0))

    for scheme in PROPERTY_SCHEMES:
        gfreq = group_counts(freq, scheme)
        gfreq.to_csv(out_dir / f"position_freq_{stem}_{scheme}.csv")
        written.append(heatmap(
            gfreq, f"{title_prefix}{stem}: {scheme}",
            out_dir / f"heatmap_{stem}_{scheme}.png",
            cbar_label="group frequency", annotate=True, vmax=1.0))

    return written


# ── Cross-length aggregation (for the length sweep) ─────────────────────────────
def mean_group_frequency(counts: pd.DataFrame, scheme: str) -> pd.Series:
    """
    Position-averaged frequency of each biochemical group across a loop.
    Collapses the positional axis so lengths with different position counts are
    comparable.  Returns a Series indexed by the scheme's group labels.
    """
    gfreq = group_counts(to_frequency(counts), scheme)
    return gfreq.mean(axis=1) if gfreq.shape[1] else gfreq.sum(axis=1)


def length_trend_matrix(counts_by_length: dict[int, pd.DataFrame],
                        scheme: str) -> pd.DataFrame:
    """
    group × length matrix of position-averaged group frequency — the compact
    "how does composition shift as the loop grows" view for one scheme.
    """
    cols = {}
    for L in sorted(counts_by_length):
        cols[f"L{L:02d}"] = mean_group_frequency(counts_by_length[L], scheme)
    out = pd.DataFrame(cols)
    out = out.reindex(list(PROPERTY_SCHEMES[scheme]))  # stable group order
    return out


def length_montage(freq_by_length: dict[int, pd.DataFrame], title: str,
                   out_path: Path, cbar_label: str = "frequency") -> Path:
    """
    One figure stacking the raw 20-AA frequency heatmaps for every swept length
    (one subplot per length, shared 20-AA y-axis).  Widths differ because longer
    loops have more positions, so each subplot draws its own position axis.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lengths = sorted(freq_by_length)
    n = len(lengths)
    max_cols = max((freq_by_length[L].shape[1] for L in lengths), default=1)

    fig_h = max(3.0, 0.34 * len(AA20) * 1.0 + 0.8)
    fig, axes = plt.subplots(
        1, n, figsize=(max(3.0, 0.5 * max_cols) * n + 1.5, fig_h),
        dpi=150, squeeze=False)
    axes = axes[0]
    im = None
    for ax, L in zip(axes, lengths):
        m = freq_by_length[L]
        im = ax.imshow(m.values.astype(float), aspect="auto", cmap=_CMAP,
                       vmin=0.0, vmax=1.0)
        ax.set_xticks(range(m.shape[1]))
        ax.set_xticklabels([c.replace("pos", "") for c in m.columns], fontsize=7)
        ax.set_xlabel(f"L={L}", fontsize=9)
        n_seq = m.attrs.get("n_sequences")
        if n_seq is not None:
            ax.set_title(f"n={n_seq}", fontsize=8)
        ax.set_yticks(range(len(AA20)))
        ax.set_yticklabels(AA20 if ax is axes[0] else [], fontsize=7,
                           family="monospace")
    fig.suptitle(title, fontsize=11)
    if im is not None:
        cbar = fig.colorbar(im, ax=list(axes), fraction=0.015, pad=0.02)
        cbar.set_label(cbar_label, fontsize=8)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ── Standalone CLI (rebuild figures without a GPU) ───────────────────────────────
def _main() -> None:
    ap = argparse.ArgumentParser(description="Rebuild loop-probe heatmaps from saved data.")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--counts", help="position_counts_<loop>.csv (AA x position)")
    src.add_argument("--sequences", help="sequences.csv with a loop_<name>_seq column")
    ap.add_argument("--loop", help="loop name (stem); required with --sequences")
    ap.add_argument("--out-dir", default=".", help="output directory for figures")
    ap.add_argument("--title-prefix", default="")
    args = ap.parse_args()

    if args.counts:
        counts = pd.read_csv(args.counts, index_col=0)
        counts = counts.reindex(AA20).fillna(0).astype(int)
        stem = Path(args.counts).stem.replace("position_counts_", "")
    else:
        if not args.loop:
            ap.error("--loop is required with --sequences")
        df = pd.read_csv(args.sequences)
        col = f"loop_{args.loop}_seq"
        if col not in df.columns:
            ap.error(f"{col} not found in {args.sequences}")
        seqs = [s for s in df[col].astype(str).tolist() if s and s != "MISSING"]
        counts = count_positions(seqs)
        counts.attrs["n_sequences"] = len(seqs)
        stem = args.loop

    paths = render_all_heatmaps(counts, Path(args.out_dir), stem,
                                title_prefix=args.title_prefix)
    print(f"Wrote {len(paths)} figures to {args.out_dir}")


if __name__ == "__main__":
    _main()
