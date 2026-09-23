"""
loop_probe_hotspots.py

Cross-cutting analysis of a completed loop-probe sweep.  GPU-free — reads the
per-run position_counts CSVs a sweep already wrote and produces:

  1. HOTSPOTS (per target x loop, at each loop's native length)
     Per position: dominant amino acid + its frequency, Shannon information
     content (bits, sequence-logo style), and the dominant biochemical group in
     each scheme (charge / hydrophobicity / size / type / polarity / aromaticity).
     A position is a "hotspot" when it is strongly constrained (high info content).

  2. CROSS-TARGET divergence (per loop x position)
     Jensen-Shannon divergence of the per-position AA distribution across the
     four targets (normalised 0-1).  High = the position is designed differently
     depending on the target (a specificity determinant); low = shared.

  3. LENGTH-TREND synthesis (per loop x scheme x group)
     Spearman rho of position-averaged group frequency vs loop length, from the
     length_trend CSVs — which compositions shift as a loop grows/shrinks.

  4. METAL-COORDINATION propensity (targets are Zn proteases)
     Per-position frequency of His/Cys/Asp/Glu — candidate zinc-ligand positions.

Outputs: <sweep>/analysis/*.csv + *.png, and a printed text summary.

    python loop_probe_hotspots.py <sweep_dir>
"""

from __future__ import annotations

import sys
import glob
import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = Path(__file__).parent.resolve()
sys.path.insert(0, str(_HERE))
import loop_probe_analysis as lpa   # AA20, PROPERTY_SCHEMES, group/freq helpers

METAL_LIGANDS = set("HCDE")   # canonical Zn-coordinating side chains
LOG2_20 = np.log2(20)


# ── metrics ─────────────────────────────────────────────────────────────────────
def freq_matrix(counts: pd.DataFrame) -> pd.DataFrame:
    return lpa.to_frequency(counts)


def info_bits(freq_col: np.ndarray) -> float:
    """Sequence-logo information content of one position, in bits (0..log2 20)."""
    p = freq_col[freq_col > 0]
    H = -np.sum(p * np.log2(p))
    return float(LOG2_20 - H)


def js_divergence(dists: list[np.ndarray]) -> float:
    """Jensen-Shannon divergence of N distributions, normalised to 0..1 (bits/log2 N)."""
    dists = [d / d.sum() for d in dists if d.sum() > 0]
    if len(dists) < 2:
        return np.nan
    M = np.mean(dists, axis=0)
    def H(p):
        p = p[p > 0]
        return -np.sum(p * np.log2(p))
    jsd = H(M) - np.mean([H(d) for d in dists])
    return float(jsd / np.log2(len(dists)))


# ── loading ──────────────────────────────────────────────────────────────────────
def native_geometry(sweep: Path) -> dict[str, tuple[int, int]]:
    """Read each loop's (pos, normal) from any run's summary.json loop_geometry."""
    geom = {}
    for s in glob.glob(str(sweep / "*" / "*" / "*" / "summary.json")):
        d = json.loads(Path(s).read_text())
        for lp, g in d.get("loop_geometry", {}).items():
            geom.setdefault(lp, (g["pos"], g["normal"]))
    return geom


def native_lengths(sweep: Path) -> dict[str, int]:
    """Each loop's native length."""
    return {lp: n for lp, (_p, n) in native_geometry(sweep).items()}


def native_loop_seqs(sweep: Path, geom: dict[str, tuple[int, int]]) -> dict[str, str]:
    """
    Native (WT) loop sub-sequence per loop, sliced from the prepared design
    construct the sweep actually used (`Local/loop_probe/inputs/*full*.pdb`,
    binder = chain A, renumbered from 1).  The binder is WT TIMP3 in every
    complex, so any cached full-length construct gives the same native residues.
    Returns {} (recovery columns are then omitted) if the construct or biotite
    is unavailable.
    """
    try:
        import numpy as _np
        import biotite.structure.io.pdb as _biopdb
        from biotite.sequence import ProteinSequence
    except Exception:
        return {}
    inputs = sorted((sweep.parent / "inputs").glob("*full*.pdb"))
    if not inputs:
        return {}
    try:
        arr = _biopdb.PDBFile.read(str(inputs[0])).get_structure()[0]
        ca = arr[(arr.chain_id == "A") & (arr.atom_name == "CA")]
        ca = ca[_np.argsort(ca.res_id)]
        seq = "".join(ProteinSequence.convert_letter_3to1(r) if r in
                      ProteinSequence._dict_3to1 else "X" for r in ca.res_name)
    except Exception:
        return {}
    out = {}
    for lp, (pos, normal) in geom.items():
        if pos + normal <= len(seq):
            out[lp] = seq[pos:pos + normal]   # residues pos+1 .. pos+normal (1-indexed pos)
    return out


def load_counts(sweep: Path, target: str, loop: str, length: int) -> pd.DataFrame | None:
    p = sweep / target / loop / f"L{length:02d}" / f"position_counts_{loop}.csv"
    if not p.exists():
        return None
    return pd.read_csv(p, index_col=0).reindex(lpa.AA20).fillna(0).astype(int)


# ── analyses ───────────────────────────────────────────────────────────────────
def hotspot_table(sweep: Path, targets: list[str], nat: dict[str, int],
                  native: dict[str, str] | None = None) -> pd.DataFrame:
    native = native or {}
    rows = []
    for loop, L in nat.items():
        nat_seq = native.get(loop)
        for t in targets:
            c = load_counts(sweep, t, loop, L)
            if c is None:
                continue
            f = freq_matrix(c)
            for j, col in enumerate(f.columns):
                fv = f[col].values
                dom_i = int(np.argmax(fv))
                row = {"target": t, "loop": loop, "length": L, "position": j + 1,
                       "dom_aa": lpa.AA20[dom_i], "dom_aa_freq": round(float(fv[dom_i]), 3),
                       "info_bits": round(info_bits(fv), 3),
                       "metal_ligand_freq": round(float(sum(
                           f.loc[a, col] for a in METAL_LIGANDS)), 3)}
                # native residue at this position (from the WT construct) + recovery
                if nat_seq is not None and j < len(nat_seq):
                    na = nat_seq[j]
                    row["native_aa"] = na
                    row["native_aa_freq"] = round(float(f.loc[na, col]) if na in
                                                   f.index else 0.0, 3)
                    row["native_recovered"] = bool(lpa.AA20[dom_i] == na)
                for scheme in lpa.PROPERTY_SCHEMES:
                    g = lpa.group_counts(f[[col]], scheme)[col]
                    gi = int(np.argmax(g.values))
                    row[f"dom_{scheme}"] = g.index[gi]
                    row[f"dom_{scheme}_freq"] = round(float(g.values[gi]), 3)
                rows.append(row)
    return pd.DataFrame(rows)


def cross_target_table(sweep: Path, targets: list[str], nat: dict[str, int]) -> pd.DataFrame:
    rows = []
    for loop, L in nat.items():
        freqs = {t: freq_matrix(load_counts(sweep, t, loop, L))
                 for t in targets if load_counts(sweep, t, loop, L) is not None}
        if len(freqs) < 2:
            continue
        L_cols = next(iter(freqs.values())).columns
        for j, col in enumerate(L_cols):
            dists = [freqs[t][col].values for t in freqs]
            jsd = js_divergence(dists)
            dom = {t: lpa.AA20[int(np.argmax(freqs[t][col].values))] for t in freqs}
            mean_ic = float(np.mean([info_bits(freqs[t][col].values) for t in freqs]))
            rows.append({"loop": loop, "position": j + 1, "length": L,
                         "cross_target_JSD": round(jsd, 3),
                         "mean_info_bits": round(mean_ic, 3),
                         "n_distinct_dom_aa": len(set(dom.values())),
                         **{f"dom_{t}": dom[t] for t in targets if t in dom}})
    return pd.DataFrame(rows)


def length_trend_table(sweep: Path, targets: list[str], loops: list[str]) -> pd.DataFrame:
    from scipy.stats import spearmanr
    rows = []
    for loop in loops:
        for t in targets:
            for scheme in lpa.PROPERTY_SCHEMES:
                p = sweep / t / loop / f"length_trend_{loop}_{scheme}.csv"
                if not p.exists():
                    continue
                df = pd.read_csv(p, index_col=0)
                lengths = [int(c.lstrip("L")) for c in df.columns]
                for grp in df.index:
                    y = df.loc[grp].values.astype(float)
                    if np.ptp(y) < 1e-9:
                        rho, pv = 0.0, 1.0
                    else:
                        rho, pv = spearmanr(lengths, y)
                    rows.append({"loop": loop, "target": t, "scheme": scheme,
                                 "group": grp, "spearman_rho": round(float(rho), 3),
                                 "p": round(float(pv), 4),
                                 "delta": round(float(y[-1] - y[0]), 3)})
    return pd.DataFrame(rows)


# ── figures ────────────────────────────────────────────────────────────────────
def fig_info_profiles(hot: pd.DataFrame, targets: list[str], loops: list[str], out: Path):
    fig, axes = plt.subplots(1, len(loops), figsize=(3.2 * len(loops), 3.2), dpi=150,
                             squeeze=False)
    for ax, loop in zip(axes[0], loops):
        sub = hot[hot.loop == loop]
        for t in targets:
            d = sub[sub.target == t].sort_values("position")
            if len(d):
                ax.plot(d.position, d.info_bits, marker="o", ms=3, label=t, lw=1.2)
        ax.set_title(loop, fontsize=10)
        ax.set_xlabel("position", fontsize=8)
        ax.set_ylim(0, LOG2_20)
        ax.tick_params(labelsize=7)
    axes[0][0].set_ylabel("information (bits)", fontsize=8)
    axes[0][-1].legend(fontsize=6, loc="upper right")
    fig.suptitle("Per-position constraint (higher = more conserved / hotspot)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def fig_specificity(cross: pd.DataFrame, out: Path):
    fig, ax = plt.subplots(figsize=(6, 5), dpi=150)
    loops = list(cross.loop.unique())
    cmap = plt.get_cmap("tab10")
    for i, loop in enumerate(loops):
        d = cross[cross.loop == loop]
        ax.scatter(d.mean_info_bits, d.cross_target_JSD, s=30, color=cmap(i),
                   label=loop, alpha=0.8, edgecolor="k", linewidth=0.3)
        for _, r in d.iterrows():
            if r.cross_target_JSD > 0.12 or r.mean_info_bits > 2.0:
                ax.annotate(f"{loop}{int(r.position)}", (r.mean_info_bits, r.cross_target_JSD),
                            fontsize=6, xytext=(2, 2), textcoords="offset points")
    ax.set_xlabel("mean information content across targets (bits)")
    ax.set_ylabel("cross-target JS divergence (0-1)")
    ax.set_title("Specificity map: constrained (right) x target-divergent (top)\n"
                 "top-right = target-specific determinant", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def fig_metal(hot: pd.DataFrame, targets: list[str], loops: list[str], out: Path):
    fig, axes = plt.subplots(len(loops), 1, figsize=(7, 1.1 * len(loops) + 1), dpi=150,
                             squeeze=False)
    for ax, loop in zip(axes[:, 0], loops):
        sub = hot[hot.loop == loop]
        L = int(sub.length.iloc[0]) if len(sub) else 0
        mat = np.zeros((len(targets), L))
        for ti, t in enumerate(targets):
            d = sub[sub.target == t].sort_values("position")
            for _, r in d.iterrows():
                mat[ti, int(r.position) - 1] = r.metal_ligand_freq
        im = ax.imshow(mat, aspect="auto", cmap="magma", vmin=0, vmax=1)
        ax.set_yticks(range(len(targets))); ax.set_yticklabels(targets, fontsize=7)
        ax.set_xticks(range(L)); ax.set_xticklabels(range(1, L + 1), fontsize=7)
        ax.set_ylabel(loop, fontsize=9)
    fig.colorbar(im, ax=list(axes[:, 0]), fraction=0.02, pad=0.02,
                 label="His+Cys+Asp+Glu freq")
    fig.suptitle("Metal-coordinating residue propensity per position", fontsize=10)
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)


# ── driver ──────────────────────────────────────────────────────────────────────
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweep_dir")
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args()

    sweep = Path(args.sweep_dir)
    man = json.loads((sweep / "sweep_manifest.json").read_text())
    targets = man["targets"]
    out = Path(args.out_dir) if args.out_dir else sweep / "analysis"
    out.mkdir(parents=True, exist_ok=True)

    geom = native_geometry(sweep)
    nat = {lp: n for lp, (_p, n) in geom.items()}
    natseq = native_loop_seqs(sweep, geom)
    loops = [l for l in ["AB", "C", "EF", "GH", "MTL"] if l in nat]
    print(f"targets={targets}  loops(native len)={ {l: nat[l] for l in loops} }")
    print("native loop sequences:", {l: natseq.get(l, '?') for l in loops})

    hot = hotspot_table(sweep, targets, nat, natseq)
    cross = cross_target_table(sweep, targets, nat)
    trends = length_trend_table(sweep, targets, loops)
    # put the native-recovery columns next to the dominant residue for readability
    if "native_aa" in hot.columns:
        lead = ["target", "loop", "length", "position", "dom_aa", "dom_aa_freq",
                "native_aa", "native_aa_freq", "native_recovered", "info_bits"]
        hot = hot[[c for c in lead if c in hot.columns] +
                  [c for c in hot.columns if c not in lead]]
    hot.to_csv(out / "hotspots_native.csv", index=False)
    cross.to_csv(out / "cross_target_divergence.csv", index=False)
    trends.to_csv(out / "length_trends.csv", index=False)

    fig_info_profiles(hot, targets, loops, out / "fig_info_content.png")
    fig_specificity(cross, out / "fig_specificity_map.png")
    fig_metal(hot, targets, loops, out / "fig_metal_propensity.png")

    # ── printed summary ─────────────────────────────────────────────────────────
    def bar():
        print("=" * 78)
    bar(); print("HOTSPOTS — strongest per-position constraints (native length)"); bar()
    for loop in loops:
        sub = hot[hot.loop == loop].sort_values("info_bits", ascending=False).head(4)
        print(f"\n[{loop}] top constrained positions (bits | dom AA | charge | hydrophobicity):")
        for _, r in sub.iterrows():
            print(f"  {r.target:6s} pos{int(r.position):2d}: {r.info_bits:4.2f}b  "
                  f"{r.dom_aa}={r.dom_aa_freq:.2f}  "
                  f"{r.dom_charge}({r.dom_charge_freq:.2f})  "
                  f"{r.dom_hydrophobicity}({r.dom_hydrophobicity_freq:.2f})")

    if "native_recovered" in hot.columns:
        bar(); print("NATIVE RECOVERY — fraction of positions where dom AA == WT residue"); bar()
        for loop in loops:
            sub = hot[hot.loop == loop]
            ns = natseq.get(loop, "")
            print(f"\n[{loop}] native={ns}")
            for t in targets:
                d = sub[sub.target == t].sort_values("position")
                if not len(d):
                    continue
                rec = "".join(r.dom_aa if r.native_recovered else "." for _, r in d.iterrows())
                frac = float(d.native_recovered.mean())
                print(f"  {t:6s} recovered {int(d.native_recovered.sum())}/{len(d)} "
                      f"({frac:.0%})  dom={''.join(d.dom_aa)}  match={rec}")

    bar(); print("CROSS-TARGET — most target-DIVERGENT positions (specificity determinants)"); bar()
    for _, r in cross.sort_values("cross_target_JSD", ascending=False).head(12).iterrows():
        doms = " ".join(f"{t[:5]}={r.get('dom_'+t,'?')}" for t in targets)
        print(f"  {r.loop:3s} pos{int(r.position):2d}  JSD={r.cross_target_JSD:.2f}  "
              f"IC={r.mean_info_bits:.2f}b  [{doms}]")
    bar(); print("CROSS-TARGET — most SHARED constrained positions (conserved across targets)"); bar()
    shared = cross[cross.mean_info_bits > 0.5].sort_values("cross_target_JSD").head(10)
    for _, r in shared.iterrows():
        doms = " ".join(f"{t[:5]}={r.get('dom_'+t,'?')}" for t in targets)
        print(f"  {r.loop:3s} pos{int(r.position):2d}  JSD={r.cross_target_JSD:.2f}  "
              f"IC={r.mean_info_bits:.2f}b  [{doms}]")

    bar(); print("LENGTH TRENDS — strongest composition shifts vs loop length (|rho|>=0.7)"); bar()
    strong = trends[(trends.spearman_rho.abs() >= 0.7)].copy()
    strong = strong.sort_values("spearman_rho", key=lambda s: s.abs(), ascending=False)
    seen = set()
    for _, r in strong.iterrows():
        k = (r.loop, r.scheme, r.group)
        if k in seen:
            continue
        seen.add(k)
        # only report a trend consistent across most targets
        same = trends[(trends.loop == r.loop) & (trends.scheme == r.scheme) &
                      (trends.group == r.group)]
        if (np.sign(same.spearman_rho) == np.sign(r.spearman_rho)).mean() >= 0.75 \
           and (same.spearman_rho.abs() >= 0.5).mean() >= 0.75:
            print(f"  {r.loop:3s} {r.scheme:14s} {r.group:22s} rho~{same.spearman_rho.mean():+.2f} "
                  f"(mean over targets), delta~{same.delta.mean():+.2f}")

    bar(); print("METAL-COORDINATION — positions with high His/Cys/Asp/Glu propensity"); bar()
    mh = hot.sort_values("metal_ligand_freq", ascending=False).head(12)
    for _, r in mh.iterrows():
        print(f"  {r.target:6s} {r.loop:3s} pos{int(r.position):2d}: "
              f"H/C/D/E={r.metal_ligand_freq:.2f}  (dom {r.dom_aa}={r.dom_aa_freq:.2f})")
    print(f"\nWrote CSVs + figures to {out}")


if __name__ == "__main__":
    main()
