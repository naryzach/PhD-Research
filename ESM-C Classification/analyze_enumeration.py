"""Analytics over the exhaustive C-loop enumeration outputs.

Reads the per-target top-K files and the selective table written by
enumerate_cloop.py and writes analyses + figures to <enum_dir>/analysis/:
  * position x amino-acid heatmaps (frequency AND log2-enrichment) for each
    target's top binders and the target-selective loops
  * predicted-probability and selectivity-margin distributions
  * top-K loop overlap (Jaccard) between targets -> pan-target "stickiness"
  * novelty vs the measured C-loops (how far the picks extrapolate)
  * motif-deduplicated, novelty-annotated shortlist (a bench list)
  * summary.json + summary.txt

Usage:
    python analyze_enumeration.py --config config.yaml
    python analyze_enumeration.py --config config.yaml --dir 0_64000000
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from esmc_utils import load_config, resolve_path

AAS = list("ACDEFGHIKLMNPQRSTVWY")
AA_IDX = {a: i for i, a in enumerate(AAS)}


def find_enum_dir(cfg, override):
    base = resolve_path(cfg, cfg["output_dir"]) / "enumeration"
    if override:
        d = Path(override)
        return d if d.is_absolute() else base / override
    subs = [p for p in base.iterdir() if p.is_dir()]
    if not subs:
        raise SystemExit(f"No enumeration subfolders under {base}")
    return max(subs, key=lambda p: p.stat().st_mtime)


def loops_to_counts(loops, L):
    counts = np.zeros((L, 20))
    for s in loops:
        for p, ch in enumerate(s):
            j = AA_IDX.get(ch)
            if j is not None:
                counts[p, j] += 1
    return counts


def heatmap(counts, title, path, enrichment=False):
    L = counts.shape[0]
    freq = counts / counts.sum(axis=1, keepdims=True).clip(min=1)
    consensus = "".join(AAS[int(np.argmax(freq[p]))] for p in range(L))
    if enrichment:
        mat = np.log2((freq.T + 1e-9) / (1.0 / 20))
        cmap, center, label = "RdBu_r", 0.0, "log2 enrichment vs uniform (1/20)"
    else:
        mat = freq.T
        cmap, center, label = "viridis", None, "frequency"
    plt.figure(figsize=(1.1 * L + 4, 9))
    sns.heatmap(mat, yticklabels=AAS, xticklabels=[f"P{p+1}" for p in range(L)],
                cmap=cmap, center=center, linewidths=0.3, cbar_kws={"label": label})
    plt.title(f"{title}\nconsensus: {consensus}")
    plt.xlabel("loop position"); plt.ylabel("amino acid")
    plt.tight_layout(); plt.savefig(path, dpi=150); plt.close()
    return consensus


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def dedup_by_motif(loops_scores, k, ham):
    reps = []
    for loop, score in loops_scores:
        if all(hamming(loop, r[0]) >= ham for r in reps):
            reps.append((loop, score))
            if len(reps) >= k:
                break
    return reps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--dir", default=None, help="enumeration subfolder (default: newest)")
    ap.add_argument("--loop-subtype", default=None,
                    help="Ligand_Subtype to use as the novelty reference (default: read from "
                         "run_meta.json written by enumerate_cloop.py, else 'C-loop' for "
                         "backward compatibility with pre-multirun enumeration runs)")
    ap.add_argument("--template-csv", default=None,
                    help="CSV to compute the novelty reference from (default: cfg data.csv_path)")
    ap.add_argument("--top", type=int, default=3000, help="rows/target to profile for heatmaps")
    ap.add_argument("--shortlist", type=int, default=60, help="deduped reps per list")
    ap.add_argument("--ham", type=int, default=2, help="min Hamming between shortlist reps")
    args = ap.parse_args()

    cfg = load_config(args.config)
    enum_dir = find_enum_dir(cfg, args.dir)
    out = enum_dir / "analysis"; out.mkdir(parents=True, exist_ok=True)
    print(f"Analyzing {enum_dir}\n -> {out}")

    loop_subtype = args.loop_subtype
    run_meta_path = enum_dir / "run_meta.json"
    if loop_subtype is None:
        if run_meta_path.exists():
            loop_subtype = json.loads(run_meta_path.read_text()).get("loop_subtype")
            print(f"loop_subtype (from run_meta.json): {loop_subtype!r}")
        else:
            loop_subtype = "C-loop"  # pre-multirun enumeration runs were always C-loop-only
            print(f"loop_subtype: no run_meta.json found, defaulting to {loop_subtype!r}")

    targets = sorted(p.stem.replace("top_", "") for p in enum_dir.glob("top_*.csv")
                     if p.stem != "top_selective")
    print("targets:", targets)
    summary = {"enum_dir": str(enum_dir), "targets": targets}

    # measured loops (of the enumerated subtype) = novelty reference
    d = cfg["data"]
    template_csv = args.template_csv or resolve_path(cfg, d["csv_path"])
    train = pd.read_csv(template_csv)
    sub = "Ligand_Subtype"
    col = train[train[sub] == loop_subtype] if (loop_subtype and sub in train.columns) else train
    train_loops = set(col[d["loop_col"]].dropna().astype(str))
    Lref = len(next(iter(train_loops))) if train_loops else 6
    train_arr = (np.array([[AA_IDX.get(c, -1) for c in s] for s in train_loops])
                 if train_loops else np.zeros((0, Lref), dtype=int))
    del_keys = set()
    for s in train_loops:
        del_keys.add(s)
        for i in range(len(s)):
            del_keys.add(s[:i] + s[i + 1:])

    def within1(s):
        return s in train_loops or any((s[:i] + s[i + 1:]) in del_keys for i in range(len(s)))

    def min_ham_to_train(s):
        if len(train_arr) == 0:
            return -1
        sc = np.array([AA_IDX.get(c, -1) for c in s])
        return int((train_arr != sc).sum(axis=1).min())

    summary["n_train_cloops"] = len(train_loops)

    consensus, prob_summary, top_sets = {}, {}, {}
    for t in targets:
        df = pd.read_csv(enum_dir / f"top_{t}.csv")
        pcol = f"prob_{t}"
        top_sets[t] = set(df["loop"].astype(str))
        loops = df["loop"].astype(str).head(args.top).tolist()
        L = len(loops[0])
        counts = loops_to_counts(loops, L)
        consensus[t] = heatmap(counts, f"{t}: top {len(loops)} binders (position x AA)",
                               out / f"aa_heatmap_top_{t}.png")
        heatmap(counts, f"{t}: top {len(loops)} binders",
                out / f"aa_heatmap_top_{t}_enrichment.png", enrichment=True)
        prob_summary[t] = {"n": int(len(df)), "max": float(df[pcol].max()),
                           "min_kept": float(df[pcol].min()),
                           "frac_ge_0.99": float((df[pcol] >= 0.99).mean())}
    summary["top_consensus"] = consensus
    summary["prob_summary"] = prob_summary

    ov = pd.DataFrame(index=targets, columns=targets, dtype=float)
    for a in targets:
        for b in targets:
            ia, ib = top_sets[a], top_sets[b]
            ov.loc[a, b] = len(ia & ib) / len(ia | ib) if (ia or ib) else np.nan
    ov = ov.astype(float)
    ov.to_csv(out / "topk_overlap_jaccard.csv")
    plt.figure(figsize=(6, 5))
    sns.heatmap(ov, annot=True, vmin=0, vmax=1, cmap="Reds")
    plt.title("Top-K loop overlap (Jaccard)\nhigh off-diagonal = pan-target 'stickiness'")
    plt.tight_layout(); plt.savefig(out / "topk_overlap_jaccard.png", dpi=150); plt.close()
    summary["topk_overlap_jaccard"] = ov.round(4).to_dict()

    plt.figure(figsize=(8, 5))
    for t in targets:
        df = pd.read_csv(enum_dir / f"top_{t}.csv")
        plt.hist(df[f"prob_{t}"], bins=60, alpha=0.5, label=t)
    plt.xlabel("predicted P(bind) among kept top-K"); plt.ylabel("count"); plt.legend()
    plt.title("Predicted-probability distribution (top-K per target)")
    plt.tight_layout(); plt.savefig(out / "prob_distributions.png", dpi=150); plt.close()

    sel_path = enum_dir / "top_selective.csv"
    if sel_path.exists():
        sel = pd.read_csv(sel_path)
        summary["selective_counts"] = {k: int(v) for k, v in sel["top_target"].value_counts().items()}
        summary["selective_margin"] = {k: float(v) for k, v in sel["margin"].describe().items()}
        plt.figure(figsize=(8, 5))
        for t, g in sel.groupby("top_target"):
            plt.hist(g["margin"], bins=60, alpha=0.5, label=f"{t} (n={len(g)})")
        plt.xlabel("selectivity margin (top_prob - runnerup_prob)"); plt.ylabel("count")
        plt.legend(); plt.title("Selectivity-margin distribution by winning target")
        plt.tight_layout(); plt.savefig(out / "margin_distribution.png", dpi=150); plt.close()
        for t in sel["top_target"].unique():
            loops = sel[sel["top_target"] == t]["loop"].astype(str).head(args.top).tolist()
            if len(loops) >= 10:
                heatmap(loops_to_counts(loops, len(loops[0])),
                        f"{t}-selective: top {len(loops)} by margin (position x AA)",
                        out / f"aa_heatmap_selective_{t}.png")

    nov = {}
    for t in targets:
        loops = pd.read_csv(enum_dir / f"top_{t}.csv")["loop"].astype(str)
        nov[t] = {"frac_in_train": float(loops.isin(train_loops).mean()),
                  "frac_within1_edit_first5k": float(np.mean([within1(s) for s in loops.head(5000)]))}
    summary["novelty_topK"] = nov

    rows = []
    for t in targets:
        df = pd.read_csv(enum_dir / f"top_{t}.csv")
        for loop, score in dedup_by_motif(list(zip(df["loop"].astype(str), df[f"prob_{t}"])),
                                          args.shortlist, args.ham):
            rows.append({"target": t, "mode": "top_binder", "loop": loop, "score": float(score),
                         "min_hamming_to_train": min_ham_to_train(loop), "in_train": loop in train_loops})
    if sel_path.exists():
        for t in sel["top_target"].unique():
            g = sel[sel["top_target"] == t]
            for loop, score in dedup_by_motif(list(zip(g["loop"].astype(str), g["margin"])),
                                              args.shortlist, args.ham):
                rows.append({"target": t, "mode": "selective", "loop": loop, "score": float(score),
                             "min_hamming_to_train": min_ham_to_train(loop), "in_train": loop in train_loops})
    sl = pd.DataFrame(rows)
    sl.to_csv(out / "shortlist_dedup.csv", index=False)
    if len(sl):
        summary["shortlist"] = {"rows": int(len(sl)),
            "min_hamming_to_train_dist": {int(k): int(v) for k, v in
                sl["min_hamming_to_train"].value_counts().sort_index().items()},
            "frac_in_train": float(sl["in_train"].mean())}

    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    lines = ["ENUMERATION ANALYSIS SUMMARY", "=" * 44,
             f"enumeration dir : {enum_dir}", f"targets         : {targets}",
             f"train C-loops   : {summary['n_train_cloops']} (novelty reference)",
             "", "Top-binder consensus motif per target:"]
    for t in targets:
        lines.append(f"  {t:<8} {consensus[t]}   (max P={prob_summary[t]['max']:.4f}, "
                     f"frac>=0.99={prob_summary[t]['frac_ge_0.99']:.3f})")
    lines += ["", "Top-K overlap (Jaccard) -- high off-diagonal = sticky/pan-target:",
              ov.round(3).to_string()]
    if sel_path.exists():
        lines += ["", f"Selective loops by winning target: {summary['selective_counts']}",
                  f"margin: max={summary['selective_margin']['max']:.3f} mean={summary['selective_margin']['mean']:.3f}"]
    lines += ["", "Novelty of top-K vs measured C-loops (extrapolation check):"]
    for t in targets:
        lines.append(f"  {t:<8} in-train={nov[t]['frac_in_train']:.4f}  "
                     f"within-1-edit(first5k)={nov[t]['frac_within1_edit_first5k']:.4f}")
    if len(sl):
        lines += ["", f"Deduped shortlist: {len(sl)} reps -> shortlist_dedup.csv "
                      f"(min-Hamming-to-train dist: {summary['shortlist']['min_hamming_to_train_dist']})"]
    (out / "summary.txt").write_text("\n".join(lines))
    print("\n" + "\n".join(lines))
    print(f"\nFigures + tables written to {out}")


if __name__ == "__main__":
    main()