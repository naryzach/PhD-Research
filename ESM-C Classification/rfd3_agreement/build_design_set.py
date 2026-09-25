"""Assemble the RFd3-pipeline design set that the ESM-C models will be asked about.

Writes to Local/esmc_rfd3_agreement/:
  design_meta.csv   one row per (design, target): structural metrics from the pipeline
  seqs.csv          one row per unique sequence to score (loop windows resolved), with the
                    ``priority`` tier so the cheap/important sets can be scored first
  controls.csv      held-out lab-assay sequences (real binder / non-binder labels) that anchor
                    what a "binder-like" ESM-C probability looks like

Sets (a sequence can belong to several; ``sets`` lists them):
  af3        AF3-Server-folded designs (strongest structural evidence)
  shortlist  Construct_Shortlist_2026-08 (manufacturing candidates)
  order      FINAL_ORDER_2026-08-27 + specificity candidates actually ordered
  hof        hall-of-fame (top 300 of the campaign)
  pool       stratified sample of the full pool across sv_pdockq quantiles, per target
  wt         WT TIMP3 and the lab-tested constructs

    python build_design_set.py --pool-per-bin 40
"""
from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "Local"
POOL_DIR = LOCAL / "iterative_refinement_20260831"
OUT = LOCAL / "esmc_rfd3_agreement"

# Same flank tripeptides the design pipeline uses to locate each redesigned loop.
LOOP_FLANKS = {"AB": ("LVK", "LVY"), "C": ("HTE", "GLK"), "EF": ("MYT", "FVE")}


def locate_loops(seq: str) -> dict:
    """{loop: (char_start, length)} located sequentially by flanks (0-indexed start)."""
    out, cursor = {}, 0
    for name, (left, right) in LOOP_FLANKS.items():
        m = re.compile(f"{left}([A-Z]*?){right}").search(seq[cursor:])
        if m:
            start = cursor + m.start() + len(left)
            out[name] = (start, len(m.group(1)))
            cursor = cursor + m.end() - len(right)
        else:
            out[name] = (-1, 0)
    return out


def load_pool() -> pd.DataFrame:
    rows = []
    for f in glob.glob(str(POOL_DIR / "it_*" / "round_summary.csv")):
        d = pd.read_csv(f, on_bad_lines="skip")
        d["iteration"] = int(re.search(r"it_(\d+)", f).group(1))
        rows.append(d)
    return pd.concat(rows, ignore_index=True)


def load_af3() -> pd.DataFrame:
    """Union of every table that carries an AF3 ipTM, keyed by design_id."""
    cal = LOCAL / "Pipeline_Calibration_2026-08" / "data"
    parts = []
    a = pd.read_csv(POOL_DIR / "af3_all.csv")
    parts.append(a.rename(columns={"binder_seq": "full_seq", "target_name": "target"})
                 [["design_id", "target", "full_seq", "af3_iptm", "af3_ptm"]])
    for f, tcol in [("FINAL_ORDER_2026-08-27.csv", "target"), ("af3_order_topup.csv", "target")]:
        d = pd.read_csv(cal / f)
        keep = ["design_id", tcol, "full_seq", "af3_iptm"] + (["af3_ptm"] if "af3_ptm" in d else [])
        parts.append(d[keep].rename(columns={tcol: "target"}))
    s = pd.read_csv(LOCAL / "Construct_Shortlist_2026-08" / "construct_shortlist.csv")
    parts.append(s.rename(columns={"target_name": "target"})
                 [["design_id", "target", "full_seq", "af3_iptm", "af3_ptm"]])
    u = pd.concat(parts, ignore_index=True).dropna(subset=["af3_iptm", "full_seq"])
    return u.sort_values("af3_iptm", ascending=False).drop_duplicates("design_id")


def stratified_pool(pool: pd.DataFrame, per_bin: int, bins: int, seed: int) -> pd.DataFrame:
    p = pool.dropna(subset=["sv_pdockq"]).drop_duplicates(["full_seq", "target_name"]).copy()
    out = []
    for tgt, g in p.groupby("target_name"):
        g = g.copy()
        g["_bin"] = pd.qcut(g["sv_pdockq"], bins, labels=False, duplicates="drop")
        for _, gb in g.groupby("_bin"):
            out.append(gb.sample(min(per_bin, len(gb)), random_state=seed))
    return pd.concat(out, ignore_index=True).drop(columns="_bin")


def controls(per_target_label: int, seed: int) -> pd.DataFrame:
    """Held-out lab-assay rows from each variant's test split (real binder labels)."""
    rows = []
    for variant_dir in sorted((LOCAL / "esmc_multirun").glob("*/data/test.parquet")):
        variant = variant_dir.parents[1].name
        t = pd.read_parquet(variant_dir)
        for tgt in [c[len("mask_"):] for c in t.columns if c.startswith("mask_")]:
            sub = t[t[f"mask_{tgt}"] == 1]
            for lab, g in sub.groupby(f"label_{tgt}"):
                g = g.sample(min(per_target_label, len(g)), random_state=seed)
                rows.append(pd.DataFrame({
                    "variant": variant, "target": tgt, "label": int(lab),
                    "full_seq": g["sequence"].values, "loop": g["loop"].values,
                    "loop_start": g["loop_start"].values,
                }))
    return pd.concat(rows, ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool-per-bin", type=int, default=40)
    ap.add_argument("--pool-bins", type=int, default=10)
    ap.add_argument("--controls-per-target-label", type=int, default=75)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    pool = load_pool()
    metric_cols = ["design_id", "target_name", "full_seq", "iteration", "bb_origin", "sv_pdockq",
                   "esm_iptm", "esm_ptm", "esm_plddt", "esm_iface_contact_density",
                   "esm_iface_n_iface_res", "composite_score", "sv_n_iface_res_binder",
                   "sv_sc_shape_complementarity", "loop_AB_seq", "loop_C_seq", "loop_EF_seq"]
    pool_m = pool[[c for c in metric_cols if c in pool.columns]].copy()
    af3 = load_af3()
    sl = pd.read_csv(LOCAL / "Construct_Shortlist_2026-08" / "construct_shortlist.csv")
    fo = pd.read_csv(LOCAL / "Pipeline_Calibration_2026-08" / "data" / "FINAL_ORDER_2026-08-27.csv")
    sp = pd.read_csv(LOCAL / "Pipeline_Calibration_2026-08" / "data" / "spec_candidates_2026-09-11.csv")
    hof = pd.read_csv(POOL_DIR / "hof_summary.csv")
    strat = stratified_pool(pool, args.pool_per_bin, args.pool_bins, args.seed)

    # (design_id, target) -> sets
    tagged = []
    def tag(df, idcol, tcol, seqcol, name):
        t = pd.DataFrame({"design_id": df[idcol].values, "target": df[tcol].values,
                          "full_seq": df[seqcol].values, "set": name})
        tagged.append(t)
    tag(af3, "design_id", "target", "full_seq", "af3")
    tag(sl, "design_id", "target_name", "full_seq", "shortlist")
    tag(fo, "design_id", "target", "full_seq", "order")
    tag(sp, "design_id", "target_name", "full_seq", "order")
    tag(hof, "design_id", "target_name", "full_seq", "hof")
    tag(strat, "design_id", "target_name", "full_seq", "pool")
    tg = pd.concat(tagged, ignore_index=True)
    sets = tg.groupby(["design_id", "target", "full_seq"])["set"].agg(
        lambda s: ",".join(sorted(set(s)))).reset_index().rename(columns={"set": "sets"})

    meta = sets.merge(pool_m.rename(columns={"target_name": "target"}).drop(columns=["full_seq"])
                      .drop_duplicates(["design_id", "target"]),
                      on=["design_id", "target"], how="left")
    meta = meta.merge(af3[["design_id", "af3_iptm", "af3_ptm"]], on="design_id", how="left")
    meta["n_cys"] = meta["full_seq"].str.count("C")
    for name in ("AB", "C", "EF"):
        meta[f"len_{name}"] = [locate_loops(s)[name][1] for s in meta["full_seq"]]

    # Sequences to score. Tier 1 = structural-evidence sets; tier 2 = pool sample.
    seqs = meta.groupby("full_seq")["sets"].agg(lambda s: ",".join(sorted(set(",".join(s).split(","))))).reset_index()
    seqs["priority"] = np.where(seqs["sets"].str.fullmatch("pool"), 2, 1)
    seqs.insert(0, "seq_id", [f"S{i:05d}" for i in range(len(seqs))])
    for name in ("AB", "C", "EF"):
        loc = [locate_loops(s)[name] for s in seqs["full_seq"]]
        seqs[f"start_{name}"] = [l[0] for l in loc]
        seqs[f"len_{name}"] = [l[1] for l in loc]

    # Sanity: flank-located loops agree with the pipeline's own recorded loop strings.
    chk = meta.drop_duplicates("design_id").dropna(subset=["loop_AB_seq", "loop_C_seq"])
    mism = sum(
        (locate_loops(r.full_seq)["AB"][1] != len(r.loop_AB_seq)) or
        (locate_loops(r.full_seq)["C"][1] != len(r.loop_C_seq)) for r in chk.itertuples())
    print(f"loop-length check vs pipeline columns: {mism}/{len(chk)} mismatches")

    # Lab-tested constructs + WT (from the sequences already in the analysis folder).
    pva = pd.read_csv(LOCAL / "Prediction_vs_Result_Analysis" / "predicted_vs_actual.csv")
    lab = pd.DataFrame({"full_seq": pva["Sequence"], "sets": "lab_construct"}).drop_duplicates("full_seq")
    lab = lab[~lab.full_seq.isin(seqs.full_seq)]
    if len(lab):
        base = len(seqs)
        lab.insert(0, "seq_id", [f"S{base + i:05d}" for i in range(len(lab))])
        lab["priority"] = 1
        for name in ("AB", "C", "EF"):
            loc = [locate_loops(s)[name] for s in lab["full_seq"]]
            lab[f"start_{name}"] = [l[0] for l in loc]
            lab[f"len_{name}"] = [l[1] for l in loc]
        seqs = pd.concat([seqs, lab], ignore_index=True)

    ctrl = controls(args.controls_per_target_label, args.seed)
    meta.to_csv(OUT / "design_meta.csv", index=False)
    seqs.to_csv(OUT / "seqs.csv", index=False)
    ctrl.to_csv(OUT / "controls.csv", index=False)

    print(f"design_meta rows (design x target): {len(meta)}")
    print(f"unique sequences to score: {len(seqs)}  (priority 1: {(seqs.priority == 1).sum()}, "
          f"priority 2: {(seqs.priority == 2).sum()})")
    print(f"controls: {len(ctrl)} lab-assay rows")
    print(meta.groupby("target")["sets"].apply(lambda s: s.str.split(",").explode().value_counts().to_dict()))


if __name__ == "__main__":
    main()
