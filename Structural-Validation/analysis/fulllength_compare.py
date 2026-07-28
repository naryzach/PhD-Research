"""
Full-length (188-aa) vs N-domain (121-aa) TIMP3 co-fold comparison.

The design/validation pipeline modelled the 121-aa N-domain in isolation, but the
ordered/tested protein is full length (N-domain + native C-domain), and the C-domain
is needed for ADAM10 binding. This scores the fresh 188-aa WT:target co-folds on the
same interface battery as the 121-aa co-folds and reports what changes.

Only WT:target co-folds exist at 188 (6 targets). A full construct-level FCS
re-correlation additionally needs 188-aa CONSTRUCT folds, which don't exist yet.

Outputs (Local/TIMP3_FullLength_Validation_2026-07/):
  wt_121_vs_188.csv

Run:  python Structural-Validation/analysis/fulllength_compare.py
"""
from __future__ import annotations

import csv
import glob
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
import config as C  # noqa: E402
import metrics as M  # noqa: E402
import structure_io as sio  # noqa: E402
from complex_metrics import _seq_identity, assign_chains  # noqa: E402

TIMP121 = C.TIMP3_WT_MATURE[:121]
TARGETS = ["ADAM10", "ADAM17", "MMP2", "MMP9", "MMP3", "MMP10"]
OUT = C.REPO_ROOT / "Local" / "TIMP3_FullLength_Validation_2026-07"

FOLD188 = C.REPO_ROOT / "Data/TIMP_Complexes/AF3_Full_Folds/timp3full_{t}/fold_timp3full_{t}_model_0.cif"
FOLD121 = C.OUT_AF3 / "results" / "timp3_wt_{t}" / "fold_timp3_wt_{t}_model_0.cif"


def _find(pattern_path: Path):
    hits = [p for p in glob.glob(str(pattern_path)) if ":" not in Path(p).name]
    return hits[0] if hits else None


def score(cif: str, target: str) -> dict:
    chains = sio.get_chains(cif)
    binder = max(chains.values(), key=lambda c: _seq_identity(c.seq[:121], TIMP121))
    tgt = max((c for c in chains.values() if c.cid != binder.cid), key=lambda c: len(c.seq))
    motif = sio.zinc_motif_resids(tgt)
    im = M.interface_summary(cif, binder.cid, tgt.cid, motif_resids_b=motif or None)
    cm = M.confidence_metrics(binder, tgt, model_path=cif, cid=binder.cid, tid=tgt.cid)
    row = {
        "binder_len": len(binder.seq), "target_len": len(tgt.seq),
        "bsa": im.bsa, "n_hbonds": im.n_hbonds, "n_salt_bridges": im.n_salt_bridges,
        "sc": im.sc_shape_complementarity, "contact_density": im.contact_density,
        "catalytic_occlusion": im.catalytic_occlusion,
        "interface_plddt": cm.get("interface_plddt", ""), "pdockq": cm.get("pdockq", ""),
        "pdockq2": cm.get("pdockq2", ""), "lis": cm.get("lis", ""),
        "interface_pae": cm.get("interface_pae", ""),
        "n_iface_res_binder": im.n_iface_res_A,
    }
    # fraction of the binder interface contributed by the C-domain (resid > 121)
    ac = np.array([r["CA"].coord for r in binder.residues])
    bc = np.array([a.coord for r in tgt.residues for a in r])
    tree = cKDTree(bc)
    ndom = sum(1 for i, r in enumerate(binder.residues)
               if r.id[1] <= 121 and tree.query_ball_point(ac[i], 8.0))
    cdom = sum(1 for i, r in enumerate(binder.residues)
               if r.id[1] > 121 and tree.query_ball_point(ac[i], 8.0))
    row["iface_res_Ndom"] = ndom
    row["iface_res_Cdom"] = cdom
    row["cdom_iface_frac"] = round(cdom / (ndom + cdom), 3) if (ndom + cdom) else 0.0

    # DockQ vs native, ADAM17 only
    if target == "ADAM17":
        ref, _ = None, None
        rp = C.DATA_DIR / "TIMP3_vs_ADAM17_X_ray.pdb"
        if rp.exists():
            try:
                rc = sio.get_chains(rp)
                rcc, rtc = assign_chains(rc, binder.seq, tgt.seq)
                dq = M.dockq(cif, rp, tgt.cid, binder.cid, rtc.cid, rcc.cid)
                row["dockq_native"] = dq.dockq
                row["lrms_native"] = dq.lrms
                row["fnat_native"] = dq.fnat
            except Exception as e:
                row["dockq_native"] = f"err:{e}"
    return row


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for t in TARGETS:
        f188 = _find(Path(str(FOLD188).format(t=t.lower())))
        f121 = _find(Path(str(FOLD121).format(t=t.lower())))
        for length, f in (("121_Ndomain", f121), ("188_fulllen", f188)):
            if not f:
                rows.append({"target": t, "length": length, "status": "missing"})
                continue
            try:
                r = {"target": t, "length": length, "status": "ok"}
                r.update(score(f, t))
                rows.append(r)
            except Exception as e:
                rows.append({"target": t, "length": length, "status": f"err:{e}"})

    fields = ["target", "length", "status", "binder_len", "target_len", "bsa",
              "n_hbonds", "n_salt_bridges", "sc", "contact_density",
              "catalytic_occlusion", "interface_plddt", "pdockq", "pdockq2", "lis",
              "interface_pae", "n_iface_res_binder", "iface_res_Ndom",
              "iface_res_Cdom", "cdom_iface_frac", "dockq_native", "lrms_native",
              "fnat_native"]
    with open(OUT / "wt_121_vs_188.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    # ── comparison table ──
    ok = {(r["target"], r["length"]): r for r in rows if r.get("status") == "ok"}
    print(f"{'target':8s}{'len':>6}{'bsa':>7}{'hbnd':>5}{'Sc':>6}{'catOcc':>8}"
          f"{'pdockq2':>8}{'iPAE':>6}{'Cdom%iface':>11}{'DockQ':>7}")
    print("-" * 72)
    for t in TARGETS:
        for L in ("121_Ndomain", "188_fulllen"):
            r = ok.get((t, L))
            if not r:
                continue
            dq = r.get("dockq_native", "")
            dq = f"{dq:.2f}" if isinstance(dq, float) else "-"
            print(f"{t:8s}{L.split('_')[0]:>6}{r['bsa']:>7.0f}{r['n_hbonds']:>5}"
                  f"{r['sc']:>6.2f}{r['catalytic_occlusion']:>8}"
                  f"{str(r['pdockq2']):>8}{str(r['interface_pae']):>6}"
                  f"{r['cdom_iface_frac']:>11}{dq:>7}")
    print(f"\n-> {(OUT / 'wt_121_vs_188.csv').relative_to(C.REPO_ROOT)}")


if __name__ == "__main__":
    main()
