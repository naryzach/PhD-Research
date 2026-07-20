"""
TEST 1 — mechanistic chelation test (crystal-independent).

TIMP inhibition works by the N-terminal Cys1-Thr2 chelating the catalytic zinc:
the reactive edge sits among the three Zn-coordinating histidines and the
catalytic glutamate of the HExxHxxGxxH motif. This test asks whether each docked /
co-folded pose reproduces that geometry — using only the motif found in the target
chain itself, with NO reference to any crystal structure.

Why it matters: our only native complex benchmark is TIMP3:ADAM17, which AF3 and
ESMFold2 have almost certainly seen in training, so DockQ 0.76 (co-folds) vs 0.05
(HADDOCK) may partly reflect memorization. This test cannot be won by memorization
— it only asks whether the chemistry is right.

For each pose:
  key atoms = heavy atoms of motif[0], motif[4], motif[10] (the three Zn-His)
              and motif[1] (the catalytic Glu)
  d_res1..d_res6 = min heavy-atom distance from each of the construct's first six
              residues to that key-atom set
  closest_res  = which of residues 1-6 is nearest
  edge_leads   = closest_res <= 2  (i.e. the C1/T2 reactive edge leads, as in the
                 real chelation mechanism, rather than some other residue)

A REFERENCE row is computed identically on the native TIMP3:ADAM17 crystal, so the
model rows can be read against the real thing on the same scale.

Output: analysis/chelation_test.csv

Run:  python Structural-Validation/analysis/chelation_test.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
import config as C  # noqa: E402
import model_registry as MR  # noqa: E402
import structure_io as sio  # noqa: E402
from complex_metrics import _resolve, _source_specs, assign_chains  # noqa: E402

N_EDGE = 6                     # construct residues profiled (1..6)
ZN_HIS_IDX = (0, 4, 10)        # H..H..H within HExxHxxGxxH
CAT_GLU_IDX = 1                # the catalytic E
FIELDS = (["construct_id", "target_id", "source", "family", "construct_src",
           "target_src", "status"]
          + [f"d_res{i}" for i in range(1, N_EDGE + 1)]
          + ["closest_res", "edge_leads", "n_motif", "model_path"])


def _res_heavy(res) -> np.ndarray:
    return np.asarray([a.coord for a in res if a.element != "H"], dtype=float)


def chelation_geometry(path, construct_seq: str, target_seq: str) -> dict:
    """Min distances from the construct's reactive edge to the catalytic key atoms."""
    chains = sio.get_chains(path)
    if len(chains) < 2:
        return {"status": "single_chain_only"}
    cchain, tchain = assign_chains(chains, construct_seq, target_seq)
    return chelation_from_chains(cchain, tchain)


def chelation_from_chains(cchain, tchain) -> dict:
    """Same measurement given already-assigned construct/target chains."""
    motif = sio.zinc_motif_resids(tchain)
    if len(motif) < 11:
        return {"status": f"no_zinc_motif(n={len(motif)})"}
    key_resids = {motif[i] for i in ZN_HIS_IDX} | {motif[CAT_GLU_IDX]}

    key_coords = [c for r in tchain.residues if r.id[1] in key_resids
                  for c in ([] if len(_res_heavy(r)) == 0 else _res_heavy(r))]
    if not key_coords:
        return {"status": "no_key_atoms"}
    tree = cKDTree(np.asarray(key_coords, dtype=float))

    row = {"n_motif": len(motif)}
    edge_ids = cchain.resids[:N_EDGE]
    by_id = {r.id[1]: r for r in cchain.residues}
    dists = []
    for k, rid in enumerate(edge_ids, start=1):
        res = by_id.get(rid)
        coords = _res_heavy(res) if res is not None else np.empty((0, 3))
        d = float(tree.query(coords, k=1)[0].min()) if len(coords) else float("nan")
        row[f"d_res{k}"] = round(d, 2) if d == d else ""
        dists.append(d)

    finite = [(i + 1, d) for i, d in enumerate(dists) if d == d]
    if finite:
        closest = min(finite, key=lambda t: t[1])[0]
        row["closest_res"] = closest
        row["edge_leads"] = closest <= 2
    row["status"] = "ok"
    return row


def reference_row(reg: dict) -> dict | None:
    """The same measurement on the native TIMP3:ADAM17 crystal."""
    ref = C.DATA_DIR / "TIMP3_vs_ADAM17_X_ray.pdb"
    if not ref.exists() or "TIMP3_WT" not in reg or "ADAM17" not in reg:
        return None
    base = {k: "" for k in FIELDS}
    base.update(construct_id="TIMP3_WT", target_id="ADAM17",
                source="CRYSTAL_REFERENCE", family="crystal",
                construct_src="Crystal", target_src="Crystal",
                model_path=str(ref))
    try:
        base.update(chelation_geometry(ref, reg["TIMP3_WT"]["sequence"],
                                       reg["ADAM17"]["sequence"]))
    except Exception as e:
        base["status"] = f"error:{e}"
    return base


def main() -> None:
    C.OUT_ANALYSIS.mkdir(parents=True, exist_ok=True)
    reg = {r["id"]: r for r in MR.load_registry()}
    constructs = [r for r in reg.values() if r["kind"] == "construct"]
    targets = [r for r in reg.values() if r["kind"] == "target"]
    specs = _source_specs()

    rows = []
    ref = reference_row(reg)
    if ref:
        rows.append(ref)

    for c in constructs:
        for t in targets:
            for spec in specs:
                base = {k: "" for k in FIELDS}
                base.update(construct_id=c["id"], target_id=t["id"], **spec)
                mp = _resolve(spec, c["id"], t["id"])
                if mp is None:
                    base["status"] = "pending_complex"
                    rows.append(base)
                    continue
                base["model_path"] = str(mp)
                try:
                    base.update(chelation_geometry(mp, c["sequence"], t["sequence"]))
                except Exception as e:
                    base["status"] = f"error:{e}"
                rows.append(base)

    out = C.OUT_ANALYSIS / "chelation_test.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    # ── per-source summary ────────────────────────────────────────────────
    def _f(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return float("nan")

    ok = [r for r in rows if r["status"] == "ok"]
    order = ["CRYSTAL_REFERENCE"] + [s["source"] for s in specs]
    print(f"chelation_test.csv: {len(rows)} rows ({len(ok)} ok)")
    print(f"\n{'source':<28}{'d_res1':>9}{'d_res2':>9}{'edge_leads':>12}{'n':>5}")
    print("-" * 63)
    for s in order:
        sub = [r for r in ok if r["source"] == s]
        if not sub:
            continue
        d1 = np.nanmean([_f(r["d_res1"]) for r in sub])
        d2 = np.nanmean([_f(r["d_res2"]) for r in sub])
        el = np.mean([1.0 if r["edge_leads"] is True else 0.0 for r in sub])
        print(f"{s:<28}{d1:>9.2f}{d2:>9.2f}{el:>12.2f}{len(sub):>5}")
    print(f"\n-> {out.relative_to(C.REPO_ROOT)}")


if __name__ == "__main__":
    main()
