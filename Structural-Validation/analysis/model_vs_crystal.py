"""
Monomer analysis: for every construct/target and every fold modeler (AF3,
ESMFold2), score the model against its experimental crystal structure and on
model-only quality, then quantify AF3-vs-ESMFold2 agreement.

Outputs (analysis/):
  master_monomer_metrics.csv   one row per (entity, folder)
  af3_vs_esmfold_monomer.csv    one row per entity (cross-modeler agreement)

Robust to not-yet-produced models: missing outputs are recorded with a status
so the table is complete and the run never fails. Runs today; fills in as
folding results arrive.

Run:  python Structural-Validation/analysis/model_vs_crystal.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
import config as C  # noqa: E402
import metrics as M  # noqa: E402
import model_registry as MR  # noqa: E402
import structure_io as sio  # noqa: E402


def _esmfold_conf() -> dict[str, dict]:
    p = C.OUT_ESM / "esmfold2_metrics.csv"
    if not p.exists():
        return {}
    with open(p, newline="") as fh:
        return {r["id"]: r for r in csv.DictReader(fh)}


def _af3_confidence(model_path: Path) -> tuple[float, float]:
    """(plddt_mean, ptm) for an AF3 model: pLDDT from CA b-factors, pTM from json."""
    chains = sio.get_chains(model_path)
    ch = sio.longest_chain(chains)
    bfac = [r["CA"].bfactor for r in ch.residues if "CA" in r]
    plddt = sum(bfac) / len(bfac) if bfac else float("nan")
    ptm = float("nan")
    for j in model_path.parent.glob("*summary_confidences*.json"):
        try:
            data = json.loads(j.read_text())
            ptm = float(data.get("ptm", ptm))
            break
        except Exception:
            pass
    return plddt, ptm


def _model_only(chain) -> dict:
    return {
        "rg_model": round(M.radius_of_gyration(chain), 2),
        "clashscore": round(M.clashscore(chain), 2),
        "rama_outlier_frac": round(M.ramachandran_outlier_frac(chain), 3),
    }


def _vs_crystal(model_chain, ref_chain) -> dict:
    sup = sio.superpose_ca(model_chain, ref_chain)
    return {
        "n_common": sup.n_common,
        "rmsd_ca": round(sup.rmsd, 3),
        "tm_score": round(M.tm_score(model_chain, ref_chain), 3),
        "gdt_ts": round(M.gdt_ts(model_chain, ref_chain), 1),
        "gdt_ha": round(M.gdt_ha(model_chain, ref_chain), 1),
        "lddt": round(M.lddt(model_chain, ref_chain), 1),
        "rg_ref": round(M.radius_of_gyration(ref_chain), 2),
    }


FIELDS = ["entity_id", "kind", "folder", "status", "n_common", "rmsd_ca",
          "tm_score", "gdt_ts", "gdt_ha", "lddt", "rg_model", "rg_ref",
          "clashscore", "rama_outlier_frac", "plddt_mean", "ptm",
          "model_path", "ref_path"]


def main() -> None:
    C.OUT_ANALYSIS.mkdir(parents=True, exist_ok=True)
    reg = MR.load_registry()
    esm_conf = _esmfold_conf()

    rows = []
    cache: dict[tuple[str, str], object] = {}   # (entity, folder) -> chain
    for r in reg:
        eid, kind = r["id"], r["kind"]
        ref_path = MR.crystal_reference(eid)
        ref_chain = None
        if ref_path:
            try:
                ref_chain = sio.longest_chain(sio.get_chains(ref_path))
            except Exception as e:
                ref_path = f"PARSE_ERROR:{e}"

        for folder in C.FOLDERS:
            row = {k: "" for k in FIELDS}
            row.update(entity_id=eid, kind=kind, folder=folder,
                       ref_path=str(ref_path) if ref_path else "")
            mp = MR.monomer_model(eid, folder)
            if mp is None:
                row["status"] = "pending_model"
                rows.append(row)
                continue
            try:
                mchain = sio.longest_chain(sio.get_chains(mp))
            except Exception as e:
                row["status"] = f"model_parse_error:{e}"
                rows.append(row)
                continue
            cache[(eid, folder)] = mchain
            row["model_path"] = str(mp)
            row.update(_model_only(mchain))

            # confidence
            if folder == "ESMFold2":
                cf = esm_conf.get(eid, {})
                row["plddt_mean"] = cf.get("plddt_mean", "")
                row["ptm"] = cf.get("ptm", "")
            else:
                pl, pt = _af3_confidence(mp)
                row["plddt_mean"] = round(pl, 2)
                row["ptm"] = "" if pt != pt else round(pt, 4)

            if ref_chain is not None:
                try:
                    row.update(_vs_crystal(mchain, ref_chain))
                    row["status"] = "ok"
                except Exception as e:
                    row["status"] = f"align_error:{e}"
            else:
                row["status"] = "no_crystal_ref"
            rows.append(row)

    with open(C.OUT_ANALYSIS / "master_monomer_metrics.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    # AF3 vs ESMFold2 agreement (structural, independent of any crystal)
    agree = []
    for r in reg:
        eid = r["id"]
        a = cache.get((eid, "AF3"))
        e = cache.get((eid, "ESMFold2"))
        rec = {"entity_id": eid, "kind": r["kind"], "status": "pending"}
        if a is not None and e is not None:
            try:
                sup = sio.superpose_ca(e, a)  # move ESMFold onto AF3
                rec.update(status="ok", n_common=sup.n_common,
                           rmsd_af3_esm=round(sup.rmsd, 3),
                           tm_af3_esm=round(M.tm_score(e, a), 3))
            except Exception as ex:
                rec["status"] = f"error:{ex}"
        agree.append(rec)

    afields = ["entity_id", "kind", "status", "n_common", "rmsd_af3_esm", "tm_af3_esm"]
    with open(C.OUT_ANALYSIS / "af3_vs_esmfold_monomer.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=afields)
        w.writeheader()
        for rec in agree:
            w.writerow({k: rec.get(k, "") for k in afields})

    done = sum(1 for r in rows if r["status"] == "ok")
    pend = sum(1 for r in rows if r["status"] == "pending_model")
    print(f"master_monomer_metrics.csv: {len(rows)} rows "
          f"({done} scored vs crystal, {pend} awaiting models)")
    print(f"af3_vs_esmfold_monomer.csv: {len(agree)} entities")
    print(f"-> {C.OUT_ANALYSIS.relative_to(C.REPO_ROOT)}")


if __name__ == "__main__":
    main()
