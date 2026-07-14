"""
Complex analysis: for every construct x target pair, score the docked/co-folded
complex on the full interface battery, catalytic-zinc-loop proximity, HADDOCK
energy components, AF3 complex confidence, and DockQ against a native reference
complex where one exists (ADAM17).

Two methods per pair:
  HADDOCK   - the docked model (docking/<fold_source>/best_models/)
  AF3       - the AF3 co-fold (af3/results/), if the optional complex batches were run

Output (analysis/):
  master_complex_metrics.csv   one row per (construct, target, method)

Robust to not-yet-produced complexes (status column). Runs today.

Run:  python Structural-Validation/analysis/complex_metrics.py
"""
from __future__ import annotations

import csv
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
import config as C  # noqa: E402
import metrics as M  # noqa: E402
import model_registry as MR  # noqa: E402
import structure_io as sio  # noqa: E402


def _seq_identity(a: str, b: str) -> float:
    pairs = sio.aligned_indices(a, b)
    if not pairs:
        return 0.0
    match = sum(1 for i, j in pairs if a[i] == b[j])
    return match / max(len(a), len(b))


def assign_chains(chains: dict, construct_seq: str, target_seq: str):
    """Return (construct_chain, target_chain) by best sequence identity."""
    best_c = max(chains.values(), key=lambda c: _seq_identity(c.seq, construct_seq))
    best_t = max(chains.values(), key=lambda c: _seq_identity(c.seq, target_seq))
    if best_c.cid == best_t.cid and len(chains) >= 2:  # tie-break
        others = [c for c in chains.values() if c.cid != best_c.cid]
        best_t = max(others, key=lambda c: _seq_identity(c.seq, target_seq))
    return best_c, best_t


def _min_ca_ca(chain_a, resids_a, chain_b, resids_b) -> float:
    ca = np.array([r["CA"].coord for r in chain_a.residues if r.id[1] in set(resids_a)])
    cb = np.array([r["CA"].coord for r in chain_b.residues if r.id[1] in set(resids_b)])
    if len(ca) == 0 or len(cb) == 0:
        return float("nan")
    return float(np.min(np.linalg.norm(ca[:, None] - cb[None], axis=2)))


def parse_haddock_remarks(path: Path) -> dict:
    """Parse HADDOCK3 emref REMARK header.

    Real format (HADDOCK3):
        REMARK HADDOCK score: -84.66
        REMARK   total,bonds,angles,improper,dihe,vdw,elec,air,cdih,...
        REMARK energies: -47.2, 40.4, 193.5, 61.5, 1537.3, -66.8, -90.5, 110.1, ...
    vdW/elec/AIR are positional (indices 5/6/7). BSA/desolv/violations are not in
    the PDB header (BSA is computed independently as the `bsa` column), so they
    stay blank unless a workflow injects them.
    """
    text = path.read_text(errors="ignore")
    out = {"haddock_score": "", "vdw": "", "elec": "", "desolv": "",
           "air": "", "haddock_bsa": "", "violations": ""}
    m = re.search(r"HADDOCK score:\s*([-\d.eE+]+)", text)
    if m:
        out["haddock_score"] = float(m.group(1))
    m = re.search(r"energies:\s*([-\d.eE+,\s]+)", text)
    if m:
        vals = [v.strip() for v in m.group(1).split(",")]
        try:
            out["vdw"], out["elec"], out["air"] = (
                float(vals[5]), float(vals[6]), float(vals[7]))
        except (IndexError, ValueError):
            pass
    b = re.search(r"[Bb]uried [Ss]urface [Aa]rea:\s*([-\d.]+)", text)
    if b:
        out["haddock_bsa"] = float(b.group(1))
    v = re.search(r"[Vv]iolations[^:]*:\s*([-\d.]+)", text)
    if v:
        out["violations"] = float(v.group(1))
    return out


def af3_complex_conf(path: Path) -> dict:
    for j in path.parent.glob("*summary_confidences*.json"):
        try:
            d = json.loads(j.read_text())
            return {"iptm": d.get("iptm", ""), "ptm": d.get("ptm", ""),
                    "pae": (np.mean(d["pae"]) if "pae" in d else "")}
        except Exception:
            pass
    return {"iptm": "", "ptm": "", "pae": ""}


FIELDS = ["construct_id", "target_id", "source", "family", "construct_src",
          "target_src", "status",
          "bsa", "n_iface_res_construct", "n_iface_res_target", "n_contacts_5A",
          "n_hbonds", "n_salt_bridges", "n_hydrophobic", "contact_density",
          "min_ca_ca_zincloop", "haddock_score", "vdw", "elec", "desolv", "air",
          "haddock_bsa", "violations", "iptm", "ptm", "pae", "esm_lplddt",
          "dockq", "fnat", "irms", "lrms", "capri", "dockq_ref_type",
          "model_path", "ref_path"]

TIMP3_ACTIVE_N = C.TIMP3_ACTIVE[1]  # first N residues of the construct = reactive edge


def _esm_complex_conf() -> dict[tuple[str, str], dict]:
    p = C.OUT_ESM / "esmfold2_complex_metrics.csv"
    if not p.exists():
        return {}
    with open(p, newline="") as fh:
        return {(r["construct_id"], r["target_id"]): r for r in csv.DictReader(fh)}


def score_complex(path, construct_seq, target_seq, family, ref_path,
                  ref_type=None, esm_conf=None) -> dict:
    chains = sio.get_chains(path)
    if len(chains) < 2:
        return {"status": "single_chain_only"}
    cchain, tchain = assign_chains(chains, construct_seq, target_seq)
    cid_pdb, tid_pdb = cchain.cid, tchain.cid

    im = M.interface_summary(path, cid_pdb, tid_pdb)
    row = {
        "bsa": im.bsa, "n_iface_res_construct": im.n_iface_res_A,
        "n_iface_res_target": im.n_iface_res_B, "n_contacts_5A": im.n_contacts_5A,
        "n_hbonds": im.n_hbonds, "n_salt_bridges": im.n_salt_bridges,
        "n_hydrophobic": im.n_hydrophobic, "contact_density": im.contact_density,
    }

    timp_edge = cchain.resids[:TIMP3_ACTIVE_N]
    motif_res = sio.zinc_motif_resids(tchain)
    row["min_ca_ca_zincloop"] = round(_min_ca_ca(cchain, timp_edge, tchain, motif_res), 2)

    if family == "HADDOCK":
        row.update(parse_haddock_remarks(Path(path)))
    elif family == "AF3_cofold":
        row.update(af3_complex_conf(Path(path)))
    elif family == "ESMFold2_cofold" and esm_conf:
        row.update(iptm=esm_conf.get("esm_iptm", ""), ptm=esm_conf.get("esm_ptm", ""),
                   pae=esm_conf.get("esm_pae", ""), esm_lplddt=esm_conf.get("esm_lplddt", ""))

    if ref_path:
        try:
            rc = sio.get_chains(ref_path)
            rcc, rtc = assign_chains(rc, construct_seq, target_seq)
            dq = M.dockq(path, ref_path, tid_pdb, cid_pdb, rtc.cid, rcc.cid)
            row.update(dockq=dq.dockq, fnat=dq.fnat, irms=dq.irms,
                       lrms=dq.lrms, capri=dq.capri, dockq_ref_type=ref_type)
        except Exception as e:
            row["capri"] = f"dockq_error:{e}"
    row["status"] = "ok"
    return row


def _source_specs() -> list[dict]:
    """Every complex model source scored: docking tracks + native co-folds."""
    specs = []
    for cs, ts in C.DOCK_MATRIX["full"]:
        specs.append({"source": f"HADDOCK:{cs}x{ts}", "family": "HADDOCK",
                      "construct_src": cs, "target_src": ts})
    specs.append({"source": "AF3_cofold", "family": "AF3_cofold",
                  "construct_src": "AF3", "target_src": "AF3"})
    specs.append({"source": "ESMFold2_cofold", "family": "ESMFold2_cofold",
                  "construct_src": "ESMFold2", "target_src": "ESMFold2"})
    return specs


def _resolve(spec, cid, tid) -> Path | None:
    if spec["family"] == "HADDOCK":
        return MR.haddock_complex(cid, tid, spec["construct_src"], spec["target_src"])
    if spec["family"] == "AF3_cofold":
        return MR.af3_complex(cid, tid)
    return MR.esmfold_complex(cid, tid)


def main() -> None:
    C.OUT_ANALYSIS.mkdir(parents=True, exist_ok=True)
    reg = {r["id"]: r for r in MR.load_registry()}
    constructs = [r for r in reg.values() if r["kind"] == "construct"]
    targets = [r for r in reg.values() if r["kind"] == "target"]
    specs = _source_specs()
    esm_conf = _esm_complex_conf()

    rows = []
    for c in constructs:
        for t in targets:
            ref, ref_type = MR.complex_reference(t["id"])
            for spec in specs:
                base = {k: "" for k in FIELDS}
                base.update(construct_id=c["id"], target_id=t["id"], **spec,
                            ref_path=str(ref) if ref else "")
                mp = _resolve(spec, c["id"], t["id"])
                if mp is None:
                    base["status"] = "pending_complex"
                    rows.append(base)
                    continue
                base["model_path"] = str(mp)
                try:
                    ec = esm_conf.get((c["id"], t["id"])) if spec["family"] == "ESMFold2_cofold" else None
                    base.update(score_complex(mp, c["sequence"], t["sequence"],
                                              spec["family"], ref, ref_type, ec))
                except Exception as e:
                    base["status"] = f"error:{e}"
                rows.append(base)

    with open(C.OUT_ANALYSIS / "master_complex_metrics.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    ok = sum(1 for r in rows if r["status"] == "ok")
    pend = sum(1 for r in rows if r["status"] == "pending_complex")
    print(f"master_complex_metrics.csv: {len(rows)} rows ({ok} scored, {pend} pending)")
    print(f"  {len(constructs)} constructs x {len(targets)} targets x {len(specs)} sources: "
          f"{[s['source'] for s in specs]}")
    print(f"-> {(C.OUT_ANALYSIS / 'master_complex_metrics.csv').relative_to(C.REPO_ROOT)}")


if __name__ == "__main__":
    main()
