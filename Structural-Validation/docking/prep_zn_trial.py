"""
Zinc-aware HADDOCK trial — input preparation.

Docks INDEPENDENTLY-folded AF3 monomers (unbound conformations), matching the prior
pipeline, so HADDOCK has to find the pose from scratch. A first version split the
AF3 CO-FOLD into two chains — that fed HADDOCK the pre-formed bound conformations
and biased the result (ADAM17 noZn DockQ 0.78 vs the pipeline's 0.05); it is
replaced here.

Protocol = legacy HADDOCK_Outputs (HADDOCK3: rigidbody 500 -> seletop 200 ->
flexref -> emref -> clustfcc -> caprieval). Two arms differ by ONE variable:

  noZn : TIMP3 monomer + target monomer, as folded          (control = prior setup)
  Zn   : identical, except the target carries its catalytic Zn(2+)

Inputs:
  TIMP3  = MR.af3_monomer("TIMP3_WT")   (independent fold, unbound)
  target = MR.af3_monomer(<target>)     (independent fold, unbound)
so the only difference between arms is the metal.

The catalytic zinc is identified in the matching target crystal as the Zn nearest
the HExxHxxGxxH histidine triad (MMPs carry a second, structural zinc), then
transplanted onto the AF target by superposing the crystal on the AF model.

Restraint policy — deliberate:
  * AIRs are the legacy ones: target zinc motif <-> TIMP3 reactive edge (1-5) /
    passive (6-10). Identical in both arms.
  * In the Zn arm we additionally restrain the Zn to its OWN three histidines so
    CNS cannot drift the metal out of the site during flexref. That is
    target-internal geometry.
  * We do NOT restrain TIMP3 Cys1 to the zinc. Whether the reactive edge finds the
    metal is the hypothesis under test; restraining it would assume the answer.

Run:  python Structural-Validation/docking/prep_zn_trial.py [--check]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO, Select

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis"))
import config as C  # noqa: E402
import metrics as M  # noqa: E402
import model_registry as MR  # noqa: E402
import structure_io as sio  # noqa: E402
from complex_metrics import _seq_identity  # noqa: E402

TRIAL = C.OUT_DOCK / "zn_trial"          # scratch inputs/configs live under Local/
ARMS = ("noZn", "Zn")
TIMP_ENTITY = "TIMP3_WT"

# Panel targets with both an independent AF3 monomer fold and a target crystal.
TARGETS = ["ADAM10", "ADAM17", "MMP2", "MMP3", "MMP9"]
CRYSTAL = {"ADAM10": "ADAM10_Xray.pdb", "ADAM17": "ADAM17_Xray.pdb",
           "MMP2": "MMP2_Xray.pdb", "MMP3": "MMP3cd_X_ray.pdb",
           "MMP9": "MMP9_Xray.pdb"}
ZN_ELEMENTS = {"ZN"}


def _force_chain(path: Path, cid: str) -> None:
    """Rewrite the chain-ID column (22) of every coordinate record."""
    out = []
    for ln in path.read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")) and len(ln) > 21:
            ln = ln[:21] + cid + ln[22:]
        out.append(ln)
    path.write_text("\n".join(out) + "\n")


class _ChainSel(Select):
    def __init__(self, cid):
        self.cid = cid

    def accept_chain(self, chain):
        return chain.id == self.cid

    def accept_residue(self, residue):
        return residue.id[0] == " "          # standard amino acids only


def _save_monomer(struct_path, src_cid: str, new_cid: str, out_pdb: Path) -> None:
    """Write one protein chain of a fold (CIF/PDB) to PDB, relabelled to new_cid."""
    st = sio.load(struct_path)
    io = PDBIO()
    io.set_structure(st)
    io.save(str(out_pdb), _ChainSel(src_cid))
    _force_chain(out_pdb, new_cid)


def _write_ensemble(conformers, cid: str, out_pdb: Path) -> int:
    """Multi-MODEL ensemble PDB from several folds of the same chain.

    Conformers are superposed onto the first before writing. AF3's own 5 seeds are
    near-identical (0.2-0.8 A), so a useful ensemble needs CROSS-METHOD members
    (AF3 + ESMFold2), which differ by ~0.9 A globally and up to 2.6 A at the
    reactive edge — the conformational freedom rigid-body docking otherwise lacks.
    """
    blocks, ref = [], None
    for i, p in enumerate(conformers):
        st = sio.load(p)
        ch = max(sio.get_chains(p).values(), key=lambda c: len(c.seq))
        if ref is None:
            ref = ch
        else:
            m, r = M._matched_ca(ch, ref)
            if len(m) >= 3:
                R, t, _ = M.kabsch(m, r)
                for a in st.get_atoms():
                    a.coord = M.apply_transform(np.asarray([a.coord], float), R, t)[0]
        tmp = out_pdb.with_name(out_pdb.name + f".m{i}.tmp")
        io = PDBIO()
        io.set_structure(st)
        io.save(str(tmp), _ChainSel(ch.cid))
        _force_chain(tmp, cid)
        blocks.append([ln for ln in tmp.read_text().splitlines()
                       if ln.startswith(("ATOM", "TER"))])
        tmp.unlink()
    with open(out_pdb, "w") as f:
        for i, b in enumerate(blocks, 1):
            f.write(f"MODEL     {i}\n")
            f.write("\n".join(b) + "\n")
            f.write("ENDMDL\n")
        f.write("END\n")
    return len(blocks)


def _zn_atoms(structure):
    """All zinc atoms.

    Cannot trust the parsed element: several of these crystals write the element
    field one column left of spec, so BioPython falls back to guessing from the
    atom name "ZN+2" and returns 'N'. Residue name (ZN / ZN1 / ZN2) and atom name
    both start with ZN in every file we have, so key off those as well.
    """
    out = []
    for a in structure.get_atoms():
        if a.get_parent().id[0] == " ":          # standard residue, not a hetero ion
            continue
        el = (a.element or "").strip().upper()
        resn = a.get_parent().get_resname().strip().upper()
        name = a.get_name().strip().upper()
        if el in ZN_ELEMENTS or resn.startswith("ZN") or name.startswith("ZN"):
            out.append(a)
    return out


def catalytic_zn(crystal_path: Path, motif_resids: list[int], chain_id: str):
    """The Zn nearest the motif histidine triad (MMPs also carry a structural Zn)."""
    st = sio.load(crystal_path)
    zns = _zn_atoms(st)
    if not zns:
        return None, None
    model = list(st)[0]
    ch = model[chain_id] if chain_id in [c.id for c in model] else None
    if ch is None:
        return None, None
    his = []
    for r in ch:
        if r.id[1] in set(motif_resids) and r.get_resname() == "HIS":
            for an in ("NE2", "ND1"):
                if an in r:
                    his.append(r[an].coord)
    if not his:
        return None, None
    hc = np.asarray(his, float)
    best, bestd = None, 1e9
    for z in zns:
        d = float(np.linalg.norm(hc - z.coord, axis=1).min())
        if d < bestd:
            best, bestd = z, d
    return best, bestd


def prep_target(target: str, out_dir: Path, timp_af, tc_ref_seq=None) -> dict:
    """Prepare independent-monomer inputs + Zn transplant for one target.

    timp_af = path to the independently-folded TIMP3 monomer (shared across targets);
    the target monomer is its own independent AF3 fold. No co-fold is used, so the
    partners arrive in UNBOUND conformations and HADDOCK must dock from scratch.
    """
    info = {"target": target, "status": "ok"}
    tgt_af = MR.af3_monomer(target)
    if tgt_af is None or timp_af is None:
        return {"target": target, "status": "no_af_monomer"}
    info["timp_input"], info["target_input"] = Path(timp_af).name, Path(tgt_af).name

    cchains, tchains = sio.get_chains(timp_af), sio.get_chains(tgt_af)
    if not cchains or not tchains:
        return {"target": target, "status": "unreadable_monomer"}
    cc = max(cchains.values(), key=lambda c: len(c.seq))   # TIMP3 (one chain)
    tc = max(tchains.values(), key=lambda c: len(c.seq))   # target (one chain)
    info["timp_len"], info["target_len"] = len(cc.seq), len(tc.seq)

    motif = sio.zinc_motif_resids(tc)
    info["n_motif"] = len(motif)
    if len(motif) < 11:
        info["status"] = "no_zinc_motif"
        return info
    info["motif_resids"] = motif

    out_dir.mkdir(parents=True, exist_ok=True)
    timp_pdb = out_dir / f"{target}_TIMP3.pdb"
    tgt_pdb = out_dir / f"{target}_target.pdb"
    # AIR convention: target = segid A, TIMP3 = segid B.
    _save_monomer(tgt_af, tc.cid, "A", tgt_pdb)
    _save_monomer(timp_af, cc.cid, "B", timp_pdb)
    info["timp_pdb"], info["target_pdb"] = str(timp_pdb), str(tgt_pdb)

    # ── transplant the catalytic Zn from the crystal onto the AF target ──
    cx = C.CRYSTAL_DIR / CRYSTAL[target]
    if not cx.exists():
        info["zn_status"] = f"no_crystal:{CRYSTAL[target]}"
        return info
    xchains = sio.get_chains(cx)
    if not xchains:
        info["zn_status"] = "crystal_unreadable"
        return info
    xt = max(xchains.values(), key=lambda c: _seq_identity(c.seq, tc.seq))
    xmotif = sio.zinc_motif_resids(xt)
    if len(xmotif) < 11:
        info["zn_status"] = "crystal_no_motif"
        return info
    zn, dhis = catalytic_zn(cx, xmotif, xt.cid)
    if zn is None:
        info["zn_status"] = "no_zn_in_crystal"
        return info
    info["zn_his_dist_crystal"] = round(dhis, 2)

    mob, ref = M._matched_ca(xt, tc)          # crystal -> AF frame
    if len(mob) < 3:
        info["zn_status"] = "superpose_failed"
        return info
    R, t, rms = M.kabsch(mob, ref)
    info["zn_superpose_rmsd"] = round(rms, 2)
    zn_af = M.apply_transform(np.asarray([zn.coord], float), R, t)[0]

    # sanity: the transplanted Zn must land in the AF model's own motif His triad
    af_his = []
    for r in tc.residues:
        if r.id[1] in set(motif) and r.get_resname() == "HIS":
            for an in ("NE2", "ND1"):
                if an in r:
                    af_his.append(r[an].coord)
    d_af = (float(np.linalg.norm(np.asarray(af_his, float) - zn_af, axis=1).min())
            if af_his else float("nan"))
    info["zn_his_dist_af"] = round(d_af, 2)
    info["zn_status"] = "ok" if d_af < 3.5 else f"suspect_geometry({d_af:.2f}A)"

    # write the Zn-bearing target: HADDOCK topology residue name ZN2 = zinc 2+
    lines = [ln for ln in tgt_pdb.read_text().splitlines()
             if not ln.startswith(("END", "MASTER", "CONECT"))]
    # Serial and resSeq must follow on from the file: a hard-coded 99999 makes CNS
    # abort with "Structure contains more than 99.999 atoms".
    serials, resids = [0], [0]
    for ln in lines:
        if ln.startswith(("ATOM", "HETATM")) and len(ln) > 26:
            try:
                serials.append(int(ln[6:11]))
                resids.append(int(ln[22:26]))
            except ValueError:
                pass
    serial, zn_resid = max(serials) + 1, max(resids) + 1
    lines.append(
        f"HETATM{serial:>5}  ZN  ZN2 A{zn_resid:>4}    "
        f"{zn_af[0]:8.3f}{zn_af[1]:8.3f}{zn_af[2]:8.3f}  1.00  0.00          ZN")
    lines.append("END")
    tgt_zn = out_dir / f"{target}_target_Zn.pdb"
    tgt_zn.write_text("\n".join(lines) + "\n")
    info["target_pdb_zn"] = str(tgt_zn)
    # what the (stdlib-only, WSL-side) runner needs to write restraints
    info["zn_resid"] = zn_resid
    info["zn_his_resids"] = [motif[i] for i in (0, 4, 10)]
    info["timp_active"] = cc.resids[:5]      # C1-T2-C3-S4-P5 reactive edge
    info["timp_passive"] = cc.resids[5:10]

    # ── ensemble inputs: AF3 seeds + the ESMFold2 fold (cross-method diversity) ──
    def _members(entity, af_path):
        out = [af_path]
        alt = Path(af_path).parent / Path(af_path).name.replace("model_0", "model_4")
        if alt.exists():
            out.append(alt)
        esm = MR.esmfold_monomer(entity)
        if esm:
            out.append(esm)
        return out
    try:
        timp_ens = out_dir / f"{target}_TIMP3_ens.pdb"
        tgt_ens = out_dir / f"{target}_target_ens.pdb"
        n_c = _write_ensemble(_members(TIMP_ENTITY, timp_af), "B", timp_ens)
        n_t = _write_ensemble(_members(target, tgt_af), "A", tgt_ens)
        info["timp_pdb_ens"], info["target_pdb_ens"] = str(timp_ens), str(tgt_ens)
        info["n_conformers"] = f"{n_c}x{n_t}"
        # semi-flexible segments for flexref: target zinc loop, TIMP3 reactive edge
        info["flex_target"] = [min(motif) - 2, max(motif) + 2]
        info["flex_timp"] = [cc.resids[0], cc.resids[min(11, len(cc.resids) - 1)]]
    except Exception as e:
        info["ens_status"] = f"error:{e}"
    return info


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="prepare and report only (no configs written)")
    args = ap.parse_args()

    timp_af = MR.af3_monomer(TIMP_ENTITY)
    if timp_af is None:
        print(f"ERROR: no independent AF3 monomer for {TIMP_ENTITY}")
        return
    print(f"TIMP3 monomer: {Path(timp_af).name}")
    out_dir = TRIAL / "inputs"

    rows = []
    for t in TARGETS:
        try:
            rows.append(prep_target(t, out_dir, timp_af))
        except Exception as e:
            rows.append({"target": t, "status": f"error:{e}"})

    hdr = (f"{'target':<9}{'status':<16}{'timp':>5}{'tgt':>5}{'motif':>7}"
           f"{'supRMSD':>9}{'Zn-His(x)':>11}{'Zn-His(AF)':>12}  zn_status")
    print(hdr); print("-" * len(hdr))
    for r in rows:
        print(f"{r.get('target',''):<9}{r.get('status',''):<16}"
              f"{r.get('timp_len',''):>5}{r.get('target_len',''):>5}"
              f"{r.get('n_motif',''):>7}{r.get('zn_superpose_rmsd',''):>9}"
              f"{r.get('zn_his_dist_crystal',''):>11}{r.get('zn_his_dist_af',''):>12}"
              f"  {r.get('zn_status','')}")
    ok = [r for r in rows if r.get("zn_status") == "ok"]
    manifest = TRIAL / "manifest.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(
        {"targets": rows, "arms": list(ARMS),
         "motif_note": "AIRs identical in both arms; Zn arm adds Zn<->its own His only"},
        indent=2))
    print(f"\n{len(ok)}/{len(rows)} targets ready with a validated catalytic Zn")
    print(f"inputs   -> {out_dir}")
    print(f"manifest -> {manifest}")


if __name__ == "__main__":
    main()
