#!/usr/bin/env python3
"""
conformation_seed.py — carry a design's LOOP CONFORMATION forward without carrying
its scaffold drift.

Why this exists
---------------
Measured 2026-08-17 on iterations 46-54: ~89% of sv_pdockq variance is set by the
loop BACKBONE conformation (RFd3) and only ~11% by the sequence on it (LigandMPNN).
The refinement loop only ever re-rolled the sequence, so nine rounds moved the best
design by 0.001. To improve we must carry conformation forward.

But the predicted structures are NOT faithful scaffolds: median CA RMSD of the
constant C-terminal 92 residues against the AF3 template is 4.3 A (up to 9.7), and
the supposedly-fixed MMP2 target chain drifts 1.6-3.8 A. Re-templating on a raw
prediction would compound that every round.

So: graft the loop coordinates onto the CANONICAL template scaffold and target.
The scaffold and target are reset to the template every round, which makes drift
impossible by construction rather than something we monitor and hope stays small.
The loops -- the only part carrying signal -- are inherited.

    seed = build_seed_template(design_cif, template_pdb, ["AB","C","EF"], out_pdb)
    scaf, tgt = drift_rmsd(design_cif, template_pdb)      # seed-eligibility gate
"""
import re
import sys
from pathlib import Path

import numpy as np
from biotite.structure import superimpose, rmsd as bio_rmsd
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile
import biotite.structure.io.pdbx as pdbx
from biotite.sequence import ProteinSequence

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

BINDER_CHAIN, TARGET_CHAIN = "A", "B"

# Flank tripeptides that bracket each redesigned loop. Kept in step with
# iterative_refinement.LOOP_CONFIGS; imported lazily so this module stays usable
# without atomworks.
_FLANKS = {
    "AB": ("LVK", "LVY"),
    "C":  ("HTE", "GLK"),
    "EF": ("MYT", "FVE"),
    "GH": ("KSC", "NEC"),
    "MTL": ("LWT", "YQS"),
}


def _chain(arr, ch):
    return arr[arr.chain_id == ch]


def _seq_and_resids(arr, ch):
    """One-letter sequence and the residue ids, in order, for one chain."""
    a = _chain(arr, ch)
    ca = a[a.atom_name == "CA"]
    seq = "".join(ProteinSequence.convert_letter_3to1(r) for r in ca.res_name)
    return seq, ca.res_id


def loop_residue_ids(arr, ch, loops):
    """
    {loop: array of res_id} for the designed loops, located by flanking tripeptides
    so it works regardless of how long each loop turned out.
    """
    seq, resids = _seq_and_resids(arr, ch)
    out, cursor = {}, 0
    for name in loops:
        left, right = _FLANKS[name]
        m = re.compile(f"{left}([A-Z]*?){right}").search(seq[cursor:])
        if not m:
            out[name] = np.array([], dtype=int)
            continue
        start = cursor + m.start() + len(left)          # 0-indexed loop start
        out[name] = resids[start:start + len(m.group(1))]
        cursor = cursor + m.end() - len(right)
    return out


def _split_const_loop(arr, ch, loops):
    """(constant res_ids, loop res_ids) for one chain."""
    loop_ids = loop_residue_ids(arr, ch, loops)
    loop_set = set(int(i) for ids in loop_ids.values() for i in ids)
    _, resids = _seq_and_resids(arr, ch)
    const = np.array([int(i) for i in resids if int(i) not in loop_set])
    return const, loop_ids


def _res_mask(arr, ch, resids):
    return (arr.chain_id == ch) & np.isin(arr.res_id, resids)


def _ca_of(arr, ch, resids):
    m = _res_mask(arr, ch, resids) & (arr.atom_name == "CA")
    sub = arr[m]
    order = np.argsort(sub.res_id)
    return sub[order]


def _load(path):
    p = str(path)
    if p.lower().endswith((".cif", ".mmcif")):
        return pdbx.get_structure(CIFFile.read(p), model=1)
    return PDBFile.read(p).get_structure()[0]


def drift_rmsd(design_path, template_path, loops=("AB", "C", "EF")):
    """
    (scaffold_rmsd, target_rmsd) in Angstrom for a predicted complex vs the template.

    scaffold = the binder's CONSTANT residues (everything outside the designed loops),
    paired 1:1 with the template's constant residues -- their counts match by
    construction even though loop lengths differ.
    target   = the fixed target chain, which should barely move at all.
    """
    des, tpl = _load(design_path), _load(template_path)
    d_const, _ = _split_const_loop(des, BINDER_CHAIN, loops)
    t_const, _ = _split_const_loop(tpl, BINDER_CHAIN, loops)
    n = min(len(d_const), len(t_const))
    dca, tca = _ca_of(des, BINDER_CHAIN, d_const[:n]), _ca_of(tpl, BINDER_CHAIN, t_const[:n])
    k = min(len(dca), len(tca))
    fit, _ = superimpose(tca[:k], dca[:k])
    scaf = float(bio_rmsd(tca[:k], fit))

    dB, tB = _chain(des, TARGET_CHAIN), _chain(tpl, TARGET_CHAIN)
    dB, tB = dB[dB.atom_name == "CA"], tB[tB.atom_name == "CA"]
    m = min(len(dB), len(tB))
    fitB, _ = superimpose(tB[:m], dB[:m])
    return scaf, float(bio_rmsd(tB[:m], fitB))


def build_seed_template(design_path, template_path, loops=("AB", "C", "EF"),
                        out_path=None):
    """
    Canonical scaffold + canonical target + THIS DESIGN'S loop conformations.

    The design is first superimposed onto the template through their paired constant
    binder residues, which puts its loops into the template's frame. The output chain
    A is then assembled residue by residue in sequence order: template coordinates
    for constant positions, design coordinates for loop positions. Chain B is copied
    from the template verbatim.

    Returns the assembled AtomArray, and writes a PDB if out_path is given.
    """
    des, tpl = _load(design_path), _load(template_path)
    d_const, d_loops = _split_const_loop(des, BINDER_CHAIN, loops)
    t_const, t_loops = _split_const_loop(tpl, BINDER_CHAIN, loops)
    n = min(len(d_const), len(t_const))
    if n < 20:
        raise ValueError(f"only {n} paired constant residues; loop detection failed")

    # Put the design into the template frame using the constant scaffold only.
    dca = _ca_of(des, BINDER_CHAIN, d_const[:n])
    tca = _ca_of(tpl, BINDER_CHAIN, t_const[:n])
    k = min(len(dca), len(tca))
    _, transform = superimpose(tca[:k], dca[:k])
    des_fit = transform.apply(des)

    d_loop_set = {int(i) for ids in d_loops.values() for i in ids}
    _, d_resids = _seq_and_resids(des, BINDER_CHAIN)
    t_const_iter = iter(t_const)

    pieces, new_id = [], 1
    for rid in d_resids:
        rid = int(rid)
        if rid in d_loop_set:                      # inherited loop: design coords
            sub = des_fit[_res_mask(des_fit, BINDER_CHAIN, [rid])]
        else:                                       # canonical scaffold: template coords
            try:
                t_rid = int(next(t_const_iter))
            except StopIteration:
                break
            sub = tpl[_res_mask(tpl, BINDER_CHAIN, [t_rid])]
        sub = sub.copy()
        sub.res_id[:] = new_id
        sub.chain_id[:] = BINDER_CHAIN
        pieces.append(sub)
        new_id += 1

    tgt = _chain(tpl, TARGET_CHAIN).copy()
    seed = pieces[0]
    for p in pieces[1:]:
        seed = seed + p
    seed = seed + tgt

    if out_path:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        f = PDBFile()
        f.set_structure(seed)
        f.write(str(out_path))
    return seed


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Graft a design's loops onto the canonical template.")
    ap.add_argument("design"); ap.add_argument("template")
    ap.add_argument("--out", default=None)
    ap.add_argument("--loops", nargs="+", default=["AB", "C", "EF"])
    a = ap.parse_args()
    scaf, tgt = drift_rmsd(a.design, a.template, a.loops)
    print(f"design drift vs template: scaffold {scaf:.2f} A | target {tgt:.2f} A")
    seed = build_seed_template(a.design, a.template, a.loops, a.out)
    print(f"seed built: {seed.array_length()} atoms, "
          f"{len(np.unique(seed.res_id[seed.chain_id == 'A']))} binder residues")
    if a.out:
        print(f"wrote {a.out}")
