"""
Extract clean, catalytic-domain-only, single-chain target structures from the
crystal files, so docking uses a consistent ground-truth target and monomer
comparisons align cleanly even when a crystal is a multi-domain / pro-form
structure (MMP2, ADAM10) or carries other chains, ions, or waters.

For each target:
  * pick the protein chain best matching the registry CD sequence,
  * keep the aligned catalytic-domain span (+ a small margin), drop the rest,
  * strip heteroatoms/waters/other chains,
  * write crystal/clean_targets/<TARGET>.pdb (chain A).

Also writes a clean TIMP3 (chain B) and, for ADAM17, a clean two-chain
reference complex (target A + TIMP3 B) for consistent DockQ.

Run:  python Structural-Validation/crystal/extract_targets.py
"""
from __future__ import annotations

import csv
import sys
import warnings
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
import config as C  # noqa: E402
import structure_io as sio  # noqa: E402

warnings.filterwarnings("ignore")
MARGIN = 5  # residues kept on each side of the aligned CD span


def _parser(path: Path):
    return MMCIFParser(QUIET=True) if path.suffix.lower() == ".cif" else PDBParser(QUIET=True)


class _KeepSpan(Select):
    def __init__(self, chain_id, keep_resids):
        self.chain_id = chain_id
        self.keep = set(keep_resids)

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return is_aa(residue, standard=True) and residue.id[1] in self.keep

    def accept_atom(self, atom):
        return atom.element != "H"


def _best_chain_for(path: Path, ref_seq: str):
    """Return (chain_obj, chain_helper) of the chain best matching ref_seq."""
    chains = sio.get_chains(path)

    def ident(ch):
        pairs = sio.aligned_indices(ch.seq, ref_seq)
        if not pairs:
            return 0.0
        return sum(1 for i, j in pairs if ch.seq[i] == ref_seq[j]) / len(ref_seq)

    best = max(chains.values(), key=ident)
    return best, ident(best)


def _aligned_span_resids(chain, ref_seq: str) -> list[int]:
    """Residue numbers matched to the CD reference — keeps the catalytic domain
    only, dropping mid-domain inserts (MMP2 fibronectin-II) and flanking domains
    (ADAM10 ectodomain) that the purchased construct does not contain."""
    pairs = sio.aligned_indices(chain.seq, ref_seq)
    if not pairs:
        return chain.resids
    return [chain.resids[i] for i, _ in pairs]


def extract_target(target_id: str, out_dir: Path) -> dict:
    meta = C.TARGETS[target_id]
    src = C.CRYSTAL_DIR / meta["crystal"]
    reg = {r["id"]: r for r in csv.DictReader(open(C.OUT_SEQ / "registry.csv"))}
    ref_seq = reg[target_id]["sequence"]
    if not src.exists():
        return {"target": target_id, "status": "missing_crystal", "file": meta["crystal"]}

    chain, identity = _best_chain_for(src, ref_seq)
    keep = _aligned_span_resids(chain, ref_seq)

    structure = _parser(src).get_structure(target_id, str(src))
    model = list(structure)[0]
    model[chain.cid].id = "A" if chain.cid != "A" else "A"
    dst = out_dir / f"{target_id}.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(dst), _KeepSpan("A" if chain.cid == "A" else chain.cid, keep))

    kept = sio.longest_chain(sio.get_chains(dst))
    return {"target": target_id, "status": "ok", "file": meta["crystal"],
            "src_chain": chain.cid, "identity": round(identity, 3),
            "kept_residues": len(kept.seq), "crystal_chain_len": len(chain.seq),
            "out": str(dst.relative_to(C.REPO_ROOT))}


def extract_timp3(out_dir: Path) -> dict:
    src = C.TIMP3_CRYSTAL
    if not src.exists():
        return {"target": "TIMP3", "status": "missing_crystal"}
    structure = PDBParser(QUIET=True).get_structure("TIMP3", str(src))
    model = list(structure)[0]
    best = sio.longest_chain(sio.get_chains(src))
    dst = out_dir / "TIMP3.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(dst), _KeepSpan(best.cid, best.resids))
    return {"target": "TIMP3", "status": "ok", "src_chain": best.cid,
            "kept_residues": len(best.seq), "out": str(dst.relative_to(C.REPO_ROOT))}


def main() -> None:
    out_dir = C.OUT_CLEAN_TARGETS
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = [extract_target(t, out_dir) for t in C.TARGET_ORDER]
    rows.append(extract_timp3(out_dir))

    keys = sorted({k for r in rows for k in r})
    with open(out_dir / "extraction_manifest.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)

    print(f"Clean targets -> {out_dir.relative_to(C.REPO_ROOT)}")
    for r in rows:
        print(f"  {r['target']:<8} {r['status']:<16} "
              f"kept {r.get('kept_residues','?')} aa "
              f"(from {r.get('crystal_chain_len','?')}), id={r.get('identity','')}")


if __name__ == "__main__":
    main()
