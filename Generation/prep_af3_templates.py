"""
Build AF3-co-fold RFd3 templates to REPLACE the HADDOCK templates.

Why: iterative_refinement.py fed RFd3 the HADDOCK complexes (Data/TIMP_Complexes/
HADDOCK_Outputs/<TARGET>_TIMP3_HADDOCK.pdb) as the design template. RFd3 holds the
target fixed and inpaints the binder loops INTO the template's relative pose — so
the loops were being designed against HADDOCK's binding orientation, which the
structural validation showed is ~30 A off-native (LRMS ~30 A, DockQ ~0.05). The
AF3/ESMFold2 co-folds recover the native mode (LRMS ~2.3 A, DockQ ~0.76), so they
are the correct template.

Each output template: chain A = TIMP3 N-terminal scaffold (resid 1-121, matching
TARGETS[*]["scaffold_len"]), chain B = target, in the AF3 co-fold pose.

Source preference per target:
  1. Local/.../af3/results/timp3_wt_<target>/  — the 121-aa scaffold co-folded with
     the target (cleanest: exactly the scaffold, no extra domain).
  2. Data/TIMP_Complexes/AlphaFold_PDB/PDB_fold_timp3_*_<target>cd_wt_model_0.pdb —
     the 188-aa TIMP3 co-fold, truncated to residues 1-121 (fallback, e.g. MMP10).

Output: Data/TIMP_Complexes/AF3_Templates/<TARGET>_TIMP3_AF3.pdb

Run:  python Generation/prep_af3_templates.py
"""
from __future__ import annotations

import glob
import sys
from pathlib import Path

from Bio.PDB import PDBIO, Select

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parent
sys.path.insert(0, str(_REPO / "Structural-Validation" / "utils"))
sys.path.insert(0, str(_REPO / "Structural-Validation" / "analysis"))
import structure_io as sio          # noqa: E402
from complex_metrics import _seq_identity  # noqa: E402

# N-domain used only to identify the TIMP3 (binder) chain by its first 121 residues.
TIMP3_NDOMAIN = ("CTCSPSHPQDAFCNSDIVIRAKVVGKKLVKEGPFGTLVYTIKQMKMYRGFTKMPHVQYIHTE"
                 "ASESLCGLKLEVNKYQYLLTGRVYDGKMYTGLCNFVERWDQLTLSQRKGLNYRYHLGCN")  # 121 aa

# FULL-LENGTH TIMP3 (188 aa). The pipeline designs the N-domain loops but ORDERS +
# tests the full protein (N-domain design + native C-domain), and the literature
# says the C-domain is needed for ADAM10 binding — so RFd3 must design in the
# full-length context. Templates keep all 188 TIMP3 residues; the C-domain (122-188)
# becomes fixed scaffold in the RFd3 contig, and scaffold_len=188 in TARGETS.
SCAFFOLD_LEN = 188

TC = _REPO / "Data" / "TIMP_Complexes"
FULL_FOLDS = TC / "AF3_Full_Folds"          # fresh full-length AF3 co-folds (5-model)
OUT_DIR = TC / "AF3_Templates"

TARGETS = ["ADAM10", "ADAM17", "MMP2", "MMP3", "MMP9", "MMP10"]
# older 188-aa TIMP3 co-fold PDBs (names inconsistent: _v_ vs _variant_) — fallback.
COFOLD188_GLOB = "PDB_fold_timp3_*_{t}cd_wt_model_0.pdb"


def _source_for(target: str):
    """Prefer the fresh full-length AF3 co-folds in Data/AF3_Full_Folds (5-model,
    confidence-scored); fall back to the older single-model AlphaFold_PDB co-folds."""
    fresh = FULL_FOLDS / f"timp3full_{target.lower()}"
    if fresh.is_dir():
        cif = [c for c in sorted(fresh.glob("*model_0.cif")) if ":" not in c.name]
        if cif:
            return cif[0], "fresh_af3_188"
    hits = glob.glob(str(TC / "AlphaFold_PDB" / COFOLD188_GLOB.format(t=target.lower())))
    if hits:
        return Path(hits[0]), "alphafold_pdb_188"
    return None, None


class _TemplateSelect(Select):
    """Keep binder chain residues 1..SCAFFOLD_LEN + all target-chain residues."""
    def __init__(self, binder_cid, target_cid):
        self.b, self.t = binder_cid, target_cid

    def accept_chain(self, chain):
        return chain.id in (self.b, self.t)

    def accept_residue(self, residue):
        if residue.id[0] != " ":
            return False
        if residue.get_parent().id == self.b:
            return residue.id[1] <= SCAFFOLD_LEN
        return True


def _force_chain(path: Path, mapping: dict) -> None:
    """Relabel chain IDs (col 22) per {old: new}."""
    out = []
    for ln in path.read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")) and len(ln) > 21 and ln[21] in mapping:
            ln = ln[:21] + mapping[ln[21]] + ln[22:]
        out.append(ln)
    path.write_text("\n".join(out) + "\n")


def build_one(target: str) -> dict:
    src, kind = _source_for(target)
    if src is None:
        return {"target": target, "status": "no_af3_source"}
    chains = sio.get_chains(src)
    binder = max(chains.values(), key=lambda c: _seq_identity(c.seq[:121], TIMP3_NDOMAIN))
    tgt = max((c for c in chains.values() if c.cid != binder.cid), key=lambda c: len(c.seq))
    bid_check = _seq_identity(binder.seq[:121], TIMP3_NDOMAIN)   # first 121 = N-domain
    if bid_check < 0.9:
        return {"target": target, "status": f"binder_id_low({bid_check:.2f})"}

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_pdb = OUT_DIR / f"{target}_TIMP3_AF3.pdb"
    st = sio.load(src)
    io = PDBIO()
    io.set_structure(st)
    io.save(str(out_pdb), _TemplateSelect(binder.cid, tgt.cid))
    # normalise to the RFd3/TARGETS convention: binder = chain A, target = chain B
    remap = {}
    if binder.cid != "A":
        remap[binder.cid] = "A"
    if tgt.cid != "B":
        remap[tgt.cid] = "B"
    if remap:
        # two-step relabel via temp letters to avoid collisions
        tmp = {binder.cid: "\x01", tgt.cid: "\x02"}
        _force_chain(out_pdb, tmp)
        _force_chain(out_pdb, {"\x01": "A", "\x02": "B"})

    # verify what we wrote
    v = sio.get_chains(out_pdb)
    a_len = len(v["A"].seq) if "A" in v else 0
    b_len = len(v["B"].seq) if "B" in v else 0
    ok = (a_len == SCAFFOLD_LEN and b_len == len(tgt.seq)
          and _seq_identity(v["A"].seq[:121], TIMP3_NDOMAIN) > 0.9)
    return {"target": target, "status": "ok" if ok else "verify_failed",
            "source": src.name, "kind": kind, "binderA_len": a_len,
            "targetB_len": b_len, "out": out_pdb.name}


def main() -> None:
    print(f"Building AF3 templates -> {OUT_DIR.relative_to(_REPO)}\n")
    hdr = f"{'target':<9}{'status':<18}{'A(TIMP3)':>9}{'B(target)':>10}  {'kind':<28}source"
    print(hdr); print("-" * len(hdr))
    rows = [build_one(t) for t in TARGETS]
    for r in rows:
        print(f"{r['target']:<9}{r['status']:<18}{r.get('binderA_len',''):>9}"
              f"{r.get('targetB_len',''):>10}  {r.get('kind',''):<28}{r.get('source','')}")
    ok = sum(1 for r in rows if r["status"] == "ok")
    print(f"\n{ok}/{len(rows)} templates built.")
    if ok == len(rows):
        print("\nNext: point iterative_refinement.py at these — set DATA_DIR to "
              "Data/TIMP_Complexes/AF3_Templates and, in TARGETS, "
              'pdb="<TARGET>_TIMP3_AF3.pdb", binder_chain="A", target_chain="B".')


if __name__ == "__main__":
    main()
