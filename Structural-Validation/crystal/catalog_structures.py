"""
Catalog the experimental structures already in the repo and map them to the
validation panel, so model-vs-crystal comparison has a clean reference set.

For each mapped file it records: path, format, chains, per-chain length, any
HETATM ligands (e.g. catalytic ZN), and resolution if present in the header.
Copies the mapped files into crystal/structures/ as a tidy working set and
writes crystal/catalog.csv.

Run:  python Structural-Validation/crystal/catalog_structures.py
"""
from __future__ import annotations

import csv
import shutil
import sys
import warnings
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import is_aa

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C  # noqa: E402

warnings.filterwarnings("ignore")


def _parser(path: Path):
    return MMCIFParser(QUIET=True) if path.suffix.lower() == ".cif" \
        else PDBParser(QUIET=True)


def inspect(path: Path) -> dict:
    """Return chain composition, ligands, and resolution for a structure file."""
    s = _parser(path).get_structure(path.stem, str(path))
    model = next(iter(s))
    chains, ligands = [], set()
    for ch in model:
        n_aa = sum(1 for r in ch if is_aa(r, standard=True))
        het = {r.resname.strip() for r in ch
               if r.id[0] != " " and r.resname.strip() not in ("HOH", "WAT")}
        ligands |= het
        if n_aa:
            chains.append(f"{ch.id}:{n_aa}aa")
    res = getattr(s.header, "get", lambda *_: None) and s.header.get("resolution")
    return {
        "chains": ";".join(chains) if chains else "none",
        "n_chains": len(chains),
        "ligands": ";".join(sorted(ligands)) if ligands else "",
        "resolution": res if res else "",
    }


def build_map() -> list[tuple[str, str, Path]]:
    """(entity, role, path) tuples for every reference structure we can use."""
    items: list[tuple[str, str, Path]] = []

    # Target monomer crystals
    for tgt in C.TARGET_ORDER:
        p = C.CRYSTAL_DIR / C.TARGETS[tgt]["crystal"]
        if p.exists():
            items.append((tgt, "target_crystal", p))
        af = C.CRYSTAL_DIR / C.TARGETS[tgt]["af_cif"]
        if af.exists():
            items.append((tgt, "target_af_cif", af))

    # TIMP3 monomer crystal
    if C.TIMP3_CRYSTAL.exists():
        items.append(("TIMP3", "timp3_crystal", C.TIMP3_CRYSTAL))

    # Reference complexes (native + approximate homologous)
    for tgt, bytype in C.COMPLEX_REFERENCES.items():
        for ref_type, paths in bytype.items():
            for p in paths:
                if p.exists():
                    items.append((tgt, f"complex_ref_{ref_type}", p))
    return items


def main() -> None:
    C.OUT_CRYSTAL.mkdir(parents=True, exist_ok=True)
    store = C.OUT_CRYSTAL / "structures"
    store.mkdir(exist_ok=True)

    rows = []
    for entity, role, path in build_map():
        try:
            info = inspect(path)
        except Exception as e:
            info = {"chains": f"PARSE_ERROR: {e}", "n_chains": 0,
                    "ligands": "", "resolution": ""}
        dest = store / path.name
        try:
            shutil.copy(path, dest)
        except Exception:
            dest = path
        rows.append({
            "entity": entity, "role": role, "file": path.name,
            "format": path.suffix.lstrip("."), **info,
            "stored_path": str(dest.relative_to(C.REPO_ROOT))
            if str(dest).startswith(str(C.REPO_ROOT)) else str(dest),
        })

    fields = ["entity", "role", "file", "format", "n_chains", "chains",
              "ligands", "resolution", "stored_path"]
    with open(C.OUT_CRYSTAL / "catalog.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    print(f"Cataloged {len(rows)} reference structures -> "
          f"{(C.OUT_CRYSTAL / 'catalog.csv').relative_to(C.REPO_ROOT)}\n")
    print(f"{'entity':<8} {'role':<18} {'chains':<24} ligands")
    for r in rows:
        print(f"{r['entity']:<8} {r['role']:<18} {r['chains']:<24} {r['ligands']}")

    have = {r["entity"] for r in rows if r["role"] == "target_crystal"}
    missing = [t for t in C.TARGET_ORDER if t not in have]
    if missing:
        print(f"\nTargets without a monomer crystal: {missing}")
    cref = {r['entity']: r['role'].replace('complex_ref_', '')
            for r in rows if r['role'].startswith('complex_ref')}
    print(f"Complex references (target:type): {cref or 'none'}")


if __name__ == "__main__":
    main()
