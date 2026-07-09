"""
Download homologous TIMP:metalloprotease co-crystals used as APPROXIMATE DockQ
references for targets that lack a native TIMP3 complex.

  1UEA  TIMP-1 : MMP-3 catalytic domain  (used for MMP2/MMP9/MMP3)

ADAM10 reuses the in-repo TIMP3:ADAM17 complex (sister ADAM, same TIMP3 ligand),
so nothing is downloaded for it. DockQ against any of these is approximate and is
labelled `dockq_ref_type=approximate` in master_complex_metrics.csv.

Run:  python Structural-Validation/crystal/fetch_homologous_refs.py
"""
from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C  # noqa: E402

PDB_IDS = {
    "1UEA": "TIMP-1 : MMP-3 catalytic domain (approx ref for MMP2/MMP9/MMP3)",
}


def main() -> None:
    C.HOMOLOG_DIR.mkdir(parents=True, exist_ok=True)
    for pdb_id, desc in PDB_IDS.items():
        dst = C.HOMOLOG_DIR / f"{pdb_id}.pdb"
        if dst.exists() and dst.stat().st_size > 1000:
            print(f"  have {pdb_id}.pdb  ({desc})")
            continue
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        try:
            urllib.request.urlretrieve(url, dst)
            print(f"  fetched {pdb_id}.pdb  ({dst.stat().st_size} bytes)  {desc}")
        except Exception as e:
            print(f"  FAILED {pdb_id}: {e}  -- MMP DockQ will stay pending")
    print(f"\nApproximate references -> {C.HOMOLOG_DIR.relative_to(C.REPO_ROOT)}")


if __name__ == "__main__":
    main()
