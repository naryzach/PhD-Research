"""
Resolve where each model / complex structure lives, tolerant of the many ways
AF3-server downloads and HADDOCK runs get named. Returns None when an output
has not been produced yet, so the analysis tables stay complete (with a
`status` column) instead of crashing.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C  # noqa: E402


def load_registry() -> list[dict]:
    with open(C.OUT_SEQ / "registry.csv", newline="") as fh:
        return list(csv.DictReader(fh))


def _sanitize(name: str) -> str:
    """Reproduce the AlphaFold-server job-name sanitization: lowercase and
    collapse runs of underscores (so 'AB_5__MMP2' -> 'ab_5_mmp2', 'AB_1'->'ab_1')."""
    return re.sub(r"_+", "_", name.lower())


def _first(patterns, root: Path) -> Path | None:
    for pat in patterns:
        hits = [p for p in sorted(root.rglob(pat)) if ":" not in p.name]  # skip ADS
        if hits:
            hits.sort(key=lambda p: ("model_0" not in p.name.lower(),
                                     "rank_1" not in p.name.lower(), len(str(p))))
            return hits[0]
    return None


def _af3_job_dir(job_name: str) -> Path | None:
    """The exact AF3 result directory for a job, matched on the sanitized name so
    a target monomer ('mmp2') never collides with a co-fold folder ('ab_1_mmp2')."""
    root = C.OUT_AF3 / "results"
    if not root.exists():
        return None
    san = _sanitize(job_name)
    d = root / san
    if d.is_dir():
        return d
    # tolerate a wrapping folder (e.g. an unzipped batch dir)
    for cand in root.rglob(san):
        if cand.is_dir():
            return cand
    return None


def af3_monomer(entity_id: str) -> Path | None:
    d = _af3_job_dir(entity_id)
    if d is None:
        return None
    return _first(["*model_0.cif", "*model*.cif", "*.cif", "*.pdb"], d)


def af3_complex(construct_id: str, target_id: str) -> Path | None:
    d = _af3_job_dir(f"{construct_id}__{target_id}")
    if d is None:
        return None
    return _first(["*model_0.cif", "*model*.cif", "*.cif"], d)


def esmfold_monomer(entity_id: str) -> Path | None:
    for ext in ("pdb", "cif"):
        p = C.OUT_ESM / "results" / f"{entity_id}.{ext}"
        if p.exists():
            return p
    return None


def esmfold_complex(construct_id: str, target_id: str) -> Path | None:
    root = C.OUT_ESM / "results_complex"
    if not root.exists():
        return None
    key = f"{construct_id}__{target_id}"
    for ext in ("pdb", "cif"):
        p = root / f"{key}.{ext}"
        if p.exists():
            return p
    return _first([f"*{key}*.pdb", f"*{key}*.cif"], root)


def crystal_target(target_id: str) -> Path | None:
    p = C.OUT_CLEAN_TARGETS / f"{target_id}.pdb"
    return p if p.exists() else None


def dock_track(construct_src: str, target_src: str) -> str:
    return f"{construct_src}__{target_src}"


def haddock_complex(construct_id: str, target_id: str,
                    construct_src: str = "AF3",
                    target_src: str = "AF3") -> Path | None:
    """Best HADDOCK model for a construct-target pair on a given docking track."""
    root = C.OUT_DOCK / dock_track(construct_src, target_src) / "best_models"
    cand = root / f"{construct_id}__{target_id}_HADDOCK.pdb"
    if cand.exists():
        return cand
    if root.exists():
        return _first([f"*{construct_id}*{target_id}*.pdb"], root)
    return None


def monomer_model(entity_id: str, folder: str) -> Path | None:
    return af3_monomer(entity_id) if folder == "AF3" else esmfold_monomer(entity_id)


def crystal_reference(entity_id: str) -> Path | None:
    if entity_id.startswith("TIMP3") or entity_id.startswith(("AB_", "C_")):
        return C.TIMP3_CRYSTAL if C.TIMP3_CRYSTAL.exists() else None
    meta = C.TARGETS.get(entity_id)
    if meta:
        p = C.CRYSTAL_DIR / meta["crystal"]
        return p if p.exists() else None
    return None


def complex_reference(target_id: str) -> tuple[Path, str] | tuple[None, None]:
    """(path, ref_type) for the best DockQ reference; native preferred over
    approximate (homologous). Returns (None, None) if none is available."""
    spec = C.COMPLEX_REFERENCES.get(target_id, {})
    for ref_type in ("native", "approximate"):
        for p in spec.get(ref_type, []):
            if p.exists():
                return p, ref_type
    return None, None
