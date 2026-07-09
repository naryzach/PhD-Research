"""
Assemble every sequence in the validation panel into one clean registry.

Inputs (all already in the repo):
  * TIMP3 constructs  -> Local/TIMP3_Redesign_2026-07/data/esmc_predict_input_fcs_constructs.csv
  * WT TIMP3          -> config.TIMP3_WT_MATURE (UniProt P35625, mature)
  * Target CDs        -> Data/Targets_w_CDs.xlsx (Catalytic_Domain column)

Outputs (Local/TIMP3_Structural_Validation_2026-07/sequences/):
  * registry.csv          one row per entity (constructs + targets) with metadata
  * all_sequences.fasta   every entity, one record each
  * fasta/<id>.fasta      per-entity FASTA (used by AF3/ESMFold input builders)

Run:  python Structural-Validation/sequences/build_sequences.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import openpyxl

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C  # noqa: E402


def _to_ndomain(full_seq: str) -> str:
    """Return the TIMP3 N-terminal binding domain from a full construct sequence.

    Restores the native CT dipeptide and truncates at the N-domain terminus so
    every construct is numbered identically to mature TIMP3 and to prior runs.
    """
    seq = full_seq.strip().upper()
    motif = C.TIMP3_NDOMAIN_END_MOTIF
    if motif not in seq:
        raise ValueError(f"N-domain end motif {motif!r} not found in sequence")
    core = seq[: seq.index(motif) + len(motif)]
    # Prepend CT only if not already present (WT already starts with CT).
    if not core.startswith(C.TIMP3_NDOMAIN_PREFIX):
        core = C.TIMP3_NDOMAIN_PREFIX + core
    return core


def load_constructs() -> list[dict]:
    """WT + the 13 loop-graft constructs, reduced to the N-terminal domain."""
    rows: list[dict] = []

    wt_nd = _to_ndomain(C.TIMP3_WT_MATURE)
    rows.append({
        "id": "TIMP3_WT", "kind": "construct", "name": "TIMP3 WT",
        "loop": "EGPFGT/ASESLC (native)", "loop_position": "native",
        "sequence": wt_nd, "length": len(wt_nd),
        "source": "UniProt P35625 (mature, N-domain)", "notes": "wild-type reference",
    })

    with open(C.CONSTRUCTS_CSV, newline="") as fh:
        for r in csv.DictReader(fh):
            name = (r.get("Construct") or "").strip()
            full = (r.get("Full Seq") or "").strip()
            if not name or not full:
                continue
            nd = _to_ndomain(full)
            rows.append({
                "id": name.replace(" ", "_"), "kind": "construct", "name": name,
                "loop": (r.get("Residues") or "").strip(),
                "loop_position": (r.get("loop_position") or "").strip(),
                "sequence": nd, "length": len(nd),
                "source": "FCS ordered construct (N-domain)", "notes": "",
            })
    return rows


def load_targets() -> list[dict]:
    """Core-5 target catalytic domains from the spreadsheet, cleaned."""
    wb = openpyxl.load_workbook(C.TARGETS_XLSX, data_only=True)
    ws = wb.active
    all_rows = list(ws.iter_rows(values_only=True))
    header = list(all_rows[0])
    by_name = {}
    for raw in all_rows[1:]:
        d = dict(zip(header, raw))
        if d.get("Target"):
            by_name[str(d["Target"]).strip()] = d

    rows: list[dict] = []
    for tgt in C.TARGET_ORDER:
        meta = C.TARGETS[tgt]
        d = by_name.get(tgt)
        if d is None:
            print(f"  WARNING: {tgt} not found in {C.TARGETS_XLSX.name}")
            continue
        cd = str(d.get("Catalytic_Domain") or "").strip().rstrip("*").strip()
        active = str(d.get("Active_site") or "").strip().rstrip("*").strip()

        notes = []
        # Flag common purification tags so downstream steps can trim if wanted.
        for tag in ("DYKDDDDK", "HHHHHH"):
            if tag in cd:
                notes.append(f"contains {('FLAG' if tag[0]=='D' else 'His')} tag")

        motif_pos = cd.find(active) + 1 if active and active in cd else -1
        rows.append({
            "id": tgt, "kind": "target", "name": tgt,
            "loop": active, "loop_position": (
                f"{motif_pos}-{motif_pos + len(active) - 1}" if motif_pos > 0 else "NA"
            ),
            "sequence": cd, "length": len(cd),
            "source": f"{meta['vendor']} ({meta['source']}); UniProt {meta['uniprot']}",
            "notes": "; ".join(notes),
        })
    return rows


def write_outputs(rows: list[dict]) -> None:
    C.ensure_out_dirs()
    fasta_dir = C.OUT_SEQ / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)

    fields = ["id", "kind", "name", "loop", "loop_position",
              "length", "source", "notes", "sequence"]
    with open(C.OUT_SEQ / "registry.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    with open(C.OUT_SEQ / "all_sequences.fasta", "w") as fh:
        for r in rows:
            fh.write(f">{r['id']} kind={r['kind']} len={r['length']}\n{r['sequence']}\n")

    for r in rows:
        with open(fasta_dir / f"{r['id']}.fasta", "w") as fh:
            fh.write(f">{r['id']}\n{r['sequence']}\n")


def main() -> None:
    constructs = load_constructs()
    targets = load_targets()
    rows = constructs + targets
    write_outputs(rows)

    print(f"Constructs: {len(constructs)}   Targets: {len(targets)}   "
          f"Total: {len(rows)}")
    print(f"\nRegistry -> {C.OUT_SEQ / 'registry.csv'}")
    print(f"{'id':<12} {'kind':<10} {'len':>4}  loop")
    for r in rows:
        print(f"{r['id']:<12} {r['kind']:<10} {r['length']:>4}  {r['loop']}")


if __name__ == "__main__":
    main()
