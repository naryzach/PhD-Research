"""
Aggregate the per-model / per-complex metric tables into tidy, analysis-ready
products and a human-readable data dictionary.

Reads (analysis/):
  master_monomer_metrics.csv, master_complex_metrics.csv,
  af3_vs_esmfold_monomer.csv
Writes (analysis/):
  monomer_wide_<metric>.csv        entity x folder pivot per key metric
  complex_matrix_<method>_<metric>.csv   construct x target matrix per metric
  metric_correlations_<method>.csv Spearman between complex metrics (when data)
  data_dictionary.csv              every column, its meaning, units, direction
  summary.json                     completion status + quick counts

Runs today; correlation/pivots populate as real outputs arrive.

Run:  python Structural-Validation/analysis/aggregate.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C  # noqa: E402

A = C.OUT_ANALYSIS

# metric -> (definition, unit, better_direction)
DATA_DICTIONARY = {
    # monomer quality
    "rmsd_ca": ("CA RMSD of model vs crystal after Kabsch superposition on "
                "sequence-matched residues", "Å", "lower"),
    "tm_score": ("TM-score vs crystal, length-normalised (fold-level similarity)",
                 "0-1", "higher"),
    "gdt_ts": ("Global Distance Test - Total Score (mean %CA within 1/2/4/8 Å)",
               "0-100", "higher"),
    "gdt_ha": ("GDT High Accuracy (0.5/1/2/4 Å)", "0-100", "higher"),
    "lddt": ("Local Distance Difference Test (superposition-free)", "0-100", "higher"),
    "rg_model": ("Radius of gyration of the model", "Å", "match-crystal"),
    "clashscore": ("Steric clashes per 1000 atoms (heavy-atom, |Δres|≥2, <2 Å)",
                   "per 1k", "lower"),
    "rama_outlier_frac": ("Fraction of residues outside broad Ramachandran basins "
                          "(approximate)", "0-1", "lower"),
    "plddt_mean": ("Mean per-residue confidence (AF3 pLDDT / ESMFold2 pLDDT)",
                   "0-100", "higher"),
    "ptm": ("Predicted TM-score (whole-model confidence)", "0-1", "higher"),
    "rmsd_af3_esm": ("CA RMSD between AF3 and ESMFold2 models of the same sequence",
                     "Å", "lower=agree"),
    "tm_af3_esm": ("TM-score between AF3 and ESMFold2 models", "0-1", "higher=agree"),
    # complex / interface
    "bsa": ("Buried surface area at the interface (ΔSASA/2)", "Å²", "higher"),
    "n_iface_res_construct": ("TIMP3-side interface residues (<5 Å)", "count", "context"),
    "n_iface_res_target": ("Target-side interface residues (<5 Å)", "count", "context"),
    "n_contacts_5A": ("Heavy-atom contact pairs within 5 Å across interface",
                      "count", "higher"),
    "n_hbonds": ("Interface hydrogen bonds (N/O donor-acceptor <3.5 Å)", "count", "higher"),
    "n_salt_bridges": ("Interface salt bridges (charged pairs <4 Å)", "count", "higher"),
    "n_hydrophobic": ("Interface hydrophobic C-C contacts <4.5 Å", "count", "higher"),
    "contact_density": ("Contacts per interface residue (packing/shape proxy)",
                        "ratio", "higher"),
    "min_ca_ca_zincloop": ("Min CA-CA distance, TIMP3 reactive edge (C1-P5) to target "
                           "HExxHxxGxxH zinc motif", "Å", "lower"),
    "haddock_score": ("HADDOCK score (weighted vdW+elec+desolv+AIR)", "a.u.", "lower"),
    "vdw": ("van der Waals energy", "kcal/mol", "lower"),
    "elec": ("electrostatic energy", "kcal/mol", "lower"),
    "desolv": ("desolvation energy", "kcal/mol", "lower"),
    "air": ("ambiguous restraints energy", "kcal/mol", "lower"),
    "haddock_bsa": ("HADDOCK-reported buried surface area", "Å²", "higher"),
    "violations": ("AIR restraint violations", "count", "lower"),
    "iptm": ("AF3 interface predicted TM-score (co-fold binding confidence)",
             "0-1", "higher"),
    "pae": ("AF3 mean predicted aligned error", "Å", "lower"),
    "dockq": ("DockQ quality vs native complex (fnat/iRMS/LRMS composite)",
              "0-1", "higher"),
    "fnat": ("fraction of native residue contacts recovered", "0-1", "higher"),
    "irms": ("interface RMSD vs native", "Å", "lower"),
    "lrms": ("ligand RMSD vs native (receptor-fit)", "Å", "lower"),
    "capri": ("CAPRI class (high/medium/acceptable/incorrect)", "class", "higher"),
}

MONOMER_KEY = ["rmsd_ca", "tm_score", "gdt_ts", "lddt", "plddt_mean"]
COMPLEX_KEY = ["bsa", "n_contacts_5A", "n_hbonds", "n_salt_bridges",
               "contact_density", "min_ca_ca_zincloop", "haddock_score",
               "iptm", "dockq"]


def _read(name):
    p = A / name
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


def _slug(s: str) -> str:
    return s.replace(":", "_").replace("x", "X").replace(" ", "")


def main() -> None:
    A.mkdir(parents=True, exist_ok=True)
    # Clear derived products so the output always reflects the current inputs
    # (master_*.csv are written by the upstream drivers and are kept).
    for pat in ("monomer_wide_*.csv", "complex_matrix_*.csv",
                "metric_correlations_*.csv"):
        for old in A.glob(pat):
            old.unlink()
    mono = _read("master_monomer_metrics.csv")
    cplx = _read("master_complex_metrics.csv")
    agree = _read("af3_vs_esmfold_monomer.csv")

    # data dictionary
    pd.DataFrame(
        [{"column": k, "definition": v[0], "unit": v[1], "better": v[2]}
         for k, v in DATA_DICTIONARY.items()]
    ).to_csv(A / "data_dictionary.csv", index=False)

    # monomer pivots (entity x folder)
    if not mono.empty:
        for m in MONOMER_KEY:
            if m in mono.columns:
                piv = mono.pivot_table(index="entity_id", columns="folder",
                                       values=m, aggfunc="first")
                piv.to_csv(A / f"monomer_wide_{m}.csv")

    # complex matrices (construct x target) per source + metric
    if not cplx.empty and "source" in cplx:
        for source in sorted(cplx["source"].dropna().unique()):
            tag = _slug(source)
            sub = cplx[cplx["source"] == source]
            for m in COMPLEX_KEY:
                if m in sub.columns and sub[m].notna().any():
                    mat = sub.pivot_table(index="construct_id", columns="target_id",
                                          values=m, aggfunc="first")
                    if mat.notna().any().any():
                        mat.to_csv(A / f"complex_matrix_{tag}_{m}.csv")
            numeric = sub[[c for c in COMPLEX_KEY if c in sub.columns]].apply(
                pd.to_numeric, errors="coerce")
            usable = numeric.dropna(axis=1, how="all")
            if usable.shape[1] >= 2 and usable.dropna().shape[0] >= 4:
                usable.corr(method="spearman").to_csv(
                    A / f"metric_correlations_{tag}.csv")

    # summary / completion
    def _status_counts(df):
        return df["status"].value_counts().to_dict() if "status" in df else {}

    reg_path = C.OUT_SEQ / "registry.csv"
    n_seq = int(len(pd.read_csv(reg_path))) if reg_path.exists() else 0
    summary = {
        "sequences": n_seq,
        "monomer_rows": int(len(mono)),
        "monomer_status": _status_counts(mono),
        "complex_rows": int(len(cplx)),
        "complex_status": _status_counts(cplx),
        "af3_vs_esmfold_status": _status_counts(agree),
        "outputs_ready": {
            "af3_monomer_models": _count_glob(C.OUT_AF3 / "results", "*.cif"),
            "esmfold_models": _count_glob(C.OUT_ESM / "results", "*.pdb"),
            "haddock_best_models": _count_glob(C.OUT_DOCK, "*_HADDOCK.pdb"),
        },
    }
    (A / "summary.json").write_text(json.dumps(summary, indent=2))

    print("Aggregation complete.")
    print(json.dumps(summary, indent=2))
    print(f"\nData dictionary: {(A / 'data_dictionary.csv').relative_to(C.REPO_ROOT)} "
          f"({len(DATA_DICTIONARY)} metrics documented)")


def _count_glob(root: Path, pat: str) -> int:
    return len(list(root.rglob(pat))) if root.exists() else 0


if __name__ == "__main__":
    main()
