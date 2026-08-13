#!/usr/bin/env python3
"""
recover_from_cifs.py — rebuild ESMFold2 scores from predicted-complex CIFs that
are ALREADY ON DISK. No GPU, no re-folding.

Why this works: the ESMFold2 stage score deliberately excludes esm_iptm (its
negative in-sample correlation with binding is a selection-bias artifact — we
never observe the non-folders). calc_composite's ESMFold2 branch reads exactly
three inputs:

    esm_plddt                    <- mean CA B-factor of the binder chain
    esm_iface_n_iface_res        <- interface geometry
    esm_iface_contact_density    <- interface geometry

...and all three are recoverable from the saved complex. So when a run wrote its
structures but lost the scores (Aug-3 2026: 17,194 CIFs, zero esm_scores.csv),
the ranking is recoverable from disk in minutes instead of days of GPU.

What is NOT recovered: esm_iptm / esm_ptm are model-head outputs, absent from the
structure. They do not enter the gate or the composite; the only consumer is the
`promising` flag, left False for recovered rows. esm_lplddt / esm_pae are
log-only and also left blank (the raw PAE matrix survives in <design>_pae.npy
if they are ever wanted).

    python Generation/recover_from_cifs.py --dry-run
    python Generation/recover_from_cifs.py
    python Generation/recover_from_cifs.py --out-base Local/specificity_refinement

Then run rescore_pool.py to fold only the designs that have no CIF at all.
Idempotent and safe to re-run; each round summary is backed up once to *.bak.
"""
import argparse
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
import calibrated_scoring as cs        # numpy-only; no atomworks/torch needed

BINDER_CHAIN, TARGET_CHAIN = "A", "B"


def binder_plddt_from_cif(cif_text: str, binder_chain: str = BINDER_CHAIN) -> float:
    """Mean per-residue pLDDT of the binder chain (CA atoms only, so it matches the
    scorer's per-residue mean rather than an atom-count-weighted one)."""
    cols, in_loop, vals = [], False, []
    for ln in cif_text.splitlines():
        s = ln.strip()
        if s.startswith("_atom_site."):
            cols.append(s.split(".")[1]); in_loop = True
        elif in_loop and s.startswith(("ATOM", "HETATM")):
            p = s.split()
            if len(p) < len(cols):
                continue
            r = dict(zip(cols, p))
            if r.get("auth_asym_id", r.get("label_asym_id")) != binder_chain:
                continue
            if r.get("label_atom_id", r.get("auth_atom_id")) != "CA":
                continue
            try:
                vals.append(float(r["B_iso_or_equiv"]))
            except (ValueError, KeyError):
                continue
        elif in_loop and s.startswith("#"):
            in_loop = False
    return float(np.mean(vals)) if vals else float("nan")


def find_cif(out_base: Path, it: str, target: str, design_id: str, side: str = None):
    """
    Locate <out_base>/<it>/esmfold2/<target>/structures/<design_id>[<side>].{cif,pdb}.

    `side` is "on"/"off" for the specificity pool, where each design is folded against
    both targets and the scorer suffixes the id with "::on" / "::off". Those land on
    disk URL-encoded (`%3A%3Aon`), so try both spellings.
    """
    d = out_base / it / "esmfold2" / target / "structures"
    if side:
        stems = [f"{design_id}%3A%3A{side}", f"{design_id}::{side}"]
    else:
        stems = [design_id]
    for stem in stems:
        for ext in (".cif", ".pdb"):
            p = d / f"{stem}{ext}"
            if p.exists():
                return p
    return None


def recover_one(cif_path: Path) -> dict:
    """{esm_plddt, esm_iface_*} from one predicted complex, or {} if unreadable."""
    try:
        txt = cif_path.read_text(errors="ignore")
        f = cs.interface_features_from_cif(txt, binder_chain=BINDER_CHAIN,
                                           target_chain=TARGET_CHAIN)
        plddt = binder_plddt_from_cif(txt)
    except Exception:
        return {}
    if not np.isfinite(plddt):
        return {}
    return {
        "esm_plddt":                 plddt,
        "esm_iface_n_iface_res":     f.get("n_iface_res"),
        "esm_iface_contact_density": f.get("contact_density"),
        "esm_iface_iface_plddt":     f.get("iface_plddt"),
        "esm_cif":                   str(cif_path),
    }


_SPEC_MOD = None


def _spec_composite(on_metrics: dict, off_metrics: dict) -> float:
    """Selectivity-aware composite, delegated to specificity_refinement so the formula
    cannot drift from the pipeline's. Imported lazily: it pulls in iterative_refinement
    (atomworks), so --pair-mode needs the pipeline env while plain recovery does not."""
    global _SPEC_MOD
    if _SPEC_MOD is None:
        try:
            import specificity_refinement as sr
            _SPEC_MOD = sr
        except Exception as exc:
            sys.exit(f"--pair-mode needs specificity_refinement, which failed to import: {exc}\n"
                     "  Run it from the pipeline env (foundry). Plain recovery (no --pair-mode)\n"
                     "  has no such dependency. Nothing was modified.")
    return _SPEC_MOD.calc_specificity_composite(on_metrics, off_metrics, model="esm")


def main():
    ap = argparse.ArgumentParser(description="Rebuild ESMFold2 scores from CIFs on disk.")
    ap.add_argument("--out-base", default=None,
                    help="Run directory holding it_*/ (default: $REFINE_OUT_BASE or "
                         "Local/iterative_refinement).")
    ap.add_argument("--dry-run", action="store_true",
                    help="Report how many rows are recoverable and exit.")
    ap.add_argument("--rebuild-hof", action="store_true",
                    help="Also rebuild hof_summary.csv + refinement_state.json from the scored "
                         "pool (needs the pipeline env).")
    ap.add_argument("--loops", nargs="+", default=["AB", "C", "EF"],
                    help="Loops the run redesigned — must match the run being rebuilt, since the "
                         "HOF's loop-length bookkeeping is per-loop (default: AB C EF).")
    ap.add_argument("--pair-mode", action="store_true",
                    help="Specificity pool: each design was folded against BOTH targets "
                         "(<id>::on / <id>::off). Recovers both sides and scores with the "
                         "selectivity-aware composite.")
    args = ap.parse_args()

    out_base = Path(args.out_base) if args.out_base else Path(
        os.environ.get("REFINE_OUT_BASE") or (_HERE / ".." / "Local" / "iterative_refinement"))
    out_base = out_base.resolve()

    csvs = sorted(out_base.glob("it_*/round_summary.csv"),
                  key=lambda p: int(p.parent.name.split("_")[1]))
    if not csvs:
        sys.exit(f"No it_*/round_summary.csv under {out_base}")
    print(f"Run directory: {out_base}\nRound summaries: {len(csvs)}")

    n_rows = n_have_cif = n_recovered = n_already = 0
    per_target = {}
    missing = {}
    seen_targets = set()   # every target in the pool, not just newly-recovered ones

    for csv in csvs:
        df = pd.read_csv(csv)
        it = csv.parent.name
        changed = False
        for i, r in df.iterrows():
            n_rows += 1
            if isinstance(r.get("target_name"), str):
                seen_targets.add(r["target_name"])
            # "Scored" means the ranking inputs are present — NOT esm_iptm, which
            # the composite never reads.
            if "esm_plddt" in df.columns and pd.notna(r.get("esm_plddt")):
                n_already += 1
                continue
            tgt, did = r.get("target_name"), r.get("design_id")
            cif = find_cif(out_base, it, str(tgt), str(did), "on" if args.pair_mode else None)
            if cif is None:
                missing[tgt] = missing.get(tgt, 0) + 1
                continue
            n_have_cif += 1
            if args.dry_run:
                per_target[tgt] = per_target.get(tgt, 0) + 1
                continue
            m = recover_one(cif)
            if not m:
                missing[tgt] = missing.get(tgt, 0) + 1
                continue
            for k, v in m.items():
                df.at[i, k] = v
            on_metrics = {
                "esm_plddt":                 m["esm_plddt"],
                "esm_iface_contact_density": m["esm_iface_contact_density"],
                "esm_iface_n_iface_res":     m["esm_iface_n_iface_res"],
            }
            score = cs.esmfold2_stage_score(on_metrics)
            if args.pair_mode:
                # Off-target side gives the contact-density gap the selectivity
                # composite ranks on. Missing off side -> neutral gap (0.5), which is
                # exactly how the pipeline treats an unknown off-target.
                off_cif = find_cif(out_base, it, str(tgt), str(did), "off")
                off = recover_one(off_cif) if off_cif else {}
                if off:
                    df.at[i, "esm_off_iface_contact_density"] = off["esm_iface_contact_density"]
                    df.at[i, "esm_off_plddt"] = off["esm_plddt"]
                    df.at[i, "selectivity_cd"] = (float(m["esm_iface_contact_density"])
                                                  - float(off["esm_iface_contact_density"]))
                score = _spec_composite(on_metrics, off)
            df.at[i, "source"] = "ESMFold2"
            df.at[i, "composite_score"] = score
            changed = True
            n_recovered += 1
            per_target[tgt] = per_target.get(tgt, 0) + 1
        if changed:
            bak = Path(str(csv) + ".bak")
            if not bak.exists():
                shutil.copy2(csv, bak)
            df.to_csv(csv, index=False)
            print(f"  {it}: recovered, summary rewritten", flush=True)

    print(f"\nRows: {n_rows} | already scored: {n_already} | "
          f"{'recoverable' if args.dry_run else 'recovered'} from CIF: "
          f"{n_have_cif if args.dry_run else n_recovered}")
    print(f"  per target: {per_target}")
    print(f"  no usable CIF (need ESMFold2): {sum(missing.values())} {missing}")
    if args.dry_run:
        return

    if args.rebuild_hof:
        # Targets come from the pool itself, not from per_target — that only counts
        # rows recovered on THIS pass, and is empty on a re-run where everything is
        # already scored (exactly when you still want the HOF rebuilt).
        targets = sorted(seen_targets)
        if not targets:
            print("HOF rebuild skipped: no target_name values found in the round summaries.")
            return
        try:
            import iterative_refinement as ir
            ir.OUT_BASE = out_base
            refiner = ir.IterativeRefiner(active_targets=targets, active_loops=args.loops,
                                          state_path=out_base / "refinement_state.json")
            scored = []
            for csv in csvs:
                d = pd.read_csv(csv)
                d = d[pd.to_numeric(d.get("esm_plddt"), errors="coerce").notna()]
                scored.extend(d.to_dict("records"))
            for t in targets:
                refiner.state["hof"][t] = []
            refiner.update_hof(scored)
            refiner._save_state()
            sizes = {t: len(refiner.state["hof"].get(t, [])) for t in targets}
            print(f"HOF rebuilt from {len(scored)} scored designs. Sizes: {sizes}")
            # The whole point of the rebuild is to seed HOF backbone reuse, so verify
            # the backbones those entries point at actually exist on disk.
            import re as _re
            ok = tot = 0
            for t in targets:
                for e in refiner.state["hof"].get(t, []):
                    m = _re.search(r"_it(\d+)_d(\d+)", str(e.get("design_id", "")))
                    if not m:
                        continue
                    tot += 1
                    if (out_base / f"it_{m.group(1)}" / "rfd3" / t /
                            f"bb_{m.group(2)}.cif").exists():
                        ok += 1
            print(f"HOF backbones resolvable for reuse: {ok}/{tot}"
                  + ("" if ok == tot else "  <-- missing CIFs will fall back to fresh RFd3"))
        except Exception as exc:
            print(f"HOF rebuild FAILED ({exc}).")
            print("  Backbone reuse needs a populated HOF — do not start the new run until "
                  "this succeeds.")
            return

    print("\nNext: python Generation/rescore_pool.py --dry-run   "
          "(should now list only the designs with no CIF).")


if __name__ == "__main__":
    main()
