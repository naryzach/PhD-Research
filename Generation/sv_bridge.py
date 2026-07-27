"""
sv_bridge.py — Structural-Validation interface battery for the binder-generation
pipeline.

Pulls the full interface-metric battery (Structural-Validation/analysis/metrics.py)
onto each generation candidate's predicted-complex CIF and returns them as `sv_*`
columns for round_summary.csv.

DELIBERATELY LOG-ONLY. None of these metrics enter calc_composite / the ranking.
The FCS correlation (Structural-Validation reports/fcs_correlation_report.html)
showed they do not predict measured binding — they anti- or zero-correlate — so
they are recorded for transparency and to accumulate the exact feature set that
correlation consumes, so future wet-lab rounds can be folded straight into it.

The ONE structural signal allowed to act on selection is a mechanistic FILTER,
not a metric-in-the-composite: `occlusion_pass` asks whether the binder's reactive
edge actually reaches the catalytic zinc cleft (buries a minimum fraction of the
HExxHxxGxxH loop). That is a sanity gate analogous to the existing ESMFold2
foldability floor, not a ranking term.

DockQ / LRMS are computed only against the wild-type TIMP3:target complex (an
APPROXIMATE reference) and are flagged `sv_dockq_ref`; they measure "does this
design bind like wild-type", so they are informative to log but would penalise
novel binding modes if ranked — never gate on them.

Robust by construction: any failure (missing CIF, unreadable, SV import problem)
returns {} / None, exactly like the pipeline's existing _esm_iface_feats, so it can
never break a generation run.

CLI (batch post-hoc enrichment of an existing round_summary):
    python sv_bridge.py --round-summary it_3/round_summary.csv \
        --cif-col esm_cif --target-col target [--out enriched.csv]
"""
from __future__ import annotations

import sys
from pathlib import Path

# Make the Structural-Validation analysis stack importable.
_SV = Path(__file__).resolve().parents[1] / "Structural-Validation"
for _p in (_SV, _SV / "analysis", _SV / "utils"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

_OK = True
try:
    import metrics as _M          # noqa: E402
    import structure_io as _sio   # noqa: E402
except Exception:                 # SV stack unavailable -> bridge degrades to no-op
    _OK = False

BINDER_CHAIN = "A"
TARGET_CHAIN = "B"
SV = "sv_"                        # column namespace
DEFAULT_MIN_OCCLUSION = 0.15     # "reactive edge reaches the catalytic cleft" gate


def available() -> bool:
    """True if the Structural-Validation metric stack imported successfully."""
    return _OK


def _resolve_chains(chains, binder_chain, target_chain):
    if binder_chain in chains and target_chain in chains:
        return chains[binder_chain], chains[target_chain]
    if len(chains) < 2:
        return None, None
    ids = sorted(chains, key=lambda k: len(chains[k].seq))   # binder = shorter
    return chains[ids[0]], chains[ids[-1]]


def sv_battery(cif_path, binder_chain: str = BINDER_CHAIN,
               target_chain: str = TARGET_CHAIN, target_id: str | None = None) -> dict:
    """Full interface battery for one predicted-complex CIF -> {sv_*: value}.

    Reference-free geometry/chemistry + pLDDT-based confidence always; DockQ/LRMS
    only when a wild-type reference for `target_id` exists (flagged approximate).
    Returns {} on any problem.
    """
    if not _OK or not cif_path or not Path(str(cif_path)).exists():
        return {}
    try:
        chains = _sio.get_chains(cif_path)
        b, t = _resolve_chains(chains, binder_chain, target_chain)
        if b is None:
            return {}
        motif = _sio.zinc_motif_resids(t)
        im = _M.interface_summary(cif_path, b.cid, t.cid,
                                  motif_resids_b=motif or None)
        out = {SV + k: v for k, v in {
            "bsa": im.bsa, "bsa_polar": im.bsa_polar, "bsa_apolar": im.bsa_apolar,
            "f_apolar_bsa": im.f_apolar_bsa, "n_iface_res_binder": im.n_iface_res_A,
            "n_iface_res_target": im.n_iface_res_B, "n_contacts_5A": im.n_contacts_5A,
            "n_hbonds": im.n_hbonds, "n_salt_bridges": im.n_salt_bridges,
            "n_hydrophobic": im.n_hydrophobic, "n_interface_clashes": im.n_interface_clashes,
            "contact_density": im.contact_density,
            "n_buried_unsat_hbond": im.n_buried_unsat_hbond,
            "sc_shape_complementarity": im.sc_shape_complementarity,
            "charge_complementarity": im.charge_complementarity,
            "catalytic_occlusion": im.catalytic_occlusion,
            "zinc_bsa_buried": im.zinc_bsa_buried,
        }.items()}

        # confidence: pLDDT-based (from CIF B-factor) always; PAE-derived only if an
        # AF3 full_data JSON sits beside the CIF (co-fold), else blank.
        cm = _M.confidence_metrics(b, t, model_path=str(cif_path),
                                   cid=b.cid, tid=t.cid)
        for k in ("interface_plddt", "pdockq", "pdockq2", "interface_pae", "lis"):
            v = cm.get(k)
            if v not in ("", None):
                out[SV + k] = v

        # optional DockQ vs the wild-type reference (approximate; LOG ONLY)
        if target_id:
            _dockq_vs_wt(cif_path, b, t, target_id, out)
        return out
    except Exception:
        return {}


def _dockq_vs_wt(cif_path, b, t, target_id, out) -> None:
    try:
        import model_registry as _MR          # noqa: E402
        from complex_metrics import assign_chains  # noqa: E402
        ref, ref_type = _MR.complex_reference(target_id)
        if not ref:
            return
        rc = _sio.get_chains(ref)
        rcc, rtc = assign_chains(rc, b.seq, t.seq)
        dq = _M.dockq(cif_path, ref, t.cid, b.cid, rtc.cid, rcc.cid)
        out.update({SV + "dockq_vs_wt": dq.dockq, SV + "lrms_vs_wt": dq.lrms,
                    SV + "fnat_vs_wt": dq.fnat, SV + "dockq_ref": ref_type})
    except Exception:
        pass


def occlusion_pass(sv_dict: dict, min_occlusion: float = DEFAULT_MIN_OCCLUSION):
    """Mechanistic sanity FILTER (not a ranking term): does the binder's reactive
    edge reach the catalytic cleft? Returns True/False, or None when the motif was
    not found / occlusion is undefined (fail-open — never drop on unknown)."""
    v = sv_dict.get(SV + "catalytic_occlusion")
    try:
        fv = float(v)
    except (TypeError, ValueError):
        return None
    if fv != fv:              # NaN: motif not found / occlusion undefined
        return None           # fail-open — never drop a design on unknown geometry
    return fv >= min_occlusion


# ── CLI: batch-enrich an existing round_summary.csv ───────────────────────────
def _cli() -> None:
    import argparse
    import pandas as pd

    ap = argparse.ArgumentParser(description="Add sv_* battery columns to a "
                                             "round_summary.csv (log-only).")
    ap.add_argument("--round-summary", required=True)
    ap.add_argument("--cif-col", default="esm_cif",
                    help="column holding the predicted-complex CIF path")
    ap.add_argument("--target-col", default="target",
                    help="column holding the target id (for optional DockQ-vs-WT)")
    ap.add_argument("--min-occlusion", type=float, default=DEFAULT_MIN_OCCLUSION)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    if not _OK:
        sys.exit("Structural-Validation metric stack not importable; aborting.")
    df = pd.read_csv(a.round_summary)
    if a.cif_col not in df.columns:
        sys.exit(f"column {a.cif_col!r} not in {a.round_summary}")

    rows = []
    for _, r in df.iterrows():
        tid = str(r[a.target_col]).upper() if a.target_col in df.columns else None
        sv = sv_battery(r.get(a.cif_col), target_id=tid)
        sv[SV + "occlusion_pass"] = occlusion_pass(sv, a.min_occlusion)
        rows.append(sv)
    enriched = pd.concat([df.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    out = a.out or a.round_summary.replace(".csv", "_sv.csv")
    enriched.to_csv(out, index=False)
    scored = sum(1 for x in rows if x.get(SV + "bsa") is not None)
    npass = sum(1 for x in rows if x.get(SV + "occlusion_pass") is True)
    print(f"enriched {scored}/{len(df)} rows with sv_* battery; "
          f"{npass} pass the occlusion gate (>= {a.min_occlusion}). -> {out}")


if __name__ == "__main__":
    _cli()
