#!/usr/bin/env python3
"""
calibrated_scoring.py  —  in-silico -> binding priors, calibrated on the
2026-07 EXACT-sequence, purchased-only calibration (Local/TIMP3_Redesign_2026-07/).

WHAT THE CALIBRATION FOUND (exact Manufacturer sequences, human/Enzo, n=38
folds / 13 constructs, expression-controlled partial Spearman + strong/weak AUC):

  * NO in-silico metric reliably predicts binding (best AUC ~= 0.68, all |rho|
    <~ 0.3, none significant). We are encoding WEAK PRIORS, not a real ranker.
  * The consistently POSITIVE AF3 signals against the honest (expression-
    normalized) binding readouts -- Norm Bind Med (Expr+), % of Pos Ctrl,
    Norm Pos Med Ratio -- are, in order:
        af3_plddt (binder-chain pLDDT)        rho up to +0.33  (p~0.05)
        af3_aptm  (binder-chain pTM)          rho ~ +0.26 across all readouts
        af3 interface residue count           AUC 0.68 (best single)
        af3 interface-residue pLDDT           rho ~ +0.17
        af3_iptm                              rho ~ +0.14 (weak, expression-linked)
    All are foldability / interface-geometry flavored and POSITIVE.
  * BpTM is DEBUNKED: 96% target-determined (within-target std ~0.03) and its
    earlier +0.48 did NOT replicate on exact folds (-> ~-0.25). DROP it.
  * Global pTM and loop-interface PAE are ~0 against the honest readouts. Excluded.

THE ESMFold2 SELECTION-BIAS TRAP (why ESMFold2 is a FILTER, not a binding ranker):
  ESMFold2's confidence metrics went NEGATIVE against binding on the tested set
  (esm_iptm rho -0.24, esm_lplddt -0.16, esm_ptm -0.22 vs Norm Bind Med (Expr+)).
  That negative slope is almost certainly a range-restriction / collider artifact:
  we only observe constructs that were actually MADE (all of which fold), so we
  never see the non-folding non-binders that would sit at even lower ESMFold2
  scores. Rewarding LOW ESMFold2 confidence would be selecting on a spurious
  in-sample slope and would likely hurt on new designs. So:
    - esm_iptm / esm_lplddt / esm_ptm are NEVER used to select or rank (neither
      direction). They are recorded for transparency only.
    - ESMFold2 is used only as a FOLDABILITY FLOOR (its intended positive sense:
      "did the binder fold and form an interface at all") on esm_plddt, which is
      ~0-correlated with binding = neutral = safe to threshold.
    - The ONE ESMFold2 feature that is POSITIVE and mechanistic is interface
      contact density (rho +0.28 / AUC 0.67). It is the ESMFold2-stage binding-ish
      ranker. Nothing else from ESMFold2 touches ranking.

SCOPE: calibrated on AF-server metrics, i.e. AF3 directly -> AF3 priors are valid.
ESMFold2 lives on a different scale; only the contact-density prior + foldability
floor are used from it, both deliberately coarse. Weights are REASONED, not fitted
(at n=13 a least-squares blend cross-validated worse than a reasoned one). Every
prior is robust to missing inputs: absent terms drop out and the remaining weights
renormalize, so the same functions work on full AF3, legacy AF3 (no interface
features), and ESMFold2.

Sources: Local/TIMP3_Redesign_2026-07/docs/09_exact_calibration_final.md,
data/exact_partial_correlation.csv, data/interface_exact_vs_binding.csv.
"""
from __future__ import annotations
import numpy as np

# ── Normalization ranges (lo, hi, higher_is_better) ──────────────────────────
# Bounds are the ~min/max of the exact human/Enzo distributions (n=38), so a raw
# metric maps to roughly [0,1] across the observed design space.
AF3_NORM = {
    "af3_plddt":               (60.0, 73.0, True),   # binder-chain pLDDT  (best positive)
    "af3_aptm":                (0.55, 0.75, True),   # binder-chain pTM    (most consistent)
    "af3_iptm":                (0.25, 0.60, True),   # interface pTM       (weak / expression)
    "af3_iface_n_iface_res":   (30.0, 55.0, True),   # interface size      (best AUC 0.68)
    "af3_iface_iface_plddt":   (42.0, 60.0, True),   # interface-res pLDDT
}
ESM_NORM = {
    "esm_iface_contact_density": (15.0, 40.0, True), # the one trustworthy positive ESM feature
    # pDockQ from the SV battery: the ONLY local metric validated against AF3.
    # 2026-08-15, ESMFold2-predicted complexes vs AF3 ipTM, target-centred Spearman:
    #   prospective  (it>=38, n=24, pre-registered) rho = +0.655  p = 0.0005
    #   retrospective(it 0-36, n=36)                rho = +0.364  p = 0.029
    #   pooled                    (n=60)            rho = +0.434  p = 0.0005
    # The ESM composite on the same designs: rho = +0.112, p = 0.40 (i.e. nothing).
    # Range = p10-p90 of the 3,325-design pool (0.119 / 0.443), rounded.
    "sv_pdockq": (0.12, 0.44, True),
}

# ── AF3 binding-prior weights (positive foldability/interface terms only) ────
# Anchored on the honest (expression-normalized) binding readouts. BpTM dropped,
# global pTM and loop-PAE excluded (~0 correlation). ipTM kept small as a
# developability/expression prior, not a binding term.
AF3_BINDING_WEIGHTS = {
    "af3_plddt":               0.30,
    "af3_aptm":                0.25,
    "af3_iface_n_iface_res":   0.25,
    "af3_iface_iface_plddt":   0.10,
    "af3_iptm":                0.10,
}

# ── ESMFold2 foldability floor (a GARBAGE filter, NOT a within-range ranker) ──
# Gate only on esm_plddt (binding-neutral, safe to threshold) and on the mere
# EXISTENCE of an interface. We deliberately do NOT gate on esm_iptm/esm_lplddt
# (the negatively-correlated confidences) to avoid the collider bias above.
# Tunable: check the esm_plddt distribution in a real round_summary.csv before
# tightening. Set FOLD_GATE["esm_plddt_min"]=None to disable the pLDDT floor.
FOLD_GATE = {
    "esm_plddt_min":  55.0,   # binder folds at all (ESMFold2 whole-binder pLDDT, 0-100)
    "min_iface_res":  5,      # an interface actually forms (>= this many contact residues)
}

# ── ESMFold2-stage composite (filter + the one positive feature) ─────────────
# No esm_iptm term on purpose. fold_base is awarded for passing the foldability
# floor; contact_density is the ESMFold2 binding-ish ranker; esm_plddt is a mild
# (positive-sense) developability prior.
COMPOSITE_ESMFOLD2 = {
    "fold_base":       0.25,   # was 0.50 -- still awarded once past the foldability gate
    "pdockq":          0.50,   # NEW dominant term: the only AF3-validated local metric
    "contact_density": 0.15,   # was 0.35 -- kept, but it is not the validated axis
    "plddt":           0.10,   # was 0.15 -- mild developability prior
}
ESM_GATE_FAIL_SCORE = 0.05    # designs failing the foldability floor rank below all folded ones


# ── Metric normalization ──────────────────────────────────────────────────────

def _isnum(v) -> bool:
    return v is not None and not (isinstance(v, float) and np.isnan(v))


def _norm(table: dict, key: str, value) -> float:
    """Map a raw metric to [0,1] via `table` (missing/out-of-table -> nan)."""
    if key not in table or not _isnum(value):
        return float("nan")
    lo, hi, up = table[key]
    v = (float(value) - lo) / (hi - lo)
    v = v if up else 1.0 - v
    return float(min(1.0, max(0.0, v)))


def _weighted(norm_table: dict, weights: dict, metrics: dict) -> float:
    """
    Weighted mean of normalized terms present in `metrics`; absent terms drop and
    the remaining weights renormalize. Returns nan if no term is available.
    """
    num = den = 0.0
    for key, w in weights.items():
        nv = _norm(norm_table, key, metrics.get(key))
        if not np.isnan(nv):
            num += w * nv
            den += w
    return float(num / den) if den > 0 else float("nan")


# ── AF3 binding prior ─────────────────────────────────────────────────────────

def af3_binding_prior(metrics: dict, weights: dict = None) -> float:
    """
    Weak AF3 binding prior in [0,1] from the calibrated positive metrics.
    `metrics` may contain any of: af3_plddt, af3_aptm, af3_iptm,
    af3_iface_n_iface_res, af3_iface_iface_plddt. Missing terms renormalize away,
    so this also scores legacy AF3 entries that only have plddt/iptm.
    """
    return _weighted(AF3_NORM, weights or AF3_BINDING_WEIGHTS, metrics)


def af3_binding_breakdown(metrics: dict, weights: dict = None) -> dict:
    """Per-term normalized values + renormalized weighted contributions (for logs/CSV)."""
    w = weights or AF3_BINDING_WEIGHTS
    present = {k: _norm(AF3_NORM, k, metrics.get(k)) for k in w}
    den = sum(w[k] for k, nv in present.items() if not np.isnan(nv)) or 1.0
    out = {}
    for k, nv in present.items():
        out[f"term_{k}"] = nv
        out[f"wcontrib_{k}"] = (w[k] * nv / den) if not np.isnan(nv) else 0.0
    out["af3_binding_prior"] = sum(out[f"wcontrib_{k}"] for k in w)
    return out


# ── ESMFold2 foldability floor + stage composite ─────────────────────────────

def esm_passes_fold_gate(metrics: dict, gate: dict = None) -> bool:
    """
    True if the ESMFold2 prediction clears the foldability floor (binder folds and
    an interface forms). Uses esm_plddt (binding-neutral) and interface size only
    -- never esm_iptm/esm_lplddt. Missing inputs are treated leniently (pass), so a
    round with no interface features saved still runs on the pLDDT floor alone.
    """
    g = gate or FOLD_GATE
    plddt = metrics.get("esm_plddt")
    if g.get("esm_plddt_min") is not None and _isnum(plddt) and float(plddt) < g["esm_plddt_min"]:
        return False
    nres = metrics.get("esm_iface_n_iface_res")
    if g.get("min_iface_res") is not None and _isnum(nres) and int(nres) < g["min_iface_res"]:
        return False
    return True


def esmfold2_stage_score(metrics: dict, weights: dict = None, gate: dict = None) -> float:
    """
    ESMFold2-stage ranking score in [0,1]: foldability floor + the one positive
    ESMFold2 feature (interface contact density) + a mild pLDDT developability
    prior. esm_iptm is intentionally excluded from ranking. Designs failing the
    foldability floor collapse to ESM_GATE_FAIL_SCORE so they rank below folded ones.
    """
    w = weights or COMPOSITE_ESMFOLD2
    if not esm_passes_fold_gate(metrics, gate):
        return ESM_GATE_FAIL_SCORE
    cd = _norm(ESM_NORM, "esm_iface_contact_density", metrics.get("esm_iface_contact_density"))
    # pDockQ is the validated term (see ESM_NORM). Designs without it -- anything
    # predating --sv-battery and not yet backfilled by recover_from_cifs.py --sv --
    # yield NaN here and the remaining weights renormalise, reproducing the old
    # composite exactly. So a partially backfilled pool degrades gracefully rather
    # than ranking un-backfilled designs as zero.
    pq = _norm(ESM_NORM, "sv_pdockq", metrics.get("sv_pdockq"))
    plddt = metrics.get("esm_plddt")
    pl = min(1.0, max(0.0, float(plddt) / 100.0)) if _isnum(plddt) else float("nan")
    num, den = w["fold_base"], w["fold_base"]        # fold_base always awarded once gated in
    for key, val in (("pdockq", pq), ("contact_density", cd), ("plddt", pl)):
        if not np.isnan(val):
            num += w[key] * val
            den += w[key]
    return float(num / den) if den > 0 else float(w["fold_base"])


# ── Interface geometry from a complex structure ──────────────────────────────
# Ported from Local/TIMP3_Redesign_2026-07/interface_exact.py so the SAME
# interface features used in calibration are computed inside the pipeline.

def _iface_from_coords(A, Ares, Ab, B, cutoff=5.0) -> dict:
    """Core interface geometry from binder(A) / target(B) coordinate arrays."""
    if len(A) == 0 or len(B) == 0:
        return {}
    contacts = 0
    near = np.zeros(len(A), dtype=bool)
    mind = 1e9
    for i in range(0, len(B), 400):                       # chunk to bound memory
        d = np.sqrt(((A[:, None, :] - B[None, i:i + 400, :]) ** 2).sum(-1))
        contacts += int((d < cutoff).sum())
        near |= (d < cutoff).any(1)
        mind = min(mind, float(d.min()))
    nres = int(np.unique(Ares[near]).size)
    ip = float(np.nanmean(Ab[near])) if near.any() and Ab.size else float("nan")
    return {
        "n_iface_res":      nres,
        "n_iface_contacts": contacts,
        "iface_plddt":      ip,
        "contact_density":  (contacts / nres) if nres else float("nan"),
        "min_gap":          mind,
    }


def interface_features_from_atoms(atoms: list, binder_chain="A", target_chain="B",
                                  cutoff: float = 5.0) -> dict:
    """
    Interface features from a list of atom dicts with keys
    {chain, x, y, z, res_id, plddt} (e.g. select_binders_to_order._parse_cif_atoms).
    Returns {n_iface_res, n_iface_contacts, iface_plddt, contact_density, min_gap}.
    """
    A = np.array([(a["x"], a["y"], a["z"]) for a in atoms if a.get("chain") == binder_chain], float)
    Ares = np.array([a["res_id"] for a in atoms if a.get("chain") == binder_chain])
    Ab = np.array([a.get("plddt", np.nan) for a in atoms if a.get("chain") == binder_chain], float)
    B = np.array([(a["x"], a["y"], a["z"]) for a in atoms if a.get("chain") == target_chain], float)
    return _iface_from_coords(A, Ares, Ab, B, cutoff)


def interface_features_from_cif(cif_text: str, binder_chain="A", target_chain="B",
                                cutoff: float = 5.0) -> dict:
    """Interface features straight from mmCIF text (parses _atom_site loop)."""
    cols, in_loop = [], False
    coords_A, res_A, b_A, coords_B = [], [], [], []
    for ln in cif_text.splitlines():
        s = ln.strip()
        if s.startswith("_atom_site."):
            cols.append(s.split(".")[1]); in_loop = True
        elif in_loop and s.startswith(("ATOM", "HETATM")):
            p = s.split()
            if len(p) < len(cols):
                continue
            r = dict(zip(cols, p))
            ch = r.get("auth_asym_id", r.get("label_asym_id"))
            try:
                xyz = (float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]))
            except (ValueError, KeyError):
                continue
            if ch == binder_chain:
                coords_A.append(xyz); b_A.append(float(r.get("B_iso_or_equiv", "nan")))
                try:
                    res_A.append(int(r.get("auth_seq_id", r.get("label_seq_id"))))
                except (ValueError, KeyError):
                    res_A.append(-1)
            elif ch == target_chain:
                coords_B.append(xyz)
        elif in_loop and s.startswith("#"):
            in_loop = False
    return _iface_from_coords(np.array(coords_A, float), np.array(res_A),
                             np.array(b_A, float), np.array(coords_B, float), cutoff)


# ── Self-test (pure python + numpy; runnable without the GPU env) ────────────
if __name__ == "__main__":
    import json
    strong = {"af3_plddt": 72.0, "af3_aptm": 0.74, "af3_iptm": 0.55,
              "af3_iface_n_iface_res": 54, "af3_iface_iface_plddt": 60.0}
    weak = {"af3_plddt": 62.0, "af3_aptm": 0.56, "af3_iptm": 0.26,
            "af3_iface_n_iface_res": 30, "af3_iface_iface_plddt": 43.0}
    legacy = {"af3_plddt": 68.0, "af3_iptm": 0.40}   # no aptm / interface features
    print("AF3 binding prior  strong :", round(af3_binding_prior(strong), 3))
    print("AF3 binding prior  weak   :", round(af3_binding_prior(weak), 3))
    print("AF3 binding prior  legacy :", round(af3_binding_prior(legacy), 3), "(renormalized)")
    print("AF3 breakdown (strong):", json.dumps(af3_binding_breakdown(strong), indent=2))
    esm_ok  = {"esm_plddt": 82.0, "esm_iface_contact_density": 40.0, "esm_iface_n_iface_res": 40}
    esm_bad = {"esm_plddt": 48.0, "esm_iface_contact_density": 40.0, "esm_iface_n_iface_res": 40}
    print("ESM stage score  folded :", round(esmfold2_stage_score(esm_ok), 3),
          "| gate:", esm_passes_fold_gate(esm_ok))
    print("ESM stage score  garbage:", round(esmfold2_stage_score(esm_bad), 3),
          "| gate:", esm_passes_fold_gate(esm_bad))
    assert af3_binding_prior(strong) > af3_binding_prior(weak) > 0
    assert esmfold2_stage_score(esm_ok) > esmfold2_stage_score(esm_bad)
    print("OK: priors ordered as expected.")
