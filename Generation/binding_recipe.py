#!/usr/bin/env python3
"""
binding_recipe.py  —  the multi-term design-scoring recipe (v1).

Canonical location: Generation/ (tracked). Calibrated in Local/Calibration/
(gitignored); see Local/Calibration/README.md for the evidence and re-fit steps.

A deliberately MULTI-FACTOR score (no lone metric), each term justified by the
Dec-2025 prediction-vs-result calibration:

  1. affinity   — interface/fold confidence on the ON-target, from the only AF
                  metrics with REAL dynamic range: LpLDDT (loop pLDDT) and PAE
                  (inverted). pTM-family (ipTM/ApTM/BpTM/pTM) is nearly flat within
                  a target → its rank is noise. Backbone term.
  2. selectivity— ON minus mean(OFF) gap on that same confidence. The project goal
                  is SELECTIVE inhibitors and raw affinity is dominated by a
                  non-specific "stickiness" factor, so we reward target
                  discrimination — but as ONE term, never alone.
  3. independent— BpTM on the ON-target: the AF metric most INDEPENDENT of the
                  LpLDDT/PAE axis, binding-leaning, but flat → low weight, flagged.
  4. expressibility — ipTM on the ON-target. In our data ipTM tracked EXPRESSION,
                  not binding; kept as a small developability prior, not a binding
                  term. (Calibration was indifferent to it for binding.)

Weights are reasoned (NOT least-squares — at n=12 a fitted blend cross-validated
WORSE than this one). The score is ROBUST to missing metrics: any term whose
input metric is absent is dropped and the remaining weights are renormalized, so
the same function works on AF3 (no BpTM) and on ESMFold2 (ipTM/pTM/pLDDT only).

CAVEAT — predictor scale. Calibrated on AF-server metrics (== AF3, same model, so
this recipe is valid for AF3 directly). ESMFold2's metrics are a DIFFERENT scale;
score_with_esmfold2.py now derives esm_lplddt / esm_pae over the redesigned loops
(so all four terms can run on ESMFold2), but re-derive NORM_RANGES — and ideally
re-run the calibration — before trusting this recipe as the ESMFold2 ranker.
"""
from __future__ import annotations
import numpy as np

# (lo, hi, higher_is_better) → maps a raw AF metric to [0, 1].
NORM_RANGES = {
    "LpLDDT": (50.0, 90.0, True),
    "PAE":    (5.0, 22.0, False),
    "ipTM":   (0.50, 0.90, True),
    "BpTM":   (0.84, 0.93, True),
    "ApTM":   (0.79, 0.84, True),
    "pTM":    (0.79, 0.90, True),
    "pLDDT":  (80.0, 90.0, True),
}

RECIPE_WEIGHTS = {
    "affinity":       0.40,   # real-range interface/fold confidence (LpLDDT, -PAE)
    "selectivity":    0.30,   # on - mean(off) gap on that confidence (NOT lone)
    "independent":    0.12,   # BpTM, independent axis, flagged low weight
    "expressibility": 0.18,   # ipTM as a developability prior, not a binding term
}


def _present(metrics: dict, key: str) -> bool:
    v = metrics.get(key) if metrics else None
    return v is not None and not (isinstance(v, float) and np.isnan(v))


def normalize(metric: str, value) -> float:
    """Map a raw AF metric to [0, 1] (missing → 0.0)."""
    if metric not in NORM_RANGES or value is None or (isinstance(value, float) and np.isnan(value)):
        return 0.0
    lo, hi, up = NORM_RANGES[metric]
    v = (float(value) - lo) / (hi - lo)
    v = v if up else 1.0 - v
    return float(min(1.0, max(0.0, v)))


def _interface_confidence(metrics: dict):
    """(value, available?) — mean of the real-range confidences present (LpLDDT, -PAE)."""
    parts = [normalize(k, metrics.get(k)) for k in ("LpLDDT", "PAE") if _present(metrics, k)]
    return (float(np.mean(parts)), True) if parts else (0.0, False)


def _terms(on_metrics: dict, off_metrics_list: list):
    """Return {term: (value, available?)} for the four recipe terms."""
    ic_on, ic_ok = _interface_confidence(on_metrics)
    ic_offs = [c for c in (_interface_confidence(o) for o in (off_metrics_list or [])) if c[1]]
    mean_off = float(np.mean([v for v, _ in ic_offs])) if ic_offs else None

    affinity = (ic_on, ic_ok)
    sel_ok = ic_ok and mean_off is not None
    selectivity = (((ic_on - mean_off + 1.0) / 2.0) if sel_ok else 0.0, sel_ok)
    independent = (normalize("BpTM", on_metrics.get("BpTM")), _present(on_metrics, "BpTM"))
    expressibility = (normalize("ipTM", on_metrics.get("ipTM")), _present(on_metrics, "ipTM"))
    return {"affinity": affinity, "selectivity": selectivity,
            "independent": independent, "expressibility": expressibility}


def score_design(on_metrics: dict, off_metrics_list: list, weights: dict = None) -> float:
    """
    Composite design score in [0, 1] for binding ON-target while avoiding OFF-targets.
    Terms whose inputs are missing are dropped and weights renormalized.

    on_metrics       : {metric: value} for the on-target complex.
    off_metrics_list : list of {metric: value}, one per OFF-target protease (use ALL
                       predicted off-targets in production, not only the tested ones).
    """
    w = weights or RECIPE_WEIGHTS
    terms = _terms(on_metrics, off_metrics_list)
    num = sum(w[k] * val for k, (val, ok) in terms.items() if ok)
    den = sum(w[k] for k, (_, ok) in terms.items() if ok)
    return float(num / den) if den > 0 else float("nan")


def score_breakdown(on_metrics: dict, off_metrics_list: list, weights: dict = None) -> dict:
    """Per-term values + weighted contributions (renormalized over available terms)."""
    w = weights or RECIPE_WEIGHTS
    terms = _terms(on_metrics, off_metrics_list)
    den = sum(w[k] for k, (_, ok) in terms.items() if ok) or 1.0
    out = {}
    for k, (val, ok) in terms.items():
        out[f"term_{k}"] = val if ok else float("nan")
        out[f"wcontrib_{k}"] = (w[k] * val / den) if ok else 0.0
    out["terms_used"] = [k for k, (_, ok) in terms.items() if ok]
    out["score"] = sum(out[f"wcontrib_{k}"] for k in terms)
    return out


if __name__ == "__main__":
    import json
    on = {"LpLDDT": 86.0, "PAE": 6.0, "BpTM": 0.91, "ipTM": 0.87}
    offs = [{"LpLDDT": 70.0, "PAE": 12.0}, {"LpLDDT": 68.0, "PAE": 13.0}]
    print("AF3-style (all terms):"); print(json.dumps(score_breakdown(on, offs), indent=2))
    esm_on = {"LpLDDT": 84.0, "ipTM": 0.6, "pTM": 0.62}   # ESMFold2: no PAE/BpTM
    esm_off = [{"LpLDDT": 80.0}, {"LpLDDT": 78.0}]
    print("ESMFold2-style (reduced terms):"); print(json.dumps(score_breakdown(esm_on, esm_off), indent=2))
