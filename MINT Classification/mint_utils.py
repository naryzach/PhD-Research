"""Shared helpers for the MINT binder-classification pipeline.

Self-contained (mirrors ../ESM-C Classification/esmc_utils.py's config/seed/
metric helpers, which are backbone-agnostic) so this folder has no import
dependency on the ESM-C pipeline. MINT-specific tokenization lives in model.py
via the `mint` package itself (mint.helpers.extract), not here.
"""
from __future__ import annotations

import os
import random
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import yaml


# --------------------------------------------------------------------------- #
# Config / IO
# --------------------------------------------------------------------------- #
def load_config(path: str | os.PathLike) -> dict:
    """Load the YAML config and resolve paths relative to the config file."""
    path = Path(path).resolve()
    with open(path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    cfg["_config_dir"] = str(path.parent)
    return cfg


def resolve_path(cfg: dict, p: str | os.PathLike) -> Path:
    """Resolve a (possibly relative) path against the config file's directory."""
    p = Path(p)
    if p.is_absolute():
        return p
    return (Path(cfg["_config_dir"]) / p).resolve()


def set_seed(seed: int = 42) -> None:
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass


def get_device():
    import torch

    if torch.cuda.is_available():
        return torch.device("cuda")
    if getattr(torch.backends, "mps", None) is not None and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def locate_loop(seq: str, loop: Optional[str]) -> int:
    """Return the 0-based char index of ``loop`` inside ``seq`` (-1 if unknown)."""
    if not isinstance(seq, str) or not isinstance(loop, str) or not loop:
        return -1
    return seq.find(loop)


# --------------------------------------------------------------------------- #
# Metrics (identical semantics to the ESM-C pipeline's esmc_utils.py)
# --------------------------------------------------------------------------- #
def compute_metrics(y_true: np.ndarray, y_score: np.ndarray,
                    threshold: float = 0.5, beta: float = 0.5) -> Dict[str, float]:
    """Binary metrics robust to single-class slices.

    ``fbeta`` (default beta=0.5) weights precision over recall -- useful when
    false positives (calling a non-binder a binder) are costlier than false
    negatives, e.g. picking candidates for wet-lab testing.
    """
    from sklearn.metrics import (
        average_precision_score,
        f1_score,
        fbeta_score,
        matthews_corrcoef,
        roc_auc_score,
    )

    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    n = int(y_true.shape[0])
    n_pos = int(y_true.sum())
    out = {
        "n": n, "n_pos": n_pos, "pos_rate": (n_pos / n) if n else float("nan"),
        "pr_auc": float("nan"), "roc_auc": float("nan"), "mcc": float("nan"),
        "f1": float("nan"), "fbeta": float("nan"),
    }
    if n == 0:
        return out
    y_pred = (y_score >= threshold).astype(int)
    if 0 < n_pos < n:
        out["pr_auc"] = float(average_precision_score(y_true, y_score))
        out["roc_auc"] = float(roc_auc_score(y_true, y_score))
    out["mcc"] = float(matthews_corrcoef(y_true, y_pred)) if len(np.unique(y_pred)) > 1 or n_pos else 0.0
    out["f1"] = float(f1_score(y_true, y_pred, zero_division=0))
    out["fbeta"] = float(fbeta_score(y_true, y_pred, beta=beta, zero_division=0))
    return out


def best_threshold(y_true: np.ndarray, y_score: np.ndarray,
                   metric: str = "mcc", beta: float = 0.5) -> float:
    """Threshold on ``y_score`` that maximises ``metric`` (falls back to 0.5).

    ``metric`` is ``"mcc"`` or ``"fbeta"`` (precision-weighted when beta < 1).
    """
    from sklearn.metrics import fbeta_score, matthews_corrcoef

    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    if y_true.size == 0 or y_true.sum() in (0, y_true.size):
        return 0.5
    candidates = np.unique(np.concatenate([[0.0, 1.0], y_score]))
    best_t, best_m = 0.5, -2.0
    for t in candidates:
        pred = (y_score >= t).astype(int)
        m = (fbeta_score(y_true, pred, beta=beta, zero_division=0) if metric == "fbeta"
             else matthews_corrcoef(y_true, pred))
        if m > best_m:
            best_m, best_t = m, float(t)
    return best_t


def mean_metric(per_target: Dict[str, Dict[str, float]], key: str) -> float:
    """Mean of a metric across per-target-name metric dicts, ignoring NaNs."""
    vals = [m[key] for m in per_target.values() if m.get(key) == m.get(key)]  # drop NaN
    return float(np.mean(vals)) if vals else float("nan")
