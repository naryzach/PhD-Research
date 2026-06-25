"""Shared helpers for the ESM-C binder-classification pipeline.

Kept dependency-light so every script (data_prep / train / predict_*) can import
from here without pulling in heavyweight modules it does not need.
"""
from __future__ import annotations

import os
import random
from pathlib import Path
from typing import Dict, List, Optional

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


# --------------------------------------------------------------------------- #
# Compatibility shims
# --------------------------------------------------------------------------- #
def patch_dynamo_config_compat() -> None:
    """Tolerate ESM++ remote code on older torch (needed for V100/sm_70 wheels).

    ESM++'s ``modeling_esm_plusplus.py`` sets ``torch._dynamo.config`` knobs such
    as ``recompile_limit`` that don't exist on torch versions old enough to still
    ship Volta (sm_70) kernels. We never use torch.compile, so make assignments of
    unknown knobs no-ops instead of letting them crash the model import.
    """
    try:
        import torch._dynamo.config as dcfg
        cls = type(dcfg)
        if getattr(cls, "_esmc_lenient", False):
            return
        _orig_setattr = cls.__setattr__

        def _lenient(self, name, value):
            try:
                _orig_setattr(self, name, value)
            except AttributeError:
                pass  # unknown dynamo knob on this torch version -> ignore

        cls.__setattr__ = _lenient
        cls._esmc_lenient = True
    except Exception:
        pass  # never let the shim itself break loading


# --------------------------------------------------------------------------- #
# Tokenizer
# --------------------------------------------------------------------------- #
def get_tokenizer(model_id: str):
    patch_dynamo_config_compat()
    """Return the ESM++/ESM-C tokenizer.

    ESM++ registers a custom tokenizer; depending on the transformers version it
    may or may not be reachable through AutoTokenizer, so we fall back to pulling
    it off the model object.
    """
    from transformers import AutoTokenizer

    try:
        return AutoTokenizer.from_pretrained(model_id, trust_remote_code=True)
    except Exception:
        from transformers import AutoModel

        model = AutoModel.from_pretrained(model_id, trust_remote_code=True)
        tok = getattr(model, "tokenizer", None)
        if tok is None:
            raise RuntimeError(
                f"Could not obtain a tokenizer for '{model_id}'. "
                "Check the model id and that trust_remote_code is allowed."
            )
        return tok


def detect_bos_offset(tokenizer) -> int:
    """How many special tokens precede residue 0 (ESM prepends a single BOS/CLS).

    Returns the index of the first real residue token so that, for residue-level
    tokenization, sequence char ``i`` maps to token ``i + offset``.
    """
    enc = tokenizer("ACDEFG", add_special_tokens=True)
    ids = enc["input_ids"]
    # ESM tokenizers map each residue to exactly one id and wrap with BOS/EOS.
    # offset = number of leading special tokens (1 for ESM/ESM-C).
    offset = 1 if len(ids) == len("ACDEFG") + 2 else 0
    return offset


# --------------------------------------------------------------------------- #
# Loop location (the 6-aa grafted region inside the constant scaffold)
# --------------------------------------------------------------------------- #
def locate_loop(seq: str, loop: Optional[str]) -> int:
    """Return the 0-based char index of ``loop`` inside ``seq`` (-1 if unknown)."""
    if not isinstance(seq, str) or not isinstance(loop, str) or not loop:
        return -1
    return seq.find(loop)


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #
def compute_target_metrics(y_true: np.ndarray, y_score: np.ndarray,
                           threshold: float = 0.5) -> Dict[str, float]:
    """Per-target binary metrics robust to single-class slices."""
    from sklearn.metrics import (
        average_precision_score,
        f1_score,
        matthews_corrcoef,
        roc_auc_score,
    )

    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    n = int(y_true.shape[0])
    n_pos = int(y_true.sum())
    out = {
        "n": n,
        "n_pos": n_pos,
        "pos_rate": (n_pos / n) if n else float("nan"),
        "pr_auc": float("nan"),
        "roc_auc": float("nan"),
        "mcc": float("nan"),
        "f1": float("nan"),
    }
    if n == 0:
        return out
    y_pred = (y_score >= threshold).astype(int)
    # AUCs need both classes present.
    if 0 < n_pos < n:
        out["pr_auc"] = float(average_precision_score(y_true, y_score))
        out["roc_auc"] = float(roc_auc_score(y_true, y_score))
    out["mcc"] = float(matthews_corrcoef(y_true, y_pred)) if len(np.unique(y_pred)) > 1 or n_pos else 0.0
    out["f1"] = float(f1_score(y_true, y_pred, zero_division=0))
    return out


def best_threshold(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Threshold on ``y_score`` that maximises MCC (falls back to 0.5)."""
    from sklearn.metrics import matthews_corrcoef

    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    if y_true.size == 0 or y_true.sum() in (0, y_true.size):
        return 0.5
    candidates = np.unique(np.concatenate([[0.0, 1.0], y_score]))
    best_t, best_m = 0.5, -2.0
    for t in candidates:
        m = matthews_corrcoef(y_true, (y_score >= t).astype(int))
        if m > best_m:
            best_m, best_t = m, float(t)
    return best_t


def mean_metric(per_target: Dict[str, Dict[str, float]], key: str) -> float:
    """Mean of a metric across targets, ignoring NaNs."""
    vals = [m[key] for m in per_target.values() if m.get(key) == m.get(key)]  # drop NaN
    return float(np.mean(vals)) if vals else float("nan")
