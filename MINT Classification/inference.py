"""Load a trained MintBinderClassifier and run predictions on (sequence,
target) pairs.

Shared by predict_seq.py and predict_csv.py.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np
import pandas as pd
import torch

from mint_utils import get_device
from model import MintBinderClassifier, PairCollator


def load_trained_model(model_dir: str | Path, device=None):
    model_dir = Path(model_dir)
    meta = json.loads((model_dir / "model_meta.json").read_text())
    device = device or get_device()
    model = MintBinderClassifier(
        mint_cfg_json=meta["mint_config_json"], checkpoint_path=meta["mint_checkpoint_path"],
        device=device, sep_chains=meta["sep_chains"], use_multimer=meta["use_multimer"],
        truncation_seq_length=meta["truncation_seq_length"],
        head_hidden=meta["head_hidden"], dropout=meta.get("dropout", 0.1),
        pos_weight=None, freeze_percent=1.0,
    )
    state = torch.load(model_dir / "model_state.pt", map_location="cpu")
    model.load_state_dict(state)
    model.to(device).eval()
    collate = PairCollator(model.mint_collate)
    return model, collate, meta, device


@torch.no_grad()
def predict_pairs(model, collate: PairCollator, seqA: Sequence[str], seqB: Sequence[str],
                  device=None, batch_size: int = 8) -> np.ndarray:
    """Score arbitrary (seqA, seqB) pairs. Returns P(bind), shape (n,)."""
    device = device or next(model.parameters()).device
    probs = []
    for i in range(0, len(seqA), batch_size):
        a_chunk, b_chunk = seqA[i:i + batch_size], seqB[i:i + batch_size]
        batch = collate([{"seqA": a, "seqB": b, "label": 0.0, "weight": 1.0, "target_name": ""}
                        for a, b in zip(a_chunk, b_chunk)])
        chains = batch["chains"].to(device)
        chain_ids = batch["chain_ids"].to(device)
        out = model(chains=chains, chain_ids=chain_ids)
        probs.append(torch.sigmoid(out["logits"].float()).cpu().numpy())
    return np.concatenate(probs)


@torch.no_grad()
def predict_against_known_targets(model, collate: PairCollator, meta: dict,
                                  sequences: List[str], targets: Optional[List[str]] = None,
                                  device=None, batch_size: int = 8) -> pd.DataFrame:
    """Score ``sequences`` against every target the model was trained on (or a
    subset of ``targets``), using the canonical target sequence captured in
    model_meta.json at train time. Wide output, one prob/call column per target
    -- mirrors the ESM-C pipeline's predict_csv.py output shape."""
    targets = targets or meta["targets"]
    tseqs = meta.get("target_sequences", {})
    missing = [t for t in targets if t not in tseqs]
    if missing:
        raise ValueError(f"No known sequence for target(s) {missing} in model_meta.json -- "
                         f"pass an explicit target sequence instead (--seqB / a target-seq column).")

    df = pd.DataFrame({"sequence": list(sequences)})
    thresholds = meta.get("thresholds", {})
    for t in targets:
        probs = predict_pairs(model, collate, list(sequences), [tseqs[t]] * len(sequences),
                              device=device, batch_size=batch_size)
        thr = float(thresholds.get(t, 0.5))
        df[f"prob_{t}"] = probs
        df[f"call_{t}"] = (probs >= thr).astype(int)
    prob_cols = [f"prob_{t}" for t in targets]
    df["max_prob"] = df[prob_cols].max(axis=1)
    df["top_target"] = df[prob_cols].idxmax(axis=1).str.replace("prob_", "", regex=False)
    return df
