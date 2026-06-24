"""Load a trained MultiTaskESMC and run predictions on raw sequences.

Shared by predict_seq.py and predict_csv.py.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader

from esmc_utils import get_device, get_tokenizer, locate_loop
from model import Collator, MultiTaskESMC


def load_trained_model(model_dir: str | Path, device=None):
    model_dir = Path(model_dir)
    meta = json.loads((model_dir / "model_meta.json").read_text())
    device = device or get_device()
    model = MultiTaskESMC(
        model_id=meta["model_id"], targets=meta["targets"],
        pooling=meta["pooling"], dropout=meta.get("dropout", 0.1),
        pos_weight=None,
    )
    state = torch.load(model_dir / "model_state.pt", map_location="cpu")
    model.load_state_dict(state)
    model.to(device).eval()
    tokenizer = get_tokenizer(meta["model_id"])
    return model, tokenizer, meta, device


def _items(sequences: Sequence[str], loops: Optional[Sequence[Optional[str]]], targets):
    n_t = len(targets)
    loops = loops if loops is not None else [None] * len(sequences)
    items, missing = [], 0
    for seq, loop in zip(sequences, loops):
        start = locate_loop(seq, loop) if loop else -1
        if loop and start < 0:
            missing += 1
        items.append({
            "sequence": str(seq),
            "labels": [-100] * n_t,
            "count_weight": [0.0] * n_t,
            "loop_start": start,
            "loop_len": len(loop) if isinstance(loop, str) and loop else 0,
        })
    return items, missing


@torch.no_grad()
def predict_sequences(model, tokenizer, meta, sequences: List[str],
                      loops: Optional[List[Optional[str]]] = None,
                      device=None, batch_size: int = 16) -> pd.DataFrame:
    targets = meta["targets"]
    device = device or next(model.parameters()).device
    if meta["pooling"] == "loop" and (loops is None or all(l in (None, "") for l in (loops or []))):
        print("[warn] model uses loop pooling but no loop motif was provided; "
              "falling back to mean pooling for these sequences (degraded accuracy).")
    items, missing = _items(sequences, loops, targets)
    if missing:
        print(f"[warn] {missing} loop motif(s) not found in their sequence; "
              "those fall back to mean pooling.")
    collate = Collator(tokenizer, bos_offset=meta["bos_offset"], max_length=meta["max_length"])
    loader = DataLoader(items, batch_size=batch_size, shuffle=False, collate_fn=collate)

    use_amp = device.type == "cuda"
    logits_all = []
    for batch in loader:
        batch = {k: (v.to(device) if torch.is_tensor(v) else v) for k, v in batch.items()}
        with torch.amp.autocast("cuda", enabled=use_amp):
            out = model(input_ids=batch["input_ids"], attention_mask=batch["attention_mask"],
                        loop_start=batch["loop_start"], loop_len=batch["loop_len"])
        logits_all.append(out["logits"].float().cpu())
    logits = torch.cat(logits_all).numpy()
    probs = 1.0 / (1.0 + np.exp(-logits))

    thresholds = meta.get("thresholds", {})
    df = pd.DataFrame({"sequence": list(sequences)})
    if loops is not None:
        df["loop"] = list(loops)
    for ti, t in enumerate(targets):
        thr = float(thresholds.get(t, 0.5))
        df[f"prob_{t}"] = probs[:, ti]
        df[f"call_{t}"] = (probs[:, ti] >= thr).astype(int)
    # convenience summary: any-target binder + the strongest target
    prob_cols = [f"prob_{t}" for t in targets]
    df["max_prob"] = df[prob_cols].max(axis=1)
    df["top_target"] = df[prob_cols].idxmax(axis=1).str.replace("prob_", "", regex=False)
    return df
