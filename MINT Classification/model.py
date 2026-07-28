"""MINT-backed binder classifier + pair dataset/collator.

Architecture
------------
mint.helpers.extract.MINTWrapper (pooled pair embedding, ESM2-650M backbone
pretrained on 96M protein-protein interactions) -> small MLP head -> ONE
bind/no-bind logit.

Unlike the ESM-C pipeline's per-target output heads, there is a single shared
head here: the target protease is chain B of the input pair, not an output
index, so one head covers every target uniformly (including ones with no
head-specific data, which was the whole AB-loop/C-loop headache in the ESM-C
multirun).

Requires the `mint` package (https://github.com/VarunUllanat/mint) installed
from a cloned checkout (`pip install -e .`) plus its checkpoint + config JSON
-- see README.md "Install". The import is deferred into `_load_mint` so that
scripts which don't need the actual backbone (nothing in this pipeline, but
kept for import-safety) still load.

NOTE: `_backbone()` assumes MINTWrapper exposes the underlying ESM2 model as
``self.model`` (inferred from `mint/helpers/extract.py`'s forward method
referencing `self.model.cls_idx` / `self.model.eos_idx`). This has not been
verified against an actual installed `mint` package (no GPU/checkpoint here)
-- if `_backbone()` raises AttributeError on the server, check the installed
`mint.helpers.extract.MINTWrapper.__init__` for the real attribute name and
fix the one line in `_backbone()`.
"""
from __future__ import annotations

import re
from typing import Optional

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset


def _load_mint(cfg_json: str, checkpoint_path: str, freeze_percent: float,
               use_multimer: bool, sep_chains: bool, device):
    from mint.helpers.extract import CollateFn, MINTWrapper
    from mint.helpers.extract import load_config as mint_load_config

    mcfg = mint_load_config(cfg_json)
    wrapper = MINTWrapper(mcfg, checkpoint_path, freeze_percent=freeze_percent,
                          use_multimer=use_multimer, sep_chains=sep_chains, device=str(device))
    return wrapper, CollateFn


class MintBinderClassifier(nn.Module):
    def __init__(self, mint_cfg_json: str, checkpoint_path: str, device,
                sep_chains: bool = True, use_multimer: bool = True,
                truncation_seq_length: int = 512,
                head_hidden: int = 256, dropout: float = 0.1,
                pos_weight: Optional[torch.Tensor] = None,
                freeze_percent: float = 1.0):
        super().__init__()
        self.wrapper, collate_cls = _load_mint(
            mint_cfg_json, checkpoint_path, freeze_percent, use_multimer, sep_chains, device)
        self.mint_collate = collate_cls(truncation_seq_length)
        emb_dim = 2560 if sep_chains else 1280
        self.head = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(emb_dim, head_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(head_hidden, 1),
        )
        if pos_weight is None:
            pos_weight = torch.ones(1)
        self.register_buffer("pos_weight", pos_weight.float().reshape(1))

    # -- backbone (un)freezing ---------------------------------------------
    # MINTWrapper only exposes a freeze_percent set at construction; we
    # instead manage requires_grad directly (same top-N-layers scheme as the
    # ESM-C pipeline's model.py) so phase1/phase2 don't need to rebuild the
    # wrapper (and re-load the checkpoint) between phases.
    def _backbone(self) -> nn.Module:
        m = getattr(self.wrapper, "model", None)
        if m is None:
            raise AttributeError(
                "MINTWrapper has no '.model' attribute in this version of the mint package -- "
                "see the module docstring; update _backbone() to match.")
        return m

    def freeze_backbone(self):
        for p in self._backbone().parameters():
            p.requires_grad = False

    def set_trainable_top_layers(self, n_last):
        """Unfreeze the top-``n_last`` transformer blocks (+ trailing norm).
        ``n_last`` may be an int, 0 (head only), or 'all'/-1 (whole backbone)."""
        if n_last in ("all", -1, "-1"):
            for p in self._backbone().parameters():
                p.requires_grad = True
            return self._count_trainable()

        self.freeze_backbone()
        n_last = int(n_last)
        if n_last <= 0:
            return self._count_trainable()

        layer_re = re.compile(r"\.(\d+)\.")
        idxs = []
        for name, _ in self._backbone().named_parameters():
            m = layer_re.search(name)
            if m:
                idxs.append(int(m.group(1)))
        if not idxs:
            for p in self._backbone().parameters():
                p.requires_grad = True
            return self._count_trainable()
        threshold = max(idxs) - n_last + 1
        for name, p in self._backbone().named_parameters():
            m = layer_re.search(name)
            if m and int(m.group(1)) >= threshold:
                p.requires_grad = True
            elif m is None and ("norm" in name.lower() and "embed" not in name.lower()):
                p.requires_grad = True  # trailing final-norm
        return self._count_trainable()

    def count_backbone_layers(self) -> int:
        layer_re = re.compile(r"\.(\d+)\.")
        idxs = [int(m.group(1)) for name, _ in self._backbone().named_parameters()
                if layer_re.search(name)]
        return (max(idxs) + 1) if idxs else 0

    def set_trainable_by_freeze_percent(self, freeze_percent: float):
        """``freeze_percent`` in [0,1]: fraction of backbone layers kept
        frozen, counted from the input side -- matches MINTWrapper's own
        ``freeze_percent`` convention (e.g. 0.5 unfreezes the final 50% of
        transformer blocks). Internally reuses ``set_trainable_top_layers``."""
        total = self.count_backbone_layers()
        if total == 0:
            return self.set_trainable_top_layers("all" if freeze_percent < 1.0 else 0)
        n_last = round((1.0 - float(freeze_percent)) * total)
        return self.set_trainable_top_layers(n_last)

    def _count_trainable(self):
        t = sum(p.numel() for p in self.parameters() if p.requires_grad)
        total = sum(p.numel() for p in self.parameters())
        return t, total

    def enable_gradient_checkpointing(self):
        fn = getattr(self._backbone(), "gradient_checkpointing_enable", None)
        if callable(fn):
            fn()
            return True
        return False

    # -- forward -------------------------------------------------------------
    def forward(self, chains, chain_ids, labels=None, count_weight=None):
        emb = self.wrapper(chains, chain_ids)     # (batch, emb_dim)
        logits = self.head(emb).squeeze(-1)        # (batch,)

        loss = None
        if labels is not None:
            bce = F.binary_cross_entropy_with_logits(
                logits, labels.float(), reduction="none", pos_weight=self.pos_weight)
            if count_weight is not None:
                loss = (bce * count_weight).sum() / count_weight.sum().clamp(min=1e-8)
            else:
                loss = bce.mean()
        return {"loss": loss, "logits": logits}

    @torch.no_grad()
    def embed(self, chains, chain_ids):
        """Pooled pair embedding the head sees (for downstream viz/SHAP)."""
        return self.wrapper(chains, chain_ids)


# --------------------------------------------------------------------------- #
# Dataset + collator
# --------------------------------------------------------------------------- #
class PairDataset(Dataset):
    """Reads a prepared parquet split (from data_prep.py) into
    (seqA, seqB, label, weight) training examples."""

    def __init__(self, source, count_weighting: str = "log"):
        self.df = (source.reset_index(drop=True) if isinstance(source, pd.DataFrame)
                   else pd.read_parquet(source).reset_index(drop=True))
        self.count_weighting = count_weighting

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        r = self.df.iloc[idx]
        c = float(r["count"])
        w = np.log1p(c) if self.count_weighting == "log" else 1.0
        return {"seqA": str(r["seqA"]), "seqB": str(r["seqB"]),
                "label": float(r["label"]), "weight": float(w),
                "target_name": str(r["target_name"])}


class PairCollator:
    """Wraps mint's own CollateFn (tokenizes the (seqA,seqB) pairs into
    chains/chain_ids per mint's Alphabet) and additionally stacks
    labels/weights in the same batch order (mint's CollateFn only returns
    chains/chain_ids -- it has no notion of a training label)."""

    def __init__(self, mint_collate_fn):
        self._mint_collate = mint_collate_fn  # instance of mint.helpers.extract.CollateFn

    def __call__(self, batch):
        pairs = [(b["seqA"], b["seqB"]) for b in batch]
        chains, chain_ids = self._mint_collate(pairs)
        labels = torch.tensor([b["label"] for b in batch], dtype=torch.float)
        weights = torch.tensor([b["weight"] for b in batch], dtype=torch.float)
        target_names = [b["target_name"] for b in batch]  # not a tensor -- per-target reporting only
        return {"chains": chains, "chain_ids": chain_ids, "labels": labels,
                "count_weight": weights, "target_name": target_names}
