"""Multi-task ESM-C model + dataset + collator.

Architecture
------------
ESM-C backbone (Synthyra ESM++)  ->  pooled representation  ->  shared MLP  ->
``Linear(hidden, n_targets)`` producing ONE binary logit per target.

Labels use a -100 sentinel for "not assayed against this target"; the loss
(masked BCE-with-logits, per-target ``pos_weight``) only counts the targets a
sequence was actually screened against. At inference, ``sigmoid(logits)`` gives
an independent binding probability for every target.

Pooling
-------
``loop`` (default) averages only the grafted-loop token embeddings — the 6-aa
loop is the entire discriminative signal, so pooling the whole 188-aa sequence
would dilute it ~30x. ``mean`` / ``cls`` are kept as ablations.
"""
from __future__ import annotations

from contextlib import nullcontext
from typing import List, Optional

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset


# --------------------------------------------------------------------------- #
# Model
# --------------------------------------------------------------------------- #
class MultiTaskESMC(nn.Module):
    def __init__(self, model_id: str, targets: List[str], pooling: str = "loop",
                 dropout: float = 0.1, pos_weight: Optional[torch.Tensor] = None,
                 trust_remote_code: bool = True):
        super().__init__()
        from transformers import AutoModel

        self.backbone = AutoModel.from_pretrained(
            model_id, trust_remote_code=trust_remote_code
        )
        self.targets = list(targets)
        self.num_targets = len(targets)
        self.pooling = pooling
        hidden = self._infer_hidden_size()

        self.head = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(hidden, hidden // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden // 2, self.num_targets),
        )
        # Always register the buffer (default = no weighting) so a model saved
        # with pos_weight loads cleanly when rebuilt for inference (pos_weight=None).
        if pos_weight is None:
            pos_weight = torch.ones(self.num_targets)
        self.register_buffer("pos_weight", pos_weight.float())
        # Set by the trainer so autocast runs INSIDE forward (correct under
        # nn.DataParallel, whose replica threads don't inherit an outer autocast).
        self._amp_dtype = None

    def enable_gradient_checkpointing(self):
        fn = getattr(self.backbone, "gradient_checkpointing_enable", None)
        if callable(fn):
            fn()
            return True
        return False

    def _infer_hidden_size(self) -> int:
        cfg = getattr(self.backbone, "config", None)
        for attr in ("hidden_size", "embed_dim", "d_model"):
            if cfg is not None and getattr(cfg, attr, None):
                return int(getattr(cfg, attr))
        # last resort: probe a forward pass
        with torch.no_grad():
            ids = torch.tensor([[0, 1, 2, 3]])
            hs = self._hidden_states(self.backbone(input_ids=ids))
            return int(hs.shape[-1])

    @staticmethod
    def _hidden_states(out) -> torch.Tensor:
        if hasattr(out, "last_hidden_state"):
            return out.last_hidden_state
        if isinstance(out, (tuple, list)):
            return out[0]
        if hasattr(out, "hidden_states") and out.hidden_states is not None:
            return out.hidden_states[-1]
        raise RuntimeError("Could not extract hidden states from backbone output.")

    def _pool(self, hs, attention_mask, loop_start, loop_len):
        if self.pooling == "cls":
            return hs[:, 0]
        if self.pooling == "mean":
            m = attention_mask.unsqueeze(-1).to(hs.dtype)
            return (hs * m).sum(1) / m.sum(1).clamp(min=1.0)
        # loop pooling (per-example slice; batch sizes are small)
        pooled = []
        for i in range(hs.size(0)):
            s = int(loop_start[i]) if loop_start is not None else -1
            ln = int(loop_len[i]) if loop_len is not None else 0
            if s >= 0 and ln > 0 and s + ln <= hs.size(1):
                pooled.append(hs[i, s:s + ln].mean(0))
            else:  # fall back to masked mean for this example
                m = attention_mask[i].unsqueeze(-1).to(hs.dtype)
                pooled.append((hs[i] * m).sum(0) / m.sum(0).clamp(min=1.0))
        return torch.stack(pooled)

    def forward(self, input_ids, attention_mask=None, labels=None,
                count_weight=None, loop_start=None, loop_len=None):
        amp = (torch.autocast("cuda", dtype=self._amp_dtype)
               if self._amp_dtype is not None else nullcontext())
        with amp:
            out = self.backbone(input_ids=input_ids, attention_mask=attention_mask)
            hs = self._hidden_states(out)
            pooled = self._pool(hs, attention_mask, loop_start, loop_len)
            logits = self.head(pooled)

            loss = None
            if labels is not None:
                target_mask = (labels != -100).float()
                safe_labels = torch.where(labels == -100,
                                          torch.zeros_like(labels), labels).float()
                bce = F.binary_cross_entropy_with_logits(
                    logits, safe_labels, reduction="none",
                    pos_weight=self.pos_weight,
                )
                bce = bce * target_mask
                if count_weight is not None:
                    bce = bce * count_weight
                denom = target_mask.sum().clamp(min=1.0)
                loss = bce.sum() / denom
        return {"loss": loss, "logits": logits}

    @torch.no_grad()
    def embed(self, input_ids, attention_mask=None, loop_start=None, loop_len=None):
        """Return the pooled representation the classifier head sees (for viz)."""
        out = self.backbone(input_ids=input_ids, attention_mask=attention_mask)
        hs = self._hidden_states(out)
        return self._pool(hs, attention_mask, loop_start, loop_len)

    # -- (un)freezing helpers --------------------------------------------- #
    def freeze_backbone(self):
        for p in self.backbone.parameters():
            p.requires_grad = False

    def set_trainable_top_layers(self, n_last):
        """Unfreeze the top-``n_last`` transformer blocks (+ trailing norm).

        ``n_last`` may be an int, 0 (head only), or 'all'/-1 (whole backbone).
        """
        import re

        if n_last in ("all", -1, "-1"):
            for p in self.backbone.parameters():
                p.requires_grad = True
            return self._count_trainable()

        self.freeze_backbone()
        n_last = int(n_last)
        if n_last <= 0:
            return self._count_trainable()

        # Discover block indices from parameter names (e.g. ...layers.7... ).
        layer_re = re.compile(r"\.(\d+)\.")
        idxs = []
        for name, _ in self.backbone.named_parameters():
            m = layer_re.search(name)
            if m:
                idxs.append(int(m.group(1)))
        if not idxs:
            # no detectable blocks -> just unfreeze everything as a safe default
            for p in self.backbone.parameters():
                p.requires_grad = True
            return self._count_trainable()
        max_idx = max(idxs)
        threshold = max_idx - n_last + 1
        for name, p in self.backbone.named_parameters():
            m = layer_re.search(name)
            if m and int(m.group(1)) >= threshold:
                p.requires_grad = True
            elif m is None and ("norm" in name.lower() and "embed" not in name.lower()):
                p.requires_grad = True  # trailing final-norm
        return self._count_trainable()

    def _count_trainable(self):
        t = sum(p.numel() for p in self.parameters() if p.requires_grad)
        total = sum(p.numel() for p in self.parameters())
        return t, total


# --------------------------------------------------------------------------- #
# Dataset
# --------------------------------------------------------------------------- #
class SeqDataset(Dataset):
    """Reads a prepared parquet split into per-sequence training examples."""

    def __init__(self, source, targets: List[str], count_weighting: str = "log"):
        # ``source`` may be a parquet path or an already-loaded DataFrame
        # (the latter lets the trainer score arbitrary subsets, e.g. novel loops).
        self.df = (source.reset_index(drop=True) if isinstance(source, pd.DataFrame)
                   else pd.read_parquet(source).reset_index(drop=True))
        self.targets = list(targets)
        self.count_weighting = count_weighting

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        r = self.df.iloc[idx]
        labels, weights = [], []
        for t in self.targets:
            mask = int(r[f"mask_{t}"])
            lab = int(r[f"label_{t}"]) if mask else -100
            labels.append(lab)
            if mask:
                c = float(r[f"count_{t}"])
                w = np.log1p(c) if self.count_weighting == "log" else 1.0
            else:
                w = 0.0
            weights.append(float(w))
        return {
            "sequence": str(r["sequence"]),
            "labels": labels,
            "count_weight": weights,
            "loop_start": int(r["loop_start"]),
            "loop_len": int(r["loop_len"]),
        }


# --------------------------------------------------------------------------- #
# Collator
# --------------------------------------------------------------------------- #
class Collator:
    """Dynamic-padding collator that tokenises a batch of sequences.

    Token loop offset = char ``loop_start`` + ``bos_offset`` (residue-level
    tokenisation maps char i -> token i + bos_offset).
    """

    def __init__(self, tokenizer, bos_offset: int = 1, max_length: int = 1024):
        self.tokenizer = tokenizer
        self.bos_offset = bos_offset
        self.max_length = max_length

    def __call__(self, batch):
        seqs = [b["sequence"] for b in batch]
        enc = self.tokenizer(
            seqs, padding=True, truncation=True, max_length=self.max_length,
            return_tensors="pt",
        )
        labels = torch.tensor([b["labels"] for b in batch], dtype=torch.float)
        count_weight = torch.tensor([b["count_weight"] for b in batch], dtype=torch.float)
        loop_start = torch.tensor(
            [b["loop_start"] + self.bos_offset if b["loop_start"] >= 0 else -1
             for b in batch], dtype=torch.long)
        loop_len = torch.tensor([b["loop_len"] for b in batch], dtype=torch.long)
        return {
            "input_ids": enc["input_ids"],
            "attention_mask": enc.get("attention_mask"),
            "labels": labels,
            "count_weight": count_weight,
            "loop_start": loop_start,
            "loop_len": loop_len,
        }
