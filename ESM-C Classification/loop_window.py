"""Locate the fixed-position grafted-loop window for a given Ligand_Subtype.

Shared by enumerate_cloop.py and shap_hotspots.py: both need to know, for a
given training CSV and loop subtype (e.g. "C-loop" or "AB-loop"), which
character window of the constant 188-aa scaffold the 6-aa loop occupies, so
they can build synthetic full-length sequences by grafting arbitrary loops
into that window.
"""
from __future__ import annotations

from collections import Counter
from typing import Callable, List

import numpy as np
import pandas as pd
import torch

from esmc_utils import resolve_path

AAS = "ACDEFGHIKLMNPQRSTVWY"


def derive_loop_template(cfg: dict, loop_subtype: str | None = None, csv_path=None):
    """Return (masked_template, start, length) for the given loop subtype.

    ``masked_template`` has the loop window replaced with underscores.
    ``loop_subtype`` filters the ``Ligand_Subtype`` column if present and a
    value is given; pass ``None`` to use the whole (unfiltered) CSV, matching
    the original C-loop-only behaviour when the column is absent.
    """
    d = cfg["data"]
    path = csv_path if csv_path is not None else resolve_path(cfg, d["csv_path"])
    df = pd.read_csv(path)
    sub_col = "Ligand_Subtype"
    if loop_subtype and sub_col in df.columns:
        c = df[df[sub_col] == loop_subtype]
        if len(c) == 0:
            raise SystemExit(f"No rows with {sub_col} == {loop_subtype!r} in {path}")
    else:
        c = df
    seqs, loops = c[d["seq_col"]].astype(str), c[d["loop_col"]].astype(str)
    starts = [s.find(l) for s, l in zip(seqs, loops)]
    valid_starts = [x for x in starts if x >= 0]
    if not valid_starts:
        raise SystemExit(f"Could not locate any loop motif in {path} (subtype={loop_subtype})")
    start = int(pd.Series(valid_starts).mode()[0])
    L = len(loops.iloc[0])
    masked = [s[:start] + "_" * L + s[start + L:] for s, st in zip(seqs, starts) if st == start]
    uniq = Counter(masked)
    print(f"distinct scaffolds outside the window (subtype={loop_subtype}): {len(uniq)} (want 1)")
    return uniq.most_common(1)[0][0], start, L


def idx_to_loop(idx: int, L: int) -> str:
    out = []
    for _ in range(L):
        out.append(AAS[idx % 20])
        idx //= 20
    return "".join(out)


def build_scorer(model, tokenizer, meta: dict, templ_masked: str, start_char: int, L: int,
                  batch_size: int = 512) -> Callable[[List[str]], np.ndarray]:
    """Return ``score(loops) -> (n, n_targets) prob array`` for a fixed scaffold window.

    Grafts each loop string into ``templ_masked`` by direct token substitution
    (no re-tokenizing the full sequence per call), batched at ``batch_size``.
    Shared by enumerate_cloop.py (exhaustive sweep) and shap_hotspots.py
    (SHAP background/foreground scoring) so the graft mechanics only live once.
    """
    device = next(model.parameters()).device
    templ_seq = templ_masked.replace("_" * L, "A" * L)
    templ_ids = tokenizer(templ_seq, return_tensors="pt")["input_ids"][0].to(device)
    Ltok = templ_ids.shape[0]
    tok_lo = start_char + meta["bos_offset"]
    aa_id = {a: int(tokenizer.convert_tokens_to_ids(a)) for a in AAS}
    cols = list(range(tok_lo, tok_lo + L))
    attn = torch.ones(Ltok, dtype=torch.long, device=device)

    @torch.no_grad()
    def score(loops: List[str]) -> np.ndarray:
        out_all = []
        for i in range(0, len(loops), batch_size):
            chunk = loops[i:i + batch_size]
            n = len(chunk)
            ids = templ_ids.unsqueeze(0).repeat(n, 1)
            for c, p in enumerate(cols):
                ids[:, p] = torch.tensor([aa_id[chunk[j][c]] for j in range(n)], device=device)
            ls = torch.full((n,), tok_lo, dtype=torch.long, device=device)
            ll = torch.full((n,), L, dtype=torch.long, device=device)
            o = model(input_ids=ids, attention_mask=attn.unsqueeze(0).expand(n, -1),
                      loop_start=ls, loop_len=ll)
            out_all.append(torch.sigmoid(o["logits"].float()).cpu().numpy())
        return np.concatenate(out_all, axis=0)

    score.tok_lo = tok_lo
    score.Ltok = Ltok
    return score
