"""Score the RFd3-pipeline designs with every trained ESM-C variant (CPU-safe, resumable).

One backbone pass per sequence, then the classifier head is applied to the loop-pooled
embedding for EACH window (AB, C, EF). The AB/C windows are the two the models were trained
on; EF is a null-window control (none of the models ever saw EF-loop variation).

The window is pooled exactly as in training (mean of the loop-token embeddings), except that
the design loops are 2-15 aa rather than the fixed 6 aa the models saw.

Outputs: Local/esmc_rfd3_agreement/scores/<variant>.csv with one row per (seq_id, window).
Held-out lab sequences (controls.csv) are scored only by the variant they came from.

    python score_designs.py --variants all --priority 1
    python score_designs.py --variants cloop_only --parity     # check vs. predict_sequences
"""
from __future__ import annotations

import argparse
import ctypes
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("HF_HUB_OFFLINE", "1")   # weights are in the local HF cache
HERE = Path(__file__).resolve().parent
ESMC_DIR = HERE.parent
sys.path.insert(0, str(ESMC_DIR))

import numpy as np
import pandas as pd
import torch

from inference import load_trained_model, predict_sequences

LOCAL = HERE.parents[1] / "Local"
RUN = LOCAL / "esmc_rfd3_agreement"
MODELS = LOCAL / "esmc_multirun"
VARIANTS = ["abloop_only", "cloop_only", "mmp9_other", "all3_original", "everything_combined"]
WINDOWS = ["AB", "C", "EF"]
ALL_TARGETS = ["ADAM17", "MMP3", "MMP9"]


def lower_priority():
    """Yield the CPU to the user's GPU sweep / interactive work (Windows only)."""
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x4000)
    except Exception:
        pass


@torch.no_grad()
def score_windows(model, tok, meta, items, batch_size):
    """items: list of dicts {seq_id, full_seq, start_<W>, len_<W>}. Returns long DataFrame."""
    bos = meta["bos_offset"]
    order = sorted(range(len(items)), key=lambda i: len(items[i]["full_seq"]))
    rows = []
    for b in range(0, len(order), batch_size):
        idx = order[b:b + batch_size]
        batch = [items[i] for i in idx]
        dev = next(model.parameters()).device
        enc = tok([x["full_seq"] for x in batch], padding=True, return_tensors="pt")
        out = model.backbone(input_ids=enc["input_ids"].to(dev), attention_mask=enc["attention_mask"].to(dev))
        hs = model._hidden_states(out).float()
        for w in WINDOWS:
            pooled, ok = [], []
            for j, x in enumerate(batch):
                s, n = int(x[f"start_{w}"]), int(x[f"len_{w}"])
                if s >= 0 and n > 0:
                    pooled.append(hs[j, s + bos:s + bos + n].mean(0)); ok.append(True)
                else:
                    pooled.append(torch.zeros(hs.size(-1), device=dev)); ok.append(False)
            probs = torch.sigmoid(model.head(torch.stack(pooled)).float()).cpu().numpy()
            for j, x in enumerate(batch):
                rec = {"seq_id": x["seq_id"], "window": w}
                for ti, t in enumerate(meta["targets"]):
                    rec[f"prob_{t}"] = float(probs[j, ti]) if ok[j] else np.nan
                rows.append(rec)
    return pd.DataFrame(rows)


def parity_check(model, tok, meta, seqs_df, device):
    """The window scorer must reproduce the repo's own predict_sequences() on the AB/C windows."""
    sub = seqs_df[(seqs_df.len_AB > 0) & (seqs_df.len_C > 0)].head(8)
    items = sub.to_dict("records")
    mine = score_windows(model, tok, meta, items, batch_size=4)
    worst = 0.0
    for w in ("AB", "C"):
        loops = [s[st:st + n] for s, st, n in zip(sub.full_seq, sub[f"start_{w}"], sub[f"len_{w}"])]
        ref = predict_sequences(model, tok, meta, sub.full_seq.tolist(), loops=loops, device=device, batch_size=4)
        m = mine[mine.window == w].set_index("seq_id").loc[sub.seq_id]
        for t in meta["targets"]:
            worst = max(worst, float(np.abs(m[f"prob_{t}"].values - ref[f"prob_{t}"].values).max()))
    print(f"  parity vs predict_sequences: max |dprob| = {worst:.2e}")
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", default="all")
    ap.add_argument("--priority", type=int, default=2, help="score sequences with priority <= this")
    ap.add_argument("--no-controls", action="store_true")
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--device", default="auto", help="auto | cuda | cpu")
    ap.add_argument("--limit", type=int, default=0, help="score only the first N sequences (benchmark)")
    ap.add_argument("--parity", action="store_true", help="only run the parity check, then exit")
    args = ap.parse_args()

    lower_priority()
    torch.set_num_threads(args.threads)
    variants = VARIANTS if args.variants == "all" else args.variants.split(",")

    seqs = pd.read_csv(RUN / "seqs.csv")
    seqs = seqs[seqs.priority <= args.priority]
    if args.limit:
        seqs = seqs.head(args.limit)
    ctrl = pd.read_csv(RUN / "controls.csv")
    ctrl["seq_id"] = [f"C{i:05d}" for i in range(len(ctrl))]
    (RUN / "scores").mkdir(exist_ok=True)
    ctrl.to_csv(RUN / "controls_ids.csv", index=False)

    for v in variants:
        out_path = RUN / "scores" / f"{v}.csv"
        t0 = time.time()
        dev = torch.device(args.device if args.device != "auto" else ("cuda" if torch.cuda.is_available() else "cpu"))
        model, tok, meta, device = load_trained_model(MODELS / v / "model", device=dev)
        print(f"[{v}] loaded in {time.time() - t0:.0f}s  targets={meta['targets']}", flush=True)
        if args.parity:
            parity_check(model, tok, meta, seqs, device)
            continue

        items = seqs.to_dict("records")
        if not args.no_controls:
            c = ctrl[ctrl.variant == v]
            from build_design_set import locate_loops
            for r in c.itertuples():
                L = locate_loops(r.full_seq)
                items.append({"seq_id": r.seq_id, "full_seq": r.full_seq,
                              **{f"start_{w}": L[w][0] for w in WINDOWS},
                              **{f"len_{w}": L[w][1] for w in WINDOWS}})
        done = set()
        if out_path.exists():
            done = set(pd.read_csv(out_path, usecols=["seq_id"]).seq_id)
        todo = [x for x in items if x["seq_id"] not in done]
        print(f"[{v}] {len(todo)} to score ({len(done)} already done)", flush=True)

        CH = 48
        t1 = time.time()
        for k in range(0, len(todo), CH):
            df = score_windows(model, tok, meta, todo[k:k + CH], args.batch_size)
            df.to_csv(out_path, mode="a", header=not out_path.exists(), index=False)
            n = min(k + CH, len(todo))
            rate = n / (time.time() - t1)
            print(f"[{v}] {n}/{len(todo)}  {rate:.2f} seq/s  eta {(len(todo) - n) / max(rate, 1e-9) / 60:.1f} min",
                  flush=True)
        del model


if __name__ == "__main__":
    main()
