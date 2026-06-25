"""Staged multi-task fine-tuning of ESM-C for per-target binder classification.

Curriculum (this is the "staged" training, reframed as a robust optimisation
schedule rather than a fragile data hierarchy):

  Phase 1  - backbone frozen, train the randomly-initialised heads (warm-up).
  Phase 2  - unfreeze the top-N transformer blocks, low LR, fine-tune with
             early stopping on the validation set.

A SINGLE held-out test set is untouched through both phases and scored once at
the end. Per-target decision thresholds are chosen on validation (max-MCC) and
saved so downstream prediction isn't a naive 0.5 cut.

Usage
-----
    python train.py --config config.yaml            # full run
    python train.py --config config.yaml --smoke     # tiny, fast sanity run
"""
from __future__ import annotations

import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader
from tqdm.auto import tqdm

from esmc_utils import (
    best_threshold,
    compute_target_metrics,
    detect_bos_offset,
    get_device,
    get_tokenizer,
    load_config,
    mean_metric,
    resolve_path,
    set_seed,
)
from model import Collator, MultiTaskESMC, SeqDataset


def to_device(batch, device):
    return {k: (v.to(device) if torch.is_tensor(v) else v) for k, v in batch.items()}


def resolve_precision(precision, device):
    """Return (amp_dtype | None, GradScaler | None) for the requested precision.

    'auto' picks bf16 on capable GPUs (more stable, no loss scaler needed) and
    falls back to fp16 otherwise. bf16/fp32 need no scaler; fp16 does.
    Autocast itself happens inside ``MultiTaskESMC.forward`` (DataParallel-safe).
    """
    if device.type != "cuda":
        return None, None
    if precision == "auto":
        # bf16 needs Ampere+ (CC >= 8.0); on Volta/Turing (e.g. V100 = 7.0) bf16
        # has no hardware support and is emulated/slow, so use fp16 there.
        cc_major = torch.cuda.get_device_capability(0)[0]
        precision = "bf16" if cc_major >= 8 else "fp16"
    if precision == "bf16":
        return torch.bfloat16, None
    if precision == "fp16":
        return torch.float16, torch.amp.GradScaler("cuda")
    return None, None  # fp32


def run_train_epoch(model, loader, optimizer, scaler, device, grad_clip, desc, log_every=0):
    model.train()
    total, n = 0.0, 0
    params = [p for p in model.parameters() if p.requires_grad]
    n_steps = len(loader)
    last_t, last_i = time.time(), 0
    # tqdm bar on a terminal; auto-hidden when output is redirected (clean cluster logs).
    for i, batch in enumerate(tqdm(loader, desc=desc, leave=False, disable=None)):
        batch = to_device(batch, device)
        optimizer.zero_grad()
        out = model(**batch)              # autocast is inside forward
        loss = out["loss"]
        if loss.dim() > 0:                # DataParallel gathers one loss per GPU
            loss = loss.mean()
        if scaler is not None:
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            torch.nn.utils.clip_grad_norm_(params, grad_clip)
            scaler.step(optimizer)
            scaler.update()
        else:
            loss.backward()
            torch.nn.utils.clip_grad_norm_(params, grad_clip)
            optimizer.step()
        bs = batch["input_ids"].size(0)
        total += loss.item() * bs
        n += bs
        if log_every and ((i + 1) % log_every == 0 or i + 1 == n_steps):
            now = time.time()
            rate = (i + 1 - last_i) / max(now - last_t, 1e-9)  # window rate, not cumulative
            last_t, last_i = now, i + 1
            print(f"  {desc}: step {i+1}/{n_steps}  loss={loss.item():.4f}  "
                  f"{rate:.1f} it/s", flush=True)
    return total / max(n, 1)


@torch.no_grad()
def evaluate(model, loader, device, targets, thresholds=None):
    model.eval()
    logits_all, labels_all = [], []
    for batch in loader:
        b = to_device(batch, device)
        out = model(**b)                  # autocast is inside forward
        logits_all.append(out["logits"].float().cpu())
        labels_all.append(batch["labels"])
    logits = torch.cat(logits_all).numpy()
    labels = torch.cat(labels_all).numpy()
    probs = 1.0 / (1.0 + np.exp(-logits))
    per = {}
    for ti, t in enumerate(targets):
        m = labels[:, ti] != -100
        if m.sum() == 0:
            continue
        thr = 0.5 if thresholds is None else thresholds.get(t, 0.5)
        per[t] = compute_target_metrics(labels[m, ti], probs[m, ti], threshold=thr)
    return per, probs, labels


def build_optimizer(model, lr, weight_decay):
    params = [p for p in model.parameters() if p.requires_grad]
    return torch.optim.AdamW(params, lr=lr, weight_decay=weight_decay)


def fmt_report(per, targets):
    lines = [f"{'target':<8} {'n':>6} {'pos':>6} {'pos_rate':>8} "
             f"{'pr_auc':>7} {'roc_auc':>7} {'mcc':>7} {'f1':>7}"]
    for t in targets:
        if t not in per:
            continue
        m = per[t]
        def g(k):
            v = m[k]
            return f"{v:.3f}" if v == v else "  nan"
        lines.append(f"{t:<8} {m['n']:>6} {m['n_pos']:>6} "
                     f"{(m['pos_rate'] if m['pos_rate']==m['pos_rate'] else 0):>8.3f} "
                     f"{g('pr_auc'):>7} {g('roc_auc'):>7} {g('mcc'):>7} {g('f1'):>7}")
    return "\n".join(lines)


def train_and_evaluate(cfg, smoke=False, pooling=None, out_name=None,
                       model_id=None, strict_test=False, log_every=None):
    """Run the full staged fine-tune + held-out test evaluation once.

    Returns a result dict (per-target test metrics, thresholds, best val score,
    output dir). ``pooling`` overrides ``cfg['model']['pooling']`` (used by the
    pooling-comparison harness); ``model_id`` overrides the backbone;
    ``out_name`` overrides the model output folder; ``strict_test`` adds a
    novel-loop (no train near-neighbour) test report.
    """
    set_seed(int(cfg["split"].get("seed", 42)))
    targets = cfg["targets"]
    device = get_device()
    pooling = pooling or cfg["model"]["pooling"]
    model_id = model_id or cfg["model"]["model_id"]
    print(f"Device: {device} | model: {model_id} | pooling: {pooling}")
    if device.type == "cpu":
        print("\n" + "!" * 70 + "\n"
              "WARNING: no GPU detected -- training the 600M model on CPU is\n"
              "impractically slow. This usually means the installed torch CUDA build\n"
              "is newer than the GPU driver supports (see 'NVIDIA driver too old').\n"
              "Fix: pip install torch --index-url https://download.pytorch.org/whl/cu128\n"
              "(match the cuXYZ build to your driver's CUDA in `nvidia-smi`).\n" + "!" * 70 + "\n")

    # --- resolve hyper-params (apply smoke overrides) -----------------------
    args = type("A", (), {"smoke": smoke})()  # shim so the body below reads args.smoke
    tr = cfg["train"]
    smoke = cfg.get("smoke", {})
    batch_size = smoke.get("batch_size", tr["batch_size"]) if args.smoke else tr["batch_size"]
    p1_epochs = smoke.get("phase1_epochs", tr["phase1"]["epochs"]) if args.smoke else tr["phase1"]["epochs"]
    p2_epochs = smoke.get("phase2_epochs", tr["phase2"]["epochs"]) if args.smoke else tr["phase2"]["epochs"]
    unfreeze = smoke.get("unfreeze_layers", tr["phase2"]["unfreeze_layers"]) if args.smoke else tr["phase2"]["unfreeze_layers"]
    grad_clip = float(tr.get("grad_clip", 1.0))
    sel_metric = tr.get("metric", "pr_auc")
    log_every = int(tr.get("log_every", 0)) if log_every is None else int(log_every)

    data_dir = resolve_path(cfg, cfg["output_dir"]) / "data"
    manifest = json.loads((data_dir / "manifest.json").read_text())
    pos_weight = torch.tensor([manifest["pos_weight"][t] for t in targets], dtype=torch.float)
    print("pos_weight:", {t: round(float(w), 3) for t, w in zip(targets, pos_weight)})

    # --- data ---------------------------------------------------------------
    tokenizer = get_tokenizer(model_id)
    bos_offset = detect_bos_offset(tokenizer)
    collate = Collator(tokenizer, bos_offset=bos_offset,
                       max_length=int(cfg["model"]["max_length"]))
    cw = cfg["preprocess"].get("count_weighting", "log")
    nw = int(tr.get("num_workers", 0))
    pin = bool(tr.get("pin_memory", False)) and device.type == "cuda"
    ds = {sp: SeqDataset(data_dir / f"{sp}.parquet", targets, count_weighting=cw)
          for sp in ("train", "val", "test")}
    print({sp: len(ds[sp]) for sp in ds})
    dl = {
        "train": DataLoader(ds["train"], batch_size=batch_size, shuffle=True,
                            collate_fn=collate, num_workers=nw, pin_memory=pin),
        "val": DataLoader(ds["val"], batch_size=batch_size, shuffle=False,
                          collate_fn=collate, num_workers=nw, pin_memory=pin),
        "test": DataLoader(ds["test"], batch_size=batch_size, shuffle=False,
                           collate_fn=collate, num_workers=nw, pin_memory=pin),
    }

    # --- model --------------------------------------------------------------
    model = MultiTaskESMC(
        model_id=model_id, targets=targets,
        pooling=pooling, dropout=float(cfg["model"]["dropout"]),
        pos_weight=pos_weight,
    ).to(device)

    # precision (autocast happens inside forward) + optional multi-GPU
    amp_dtype, scaler = resolve_precision(str(tr.get("precision", "auto")), device)
    model._amp_dtype = amp_dtype
    print(f"precision: {('fp32' if amp_dtype is None else str(amp_dtype).split('.')[-1])} "
          f"| loss-scaler: {scaler is not None}")
    if tr.get("gradient_checkpointing", False):
        ok = model.enable_gradient_checkpointing()
        print(f"gradient checkpointing: {'enabled' if ok else 'not supported by backbone'}")
    n_gpu = torch.cuda.device_count() if device.type == "cuda" else 0
    train_model = model
    if str(tr.get("multi_gpu", "false")).lower() == "auto" and n_gpu > 1:
        train_model = torch.nn.DataParallel(model)
        print(f"DataParallel across {n_gpu} GPUs")

    # === Phase 1: frozen backbone, train heads =============================
    model.freeze_backbone()
    t, total = model._count_trainable()
    print(f"\n[Phase 1] heads only - trainable {t:,}/{total:,} params")
    opt = build_optimizer(model, float(tr["phase1"]["lr"]), float(tr["weight_decay"]))
    for ep in range(p1_epochs):
        loss = run_train_epoch(train_model, dl["train"], opt, scaler, device, grad_clip,
                               desc=f"P1 epoch {ep+1}/{p1_epochs}", log_every=log_every)
        per, *_ = evaluate(model, dl["val"], device, targets)
        print(f"  P1 ep{ep+1}: train_loss={loss:.4f} val_{sel_metric}={mean_metric(per, sel_metric):.4f}")

    # === Phase 2: unfreeze top-N layers, fine-tune with early stopping =====
    t, total = model.set_trainable_top_layers(unfreeze)
    model.to(device)
    print(f"\n[Phase 2] unfreeze={unfreeze} - trainable {t:,}/{total:,} params")
    opt = build_optimizer(model, float(tr["phase2"]["lr"]), float(tr["weight_decay"]))
    patience = int(tr["phase2"].get("early_stop_patience", 3))
    best_score, best_state, bad = -np.inf, None, 0
    for ep in range(p2_epochs):
        loss = run_train_epoch(train_model, dl["train"], opt, scaler, device, grad_clip,
                               desc=f"P2 epoch {ep+1}/{p2_epochs}", log_every=log_every)
        per, *_ = evaluate(model, dl["val"], device, targets)
        score = mean_metric(per, sel_metric)
        print(f"  P2 ep{ep+1}: train_loss={loss:.4f} val_{sel_metric}={score:.4f}")
        if score > best_score:
            best_score = score
            best_state = copy.deepcopy(model.state_dict())
            bad = 0
        else:
            bad += 1
            if bad >= patience:
                print(f"  early stop (no val improvement for {patience} epochs)")
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    print(f"\nBest val {sel_metric}: {best_score:.4f}")

    # --- choose thresholds on validation -----------------------------------
    per_val, probs_val, labels_val = evaluate(model, dl["val"], device, targets)
    thresholds = {}
    for ti, t_ in enumerate(targets):
        m = labels_val[:, ti] != -100
        thresholds[t_] = best_threshold(labels_val[m, ti], probs_val[m, ti]) if m.sum() else 0.5

    # --- final, single evaluation on the held-out test set -----------------
    per_test, *_ = evaluate(model, dl["test"], device, targets, thresholds=thresholds)
    report = fmt_report(per_test, targets)
    print("\n=== HELD-OUT TEST REPORT ===")
    print(report)
    print("thresholds:", {k: round(v, 3) for k, v in thresholds.items()})

    # --- strict novel-loop subset (no train near-neighbour) ----------------
    per_novel, novel_report = None, None
    if strict_test:
        test_df = SeqDataset(data_dir / "test.parquet", targets, count_weighting=cw).df
        if "near_train_h1" in test_df.columns:
            novel = test_df[~test_df["near_train_h1"].astype(bool)]
            print(f"\n=== STRICT NOVEL-LOOP TEST ({len(novel)}/{len(test_df)} seqs, "
                  f"no train neighbour within 1 edit) ===")
            if len(novel):
                ndl = DataLoader(SeqDataset(novel, targets, count_weighting=cw),
                                 batch_size=batch_size, shuffle=False, collate_fn=collate,
                                 num_workers=nw, pin_memory=pin)
                per_novel, *_ = evaluate(model, ndl, device, targets, thresholds=thresholds)
                novel_report = fmt_report(per_novel, targets)
                print(novel_report)
        else:
            print("\n[info] --strict-test: re-run data_prep.py to add the near_train_h1 column.")

    # --- save ---------------------------------------------------------------
    if out_name is None:
        out_name = "model_smoke" if args.smoke else "model"
    out_dir = resolve_path(cfg, cfg["output_dir"]) / out_name
    out_dir.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), out_dir / "model_state.pt")
    meta = {
        "model_id": model_id,
        "targets": targets,
        "pooling": pooling,
        "dropout": float(cfg["model"]["dropout"]),
        "bos_offset": int(bos_offset),
        "max_length": int(cfg["model"]["max_length"]),
        "thresholds": thresholds,
        "pos_weight": {t: float(w) for t, w in zip(targets, pos_weight)},
        "val_best_metric": {sel_metric: float(best_score)},
    }
    (out_dir / "model_meta.json").write_text(json.dumps(meta, indent=2))
    report_txt = report + "\n\nthresholds: " + str(thresholds) + "\n"
    if novel_report:
        report_txt += "\n=== STRICT NOVEL-LOOP TEST ===\n" + novel_report + "\n"
    (out_dir / "test_report.txt").write_text(report_txt)
    (out_dir / "test_report.json").write_text(json.dumps(
        {"per_target": per_test, "per_target_novel": per_novel,
         "thresholds": thresholds}, indent=2))
    print(f"\nSaved model + meta + test report to {out_dir}")
    return {"per_test": per_test, "per_novel": per_novel, "thresholds": thresholds,
            "best_val": float(best_score), "out_dir": str(out_dir), "pooling": pooling}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--smoke", action="store_true", help="tiny fast run")
    ap.add_argument("--pooling", default=None, choices=["loop", "mean", "cls"],
                    help="override config pooling for this run")
    ap.add_argument("--model-id", default=None, help="override config model_id (e.g. ESMplusplus_small)")
    ap.add_argument("--out-name", default=None, help="override model output folder name")
    ap.add_argument("--strict-test", action="store_true",
                    help="also report on the novel-loop test subset (no train near-neighbour)")
    ap.add_argument("--log-every", type=int, default=None,
                    help="print a plain progress line every N steps (good for cluster logs)")
    args = ap.parse_args()
    cfg = load_config(args.config)
    train_and_evaluate(cfg, smoke=args.smoke, pooling=args.pooling, model_id=args.model_id,
                       out_name=args.out_name, strict_test=args.strict_test, log_every=args.log_every)


if __name__ == "__main__":
    main()
