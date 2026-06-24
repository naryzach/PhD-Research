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


def run_train_epoch(model, loader, optimizer, scaler, device, grad_clip, desc):
    model.train()
    total, n = 0.0, 0
    for batch in tqdm(loader, desc=desc, leave=False):
        batch = to_device(batch, device)
        optimizer.zero_grad()
        with torch.amp.autocast("cuda", enabled=scaler is not None):
            out = model(**batch)
            loss = out["loss"]
        if scaler is not None:
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            torch.nn.utils.clip_grad_norm_(
                [p for p in model.parameters() if p.requires_grad], grad_clip)
            scaler.step(optimizer)
            scaler.update()
        else:
            loss.backward()
            torch.nn.utils.clip_grad_norm_(
                [p for p in model.parameters() if p.requires_grad], grad_clip)
            optimizer.step()
        bs = batch["input_ids"].size(0)
        total += loss.item() * bs
        n += bs
    return total / max(n, 1)


@torch.no_grad()
def evaluate(model, loader, device, targets, thresholds=None):
    model.eval()
    logits_all, labels_all = [], []
    use_amp = device.type == "cuda"
    for batch in loader:
        b = to_device(batch, device)
        with torch.amp.autocast("cuda", enabled=use_amp):
            out = model(**b)
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


def train_and_evaluate(cfg, smoke=False, pooling=None, out_name=None):
    """Run the full staged fine-tune + held-out test evaluation once.

    Returns a result dict (per-target test metrics, thresholds, best val score,
    output dir). ``pooling`` overrides ``cfg['model']['pooling']`` (used by the
    pooling-comparison harness); ``out_name`` overrides the model output folder.
    """
    set_seed(int(cfg["split"].get("seed", 42)))
    targets = cfg["targets"]
    device = get_device()
    pooling = pooling or cfg["model"]["pooling"]
    print(f"Device: {device} | model: {cfg['model']['model_id']} | pooling: {pooling}")

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

    data_dir = resolve_path(cfg, cfg["output_dir"]) / "data"
    manifest = json.loads((data_dir / "manifest.json").read_text())
    pos_weight = torch.tensor([manifest["pos_weight"][t] for t in targets], dtype=torch.float)
    print("pos_weight:", {t: round(float(w), 3) for t, w in zip(targets, pos_weight)})

    # --- data ---------------------------------------------------------------
    tokenizer = get_tokenizer(cfg["model"]["model_id"])
    bos_offset = detect_bos_offset(tokenizer)
    collate = Collator(tokenizer, bos_offset=bos_offset,
                       max_length=int(cfg["model"]["max_length"]))
    cw = cfg["preprocess"].get("count_weighting", "log")
    ds = {sp: SeqDataset(data_dir / f"{sp}.parquet", targets, count_weighting=cw)
          for sp in ("train", "val", "test")}
    print({sp: len(ds[sp]) for sp in ds})
    dl = {
        "train": DataLoader(ds["train"], batch_size=batch_size, shuffle=True, collate_fn=collate),
        "val": DataLoader(ds["val"], batch_size=batch_size, shuffle=False, collate_fn=collate),
        "test": DataLoader(ds["test"], batch_size=batch_size, shuffle=False, collate_fn=collate),
    }

    # --- model --------------------------------------------------------------
    model = MultiTaskESMC(
        model_id=cfg["model"]["model_id"], targets=targets,
        pooling=pooling, dropout=float(cfg["model"]["dropout"]),
        pos_weight=pos_weight,
    ).to(device)
    scaler = torch.amp.GradScaler("cuda") if (device.type == "cuda" and tr.get("fp16", True)) else None

    # === Phase 1: frozen backbone, train heads =============================
    model.freeze_backbone()
    t, total = model._count_trainable()
    print(f"\n[Phase 1] heads only - trainable {t:,}/{total:,} params")
    opt = build_optimizer(model, float(tr["phase1"]["lr"]), float(tr["weight_decay"]))
    for ep in range(p1_epochs):
        loss = run_train_epoch(model, dl["train"], opt, scaler, device, grad_clip,
                               desc=f"P1 epoch {ep+1}/{p1_epochs}")
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
        loss = run_train_epoch(model, dl["train"], opt, scaler, device, grad_clip,
                               desc=f"P2 epoch {ep+1}/{p2_epochs}")
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

    # --- save ---------------------------------------------------------------
    if out_name is None:
        out_name = "model_smoke" if args.smoke else "model"
    out_dir = resolve_path(cfg, cfg["output_dir"]) / out_name
    out_dir.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), out_dir / "model_state.pt")
    meta = {
        "model_id": cfg["model"]["model_id"],
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
    (out_dir / "test_report.txt").write_text(report + "\n\nthresholds: " + str(thresholds) + "\n")
    (out_dir / "test_report.json").write_text(json.dumps(
        {"per_target": per_test, "thresholds": thresholds}, indent=2))
    print(f"\nSaved model + meta + test report to {out_dir}")
    return {"per_test": per_test, "thresholds": thresholds,
            "best_val": float(best_score), "out_dir": str(out_dir), "pooling": pooling}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--smoke", action="store_true", help="tiny fast run")
    ap.add_argument("--pooling", default=None, choices=["loop", "mean", "cls"],
                    help="override config pooling for this run")
    ap.add_argument("--out-name", default=None, help="override model output folder name")
    args = ap.parse_args()
    cfg = load_config(args.config)
    train_and_evaluate(cfg, smoke=args.smoke, pooling=args.pooling, out_name=args.out_name)


if __name__ == "__main__":
    main()
