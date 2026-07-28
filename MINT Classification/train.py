"""Staged fine-tuning of MINT for pair (TIMP3-variant, target) binder
classification.

Curriculum (identical shape to the ESM-C pipeline):
  Phase 1 - backbone frozen, train the randomly-initialised head (warm-up).
  Phase 2 - unfreeze a fraction of the backbone (train.phase2.freeze_percent),
            low LR, fine-tune with early stopping on the validation set.

A SINGLE held-out test set is untouched through both phases and scored once
at the end. Per-TARGET decision thresholds are chosen on validation (Fbeta or
MCC, see config) even though there's only one shared head, since different
targets can have different optimal cutoffs off the same score.

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
from contextlib import nullcontext
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader
from tqdm.auto import tqdm

from mint_utils import best_threshold, compute_metrics, get_device, load_config, mean_metric, resolve_path, set_seed
from model import MintBinderClassifier, PairCollator, PairDataset


def to_device(batch, device):
    return {k: (v.to(device) if torch.is_tensor(v) else v) for k, v in batch.items()}


def resolve_precision(precision, device):
    if device.type != "cuda":
        return None, None
    if precision == "auto":
        cc_major = torch.cuda.get_device_capability(0)[0]
        precision = "bf16" if cc_major >= 8 else "fp16"
    if precision == "bf16":
        return torch.bfloat16, None
    if precision == "fp16":
        return torch.float16, torch.amp.GradScaler("cuda")
    return None, None  # fp32


def run_train_epoch(model, loader, optimizer, scaler, device, grad_clip, amp_dtype, desc, log_every=0):
    model.train()
    total, n = 0.0, 0
    params = [p for p in model.parameters() if p.requires_grad]
    n_steps = len(loader)
    last_t, last_i = time.time(), 0
    for i, batch in enumerate(tqdm(loader, desc=desc, leave=False, disable=None)):
        b = to_device(batch, device)
        optimizer.zero_grad()
        amp_ctx = (torch.autocast(device.type, dtype=amp_dtype)
                  if (amp_dtype is not None and device.type == "cuda") else nullcontext())
        with amp_ctx:
            out = model(chains=b["chains"], chain_ids=b["chain_ids"],
                       labels=b["labels"], count_weight=b["count_weight"])
            loss = out["loss"]
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
        bs = b["labels"].size(0)
        total += loss.item() * bs
        n += bs
        if log_every and ((i + 1) % log_every == 0 or i + 1 == n_steps):
            now = time.time()
            rate = (i + 1 - last_i) / max(now - last_t, 1e-9)
            last_t, last_i = now, i + 1
            print(f"  {desc}: step {i+1}/{n_steps}  loss={loss.item():.4f}  "
                  f"{rate:.1f} it/s", flush=True)
    return total / max(n, 1)


@torch.no_grad()
def evaluate(model, loader, device, thresholds=None, fbeta_beta=0.5, amp_dtype=None):
    model.eval()
    amp_ctx = (torch.autocast(device.type, dtype=amp_dtype)
              if (amp_dtype is not None and device.type == "cuda") else nullcontext())
    logits_all, labels_all, target_names = [], [], []
    for batch in loader:
        b = to_device(batch, device)
        with amp_ctx:
            out = model(chains=b["chains"], chain_ids=b["chain_ids"])
        logits_all.append(out["logits"].float().cpu())
        labels_all.append(batch["labels"])
        target_names.extend(batch["target_name"])
    logits = torch.cat(logits_all).numpy()
    labels = torch.cat(labels_all).numpy()
    probs = 1.0 / (1.0 + np.exp(-logits))
    target_names = np.array(target_names)

    per = {}
    for t in sorted(set(target_names.tolist())):
        m = target_names == t
        if m.sum() == 0:
            continue
        thr = 0.5 if thresholds is None else thresholds.get(t, 0.5)
        per[t] = compute_metrics(labels[m], probs[m], threshold=thr, beta=fbeta_beta)
    return per, probs, labels, target_names


def build_optimizer(model, lr, weight_decay):
    params = [p for p in model.parameters() if p.requires_grad]
    return torch.optim.AdamW(params, lr=lr, weight_decay=weight_decay)


def fmt_report(per, targets, beta=0.5):
    lines = [f"{'target':<8} {'n':>6} {'pos':>6} {'pos_rate':>8} "
             f"{'pr_auc':>7} {'roc_auc':>7} {'mcc':>7} {'f1':>7} {'f'+str(beta):>7}"]
    for t in targets:
        if t not in per:
            continue
        m = per[t]
        def g(k):
            v = m[k]
            return f"{v:.3f}" if v == v else "  nan"
        lines.append(f"{t:<8} {m['n']:>6} {m['n_pos']:>6} "
                     f"{(m['pos_rate'] if m['pos_rate']==m['pos_rate'] else 0):>8.3f} "
                     f"{g('pr_auc'):>7} {g('roc_auc'):>7} {g('mcc'):>7} {g('f1'):>7} {g('fbeta'):>7}")
    return "\n".join(lines)


def train_and_evaluate(cfg, smoke=False, out_name=None, strict_test=False, log_every=None):
    set_seed(int(cfg["split"].get("seed", 42)))
    device = get_device()
    print(f"Device: {device}")
    if device.type == "cpu":
        print("\n" + "!" * 70 + "\n"
              "WARNING: no GPU detected -- fine-tuning the 650M MINT backbone on CPU\n"
              "is impractically slow. Check torch's CUDA build vs the node's driver.\n" + "!" * 70 + "\n")

    args = type("A", (), {"smoke": smoke})()
    tr = cfg["train"]
    smoke_cfg = cfg.get("smoke", {})
    mc = cfg["mint"]
    batch_size = smoke_cfg.get("batch_size", tr["batch_size"]) if args.smoke else tr["batch_size"]
    p1_epochs = smoke_cfg.get("phase1_epochs", tr["phase1"]["epochs"]) if args.smoke else tr["phase1"]["epochs"]
    p2_epochs = smoke_cfg.get("phase2_epochs", tr["phase2"]["epochs"]) if args.smoke else tr["phase2"]["epochs"]
    freeze_pct = (smoke_cfg.get("freeze_percent", tr["phase2"]["freeze_percent"])
                 if args.smoke else tr["phase2"]["freeze_percent"])
    grad_clip = float(tr.get("grad_clip", 1.0))
    sel_metric = tr.get("metric", "pr_auc")
    threshold_metric = tr.get("threshold_metric", "mcc")
    fbeta_beta = float(tr.get("fbeta_beta", 0.5))
    log_every = int(tr.get("log_every", 0)) if log_every is None else int(log_every)

    data_dir = resolve_path(cfg, cfg["output_dir"]) / "data"
    manifest = json.loads((data_dir / "manifest.json").read_text())
    pos_weight = torch.tensor([manifest["pos_weight"]["overall"]], dtype=torch.float)
    print("pos_weight (overall):", round(float(pos_weight.item()), 3))

    # --- model (builds mint's tokenizer/collate internally) -----------------
    model = MintBinderClassifier(
        mint_cfg_json=str(resolve_path(cfg, mc["config_json"])),
        checkpoint_path=str(resolve_path(cfg, mc["checkpoint_path"])),
        device=device, sep_chains=bool(mc.get("sep_chains", True)),
        use_multimer=bool(mc.get("use_multimer", True)),
        truncation_seq_length=int(mc.get("truncation_seq_length", 512)),
        head_hidden=int(cfg["model"]["head_hidden"]), dropout=float(cfg["model"]["dropout"]),
        pos_weight=pos_weight, freeze_percent=1.0,
    ).to(device)
    collate = PairCollator(model.mint_collate)

    # --- data -----------------------------------------------------------
    cw = cfg["preprocess"].get("count_weighting", "log")
    nw = int(tr.get("num_workers", 0))
    pin = bool(tr.get("pin_memory", False)) and device.type == "cuda"
    ds = {sp: PairDataset(data_dir / f"{sp}.parquet", count_weighting=cw)
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

    amp_dtype, scaler = resolve_precision(str(tr.get("precision", "auto")), device)
    print(f"precision: {('fp32' if amp_dtype is None else str(amp_dtype).split('.')[-1])} "
          f"| loss-scaler: {scaler is not None}")
    if tr.get("gradient_checkpointing", False):
        ok = model.enable_gradient_checkpointing()
        print(f"gradient checkpointing: {'enabled' if ok else 'not supported by backbone'}")

    # === Phase 1: frozen backbone, train head ===============================
    model.freeze_backbone()
    t, total = model._count_trainable()
    print(f"\n[Phase 1] heads only - trainable {t:,}/{total:,} params")
    opt = build_optimizer(model, float(tr["phase1"]["lr"]), float(tr["weight_decay"]))
    for ep in range(p1_epochs):
        loss = run_train_epoch(model, dl["train"], opt, scaler, device, grad_clip, amp_dtype,
                               desc=f"P1 epoch {ep+1}/{p1_epochs}", log_every=log_every)
        per, *_ = evaluate(model, dl["val"], device, fbeta_beta=fbeta_beta, amp_dtype=amp_dtype)
        print(f"  P1 ep{ep+1}: train_loss={loss:.4f} val_{sel_metric}={mean_metric(per, sel_metric):.4f}")

    # === Phase 2: unfreeze a fraction of the backbone, fine-tune ============
    t, total = model.set_trainable_by_freeze_percent(freeze_pct)
    model.to(device)
    print(f"\n[Phase 2] freeze_percent={freeze_pct} "
          f"({model.count_backbone_layers()} backbone layers) - trainable {t:,}/{total:,} params")
    opt = build_optimizer(model, float(tr["phase2"]["lr"]), float(tr["weight_decay"]))
    patience = int(tr["phase2"].get("early_stop_patience", 3))
    best_score, best_state, bad = -np.inf, None, 0
    for ep in range(p2_epochs):
        loss = run_train_epoch(model, dl["train"], opt, scaler, device, grad_clip, amp_dtype,
                               desc=f"P2 epoch {ep+1}/{p2_epochs}", log_every=log_every)
        per, *_ = evaluate(model, dl["val"], device, fbeta_beta=fbeta_beta, amp_dtype=amp_dtype)
        score = mean_metric(per, sel_metric)
        print(f"  P2 ep{ep+1}: train_loss={loss:.4f} val_{sel_metric}={score:.4f}")
        if score > best_score:
            best_score, best_state, bad = score, copy.deepcopy(model.state_dict()), 0
        else:
            bad += 1
            if bad >= patience:
                print(f"  early stop (no val improvement for {patience} epochs)")
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    print(f"\nBest val {sel_metric}: {best_score:.4f}")

    # --- choose per-target thresholds on validation -------------------------
    per_val, probs_val, labels_val, tnames_val = evaluate(model, dl["val"], device, fbeta_beta=fbeta_beta, amp_dtype=amp_dtype)
    targets = sorted(set(tnames_val.tolist()))
    thresholds = {}
    for t in targets:
        m = tnames_val == t
        thresholds[t] = (best_threshold(labels_val[m], probs_val[m], metric=threshold_metric, beta=fbeta_beta)
                         if m.sum() else 0.5)
    print(f"threshold_metric: {threshold_metric}" +
         (f" (beta={fbeta_beta})" if threshold_metric == "fbeta" else ""))

    # --- final, single evaluation on the held-out test set -------------------
    per_test, *_ = evaluate(model, dl["test"], device, thresholds=thresholds, fbeta_beta=fbeta_beta, amp_dtype=amp_dtype)
    report = fmt_report(per_test, targets, beta=fbeta_beta)
    print("\n=== HELD-OUT TEST REPORT ===")
    print(report)
    print("thresholds:", {k: round(v, 3) for k, v in thresholds.items()})

    # --- strict novel-loop subset (no train near-neighbour) ------------------
    per_novel, novel_report = None, None
    if strict_test:
        test_df = PairDataset(data_dir / "test.parquet", count_weighting=cw).df
        if "near_train_h1" in test_df.columns:
            novel = test_df[~test_df["near_train_h1"].astype(bool)]
            print(f"\n=== STRICT NOVEL-LOOP TEST ({len(novel)}/{len(test_df)} rows, "
                  f"no train neighbour within 1 edit) ===")
            if len(novel):
                ndl = DataLoader(PairDataset(novel, count_weighting=cw), batch_size=batch_size,
                                 shuffle=False, collate_fn=collate, num_workers=nw, pin_memory=pin)
                per_novel, *_ = evaluate(model, ndl, device, thresholds=thresholds, fbeta_beta=fbeta_beta, amp_dtype=amp_dtype)
                novel_report = fmt_report(per_novel, targets, beta=fbeta_beta)
                print(novel_report)
        else:
            print("\n[info] --strict-test: re-run data_prep.py to add the near_train_h1 column.")

    # --- save -----------------------------------------------------------
    if out_name is None:
        out_name = "model_smoke" if args.smoke else "model"
    out_dir = resolve_path(cfg, cfg["output_dir"]) / out_name
    out_dir.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), out_dir / "model_state.pt")
    meta = {
        "mint_config_json": str(resolve_path(cfg, mc["config_json"])),
        "mint_checkpoint_path": str(resolve_path(cfg, mc["checkpoint_path"])),
        "sep_chains": bool(mc.get("sep_chains", True)),
        "use_multimer": bool(mc.get("use_multimer", True)),
        "truncation_seq_length": int(mc.get("truncation_seq_length", 512)),
        "head_hidden": int(cfg["model"]["head_hidden"]),
        "dropout": float(cfg["model"]["dropout"]),
        "targets": targets,
        "target_sequences": manifest.get("target_sequences", {}),
        "thresholds": thresholds,
        "threshold_metric": threshold_metric,
        "fbeta_beta": fbeta_beta,
        "pos_weight_overall": float(pos_weight.item()),
        "val_best_metric": {sel_metric: float(best_score)},
    }
    (out_dir / "model_meta.json").write_text(json.dumps(meta, indent=2))
    report_txt = report + "\n\nthresholds: " + str(thresholds) + "\n"
    if novel_report:
        report_txt += "\n=== STRICT NOVEL-LOOP TEST ===\n" + novel_report + "\n"
    (out_dir / "test_report.txt").write_text(report_txt)
    (out_dir / "test_report.json").write_text(json.dumps(
        {"per_target": per_test, "per_target_novel": per_novel, "thresholds": thresholds}, indent=2))
    print(f"\nSaved model + meta + test report to {out_dir}")
    return {"per_test": per_test, "per_novel": per_novel, "thresholds": thresholds,
            "best_val": float(best_score), "out_dir": str(out_dir)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--smoke", action="store_true", help="tiny fast run")
    ap.add_argument("--out-name", default=None, help="override model output folder name")
    ap.add_argument("--strict-test", action="store_true",
                    help="also report on the novel-loop test subset (no train near-neighbour)")
    ap.add_argument("--log-every", type=int, default=None,
                    help="print a plain progress line every N steps (good for cluster logs)")
    args = ap.parse_args()
    cfg = load_config(args.config)
    train_and_evaluate(cfg, smoke=args.smoke, out_name=args.out_name,
                       strict_test=args.strict_test, log_every=args.log_every)


if __name__ == "__main__":
    main()
