"""Exhaustively score every 20^6 6-mer loop with the trained model, for a given
loop window (C-loop or AB-loop -- pass --loop-subtype).
Keeps top-K per target + the most target-SELECTIVE loops. Inference-only, batched.
Shard with --start/--end to split across GPUs/jobs.

Checkpointing: progress (current index + the running top-K heaps) is saved to
<out_dir>/checkpoint.json every --checkpoint-every sequences (atomic write --
temp file + rename, so a crash mid-write can't corrupt it). Re-running the
same command automatically resumes from the last checkpoint instead of
restarting the whole sweep; pass --no-resume to force a clean restart."""
import argparse, heapq, json, os, time
from pathlib import Path
import numpy as np, pandas as pd, torch
from esmc_utils import load_config, resolve_path
from inference import load_trained_model
from loop_window import AAS, build_scorer, derive_loop_template, idx_to_loop

def save_checkpoint(path, start, end, idx, heaps, sel):
    tmp = path.with_suffix(".json.tmp")
    payload = {"start": start, "end": end, "idx": idx,
              "heaps": {t: h for t, h in heaps.items()}, "sel": sel}
    tmp.write_text(json.dumps(payload))
    os.replace(tmp, path)  # atomic on POSIX and Windows

def load_checkpoint(path, start, end, targets):
    if not path.exists():
        return None
    try:
        d = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as e:
        print(f"[warn] checkpoint at {path} unreadable ({e}); starting fresh")
        return None
    if d.get("start") != start or d.get("end") != end:
        print(f"[warn] checkpoint range [{d.get('start')},{d.get('end')}) doesn't match "
              f"requested [{start},{end}); starting fresh")
        return None
    heaps = {t: [tuple(x) for x in d["heaps"].get(t, [])] for t in targets}
    sel = [tuple(x) for x in d["sel"]]
    for h in heaps.values():
        heapq.heapify(h)
    heapq.heapify(sel)
    return d["idx"], heaps, sel

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--model-dir", default=None)
    ap.add_argument("--loop-subtype", default="C-loop",
                    help="Ligand_Subtype value that defines the window to enumerate "
                         "(e.g. 'C-loop' or 'AB-loop'); pass '' to use the whole CSV unfiltered")
    ap.add_argument("--template-csv", default=None,
                    help="CSV to derive the loop window from (default: cfg data.csv_path)")
    ap.add_argument("--batch-size", type=int, default=512)
    ap.add_argument("--topk", type=int, default=50000)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--end", type=int, default=20**6)
    ap.add_argument("--out", default=None)
    ap.add_argument("--log-every", type=int, default=200)
    ap.add_argument("--checkpoint-every", type=int, default=1_000_000,
                    help="save progress every N sequences (0 disables checkpointing)")
    ap.add_argument("--no-resume", action="store_true",
                    help="ignore any existing checkpoint and start this range from scratch")
    args = ap.parse_args()

    cfg = load_config(args.config)
    model_dir = args.model_dir or (resolve_path(cfg, cfg["output_dir"]) / "model")
    model, tokenizer, meta, device = load_trained_model(model_dir)
    targets = meta["targets"]
    if device.type == "cuda":
        cc = torch.cuda.get_device_capability(0)[0]
        model._amp_dtype = torch.bfloat16 if cc >= 8 else torch.float16

    subtype = args.loop_subtype or None
    templ_masked, start_char, L = derive_loop_template(cfg, loop_subtype=subtype,
                                                        csv_path=args.template_csv)
    score = build_scorer(model, tokenizer, meta, templ_masked, start_char, L,
                         batch_size=args.batch_size)
    print(f"window: char {start_char}..{start_char+L-1} | token {score.tok_lo}..{score.tok_lo+L-1} "
          f"| {score.Ltok} tokens")
    print(f"precision {model._amp_dtype} | range [{args.start},{args.end}) of {20**L}")

    tag = (subtype or "all").replace("-", "").lower()
    out_dir = resolve_path(cfg, cfg["output_dir"]) / "enumeration" / (args.out or f"{tag}_{args.start}_{args.end}")
    out_dir.mkdir(parents=True, exist_ok=True)
    ckpt_path = out_dir / "checkpoint.json"

    K = args.topk
    resumed = None if args.no_resume else load_checkpoint(ckpt_path, args.start, args.end, targets)
    if resumed:
        idx, heaps, sel = resumed
        done = idx - args.start
        print(f"[resume] picking up at {idx:,} ({done:,}/{args.end-args.start:,} already done "
              f"from a previous checkpoint)")
    else:
        idx = args.start
        heaps = {t: [] for t in targets}
        sel = []
        done = 0

    t0 = time.time(); b = 0; since_ckpt = 0
    while idx < args.end:
        n = min(args.batch_size, args.end - idx)
        loops = [idx_to_loop(idx + j, L) for j in range(n)]
        probs = score(loops)
        for j in range(n):
            pj = probs[j]
            for ti, t in enumerate(targets):
                h = heaps[t]; v = float(pj[ti])
                if len(h) < K: heapq.heappush(h, (v, loops[j]))
                elif v > h[0][0]: heapq.heapreplace(h, (v, loops[j]))
            o = np.argsort(pj)[::-1]; m = float(pj[o[0]] - pj[o[1]])
            rec = (m, loops[j], targets[o[0]], float(pj[o[0]]), float(pj[o[1]]))
            if len(sel) < K: heapq.heappush(sel, rec)
            elif m > sel[0][0]: heapq.heapreplace(sel, rec)
        idx += n; done += n; b += 1; since_ckpt += n
        if args.log_every and b % args.log_every == 0:
            el = time.time() - t0; r = done / el
            print(f"  {done:,}/{args.end-args.start:,}  {r:,.0f} seq/s  "
                  f"ETA {(args.end-idx)/max(r,1)/3600:.1f} h", flush=True)
        if args.checkpoint_every and since_ckpt >= args.checkpoint_every:
            save_checkpoint(ckpt_path, args.start, args.end, idx, heaps, sel)
            since_ckpt = 0
            print(f"  [checkpoint] {idx:,}/{args.end:,} saved -> {ckpt_path}", flush=True)

    (out_dir / "run_meta.json").write_text(
        json.dumps({"loop_subtype": subtype, "start_char": start_char, "loop_len": L,
                    "start": args.start, "end": args.end}, indent=2))
    for t in targets:
        pd.DataFrame(sorted(heaps[t], reverse=True), columns=[f"prob_{t}", "loop"]
                     ).to_csv(out_dir / f"top_{t}.csv", index=False)
    pd.DataFrame(sorted(sel, reverse=True),
                 columns=["margin", "loop", "top_target", "top_prob", "runnerup_prob"]
                 ).to_csv(out_dir / "top_selective.csv", index=False)
    if ckpt_path.exists():
        ckpt_path.unlink()  # sweep finished cleanly -- checkpoint no longer needed
    print(f"\nDone {done:,} in {(time.time()-t0)/3600:.2f} h -> {out_dir}")

if __name__ == "__main__":
    main()
