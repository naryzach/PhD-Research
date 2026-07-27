"""Exhaustively score every 20^6 6-mer loop with the trained model, for a given
loop window (C-loop or AB-loop -- pass --loop-subtype).
Keeps top-K per target + the most target-SELECTIVE loops. Inference-only, batched.
Shard with --start/--end to split across GPUs/jobs."""
import argparse, heapq, json, time
from pathlib import Path
import numpy as np, pandas as pd, torch
from esmc_utils import load_config, resolve_path
from inference import load_trained_model
from loop_window import AAS, build_scorer, derive_loop_template, idx_to_loop

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

    heaps = {t: [] for t in targets}; sel = []; K = args.topk
    t0 = time.time(); done = 0; idx = args.start; b = 0
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
        idx += n; done += n; b += 1
        if args.log_every and b % args.log_every == 0:
            el = time.time() - t0; r = done / el
            print(f"  {done:,}/{args.end-args.start:,}  {r:,.0f} seq/s  "
                  f"ETA {(args.end-idx)/max(r,1)/3600:.1f} h", flush=True)

    tag = (subtype or "all").replace("-", "").lower()
    out_dir = resolve_path(cfg, cfg["output_dir"]) / "enumeration" / (args.out or f"{tag}_{args.start}_{args.end}")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "run_meta.json").write_text(
        json.dumps({"loop_subtype": subtype, "start_char": start_char, "loop_len": L,
                    "start": args.start, "end": args.end}, indent=2))
    for t in targets:
        pd.DataFrame(sorted(heaps[t], reverse=True), columns=[f"prob_{t}", "loop"]
                     ).to_csv(out_dir / f"top_{t}.csv", index=False)
    pd.DataFrame(sorted(sel, reverse=True),
                 columns=["margin", "loop", "top_target", "top_prob", "runnerup_prob"]
                 ).to_csv(out_dir / "top_selective.csv", index=False)
    print(f"\nDone {done:,} in {(time.time()-t0)/3600:.2f} h -> {out_dir}")

if __name__ == "__main__":
    main()