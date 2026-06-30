"""Exhaustively score every 20^6 C-loop 6-mer with the trained model.
Keeps top-K per target + the most target-SELECTIVE loops. Inference-only, batched.
Shard with --start/--end to split across GPUs/jobs."""
import argparse, heapq, time
from pathlib import Path
import numpy as np, pandas as pd, torch
from esmc_utils import load_config, resolve_path
from inference import load_trained_model

AAS = "ACDEFGHIKLMNPQRSTVWY"

def derive_cloop_template(cfg):
    d = cfg["data"]; df = pd.read_csv(resolve_path(cfg, d["csv_path"]))
    sub = "Ligand_Subtype"
    c = df[df[sub] == "C-loop"] if sub in df.columns else df
    seqs, loops = c[d["seq_col"]].astype(str), c[d["loop_col"]].astype(str)
    starts = [s.find(l) for s, l in zip(seqs, loops)]
    start = int(pd.Series([x for x in starts if x >= 0]).mode()[0])
    L = len(loops.iloc[0])
    masked = [s[:start] + "_"*L + s[start+L:] for s, st in zip(seqs, starts) if st == start]
    from collections import Counter
    uniq = Counter(masked)
    print(f"distinct C-loop scaffolds outside the window: {len(uniq)} (want 1)")
    return uniq.most_common(1)[0][0], start, L

def idx_to_loop(idx, L):
    out = []
    for _ in range(L):
        out.append(AAS[idx % 20]); idx //= 20
    return "".join(out)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--model-dir", default=None)
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

    templ_masked, start_char, L = derive_cloop_template(cfg)
    templ_seq = templ_masked.replace("_"*L, "A"*L)
    templ_ids = tokenizer(templ_seq, return_tensors="pt")["input_ids"][0].to(device)
    Ltok = templ_ids.shape[0]
    tok_lo = start_char + meta["bos_offset"]
    aa_id = {a: int(tokenizer.convert_tokens_to_ids(a)) for a in AAS}
    print(f"window: char {start_char}..{start_char+L-1} | token {tok_lo}..{tok_lo+L-1} | {Ltok} tokens")
    print(f"precision {model._amp_dtype} | range [{args.start},{args.end}) of {20**L}")

    heaps = {t: [] for t in targets}; sel = []; K = args.topk
    attn = torch.ones(Ltok, dtype=torch.long, device=device)
    cols = list(range(tok_lo, tok_lo + L))
    t0 = time.time(); done = 0; idx = args.start; b = 0
    with torch.no_grad():
        while idx < args.end:
            n = min(args.batch_size, args.end - idx)
            loops = [idx_to_loop(idx + j, L) for j in range(n)]
            ids = templ_ids.unsqueeze(0).repeat(n, 1)
            for c, p in enumerate(cols):
                ids[:, p] = torch.tensor([aa_id[loops[j][c]] for j in range(n)], device=device)
            ls = torch.full((n,), tok_lo, dtype=torch.long, device=device)
            ll = torch.full((n,), L, dtype=torch.long, device=device)
            out = model(input_ids=ids, attention_mask=attn.unsqueeze(0).expand(n, -1),
                        loop_start=ls, loop_len=ll)
            probs = torch.sigmoid(out["logits"].float()).cpu().numpy()
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

    out_dir = resolve_path(cfg, cfg["output_dir"]) / "enumeration" / (args.out or f"{args.start}_{args.end}")
    out_dir.mkdir(parents=True, exist_ok=True)
    for t in targets:
        pd.DataFrame(sorted(heaps[t], reverse=True), columns=[f"prob_{t}", "loop"]
                     ).to_csv(out_dir / f"top_{t}.csv", index=False)
    pd.DataFrame(sorted(sel, reverse=True),
                 columns=["margin", "loop", "top_target", "top_prob", "runnerup_prob"]
                 ).to_csv(out_dir / "top_selective.csv", index=False)
    print(f"\nDone {done:,} in {(time.time()-t0)/3600:.2f} h -> {out_dir}")

if __name__ == "__main__":
    main()