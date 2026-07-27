"""Wrapper that runs the full ESM-C pipeline (fine-tune -> SHAP hotspots ->
64M-combination enumeration -> figures) across every dataset variant.

Variants (see configs/*.yaml for exact data/target wiring; run --list for a
one-line summary of each):

  all3_original        TIMP_binder_all.csv, 3 targets, AB-loop + C-loop mixed
                        (replicates the original run)
  abloop_only           ... filtered to AB-loop rows only (MMP3, MMP9)
  cloop_only             ... filtered to C-loop rows only (ADAM17, MMP9)
  mmp9_other            TIMP_binder_Cloop_IHTE_CGLK.csv (separate MMP9-only,
                        C-loop dataset)
  everything_combined   TIMP_binder_all.csv + TIMP_binder_Cloop_IHTE_CGLK.csv
                        unioned, 3 targets, AB-loop + C-loop mixed

The two "mixed" variants (all3_original, everything_combined) train on both
loop windows at once, so their SHAP + enumeration steps run TWICE -- once per
window (AB-loop, C-loop) -- since a single 6-position sweep/explanation can
only cover one graft window at a time.

Steps per (variant, loop_subtype):
  data_prep -> train (once per variant, not per subtype) -> shap_hotspots ->
  enumerate_cloop (full 20^6 sweep) -> analyze_enumeration

Usage
-----
    python multirun/prepare_variant_datasets.py     # one-time, or let run_all do it
    python multirun/run_all.py                      # everything, all variants
    python multirun/run_all.py --list                # show the plan, run nothing
    python multirun/run_all.py --variants abloop_only,cloop_only
    python multirun/run_all.py --steps data,train                 # just fine-tune
    python multirun/run_all.py --steps shap,enumerate,analyze     # redo analysis only
    python multirun/run_all.py --smoke                            # tiny sanity pass, all steps
    python multirun/run_all.py --dry-run                          # print commands, run nothing
    python multirun/run_all.py --force                            # ignore existing outputs, rerun

Run from anywhere -- paths are resolved relative to this repo. Each step logs
to <variant_output_dir>/logs/<step>[_<subtype>].log as well as stdout, so a
crashed multi-day cluster job leaves a readable trail per step.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

ESMC_DIR = Path(__file__).resolve().parent.parent   # ESM-C Classification/
MULTIRUN_DIR = Path(__file__).resolve().parent       # ESM-C Classification/multirun/
sys.path.insert(0, str(ESMC_DIR))
from esmc_utils import load_config, resolve_path  # noqa: E402

VARIANTS = [
    {"name": "all3_original", "loop_subtypes": ["AB-loop", "C-loop"]},
    {"name": "abloop_only", "loop_subtypes": ["AB-loop"]},
    {"name": "cloop_only", "loop_subtypes": ["C-loop"]},
    {"name": "mmp9_other", "loop_subtypes": ["C-loop"]},
    {"name": "everything_combined", "loop_subtypes": ["AB-loop", "C-loop"]},
]
VARIANT_NAMES = [v["name"] for v in VARIANTS]
DERIVED_CSVS = ["abloop_only.csv", "cloop_only.csv", "everything_combined.csv"]


def tag_of(subtype: str) -> str:
    return subtype.replace("-", "").lower()


def run(cmd: list[str], log_path: Path, dry_run: bool) -> None:
    print(f"\n$ {' '.join(cmd)}\n  (cwd={ESMC_DIR}, log={log_path})")
    if dry_run:
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with open(log_path, "w", encoding="utf-8") as fh:
        proc = subprocess.Popen(cmd, cwd=ESMC_DIR, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in proc.stdout:
            print(line, end="")
            fh.write(line)
        rc = proc.wait()
    el = time.time() - t0
    if rc != 0:
        raise SystemExit(f"[FAILED] {' '.join(cmd)} (exit {rc}, {el/60:.1f} min) -- see {log_path}")
    print(f"  [ok] {el/60:.1f} min")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", default=None,
                    help=f"comma-separated subset of: {','.join(VARIANT_NAMES)} (default: all)")
    ap.add_argument("--steps", default="prepare,data,train,shap,enumerate,analyze",
                    help="comma-separated subset of: prepare,data,train,shap,enumerate,analyze")
    ap.add_argument("--list", action="store_true", help="print the run plan and exit")
    ap.add_argument("--dry-run", action="store_true", help="print commands without running them")
    ap.add_argument("--force", action="store_true", help="rerun steps even if outputs already exist")
    ap.add_argument("--smoke", action="store_true",
                    help="tiny fast pass through every step (data_prep --limit, train --smoke)")
    ap.add_argument("--strict-test", action="store_true", help="passed through to train.py")
    ap.add_argument("--log-every", type=int, default=None, help="passed through to train.py")
    # enumerate_cloop.py passthrough
    ap.add_argument("--topk", type=int, default=50000, help="passed through to enumerate_cloop.py")
    ap.add_argument("--enum-batch-size", type=int, default=512)
    ap.add_argument("--enum-start", type=int, default=0, help="shard start (passed to enumerate_cloop.py)")
    ap.add_argument("--enum-end", type=int, default=20**6, help="shard end (passed to enumerate_cloop.py)")
    # shap_hotspots.py passthrough
    ap.add_argument("--n-explain", type=int, default=300)
    ap.add_argument("--background-size", type=int, default=50)
    args = ap.parse_args()

    variants = args.variants.split(",") if args.variants else VARIANT_NAMES
    unknown = set(variants) - set(VARIANT_NAMES)
    if unknown:
        raise SystemExit(f"Unknown variant(s): {unknown} -- choices: {VARIANT_NAMES}")
    steps = set(args.steps.split(","))

    plan = [v for v in VARIANTS if v["name"] in variants]
    print("=== Multirun plan ===")
    for v in plan:
        print(f"  {v['name']:<22} loop_subtypes={v['loop_subtypes']}  steps={sorted(steps)}")
    if args.list:
        return

    if "prepare" in steps:
        derived_dir = ESMC_DIR.parent / "Local" / "esmc_multirun" / "derived_data"
        missing = [f for f in DERIVED_CSVS if not (derived_dir / f).exists()]
        if missing or args.force:
            run([sys.executable, str(MULTIRUN_DIR / "prepare_variant_datasets.py")],
                MULTIRUN_DIR / "logs" / "prepare_variant_datasets.log", args.dry_run)
        else:
            print("[skip] derived_data/*.csv already present (--force to regenerate)")

    for v in plan:
        name = v["name"]
        cfg_rel = f"multirun/configs/{name}.yaml"
        cfg_abs = MULTIRUN_DIR / "configs" / f"{name}.yaml"
        cfg = load_config(cfg_abs)
        out_dir = resolve_path(cfg, cfg["output_dir"])
        log_dir = out_dir / "logs"
        print(f"\n{'='*70}\nVariant: {name}  ->  {out_dir}\n{'='*70}")

        if "data" in steps:
            manifest = out_dir / "data" / "manifest.json"
            if manifest.exists() and not args.force:
                print(f"[skip] data_prep already done ({manifest})")
            else:
                cmd = [sys.executable, "data_prep.py", "--config", cfg_rel]
                if args.smoke:
                    cmd += ["--limit", "1500"]
                run(cmd, log_dir / "data_prep.log", args.dry_run)

        if "train" in steps:
            model_file = out_dir / ("model_smoke" if args.smoke else "model") / "model_state.pt"
            if model_file.exists() and not args.force:
                print(f"[skip] train already done ({model_file})")
            else:
                cmd = [sys.executable, "train.py", "--config", cfg_rel]
                if args.smoke:
                    cmd += ["--smoke"]
                if args.strict_test:
                    cmd += ["--strict-test"]
                if args.log_every is not None:
                    cmd += ["--log-every", str(args.log_every)]
                run(cmd, log_dir / "train.log", args.dry_run)

        for subtype in v["loop_subtypes"]:
            tag = tag_of(subtype)
            is_full_sweep = (args.enum_start == 0 and args.enum_end == 20**6)
            out_name = f"{tag}_full" if is_full_sweep else f"{tag}_{args.enum_start}_{args.enum_end}"

            if "shap" in steps:
                shap_summary = out_dir / "shap" / tag / "summary.json"
                if shap_summary.exists() and not args.force:
                    print(f"[skip] shap[{subtype}] already done ({shap_summary})")
                else:
                    cmd = [sys.executable, "shap_hotspots.py", "--config", cfg_rel,
                          "--loop-subtype", subtype, "--n-explain", str(args.n_explain),
                          "--background-size", str(args.background_size)]
                    run(cmd, log_dir / f"shap_{tag}.log", args.dry_run)

            if "enumerate" in steps:
                enum_done = out_dir / "enumeration" / out_name / "top_selective.csv"
                if enum_done.exists() and not args.force:
                    print(f"[skip] enumerate[{subtype}] already done ({enum_done})")
                else:
                    cmd = [sys.executable, "enumerate_cloop.py", "--config", cfg_rel,
                          "--loop-subtype", subtype, "--topk", str(args.topk),
                          "--batch-size", str(args.enum_batch_size),
                          "--start", str(args.enum_start), "--end", str(args.enum_end),
                          "--out", out_name]
                    run(cmd, log_dir / f"enumerate_{tag}.log", args.dry_run)

            if "analyze" in steps:
                analysis_done = out_dir / "enumeration" / out_name / "analysis" / "summary.json"
                if analysis_done.exists() and not args.force:
                    print(f"[skip] analyze[{subtype}] already done ({analysis_done})")
                else:
                    cmd = [sys.executable, "analyze_enumeration.py", "--config", cfg_rel,
                          "--dir", out_name, "--loop-subtype", subtype]
                    run(cmd, log_dir / f"analyze_{tag}.log", args.dry_run)

    print("\nAll requested steps complete.")


if __name__ == "__main__":
    main()
