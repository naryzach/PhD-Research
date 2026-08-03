#!/usr/bin/env bash
# run_generation.sh — one long, unattended TIMP3 binder-generation run.
#
# A single stretched simulated-annealing pass: LigandMPNN sampling temperature
# cools from HOT (broad exploration) to COLD (confident exploitation) over many
# iterations, ranked in-loop by ESMFold2 (see calibrated_scoring.py). AF3 is NOT
# interleaved — it's a weak, rate-limited prior, so we fold it in ONCE at the end
# via the stratified calibration tranche (see the tail of this file).
#
# Robust for "let it run on the server for a while":
#   * logs to Local/iterative_refinement/logs/gen_<timestamp>.log
#   * per-target failure isolation is built into the pipeline; state is saved every
#     iteration, so a crash loses at most the current iteration.
#   * auto-resumes on a non-zero exit (up to MAX_RETRIES). Because state persists,
#     re-invoking the SAME command continues the anneal from where it stopped
#     (--init-temperature is ignored on resume — the saved temperature wins).
#
# Usage:
#   conda activate foundry
#   export ESMFOLD2_PYTHON=/path/to/envs/esmfold2/bin/python
#   bash Generation/run_generation.sh              # start (or resume if state exists)
#   FRESH=1 bash Generation/run_generation.sh      # archive old state, start clean
#
# On a SLURM cluster, submit this with sbatch (request the GPUs + a long walltime).

set -uo pipefail

# Reduce CUDA fragmentation OOMs (the failure mode when GPUs are shared/tight).
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"
# If you hit the libcue_ops.so / cublas "undefined symbol" error (cuequivariance
# mismatch — it crippled the 2026-08-03 run: RFd3 crawled, ESMFold2 folds failed),
# export DISABLE_CUEQUIVARIANCE=1 before launching. The code auto-sets it for V100
# only, so set it yourself if your GPU has the mismatch but isn't a V100. Passed
# through to the python child here.
[[ -n "${DISABLE_CUEQUIVARIANCE:-}" ]] && export DISABLE_CUEQUIVARIANCE
# IMPORTANT: do NOT co-schedule this with run_specificity.sh on the SAME GPUs — they
# will starve each other's ESMFold2 (CUDA OOM) and silently leave designs unscored,
# hitting the last-processed targets hardest. Pin each run to its own device(s):
#   CUDA_VISIBLE_DEVICES=0,1 bash Generation/run_generation.sh
#   CUDA_VISIBLE_DEVICES=2,3 bash Generation/run_specificity.sh
# If you only have one GPU set, run them SEQUENTIALLY, not at once.

# ── Config (edit here) ───────────────────────────────────────────────────────
TARGETS="MMP2 MMP9 ADAM10 ADAM17"    # purchased/human = the calibratable set
LOOPS="AB C EF"
BACKBONES_PER_TARGET=50              # RFd3 is the slow stage; scales run time
SEQS_PER_BACKBONE=3                  # cheap LigandMPNN diversity + more ESMFold2 folds
INIT_TEMPERATURE=0.60               # HOT start (fresh run only)
MIN_TEMPERATURE=0.10                # COLD confident floor
TEMP_DECAY=0.94                     # slow cool -> ~29 iters to reach the floor
MAX_ITERATIONS=40                   # full anneal + ~11 exploit iterations at the floor
ESMFOLD2_GPUS="${ESMFOLD2_GPUS:-auto}"  # shard ESMFold2 across free GPUs; set =1 (env) when co-
                                    # running with another campaign pinned to a different GPU, since
                                    # 'auto' detects via nvidia-smi and ignores CUDA_VISIBLE_DEVICES
MAX_RETRIES=5                       # auto-restart budget on hard crashes
# ── SV structural battery (sv_bridge) ──
SV_BATTERY=1                        # 1 = log the full Structural-Validation interface battery (sv_* cols); 0 = off
SV_OCCLUSION_FILTER=1               # 1 = zero the composite for cleft-missing designs (mechanistic sanity GATE)
SV_OCCLUSION_MIN=""                 # optional: override catalytic-occlusion threshold (blank = default 0.15)

# ── Paths ────────────────────────────────────────────────────────────────────
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Output root — overridable via REFINE_OUT_BASE so a fresh anneal can run in its own
# directory without clobbering a preserved/salvaged pool. Exported so the python
# child (iterative_refinement.py) uses the SAME root.
export REFINE_OUT_BASE="${REFINE_OUT_BASE:-$HERE/../Local/iterative_refinement}"
OUT_BASE="$REFINE_OUT_BASE"
STATE="$OUT_BASE/refinement_state.json"
LOG_DIR="$OUT_BASE/logs"
mkdir -p "$LOG_DIR"
LOG="$LOG_DIR/gen_$(date +%Y%m%d_%H%M%S).log"

# ── Preflight ────────────────────────────────────────────────────────────────
# ESMFold2 runs as a subprocess under $ESMFOLD2_PYTHON. It does NOT need a separate
# env — if your ACTIVE env (e.g. foundry) already has ESMFold2, just let it default
# to the active interpreter. Only set ESMFOLD2_PYTHON yourself if ESMFold2 lives in
# a different env than the one you launch this from.
export ESMFOLD2_PYTHON="${ESMFOLD2_PYTHON:-$(command -v python)}"
echo "ESMFOLD2_PYTHON=$ESMFOLD2_PYTHON" | tee -a "$LOG"
if [[ ! -x "$ESMFOLD2_PYTHON" ]]; then
  echo "WARNING: '$ESMFOLD2_PYTHON' is not an executable python — ESMFold2 ranking will be SKIPPED." | tee -a "$LOG"
elif ! "$ESMFOLD2_PYTHON" -c "import esm" >/dev/null 2>&1; then
  echo "WARNING: '$ESMFOLD2_PYTHON' cannot 'import esm' — ESMFold2 ranking will be SKIPPED" | tee -a "$LOG"
  echo "         (designs would fall back to the foldability floor only)." | tee -a "$LOG"
  echo "         Fix: run with your ESMFold2-capable env active, or set" | tee -a "$LOG"
  echo "           export ESMFOLD2_PYTHON=/path/to/that/env/bin/python" | tee -a "$LOG"
else
  echo "ESMFold2 backend (esm) import OK." | tee -a "$LOG"
fi

if [[ "${FRESH:-0}" == "1" ]]; then
  # Archive the WHOLE prior run — state AND it_* dirs + hof_summary — so a fresh
  # anneal can't inherit stale iteration dirs (which would contaminate the pool
  # and the stratified export). Moved, not deleted.
  if compgen -G "$OUT_BASE/it_*" > /dev/null || [[ -f "$STATE" ]]; then
    arch="$OUT_BASE/_archived_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$arch"
    [[ -f "$STATE" ]] && mv "$STATE" "$arch"/
    [[ -f "$OUT_BASE/hof_summary.csv" ]] && mv "$OUT_BASE/hof_summary.csv" "$arch"/
    for d in "$OUT_BASE"/it_*; do [[ -e "$d" ]] && mv "$d" "$arch"/; done
    echo "FRESH: archived prior run -> $arch (anneal starts hot at T=$INIT_TEMPERATURE)" | tee -a "$LOG"
  fi
fi
if [[ -f "$STATE" ]]; then
  echo "RESUMING: $STATE exists — the anneal continues from the saved iteration/temperature." | tee -a "$LOG"
else
  echo "STARTING FRESH: no state file — anneal begins at T=$INIT_TEMPERATURE." | tee -a "$LOG"
fi

# Assemble optional SV-battery flags (log-only metrics + optional occlusion gate).
sv_flags=()
(( SV_BATTERY )) && sv_flags+=(--sv-battery)
(( SV_OCCLUSION_FILTER )) && sv_flags+=(--sv-occlusion-filter)
[[ -n "$SV_OCCLUSION_MIN" ]] && sv_flags+=(--sv-occlusion-min "$SV_OCCLUSION_MIN")

echo "Log: $LOG"
echo "Targets: $TARGETS | Loops: $LOOPS | ${BACKBONES_PER_TARGET}bb x ${SEQS_PER_BACKBONE}seq | " \
     "T:${INIT_TEMPERATURE}->${MIN_TEMPERATURE} decay ${TEMP_DECAY} | ${MAX_ITERATIONS} iters | " \
     "SV: ${sv_flags[*]:-none}" | tee -a "$LOG"

# ── Run with crash-resume ────────────────────────────────────────────────────
attempt=0
while (( attempt <= MAX_RETRIES )); do
  echo "=== launch attempt $attempt at $(date) ===" | tee -a "$LOG"
  python "$HERE/iterative_refinement.py" \
    --targets $TARGETS \
    --loops $LOOPS \
    --backbones-per-target "$BACKBONES_PER_TARGET" \
    --seqs-per-backbone "$SEQS_PER_BACKBONE" \
    --init-temperature "$INIT_TEMPERATURE" \
    --min-temperature "$MIN_TEMPERATURE" \
    --temp-decay "$TEMP_DECAY" \
    --esmfold2-gpus "$ESMFOLD2_GPUS" \
    --max-iterations "$MAX_ITERATIONS" \
    ${sv_flags[@]+"${sv_flags[@]}"} \
    >> "$LOG" 2>&1
  rc=$?
  if (( rc == 0 )); then
    echo "=== generation completed cleanly at $(date) ===" | tee -a "$LOG"
    break
  fi
  attempt=$((attempt+1))
  echo "!!! exited rc=$rc — resuming from saved state ($attempt/$MAX_RETRIES) after 30s" | tee -a "$LOG"
  sleep 30
done

if (( attempt > MAX_RETRIES )); then
  echo "GAVE UP after $MAX_RETRIES retries — check $LOG" | tee -a "$LOG"
  exit 1
fi

# ── Next steps (printed, not run — these need the AF3 server in the loop) ─────
cat <<EOF | tee -a "$LOG"

Generation done. HOF: $OUT_BASE/hof_summary.csv

STAGE 3 — AF3 calibration tranche (spans predicted strong->weak across the WHOLE pool):
  python $HERE/iterative_refinement.py --targets $TARGETS --stratified-export 30
  # upload af3_submission_stratified.json, download results, save as folds_dayN.zip here,
  # then: python $HERE/iterative_refinement.py --targets $TARGETS --import-af3 folds_dayN.zip
  # (60-design tranche = run the export/upload/import on two separate days)

STAGE 4 — order selection + calibrated-vs-legacy rerank report:
  python $HERE/select_binders_to_order.py --criteria all --n 15
EOF
