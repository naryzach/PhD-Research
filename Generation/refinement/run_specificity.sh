#!/usr/bin/env bash
# run_specificity.sh — one long, unattended TIMP3 SELECTIVITY-generation run.
#
# Same hot->cold anneal and crash-resume machinery as run_generation.sh, but for
# specificity_refinement.py: each design is folded against its ON- and OFF-target
# and ranked by (calibrated on-target quality) + (interface contact-density gap).
# AF3 is NOT interleaved — validate at the end.
#
# Usage:
#   conda activate foundry
#   bash Generation/refinement/run_specificity.sh              # start (or resume if state exists)
#   FRESH=1 bash Generation/refinement/run_specificity.sh      # archive old state, start clean
#
# Run it on different GPUs than run_generation.sh if you want both campaigns at once
# (ESMFOLD2_GPUS=auto skips GPUs already busy, but RFd3/LMPNN are single-GPU, so pin
# each launcher to its own device with CUDA_VISIBLE_DEVICES if you co-schedule).

set -uo pipefail

# Reduce CUDA fragmentation OOMs. And do NOT share GPUs with run_generation.sh —
# pin this to its own device(s) (CUDA_VISIBLE_DEVICES=...) or run it sequentially.
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

# ── Config (edit here) ───────────────────────────────────────────────────────
PAIRS="MMP ADAM"                     # MMP2-vs-MMP9, ADAM10-vs-ADAM17
LOOPS="AB C EF"
BACKBONES_PER_TARGET=50
SEQS_PER_BACKBONE=3
INIT_TEMPERATURE=0.60
MIN_TEMPERATURE=0.10
TEMP_DECAY=0.94
MAX_ITERATIONS=40
ESMFOLD2_GPUS="${ESMFOLD2_GPUS:-auto}"  # =1 (env) when co-running pinned to a different GPU
MAX_RETRIES=5

# ── Paths ────────────────────────────────────────────────────────────────────
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Output root — overridable via SPEC_OUT_BASE. Exported so specificity_refinement.py
# uses the SAME root (it redirects the inherited iterative OUT_BASE to this).
export SPEC_OUT_BASE="${SPEC_OUT_BASE:-$HERE/../../Local/specificity_refinement}"
OUT_BASE="$SPEC_OUT_BASE"
STATE="$OUT_BASE/specificity_state.json"
LOG_DIR="$OUT_BASE/logs"
mkdir -p "$LOG_DIR"
LOG="$LOG_DIR/spec_$(date +%Y%m%d_%H%M%S).log"

# ── Preflight (ESMFold2 subprocess python) ───────────────────────────────────
export ESMFOLD2_PYTHON="${ESMFOLD2_PYTHON:-$(command -v python)}"
export PYTHONUNBUFFERED=1   # see run_generation.sh — keeps the scoring log lines
                            # from being lost in the stderr buffer on a hard kill
echo "ESMFOLD2_PYTHON=$ESMFOLD2_PYTHON" | tee -a "$LOG"

# Reflexive fold test, hard gate — see run_generation.sh for the rationale.
# `import esm` is not evidence that a fold will succeed.
if [[ "${SKIP_ESM_PREFLIGHT:-0}" == "1" ]]; then
  echo "WARNING: SKIP_ESM_PREFLIGHT=1 — launching WITHOUT proving ESMFold2 can fold." | tee -a "$LOG"
else
  echo "Preflight: folding a test complex through the real scorer path (~1-2 min)..." | tee -a "$LOG"
  if python "$HERE/esmfold2_smoketest.py" >>"$LOG" 2>&1; then
    echo "Preflight OK: ESMFold2 folded a test complex." | tee -a "$LOG"
  else
    {
      echo "FATAL: ESMFold2 preflight fold FAILED — aborting before any GPU time is spent."
      echo "       interpreter: $ESMFOLD2_PYTHON"
      echo "       Diagnostics in: $LOG    Re-check: python $HERE/esmfold2_smoketest.py"
    } | tee -a "$LOG"
    exit 1
  fi
fi

if [[ "${FRESH:-0}" == "1" && -f "$STATE" ]]; then
  bak="$STATE.$(date +%Y%m%d_%H%M%S).bak"; mv "$STATE" "$bak"
  echo "FRESH: archived existing state -> $bak (anneal starts hot at T=$INIT_TEMPERATURE)" | tee -a "$LOG"
fi
if [[ -f "$STATE" ]]; then
  echo "RESUMING: $STATE exists — anneal continues from the saved iteration/temperature." | tee -a "$LOG"
else
  echo "STARTING FRESH: no state file — anneal begins at T=$INIT_TEMPERATURE." | tee -a "$LOG"
fi

echo "Log: $LOG"
echo "Pairs: $PAIRS | Loops: $LOOPS | ${BACKBONES_PER_TARGET}bb x ${SEQS_PER_BACKBONE}seq | " \
     "T:${INIT_TEMPERATURE}->${MIN_TEMPERATURE} decay ${TEMP_DECAY} | ${MAX_ITERATIONS} iters" | tee -a "$LOG"

# ── Run with crash-resume ────────────────────────────────────────────────────
attempt=0
while (( attempt <= MAX_RETRIES )); do
  echo "=== launch attempt $attempt at $(date) ===" | tee -a "$LOG"
  python "$HERE/specificity_refinement.py" \
    --pair $PAIRS \
    --loops $LOOPS \
    --backbones-per-target "$BACKBONES_PER_TARGET" \
    --seqs-per-backbone "$SEQS_PER_BACKBONE" \
    --init-temperature "$INIT_TEMPERATURE" \
    --min-temperature "$MIN_TEMPERATURE" \
    --temp-decay "$TEMP_DECAY" \
    --esmfold2-gpus "$ESMFOLD2_GPUS" \
    --max-iterations "$MAX_ITERATIONS" \
    >> "$LOG" 2>&1
  rc=$?
  if (( rc == 0 )); then
    echo "=== specificity generation completed cleanly at $(date) ===" | tee -a "$LOG"
    break
  fi
  if (( rc == 3 )); then   # EXIT_ESMFOLD2_DEAD — retrying cannot fix a dead ranker
    {
      echo "=== ABORTED at $(date): ESMFold2 scored 0 designs for a whole target ==="
      echo "    Not retrying — environment fault. Check: python $HERE/esmfold2_smoketest.py"
      echo "    Saved state is intact; relaunch (without FRESH=1) to resume."
    } | tee -a "$LOG"
    exit 3
  fi
  attempt=$((attempt+1))
  echo "!!! exited rc=$rc — resuming from saved state ($attempt/$MAX_RETRIES) after 30s" | tee -a "$LOG"
  sleep 30
done

if (( attempt > MAX_RETRIES )); then
  echo "GAVE UP after $MAX_RETRIES retries — check $LOG" | tee -a "$LOG"
  exit 1
fi

cat <<EOF | tee -a "$LOG"

Specificity generation done. HOF: $OUT_BASE/specificity_hof_summary.csv

AF3 selectivity validation (on+off jobs are auto-exported every 2 iters):
  # upload the latest $OUT_BASE/af3_specificity_it*.json, download, then:
  python $HERE/specificity_refinement.py --pair $PAIRS --import-af3 <results.zip>

Ordering the selective designs (cross-fold each candidate against all targets):
  python $HERE/select_binders_to_order.py --emit-crossfold-input
  # fold ordering/crossfold_input.csv with ESMFold2, then:
  python $HERE/select_binders_to_order.py --criteria best_specificity \\
      --specificity-scores ordering/crossfold_scores.csv
EOF
