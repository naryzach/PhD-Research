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
#   bash Generation/run_specificity.sh              # start (or resume if state exists)
#   FRESH=1 bash Generation/run_specificity.sh      # archive old state, start clean
#
# Run it on different GPUs than run_generation.sh if you want both campaigns at once
# (ESMFOLD2_GPUS=auto skips GPUs already busy, but RFd3/LMPNN are single-GPU, so pin
# each launcher to its own device with CUDA_VISIBLE_DEVICES if you co-schedule).

set -uo pipefail

# ── Config (edit here) ───────────────────────────────────────────────────────
PAIRS="MMP ADAM"                     # MMP2-vs-MMP9, ADAM10-vs-ADAM17
LOOPS="AB C EF"
BACKBONES_PER_TARGET=50
SEQS_PER_BACKBONE=3
INIT_TEMPERATURE=0.60
MIN_TEMPERATURE=0.10
TEMP_DECAY=0.94
MAX_ITERATIONS=40
ESMFOLD2_GPUS=auto
MAX_RETRIES=5

# ── Paths ────────────────────────────────────────────────────────────────────
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_BASE="$HERE/../Local/specificity_refinement"
STATE="$OUT_BASE/specificity_state.json"
LOG_DIR="$OUT_BASE/logs"
mkdir -p "$LOG_DIR"
LOG="$LOG_DIR/spec_$(date +%Y%m%d_%H%M%S).log"

# ── Preflight (ESMFold2 subprocess python) ───────────────────────────────────
export ESMFOLD2_PYTHON="${ESMFOLD2_PYTHON:-$(command -v python)}"
echo "ESMFOLD2_PYTHON=$ESMFOLD2_PYTHON" | tee -a "$LOG"
if [[ ! -x "$ESMFOLD2_PYTHON" ]]; then
  echo "WARNING: '$ESMFOLD2_PYTHON' is not an executable python — ESMFold2 ranking will be SKIPPED." | tee -a "$LOG"
elif ! "$ESMFOLD2_PYTHON" -c "import esm" >/dev/null 2>&1; then
  echo "WARNING: '$ESMFOLD2_PYTHON' cannot 'import esm' — ESMFold2 ranking will be SKIPPED." | tee -a "$LOG"
  echo "         Fix: run with your ESMFold2-capable env active, or export ESMFOLD2_PYTHON=/path/to/env/bin/python" | tee -a "$LOG"
else
  echo "ESMFold2 backend (esm) import OK." | tee -a "$LOG"
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
