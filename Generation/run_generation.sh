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

# ── Config (edit here) ───────────────────────────────────────────────────────
TARGETS="MMP2 MMP9 ADAM10 ADAM17"    # purchased/human = the calibratable set
LOOPS="AB C EF"
BACKBONES_PER_TARGET=50              # RFd3 is the slow stage; scales run time
SEQS_PER_BACKBONE=3                  # cheap LigandMPNN diversity + more ESMFold2 folds
INIT_TEMPERATURE=0.60               # HOT start (fresh run only)
MIN_TEMPERATURE=0.10                # COLD confident floor
TEMP_DECAY=0.94                     # slow cool -> ~29 iters to reach the floor
MAX_ITERATIONS=40                   # full anneal + ~11 exploit iterations at the floor
ESMFOLD2_GPUS=auto                  # shard ESMFold2 across all free GPUs
MAX_RETRIES=5                       # auto-restart budget on hard crashes

# ── Paths ────────────────────────────────────────────────────────────────────
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_BASE="$HERE/../Local/iterative_refinement"
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

if [[ "${FRESH:-0}" == "1" && -f "$STATE" ]]; then
  bak="$STATE.$(date +%Y%m%d_%H%M%S).bak"
  mv "$STATE" "$bak"
  echo "FRESH: archived existing state -> $bak (anneal will start hot at T=$INIT_TEMPERATURE)" | tee -a "$LOG"
fi
if [[ -f "$STATE" ]]; then
  echo "RESUMING: $STATE exists — the anneal continues from the saved iteration/temperature." | tee -a "$LOG"
else
  echo "STARTING FRESH: no state file — anneal begins at T=$INIT_TEMPERATURE." | tee -a "$LOG"
fi

echo "Log: $LOG"
echo "Targets: $TARGETS | Loops: $LOOPS | ${BACKBONES_PER_TARGET}bb x ${SEQS_PER_BACKBONE}seq | " \
     "T:${INIT_TEMPERATURE}->${MIN_TEMPERATURE} decay ${TEMP_DECAY} | ${MAX_ITERATIONS} iters" | tee -a "$LOG"

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
