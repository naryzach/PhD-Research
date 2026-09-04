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
PAIRS="MMP ADAM"                     # MMP9-vs-MMP2, ADAM17-vs-ADAM10 (flipped 2026-08-27)
# Conformation refinement, using the settings the generation campaign measured.
# partial_t 3.0->1.0 comes from partial_t_ladder.py: parent->child transmission is
# +0.87/+0.80/+0.69 at t=1/2/3 and collapses to +0.18/-0.31 at t=5/8.
CONFORMATION_MODE="${CONFORMATION_MODE:-1}"
BEAM_WIDTH="${BEAM_WIDTH:-10}"
PARTIAL_T="${PARTIAL_T:-3.0}"
PARTIAL_T_MIN="${PARTIAL_T_MIN:-1.0}"
BEAM_FRESH_FRACTION="${BEAM_FRESH_FRACTION:-0.25}"
SV_BATTERY="${SV_BATTERY:-1}"        # sv_* columns; the first campaign ran without them
LOOPS="AB C EF"
BACKBONES_PER_TARGET=50
# 1, not 3. This path folds the on- AND off-target, so 50x3 designs = 300 complexes per
# round -- double the generation campaign's batch, on the two LARGEST complexes (ADAM17
# 188+256 aa). Measured 2026-09-02/03: that batch does not finish inside three hours and
# every round was being cut off by the timeout.
# Cutting SEQUENCES rather than backbones is the right lever: the variance decomposition
# says loop CONFORMATION carries ~89% of sv_pdockq and loop SEQUENCE ~11%, so sampling 3
# sequences per backbone spent 3x the fold budget on the axis that matters least.
# 50 designs = 100 complexes/round, which fits comfortably.
SEQS_PER_BACKBONE="${SEQS_PER_BACKBONE:-1}"   # SEQS_PER_BACKBONE=3 restores the old batch
INIT_TEMPERATURE=0.60
MIN_TEMPERATURE=0.10
TEMP_DECAY=0.94
MAX_ITERATIONS=40
ESMFOLD2_GPUS="${ESMFOLD2_GPUS:-auto}"  # =1 (env) when co-running pinned to a different GPU
MAX_RETRIES=5

# ── Paths ────────────────────────────────────────────────────────────────────
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Whole-batch ESMFold2 timeout. This is a throughput budget, not a safety
# limit: 300 complexes/round (on+off per design) and the ADAM pair is the largest.
# A near-constant scored count each round means the batch is being cut off here.
export ESMFOLD2_TIMEOUT_S="${ESMFOLD2_TIMEOUT_S:-10800}"
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

if [[ "${FRESH:-0}" == "1" ]]; then
  # Archive the WHOLE prior run — state AND it_* dirs + hof summaries — not just the
  # state file. Moving only the state (the old behaviour) leaves the round dirs in
  # place, so the fresh run restarts at it_0 and OVERWRITES it_0/round_summary.csv
  # while it_1..it_N survive from the previous objective. Measured 2026-09-04: the
  # pool then reported 10 rounds with the state at iteration 1, blending the
  # contact-density and pdockq objectives in one directory with no way to tell them
  # apart. Same archiving as run_generation.sh. Moved, not deleted.
  if compgen -G "$OUT_BASE/it_*" > /dev/null || [[ -f "$STATE" ]]; then
    arch="$OUT_BASE/_archived_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$arch"
    [[ -f "$STATE" ]] && mv "$STATE" "$arch"/
    for s in specificity_hof_summary.csv hof_summary.csv; do
      [[ -f "$OUT_BASE/$s" ]] && mv "$OUT_BASE/$s" "$arch"/
    done
    for d in "$OUT_BASE"/it_*; do [[ -e "$d" ]] && mv "$d" "$arch"/; done
    echo "FRESH: archived prior run -> $arch (anneal starts hot at T=$INIT_TEMPERATURE)" | tee -a "$LOG"
  fi
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
# See run_generation.sh: a signal is an operator decision, not a transient crash.
_child=""
_shutdown() {
  echo "=== stopped by signal at $(date) ===" | tee -a "$LOG"
  [[ -n "$_child" ]] && kill "$_child" 2>/dev/null
  exit 130
}
trap _shutdown INT TERM

attempt=0
while (( attempt <= MAX_RETRIES )); do
  echo "=== launch attempt $attempt at $(date) ===" | tee -a "$LOG"
  conf_flags=()
  if [[ "$CONFORMATION_MODE" == "1" ]]; then
    conf_flags=(--conformation-mode --beam-width "$BEAM_WIDTH" --partial-t "$PARTIAL_T"
                --partial-t-min "$PARTIAL_T_MIN" --beam-fresh-fraction "$BEAM_FRESH_FRACTION")
  fi
  sv_flags=()
  (( SV_BATTERY )) && sv_flags+=(--sv-battery)
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
    ${conf_flags[@]+"${conf_flags[@]}"} \
    ${sv_flags[@]+"${sv_flags[@]}"} \
    >> "$LOG" 2>&1 &
  _child=$!
  wait "$_child"
  rc=$?
  if (( rc == 0 )); then
    echo "=== specificity generation completed cleanly at $(date) ===" | tee -a "$LOG"
    break
  fi
  if (( rc >= 128 )); then  # killed by a signal: deliberate, so stop
    echo "=== terminated by signal (rc=$rc) at $(date); not retrying ===" | tee -a "$LOG"
    exit "$rc"
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
