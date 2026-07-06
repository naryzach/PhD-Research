#!/bin/bash
# Pre-flight check for MetalBinder pipeline dependencies.
# Run with: bash MetalBinder/check_pipeline_env.sh
# All checks must pass before submitting the SLURM job.

PASS=0; FAIL=0

ok()   { echo "  [OK]   $1"; ((PASS++)); }
fail() { echo "  [FAIL] $1"; ((FAIL++)); }
info() { echo ""; echo "=== $1 ==="; }

# ── Active env ─────────────────────────────────────────────────────────────────
info "Active conda environment"
if [[ "$CONDA_DEFAULT_ENV" == "protdesign" ]]; then
    ok "conda env: protdesign"
else
    fail "expected 'protdesign', got '${CONDA_DEFAULT_ENV:-none}' — run: conda activate protdesign"
fi

# ── Core Python packages ────────────────────────────────────────────────────────
info "Core packages (protdesign env)"
for pkg in torch torchvision numpy scipy lightning rfd3 biotite; do
    python -c "import $pkg; print('  [OK]   $pkg', getattr($pkg, '__version__', '?'))" 2>/dev/null \
        || fail "$pkg not importable"
    ((PASS++)) 2>/dev/null; true  # suppress arithmetic error; ok() already increments
done

# torch/torchvision version match
TV=$(python -c "import torchvision; print(torchvision.__version__.split('+')[0])" 2>/dev/null)
T=$(python -c "import torch; print(torch.__version__.split('+')[0])" 2>/dev/null)
T_MINOR=$(echo "$T" | cut -d. -f2)
TV_MINOR=$(echo "$TV" | cut -d. -f2)
EXPECTED_TV_MINOR=$((T_MINOR + 15))
if [[ "$TV_MINOR" == "$EXPECTED_TV_MINOR" ]]; then
    ok "torch $T <-> torchvision $TV version match"
else
    fail "torch $T expects torchvision 0.${EXPECTED_TV_MINOR}.x, got $TV — fix: pip install torchvision==0.${EXPECTED_TV_MINOR}.0"
fi

# numpy >= 2.0
NP_MAJOR=$(python -c "import numpy; print(numpy.__version__.split('.')[0])" 2>/dev/null)
if [[ "$NP_MAJOR" -ge 2 ]]; then
    NP_VER=$(python -c "import numpy; print(numpy.__version__)" 2>/dev/null)
    ok "numpy $NP_VER (>= 2.0)"
else
    fail "numpy < 2.0 — fix: pip install 'numpy>=2.0'"
fi

# GPU accessible
python -c "import torch; assert torch.cuda.is_available(), 'no GPU'" 2>/dev/null \
    && ok "CUDA GPU available ($(python -c 'import torch; print(torch.cuda.get_device_name(0))' 2>/dev/null))" \
    || fail "CUDA GPU not available — ensure SLURM gres=gpu:1 is set"

# ── Chai-1 env ──────────────────────────────────────────────────────────────────
info "Chai-1 subprocess env"
CHAI1_PY="${CHAI1_PYTHON:-$HOME/miniconda3/envs/chai1/bin/python}"
if [[ -x "$CHAI1_PY" ]]; then
    ok "Chai-1 python found: $CHAI1_PY"
    CHAI_VER=$("$CHAI1_PY" -c "import chai_lab; print(chai_lab.__version__)" 2>/dev/null)
    if [[ -n "$CHAI_VER" ]]; then
        ok "chai_lab $CHAI_VER importable"
    else
        fail "chai_lab not importable in $CHAI1_PY — fix: $CHAI1_PY -m pip install chai-lab"
    fi
else
    fail "Chai-1 python not found at: $CHAI1_PY"
    echo "         Fix:"
    echo "           conda create -n chai1 python=3.10 -y"
    echo "           conda activate chai1 && pip install chai-lab && conda deactivate"
    echo "         Then add to SLURM script: export CHAI1_PYTHON=$CHAI1_PY"
fi

# ── Input files ─────────────────────────────────────────────────────────────────
info "Input files"
INPUT_CIF="$(dirname "$0")/../Tools/EF_Hand/calmodulin_ef_hand.cif"
if [[ -f "$INPUT_CIF" ]]; then
    ok "Input CIF: $INPUT_CIF"
else
    fail "Input CIF not found: $INPUT_CIF"
fi

# ── Summary ─────────────────────────────────────────────────────────────────────
echo ""
echo "────────────────────────────────────────"
echo "  Passed: $PASS   Failed: $FAIL"
if [[ "$FAIL" -eq 0 ]]; then
    echo "  All checks passed — safe to submit SLURM job."
else
    echo "  Fix the FAIL items above before submitting."
fi
echo "────────────────────────────────────────"
