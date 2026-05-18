#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_pdb_to_stl.sh  —  One-stop runner for pdb_to_stl.py
#
# Uses the "haddock" conda env (has biopython, numpy, scipy).
# Installs scikit-image into that env on first run if missing.
#
# Usage (run from project root):
#   bash 3D-Print/run_pdb_to_stl.sh                      # all PDB + CIF in Data/
#   bash 3D-Print/run_pdb_to_stl.sh --resolution 0.8     # finer mesh
#   bash 3D-Print/run_pdb_to_stl.sh Data/Target_Crystal_Structures/MMP2_Xray.pdb
# ---------------------------------------------------------------------------

CONDA_ENV="haddock"
CONDA_BASE="$HOME/miniconda3"
PYTHON="$CONDA_BASE/envs/$CONDA_ENV/bin/python"
PIP="$CONDA_BASE/envs/$CONDA_ENV/bin/pip"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo "=== pdb_to_stl  (env: $CONDA_ENV) ==="

# ---------- install scikit-image if missing ----------
if ! "$PYTHON" -c "import skimage" 2>/dev/null; then
    echo "[setup] Installing scikit-image into '$CONDA_ENV' env ..."
    "$PIP" install scikit-image --quiet
fi

# ---------- build file list ----------
if [ "$#" -eq 0 ] || [[ "$1" == --* ]]; then
    # Default: all PDB and CIF files under Data/
    INPUTS=(
        "$PROJECT_ROOT/Data"
    )
    FIND_ARGS=(-type f \( -name "*.pdb" -o -name "*.cif" \))
    FILES=()
    while IFS= read -r f; do FILES+=("$f"); done < <(find "${INPUTS[@]}" "${FIND_ARGS[@]}")
    echo "Found ${#FILES[@]} structure file(s) in Data/"
    EXTRA_ARGS=("$@")
else
    FILES=("$@")
    EXTRA_ARGS=()
fi

if [ "${#FILES[@]}" -eq 0 ]; then
    echo "No PDB/CIF files found. Exiting."
    exit 1
fi

# ---------- run converter ----------
"$PYTHON" "$SCRIPT_DIR/pdb_to_stl.py" \
    "${FILES[@]}" \
    -o "$PROJECT_ROOT/STL_Output" \
    "${EXTRA_ARGS[@]}"
